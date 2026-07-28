"""Bounded-memory, resumable search of the indexed catalyst population."""

from __future__ import annotations

import hashlib
import json
import multiprocessing as mp
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional

import numpy as np

from pipeline.search.discovery import candidate_id, canonicalize_genome, discovery_region
from pipeline.search.indexed_space import (
    ADMISSIBILITY_POLICY_VERSION, TOTAL_SIZE, candidate_at,
    is_physically_admissible,
)


@dataclass
class ScanConfig:
    application: str
    database: str
    start: int = 0
    stop: int = TOTAL_SIZE
    batch_size: int = 65536
    worker_id: int = 0
    num_workers: int = 1
    global_archive_size: int = 10000
    max_batches: Optional[int] = None
    state_id: Optional[str] = None  # independent resume cursor for branch leaves
    deadline_epoch_s: Optional[float] = None


def _connect(path: str) -> sqlite3.Connection:
    db_path = Path(path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path), timeout=60)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    conn.executescript("""
        CREATE TABLE IF NOT EXISTS scan_state (
            application TEXT NOT NULL, worker_id INTEGER NOT NULL,
            num_workers INTEGER NOT NULL, next_index INTEGER NOT NULL,
            stop_index INTEGER NOT NULL, processed INTEGER NOT NULL DEFAULT 0,
            accepted INTEGER NOT NULL DEFAULT 0, rejected INTEGER NOT NULL DEFAULT 0,
            updated_at REAL NOT NULL, PRIMARY KEY(application, worker_id, num_workers)
        );
        CREATE TABLE IF NOT EXISTS chunks (
            application TEXT NOT NULL, worker_id INTEGER NOT NULL,
            first_index INTEGER NOT NULL, last_index INTEGER NOT NULL,
            processed INTEGER NOT NULL, accepted INTEGER NOT NULL,
            rejected INTEGER NOT NULL, identity_digest TEXT NOT NULL,
            elapsed_s REAL NOT NULL,
            PRIMARY KEY(application, worker_id, first_index)
        );
        CREATE TABLE IF NOT EXISTS region_champions (
            application TEXT NOT NULL, region TEXT NOT NULL,
            candidate_id TEXT NOT NULL, global_index INTEGER NOT NULL,
            primary_score REAL NOT NULL, objectives TEXT NOT NULL,
            genome TEXT NOT NULL, PRIMARY KEY(application, region)
        );
        CREATE TABLE IF NOT EXISTS global_archive (
            application TEXT NOT NULL, candidate_id TEXT NOT NULL,
            global_index INTEGER NOT NULL, primary_score REAL NOT NULL,
            objectives TEXT NOT NULL, genome TEXT NOT NULL,
            PRIMARY KEY(application, candidate_id)
        );
        CREATE INDEX IF NOT EXISTS global_score_idx
            ON global_archive(application, primary_score);
        CREATE TABLE IF NOT EXISTS objective_archive (
            application TEXT NOT NULL, objective_index INTEGER NOT NULL,
            candidate_id TEXT NOT NULL, global_index INTEGER NOT NULL,
            objective_score REAL NOT NULL, objectives TEXT NOT NULL,
            genome TEXT NOT NULL,
            PRIMARY KEY(application, objective_index, candidate_id)
        );
        CREATE INDEX IF NOT EXISTS objective_score_idx
            ON objective_archive(application, objective_index, objective_score);
        CREATE TABLE IF NOT EXISTS regional_objective_champions (
            application TEXT NOT NULL, region TEXT NOT NULL,
            objective_index INTEGER NOT NULL, candidate_id TEXT NOT NULL,
            global_index INTEGER NOT NULL, objective_score REAL NOT NULL,
            objectives TEXT NOT NULL, genome TEXT NOT NULL,
            PRIMARY KEY(application, region, objective_index)
        );
        CREATE TABLE IF NOT EXISTS scan_progress (
            application TEXT NOT NULL, state_id TEXT NOT NULL,
            next_index INTEGER NOT NULL, stop_index INTEGER NOT NULL,
            processed INTEGER NOT NULL DEFAULT 0, accepted INTEGER NOT NULL DEFAULT 0,
            rejected INTEGER NOT NULL DEFAULT 0, updated_at REAL NOT NULL,
            PRIMARY KEY(application, state_id)
        );
        CREATE TABLE IF NOT EXISTS scan_chunks (
            application TEXT NOT NULL, state_id TEXT NOT NULL,
            first_index INTEGER NOT NULL, last_index INTEGER NOT NULL,
            processed INTEGER NOT NULL, accepted INTEGER NOT NULL,
            rejected INTEGER NOT NULL, identity_digest TEXT NOT NULL,
            elapsed_s REAL NOT NULL,
            PRIMARY KEY(application, state_id, first_index)
        );
        CREATE TABLE IF NOT EXISTS scan_metadata (
            key TEXT PRIMARY KEY, value TEXT NOT NULL
        );
    """)
    return conn


def _verify_policy(conn: sqlite3.Connection) -> None:
    row = conn.execute(
        "SELECT value FROM scan_metadata WHERE key='admissibility_policy_version'"
    ).fetchone()
    if row is None:
        populated = any(conn.execute(f"SELECT EXISTS(SELECT 1 FROM {table} LIMIT 1)").fetchone()[0]
                        for table in ('scan_progress', 'scan_chunks', 'global_archive'))
        if populated:
            raise RuntimeError(
                "Legacy scan database has no admissibility-policy version; use a "
                "fresh database so old and canonicalized evidence are not mixed."
            )
        conn.execute(
            "INSERT INTO scan_metadata VALUES ('admissibility_policy_version', ?)",
            (ADMISSIBILITY_POLICY_VERSION,))
        conn.commit()
    elif row[0] != ADMISSIBILITY_POLICY_VERSION:
        raise RuntimeError(
            f"Scan database uses admissibility policy {row[0]!r}, expected "
            f"{ADMISSIBILITY_POLICY_VERSION!r}; use a fresh database."
        )


def _resume_index(conn: sqlite3.Connection, cfg: ScanConfig) -> int:
    state_id = cfg.state_id or f"worker:{cfg.worker_id}/{cfg.num_workers}"
    row = conn.execute(
        "SELECT next_index FROM scan_progress WHERE application=? AND state_id=?",
        (cfg.application, state_id),
    ).fetchone()
    return max(cfg.start + cfg.worker_id, int(row[0])) if row else cfg.start + cfg.worker_id


def _candidate_record(global_index: int, genome: tuple,
                      objectives: np.ndarray) -> tuple:
    region = "|".join(discovery_region(genome))
    primary = float(objectives[0])
    return canonicalize_genome(genome), region, global_index, primary, objectives, genome


def _upsert_region(conn, cfg: ScanConfig, record: tuple) -> None:
    _, region, global_index, primary, objectives, genome = record
    cid = candidate_id(genome)
    obj_json = json.dumps([float(x) for x in objectives], separators=(",", ":"))
    genome_json = json.dumps(genome, separators=(",", ":"))
    conn.execute("""
        INSERT INTO region_champions VALUES (?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(application, region) DO UPDATE SET
          candidate_id=excluded.candidate_id, global_index=excluded.global_index,
          primary_score=excluded.primary_score, objectives=excluded.objectives,
          genome=excluded.genome
        WHERE excluded.primary_score < region_champions.primary_score
    """, (cfg.application, region, cid, global_index, primary, obj_json, genome_json))


def _upsert_global(conn, cfg: ScanConfig, record: tuple) -> None:
    _, _, global_index, primary, objectives, genome = record
    cid = candidate_id(genome)
    obj_json = json.dumps([float(x) for x in objectives], separators=(",", ":"))
    genome_json = json.dumps(genome, separators=(",", ":"))
    conn.execute("""
        INSERT INTO global_archive VALUES (?, ?, ?, ?, ?, ?)
        ON CONFLICT(application, candidate_id) DO UPDATE SET
          global_index=excluded.global_index, primary_score=excluded.primary_score,
          objectives=excluded.objectives, genome=excluded.genome
        WHERE excluded.primary_score < global_archive.primary_score
    """, (cfg.application, cid, global_index, primary, obj_json, genome_json))


def _upsert_objective(conn, cfg: ScanConfig, record: tuple,
                      objective_index: int, regional: bool) -> None:
    _, region, global_index, _, objectives, genome = record
    cid = candidate_id(genome)
    score = float(objectives[objective_index])
    obj_json = json.dumps([float(x) for x in objectives], separators=(",", ":"))
    genome_json = json.dumps(genome, separators=(",", ":"))
    if regional:
        conn.execute("""INSERT INTO regional_objective_champions VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, region, objective_index) DO UPDATE SET
              candidate_id=excluded.candidate_id, global_index=excluded.global_index,
              objective_score=excluded.objective_score, objectives=excluded.objectives,
              genome=excluded.genome
            WHERE excluded.objective_score < regional_objective_champions.objective_score
        """, (cfg.application, region, objective_index, cid, global_index,
                score, obj_json, genome_json))
    else:
        conn.execute("""INSERT INTO objective_archive VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, objective_index, candidate_id) DO UPDATE SET
              global_index=excluded.global_index, objective_score=excluded.objective_score,
              objectives=excluded.objectives, genome=excluded.genome
            WHERE excluded.objective_score < objective_archive.objective_score
        """, (cfg.application, objective_index, cid, global_index,
                score, obj_json, genome_json))


def run_streaming_scan(config: ScanConfig,
                       scorer: Callable[[List[tuple]], np.ndarray]) -> dict:
    """Score a complete or partial global range and resume safely after exits.

    ``scorer`` must return an ``(N, M)`` minimization-objective array.  The
    database stores only global and per-region champions plus a cryptographic
    digest for each processed chunk, keeping storage independent of 21.1B N.
    """
    if not 0 <= config.start <= config.stop <= TOTAL_SIZE:
        raise ValueError("scan bounds outside indexed space")
    if config.batch_size <= 0:
        raise ValueError("batch_size must be positive")
    conn = _connect(config.database)
    _verify_policy(conn)
    next_index = _resume_index(conn, config)
    batches = processed_total = accepted_total = rejected_total = 0
    started = time.time()
    # Exact bounded top-K state lives in memory; only entrants are written.
    global_best = {}
    for row in conn.execute(
            "SELECT candidate_id, global_index, primary_score, objectives, genome "
            "FROM global_archive WHERE application=?", (config.application,)
        ):
        genome = json.loads(row[4])
        if genome[0] == 'SolidCatalyst': genome[5] = tuple(genome[5])
        if genome[0] == 'HEA': genome[1] = tuple(genome[1])
        genome = tuple(genome)
        key = canonicalize_genome(genome)
        global_best[key] = (key, "", row[1], row[2], np.asarray(json.loads(row[3])), genome)

    while next_index < config.stop:
        if config.deadline_epoch_s is not None and time.time() >= config.deadline_epoch_s:
            break
        if config.max_batches is not None and batches >= config.max_batches:
            break
        indices = list(range(next_index, min(config.stop, next_index + config.batch_size * config.num_workers),
                             config.num_workers))
        genomes, accepted_indices, rejected = [], [], 0
        digest = hashlib.sha256()
        for index in indices:
            genome = candidate_at(index)
            # The index-to-genome mapping is deterministic and versioned in
            # code; hashing compact indices avoids JSON+SHA work for billions
            # of candidates while still detecting missing/reordered chunks.
            digest.update(index.to_bytes(8, byteorder="little", signed=False))
            admissible, _ = is_physically_admissible(genome)
            if admissible:
                genomes.append(genome)
                accepted_indices.append(index)
            else:
                rejected += 1

        t_batch = time.time()
        if genomes:
            objectives = np.asarray(scorer(genomes), dtype=float)
            if objectives.ndim != 2 or len(objectives) != len(genomes):
                raise ValueError("scorer must return (N, M) objectives")
            if not np.all(np.isfinite(objectives)):
                raise ValueError("scorer returned NaN or infinite objectives")
            region_best = {}
            regional_objective_best = {}
            changed_global = set()
            for index, genome, obj in zip(accepted_indices, genomes, objectives):
                record = _candidate_record(index, genome, obj)
                key, region, _, primary, _, _ = record
                previous_region = region_best.get(region)
                if previous_region is None or primary < previous_region[3]:
                    region_best[region] = record
                for objective_index, score in enumerate(obj):
                    key_region = (region, objective_index)
                    prior = regional_objective_best.get(key_region)
                    if prior is None or float(score) < float(prior[4][objective_index]):
                        regional_objective_best[key_region] = record
                previous_global = global_best.get(key)
                if previous_global is None or primary < previous_global[3]:
                    global_best[key] = record
                    changed_global.add(key)

            # Reduce in memory before touching SQLite. This is exact for the
            # primary-objective top-K and turns billions of potential writes
            # into at most region winners plus actual archive entrants.
            if len(global_best) > config.global_archive_size:
                # Candidate ID is the stable total-order tie break used by the
                # SQLite archive as well; results must not depend on batch size
                # or whether indices were processed serially or by shards.
                keep = sorted(
                    global_best.values(),
                    key=lambda record: (record[3], candidate_id(record[5]))
                )[:config.global_archive_size]
                global_best = {r[0]: r for r in keep}
            for record in region_best.values():
                _upsert_region(conn, config, record)
            for (_, objective_index), record in regional_objective_best.items():
                _upsert_objective(conn, config, record, objective_index, regional=True)
            for key in changed_global & global_best.keys():
                _upsert_global(conn, config, global_best[key])
            # Preserve strong candidates for every objective, not only objective
            # zero. Batch reduction keeps database traffic bounded.
            per_objective_limit = max(1, config.global_archive_size // objectives.shape[1])
            for objective_index in range(objectives.shape[1]):
                # Reduce every batch with the same total ordering used by the
                # database. This makes tied objectives invariant to batches and
                # worker shards while keeping writes bounded by archive size.
                best_local = sorted(
                    range(len(genomes)),
                    key=lambda i: (float(objectives[i, objective_index]),
                                   candidate_id(genomes[i]))
                )[:per_objective_limit]
                archive_count = conn.execute(
                    "SELECT COUNT(*) FROM objective_archive "
                    "WHERE application=? AND objective_index=?",
                    (config.application, objective_index)).fetchone()[0]
                if archive_count >= per_objective_limit:
                    worst_score, worst_id = conn.execute(
                        "SELECT objective_score, candidate_id FROM objective_archive "
                        "WHERE application=? AND objective_index=? "
                        "ORDER BY objective_score DESC, candidate_id DESC LIMIT 1",
                        (config.application, objective_index)).fetchone()
                    best_local = [
                        i for i in best_local
                        if (float(objectives[i, objective_index]),
                            candidate_id(genomes[i])) < (worst_score, worst_id)
                    ]
                for local_index in best_local:
                    record = _candidate_record(accepted_indices[local_index], genomes[local_index], objectives[local_index])
                    _upsert_objective(conn, config, record, objective_index, regional=False)
                excess_obj = conn.execute(
                    "SELECT COUNT(*)-? FROM objective_archive WHERE application=? AND objective_index=?",
                    (per_objective_limit, config.application, objective_index)).fetchone()[0]
                if excess_obj > 0:
                    conn.execute("""DELETE FROM objective_archive WHERE rowid IN (
                        SELECT rowid FROM objective_archive WHERE application=? AND objective_index=?
                        ORDER BY objective_score DESC, candidate_id DESC LIMIT ?)
                    """, (config.application, objective_index, excess_obj))

        # Bound the global archive. Region champions remain independent, so
        # unfamiliar families are retained even when absent from global top-K.
        excess = conn.execute(
            "SELECT COUNT(*)-? FROM global_archive WHERE application=?",
            (config.global_archive_size, config.application),
        ).fetchone()[0]
        if excess > 0:
            conn.execute("""DELETE FROM global_archive WHERE rowid IN (
                SELECT rowid FROM global_archive WHERE application=?
                ORDER BY primary_score DESC, candidate_id DESC LIMIT ?)
            """, (config.application, excess))

        processed, accepted = len(indices), len(genomes)
        following = indices[-1] + config.num_workers if indices else config.stop
        state_id = config.state_id or f"worker:{config.worker_id}/{config.num_workers}"
        conn.execute("INSERT OR REPLACE INTO scan_chunks VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)", (
            config.application, state_id, indices[0], indices[-1], processed,
            accepted, rejected, digest.hexdigest(), time.time() - t_batch))
        conn.execute("""INSERT INTO scan_progress VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, state_id) DO UPDATE SET
              next_index=excluded.next_index, stop_index=excluded.stop_index,
              processed=scan_progress.processed+excluded.processed,
              accepted=scan_progress.accepted+excluded.accepted,
              rejected=scan_progress.rejected+excluded.rejected,
              updated_at=excluded.updated_at
        """, (config.application, state_id, following, config.stop, processed,
                accepted, rejected, time.time()))
        conn.commit()
        next_index = following
        processed_total += processed
        accepted_total += accepted
        rejected_total += rejected
        batches += 1

    region_count = conn.execute(
        "SELECT COUNT(*) FROM region_champions WHERE application=?", (config.application,)
    ).fetchone()[0]
    global_count = conn.execute(
        "SELECT COUNT(*) FROM global_archive WHERE application=?", (config.application,)
    ).fetchone()[0]
    conn.close()
    return {
        "application": config.application,
        "processed_this_run": processed_total,
        "accepted_this_run": accepted_total,
        "rejected_this_run": rejected_total,
        "next_index": next_index,
        "stop_index": config.stop,
        "complete": next_index >= config.stop,
        "region_champions": region_count,
        "global_archive": global_count,
        "elapsed_s": time.time() - started,
    }


_SHARD_SCORER = None


def _run_scan_shard(config: ScanConfig) -> dict:
    if _SHARD_SCORER is None:
        raise RuntimeError('sharded scanner has no inherited scorer')
    return run_streaming_scan(config, _SHARD_SCORER)


def _merge_shard_archives(config: ScanConfig, shard_databases: list[Path],
                          aggregate_state_id: str) -> None:
    """Merge complete independent shard archives into the campaign database."""
    destination = _connect(config.database)
    _verify_policy(destination)
    table_specs = {
        'region_champions': (
            7, """INSERT INTO region_champions VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, region) DO UPDATE SET
              candidate_id=excluded.candidate_id, global_index=excluded.global_index,
              primary_score=excluded.primary_score, objectives=excluded.objectives,
              genome=excluded.genome
            WHERE excluded.primary_score < region_champions.primary_score OR
              (excluded.primary_score = region_champions.primary_score AND
               excluded.global_index < region_champions.global_index)"""),
        'global_archive': (
            6, """INSERT INTO global_archive VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, candidate_id) DO UPDATE SET
              global_index=excluded.global_index, primary_score=excluded.primary_score,
              objectives=excluded.objectives, genome=excluded.genome
            WHERE excluded.primary_score < global_archive.primary_score OR
              (excluded.primary_score = global_archive.primary_score AND
               excluded.global_index < global_archive.global_index)"""),
        'objective_archive': (
            7, """INSERT INTO objective_archive VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, objective_index, candidate_id) DO UPDATE SET
              global_index=excluded.global_index, objective_score=excluded.objective_score,
              objectives=excluded.objectives, genome=excluded.genome
            WHERE excluded.objective_score < objective_archive.objective_score OR
              (excluded.objective_score = objective_archive.objective_score AND
               excluded.global_index < objective_archive.global_index)"""),
        'regional_objective_champions': (
            8, """INSERT INTO regional_objective_champions VALUES (?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(application, region, objective_index) DO UPDATE SET
              candidate_id=excluded.candidate_id, global_index=excluded.global_index,
              objective_score=excluded.objective_score, objectives=excluded.objectives,
              genome=excluded.genome
            WHERE excluded.objective_score < regional_objective_champions.objective_score OR
              (excluded.objective_score = regional_objective_champions.objective_score AND
               excluded.global_index < regional_objective_champions.global_index)"""),
    }
    processed = accepted = rejected = 0
    digest = hashlib.sha256()
    for shard_path in sorted(shard_databases, key=str):
        with _connect(str(shard_path)) as source:
            _verify_policy(source)
            progress = source.execute(
                "SELECT processed, accepted, rejected FROM scan_progress "
                "WHERE application=?", (config.application,)).fetchone()
            if progress is None:
                destination.close()
                raise RuntimeError(f'shard lacks progress: {shard_path}')
            processed += int(progress[0])
            accepted += int(progress[1])
            rejected += int(progress[2])
            for row in source.execute(
                    "SELECT identity_digest FROM scan_chunks WHERE application=? "
                    "ORDER BY state_id, first_index", (config.application,)):
                digest.update(row[0].encode())
            for table, (width, statement) in table_specs.items():
                rows = source.execute(
                    f"SELECT * FROM {table} WHERE application=?",
                    (config.application,)).fetchall()
                if any(len(row) != width for row in rows):
                    destination.close()
                    raise RuntimeError(f'invalid {table} row width in {shard_path}')
                destination.executemany(statement, rows)

    expected = config.stop - config.start
    if processed != expected or accepted + rejected != processed:
        destination.close()
        raise RuntimeError(
            f'shard coverage mismatch: processed={processed}, expected={expected}, '
            f'accepted={accepted}, rejected={rejected}')

    excess = destination.execute(
        "SELECT COUNT(*)-? FROM global_archive WHERE application=?",
        (config.global_archive_size, config.application)).fetchone()[0]
    if excess > 0:
        destination.execute("""DELETE FROM global_archive WHERE rowid IN (
            SELECT rowid FROM global_archive WHERE application=?
            ORDER BY primary_score DESC, candidate_id DESC LIMIT ?)""",
            (config.application, excess))
    objective_indices = [
        row[0] for row in destination.execute(
            "SELECT DISTINCT objective_index FROM objective_archive WHERE application=?",
            (config.application,))]
    for objective_index in objective_indices:
        objective_count = destination.execute(
            "SELECT COUNT(*) FROM objective_archive WHERE application=? "
            "AND objective_index=?", (config.application, objective_index)).fetchone()[0]
        objective_limit = max(1, config.global_archive_size // max(1, len(objective_indices)))
        excess = objective_count - objective_limit
        if excess > 0:
            destination.execute("""DELETE FROM objective_archive WHERE rowid IN (
                SELECT rowid FROM objective_archive WHERE application=? AND objective_index=?
                ORDER BY objective_score DESC, candidate_id DESC LIMIT ?)""",
                (config.application, objective_index, excess))

    destination.execute(
        "INSERT OR REPLACE INTO scan_chunks VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
        (config.application, aggregate_state_id, config.start, config.stop - 1,
         processed, accepted, rejected, digest.hexdigest(), 0.0))
    destination.execute("""INSERT INTO scan_progress VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(application, state_id) DO UPDATE SET
          next_index=excluded.next_index, stop_index=excluded.stop_index,
          processed=excluded.processed, accepted=excluded.accepted,
          rejected=excluded.rejected, updated_at=excluded.updated_at""",
        (config.application, aggregate_state_id, config.stop, config.stop,
         processed, accepted, rejected, time.time()))
    destination.commit()
    destination.close()


def run_sharded_scan(config: ScanConfig,
                     scorer: Callable[[List[tuple]], np.ndarray],
                     workers: int) -> dict:
    """Run one range through independent deterministic process/database shards."""
    if workers < 1:
        raise ValueError('workers must be positive')
    range_size = config.stop - config.start
    effective_workers = min(workers, range_size)
    if effective_workers <= 1:
        return run_streaming_scan(config, scorer)
    if config.worker_id != 0 or config.num_workers != 1:
        raise ValueError('sharded scan expects an unsharded base configuration')
    if 'fork' not in mp.get_all_start_methods():
        raise RuntimeError('process sharding requires fork semantics for the fitted scorer')

    state_id = config.state_id or f'range:{config.start}:{config.stop}'
    shard_key = hashlib.sha256(
        f'{config.application}:{state_id}:{config.start}:{config.stop}:{effective_workers}'.encode()
    ).hexdigest()[:20]
    shard_root = Path(f'{config.database}.shards') / shard_key
    shard_root.mkdir(parents=True, exist_ok=True)
    shard_configs, shard_databases = [], []
    for worker_id in range(effective_workers):
        database = shard_root / f'worker_{worker_id:03d}.sqlite'
        shard_databases.append(database)
        shard_configs.append(ScanConfig(
            application=config.application, database=str(database),
            start=config.start, stop=config.stop, batch_size=config.batch_size,
            worker_id=worker_id, num_workers=effective_workers,
            global_archive_size=config.global_archive_size,
            max_batches=config.max_batches,
            state_id=f'{state_id}:shard:{worker_id}/{effective_workers}',
            deadline_epoch_s=config.deadline_epoch_s,
        ))

    global _SHARD_SCORER
    _SHARD_SCORER = scorer
    context = mp.get_context('fork')
    try:
        with context.Pool(processes=effective_workers) as pool:
            summaries = pool.map(_run_scan_shard, shard_configs)
    finally:
        _SHARD_SCORER = None

    complete = all(summary['complete'] for summary in summaries)
    processed = sum(summary['processed_this_run'] for summary in summaries)
    if complete:
        _merge_shard_archives(config, shard_databases, state_id)
    return {
        'application': config.application,
        'workers': effective_workers,
        'complete': complete,
        'processed_this_run': processed,
        'shards': summaries,
        'shard_root': str(shard_root),
        'next_index': config.stop if complete else None,
        'stop_index': config.stop,
    }


def load_archive_genomes(database: str, application: str,
                         limit: int = 10000) -> List[tuple]:
    """Load deduplicated global and region champions for downstream search."""
    conn = _connect(database)
    rows = conn.execute("""
        SELECT genome, primary_score FROM global_archive WHERE application=?
        UNION ALL
        SELECT genome, primary_score FROM region_champions WHERE application=?
        UNION ALL
        SELECT genome, objective_score FROM objective_archive WHERE application=?
        UNION ALL
        SELECT genome, objective_score FROM regional_objective_champions WHERE application=?
        ORDER BY primary_score ASC LIMIT ?
    """, (application, application, application, application, limit * 8)).fetchall()
    conn.close()
    result, seen = [], set()
    for raw, _ in rows:
        genome = tuple(json.loads(raw))
        # Restore nested tuple fields used by encoders.
        if genome[0] in ("SolidCatalyst", "HEA"):
            genome = tuple(list(genome[:1]) + [tuple(genome[1])] + list(genome[2:])) if genome[0] == "HEA" else tuple(list(genome[:5]) + [tuple(genome[5])] + list(genome[6:]))
        cid = candidate_id(genome)
        if cid not in seen:
            result.append(genome)
            seen.add(cid)
            if len(result) >= limit:
                break
    return result
