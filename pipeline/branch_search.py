"""Safe deterministic divide-and-conquer over the indexed catalyst space.

Surrogate probes determine processing order, never exclusion.  A node is pruned
only when every encoded member has been checked by conservative hard constraints.
All other leaves are delegated to the exhaustive scanner.
"""

from __future__ import annotations

import hashlib
import json
import math
import sqlite3
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, List, Optional

import numpy as np

from pipeline.exhaustive_search import ScanConfig, run_streaming_scan
from pipeline.indexed_space import (CLASS_OFFSETS, CLASS_ORDER, CLASS_SIZES,
                                    TOTAL_SIZE, candidate_at, is_physically_admissible)
from pipeline.ood_detector import CLASS_CONFIDENCE


@dataclass
class BranchConfig:
    application: str
    database: str
    leaf_size: int = 1_000_000
    probe_count: int = 9
    scan_batch_size: int = 65536
    hard_prune_limit: int = 4096
    max_leaves: Optional[int] = None
    global_archive_size: int = 10000
    material_classes: Optional[tuple] = None
    expected_population: Optional[int] = None
    certificate_path: Optional[str] = None
    max_runtime_s: Optional[float] = None


def _node_id(application: str, start: int, stop: int) -> str:
    return hashlib.sha1(f"{application}:{start}:{stop}".encode()).hexdigest()[:20]


def _open(path: str) -> sqlite3.Connection:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(path, timeout=60)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("""CREATE TABLE IF NOT EXISTS branch_nodes (
        application TEXT NOT NULL, node_id TEXT NOT NULL, material_class TEXT NOT NULL,
        start_index INTEGER NOT NULL, stop_index INTEGER NOT NULL, depth INTEGER NOT NULL,
        status TEXT NOT NULL, priority REAL, probe_best REAL, probe_spread REAL,
        parent_id TEXT, reason TEXT, updated_at REAL NOT NULL,
        PRIMARY KEY(application, node_id)
    )""")
    return conn


def _probe_indices(start: int, stop: int, count: int) -> List[int]:
    size = stop - start
    if size <= count:
        return list(range(start, stop))
    # Integer linspace includes both boundaries and is deterministic.
    return sorted({start + (i * (size - 1)) // (count - 1) for i in range(count)})


def _node_priority(material_class: str, start: int, stop: int,
                   scorer: Callable[[List[tuple]], np.ndarray], probe_count: int):
    probes = []
    for index in _probe_indices(start, stop, max(2, probe_count)):
        genome = candidate_at(index)
        admissible, _ = is_physically_admissible(genome)
        if admissible:
            probes.append(genome)
    if not probes:
        return -1e9, None, None
    objectives = np.asarray(scorer(probes), float)
    if objectives.ndim != 2 or len(objectives) != len(probes):
        raise ValueError("branch scorer must return (N, M) objectives")
    primary = objectives[:, 0]
    best = float(primary.min())
    spread = float(np.ptp(primary))
    # Lower is processed first. Spread and low class confidence increase
    # exploration priority, but are never used as a pruning certificate.
    ood_bonus = 1.0 - CLASS_CONFIDENCE.get(material_class, 0.5)
    size_bonus = 0.01 * math.log10(max(1, stop - start))
    priority = best - 0.5 * spread - 0.2 * ood_bonus - size_bonus
    return priority, best, spread


def _insert_node(conn, cfg, material_class, start, stop, depth, scorer,
                 parent_id=None):
    node_id = _node_id(cfg.application, start, stop)
    exists = conn.execute(
        "SELECT 1 FROM branch_nodes WHERE application=? AND node_id=?",
        (cfg.application, node_id)).fetchone()
    if exists:
        return node_id
    priority, best, spread = _node_priority(
        material_class, start, stop, scorer, cfg.probe_count)
    conn.execute("INSERT INTO branch_nodes VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", (
        cfg.application, node_id, material_class, start, stop, depth, 'pending',
        priority, best, spread, parent_id, None, time.time()))
    return node_id


def run_branch_and_bound(config: BranchConfig,
                         scorer: Callable[[List[tuple]], np.ndarray]) -> dict:
    """Recursively schedule and exhaustively resolve catalyst-space leaves."""
    if config.leaf_size <= 0 or config.probe_count < 2:
        raise ValueError("invalid branch configuration")
    conn = _open(config.database)
    classes = config.material_classes or CLASS_ORDER
    unknown = set(classes) - set(CLASS_ORDER)
    if unknown:
        raise ValueError(f"unknown material classes: {sorted(unknown)}")
    declared_population = sum(CLASS_SIZES[c] for c in classes)
    expected = config.expected_population
    if expected is None and tuple(classes) == tuple(CLASS_ORDER):
        expected = TOTAL_SIZE
    if expected is not None and expected != declared_population:
        raise ValueError(
            f"population denominator mismatch: expected {expected:,}, "
            f"indexed space contains {declared_population:,}"
        )
    for cls in classes:
        start = CLASS_OFFSETS[cls]
        _insert_node(conn, config, cls, start, start + CLASS_SIZES[cls], 0, scorer)
    conn.commit()

    leaves_scanned = expanded = pruned = 0
    deadline = None if config.max_runtime_s is None else time.time() + config.max_runtime_s
    while config.max_leaves is None or leaves_scanned < config.max_leaves:
        if deadline is not None and time.time() >= deadline:
            break
        row = conn.execute("""SELECT node_id, material_class, start_index, stop_index, depth
            FROM branch_nodes WHERE application=? AND status='pending'
            ORDER BY priority ASC, start_index ASC LIMIT 1""",
            (config.application,)).fetchone()
        if row is None:
            break
        node_id, cls, start, stop, depth = row
        size = stop - start

        # A hard prune is a proof over every member, never a surrogate guess.
        if size <= config.hard_prune_limit:
            any_admissible = any(is_physically_admissible(candidate_at(i))[0]
                                 for i in range(start, stop))
            if not any_admissible:
                conn.execute("UPDATE branch_nodes SET status='pruned', reason=?, updated_at=? "
                             "WHERE application=? AND node_id=?",
                             ('all_members_fail_hard_constraints', time.time(),
                              config.application, node_id))
                conn.commit()
                pruned += 1
                continue

        if size > config.leaf_size:
            mid = start + size // 2
            _insert_node(conn, config, cls, start, mid, depth + 1, scorer, node_id)
            _insert_node(conn, config, cls, mid, stop, depth + 1, scorer, node_id)
            conn.execute("UPDATE branch_nodes SET status='expanded', updated_at=? "
                         "WHERE application=? AND node_id=?",
                         (time.time(), config.application, node_id))
            conn.commit()
            expanded += 1
            continue

        # Close our planning transaction before the scanner writes archives.
        conn.commit()
        summary = run_streaming_scan(ScanConfig(
            application=config.application, database=config.database,
            start=start, stop=stop, batch_size=config.scan_batch_size,
            global_archive_size=config.global_archive_size,
            state_id=f"branch:{node_id}",
            deadline_epoch_s=deadline,
        ), scorer)
        status = 'scanned' if summary['complete'] else 'pending'
        conn.execute("UPDATE branch_nodes SET status=?, reason=?, updated_at=? "
                     "WHERE application=? AND node_id=?",
                     (status, json.dumps(summary, sort_keys=True), time.time(),
                      config.application, node_id))
        conn.commit()
        leaves_scanned += 1

    counts = dict(conn.execute(
        "SELECT status, COUNT(*) FROM branch_nodes WHERE application=? GROUP BY status",
        (config.application,)).fetchall())
    unresolved_population = conn.execute(
        "SELECT COALESCE(SUM(stop_index-start_index),0) FROM branch_nodes "
        "WHERE application=? AND status='pending'", (config.application,)).fetchone()[0]
    conn.close()
    certificate = verify_branch_coverage(
        config.database, config.application, tuple(classes), config.certificate_path)
    return {
        "application": config.application,
        "leaves_scanned_this_run": leaves_scanned,
        "nodes_expanded_this_run": expanded,
        "nodes_pruned_this_run": pruned,
        "node_status_counts": counts,
        "unresolved_encoded_population": int(unresolved_population),
        "complete": certificate['complete'],
        "coverage_certificate": certificate,
    }


def verify_branch_coverage(database: str, application: str,
                           material_classes: Optional[tuple] = None,
                           certificate_path: Optional[str] = None) -> dict:
    """Prove that terminal nodes partition the declared space exactly.

    This verifies address coverage, scanner completion, and pruning provenance;
    it does not assert that surrogate predictions equal experimental truth.
    """
    classes = material_classes or CLASS_ORDER
    conn = _open(database)
    errors = []
    terminal_count = resolved_count = unresolved_count = covered = 0
    digest = hashlib.sha256()

    for cls in classes:
        expected_start = CLASS_OFFSETS[cls]
        class_stop = expected_start + CLASS_SIZES[cls]
        rows = conn.execute("""SELECT node_id, start_index, stop_index, status, reason
            FROM branch_nodes WHERE application=? AND material_class=?
              AND status!='expanded' ORDER BY start_index, stop_index""",
            (application, cls)).fetchall()
        cursor = expected_start
        for node_id, start, stop, status, reason in rows:
            terminal_count += 1
            if start != cursor:
                kind = 'overlap' if start < cursor else 'gap'
                errors.append(f"{cls}: {kind} before [{start},{stop}); expected {cursor}")
            if stop <= start or stop > class_stop:
                errors.append(f"{cls}: invalid terminal interval [{start},{stop})")
            cursor = max(cursor, stop)
            covered += max(0, stop - start)
            digest.update(f"{cls}:{start}:{stop}:{status}".encode())

            if status == 'scanned':
                progress = conn.execute("""SELECT next_index, stop_index FROM scan_progress
                    WHERE application=? AND state_id=?""",
                    (application, f"branch:{node_id}")).fetchone()
                if not progress or progress[0] < stop or progress[1] != stop:
                    errors.append(f"{cls}: scanned leaf {node_id} lacks complete scan cursor")
                else:
                    resolved_count += 1
            elif status == 'pruned':
                if reason != 'all_members_fail_hard_constraints':
                    errors.append(f"{cls}: pruned leaf {node_id} lacks hard proof")
                else:
                    # Recheck the proof when the certificate is generated.
                    if any(is_physically_admissible(candidate_at(i))[0]
                           for i in range(start, stop)):
                        errors.append(f"{cls}: invalid hard-prune proof for {node_id}")
                    else:
                        resolved_count += 1
            elif status == 'pending':
                unresolved_count += 1
            else:
                errors.append(f"{cls}: unknown terminal status {status}")
        if cursor != class_stop:
            errors.append(f"{cls}: terminal coverage ends at {cursor}, expected {class_stop}")

    declared = sum(CLASS_SIZES[c] for c in classes)
    if covered != declared:
        errors.append(f"terminal interval sum {covered} != declared population {declared}")
    certificate = {
        'application': application,
        'material_classes': list(classes),
        'declared_encoded_population': declared,
        'terminal_interval_population': covered,
        'terminal_nodes': terminal_count,
        'resolved_terminal_nodes': resolved_count,
        'unresolved_terminal_nodes': unresolved_count,
        'partition_digest_sha256': digest.hexdigest(),
        'gap_free': not any('gap' in e for e in errors),
        'overlap_free': not any('overlap' in e for e in errors),
        'complete': not errors and unresolved_count == 0,
        'errors': errors,
        'generated_at': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
    }
    conn.close()
    if certificate_path:
        path = Path(certificate_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(certificate, indent=2, sort_keys=True) + '\n')
    return certificate
