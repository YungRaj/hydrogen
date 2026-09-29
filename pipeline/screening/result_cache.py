"""Content-addressed reuse of completed atomistic screening calculations."""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import json
import sqlite3
import sys
import time
from pathlib import Path
from typing import Callable, Sequence

from pipeline.search.discovery import candidate_id


_SCIENTIFIC_MODULES = (
    'pipeline.screening.surface_screener',
    'pipeline.screening.protocols',
    'pipeline.screening.relaxation',
    'pipeline.screening.surface_calculator',
    'pipeline.screening.batched_calculator',
)


def implementation_digest(worker_target: Callable) -> str:
    """Hash scientific source and dependency versions affecting a result.

    Args:
        worker_target: Application worker whose complete module contains the
            structure generation, reference energies, and descriptor equations.

    Returns:
        SHA-256 digest identifying the exact reusable implementation.
    """
    import importlib

    digest = hashlib.sha256()
    module_names = tuple(dict.fromkeys(
        (worker_target.__module__, *_SCIENTIFIC_MODULES)))
    for name in module_names:
        module = importlib.import_module(name)
        path = inspect.getsourcefile(module)
        if path is None:
            raise RuntimeError(f'cannot fingerprint scientific module {name}')
        digest.update(name.encode())
        digest.update(Path(path).read_bytes())
    digest.update(sys.version.encode())
    for distribution in ('ase', 'fairchem-core', 'torch'):
        try:
            version = importlib.metadata.version(distribution)
        except importlib.metadata.PackageNotFoundError:
            version = 'missing'
        digest.update(f'{distribution}={version}'.encode())
    return digest.hexdigest()


def _connect(path: Path) -> sqlite3.Connection:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path, timeout=60)
    connection.execute('PRAGMA journal_mode=WAL')
    connection.execute('PRAGMA synchronous=NORMAL')
    connection.execute('''CREATE TABLE IF NOT EXISTS screening_results (
        application TEXT NOT NULL,
        protocol_id TEXT NOT NULL,
        implementation_digest TEXT NOT NULL,
        candidate_id TEXT NOT NULL,
        result_json TEXT NOT NULL,
        created_at REAL NOT NULL,
        PRIMARY KEY(application, protocol_id, implementation_digest, candidate_id)
    )''')
    return connection


def _reusable(record: dict, protocol_id: str) -> bool:
    """Accept only complete successful records from the exact protocol."""
    if record.get('screening_protocol') != protocol_id or record.get('valid') is not True:
        return False
    convergence = [value for key, value in record.items()
                   if key.startswith('relax_') and key.endswith('_converged')]
    return bool(convergence) and all(bool(value) for value in convergence)


def load_cached_results(path: Path, application: str, protocol_id: str,
                        digest: str, genomes: Sequence[tuple]) -> dict[int, dict]:
    """Load reusable rows indexed by their requested input position.

    Args:
        path: SQLite cache location.
        application: Scientific application namespace.
        protocol_id: Immutable screening-protocol identifier.
        digest: Scientific implementation and dependency fingerprint.
        genomes: Ordered candidate requests to look up.

    Returns:
        Mapping from input positions to independently copied cached records.
    """
    if not path.exists() or not genomes:
        return {}
    ids = [candidate_id(genome) for genome in genomes]
    unique = sorted(set(ids))
    rows = []
    with _connect(path) as connection:
        for offset in range(0, len(unique), 500):
            chunk = unique[offset:offset + 500]
            placeholders = ','.join('?' for _ in chunk)
            rows.extend(connection.execute(
                f'''SELECT candidate_id, result_json FROM screening_results
                    WHERE application=? AND protocol_id=? AND implementation_digest=?
                      AND candidate_id IN ({placeholders})''',
                (application, protocol_id, digest, *chunk)).fetchall())
    records = {}
    for cid, payload in rows:
        try:
            record = json.loads(payload)
        except (json.JSONDecodeError, TypeError):
            continue
        if isinstance(record, dict):
            records[cid] = record
    return {index: dict(records[cid]) for index, cid in enumerate(ids)
            if cid in records and _reusable(records[cid], protocol_id)}


def _json_default(value):
    if hasattr(value, 'item'):
        return value.item()
    raise TypeError(f'{type(value).__name__} is not JSON serializable')


def store_completed_results(path: Path, application: str, protocol_id: str,
                            digest: str, genomes: Sequence[tuple],
                            results: Sequence[dict]) -> int:
    """Persist fully converged results and return the number accepted.

    Args:
        path: SQLite cache location.
        application: Scientific application namespace.
        protocol_id: Immutable screening-protocol identifier.
        digest: Scientific implementation and dependency fingerprint.
        genomes: Candidates corresponding positionally to ``results``.
        results: Screening records considered for fail-closed reuse.

    Returns:
        Number of valid, fully converged records written to the cache.
    """
    if len(genomes) != len(results):
        raise ValueError('cache genomes and results must have identical lengths')
    rows = []
    for genome, result in zip(genomes, results):
        if _reusable(result, protocol_id):
            rows.append((application, protocol_id, digest, candidate_id(genome),
                         json.dumps(result, sort_keys=True, separators=(',', ':'),
                                    default=_json_default), time.time()))
    if rows:
        with _connect(path) as connection:
            connection.executemany(
                'INSERT OR REPLACE INTO screening_results VALUES (?, ?, ?, ?, ?, ?)',
                rows)
    return len(rows)
