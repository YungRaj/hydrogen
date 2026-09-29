"""Typed, content-addressed reuse of completed atomistic screening results."""

from __future__ import annotations

import hashlib
import importlib.metadata
import inspect
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Callable, Mapping, Sequence, TypeAlias

from pipeline.search.discovery import candidate_id


_SCIENTIFIC_MODULES = (
    'pipeline.screening.surface_screener',
    'pipeline.screening.protocols',
    'pipeline.screening.relaxation',
    'pipeline.screening.surface_calculator',
    'pipeline.screening.batched_calculator',
)
_FIELD_TYPES = frozenset({'none', 'bool', 'int', 'float', 'str'})
ScalarValue: TypeAlias = None | bool | int | float | str
CachedRecord: TypeAlias = dict[str, ScalarValue]


def implementation_digest(worker_target: Callable) -> str:
    """Hash scientific source and dependency versions affecting a result.

    Args:
        worker_target: Application worker containing the descriptor equations.

    Returns:
        SHA-256 digest identifying the exact reusable implementation.
    """
    import importlib

    digest = hashlib.sha256()
    modules = tuple(dict.fromkeys(
        (worker_target.__module__, *_SCIENTIFIC_MODULES)))
    for name in modules:
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
    connection.executescript('''
        CREATE TABLE IF NOT EXISTS screening_cache_entries (
            application TEXT NOT NULL,
            protocol_id TEXT NOT NULL,
            implementation_digest TEXT NOT NULL,
            candidate_id TEXT NOT NULL,
            created_at REAL NOT NULL,
            PRIMARY KEY(application, protocol_id, implementation_digest, candidate_id)
        );
        CREATE TABLE IF NOT EXISTS screening_cache_fields (
            application TEXT NOT NULL,
            protocol_id TEXT NOT NULL,
            implementation_digest TEXT NOT NULL,
            candidate_id TEXT NOT NULL,
            field_name TEXT NOT NULL,
            field_type TEXT NOT NULL CHECK(field_type IN ('none','bool','int','float','str')),
            bool_value INTEGER,
            int_value INTEGER,
            float_value TEXT,
            str_value TEXT,
            PRIMARY KEY(application, protocol_id, implementation_digest,
                        candidate_id, field_name),
            FOREIGN KEY(application, protocol_id, implementation_digest, candidate_id)
                REFERENCES screening_cache_entries(
                    application, protocol_id, implementation_digest, candidate_id)
                ON DELETE CASCADE
        );
    ''')
    connection.execute('PRAGMA foreign_keys=ON')
    return connection


def _reusable(record: Mapping[str, object], protocol_id: str) -> bool:
    if record.get('screening_protocol') != protocol_id or record.get('valid') is not True:
        return False
    convergence = [value for key, value in record.items()
                   if key.startswith('relax_') and key.endswith('_converged')]
    return bool(convergence) and all(bool(value) for value in convergence)


def _encode_field(
        value: object,
        ) -> tuple[str, int | None, int | None, str | None, str | None]:
    if hasattr(value, 'item'):
        value = value.item()
    if value is None:
        return 'none', None, None, None, None
    if type(value) is bool:
        return 'bool', int(value), None, None, None
    if type(value) is int:
        return 'int', None, value, None, None
    if type(value) is float:
        return 'float', None, None, value.hex(), None
    if type(value) is str:
        return 'str', None, None, None, value
    raise TypeError(f'unsupported cache field type: {type(value).__name__}')


def _decode_field(row: tuple[object, ...]) -> ScalarValue:
    field_type, bool_value, int_value, float_value, str_value = row
    if field_type not in _FIELD_TYPES:
        raise ValueError(f'unknown cache field type: {field_type}')
    if field_type == 'none':
        return None
    if field_type == 'bool':
        return bool(bool_value)
    if field_type == 'int':
        return int(int_value)
    if field_type == 'float':
        return float.fromhex(float_value)
    return str(str_value)


def load_cached_results(path: Path, application: str, protocol_id: str,
                        digest: str,
                        genomes: Sequence[tuple]) -> dict[int, CachedRecord]:
    """Load typed reusable rows indexed by their requested input position.

    Args:
        path: SQLite cache location.
        application: Scientific application namespace.
        protocol_id: Immutable screening-protocol identifier.
        digest: Scientific implementation and dependency fingerprint.
        genomes: Ordered candidate requests to look up.

    Returns:
        Mapping from input positions to independently reconstructed records.
    """
    if not path.exists() or not genomes:
        return {}
    ids = [candidate_id(genome) for genome in genomes]
    records: dict[str, CachedRecord] = defaultdict(dict)
    unique = sorted(set(ids))
    with _connect(path) as connection:
        for offset in range(0, len(unique), 400):
            chunk = unique[offset:offset + 400]
            placeholders = ','.join('?' for _ in chunk)
            rows = connection.execute(
                f'''SELECT candidate_id, field_name, field_type, bool_value,
                           int_value, float_value, str_value
                    FROM screening_cache_fields
                    WHERE application=? AND protocol_id=?
                      AND implementation_digest=?
                      AND candidate_id IN ({placeholders})''',
                (application, protocol_id, digest, *chunk)).fetchall()
            for cid, name, *encoded in rows:
                records[cid][name] = _decode_field(tuple(encoded))
    return {index: dict(records[cid]) for index, cid in enumerate(ids)
            if cid in records and _reusable(records[cid], protocol_id)}


def store_completed_results(path: Path, application: str, protocol_id: str,
                            digest: str, genomes: Sequence[tuple],
                            results: Sequence[Mapping[str, object]]) -> int:
    """Persist fully converged scalar records in typed SQLite columns.

    Args:
        path: SQLite cache location.
        application: Scientific application namespace.
        protocol_id: Immutable screening-protocol identifier.
        digest: Scientific implementation and dependency fingerprint.
        genomes: Candidates corresponding positionally to ``results``.
        results: Screening records considered for fail-closed reuse.

    Returns:
        Number of valid, fully converged scalar records written.
    """
    if len(genomes) != len(results):
        raise ValueError('cache genomes and results must have identical lengths')
    accepted = []
    for genome, result in zip(genomes, results):
        if not _reusable(result, protocol_id):
            continue
        try:
            fields = [(name, *_encode_field(value))
                      for name, value in result.items()]
        except TypeError:
            continue
        accepted.append((candidate_id(genome), fields))
    if not accepted:
        return 0
    with _connect(path) as connection:
        for cid, fields in accepted:
            key = (application, protocol_id, digest, cid)
            connection.execute(
                'INSERT OR REPLACE INTO screening_cache_entries VALUES (?, ?, ?, ?, ?)',
                (*key, time.time()))
            connection.execute(
                '''DELETE FROM screening_cache_fields WHERE application=?
                   AND protocol_id=? AND implementation_digest=? AND candidate_id=?''', key)
            connection.executemany(
                '''INSERT INTO screening_cache_fields VALUES
                   (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)''',
                [(*key, *field) for field in fields])
    return len(accepted)
