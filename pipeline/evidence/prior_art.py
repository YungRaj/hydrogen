"""Versioned evidence registry for literature, patents, and experiments."""

import ast
import csv
import hashlib
import json
import sqlite3
import time
from pathlib import Path

from pipeline.search.discovery import candidate_id, discovery_region


class PriorArtRegistry:
    def __init__(self, database: str):
        Path(database).parent.mkdir(parents=True, exist_ok=True)
        self.database = database
        with self._connect() as conn:
            conn.execute("""CREATE TABLE IF NOT EXISTS prior_art_records (
                candidate_id TEXT NOT NULL, source_id TEXT NOT NULL,
                region TEXT NOT NULL,
                genome TEXT NOT NULL, source_type TEXT NOT NULL,
                citation TEXT, evidence_level TEXT NOT NULL,
                publication_year INTEGER, record_hash TEXT NOT NULL,
                imported_at REAL NOT NULL,
                PRIMARY KEY(candidate_id, source_id)
            )""")
            conn.execute("CREATE INDEX IF NOT EXISTS prior_region_idx "
                         "ON prior_art_records(region)")
            # One-way compatibility migration from the original single-source
            # table. Multiple citations can now coexist for one candidate.
            legacy = conn.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name='prior_art'"
            ).fetchone()
            if legacy:
                for row in conn.execute(
                        "SELECT candidate_id, region, genome, source_type, source_id, "
                        "citation, evidence_level, imported_at FROM prior_art"):
                    record_hash = hashlib.sha256(
                        json.dumps(list(row[:-1]), sort_keys=True).encode()).hexdigest()
                    conn.execute(
                        "INSERT OR IGNORE INTO prior_art_records VALUES "
                        "(?,?,?,?,?,?,?,?,?,?)",
                        (row[0], row[4], row[1], row[2], row[3], row[5],
                         row[6], None, record_hash, row[7]))

    def _connect(self):
        return sqlite3.connect(self.database, timeout=60)

    def add(self, genome: tuple, source_type: str, source_id: str,
            citation: str = '', evidence_level: str = 'reported',
            publication_year: int | None = None):
        if not str(source_id).strip():
            raise ValueError('prior-art source_id is required')
        if publication_year is not None and not 1800 <= int(publication_year) <= 2200:
            raise ValueError('invalid publication_year')
        cid = candidate_id(genome)
        region = '|'.join(discovery_region(genome))
        payload = {
            'candidate_id': cid, 'source_id': str(source_id),
            'region': region, 'genome': repr(genome),
            'source_type': str(source_type), 'citation': str(citation),
            'evidence_level': str(evidence_level),
            'publication_year': publication_year,
        }
        record_hash = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()
        ).hexdigest()
        with self._connect() as conn:
            conn.execute(
                "INSERT OR REPLACE INTO prior_art_records VALUES "
                "(?,?,?,?,?,?,?,?,?,?)",
                (cid, str(source_id), region, repr(genome), str(source_type),
                 str(citation), str(evidence_level), publication_year,
                 record_hash, time.time()))

    def import_csv(self, path: str) -> int:
        count = 0
        with open(path, newline='') as handle:
            for row in csv.DictReader(handle):
                genome = ast.literal_eval(row['genome'])
                self.add(genome, row.get('source_type', 'literature'),
                         row.get('source_id', ''), row.get('citation', ''),
                         row.get('evidence_level', 'reported'),
                         int(row['publication_year'])
                         if row.get('publication_year') else None)
                count += 1
        return count

    def classify(self, genome: tuple) -> dict:
        cid = candidate_id(genome)
        region = '|'.join(discovery_region(genome))
        with self._connect() as conn:
            exact = conn.execute(
                "SELECT source_type, source_id, citation, evidence_level, "
                "publication_year, record_hash FROM prior_art_records "
                "WHERE candidate_id=? ORDER BY source_id", (cid,)).fetchall()
            related = conn.execute(
                "SELECT COUNT(DISTINCT candidate_id) FROM prior_art_records "
                "WHERE region=?", (region,)).fetchone()[0]
        return {'exact_prior_art': bool(exact), 'exact_records': exact,
                'region_prior_art_count': int(related),
                'novelty_status': 'known' if exact else ('region_known' if related else 'unseen')}

    def count(self) -> int:
        with self._connect() as conn:
            return int(conn.execute(
                "SELECT COUNT(*) FROM prior_art_records").fetchone()[0])


def annotate_prior_art(frame, database: str):
    if frame is None or 'genome' not in frame.columns:
        return frame
    registry = PriorArtRegistry(database)
    statuses, exact, related = [], [], []
    for raw in frame['genome']:
        try:
            genome = ast.literal_eval(raw) if isinstance(raw, str) else tuple(raw)
            result = registry.classify(genome)
            statuses.append(result['novelty_status'])
            exact.append(result['exact_prior_art'])
            related.append(result['region_prior_art_count'])
        except (ValueError, SyntaxError, TypeError):
            statuses.append('unknown'); exact.append(False); related.append(0)
    frame = frame.copy()
    frame['prior_art_status'] = statuses
    frame['exact_prior_art'] = exact
    frame['region_prior_art_count'] = related
    return frame
