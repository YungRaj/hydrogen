"""Versioned evidence registry for literature, patents, and experiments."""

import ast
import csv
import sqlite3
import time
from pathlib import Path

from pipeline.search.discovery import candidate_id, discovery_region


class PriorArtRegistry:
    def __init__(self, database: str):
        Path(database).parent.mkdir(parents=True, exist_ok=True)
        self.database = database
        with self._connect() as conn:
            conn.execute("""CREATE TABLE IF NOT EXISTS prior_art (
                candidate_id TEXT PRIMARY KEY, region TEXT NOT NULL,
                genome TEXT NOT NULL, source_type TEXT NOT NULL,
                source_id TEXT NOT NULL, citation TEXT,
                evidence_level TEXT NOT NULL, imported_at REAL NOT NULL
            )""")
            conn.execute("CREATE INDEX IF NOT EXISTS prior_region_idx ON prior_art(region)")

    def _connect(self):
        return sqlite3.connect(self.database, timeout=60)

    def add(self, genome: tuple, source_type: str, source_id: str,
            citation: str = '', evidence_level: str = 'reported'):
        cid = candidate_id(genome)
        region = '|'.join(discovery_region(genome))
        with self._connect() as conn:
            conn.execute("INSERT OR REPLACE INTO prior_art VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
                         (cid, region, repr(genome), source_type, source_id,
                          citation, evidence_level, time.time()))

    def import_csv(self, path: str) -> int:
        count = 0
        with open(path, newline='') as handle:
            for row in csv.DictReader(handle):
                genome = ast.literal_eval(row['genome'])
                self.add(genome, row.get('source_type', 'literature'),
                         row.get('source_id', ''), row.get('citation', ''),
                         row.get('evidence_level', 'reported'))
                count += 1
        return count

    def classify(self, genome: tuple) -> dict:
        cid = candidate_id(genome)
        region = '|'.join(discovery_region(genome))
        with self._connect() as conn:
            exact = conn.execute("SELECT source_type, source_id, citation, evidence_level "
                                 "FROM prior_art WHERE candidate_id=?", (cid,)).fetchall()
            related = conn.execute("SELECT COUNT(*) FROM prior_art WHERE region=?", (region,)).fetchone()[0]
        return {'exact_prior_art': bool(exact), 'exact_records': exact,
                'region_prior_art_count': int(related),
                'novelty_status': 'known' if exact else ('region_known' if related else 'unseen')}

    def count(self) -> int:
        with self._connect() as conn:
            return int(conn.execute("SELECT COUNT(*) FROM prior_art").fetchone()[0])


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
