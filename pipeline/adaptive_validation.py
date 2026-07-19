"""Adaptive, coverage-safe allocation of expensive validation calculations."""

from __future__ import annotations

import json
import ast
import sqlite3
import time
from collections import defaultdict
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np

from pipeline.discovery import candidate_id, discovery_region, _quality_score


def _connect(database: str):
    Path(database).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(database, timeout=60)
    conn.execute("""CREATE TABLE IF NOT EXISTS validation_observations (
        application TEXT NOT NULL, candidate_id TEXT NOT NULL,
        material_class TEXT NOT NULL, region TEXT NOT NULL, fidelity TEXT NOT NULL,
        predicted REAL NOT NULL, observed REAL NOT NULL, absolute_error REAL NOT NULL,
        productive INTEGER NOT NULL, provenance TEXT NOT NULL, recorded_at REAL NOT NULL,
        PRIMARY KEY(application, candidate_id, fidelity)
    )""")
    conn.execute("CREATE INDEX IF NOT EXISTS validation_region_idx ON "
                 "validation_observations(application, region)")
    return conn


def record_validation(database: str, application: str, genome: tuple,
                      predicted: float, observed: float, fidelity: str,
                      productive: bool, provenance: Mapping) -> None:
    """Persist a paired prediction/observation; provenance is mandatory."""
    if not provenance or not provenance.get('source_id'):
        raise ValueError('validation provenance requires source_id')
    values = (predicted, observed)
    if not all(np.isfinite(float(x)) for x in values):
        raise ValueError('validation values must be finite')
    region = '|'.join(discovery_region(genome))
    with _connect(database) as conn:
        conn.execute("INSERT OR REPLACE INTO validation_observations VALUES "
                     "(?,?,?,?,?,?,?,?,?,?,?)", (
            application, candidate_id(genome), genome[0], region, fidelity,
            float(predicted), float(observed), abs(float(predicted) - float(observed)),
            int(bool(productive)), json.dumps(dict(provenance), sort_keys=True), time.time()))


def regional_calibration(database: str, application: str) -> dict:
    """Return MAE, bias, disagreement, and productivity for each chemistry region."""
    with _connect(database) as conn:
        rows = conn.execute("""SELECT region, COUNT(*), AVG(absolute_error),
            AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
            WHERE application=? GROUP BY region""", (application,)).fetchall()
    return {region: {'n': int(n), 'mae': float(mae), 'bias': float(bias),
                     'productivity': float(productivity), 'best_observed': float(best)}
            for region, n, mae, bias, productivity, best in rows}


def class_calibration(database: str, application: str) -> dict:
    """Aggregate calibration when an exact chemistry region is still unseen."""
    with _connect(database) as conn:
        rows = conn.execute("""SELECT material_class, COUNT(*), AVG(absolute_error),
            AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
            WHERE application=? GROUP BY material_class""", (application,)).fetchall()
    return {material_class: {'n': int(n), 'mae': float(mae), 'bias': float(bias),
                             'productivity': float(productivity), 'best_observed': float(best)}
            for material_class, n, mae, bias, productivity, best in rows}


def _calibration_for(genome: tuple, regional: dict, by_class: dict) -> dict:
    """Use the most local available calibration without inventing evidence."""
    return regional.get('|'.join(discovery_region(genome)), by_class.get(genome[0], {}))


def priority_adjustment(database: str, application: str, genomes: Sequence[tuple]) -> float:
    """Lower disagreement priority; defer repeatedly unproductive regions."""
    stats = regional_calibration(database, application)
    classes = class_calibration(database, application)
    matched = [_calibration_for(g, stats, classes) for g in genomes]
    matched = [x for x in matched if x]
    if not matched:
        return 0.0
    disagreement_bonus = np.mean([min(x['mae'], 2.0) for x in matched])
    productivity = np.mean([x['productivity'] for x in matched])
    evidence = np.mean([min(x['n'] / 3.0, 1.0) for x in matched])
    # Branch scheduler processes lower numbers first. This never prunes.
    return float(-0.5 * disagreement_bonus + 0.5 * evidence * (1.0 - productivity))


def allocate_validation_batch(candidates: Sequence[tuple], objectives: np.ndarray,
                              n_select: int, database: str, application: str,
                              min_per_class: int = 1,
                              uncertainties=None) -> list[int]:
    """Guarantee class quotas, then allocate by improvement/error/uncertainty."""
    if n_select <= 0 or not candidates:
        return []
    n_select = min(n_select, len(candidates))
    objectives = np.asarray(objectives, float)
    quality = _quality_score(objectives)
    uncertainty = np.zeros(len(candidates)) if uncertainties is None else np.asarray(uncertainties, float)
    if len(uncertainty) != len(candidates):
        raise ValueError('uncertainties must match candidates')
    if np.ptp(uncertainty) > 0:
        uncertainty = (uncertainty - uncertainty.min()) / np.ptp(uncertainty)
    stats = regional_calibration(database, application)
    classes = class_calibration(database, application)
    candidate_stats = [_calibration_for(g, stats, classes) for g in candidates]
    error = np.array([item.get('mae', 0.0) for item in candidate_stats])
    if np.ptp(error) > 0:
        error = (error - error.min()) / np.ptp(error)
    productivity = np.array([item.get('productivity', 0.5) for item in candidate_stats])
    primary = objectives[:, 0]
    improvement = np.array([
        max(candidate_stats[i].get('best_observed', primary[i]) - primary[i], 0.0)
        for i in range(len(candidates))])
    if np.ptp(improvement) > 0:
        improvement = (improvement - improvement.min()) / np.ptp(improvement)
    else:
        improvement = quality
    score = (0.30 * quality + 0.20 * improvement + 0.25 * uncertainty +
             0.15 * error + 0.10 * productivity)
    ids = [candidate_id(g) for g in candidates]

    by_class = defaultdict(list)
    for i, genome in enumerate(candidates):
        by_class[genome[0]].append(i)
    quota_required = min_per_class * len(by_class)
    if quota_required > n_select:
        raise ValueError(f'validation budget {n_select} cannot satisfy class quota {quota_required}')
    selected = []
    for material_class in sorted(by_class):
        ranked = sorted(by_class[material_class], key=lambda i: (-score[i], ids[i]))
        selected.extend(ranked[:min(min_per_class, len(ranked))])
    selected = sorted(set(selected), key=lambda i: (-score[i], ids[i]))
    chosen = set(selected)
    remainder = sorted((i for i in range(len(candidates)) if i not in chosen),
                       key=lambda i: (-score[i], ids[i]))
    selected.extend(remainder[:n_select - len(selected)])
    return selected


def experimental_slate(candidates: Sequence[tuple], objectives: np.ndarray,
                       n_select: int) -> list[int]:
    """Preserve chemistry diversity: one regional champion before repeats."""
    quality = _quality_score(np.asarray(objectives, float))
    ids = [candidate_id(g) for g in candidates]
    by_region = defaultdict(list)
    for i, genome in enumerate(candidates):
        by_region[discovery_region(genome)].append(i)
    champions = [max(indices, key=lambda i: (quality[i], ids[i]))
                 for _, indices in sorted(by_region.items())]
    champions.sort(key=lambda i: (-quality[i], ids[i]))
    selected = champions[:min(n_select, len(champions))]
    chosen = set(selected)
    remainder = sorted((i for i in range(len(candidates)) if i not in chosen),
                       key=lambda i: (-quality[i], ids[i]))
    selected.extend(remainder[:max(0, n_select - len(selected))])
    return selected


def record_screening_frame(database: str, application: str,
                           predictions: Mapping[str, float], frame,
                           observed_column: str, fidelity: str,
                           productive_threshold: float) -> int:
    """Ingest paired surrogate/high-fidelity rows from a screening DataFrame."""
    recorded = 0
    if frame is None or 'genome' not in frame.columns or observed_column not in frame.columns:
        return recorded
    for _, row in frame.iterrows():
        try:
            genome = ast.literal_eval(row['genome']) if isinstance(row['genome'], str) else tuple(row['genome'])
            cid = candidate_id(genome)
            observed = float(row[observed_column])
            if cid not in predictions or not np.isfinite(observed) or not bool(row.get('valid', True)):
                continue
            record_validation(database, application, genome, predictions[cid], observed,
                              fidelity, observed <= productive_threshold,
                              {'source_id': f'{fidelity}:{cid}'})
            recorded += 1
        except (ValueError, TypeError, SyntaxError):
            continue
    return recorded


def persist_experimental_slate(database: str, application: str,
                               candidates: Sequence[tuple], objectives: np.ndarray,
                               indices: Sequence[int]) -> None:
    """Persist the diverse shortlist for synthesis/experimental handoff."""
    with _connect(database) as conn:
        conn.execute("""CREATE TABLE IF NOT EXISTS experimental_slate (
            application TEXT NOT NULL, rank INTEGER NOT NULL, candidate_id TEXT NOT NULL,
            material_class TEXT NOT NULL, region TEXT NOT NULL, genome TEXT NOT NULL,
            objectives TEXT NOT NULL, created_at REAL NOT NULL,
            PRIMARY KEY(application, rank))""")
        conn.execute('DELETE FROM experimental_slate WHERE application=?', (application,))
        for rank, i in enumerate(indices, 1):
            genome = candidates[i]
            conn.execute('INSERT INTO experimental_slate VALUES (?,?,?,?,?,?,?,?)', (
                application, rank, candidate_id(genome), genome[0],
                '|'.join(discovery_region(genome)), repr(genome),
                json.dumps([float(x) for x in objectives[i]]), time.time()))
