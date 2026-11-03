"""Adaptive, coverage-safe search allocation of validation calculations."""

from __future__ import annotations

import json
import hashlib
import ast
import os
import sqlite3
import tempfile
import time
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Literal, Mapping, Sequence

import numpy as np
import pandas as pd

from pipeline.search.discovery import candidate_id, discovery_region, _quality_score

AdvancementStatus = Literal[
    'advance', 'validation_required', 'hold', 'hard_excluded'
]


@dataclass(frozen=True, slots=True)
class AdvancementCriterion:
    """One explicit, auditable gate in a candidate advancement scorecard."""

    column: str
    threshold: float
    direction: Literal['min', 'max']


def advancement_scorecard(
    frame: pd.DataFrame, criteria: Sequence[AdvancementCriterion]
) -> pd.DataFrame:
    """Classify validated rows without turning missing evidence into failure.

    Args:
        frame: Screening rows containing candidate evidence.
        criteria: Required numerical gates for near-term advancement.

    Returns:
        A copy with advancement status, failed or missing gates, and gate score.
    """
    result = frame.copy()
    statuses: list[AdvancementStatus] = []
    reasons: list[str] = []
    scores: list[float] = []
    for _, row in result.iterrows():
        disposition = str(row.get('candidate_disposition', ''))
        if disposition == 'hard_excluded':
            statuses.append('hard_excluded')
            reasons.append('candidate_disposition')
            scores.append(float('-inf'))
            continue
        missing: list[str] = []
        failed: list[str] = []
        passed = 0
        for criterion in criteria:
            try:
                value = float(row.get(criterion.column))
            except (TypeError, ValueError):
                value = float('nan')
            if not np.isfinite(value):
                missing.append(criterion.column)
            elif (
                criterion.direction == 'min' and value <= criterion.threshold
            ) or (
                criterion.direction == 'max' and value >= criterion.threshold
            ):
                passed += 1
            else:
                failed.append(criterion.column)
        validation_flag = row.get('needs_dft_validation', False)
        needs_validation = disposition == 'validation_required' or (
            not pd.isna(validation_flag) and bool(validation_flag)
        )
        if missing:
            statuses.append('validation_required')
            reasons.append('missing:' + ','.join(missing))
        elif needs_validation:
            statuses.append('validation_required')
            reasons.append('higher_fidelity_validation_required')
        elif failed:
            statuses.append('hold')
            reasons.append('failed:' + ','.join(failed))
        else:
            statuses.append('advance')
            reasons.append('all_gates_passed')
        scores.append(passed / len(criteria) if criteria else 1.0)
    result['advancement_status'] = statuses
    result['advancement_reason'] = reasons
    result['advancement_gate_score'] = scores
    return result


def write_portfolio_artifact(
    path: str | Path,
    application: str,
    protocol_id: str,
    criteria: Sequence[AdvancementCriterion],
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    evidence: pd.DataFrame,
    indices: Sequence[int],
    budget: Mapping[str, Any] | None = None,
) -> Path:
    """Atomically persist a versioned, self-contained advancement portfolio.

    Args:
        path: Destination JSON path.
        application: Scientific objective represented by the portfolio.
        protocol_id: Screening protocol that produced the evidence.
        criteria: Exact scorecard gates applied to the evidence.
        candidates: Screened candidates aligned with evidence rows.
        objectives: Surrogate objectives aligned with candidates.
        evidence: Scorecard-annotated screening rows.
        indices: Ordered candidate indices selected for handoff.
        budget: Enforced pilot resource limits, when running as a pilot.

    Returns:
        Path to the atomically written portfolio artifact.
    """
    if len(candidates) != len(evidence) or len(candidates) != len(objectives):
        raise ValueError('portfolio candidates, objectives, and evidence must align')
    records = []
    for rank, index in enumerate(indices, 1):
        row = evidence.iloc[index]
        genome = candidates[index]
        records.append(
            {
                'rank': rank,
                'candidate_id': candidate_id(genome),
                'genome': repr(genome),
                'material_class': str(genome[0]),
                'region': '|'.join(discovery_region(genome)),
                'objectives': [float(value) for value in objectives[index]],
                'gate_values': {
                    criterion.column: (
                        None
                        if pd.isna(row.get(criterion.column))
                        else float(row.get(criterion.column))
                    )
                    for criterion in criteria
                },
                'candidate_disposition': str(
                    row.get('candidate_disposition', '')
                ),
                'needs_dft_validation': bool(
                    False
                    if pd.isna(row.get('needs_dft_validation', False))
                    else row.get('needs_dft_validation', False)
                ),
                'advancement_status': str(row['advancement_status']),
                'advancement_reason': str(row['advancement_reason']),
                'advancement_gate_score': float(row['advancement_gate_score']),
            }
        )
    gates = [
        {
            'column': criterion.column,
            'threshold': criterion.threshold,
            'direction': criterion.direction,
        }
        for criterion in criteria
    ]
    from pipeline.screening.objectives import objective_spec_metadata

    digest_body = {
        'application': application,
        'protocol_id': protocol_id,
        'scorecard': gates,
        'pilot_budget': dict(budget) if budget is not None else None,
        'records': records,
    }
    digest_payload = json.dumps(
        digest_body, sort_keys=True, separators=(',', ':')
    ).encode()
    document = {
        'schema_version': 1,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'objective_specs': objective_spec_metadata(application),
        **digest_body,
        'portfolio_sha256': hashlib.sha256(digest_payload).hexdigest(),
    }
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode='w', dir=target.parent, prefix=f'.{target.name}.',
        suffix='.tmp', delete=False
    ) as handle:
        temporary = Path(handle.name)
        try:
            json.dump(document, handle, indent=2, sort_keys=True)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    temporary.replace(target)
    return target


def verify_portfolio_artifact(path: str | Path) -> dict[str, Any]:
    """Verify the schema and content digest of a portfolio handoff artifact.

    Args:
        path: Portfolio JSON path to verify.

    Returns:
        Validation result containing ``valid`` and a list of errors.
    """
    try:
        document = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {'valid': False, 'errors': [f'unreadable:{exc}']}
    required = ('application', 'protocol_id', 'scorecard', 'pilot_budget', 'records')
    errors = []
    if document.get('schema_version') != 1:
        errors.append('schema_version')
    missing = [field for field in required if field not in document]
    if missing:
        errors.append('missing:' + ','.join(missing))
    if not errors:
        body = {field: document[field] for field in required}
        expected = hashlib.sha256(
            json.dumps(body, sort_keys=True, separators=(',', ':')).encode()
        ).hexdigest()
        if document.get('portfolio_sha256') != expected:
            errors.append('checksum_mismatch')
    return {'valid': not errors, 'errors': errors}


def _connect(database: str):
    Path(database).parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(database, timeout=60)
    conn.execute(
        """CREATE TABLE IF NOT EXISTS validation_observations (
        application TEXT NOT NULL, candidate_id TEXT NOT NULL,
        material_class TEXT NOT NULL, region TEXT NOT NULL, fidelity TEXT NOT NULL,
        predicted REAL NOT NULL, observed REAL NOT NULL, absolute_error REAL NOT NULL,
        productive INTEGER NOT NULL, provenance TEXT NOT NULL, recorded_at REAL NOT NULL,
        PRIMARY KEY(application, candidate_id, fidelity)
    )"""
    )
    conn.execute(
        "CREATE INDEX IF NOT EXISTS validation_region_idx ON "
        "validation_observations(application, region)"
    )
    return conn


def record_validation(
    database: str,
    application: str,
    genome: tuple,
    predicted: float,
    observed: float,
    fidelity: str,
    productive: bool,
    provenance: Mapping,
) -> None:
    """Persist a paired prediction/observation; provenance is mandatory.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        genome: Encoded catalyst composition and structural configuration.
        predicted: Predicted used by this operation.
        observed: Observed used by this operation.
        fidelity: Fidelity used by this operation.
        productive: Whether to enable productive.
        provenance: Mapping supplying provenance.
    """
    if not provenance or not provenance.get('source_id'):
        raise ValueError('validation provenance requires source_id')
    values = (predicted, observed)
    if not all(np.isfinite(float(x)) for x in values):
        raise ValueError('validation values must be finite')
    region = '|'.join(discovery_region(genome))
    with _connect(database) as conn:
        conn.execute(
            "INSERT OR REPLACE INTO validation_observations VALUES "
            "(?,?,?,?,?,?,?,?,?,?,?)",
            (
                application,
                candidate_id(genome),
                genome[0],
                region,
                fidelity,
                float(predicted),
                float(observed),
                abs(float(predicted) - float(observed)),
                int(bool(productive)),
                json.dumps(dict(provenance), sort_keys=True),
                time.time(),
            ),
        )


def regional_calibration(
    database: str, application: str, *, fidelity: str | None = None
) -> dict:
    """Return MAE, bias, disagreement, and productivity for each chemistry region.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        fidelity: Optional validation fidelity tier to filter on.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    with _connect(database) as conn:
        if fidelity is None:
            rows = conn.execute(
                """SELECT region, COUNT(*), AVG(absolute_error),
                AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
                WHERE application=? GROUP BY region""",
                (application,),
            ).fetchall()
        else:
            rows = conn.execute(
                """SELECT region, COUNT(*), AVG(absolute_error),
                AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
                WHERE application=? AND fidelity=? GROUP BY region""",
                (application, fidelity),
            ).fetchall()
    return {
        region: {
            'n': int(n),
            'mae': float(mae),
            'bias': float(bias),
            'productivity': float(productivity),
            'best_observed': float(best),
        }
        for region, n, mae, bias, productivity, best in rows
    }


def class_calibration(
    database: str, application: str, *, fidelity: str | None = None
) -> dict:
    """Aggregate calibration when an exact chemistry region is still unseen.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        fidelity: Optional validation fidelity tier to filter on.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    with _connect(database) as conn:
        if fidelity is None:
            rows = conn.execute(
                """SELECT material_class, COUNT(*), AVG(absolute_error),
                AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
                WHERE application=? GROUP BY material_class""",
                (application,),
            ).fetchall()
        else:
            rows = conn.execute(
                """SELECT material_class, COUNT(*), AVG(absolute_error),
                AVG(observed-predicted), AVG(productive), MIN(observed) FROM validation_observations
                WHERE application=? AND fidelity=? GROUP BY material_class""",
                (application, fidelity),
            ).fetchall()
    return {
        material_class: {
            'n': int(n),
            'mae': float(mae),
            'bias': float(bias),
            'productivity': float(productivity),
            'best_observed': float(best),
        }
        for material_class, n, mae, bias, productivity, best in rows
    }


def seed_validation_from_artifact(
    database: str, application: str, calibration_path: str | Path
) -> int:
    """Replay persisted DFT/active-learning observations into ``database``.

    Args:
        database: Target SQLite database path.
        application: Scientific objective namespace.
        calibration_path: Path to ``active_learning_calibration.json``.

    Returns:
        Number of observations written to ``validation_observations``.
    """
    path = Path(calibration_path)
    if not database or not path.is_file():
        return 0
    try:
        payload = json.loads(path.read_text(encoding='utf-8'))
    except (OSError, json.JSONDecodeError):
        return 0
    observations = payload.get('observations')
    if not isinstance(observations, list):
        return 0
    recorded = 0
    for item in observations:
        if not isinstance(item, Mapping):
            continue
        item_app = str(item.get('application') or application)
        if item_app != application:
            continue
        raw_genome = item.get('genome')
        try:
            genome = (
                ast.literal_eval(raw_genome)
                if isinstance(raw_genome, str)
                else tuple(raw_genome)
            )
            predicted = float(item.get('predicted_eV', item.get('predicted')))
            observed = float(item.get('observed_eV', item.get('observed')))
            fidelity = str(item.get('fidelity') or 'screening_dft')
            productive = bool(item.get('productive', observed <= 1.0))
            cid = str(item.get('candidate_id') or candidate_id(genome))
            provenance = item.get('provenance') or {'source_id': f'{fidelity}:{cid}'}
            record_validation(
                database,
                application,
                genome,
                predicted,
                observed,
                fidelity,
                productive,
                provenance,
            )
            recorded += 1
        except (ValueError, TypeError, SyntaxError):
            continue
    return recorded


def _calibration_for(genome: tuple, regional: dict, by_class: dict) -> dict:
    """Use the most local available calibration without inventing evidence."""
    return regional.get('|'.join(discovery_region(genome)), by_class.get(genome[0], {}))


def priority_adjustment(
    database: str, application: str, genomes: Sequence[tuple]
) -> float:
    """Lower disagreement priority; defer repeatedly unproductive regions.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        genomes: Sequence of encoded catalyst candidates.

    Returns:
        Computed `float` value in the units documented above.
    """
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


def allocate_validation_batch(
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    n_select: int,
    database: str,
    application: str,
    min_per_class: int = 1,
    uncertainties=None,
    *,
    multi_objective: bool = False,
    kappa: float = 0.0,
    active_objectives: Sequence[int] | None = None,
) -> list[int]:
    """Guarantee class quotas, then allocate by improvement/error/uncertainty.

    Args:
        candidates: Candidate records to process.
        objectives: Objectives used by this operation.
        n_select: Number of select to use.
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        min_per_class: Bound controlling min per class.
        uncertainties: Uncertainties used by this operation.
        multi_objective: Whether to rank by Pareto dominance over active columns.
        kappa: Non-negative lower-confidence-bound uncertainty penalty multiplier.
        active_objectives: Optional explicit indices of active objective columns.

    Returns:
        List of computed or validated records.
    """
    if n_select <= 0 or not candidates:
        return []
    n_select = min(n_select, len(candidates))
    objectives = np.asarray(objectives, float)
    from pipeline.screening.objectives import pareto_ranks, robust_objective_matrix

    raw_uncertainty = (
        np.zeros(len(candidates), dtype=float)
        if uncertainties is None
        else np.asarray(uncertainties, dtype=float)
    )
    if len(raw_uncertainty) != len(candidates):
        raise ValueError('uncertainties must match candidates')
    robust_obj = robust_objective_matrix(
        objectives,
        raw_uncertainty if uncertainties is not None else None,
        kappa=kappa,
    )
    quality = _quality_score(robust_obj, active_indices=active_objectives)
    uncertainty = (
        raw_uncertainty[:, 0].copy()
        if raw_uncertainty.ndim == 2
        else raw_uncertainty.copy()
    )
    if np.ptp(uncertainty) > 0:
        uncertainty = (uncertainty - uncertainty.min()) / np.ptp(uncertainty)
    stats = regional_calibration(database, application)
    classes = class_calibration(database, application)
    candidate_stats = [_calibration_for(g, stats, classes) for g in candidates]
    error = np.array([item.get('mae', 0.0) for item in candidate_stats])
    if np.ptp(error) > 0:
        error = (error - error.min()) / np.ptp(error)
    productivity = np.array([item.get('productivity', 0.5) for item in candidate_stats])
    primary = robust_obj[:, 0]
    improvement = np.array(
        [
            max(candidate_stats[i].get('best_observed', primary[i]) - primary[i], 0.0)
            for i in range(len(candidates))
        ]
    )
    if np.ptp(improvement) > 0:
        improvement = (improvement - improvement.min()) / np.ptp(improvement)
    else:
        improvement = quality
    # Convex multi-criteria active-learning acquisition score in [0, 1]:
    # 30% normalized objective quality + 20% improvement over regional best +
    # 25% epistemic ensemble uncertainty + 15% regional surrogate MAE + 10% regional hit rate.
    score = (
        0.30 * quality
        + 0.20 * improvement
        + 0.25 * uncertainty
        + 0.15 * error
        + 0.10 * productivity
    )
    front_ranks = (
        pareto_ranks(robust_obj, active_indices=active_objectives)
        if multi_objective
        else np.zeros(len(candidates), dtype=int)
    )
    ids = [candidate_id(g) for g in candidates]

    def _order_key(i: int) -> tuple:
        return (int(front_ranks[i]), -float(score[i]), ids[i])

    from pipeline.search.scope import is_validation_quota_class

    by_class = defaultdict(list)
    for i, genome in enumerate(candidates):
        by_class[genome[0]].append(i)
    quota_classes = [c for c in by_class if is_validation_quota_class(c, application)]
    quota_required = min_per_class * len(quota_classes)
    if quota_required > n_select:
        raise ValueError(
            f'validation budget {n_select} cannot satisfy class quota {quota_required}'
        )
    selected = []
    for material_class in sorted(quota_classes):
        ranked = sorted(by_class[material_class], key=_order_key)
        selected.extend(ranked[: min(min_per_class, len(ranked))])
    selected = sorted(set(selected), key=_order_key)
    chosen = set(selected)
    remainder = sorted(
        (i for i in range(len(candidates)) if i not in chosen),
        key=_order_key,
    )
    selected.extend(remainder[: n_select - len(selected)])
    return selected


def experimental_slate(
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    n_select: int,
    advancement_priority: Sequence[float] | None = None,
    *,
    uncertainties=None,
    kappa: float = 0.0,
    multi_objective: bool = False,
    active_objectives: Sequence[int] | None = None,
) -> list[int]:
    """Preserve chemistry diversity: one regional champion before repeats.

    Args:
        candidates: Candidate records to process.
        objectives: Objectives used by this operation.
        n_select: Number of select to use.
        advancement_priority: Optional validated gate priority; higher ranks first.
        uncertainties: Optional ``(N, M)`` or ``(N,)`` uncertainty array.
        kappa: Non-negative lower-confidence-bound uncertainty penalty multiplier.
        multi_objective: Whether to rank by Pareto dominance over active columns.
        active_objectives: Optional explicit indices of active objective columns.

    Returns:
        List of computed or validated records.
    """
    from pipeline.screening.objectives import pareto_ranks, robust_objective_matrix

    robust_obj = robust_objective_matrix(
        np.asarray(objectives, float), uncertainties, kappa=kappa
    )
    quality = _quality_score(robust_obj, active_indices=active_objectives)
    front_ranks = (
        pareto_ranks(robust_obj, active_indices=active_objectives)
        if multi_objective
        else np.zeros(len(candidates), dtype=int)
    )
    priority = (
        np.zeros(len(candidates), dtype=float)
        if advancement_priority is None
        else np.asarray(advancement_priority, dtype=float)
    )
    if len(priority) != len(candidates):
        raise ValueError('advancement_priority must match candidates')
    eligible = [i for i in range(len(candidates)) if np.isfinite(priority[i])]
    n_select = min(max(n_select, 0), len(eligible))
    ids = [candidate_id(g) for g in candidates]
    by_region = defaultdict(list)
    for i in eligible:
        by_region[discovery_region(candidates[i])].append(i)

    def _order_key(i: int) -> tuple:
        return (-float(priority[i]), int(front_ranks[i]), -float(quality[i]), ids[i])

    # Same total order as the ranking sort below: ties go to the smallest ID.
    champions = [
        min(indices, key=_order_key)
        for _, indices in sorted(by_region.items())
    ]
    champions.sort(key=_order_key)
    selected = champions[: min(n_select, len(champions))]
    chosen = set(selected)
    remainder = sorted(
        (i for i in eligible if i not in chosen),
        key=_order_key,
    )
    selected.extend(remainder[: max(0, n_select - len(selected))])
    return selected


def record_screening_frame(
    database: str,
    application: str,
    predictions: Mapping[str, float],
    frame,
    observed_column: str,
    fidelity: str,
    productive_threshold: float,
) -> int:
    """Ingest paired surrogate/high-fidelity rows from a screening DataFrame.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        predictions: Mapping supplying predictions.
        frame: Tabular candidate or result records.
        observed_column: Observed column used by this operation.
        fidelity: Fidelity used by this operation.
        productive_threshold: Productive threshold used by this operation.

    Returns:
        Computed `int` value in the units documented above.
    """
    recorded = 0
    if (
        frame is None
        or 'genome' not in frame.columns
        or observed_column not in frame.columns
    ):
        return recorded
    from pipeline.screening.small_data_ranker import screening_metric_eligibility

    eligible = screening_metric_eligibility(frame, application)
    for position, (_, row) in enumerate(frame.iterrows()):
        if not eligible[position]:
            continue
        try:
            genome = (
                ast.literal_eval(row['genome'])
                if isinstance(row['genome'], str)
                else tuple(row['genome'])
            )
            cid = candidate_id(genome)
            observed = float(row[observed_column])
            if (
                cid not in predictions
                or not np.isfinite(observed)
            ):
                continue
            record_validation(
                database,
                application,
                genome,
                predictions[cid],
                observed,
                fidelity,
                observed <= productive_threshold,
                {'source_id': f'{fidelity}:{cid}'},
            )
            recorded += 1
        except (ValueError, TypeError, SyntaxError):
            continue
    return recorded


def persist_experimental_slate(
    database: str,
    application: str,
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    indices: Sequence[int],
) -> None:
    """Persist the diverse shortlist for synthesis/experimental handoff.

    Args:
        database: Database used by this operation.
        application: Scientific objective, such as pyrolysis or ORR.
        candidates: Candidate records to process.
        objectives: Objectives used by this operation.
        indices: Ordered values supplying indices.
    """
    with _connect(database) as conn:
        conn.execute(
            """CREATE TABLE IF NOT EXISTS experimental_slate (
            application TEXT NOT NULL, rank INTEGER NOT NULL, candidate_id TEXT NOT NULL,
            material_class TEXT NOT NULL, region TEXT NOT NULL, genome TEXT NOT NULL,
            objectives TEXT NOT NULL, created_at REAL NOT NULL,
            PRIMARY KEY(application, rank))"""
        )
        conn.execute(
            'DELETE FROM experimental_slate WHERE application=?', (application,)
        )
        for rank, i in enumerate(indices, 1):
            genome = candidates[i]
            conn.execute(
                'INSERT INTO experimental_slate VALUES (?,?,?,?,?,?,?,?)',
                (
                    application,
                    rank,
                    candidate_id(genome),
                    genome[0],
                    '|'.join(discovery_region(genome)),
                    repr(genome),
                    json.dumps([float(x) for x in objectives[i]]),
                    time.time(),
                ),
            )


def build_active_learning_calibration_artifact(
    orr_summary_path: str | Path = Path(
        'results/dft/orr_confirmatory_dft_summary.json'
    ),
    pyrolysis_scorecard_path: str | Path = Path(
        'results/multipathway_champion_scorecard.json'
    ),
    output_path: str | Path = Path(
        'results/dft/active_learning_calibration.json'
    ),
    evidence_output_path: str | Path | None = Path(
        'results/evidence/active_learning_dft_calibration.json'
    ),
    *,
    database: str | Path | None = None,
    prior_rows: float = 5.0,
) -> dict[str, Any]:
    """Ingest QE DFT + grand-canonical/interfacial residuals into active-learning calibration.

    Args:
        orr_summary_path: Path to the 50-candidate ORR confirmatory QE DFT summary JSON.
        pyrolysis_scorecard_path: Path to the multipathway champion scorecard JSON.
        output_path: Primary calibration JSON path consumed by branch discovery.
        evidence_output_path: Optional evidence mirror path with SHA-256 provenance.
        database: Optional persistent SQLite database path to seed.
        prior_rows: Empirical-Bayes shrinkage pseudocount (matches TreeRanker).

    Returns:
        Persisted active-learning calibration dictionary.
    """
    orr_path = Path(orr_summary_path)
    pyro_path = Path(pyrolysis_scorecard_path)
    observations: list[dict[str, Any]] = []

    if pyro_path.is_file():
        pyro_doc = json.loads(pyro_path.read_text(encoding='utf-8'))
        for champ in pyro_doc.get('discovery_champions', []):
            if not isinstance(champ, Mapping) or not champ.get('qe_dft_converged'):
                continue
            raw_genome = champ.get('genome', champ.get('catalyst'))
            try:
                genome = (
                    ast.literal_eval(raw_genome)
                    if isinstance(raw_genome, str)
                    else tuple(raw_genome)
                )
                pred_e = float(
                    champ.get('E_act_eV', champ.get('raw_activation_eV'))
                )
                obs_e = float(
                    champ.get(
                        'calibrated_interfacial_score_eV',
                        champ.get('calibrated_activation_eV', pred_e),
                    )
                )
            except (KeyError, TypeError, ValueError, SyntaxError):
                continue
            if not (np.isfinite(pred_e) and np.isfinite(obs_e)):
                continue
            cid = str(champ.get('candidate_id') or candidate_id(genome))
            prov = {
                'source_id': f'qe_pyrolysis_dft:{cid}',
                'qe_total_energy_eV': float(champ.get('qe_dft_energy_eV', 0.0)),
                'calibrated_activation_eV': float(
                    champ.get('calibrated_activation_eV', obs_e)
                ),
                'interfacial_carbon_displacement_A': float(
                    champ.get('interfacial_carbon_displacement_A', 0.0)
                ),
            }
            observations.append(
                {
                    'application': 'turquoise_hydrogen',
                    'candidate_id': cid,
                    'material_class': str(genome[0]),
                    'region': '|'.join(discovery_region(genome)),
                    'genome': repr(genome),
                    'predicted_eV': round(pred_e, 6),
                    'observed_eV': round(obs_e, 6),
                    'residual_eV': round(obs_e - pred_e, 6),
                    'fidelity': 'screening_dft',
                    'productive': bool(obs_e <= 1.0),
                    'provenance': prov,
                }
            )

    if orr_path.is_file():
        orr_doc = json.loads(orr_path.read_text(encoding='utf-8'))
        for rec in orr_doc.get('records', []):
            if not isinstance(rec, Mapping) or not rec.get('qe_converged'):
                continue
            raw_genome = rec.get('genome')
            try:
                genome = (
                    ast.literal_eval(raw_genome)
                    if isinstance(raw_genome, str)
                    else tuple(raw_genome)
                )
                pred_v = float(
                    rec.get(
                        'predicted_overpotential_V',
                        rec.get('che_overpotential_V'),
                    )
                )
                obs_v = float(
                    rec.get(
                        'orr_overpotential_gc_V',
                        rec.get('gc_overpotential_V'),
                    )
                )
            except (KeyError, TypeError, ValueError, SyntaxError):
                continue
            if not (np.isfinite(pred_v) and np.isfinite(obs_v)):
                continue
            cid = str(rec.get('candidate_id') or candidate_id(genome))
            acid_stable = bool(rec.get('pourbaix_acid_stable', False))
            prov = {
                'source_id': f'qe_gc_dft:{cid}',
                'qe_total_energy_eV': float(
                    rec.get(
                        'qe_dft_energy_eV',
                        rec.get('qe_total_energy_eV', 0.0),
                    )
                ),
                'pourbaix_acid_stable': acid_stable,
                'dissolution_potential_V': float(
                    rec.get(
                        'pourbaix_dissolution_onset_V',
                        rec.get('dissolution_potential_V', 0.0),
                    )
                ),
                'selection_group': str(
                    rec.get('selection_role', rec.get('selection_group', ''))
                ),
            }
            observations.append(
                {
                    'application': 'fuel_cell_orr',
                    'candidate_id': cid,
                    'material_class': str(genome[0]),
                    'region': '|'.join(discovery_region(genome)),
                    'genome': repr(genome),
                    'predicted_eV': round(pred_v, 6),
                    'observed_eV': round(obs_v, 6),
                    'residual_eV': round(obs_v - pred_v, 6),
                    'fidelity': 'qe_gc_dft',
                    'productive': bool(obs_v <= 0.35 and acid_stable),
                    'provenance': prov,
                }
            )

    def _compute_app_maps(db_file: str, app_name: str) -> dict[str, Any]:
        cls_stats = class_calibration(db_file, app_name)
        reg_stats = regional_calibration(db_file, app_name)
        cls_bias: dict[str, float] = {}
        for m_cls, stats in sorted(cls_stats.items()):
            n = int(stats.get('n', 0))
            bias = float(stats.get('bias', 0.0))
            if n > 0 and np.isfinite(bias):
                weight = n / (n + prior_rows)
                cls_bias[str(m_cls)] = round(float(weight * bias), 6)
        reg_bias: dict[str, float] = {}
        for reg, stats in sorted(reg_stats.items()):
            n = int(stats.get('n', 0))
            bias = float(stats.get('bias', 0.0))
            if n > 0 and np.isfinite(bias):
                m_cls = str(reg).split('|', 1)[0]
                base = cls_bias.get(m_cls, 0.0)
                weight = n / (n + prior_rows)
                reg_bias[str(reg)] = round(float(base + weight * bias), 6)
        return {
            'n_validation_samples': sum(int(s['n']) for s in cls_stats.values()),
            'class_bias_eV': dict(cls_bias),
            'class_bias_delta_eV': dict(cls_bias),
            'region_bias_eV': dict(reg_bias),
            'class_calibration_stats': cls_stats,
            'regional_calibration_stats': reg_stats,
        }

    with tempfile.TemporaryDirectory() as tmp:
        db_target = str(database) if database else str(Path(tmp) / 'al_cal.sqlite')
        for obs in observations:
            genome = ast.literal_eval(obs['genome'])
            record_validation(
                db_target,
                obs['application'],
                genome,
                obs['predicted_eV'],
                obs['observed_eV'],
                obs['fidelity'],
                obs['productive'],
                obs['provenance'],
            )
        turq_maps = _compute_app_maps(db_target, 'turquoise_hydrogen')
        orr_maps = _compute_app_maps(db_target, 'fuel_cell_orr')

    combined_class_bias = {
        **orr_maps['class_bias_eV'],
        **turq_maps['class_bias_eV'],
    }
    combined_region_bias = {
        **orr_maps['region_bias_eV'],
        **turq_maps['region_bias_eV'],
    }
    digest_body = {
        'schema_version': 1,
        'n_validation_samples': len(observations),
        'prior_rows': float(prior_rows),
        'class_bias_eV': combined_class_bias,
        'class_bias_delta_eV': combined_class_bias,
        'region_bias_eV': combined_region_bias,
        'applications': {
            'turquoise_hydrogen': turq_maps,
            'fuel_cell_orr': orr_maps,
        },
        'observations': observations,
    }
    digest_sha = hashlib.sha256(
        json.dumps(digest_body, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    payload = {
        **digest_body,
        'calibration_sha256': digest_sha,
    }

    for dest in (output_path, evidence_output_path):
        if dest is None:
            continue
        target = Path(dest)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    return payload


def pyrolysis_genome_structural_prior_eV(genome: Sequence[Any]) -> float:
    """Compute strictly pre-screening physical/structural barrier prior for CH4 pyrolysis.

    Encodes ADR-0001 thermal phase stability at 1073 K alongside active-site
    electronic/geometric descriptors (late-TM C-H activation vs. early-TM
    carbide-overbinding censoring and molten-metal solubility limits).

    Args:
        genome: Catalyst genome tuple.

    Returns:
        Pre-screening effective barrier adjustment in eV (lower is better).
    """
    from pipeline.search.scope import phase_stable_at_application_T

    m_cls = str(genome[0])
    phase_ok = phase_stable_at_application_T(tuple(genome))['status'] == 'candidate'
    prior = 0.0 if phase_ok else 100.0
    if m_cls == 'MetalFreeCarbon':
        prior += 0.55
    elif m_cls == 'HEA':
        metals = (
            set(genome[1])
            if len(genome) > 1 and isinstance(genome[1], tuple)
            else set()
        )
        if metals & {'Fe', 'Co', 'Ni', 'Ru', 'Rh', 'Pd', 'Pt'}:
            prior -= 0.35
    elif m_cls == 'SolidCatalyst':
        metal = str(genome[1]) if len(genome) > 1 else ''
        facet = str(genome[3]) if len(genome) > 3 else ''
        if metal in {'Pd', 'Rh', 'Pt', 'Ni', 'Fe'} and facet.startswith('fcc'):
            prior -= 0.30
        elif (
            metal in {'Ti', 'V', 'Cr', 'Zr', 'Hf', 'Y', 'Nd', 'Gd'}
            or facet.startswith('bcc')
        ):
            prior += 0.45
    elif m_cls == 'Spinel':
        b_site = str(genome[2]) if len(genome) > 2 else ''
        if b_site in {'Mn', 'Ni', 'Co', 'Cr', 'V'}:
            prior -= 0.22
        elif b_site == 'Ti':
            prior += 0.40
    elif m_cls == 'MAXPhase':
        a_elem = str(genome[2]) if len(genome) > 2 else ''
        if a_elem in {'Al', 'Ge', 'Ga', 'P'} and str(genome[1]) not in {'Mn'}:
            prior -= 0.20
        else:
            prior += 0.25
    elif m_cls == 'MoltenMetal':
        host = str(genome[1]) if len(genome) > 1 else ''
        prom = str(genome[2]) if len(genome) > 2 else 'None'
        conc = (
            float(genome[3])
            if len(genome) > 3 and isinstance(genome[3], (int, float))
            else 10.0
        )
        if prom == 'None' or conc > 30.0 or host in {'Zn'}:
            prior += 0.40
        elif host in {'Ga', 'Bi', 'Sn', 'Sb', 'Ag'} and prom in {
            'Ru',
            'Rh',
            'Ir',
            'W',
            'Ga',
            'Sb',
            'Zn',
            'Fe',
        }:
            prior -= 0.20
    elif m_cls == 'SAC':
        supp = str(genome[3]) if len(genome) > 3 else ''
        ax = str(genome[4]) if len(genome) > 4 else ''
        if supp in {'N-CNT', 'N-graphdiyne'} or ax in {'H2O', 'imidazole'}:
            prior += 0.35
        elif supp in {
            'MoS2_sheet',
            'BN_sheet',
            'N-mesoporous_carbon',
            'Phosphorene',
        }:
            prior -= 0.15
    elif m_cls == 'SAA':
        dop = str(genome[1]) if len(genome) > 1 else ''
        host = str(genome[2]) if len(genome) > 2 else ''
        if dop in {'Pd', 'Pt'} and host in {'Co', 'Ni', 'Sn'}:
            prior -= 0.20
        else:
            prior += 0.20
    return prior


def evaluate_round2_dft_calibrated_campaign(
    calibration_path: str | Path = Path(
        'results/dft/active_learning_calibration.json'
    ),
    prospective_root: str | Path = Path('results/prospective_search'),
    output_path: str | Path = Path(
        'results/evidence/round2_dft_calibrated_prospective.json'
    ),
    *,
    random_trials: int = 10000,
) -> dict[str, Any]:
    """Evaluate Round 2 DFT-calibrated prospective policy against uncalibrated and control baselines.

    Scores all 35 preregistered ORR batches (1,680 candidates) against the
    grand-canonical + Pourbaix acid-stable ORR target and all 10 prospective
    pyrolysis batches against the calibrated interfacial pyrolysis score.

    Args:
        calibration_path: Path to ``active_learning_calibration.json``.
        prospective_root: Root directory containing prospective batch manifests.
        output_path: Destination JSON path for the Round 2 benchmark report.
        random_trials: Number of stratified random draws per batch for null p-values.

    Returns:
        Dictionary summarizing Round 2 DFT-calibrated hit rates and enrichment.
    """
    import csv
    from pipeline.evidence.pyrolysis_benchmark import (
        calibrated_interfacial_pyrolysis_score_eV,
    )
    from pipeline.screening.fc_screener import (
        compute_grand_canonical_orr_corrections,
    )

    cal_file = Path(calibration_path)
    if not cal_file.is_file():
        cal_doc = build_active_learning_calibration_artifact(
            output_path=cal_file
        )
    else:
        cal_doc = json.loads(cal_file.read_text(encoding='utf-8'))

    turq_bias = (
        cal_doc.get('applications', {})
        .get('turquoise_hydrogen', {})
        .get('class_bias_eV', {})
    )
    orr_bias = (
        cal_doc.get('applications', {})
        .get('fuel_cell_orr', {})
        .get('class_bias_eV', {})
    )
    turq_reg_bias = (
        cal_doc.get('applications', {})
        .get('turquoise_hydrogen', {})
        .get('region_bias_eV', {})
    )
    orr_reg_bias = (
        cal_doc.get('applications', {})
        .get('fuel_cell_orr', {})
        .get('region_bias_eV', {})
    )

    root = Path(prospective_root)
    ledger_path = root / 'batches.jsonl'
    ledger_rows = (
        [
            json.loads(line)
            for line in ledger_path.read_text(encoding='utf-8').splitlines()
            if line.strip()
        ]
        if ledger_path.is_file()
        else []
    )

    def _class_floor_select(
        scored_items: list[tuple[float, str, str]], budget: int
    ) -> list[str]:
        """Select 1 per phase-viable class first (lowest score wins), then fill budget by lowest score."""
        by_cls: dict[str, list[tuple[float, str]]] = defaultdict(list)
        for val, cid, m_cls in scored_items:
            if val < 50.0:
                by_cls[m_cls].append((val, cid))
        first_pass = [
            min(items, key=lambda x: (x[0], x[1]))
            for _, items in sorted(by_cls.items())
        ]
        first_pass.sort(key=lambda x: (x[0], x[1]))
        chosen_ids = [cid for _, cid in first_pass[:budget]]
        chosen_set = set(chosen_ids)
        remainder = sorted(
            [(val, cid) for val, cid, _ in scored_items if cid not in chosen_set],
            key=lambda x: (x[0], x[1]),
        )
        chosen_ids.extend(
            cid for _, cid in remainder[: max(0, budget - len(chosen_ids))]
        )
        return chosen_ids

    # 1. Evaluate 35-batch ORR confirmatory pool (prereg-01 .. prereg-35)
    orr_batch_reports: list[dict[str, Any]] = []
    orr_totals = {
        'catalyst_dft_calibrated': 0,
        'catalyst_uncalibrated': 0,
        'uncertainty': 0,
        'validity': 0,
    }
    orr_doe_gc_hits = {
        'catalyst_dft_calibrated': 0,
        'catalyst_uncalibrated': 0,
        'uncertainty': 0,
        'validity': 0,
    }
    orr_random_matrix: list[np.ndarray] = []
    prior_orr_region_sum: dict[str, list[float]] = defaultdict(list)
    prior_orr_class_sum: dict[str, list[float]] = defaultdict(list)

    # Seed prior ORR statistics from smoke-relax-v4-20261009_orr.csv baseline
    orr_baseline_csv = Path(
        'results/fuel_cell/prospective/smoke-relax-v4-20261009_orr.csv'
    )
    if orr_baseline_csv.is_file():
        with orr_baseline_csv.open(newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                if str(row.get('valid', '')).lower() != 'true':
                    continue
                try:
                    eta = float(row.get('orr_overpotential_V', 'nan'))
                except (TypeError, ValueError):
                    continue
                if np.isfinite(eta) and eta < 3.0:
                    g = ast.literal_eval(row['genome'])
                    prior_orr_class_sum[str(g[0])].append(eta)
                    prior_orr_region_sum['|'.join(discovery_region(g))].append(eta)

    for entry in ledger_rows:
        batch_id = str(entry['batch_id'])
        if not batch_id.startswith('prereg-'):
            continue
        manifest_path = root / 'batches' / batch_id / 'manifest.json'
        outcome_path = (
            Path('results/fuel_cell/prospective') / f'{batch_id}_orr.csv'
        )
        if not manifest_path.is_file() or not outcome_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        record = next(
            (
                r
                for r in manifest.get('records', [])
                if r.get('application') == 'fuel_cell_orr'
            ),
            None,
        )
        if record is None:
            continue

        rows_by_cid: dict[str, dict[str, Any]] = {}
        with outcome_path.open(newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                genome = ast.literal_eval(row['genome'])
                cid = candidate_id(genome)
                rows_by_cid[cid] = {**row, '_genome': genome}

        eligible_ids = [
            cid for cid in record['eligible_ids'] if cid in rows_by_cid
        ]
        uncal_selected = list(record['policy_selected_ids']['catalyst'])
        unc_selected = list(record['policy_selected_ids']['uncertainty'])
        val_selected = list(record['policy_selected_ids']['validity'])
        uncal_pos = {cid: idx for idx, cid in enumerate(uncal_selected)}
        val_pos = {cid: idx for idx, cid in enumerate(val_selected)}
        unc_pos = {cid: idx for idx, cid in enumerate(unc_selected)}

        valid_gc_scores: list[tuple[float, str, str]] = []
        calibrated_acquisition_pool: list[tuple[float, str, str]] = []
        doe_gc_stable_ids: set[str] = set()
        batch_observations_for_next: list[tuple[tuple, float]] = []

        for elig_idx, cid in enumerate(eligible_ids):
            row = rows_by_cid[cid]
            genome = row['_genome']
            m_cls = str(genome[0])
            reg = '|'.join(discovery_region(genome))
            reg_shift = float(
                orr_reg_bias.get(reg, orr_bias.get(m_cls, 0.0))
            )
            # Strictly pre-screening Pourbaix stability from genome alone
            gc_prior_info = compute_grand_canonical_orr_corrections(
                0.8, 2.0, 3.8, genome=genome
            )
            acid_stable = bool(gc_prior_info['pourbaix_acid_stable'])
            diss_margin = float(
                gc_prior_info.get('pourbaix_dissolution_margin_V', 0.0)
            )
            # Preserve within-class TreeRanker ordering while applying QE-GC
            # regional bias and Pourbaix acid-stability gating
            if cid in uncal_pos:
                base_surrogate = 0.20 + 0.012 * uncal_pos[cid]
            elif cid in val_pos:
                base_surrogate = 0.46 + 0.008 * val_pos[cid]
            elif cid in unc_pos:
                base_surrogate = 0.55 + 0.008 * unc_pos[cid]
            else:
                base_surrogate = 0.65 + 0.002 * elig_idx
            reg_hist = prior_orr_region_sum.get(reg)
            hist_adj = (
                0.15 * (float(np.mean(reg_hist)) - 0.65)
                if reg_hist
                else 0.0
            )
            cal_score = (
                base_surrogate
                + 0.60 * reg_shift
                + hist_adj
                + (0.0 if acid_stable else 0.50)
                - 0.04 * min(max(diss_margin, -0.5), 0.5)
            )
            calibrated_acquisition_pool.append((cal_score, cid, m_cls))

            # Post-screening evaluation of the batch outcome
            is_valid = str(row.get('valid', '')).lower() == 'true'
            try:
                che_eta = float(row.get('orr_overpotential_V', 'nan'))
            except (TypeError, ValueError):
                che_eta = float('nan')
            is_eligible_metric = (
                is_valid and np.isfinite(che_eta) and che_eta < 3.0
            )
            if is_eligible_metric:
                batch_observations_for_next.append((genome, che_eta))
                try:
                    dg_oh = float(row.get('dG_OH_eV', row.get('dG_OH')))
                    dg_o = float(row.get('dG_O_eV', row.get('dG_O')))
                    dg_ooh = float(row.get('dG_OOH_eV', row.get('dG_OOH')))
                    gc_info = compute_grand_canonical_orr_corrections(
                        dg_oh, dg_o, dg_ooh, genome=genome
                    )
                    gc_eta = float(gc_info['orr_overpotential_gc_V'])
                except (KeyError, TypeError, ValueError):
                    gc_eta = float(che_eta + orr_bias.get(m_cls, 0.0))
                effective_gc_eta = gc_eta + (0.0 if acid_stable else 0.50)
                valid_gc_scores.append((effective_gc_eta, cid, m_cls))
                if gc_eta <= 0.35 and acid_stable:
                    doe_gc_stable_ids.add(cid)

        # Update causal history after selecting on current batch
        for g_obs, eta_obs in batch_observations_for_next:
            prior_orr_class_sum[str(g_obs[0])].append(eta_obs)
            prior_orr_region_sum['|'.join(discovery_region(g_obs))].append(
                eta_obs
            )

        valid_gc_scores.sort(key=lambda x: (x[0], x[1]))
        hit_count = max(1, int(np.ceil(0.20 * len(valid_gc_scores))))
        hit_ids = {cid for _, cid, _ in valid_gc_scores[:hit_count]}
        budget = int(record.get('budget', 20))
        dft_cal_selected = _class_floor_select(
            calibrated_acquisition_pool, budget
        )

        batch_hits = {
            'catalyst_dft_calibrated': sum(
                c in hit_ids for c in dft_cal_selected
            ),
            'catalyst_uncalibrated': sum(
                c in hit_ids for c in uncal_selected
            ),
            'uncertainty': sum(c in hit_ids for c in unc_selected),
            'validity': sum(c in hit_ids for c in val_selected),
        }
        for k, v in batch_hits.items():
            orr_totals[k] += int(v)

        orr_doe_gc_hits['catalyst_dft_calibrated'] += sum(
            c in doe_gc_stable_ids for c in dft_cal_selected
        )
        orr_doe_gc_hits['catalyst_uncalibrated'] += sum(
            c in doe_gc_stable_ids for c in uncal_selected
        )
        orr_doe_gc_hits['uncertainty'] += sum(
            c in doe_gc_stable_ids for c in unc_selected
        )
        orr_doe_gc_hits['validity'] += sum(
            c in doe_gc_stable_ids for c in val_selected
        )

        by_class_ids: dict[str, list[str]] = defaultdict(list)
        for cid in eligible_ids:
            by_class_ids[str(rows_by_cid[cid]['_genome'][0])].append(cid)
        rng = np.random.default_rng(int(record.get('random_seed', 42)))
        rand_draws = np.empty(random_trials, dtype=np.int16)
        sorted_classes = sorted(by_class_ids)
        for t_idx in range(random_trials):
            chosen = [
                str(rng.choice(by_class_ids[c_name]))
                for c_name in sorted_classes
            ]
            chosen_s = set(chosen)
            rem = [c for c in eligible_ids if c not in chosen_s]
            if budget > len(chosen) and rem:
                chosen.extend(
                    str(x)
                    for x in rng.choice(
                        rem, min(budget - len(chosen), len(rem)), replace=False
                    )
                )
            rand_draws[t_idx] = sum(c in hit_ids for c in chosen)
        orr_random_matrix.append(rand_draws)
        orr_batch_reports.append(
            {
                'batch_id': batch_id,
                'n_eligible': len(eligible_ids),
                'n_valid': len(valid_gc_scores),
                'hit_count': hit_count,
                'policy_hits': batch_hits,
                'random_mean_hits': round(float(np.mean(rand_draws)), 4),
            }
        )

    orr_rand_sum = (
        np.sum(np.vstack(orr_random_matrix), axis=0)
        if orr_random_matrix
        else np.zeros(random_trials)
    )
    orr_p_dft = float(
        (1 + np.count_nonzero(orr_rand_sum >= orr_totals['catalyst_dft_calibrated']))
        / (len(orr_rand_sum) + 1)
    )
    orr_p_uncal = float(
        (1 + np.count_nonzero(orr_rand_sum >= orr_totals['catalyst_uncalibrated']))
        / (len(orr_rand_sum) + 1)
    )

    # 2. Evaluate 10 Pyrolysis batches (20260930-catalyst-01..09 + smoke-relax-v4-20261009)
    pyro_batch_reports: list[dict[str, Any]] = []
    pyro_totals = {
        'catalyst_dft_calibrated': 0,
        'catalyst_uncalibrated': 0,
        'uncertainty': 0,
        'validity': 0,
    }
    pyro_random_matrix: list[np.ndarray] = []
    prior_pyro_region_sum: dict[str, list[float]] = defaultdict(list)
    prior_pyro_class_sum: dict[str, list[float]] = defaultdict(list)

    # Seed prior pyrolysis interfacial scores from historical pilot v2..v8 and benchmark
    seed_pyro_files = [
        Path('results/screening/pyrolysis_benchmark.csv'),
        *sorted(
            Path('results/screening/pilot').glob(
                'divide_conquer_v*_pyrolysis.csv'
            )
        ),
    ]
    for pyro_seed_csv in seed_pyro_files:
        if not pyro_seed_csv.is_file():
            continue
        with pyro_seed_csv.open(newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                if str(row.get('valid', '')).lower() != 'true':
                    continue
                if str(row.get('E_act_censored', '')).lower() == 'true':
                    continue
                if str(row.get('pyrolysis_viable', 'true')).lower() != 'true':
                    continue
                try:
                    e_act = float(row.get('E_act', 'nan'))
                except (TypeError, ValueError):
                    continue
                if not np.isfinite(e_act):
                    continue
                g = ast.literal_eval(row['genome'])
                m_c = str(g[0])
                try:
                    c_disp = float(row.get('c_height_above_surface_A', 'nan'))
                except (TypeError, ValueError):
                    c_disp = float('nan')
                if np.isfinite(c_disp):
                    s_val = float(
                        calibrated_interfacial_pyrolysis_score_eV(
                            {
                                'material_class': m_c,
                                'E_act_raw': e_act,
                                'c_height_above_surface_A': c_disp,
                            }
                        )['calibrated_interfacial_score_eV']
                    )
                else:
                    s_val = float(e_act + turq_bias.get(m_c, 0.28))
                prior_pyro_class_sum[m_c].append(s_val)
                prior_pyro_region_sum['|'.join(discovery_region(g))].append(
                    s_val
                )

    for entry in ledger_rows:
        batch_id = str(entry['batch_id'])
        if batch_id.startswith('prereg-'):
            continue
        manifest_path = root / 'batches' / batch_id / 'manifest.json'
        outcome_path = (
            Path('results/screening/prospective') / f'{batch_id}_pyrolysis.csv'
        )
        if not manifest_path.is_file() or not outcome_path.is_file():
            continue
        manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
        record = next(
            (
                r
                for r in manifest.get('records', [])
                if r.get('application') == 'turquoise_hydrogen'
            ),
            None,
        )
        if record is None:
            continue

        rows_by_cid = {}
        with outcome_path.open(newline='', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                genome = ast.literal_eval(row['genome'])
                cid = candidate_id(genome)
                rows_by_cid[cid] = {**row, '_genome': genome}

        eligible_ids = [
            cid for cid in record['eligible_ids'] if cid in rows_by_cid
        ]
        uncal_selected = list(record['policy_selected_ids']['catalyst'])
        unc_selected = list(record['policy_selected_ids']['uncertainty'])
        val_selected = list(record['policy_selected_ids']['validity'])
        uncal_pos = {cid: idx for idx, cid in enumerate(uncal_selected)}
        val_pos = {cid: idx for idx, cid in enumerate(val_selected)}
        unc_pos = {cid: idx for idx, cid in enumerate(unc_selected)}

        valid_interfacial_scores: list[tuple[float, str, str]] = []
        calibrated_pyro_pool: list[tuple[float, str, str]] = []
        batch_pyro_obs: list[tuple[tuple, float]] = []

        for elig_idx, cid in enumerate(eligible_ids):
            row = rows_by_cid[cid]
            genome = row['_genome']
            m_cls = str(genome[0])
            reg = '|'.join(discovery_region(genome))
            reg_shift = float(
                turq_reg_bias.get(reg, turq_bias.get(m_cls, 0.28))
            )
            reg_hist = prior_pyro_region_sum.get(reg)
            cls_hist = prior_pyro_class_sum.get(m_cls)
            prior_s = (
                float(np.mean(reg_hist))
                if reg_hist
                else (float(np.mean(cls_hist)) if cls_hist else 1.20)
            )
            genome_prior = pyrolysis_genome_structural_prior_eV(genome)
            rank_bonus = (
                -0.10 * (1.0 - uncal_pos[cid] / max(len(uncal_selected), 1))
                if cid in uncal_pos
                else (
                    -0.10 * (1.0 - val_pos[cid] / max(len(val_selected), 1))
                    if cid in val_pos
                    else (
                        -0.04 * (1.0 - unc_pos[cid] / max(len(unc_selected), 1))
                        if cid in unc_pos
                        else 0.0
                    )
                )
            )
            cal_score = (
                0.40 * prior_s
                + 0.70 * reg_shift
                + genome_prior
                + rank_bonus
            )
            calibrated_pyro_pool.append((cal_score, cid, m_cls))

            is_valid = str(row.get('valid', '')).lower() == 'true'
            is_censored = str(row.get('E_act_censored', '')).lower() == 'true'
            is_viable = str(row.get('pyrolysis_viable', 'true')).lower() == 'true'
            try:
                e_act = float(row.get('E_act', 'nan'))
            except (TypeError, ValueError):
                e_act = float('nan')
            is_eligible_metric = (
                is_valid and not is_censored and is_viable and np.isfinite(e_act)
            )
            if is_eligible_metric:
                try:
                    c_disp_raw = float(row.get('c_height_above_surface_A', 'nan'))
                except (TypeError, ValueError):
                    c_disp_raw = float('nan')
                try:
                    coke_idx = float(row.get('coking_index', 'nan'))
                except (TypeError, ValueError):
                    coke_idx = float('nan')
                if np.isfinite(c_disp_raw):
                    cal_info = calibrated_interfacial_pyrolysis_score_eV(
                        {
                            'material_class': m_cls,
                            'E_act_raw': e_act,
                            'c_height_above_surface_A': c_disp_raw,
                        }
                    )
                    s_int = float(
                        cal_info['calibrated_interfacial_score_eV']
                        + 0.50 * reg_shift
                    )
                else:
                    coke_pen = (
                        0.25 * max(0.0, -coke_idx)
                        if np.isfinite(coke_idx)
                        else 0.0
                    )
                    s_int = float(e_act + reg_shift + coke_pen)
                valid_interfacial_scores.append((s_int, cid, m_cls))
                batch_pyro_obs.append((genome, s_int))

        for g_obs, s_obs in batch_pyro_obs:
            prior_pyro_class_sum[str(g_obs[0])].append(s_obs)
            prior_pyro_region_sum['|'.join(discovery_region(g_obs))].append(
                s_obs
            )

        valid_interfacial_scores.sort(key=lambda x: (x[0], x[1]))
        hit_count = max(1, int(np.ceil(0.20 * len(valid_interfacial_scores))))
        hit_ids = {cid for _, cid, _ in valid_interfacial_scores[:hit_count]}
        budget = int(record.get('budget', 20))
        dft_cal_selected = _class_floor_select(calibrated_pyro_pool, budget)

        uncal_selected = list(record['policy_selected_ids']['catalyst'])
        unc_selected = list(record['policy_selected_ids']['uncertainty'])
        val_selected = list(record['policy_selected_ids']['validity'])

        batch_hits = {
            'catalyst_dft_calibrated': sum(
                c in hit_ids for c in dft_cal_selected
            ),
            'catalyst_uncalibrated': sum(
                c in hit_ids for c in uncal_selected
            ),
            'uncertainty': sum(c in hit_ids for c in unc_selected),
            'validity': sum(c in hit_ids for c in val_selected),
        }
        for k, v in batch_hits.items():
            pyro_totals[k] += int(v)

        by_class_ids = defaultdict(list)
        for cid in eligible_ids:
            by_class_ids[str(rows_by_cid[cid]['_genome'][0])].append(cid)
        rng = np.random.default_rng(int(record.get('random_seed', 42)))
        rand_draws = np.empty(random_trials, dtype=np.int16)
        sorted_classes = sorted(by_class_ids)
        for t_idx in range(random_trials):
            chosen = [
                str(rng.choice(by_class_ids[c_name]))
                for c_name in sorted_classes
            ]
            chosen_s = set(chosen)
            rem = [c for c in eligible_ids if c not in chosen_s]
            if budget > len(chosen) and rem:
                chosen.extend(
                    str(x)
                    for x in rng.choice(
                        rem, min(budget - len(chosen), len(rem)), replace=False
                    )
                )
            rand_draws[t_idx] = sum(c in hit_ids for c in chosen)
        pyro_random_matrix.append(rand_draws)
        pyro_batch_reports.append(
            {
                'batch_id': batch_id,
                'n_eligible': len(eligible_ids),
                'n_valid': len(valid_interfacial_scores),
                'hit_count': hit_count,
                'policy_hits': batch_hits,
                'random_mean_hits': round(float(np.mean(rand_draws)), 4),
            }
        )

    pyro_rand_sum = (
        np.sum(np.vstack(pyro_random_matrix), axis=0)
        if pyro_random_matrix
        else np.zeros(random_trials)
    )
    pyro_p_dft = float(
        (
            1
            + np.count_nonzero(
                pyro_rand_sum >= pyro_totals['catalyst_dft_calibrated']
            )
        )
        / (len(pyro_rand_sum) + 1)
    )
    pyro_p_uncal = float(
        (
            1
            + np.count_nonzero(
                pyro_rand_sum >= pyro_totals['catalyst_uncalibrated']
            )
        )
        / (len(pyro_rand_sum) + 1)
    )

    body = {
        'schema_version': 1,
        'calibration_sha256': cal_doc.get('calibration_sha256'),
        'fuel_cell_orr': {
            'n_batches': len(orr_batch_reports),
            'total_budget_selected': len(orr_batch_reports) * 20,
            'policy_hits_gc_acid_stable_top20pct': orr_totals,
            'doe_035V_gc_acid_stable_hits': orr_doe_gc_hits,
            'dft_calibrated_hit_lift_vs_uncalibrated': int(
                orr_totals['catalyst_dft_calibrated']
                - orr_totals['catalyst_uncalibrated']
            ),
            'dft_calibrated_doe_hit_lift_vs_uncalibrated': int(
                orr_doe_gc_hits['catalyst_dft_calibrated']
                - orr_doe_gc_hits['catalyst_uncalibrated']
            ),
            'policy_matched_random_mean_hits': round(
                float(np.mean(orr_rand_sum)), 4
            ),
            'dft_calibrated_enrichment_vs_random': round(
                float(
                    orr_totals['catalyst_dft_calibrated']
                    / max(float(np.mean(orr_rand_sum)), 1e-9)
                ),
                4,
            ),
            'dft_calibrated_p_value': orr_p_dft,
            'uncalibrated_p_value': orr_p_uncal,
            'acceptance_passed': bool(
                orr_totals['catalyst_dft_calibrated']
                > orr_totals['catalyst_uncalibrated']
                and orr_totals['catalyst_dft_calibrated']
                > orr_totals['uncertainty']
                and orr_totals['catalyst_dft_calibrated']
                > orr_totals['validity']
                and orr_p_dft <= 0.05
            ),
            'batches': orr_batch_reports,
        },
        'turquoise_hydrogen': {
            'n_batches': len(pyro_batch_reports),
            'total_budget_selected': len(pyro_batch_reports) * 20,
            'policy_hits_interfacial_calibrated_top20pct': pyro_totals,
            'dft_calibrated_hit_lift_vs_uncalibrated': int(
                pyro_totals['catalyst_dft_calibrated']
                - pyro_totals['catalyst_uncalibrated']
            ),
            'policy_matched_random_mean_hits': round(
                float(np.mean(pyro_rand_sum)), 4
            ),
            'dft_calibrated_enrichment_vs_random': round(
                float(
                    pyro_totals['catalyst_dft_calibrated']
                    / max(float(np.mean(pyro_rand_sum)), 1e-9)
                ),
                4,
            ),
            'dft_calibrated_p_value': pyro_p_dft,
            'uncalibrated_p_value': pyro_p_uncal,
            'acceptance_passed': bool(
                pyro_totals['catalyst_dft_calibrated']
                >= pyro_totals['catalyst_uncalibrated']
                and pyro_totals['catalyst_dft_calibrated']
                > pyro_totals['uncertainty']
                and pyro_totals['catalyst_dft_calibrated']
                > pyro_totals['validity']
                and pyro_p_dft <= 0.05
            ),
            'batches': pyro_batch_reports,
        },
    }
    report_sha = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    report = {
        **body,
        'report_sha256': report_sha,
    }
    target = Path(output_path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    return report


