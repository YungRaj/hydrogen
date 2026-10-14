"""Fail-closed contracts for the external pyrolysis descriptor benchmark."""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

from pipeline.search.discovery import candidate_id

MIN_PER_LABEL = 12
LABELS = ('active', 'poor_or_poisoning')
REQUIRED_TEXT_FIELDS = (
    'condition_id',
    'catalyst_identity',
    'reactor_context',
    'source_id',
    'source_family',
    'label_basis',
    'source_url',
)
TRUE_VALUES = (True, 'True', 'true', 1, '1')
EV_PER_KJ_MOL = 1.0 / 96.4853321233
K_B_EV_K = 8.617333262145e-5
H2_H_298_MINUS_H_0_KJ_MOL = 8.468
H2_SHOMATE = (
    (298.0, 1000.0, (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974)),
    (1000.0, 2500.0, (18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133)),
)


def effective_barrier_eV(barriers: list[float], temperature_K: float) -> float:
    """Return the barrier whose Arrhenius factor equals the mean site rate.

    Args:
        barriers: Non-negative replica barriers in eV.
        temperature_K: Absolute condition temperature.

    Returns:
        Numerically stable rate-averaged effective barrier in eV.
    """
    if not barriers or temperature_K <= 0.0:
        raise ValueError('rate averaging requires barriers and positive temperature')
    values = [float(value) for value in barriers]
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError('rate averaging requires finite non-negative barriers')
    minimum = min(values)
    mean_relative_rate = sum(
        math.exp(-(value - minimum) / (K_B_EV_K * temperature_K))
        for value in values
    ) / len(values)
    return minimum - K_B_EV_K * temperature_K * math.log(mean_relative_rate)


def h2_chemical_potential_eV(
    temperature_K: float, pressure_bar: float = 1.0
) -> float:
    """Return H2 ideal-gas chemical potential relative to 0 K at 1 bar.

    Args:
        temperature_K: Absolute gas temperature within 298--2500 K.
        pressure_bar: Ideal-gas H2 partial pressure in bar.

    Returns:
        Thermal enthalpy, entropy, and pressure correction in eV per H2.
    """
    if pressure_bar <= 0.0:
        raise ValueError('H2 pressure must be positive')
    selected = next(
        (coefficients for low, high, coefficients in H2_SHOMATE
         if low <= temperature_K <= high),
        None,
    )
    if selected is None:
        raise ValueError('H2 Shomate correction requires 298 <= T <= 2500 K')
    a, b, c, d, e, f, g = selected
    t = temperature_K / 1000.0
    enthalpy_kj_mol = (
        a * t + b * t**2 / 2.0 + c * t**3 / 3.0 + d * t**4 / 4.0
        - e / t + f + H2_H_298_MINUS_H_0_KJ_MOL
    )
    entropy_j_mol_k = (
        a * math.log(t) + b * t + c * t**2 / 2.0 + d * t**3 / 3.0
        - e / (2.0 * t**2) + g
    )
    standard_eV = (
        enthalpy_kj_mol - temperature_K * entropy_j_mol_k / 1000.0
    ) * EV_PER_KJ_MOL
    return standard_eV + K_B_EV_K * temperature_K * math.log(pressure_bar)


def rate_free_energy_span_eV(row: dict[str, Any]) -> float:
    """Return the rate-averaged, H2-free-energy-corrected molten score.

    Args:
        row: Completed ensemble screening row for a promoted molten alloy.

    Returns:
        Effective turnover-limiting free-energy span in eV.

    Raises:
        ValueError: If ensemble data are absent or the BEP scope is unsupported.
    """
    try:
        genome = ast.literal_eval(str(row['genome']))
        temperature_K = float(row['temperature_K'])
        replicas = json.loads(str(row['molten_replicas_json']))
    except (KeyError, TypeError, ValueError, SyntaxError, json.JSONDecodeError) as exc:
        raise ValueError('missing molten rate/free-energy input') from exc
    if genome[0] != 'MoltenMetal' or genome[2] == 'None':
        raise ValueError('pure-host molten BEP barrier is outside validated scope')
    h2_mu = h2_chemical_potential_eV(temperature_K)
    spans: list[float] = []
    for replica in replicas:
        try:
            activation = max(0.0, float(replica['E_act_raw']))
            hydrogen_release = -2.0 * float(replica['dE_H']) + h2_mu
            carbon_release = -float(replica['dE_C'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('missing molten rate/free-energy input') from exc
        span = max(activation, hydrogen_release, carbon_release, 0.0)
        if not math.isfinite(span):
            raise ValueError('missing molten rate/free-energy input')
        spans.append(span)
    return effective_barrier_eV(spans, temperature_K)


def energetic_span_eV(row: dict[str, Any]) -> float:
    """Return the largest uphill step in the minimal CH4 -> C + 2 H2 cycle.

    The score has no fitted coefficients.  It compares methane activation with
    the thermodynamic costs of releasing two adsorbed H atoms as H2 and one
    adsorbed C atom as solid carbon.  Lower values are better.

    Args:
        row: Completed pyrolysis screening result.

    Returns:
        Two-sided energetic-span score in eV.

    Raises:
        ValueError: If a required value is missing, non-finite, or censored.
    """
    names = ('E_act_raw', 'dE_H', 'dE_C')
    try:
        e_act, d_e_h, d_e_c = (float(row[name]) for name in names)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('missing finite energetic-span input') from exc
    if not all(math.isfinite(value) for value in (e_act, d_e_h, d_e_c)):
        raise ValueError('missing finite energetic-span input')
    return max(e_act, -2.0 * d_e_h, -d_e_c, 0.0)


def roc_auc(labels: list[bool], scores: list[float]) -> float:
    """Compute tie-aware ROC AUC where larger scores predict ``True``.

    Args:
        labels: Binary ground-truth values.
        scores: Aligned prediction scores.

    Returns:
        Area under the receiver-operating-characteristic curve.
    """
    if len(labels) != len(scores):
        raise ValueError('ROC AUC requires aligned positive and negative examples')
    positives = [score for label, score in zip(labels, scores) if label]
    negatives = [score for label, score in zip(labels, scores) if not label]
    if not positives or not negatives:
        raise ValueError('ROC AUC requires aligned positive and negative examples')
    wins = sum(
        float(positive > negative) + 0.5 * float(positive == negative)
        for positive in positives
        for negative in negatives
    )
    return wins / (len(positives) * len(negatives))


def grouped_bootstrap_auc(
    rows: list[dict[str, Any]],
    score_key: str,
    resamples: int = 10_000,
    seed: int = 20261010,
) -> dict[str, float]:
    """Estimate AUC and a percentile interval by resampling source families.

    Args:
        rows: Labelled score rows containing source family and label.
        score_key: Lower-is-better score to evaluate.
        resamples: Number of deterministic bootstrap resamples.
        seed: Reproducible pseudo-random seed.

    Returns:
        Point AUC and 95% percentile bounds.
    """
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(str(row['source_family']), []).append(row)
    group_rows = list(groups.values())

    def auc(sample: list[dict[str, Any]]) -> float:
        labels = [row['label'] == 'active' for row in sample]
        scores = [-float(row[score_key]) for row in sample]
        return roc_auc(labels, scores)

    point = auc(rows)
    generator = random.Random(seed)
    estimates: list[float] = []
    while len(estimates) < resamples:
        sample = [
            row
            for _ in group_rows
            for row in generator.choice(group_rows)
        ]
        try:
            estimates.append(auc(sample))
        except ValueError:
            continue
    estimates.sort()
    return {
        'auc': point,
        'ci95_low': estimates[int(0.025 * resamples)],
        'ci95_high': estimates[int(0.975 * resamples) - 1],
    }


def dataset_digest(document: dict[str, Any]) -> str:
    """Return the content identity of a benchmark document without its digest.

    Args:
        document: Parsed benchmark dataset.

    Returns:
        Lowercase SHA-256 digest of canonical JSON content.
    """
    body = {key: value for key, value in document.items() if key != 'dataset_sha256'}
    return hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()


def validate_benchmark(document: dict[str, Any], require_locked: bool = True) -> list[str]:
    """Validate labels, provenance, independence, and exact genome mappings.

    Args:
        document: Parsed external benchmark dataset.
        require_locked: Whether readiness requires a matching immutable digest.

    Returns:
        Human-readable contract violations; an empty list means scoreable.
    """
    errors: list[str] = []
    if document.get('schema_version') != 1:
        errors.append('schema_version must be 1')
    records = document.get('records')
    if not isinstance(records, list):
        return errors + ['records must be a list']
    condition_ids: set[str] = set()
    sources: dict[str, set[str]] = {label: set() for label in LABELS}
    for index, row in enumerate(records):
        prefix = f'record {index}'
        if not isinstance(row, dict):
            errors.append(f'{prefix} must be an object')
            continue
        label = str(row.get('label', ''))
        if label not in sources:
            errors.append(f'{prefix} has invalid label')
        if row.get('exact_encoding') is not True:
            errors.append(f'{prefix} is not an exact encoding')
        for field in REQUIRED_TEXT_FIELDS:
            if not isinstance(row.get(field), str) or not row[field].strip():
                errors.append(f'{prefix} missing {field}')
        condition_id = row.get('condition_id')
        if isinstance(condition_id, str) and condition_id.strip():
            if condition_id in condition_ids:
                errors.append(f'{prefix} duplicates condition_id')
            condition_ids.add(condition_id)
        try:
            temperature_K = float(row['temperature_K'])
            if not math.isfinite(temperature_K) or temperature_K < 700.0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f'{prefix} temperature_K must be finite and >= 700 K')
        if not str(row.get('source_url', '')).startswith(('https://doi.org/', 'https://')):
            errors.append(f'{prefix} source_url is not stable HTTPS provenance')
        source_id = str(row.get('source_id', ''))
        if source_id.startswith('doi:'):
            expected_url = 'https://doi.org/' + source_id.removeprefix('doi:')
            if row.get('source_url') != expected_url:
                errors.append(f'{prefix} source_id and source_url disagree')
        try:
            genome = ast.literal_eval(str(row['genome']))
            candidate_id(genome)
        except (KeyError, SyntaxError, ValueError, TypeError):
            errors.append(f'{prefix} has invalid genome')
            continue
        if genome[0] == 'MoltenMetal':
            from pipeline.search.design_space import validate_molten_metal

            if not validate_molten_metal(genome):
                errors.append(f'{prefix} molten genome is outside the design space')
        if label in sources and isinstance(row.get('source_id'), str):
            sources[label].add(row['source_id'])
    for label, independent_sources in sources.items():
        if len(independent_sources) < MIN_PER_LABEL:
            errors.append(
                f'{label} has {len(independent_sources)} independent primary sources; '
                f'{MIN_PER_LABEL} required'
            )
    if require_locked:
        if document.get('locked') is not True:
            errors.append('dataset is not locked')
        if document.get('dataset_sha256') != dataset_digest(document):
            errors.append('dataset digest mismatch')
    return errors


def load_scoreable_benchmark(path: str | Path) -> dict[str, Any]:
    """Load a benchmark only after every preregistered curation gate passes.

    Args:
        path: JSON benchmark dataset to load.

    Returns:
        Parsed, checksum-verified benchmark document.

    Raises:
        ValueError: If any provenance, independence, encoding, or lock gate fails.
    """
    document = json.loads(Path(path).read_text())
    errors = validate_benchmark(document)
    if errors:
        raise ValueError('pyrolysis benchmark is not scoreable: ' + '; '.join(errors))
    return document


def verify_analysis_lock(path: str | Path, root: str | Path = '.') -> list[str]:
    """Return mismatches between an analysis lock and its frozen files.

    Args:
        path: JSON lock document.
        root: Repository root against which locked paths are resolved.

    Returns:
        Missing or content-mismatched locked-file errors.
    """
    lock = json.loads(Path(path).read_text())
    errors: list[str] = []
    for relative, expected in lock.get('sha256', {}).items():
        source = Path(root) / relative
        if not source.is_file():
            errors.append(f'missing locked file: {relative}')
            continue
        actual = hashlib.sha256(source.read_bytes()).hexdigest()
        if actual != expected:
            errors.append(f'locked file changed: {relative}')
    if not lock.get('sha256'):
        errors.append('lock contains no file hashes')
    return errors


def analyze_screening_results(
    benchmark: dict[str, Any], screening_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    """Join locked conditions to screening results and run the frozen analysis.

    Args:
        benchmark: Checksum-verified benchmark document.
        screening_rows: Screening results keyed by their serialized genome.

    Returns:
        Complete benchmark report including the preregistered pass verdict.

    Raises:
        ValueError: If any condition lacks an eligible screening result.
    """
    by_genome = {str(row.get('genome')): row for row in screening_rows}
    scored: list[dict[str, Any]] = []
    failures: list[str] = []
    for condition in benchmark['records']:
        result = by_genome.get(condition['genome'])
        if result is None:
            failures.append(f"{condition['condition_id']}: missing screening result")
            continue
        if result.get('valid') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: invalid screening result")
            continue
        if result.get('pyrolysis_viable') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: non-viable screening result")
            continue
        try:
            span = energetic_span_eV(result)
            current = float(result['E_act_raw'])
        except (KeyError, TypeError, ValueError) as exc:
            failures.append(f"{condition['condition_id']}: {exc}")
            continue
        scored.append({
            'condition_id': condition['condition_id'],
            'source_family': condition['source_family'],
            'label': condition['label'],
            'genome': condition['genome'],
            'E_act_eV': current,
            'energetic_span_eV': span,
            'limiting_step': max(
                (
                    ('methane_activation', current),
                    ('hydrogen_release', -2.0 * float(result['dE_H'])),
                    ('carbon_release', -float(result['dE_C'])),
                ),
                key=lambda item: item[1],
            )[0],
        })
    if failures:
        raise ValueError('benchmark screening is incomplete: ' + '; '.join(failures))
    current_stats = grouped_bootstrap_auc(scored, 'E_act_eV')
    span_stats = grouped_bootstrap_auc(scored, 'energetic_span_eV')
    improvement = span_stats['auc'] - current_stats['auc']
    passed = (
        span_stats['auc'] >= 0.70
        and span_stats['ci95_low'] > 0.50
        and improvement >= 0.10
    )
    return {
        'schema_version': 1,
        'dataset_sha256': benchmark['dataset_sha256'],
        'conditions_scored': len(scored),
        'current_E_act': current_stats,
        'energetic_span': span_stats,
        'auc_improvement': improvement,
        'acceptance_passed': passed,
        'records': scored,
    }


def main() -> None:
    """Screen or analyze the locked external benchmark."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('screen', 'analyze'))
    parser.add_argument(
        '--dataset', default='docs/evidence/pyrolysis_benchmark_dataset.json'
    )
    parser.add_argument(
        '--screening-csv', default='results/screening/pyrolysis_benchmark.csv'
    )
    parser.add_argument(
        '--report', default='results/evidence/pyrolysis_descriptor_benchmark.json'
    )
    parser.add_argument(
        '--lock', default='docs/evidence/pyrolysis_descriptor_benchmark_lock.json'
    )
    parser.add_argument(
        '--exploratory', action='store_true',
        help='allow code that differs from the lock (e.g. another screening '
             'model); the report is labelled exploratory and cannot pass')
    args = parser.parse_args()
    lock_errors = verify_analysis_lock(args.lock)
    if lock_errors and not args.exploratory:
        parser.error(
            '; '.join(lock_errors) + ' -- this analysis is frozen; re-run it from '
            'the commit recorded in docs/evidence/lock_provenance.json'
        )
    benchmark = load_scoreable_benchmark(args.dataset)
    if args.action == 'screen':
        from pipeline.screening.surface_screener import run_screening

        genomes = list(dict.fromkeys(
            ast.literal_eval(row['genome']) for row in benchmark['records']
        ))
        run_screening(genomes, Path(args.screening_csv).name)
        return

    import pandas as pd

    rows = pd.read_csv(args.screening_csv).to_dict(orient='records')
    report = analyze_screening_results(benchmark, rows)
    report['screening_models'] = sorted(
        {str(row['screening_model']) for row in rows if row.get('screening_model')})
    if lock_errors:
        # Not the preregistered analysis: record why and never claim a pass.
        report['exploratory_modified_implementation'] = True
        report['lock_errors'] = lock_errors
        report['acceptance_passed'] = False
    target = Path(args.report)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(target)
    print(json.dumps({key: value for key, value in report.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
