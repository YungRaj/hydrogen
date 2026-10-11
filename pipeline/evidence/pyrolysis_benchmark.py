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
    values: dict[str, float] = {}
    for name in ('E_act_raw', 'dE_H', 'dE_C'):
        try:
            value = float(row[name])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f'missing finite {name}') from exc
        if not math.isfinite(value):
            raise ValueError(f'missing finite {name}')
        values[name] = value
    return max(values['E_act_raw'], -2.0 * values['dE_H'], -values['dE_C'], 0.0)


def roc_auc(labels: list[bool], scores: list[float]) -> float:
    """Compute tie-aware ROC AUC where larger scores predict ``True``.

    Args:
        labels: Binary ground-truth values.
        scores: Aligned prediction scores.

    Returns:
        Area under the receiver-operating-characteristic curve.
    """
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0 or len(labels) != len(scores):
        raise ValueError('ROC AUC requires aligned positive and negative examples')
    wins = 0.0
    for index, positive in enumerate(labels):
        if not positive:
            continue
        for other, negative in enumerate(labels):
            if negative:
                continue
            if scores[index] > scores[other]:
                wins += 1.0
            elif scores[index] == scores[other]:
                wins += 0.5
    return wins / (positives * negatives)


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
    sources: dict[str, set[str]] = {'active': set(), 'poor_or_poisoning': set()}
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
        for field in (
            'condition_id',
            'catalyst_identity',
            'reactor_context',
            'source_id',
            'source_family',
            'label_basis',
            'source_url',
        ):
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
        if result.get('valid') not in (True, 'True', 'true', 1, '1'):
            failures.append(f"{condition['condition_id']}: invalid screening result")
            continue
        if result.get('pyrolysis_viable') not in (True, 'True', 'true', 1, '1'):
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
    args = parser.parse_args()
    lock_errors = verify_analysis_lock(args.lock)
    if lock_errors:
        parser.error('; '.join(lock_errors))
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
    target = Path(args.report)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(target)
    print(json.dumps({key: value for key, value in report.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
