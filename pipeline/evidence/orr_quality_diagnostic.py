"""Diagnose ORR censoring and screening failures without changing policy."""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path('results/prospective_search')
OUTPUT = ROOT / 'orr_quality_diagnostic.json'


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _failure_family(error: str) -> str:
    lowered = error.lower()
    for stage in ('ooh', 'oh', 'o', 'clean'):
        if f'invalid {stage} initial geometry: atomic_overlap' in lowered:
            return f'{stage}_atomic_overlap'
        if f'unconverged {stage} relaxation' in lowered:
            return f'{stage}_unconverged'
    if 'unphysical' in lowered:
        return 'unphysical_energy'
    if 'timeout' in lowered:
        return 'timeout'
    return 'other' if error else 'unspecified'


def _binding_bin(raw: str) -> str:
    value = abs(float(raw))
    if value <= 1.0:
        return 'abs_le_1_eV'
    if value <= 2.0:
        return 'abs_1_to_2_eV'
    if value <= 3.0:
        return 'abs_2_to_3_eV'
    if value <= 5.0:
        return 'abs_3_to_5_eV'
    return 'abs_gt_5_eV'


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f'.{path.name}.', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w') as stream:
            json.dump(payload, stream, indent=2, sort_keys=True)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def build_diagnostic(output: Path = OUTPUT) -> Path:
    """Summarize failure and metric-floor rates in confirmed ORR evidence.

    Args:
        output: Destination for the immutable diagnostic JSON.

    Returns:
        Path to the written diagnostic artifact.
    """
    final_path = ROOT / 'confirmatory_analysis.json'
    final = json.loads(final_path.read_text())
    if final.get('status') != 'complete':
        raise RuntimeError('completed confirmatory evidence is required')
    sources = {str(final_path): _sha256(final_path)}
    totals: Counter[str] = Counter()
    by_class: dict[str, Counter[str]] = defaultdict(Counter)
    failures: Counter[str] = Counter()
    censored_steps: Counter[str] = Counter()
    censored_binding: Counter[str] = Counter()
    for batch in final['batch_ids']:
        analysis_path = ROOT / 'batches' / str(batch) / 'analysis.json'
        analysis = json.loads(analysis_path.read_text())
        outcome_path = Path('results/fuel_cell/prospective') / f'{batch}_orr.csv'
        if analysis['outcome_sha256']['fuel_cell_orr'] != _sha256(outcome_path):
            raise RuntimeError(f'outcome checksum mismatch: {batch}')
        sources[str(analysis_path)] = _sha256(analysis_path)
        sources[str(outcome_path)] = _sha256(outcome_path)
        with outcome_path.open(newline='') as stream:
            for row in csv.DictReader(stream):
                material_class = str(ast.literal_eval(row['genome'])[0])
                valid = row['valid'].lower() == 'true'
                totals['screened'] += 1
                by_class[material_class]['screened'] += 1
                if not valid:
                    totals['invalid'] += 1
                    by_class[material_class]['invalid'] += 1
                    failures[_failure_family(row.get('error', ''))] += 1
                    continue
                totals['valid'] += 1
                by_class[material_class]['valid'] += 1
                value = float(row['orr_overpotential_V'])
                censored = math.isfinite(value) and value == 3.0
                if censored:
                    totals['censored_valid'] += 1
                    by_class[material_class]['censored_valid'] += 1
                    censored_steps[str(row.get('rate_determining_step') or 'unknown')] += 1
                    raw_binding = str(row.get('binding_strength') or '')
                    try:
                        censored_binding[_binding_bin(raw_binding)] += 1
                    except ValueError:
                        censored_binding['unknown'] += 1
                else:
                    totals['uncensored_valid'] += 1
                    by_class[material_class]['uncensored_valid'] += 1
    screened = totals['screened']
    valid = totals['valid']
    report = {
        'schema_version': 1,
        'interpretation': 'descriptive_only_no_policy_or_verdict_change',
        'metric_floor_definition': 'valid and orr_overpotential_V == 3.0',
        'totals': {
            **dict(totals),
            'invalid_fraction': totals['invalid'] / screened,
            'censored_fraction_of_valid': totals['censored_valid'] / valid,
            'scoreable_fraction': totals['uncensored_valid'] / screened,
        },
        'by_material_class': {
            name: {
                **dict(counts),
                'invalid_fraction': counts['invalid'] / counts['screened'],
                'censored_fraction_of_valid': (
                    counts['censored_valid'] / counts['valid'] if counts['valid'] else None
                ),
            }
            for name, counts in sorted(by_class.items())
        },
        'invalid_failure_families': dict(failures.most_common()),
        'censored_rate_determining_steps': dict(censored_steps.most_common()),
        'censored_binding_strength': dict(censored_binding.most_common()),
        'source_sha256': sources,
    }
    body = json.dumps(report, sort_keys=True, separators=(',', ':')).encode()
    report['diagnostic_sha256'] = hashlib.sha256(body).hexdigest()
    _atomic_json(output, report)
    return output


def verify_diagnostic(path: Path = OUTPUT) -> bool:
    """Verify the diagnostic digest and every source checksum.

    Args:
        path: Diagnostic artifact to check.

    Returns:
        True only when the diagnostic and all evidence inputs are unchanged.
    """
    report = json.loads(path.read_text())
    expected = report.pop('diagnostic_sha256')
    actual = hashlib.sha256(
        json.dumps(report, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    return expected == actual and all(
        Path(source).is_file() and _sha256(Path(source)) == digest
        for source, digest in report['source_sha256'].items()
    )


if __name__ == '__main__':
    target = build_diagnostic()
    if not verify_diagnostic(target):
        raise SystemExit('generated ORR diagnostic failed verification')
    print(target)
