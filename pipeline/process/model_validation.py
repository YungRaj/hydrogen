"""Independent calibration/holdout scoring for external reactor models."""

from __future__ import annotations

import json
import math
from pathlib import Path


def score_holdout(path: str | Path, calibration: dict) -> dict:
    """Calculate held-out error from raw prediction/observation records."""
    try:
        payload = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'validation records are unreadable: {exc}') from exc
    training_ids = set(calibration.get('training_ids', []))
    validation_ids = set(calibration.get('validation_ids', []))
    if not training_ids or not validation_ids or training_ids & validation_ids:
        raise ValueError('calibration and holdout IDs must be nonempty and disjoint')
    records = payload.get('records', [])
    if not payload.get('source') or payload.get('source') != calibration.get('source'):
        raise ValueError('validation-record source must match calibration source')
    by_id = {str(row.get('id')): row for row in records if isinstance(row, dict)}
    missing = sorted((training_ids | validation_ids) - set(by_id))
    if missing:
        raise ValueError(f'validation records missing declared IDs: {missing}')
    metric = calibration.get('metric', 'relative_rmse')
    threshold = float(calibration.get('acceptance_threshold'))
    if metric not in {'rmse', 'mae', 'relative_rmse'} or not math.isfinite(
            threshold) or threshold < 0:
        raise ValueError('invalid validation metric or acceptance threshold')
    errors = []
    relative = []
    for record_id in sorted(validation_ids):
        row = by_id[record_id]
        predicted, observed = float(row['predicted']), float(row['observed'])
        if not math.isfinite(predicted) or not math.isfinite(observed):
            raise ValueError(f'non-finite holdout record: {record_id}')
        errors.append(predicted - observed)
        relative.append((predicted - observed) / max(abs(observed), 1e-12))
    if metric == 'mae':
        error = sum(abs(value) for value in errors) / len(errors)
    elif metric == 'rmse':
        error = math.sqrt(sum(value * value for value in errors) / len(errors))
    else:
        error = math.sqrt(sum(value * value for value in relative) / len(relative))
    return {
        'metric': metric, 'holdout_error': error,
        'acceptance_threshold': threshold, 'passed': error <= threshold,
        'training_count': len(training_ids), 'holdout_count': len(validation_ids),
        'record_source': payload.get('source'),
    }
