"""Independent calibration/holdout scoring for external reactor models."""

from __future__ import annotations

from collections.abc import Mapping
import json
import math
from pathlib import Path
from typing import Any, cast


def score_holdout(path: str | Path, calibration: Mapping[str, Any]) -> dict[str, Any]:
    """Calculate held-out error from raw prediction/observation records.

    Args:
        path: Filesystem path to the input or output artifact.
        calibration: Mapping supplying calibration.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    try:
        payload = json.loads(Path(path).read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f'validation records are unreadable: {exc}') from exc
    training_ids = {str(x) for x in calibration.get('training_ids', [])}
    validation_ids = {str(x) for x in calibration.get('validation_ids', [])}
    if not training_ids or not validation_ids or training_ids & validation_ids:
        raise ValueError('calibration and holdout IDs must be nonempty and disjoint')
    raw_records = payload.get('records', [])
    records = cast(list[Any], raw_records) if isinstance(raw_records, list) else []
    if not payload.get('source') or payload.get('source') != calibration.get('source'):
        raise ValueError('validation-record source must match calibration source')
    by_id: dict[str, dict[str, Any]] = {}
    for raw_row in records:
        if isinstance(raw_row, dict):
            typed_row = cast(dict[str, Any], raw_row)
            by_id[str(typed_row.get('id'))] = typed_row
    missing = sorted((training_ids | validation_ids) - set(by_id.keys()))
    if missing:
        raise ValueError(f'validation records missing declared IDs: {missing}')
    metric = str(calibration.get('metric', 'relative_rmse'))
    raw_threshold = calibration.get('acceptance_threshold')
    if raw_threshold is None:
        raise ValueError('invalid validation metric or acceptance threshold')
    threshold = float(raw_threshold)
    if (
        metric not in {'rmse', 'mae', 'relative_rmse'}
        or not math.isfinite(threshold)
        or threshold < 0
    ):
        raise ValueError('invalid validation metric or acceptance threshold')
    errors: list[float] = []
    relative: list[float] = []
    for record_id in sorted(validation_ids):
        row = by_id[record_id]
        predicted, observed = float(row['predicted']), float(row['observed'])
        if not math.isfinite(predicted) or not math.isfinite(observed):
            raise ValueError(f'non-finite holdout record: {record_id}')
        errors.append(predicted - observed)
        # 1e-12 denominator floor guards relative_rmse against division by zero on zero-valued observations
        relative.append((predicted - observed) / max(abs(observed), 1e-12))
    if metric == 'mae':
        error = sum(abs(value) for value in errors) / len(errors)
    elif metric == 'rmse':
        error = math.sqrt(sum(value * value for value in errors) / len(errors))
    else:
        error = math.sqrt(sum(value * value for value in relative) / len(relative))
    return {
        'metric': metric,
        'holdout_error': error,
        'acceptance_threshold': threshold,
        'passed': error <= threshold,
        'training_count': len(training_ids),
        'holdout_count': len(validation_ids),
        'record_source': payload.get('source'),
    }
