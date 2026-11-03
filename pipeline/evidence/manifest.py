"""Hash-verified scientific evidence manifests."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, cast

from pipeline.io.hashing import file_sha256

# Schema version 2 requires per-record SHA-256 content digests and explicit status gates.
MANIFEST_SCHEMA_VERSION = 2

EVIDENCE_KEYS = (
    'converged_dft',
    'converged_orr_dft',
    'ntec_control_pair',
    'measured_reactor',
    'measured_deactivation',
    'measured_mea',
    'measured_durability',
    'hydrogen_impurity_test',
    'time_split_benchmark',
    'curated_prior_art_source',
)

REQUIRED_STATUS = {
    'converged_dft': 'converged',
    'converged_orr_dft': 'converged',
    'ntec_control_pair': 'measured',
    'measured_reactor': 'measured',
    'measured_deactivation': 'measured',
    'measured_mea': 'measured',
    'measured_durability': 'measured',
    'hydrogen_impurity_test': 'measured',
    'time_split_benchmark': 'valid',
    'curated_prior_art_source': 'curated',
}


def verify_evidence_manifest(path: str | Path) -> dict[str, Any]:
    """Verify every evidence record against an immutable source artifact.

    Args:
        path: Filesystem path to the input or output artifact.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    manifest_path = Path(path)
    if not manifest_path.is_file():
        return {
            'valid': False,
            'errors': ['evidence_manifest_missing'],
            'counts': {key: 0 for key in EVIDENCE_KEYS},
        }
    try:
        payload = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {
            'valid': False,
            'errors': [f'evidence_manifest_invalid:{exc}'],
            'counts': {key: 0 for key in EVIDENCE_KEYS},
        }
    if payload.get('schema_version') != MANIFEST_SCHEMA_VERSION or not isinstance(
        payload.get('records'), dict
    ):
        return {
            'valid': False,
            'errors': ['evidence_manifest_schema_mismatch'],
            'counts': {key: 0 for key in EVIDENCE_KEYS},
        }

    errors: list[str] = []
    verified: dict[str, list[dict[str, Any]]] = {}
    counts: dict[str, int] = {}
    records_map = cast(dict[str, Any], payload['records'])
    for key in EVIDENCE_KEYS:
        raw_records = records_map.get(key, [])
        if not isinstance(raw_records, list):
            errors.append(f'{key}:records_not_list')
            records: list[Any] = []
        else:
            records = cast(list[Any], raw_records)
        accepted: list[dict[str, Any]] = []
        seen: set[tuple[str, str, str]] = set()
        for index, raw_record in enumerate(records):
            label = f'{key}[{index}]'
            if not isinstance(raw_record, dict):
                errors.append(f'{label}:not_object')
                continue
            record = cast(dict[str, Any], raw_record)
            required = (
                'candidate_id',
                'source_path',
                'sha256',
                'protocol_id',
                'status',
            )
            missing = [field for field in required if not record.get(field)]
            if missing:
                errors.append(f'{label}:missing:{",".join(missing)}')
                continue
            if record['status'] != REQUIRED_STATUS[key]:
                errors.append(f'{label}:status_must_be:{REQUIRED_STATUS[key]}')
                continue
            identity = (
                str(record['candidate_id']),
                str(record['source_path']),
                str(record['protocol_id']),
            )
            if identity in seen:
                errors.append(f'{label}:duplicate')
                continue
            seen.add(identity)
            source = Path(str(record['source_path']))
            if not source.is_absolute():
                source = manifest_path.parent / source
            if not source.is_file():
                errors.append(f'{label}:source_missing')
                continue
            if file_sha256(source) != str(record['sha256']).lower():
                errors.append(f'{label}:checksum_mismatch')
                continue
            accepted.append({**record, 'resolved_source': str(source.resolve())})
        verified[key] = accepted
        counts[key] = len(accepted)
    return {
        'valid': not errors,
        'errors': errors,
        'counts': counts,
        'verified_records': verified,
        'schema_version': MANIFEST_SCHEMA_VERSION,
    }
