"""Hash-verified scientific evidence manifests."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


EVIDENCE_KEYS = (
    'converged_dft', 'converged_orr_dft', 'ntec_control_pair',
    'measured_reactor', 'measured_deactivation', 'measured_mea',
    'measured_durability', 'hydrogen_impurity_test',
    'time_split_benchmark', 'curated_prior_art_source',
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


def file_sha256(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def verify_evidence_manifest(path: str | Path) -> dict:
    """Verify every evidence record against an immutable source artifact."""
    manifest_path = Path(path)
    if not manifest_path.is_file():
        return {'valid': False, 'errors': ['evidence_manifest_missing'],
                'counts': {key: 0 for key in EVIDENCE_KEYS}}
    try:
        payload = json.loads(manifest_path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        return {'valid': False, 'errors': [f'evidence_manifest_invalid:{exc}'],
                'counts': {key: 0 for key in EVIDENCE_KEYS}}
    if payload.get('schema_version') != 2 or not isinstance(payload.get('records'), dict):
        return {'valid': False, 'errors': ['evidence_manifest_schema_mismatch'],
                'counts': {key: 0 for key in EVIDENCE_KEYS}}

    errors, verified, counts = [], {}, {}
    for key in EVIDENCE_KEYS:
        records = payload['records'].get(key, [])
        if not isinstance(records, list):
            errors.append(f'{key}:records_not_list')
            records = []
        accepted = []
        seen = set()
        for index, record in enumerate(records):
            label = f'{key}[{index}]'
            if not isinstance(record, dict):
                errors.append(f'{label}:not_object')
                continue
            required = ('candidate_id', 'source_path', 'sha256',
                        'protocol_id', 'status')
            missing = [field for field in required if not record.get(field)]
            if missing:
                errors.append(f'{label}:missing:{",".join(missing)}')
                continue
            if record['status'] != REQUIRED_STATUS[key]:
                errors.append(
                    f'{label}:status_must_be:{REQUIRED_STATUS[key]}')
                continue
            identity = (record['candidate_id'], record['source_path'],
                        record['protocol_id'])
            if identity in seen:
                errors.append(f'{label}:duplicate')
                continue
            seen.add(identity)
            source = Path(record['source_path'])
            if not source.is_absolute():
                source = manifest_path.parent / source
            if not source.is_file():
                errors.append(f'{label}:source_missing')
                continue
            if file_sha256(source) != str(record['sha256']).lower():
                errors.append(f'{label}:checksum_mismatch')
                continue
            accepted.append(dict(record, resolved_source=str(source.resolve())))
        verified[key] = accepted
        counts[key] = len(accepted)
    return {'valid': not errors, 'errors': errors, 'counts': counts,
            'verified_records': verified, 'schema_version': 2}
