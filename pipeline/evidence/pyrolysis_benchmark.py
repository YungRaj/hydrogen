"""Fail-closed contracts for the external pyrolysis descriptor benchmark."""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
from typing import Any

from pipeline.search.discovery import candidate_id

MIN_PER_LABEL = 12


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
    groups: dict[str, set[str]] = {'active': set(), 'poor_or_poisoning': set()}
    for index, row in enumerate(records):
        prefix = f'record {index}'
        if not isinstance(row, dict):
            errors.append(f'{prefix} must be an object')
            continue
        label = str(row.get('label', ''))
        if label not in groups:
            errors.append(f'{prefix} has invalid label')
        if row.get('exact_encoding') is not True:
            errors.append(f'{prefix} is not an exact encoding')
        for field in (
            'condition_id',
            'catalyst_identity',
            'reactor_context',
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
        try:
            genome = ast.literal_eval(str(row['genome']))
            candidate_id(genome)
        except (KeyError, SyntaxError, ValueError, TypeError):
            errors.append(f'{prefix} has invalid genome')
            continue
        if label in groups and isinstance(row.get('source_family'), str):
            groups[label].add(row['source_family'])
    for label, source_groups in groups.items():
        if len(source_groups) < MIN_PER_LABEL:
            errors.append(
                f'{label} has {len(source_groups)} independent source families; '
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
