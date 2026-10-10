"""Build a checksum-locked ORR DFT slate from confirmatory evidence."""

from __future__ import annotations

import ast
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

from pipeline.search.discovery import candidate_id

ROOT = Path('results/prospective_search')
OUTPUT = Path('results/dft/orr_confirmatory_slate.json')
SLATE_SIZE = 50


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def _rank(rows: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    return sorted(rows, key=lambda row: (row['predicted_overpotential_V'], row['candidate_id']))


def _take(
    candidates: Iterable[dict[str, Any]], count: int, selected: dict[str, dict[str, Any]], role: str
) -> None:
    for row in candidates:
        if len([item for item in selected.values() if item['selection_role'] == role]) >= count:
            return
        identity = str(row['candidate_id'])
        if identity not in selected:
            selected[identity] = {**row, 'selection_role': role}


def build_slate(output: Path = OUTPUT) -> Path:
    """Verify the completed campaign and write its diverse ORR DFT handoff.

    Args:
        output: Immutable JSON artifact destination.

    Returns:
        Path to the checksum-locked DFT slate.
    """
    final_path = ROOT / 'confirmatory_analysis.json'
    final = json.loads(final_path.read_text())
    if final.get('status') != 'complete' or final.get('acceptance_passed') is not True:
        raise RuntimeError('a passing completed confirmatory campaign is required')
    batch_ids = [str(value) for value in final['batch_ids']]
    if len(batch_ids) != int(final['planned_batches']):
        raise RuntimeError('confirmatory batch count does not match the locked plan')

    source_hashes = {str(final_path): _sha256(final_path)}
    rows: list[dict[str, Any]] = []
    for batch in batch_ids:
        folder = ROOT / 'batches' / batch
        manifest_path = folder / 'manifest.json'
        analysis_path = folder / 'analysis.json'
        outcome_path = Path('results/fuel_cell/prospective') / f'{batch}_orr.csv'
        manifest = json.loads(manifest_path.read_text())
        analysis = json.loads(analysis_path.read_text())
        if analysis['manifest_sha256'] != _sha256(manifest_path):
            raise RuntimeError(f'manifest checksum mismatch: {batch}')
        if analysis['outcome_sha256']['fuel_cell_orr'] != _sha256(outcome_path):
            raise RuntimeError(f'outcome checksum mismatch: {batch}')
        records = [item for item in manifest['records'] if item['application'] == 'fuel_cell_orr']
        if len(records) != 1:
            raise RuntimeError(f'exactly one ORR record is required: {batch}')
        record = records[0]
        memberships = {
            name: set(identities) for name, identities in record['policy_selected_ids'].items()
        }
        source_hashes[str(manifest_path)] = _sha256(manifest_path)
        source_hashes[str(analysis_path)] = _sha256(analysis_path)
        source_hashes[str(outcome_path)] = _sha256(outcome_path)
        with outcome_path.open(newline='') as stream:
            for raw in csv.DictReader(stream):
                if raw['valid'].lower() != 'true':
                    continue
                metric = float(raw['orr_overpotential_V'])
                if not math.isfinite(metric):
                    continue
                genome = ast.literal_eval(raw['genome'])
                identity = candidate_id(genome)
                if identity != raw['candidate_id']:
                    raise RuntimeError(f'candidate identity mismatch: {batch}')
                policies = sorted(name for name, ids in memberships.items() if identity in ids)
                rows.append({
                    'batch_id': batch,
                    'candidate_id': identity,
                    'genome': repr(genome),
                    'material_class': str(genome[0]),
                    'predicted_overpotential_V': metric,
                    'dG_OH_eV': float(raw['dG_OH_eV']),
                    'dG_O_eV': float(raw['dG_O_eV']),
                    'dG_OOH_eV': float(raw['dG_OOH_eV']),
                    'model_confidence': float(raw['model_confidence']),
                    'policy_membership': policies,
                })

    selected: dict[str, dict[str, Any]] = {}
    by_class: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        by_class[str(row['material_class'])].append(row)
    for material_class in sorted(by_class):
        champion = _rank(by_class[material_class])[0]
        selected[str(champion['candidate_id'])] = {
            **champion, 'selection_role': 'class_champion'
        }
    _take(_rank(row for row in rows if 'catalyst' in row['policy_membership']), 20, selected, 'guided_top')
    _take(_rank(row for row in rows if row['policy_membership'] == ['uncertainty']), 6, selected, 'uncertainty_control')
    _take(_rank(row for row in rows if row['policy_membership'] == ['validity']), 6, selected, 'validity_control')
    controls = sorted(
        (row for row in rows if not row['policy_membership']),
        key=lambda row: hashlib.sha256(str(row['candidate_id']).encode()).hexdigest(),
    )
    _take(controls, 6, selected, 'coverage_control')
    if len(selected) != SLATE_SIZE:
        raise RuntimeError(f'DFT slate has {len(selected)} candidates, expected {SLATE_SIZE}')

    records = []
    for rank, row in enumerate(selected.values(), 1):
        records.append({'rank': rank, **row})
    body = {
        'application': 'fuel_cell_orr',
        'screening_evidence': 'esen_gnn_not_dft',
        'validation_goal': 'measure GNN-to-DFT error, rank retention, and policy enrichment',
        'selection_counts': dict(sorted(
            (role, sum(item['selection_role'] == role for item in records))
            for role in {item['selection_role'] for item in records}
        )),
        'source_sha256': source_hashes,
        'records': records,
    }
    digest = hashlib.sha256(
        json.dumps(body, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    _atomic_json(output, {'schema_version': 1, **body, 'slate_sha256': digest})
    return output


def verify_slate(path: Path = OUTPUT) -> bool:
    """Return whether a slate's content digest and all source hashes match.

    Args:
        path: Slate artifact to verify.

    Returns:
        True only when the content digest and every evidence source match.
    """
    document = json.loads(path.read_text())
    digest = document.pop('slate_sha256')
    document.pop('schema_version')
    actual = hashlib.sha256(
        json.dumps(document, sort_keys=True, separators=(',', ':')).encode()
    ).hexdigest()
    return digest == actual and all(
        Path(source).is_file() and _sha256(Path(source)) == expected
        for source, expected in document['source_sha256'].items()
    )


if __name__ == '__main__':
    artifact = build_slate()
    if not verify_slate(artifact):
        raise SystemExit('generated ORR DFT slate failed verification')
    print(artifact)
