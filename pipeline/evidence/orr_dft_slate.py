"""Build a checksum-locked ORR DFT slate from confirmatory evidence."""

from __future__ import annotations

import ast
from collections import defaultdict
from collections.abc import Iterable
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import cast

from pipeline.data_models.evidence import (
    ORRConfirmatoryDFTSummary,
    ORRPolicyEnrichmentSummary,
    ORRRoleSummary,
    ORRSlateSelectedRow,
    ORRValidatedSlateRecord,
)
from pipeline.io.hashing import file_sha256
from pipeline.io.json import write_atomic_json
from pipeline.search.discovery import candidate_id

ROOT = Path('results/prospective_search')
OUTPUT = Path('results/dft/orr_confirmatory_slate.json')
# 50-candidate DFT slate: 14 class champions + 20 guided_top (with 2 champion overlaps)
# + 6 uncertainty controls + 6 validity controls + 6 coverage controls = 50 unique genomes
SLATE_SIZE = 50


def _rank(rows: Iterable[ORRSlateSelectedRow]) -> list[ORRSlateSelectedRow]:
    return sorted(
        rows,
        key=lambda row: (row['predicted_overpotential_V'], row['candidate_id']),
    )


def _take(
    candidates: Iterable[ORRSlateSelectedRow],
    count: int,
    selected: dict[str, ORRSlateSelectedRow],
    role: str,
) -> None:
    for row in candidates:
        if (
            len(
                [
                    item
                    for item in selected.values()
                    if item.get('selection_role') == role
                ]
            )
            >= count
        ):
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

    source_hashes = {str(final_path): file_sha256(final_path)}
    rows: list[ORRSlateSelectedRow] = []
    for batch in batch_ids:
        folder = ROOT / 'batches' / batch
        manifest_path = folder / 'manifest.json'
        analysis_path = folder / 'analysis.json'
        outcome_path = Path('results/fuel_cell/prospective') / f'{batch}_orr.csv'
        manifest = json.loads(manifest_path.read_text())
        analysis = json.loads(analysis_path.read_text())
        if analysis['manifest_sha256'] != file_sha256(manifest_path):
            raise RuntimeError(f'manifest checksum mismatch: {batch}')
        if analysis['outcome_sha256']['fuel_cell_orr'] != file_sha256(outcome_path):
            raise RuntimeError(f'outcome checksum mismatch: {batch}')
        records = [item for item in manifest['records'] if item['application'] == 'fuel_cell_orr']
        if len(records) != 1:
            raise RuntimeError(f'exactly one ORR record is required: {batch}')
        record = records[0]
        memberships = {
            name: set(identities) for name, identities in record['policy_selected_ids'].items()
        }
        source_hashes[str(manifest_path)] = file_sha256(manifest_path)
        source_hashes[str(analysis_path)] = file_sha256(analysis_path)
        source_hashes[str(outcome_path)] = file_sha256(outcome_path)
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

    selected: dict[str, ORRSlateSelectedRow] = {}
    by_class: dict[str, list[ORRSlateSelectedRow]] = defaultdict(list)
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

    slate_records: list[ORRSlateSelectedRow] = []
    for rank, row in enumerate(selected.values(), 1):
        slate_records.append({'rank': rank, **row})
    body = {
        'application': 'fuel_cell_orr',
        'screening_evidence': 'esen_gnn_not_dft',
        'validation_goal': 'measure GNN-to-DFT error, rank retention, and policy enrichment',
        'selection_counts': dict(sorted(
            (
                role,
                sum(item.get('selection_role') == role for item in slate_records),
            )
            for role in {
                str(item.get('selection_role', '')) for item in slate_records
            }
        )),
        'source_sha256': source_hashes,
        'records': slate_records,
    }
    from pipeline.io.hashing import canonical_json_sha256

    digest = canonical_json_sha256(body)
    write_atomic_json(output, {'schema_version': 1, **body, 'slate_sha256': digest})
    return output


def verify_slate(path: Path = OUTPUT) -> bool:
    """Return whether a slate's content digest and all source hashes match.

    Args:
        path: Slate artifact to verify.

    Returns:
        True only when the content digest and every evidence source match.
    """
    from pipeline.io.hashing import canonical_json_sha256

    document = json.loads(path.read_text())
    digest = document.pop('slate_sha256')
    document.pop('schema_version')
    actual = canonical_json_sha256(document)
    return digest == actual and all(
        Path(source).is_file() and file_sha256(Path(source)) == expected
        for source, expected in document['source_sha256'].items()
    )


SUMMARY_OUTPUT = Path('results/dft/orr_confirmatory_dft_summary.json')
# DOE high-activity ORR cathode screening threshold (eta_ORR <= 0.35 V vs. RHE)
DOE_ORR_OVERPOTENTIAL_TARGET_V = 0.35


def summarize_orr_slate_dft_validation(
    slate_path: Path = OUTPUT,
    dft_dir: Path = Path('results/dft'),
    output_path: Path = SUMMARY_OUTPUT,
) -> ORRConfirmatoryDFTSummary:
    """Aggregate QE pw.x DFT validation, GC-DFT shifts, and policy enrichment for the locked ORR slate.

    Args:
        slate_path: Path to the checksum-locked 50-candidate ORR confirmatory slate.
        dft_dir: Directory containing per-candidate QE ``*_dft.json`` outputs.
        output_path: Destination path for the canonical summary artifact.

    Returns:
        JSON-serializable summary dictionary written to ``output_path``.
    """
    if not verify_slate(slate_path):
        raise RuntimeError(f'ORR DFT slate failed checksum verification: {slate_path}')
    from pipeline.io.hashing import canonical_json_sha256
    from pipeline.screening.fc_screener import (
        compute_grand_canonical_orr_corrections,
    )

    slate = json.loads(slate_path.read_text(encoding='utf-8'))
    validated_records: list[ORRValidatedSlateRecord] = []
    by_role: dict[str, list[ORRValidatedSlateRecord]] = defaultdict(list)

    for row in slate['records']:
        cid = str(row['candidate_id'])
        genome = ast.literal_eval(str(row['genome']))
        dft_path = dft_dir / f'orr_slate_{cid}_dft.json'
        if not dft_path.is_file():
            alt_path = dft_dir / f'campaign_{cid}_dft.json'
            if alt_path.is_file():
                dft_path = alt_path
        dft_rec: dict[str, object] = (
            json.loads(dft_path.read_text(encoding='utf-8'))
            if dft_path.is_file()
            else {}
        )
        dg_oh = float(row['dG_OH_eV'])
        dg_o = float(row['dG_O_eV'])
        dg_ooh = float(row['dG_OOH_eV'])
        gc_corr = compute_grand_canonical_orr_corrections(
            dg_oh,
            dg_o,
            dg_ooh,
            genome=genome,
        )
        acid_stable = bool(gc_corr['pourbaix_acid_stable'])
        eta_che = float(row['predicted_overpotential_V'])
        eta_gc = float(gc_corr['orr_overpotential_gc_V'])
        qe_conv = bool(dft_rec.get('converged', False))
        dft_e_ry = dft_rec.get('dft_energy_Ry')
        dft_e_ev = dft_rec.get('dft_energy_eV')
        max_f = dft_rec.get('max_force_Ry_bohr')
        enriched = cast(
            ORRValidatedSlateRecord,
            {
                **row,
                'orr_overpotential_gc_V': round(eta_gc, 6),
                'gc_field_shift_OH_eV': round(
                    float(gc_corr['dG_OH_gc_eV']) - dg_oh, 6
                ),
                'pourbaix_dissolution_onset_V': round(
                    float(gc_corr['pourbaix_dissolution_onset_V']), 4
                ),
                'pourbaix_dissolution_margin_V': round(
                    float(gc_corr['pourbaix_dissolution_margin_V']), 4
                ),
                'pourbaix_acid_stable': acid_stable,
                'qe_converged': qe_conv,
                'qe_dft_energy_Ry': (
                    round(float(dft_e_ry), 6)
                    if isinstance(dft_e_ry, (int, float))
                    else None
                ),
                'qe_dft_energy_eV': (
                    round(float(dft_e_ev), 4)
                    if isinstance(dft_e_ev, (int, float))
                    else None
                ),
                'qe_max_force_Ry_bohr': (
                    round(float(max_f), 6)
                    if isinstance(max_f, (int, float))
                    else None
                ),
                'doe_target_hit': eta_che <= DOE_ORR_OVERPOTENTIAL_TARGET_V,
                'doe_target_hit_gc': eta_gc <= DOE_ORR_OVERPOTENTIAL_TARGET_V,
            },
        )
        validated_records.append(enriched)
        by_role[str(row['selection_role'])].append(enriched)

    role_summary: dict[str, ORRRoleSummary] = {}
    for role, items in sorted(by_role.items()):
        etas = [float(x['predicted_overpotential_V']) for x in items]
        etas_gc = [float(x['orr_overpotential_gc_V']) for x in items]
        role_summary[role] = {
            'count': len(items),
            'n_qe_converged': sum(1 for x in items if x['qe_converged']),
            'n_pourbaix_acid_stable': sum(1 for x in items if x['pourbaix_acid_stable']),
            'mean_overpotential_V': round(sum(etas) / len(etas), 4),
            'min_overpotential_V': round(min(etas), 4),
            'mean_overpotential_gc_V': round(sum(etas_gc) / len(etas_gc), 4),
            'min_overpotential_gc_V': round(min(etas_gc), 4),
            'doe_target_hits_le_0_35V': sum(1 for x in items if x['doe_target_hit']),
            'doe_target_hits_gc_le_0_35V': sum(
                1 for x in items if x['doe_target_hit_gc']
            ),
        }

    guided_pool = [
        x
        for x in validated_records
        if x.get('selection_role') in ('class_champion', 'guided_top')
    ]
    control_pool = [
        x
        for x in validated_records
        if x.get('selection_role')
        in ('uncertainty_control', 'validity_control', 'coverage_control')
    ]
    guided_hit_rate = (
        sum(1 for x in guided_pool if x['doe_target_hit']) / len(guided_pool)
        if guided_pool
        else 0.0
    )
    control_hit_rate = (
        sum(1 for x in control_pool if x['doe_target_hit']) / len(control_pool)
        if control_pool
        else 0.0
    )
    n_conv = sum(1 for x in validated_records if x['qe_converged'])
    enrichment: ORRPolicyEnrichmentSummary = {
        'guided_and_champion_count': len(guided_pool),
        'control_count': len(control_pool),
        'guided_doe_hit_rate': round(guided_hit_rate, 4),
        'control_doe_hit_rate': round(control_hit_rate, 4),
        'enrichment_ratio': (
            round(guided_hit_rate / control_hit_rate, 4)
            if control_hit_rate > 0
            else None
        ),
    }
    body: dict[str, object] = {
        'application': 'fuel_cell_orr',
        'slate_sha256': str(slate['slate_sha256']),
        'n_slate_candidates': len(validated_records),
        'n_qe_converged': n_conv,
        'qe_convergence_rate': round(n_conv / max(len(validated_records), 1), 4),
        'n_pourbaix_acid_stable': sum(
            1 for x in validated_records if x['pourbaix_acid_stable']
        ),
        'role_summary': role_summary,
        'policy_enrichment': enrichment,
        'records': validated_records,
    }
    digest = canonical_json_sha256(body)
    payload: ORRConfirmatoryDFTSummary = {
        'schema_version': 1,
        'application': 'fuel_cell_orr',
        'slate_sha256': str(slate['slate_sha256']),
        'n_slate_candidates': len(validated_records),
        'n_qe_converged': n_conv,
        'qe_convergence_rate': round(n_conv / max(len(validated_records), 1), 4),
        'n_pourbaix_acid_stable': sum(
            1 for x in validated_records if x['pourbaix_acid_stable']
        ),
        'role_summary': role_summary,
        'policy_enrichment': enrichment,
        'records': validated_records,
        'summary_sha256': digest,
    }
    write_atomic_json(output_path, dict(payload))
    return payload


if __name__ == '__main__':
    artifact = build_slate()
    if not verify_slate(artifact):
        raise SystemExit('generated ORR DFT slate failed verification')
    print(artifact)
