"""Candidate-to-pathway stage shared by interactive and production campaigns."""

from __future__ import annotations

from typing import Mapping, Sequence

from pipeline.process.pathway_modes import DEFAULT_MODE, resolve_pathway_mode


def simulate_candidate(row: Mapping, catalyst_name: str,
                       temperatures: Sequence[float],
                       reactor_types: Sequence[str] | None = None,
                       forbid_mock: bool = True,
                       kinetics_validation: Mapping | None = None,
                       pathway_mode: str = DEFAULT_MODE,
                       multiphysics_results_dir: str | None = None) -> dict:
    """Build only the mechanism appropriate to the selected pathway and run it."""
    from pipeline.process.reactor_mechanisms import (
        CandidateKinetics, write_full_mechanism)
    from pipeline.process.reactor_models import run_reactor_sweep

    mode = resolve_pathway_mode(pathway_mode)
    candidate_id = str(row.get('candidate_id', catalyst_name))
    mechanism = None
    if not mode.requires_specialized_validation:
        kinetics = CandidateKinetics.from_screening_row(
            row, candidate_id=candidate_id, validation=kinetics_validation)
        mechanism = write_full_mechanism(catalyst_name, kinetics=kinetics)
    barrier = float(row.get('E_act'))
    sweep = run_reactor_sweep(
        catalyst_name, str(mechanism or ''), temperatures=list(temperatures),
        reactor_types=None if reactor_types is None else list(reactor_types),
        catalyst_E_act_eV=barrier, pathway_mode=pathway_mode,
        material_class=row.get('material_class'), candidate_id=candidate_id,
        multiphysics_results_dir=multiphysics_results_dir)
    if forbid_mock and any(result.get('mock') for result in sweep):
        raise RuntimeError('mock reactor output is forbidden in production')
    by_status = {}
    for result in sweep:
        by_status.setdefault(result.get('status', 'complete'), []).append(result)
    completed = by_status.get('complete', [])
    failed = by_status.get('failed', [])
    pending = by_status.get('validation_required', [])
    not_applicable = by_status.get('not_applicable', [])
    best = max(completed, key=lambda result: result.get('CH4_conversion', 0.0)) \
        if completed else {}
    return {
        'catalyst': catalyst_name,
        'candidate_id': candidate_id,
        'material_class': row.get('material_class'),
        'pathway_mode': pathway_mode,
        'E_act': barrier,
        'mechanism_file': str(mechanism) if mechanism is not None else None,
        'sweep': sweep,
        'best_condition': best,
        'sweep_status': (
            'complete' if len(completed) == len(sweep) else
            'partial' if completed else
            'validation_required' if pending else
            'not_applicable' if not_applicable else 'failed'),
        'completed_conditions': len(completed),
        'failed_conditions': len(failed),
        'pending_conditions': len(pending),
        'not_applicable_conditions': len(not_applicable),
        'unresolved_conditions': len(sweep) - len(completed),
        'can_exclude_candidate': bool(completed) and not failed and all(
            result.get('can_exclude_candidate', False) for result in completed),
    }
