"""Candidate-to-pathway stage shared by interactive and production campaigns."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Mapping, Sequence

from pipeline.utils import repo_relative
from pipeline.reactors.modes import DEFAULT_MODE, resolve_pathway_mode
from pipeline.reactors.eligibility import rankable_results, usable_results


@dataclass(frozen=True)
class ReactorStageServices:
    """Injectable reactor-stage boundaries for tests and alternate workflows."""

    resolve_mode: Callable
    build_kinetics: Callable
    write_mechanism: Callable
    run_sweep: Callable


def default_reactor_services() -> ReactorStageServices:
    """Load production dependencies lazily so importing the stage stays cheap.

    Returns:
        Computed `ReactorStageServices` result.
    """
    from pipeline.reactors.mechanisms import (
        CandidateKinetics, write_full_mechanism)
    from pipeline.reactors.models import run_reactor_sweep
    return ReactorStageServices(
        resolve_mode=resolve_pathway_mode,
        build_kinetics=CandidateKinetics.from_screening_row,
        write_mechanism=write_full_mechanism,
        run_sweep=run_reactor_sweep)


def summarize_reactor_sweep(sweep: Sequence[Mapping]) -> dict:
    """Summarize condition-level evidence without executing reactor software.

    Args:
        sweep: Mapping supplying sweep.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    by_status = {}
    for result in sweep:
        by_status.setdefault(result.get('status', 'complete'), []).append(result)
    completed = by_status.get('complete', [])
    failed = by_status.get('failed', [])
    pending = by_status.get('validation_required', [])
    not_applicable = by_status.get('not_applicable', [])
    usable = usable_results(completed)
    rankable = rankable_results(usable)
    best = max(rankable, key=lambda result: result['CH4_conversion']) \
        if rankable else {}
    return {
        'best_condition': best,
        'sweep_status': (
            'complete' if len(completed) == len(sweep) else
            'partial' if completed else
            'validation_required' if pending else
            'not_applicable' if not_applicable else 'failed'),
        'completed_conditions': len(completed),
        'usable_conditions': len(usable),
        'failed_conditions': len(failed),
        'pending_conditions': len(pending),
        'not_applicable_conditions': len(not_applicable),
        'unresolved_conditions': len(sweep) - len(completed),
        'can_exclude_candidate': bool(completed) and not failed and all(
            result.get('can_exclude_candidate', False) for result in completed),
    }


def simulate_candidate(row: Mapping, catalyst_name: str,
                       temperatures: Sequence[float],
                       reactor_types: Sequence[str] | None = None,
                       forbid_mock: bool = True,
                       kinetics_validation: Mapping | None = None,
                       pathway_mode: str = DEFAULT_MODE,
                       multiphysics_results_dir: str | None = None,
                       services: ReactorStageServices | None = None) -> dict:
    """Build only the mechanism appropriate to the selected pathway and run it.

    Args:
        row: Mapping supplying row.
        catalyst_name: Catalyst name used by this operation.
        temperatures: Ordered values supplying temperatures.
        reactor_types: Ordered values supplying reactor types.
        forbid_mock: Whether to enable forbid mock.
        kinetics_validation: Mapping supplying kinetics validation.
        pathway_mode: Configured methane-conversion pathway.
        multiphysics_results_dir: Directory used for multiphysics results dir.
        services: Services used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    services = services or default_reactor_services()
    mode = services.resolve_mode(pathway_mode)
    candidate_id = str(row.get('candidate_id', catalyst_name))
    mechanism = None
    if not mode.requires_specialized_validation:
        kinetics = services.build_kinetics(
            row, candidate_id=candidate_id, validation=kinetics_validation)
        mechanism = services.write_mechanism(catalyst_name, kinetics=kinetics)
    barrier = float(row.get('E_act'))
    try:
        dE_H = float(row.get('dE_H'))
    except (TypeError, ValueError):
        dE_H = 0.0
    # Turquoise carbon policy for every production sweep: no CO2, circulating
    # fluidized solids, bounded mechanical regen. See docs/backlog/B2.
    carbon_policy = {
        'co2_permitted': False,
        'fluidized_mode': 'circulating',
        'max_regen_cycles': 3,
        'regen_mechanism': 'mechanical',
    }
    sweep = services.run_sweep(
        catalyst_name, str(mechanism or ''), temperatures=list(temperatures),
        reactor_types=None if reactor_types is None else list(reactor_types),
        catalyst_E_act_eV=barrier, pathway_mode=pathway_mode,
        material_class=row.get('material_class'), candidate_id=candidate_id,
        multiphysics_results_dir=multiphysics_results_dir,
        catalyst_dE_H_eV=dE_H, reactor_config_kwargs=carbon_policy)
    if forbid_mock and any(result.get('mock') for result in sweep):
        raise RuntimeError('mock reactor output is forbidden in production')
    summary = summarize_reactor_sweep(sweep)
    return {
        'catalyst': catalyst_name,
        'candidate_id': candidate_id,
        'material_class': row.get('material_class'),
        'pathway_mode': pathway_mode,
        'E_act': barrier,
        'mechanism_file': repo_relative(mechanism) if mechanism is not None else None,
        'sweep': sweep,
        **summary,
    }
