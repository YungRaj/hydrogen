"""Candidate-to-pathway stage shared by interactive and production campaigns."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Mapping, Sequence, cast

from pipeline.utils import repo_relative
from pipeline.reactors.modes import DEFAULT_MODE, PathwayMode, resolve_pathway_mode
from pipeline.reactors.eligibility import rankable_results, usable_results
from pipeline.data_models.reactors import (
    CandidateReactorResult,
    CarbonPolicyConfig,
    PyrolysisValidationArtifact,
    ReactorResult,
    ReactorStatus,
    ReactorSweepSummary,
)


@dataclass(frozen=True)
class ReactorStageServices:
    """Injectable reactor-stage boundaries for tests and alternate workflows."""

    resolve_mode: Callable[[str], PathwayMode]
    build_kinetics: Callable[..., object]
    write_mechanism: Callable[..., Path | str | None]
    run_sweep: Callable[..., list[ReactorResult]]


def default_reactor_services() -> ReactorStageServices:
    """Load production dependencies lazily so importing the stage stays cheap.

    Returns:
        Computed `ReactorStageServices` result.
    """
    from pipeline.reactors.mechanisms import CandidateKinetics, write_full_mechanism
    from pipeline.reactors.models import run_reactor_sweep

    return ReactorStageServices(
        resolve_mode=resolve_pathway_mode,
        build_kinetics=CandidateKinetics.from_screening_row,
        write_mechanism=write_full_mechanism,
        run_sweep=run_reactor_sweep,
    )


def summarize_reactor_sweep(sweep: Sequence[ReactorResult]) -> ReactorSweepSummary:
    """Summarize condition-level evidence without executing reactor software.

    Args:
        sweep: Mapping supplying sweep.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    by_status: dict[str, list[ReactorResult]] = {}
    for result in sweep:
        by_status.setdefault(str(result.get('status', 'complete')), []).append(result)
    completed = by_status.get('complete', [])
    failed = by_status.get('failed', [])
    pending = by_status.get('validation_required', [])
    not_applicable = by_status.get('not_applicable', [])
    usable = usable_results(completed)
    rankable = rankable_results(usable)
    best: ReactorResult = (
        max(rankable, key=lambda result: float(result.get('CH4_conversion', 0.0)))
        if rankable
        else {}
    )
    sweep_status: ReactorStatus = (
        'complete'
        if len(completed) == len(sweep)
        else (
            'partial'
            if completed
            else (
                'validation_required'
                if pending
                else 'not_applicable' if not_applicable else 'failed'
            )
        )
    )
    return {
        'best_condition': best,
        'sweep_status': sweep_status,
        'completed_conditions': len(completed),
        'usable_conditions': len(usable),
        'failed_conditions': len(failed),
        'pending_conditions': len(pending),
        'not_applicable_conditions': len(not_applicable),
        'unresolved_conditions': len(sweep) - len(completed),
        'can_exclude_candidate': bool(completed)
        and not failed
        and all(bool(result.get('can_exclude_candidate', False)) for result in completed),
    }


def simulate_candidate(
    row: Mapping[str, object],
    catalyst_name: str,
    temperatures: Sequence[float],
    reactor_types: Sequence[str] | None = None,
    forbid_mock: bool = True,
    kinetics_validation: PyrolysisValidationArtifact | Mapping[str, object] | None = None,
    pathway_mode: str = DEFAULT_MODE,
    multiphysics_results_dir: str | None = None,
    services: ReactorStageServices | None = None,
) -> CandidateReactorResult:
    """Build the candidate mechanism and run the selected pathway's reactors.

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
    services.resolve_mode(pathway_mode)
    candidate_id = str(row.get('candidate_id', catalyst_name))
    # Every mode gets the candidate's Cantera mechanism: thermal reactors run
    # it directly, and NTEC/electrochemical copy it into their (calibrated
    # artifact or uncalibrated screening) MOOSE <-> Cantera interface case.
    kinetics = services.build_kinetics(
        row, candidate_id=candidate_id, validation=kinetics_validation
    )
    mechanism = services.write_mechanism(catalyst_name, kinetics=kinetics)
    raw_barrier = row.get('E_act')
    if not isinstance(raw_barrier, (int, float)):
        barrier = float(str(raw_barrier))
    else:
        barrier = float(raw_barrier)
    raw_deh = row.get('dE_H')
    try:
        dE_H = float(raw_deh) if isinstance(raw_deh, (int, float, str)) else 0.0
    except (TypeError, ValueError):
        dE_H = 0.0
    # Turquoise carbon policy for every production sweep: no CO2, circulating
    # fluidized solids, bounded mechanical regen. See docs/backlog/B2.
    carbon_policy: CarbonPolicyConfig = {
        'co2_permitted': False,
        'fluidized_mode': 'circulating',
        'max_regen_cycles': 3,
        'regen_mechanism': 'mechanical',
    }
    raw_class = row.get('material_class')
    material_class = str(raw_class) if raw_class is not None else None
    sweep = services.run_sweep(
        catalyst_name,
        str(mechanism or ''),
        temperatures=list(temperatures),
        reactor_types=None if reactor_types is None else list(reactor_types),
        catalyst_E_act_eV=barrier,
        pathway_mode=pathway_mode,
        material_class=material_class,
        candidate_id=candidate_id,
        multiphysics_results_dir=multiphysics_results_dir,
        catalyst_dE_H_eV=dE_H,
        reactor_config_kwargs=carbon_policy,
    )
    if forbid_mock and any(result.get('mock') for result in sweep):
        raise RuntimeError('mock reactor output is forbidden in production')
    summary = summarize_reactor_sweep(sweep)
    return cast(
        CandidateReactorResult,
        {
            'catalyst': catalyst_name,
            'candidate_id': candidate_id,
            'material_class': material_class,
            'pathway_mode': pathway_mode,
            'E_act': barrier,
            'mechanism_file': repo_relative(mechanism) if mechanism is not None else None,
            'sweep': sweep,
            **summary,
        },
    )

