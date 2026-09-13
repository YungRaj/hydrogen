"""Independent coordination of a batch of candidate reactor simulations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Sequence

from pipeline.stages.contracts import StageOutcome


@dataclass(frozen=True)
class ReactorBatchServices:
    prepare_gas_mechanism: Callable
    simulate_candidate: Callable
    write_mock_mechanism: Callable
    run_mock_sweep: Callable


def default_reactor_batch_services() -> ReactorBatchServices:
    from pipeline.process.reactor_mechanisms import (
        write_full_mechanism, write_gri30_subset)
    from pipeline.process.reactor_models import run_reactor_sweep
    from pipeline.stages.reactor import simulate_candidate
    return ReactorBatchServices(
        write_gri30_subset, simulate_candidate,
        write_full_mechanism, run_reactor_sweep)


def run_reactor_batch_stage(
        candidates, *, temperatures: Sequence[float],
        reactor_types: Sequence[str], pathway_mode: str,
        multiphysics_results_dir: str, allow_mock_inputs: bool,
        services: ReactorBatchServices | None = None) -> StageOutcome:
    """Run candidate sweeps; mock inputs require the existing explicit opt-in."""
    services = services or default_reactor_batch_services()
    services.prepare_gas_mechanism()
    results = []
    if candidates is not None:
        for index, row in candidates.iterrows():
            candidate = services.simulate_candidate(
                row, f'cat_{index}', temperatures, reactor_types,
                forbid_mock=not allow_mock_inputs, pathway_mode=pathway_mode,
                multiphysics_results_dir=multiphysics_results_dir)
            results.extend(candidate['sweep'])
    else:
        if not allow_mock_inputs:
            raise RuntimeError('mock catalyst fallback is disabled')
        material_class = ('MoltenMetal' if pathway_mode == 'mmbcr'
                          else 'SolidCatalyst')
        for name, barrier in (
                ('NiBi_10', 0.85), ('FeC_supported', 0.65), ('CuSn_20', 1.1)):
            mechanism = services.write_mock_mechanism(
                name, E_act_CH4=barrier)
            results.extend(services.run_mock_sweep(
                name, str(mechanism), temperatures=list(temperatures),
                reactor_types=list(reactor_types), pathway_mode=pathway_mode,
                material_class=material_class,
                multiphysics_results_dir=multiphysics_results_dir))
    state = {'n_simulations': len(results)}
    if results:
        state['best_conversion'] = max(
            row.get('CH4_conversion', 0) for row in results)
    return StageOutcome(state=state, products={'reactor_results': results})
