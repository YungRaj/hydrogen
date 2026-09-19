"""Independent coordination of a batch of candidate reactor simulations."""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Callable, Optional, Sequence

from pipeline.stages.contracts import StageOutcome

logger = logging.getLogger('reactor_batch')


@dataclass(frozen=True)
class ReactorBatchServices:
    """Bundle replaceable reactor-batch operations.

    Attributes:
        prepare_gas_mechanism: Configured prepare gas mechanism value.
        simulate_candidate: Configured simulate candidate value.
        write_mock_mechanism: Configured write mock mechanism value.
        run_mock_sweep: Configured run mock sweep value.
        check_equilibrium: Optional gas-phase X_eq tolerance check run once
            before the batch; warn-only.
        build_scorecard: Optional solids scorecard builder
            ``(results, judge_catalyst, headline_t_min) -> dict``. When
            present it replaces the raw ``max(CH4_conversion)`` headline,
            which would otherwise rank MMBCR X_eq as "best".
        persist_scorecard: Optional ``(scorecard) -> None`` writer.
    """
    prepare_gas_mechanism: Callable
    simulate_candidate: Callable
    write_mock_mechanism: Callable
    run_mock_sweep: Callable
    check_equilibrium: Optional[Callable] = None
    build_scorecard: Optional[Callable] = None
    persist_scorecard: Optional[Callable] = None


def default_reactor_batch_services() -> ReactorBatchServices:
    """Construct production reactor-batch dependencies.

    Returns:
        A `ReactorBatchServices` containing the default reactor batch services result.
    """
    from pipeline.common.utils import save_json
    from pipeline.process.equilibrium_check import run_equilibrium_sweep
    from pipeline.process.phase2_scorecard import build_solids_scorecard
    from pipeline.process.reactor_mechanisms import (
        write_full_mechanism, write_gri30_subset)
    from pipeline.process.reactor_models import run_reactor_sweep
    from pipeline.stages.reactor import simulate_candidate
    return ReactorBatchServices(
        write_gri30_subset, simulate_candidate,
        write_full_mechanism, run_reactor_sweep,
        check_equilibrium=run_equilibrium_sweep,
        build_scorecard=build_solids_scorecard,
        persist_scorecard=lambda scorecard: save_json(
            scorecard, 'phase2_solids_scorecard.json', subdir='reactor'))


def run_reactor_batch_stage(
        candidates, *, temperatures: Sequence[float],
        reactor_types: Sequence[str], pathway_mode: str,
        multiphysics_results_dir: str, allow_mock_inputs: bool,
        judge_catalyst: Optional[str] = None,
        headline_t_min: float = 1200.0,
        services: ReactorBatchServices | None = None) -> StageOutcome:
    """Run candidate sweeps; mock inputs require the existing explicit opt-in.

    Args:
        candidates: Candidate records to process.
        temperatures: Ordered values supplying temperatures.
        reactor_types: Ordered values supplying reactor types.
        pathway_mode: Configured methane-conversion pathway.
        multiphysics_results_dir: Directory used for multiphysics results dir.
        allow_mock_inputs: Whether to enable allow mock inputs.
        judge_catalyst: Named solids judge for the scorecard (campaign
            choice; None = best non-H-parked solids at the headline band).
        headline_t_min: Lowest T (K) counted as the headline band.
        services: Services used by this operation.

    Returns:
        Computed `StageOutcome` result.
    """
    services = services or default_reactor_batch_services()
    services.prepare_gas_mechanism()
    equilibrium = None
    if services.check_equilibrium is not None:
        equilibrium = services.check_equilibrium()
        if not equilibrium.get('within_tolerance', False):
            logger.warning(
                'Equilibrium check outside tolerance '
                f"(worst_abs_error={equilibrium.get('worst_abs_error')})")
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
    if equilibrium is not None:
        state['equilibrium_check'] = {
            'within_tolerance': equilibrium.get('within_tolerance'),
            'worst_abs_error': equilibrium.get('worst_abs_error'),
        }
    if services.build_scorecard is not None:
        # Solids are judged on single-pass X of a named (or best non-H-parked)
        # solids catalyst. MMBCR X_eq is reported separately and never ranks.
        scorecard = services.build_scorecard(
            results, judge_catalyst=judge_catalyst,
            headline_t_min=headline_t_min)
        if services.persist_scorecard is not None:
            services.persist_scorecard(scorecard)
        judge_name = (scorecard.get('judge_catalyst')
                      or scorecard.get('headline_catalyst')
                      or 'best_non_h_parked')
        state.update({
            'solids_scorecard': scorecard,
            'best_conversion': scorecard.get('headline_solids_conversion'),
            'best_conversion_scope': f'solids_single_pass_judge_{judge_name}',
            'mmbcr_max_conversion': scorecard.get('mmbcr_max_conversion'),
        })
    elif results:
        state['best_conversion'] = max(
            row.get('CH4_conversion', 0) for row in results)
    return StageOutcome(state=state, products={'reactor_results': results})
