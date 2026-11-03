"""Independent coordination of a batch of candidate reactor simulations."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
import logging
from typing import Any, cast

from pipeline.data_models.reactors import ReactorResult, SolidsScorecard
from pipeline.data_models.stages import ReactorBatchProducts, ReactorBatchState
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
        load_reference_candidate: Optional loader for a missing named solids
            judge. Its row is added only to this batch, not the discovery slate.
    """

    prepare_gas_mechanism: Callable[..., Any]
    simulate_candidate: Callable[..., dict[str, Any]]
    write_mock_mechanism: Callable[..., Any]
    run_mock_sweep: Callable[..., list[ReactorResult]]
    check_equilibrium: Callable[..., dict[str, Any]] | None = None
    build_scorecard: Callable[..., SolidsScorecard] | None = None
    persist_scorecard: Callable[[SolidsScorecard], Any] | None = None
    load_reference_candidate: Callable[[str], Any] | None = None


def default_reactor_batch_services() -> ReactorBatchServices:
    """Construct production reactor-batch dependencies.

    Returns:
        A `ReactorBatchServices` containing the default reactor batch services result.
    """
    from pipeline.reactors.equilibrium import (
        run_equilibrium_sweep,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.reactors.mechanisms import (
        write_full_mechanism,  # pyright: ignore[reportUnknownVariableType]
        write_gri30_subset,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.reactors.models import (
        run_reactor_sweep,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.reactors.scorecard import build_solids_scorecard
    from pipeline.stages.candidate_io import load_reactor_reference
    from pipeline.stages.reactor import simulate_candidate
    from pipeline.utils import save_json  # pyright: ignore[reportUnknownVariableType]

    return ReactorBatchServices(
        cast(Callable[..., Any], write_gri30_subset),
        cast(Callable[..., dict[str, Any]], simulate_candidate),
        cast(Callable[..., Any], write_full_mechanism),
        cast(Callable[..., list[ReactorResult]], run_reactor_sweep),
        check_equilibrium=cast(Callable[..., dict[str, Any]], run_equilibrium_sweep),
        build_scorecard=build_solids_scorecard,
        persist_scorecard=lambda scorecard: save_json(
            dict(scorecard), 'phase2_solids_scorecard.json', subdir='reactor'
        ),
        load_reference_candidate=load_reactor_reference,
    )


def run_reactor_batch_stage(
    candidates: Any,
    *,
    temperatures: Sequence[float],
    reactor_types: Sequence[str],
    pathway_mode: str,
    multiphysics_results_dir: str,
    allow_mock_inputs: bool,
    judge_catalyst: str | None = None,
    headline_t_min: float = 1200.0,
    headline_t_max: float | None = None,
    services: ReactorBatchServices | None = None,
) -> StageOutcome[ReactorBatchState, ReactorBatchProducts]:
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
        headline_t_max: Highest T (K) counted as the headline band, or
            None for no upper bound (legacy 4-point sweep uses 1200+).
        services: Services used by this operation.

    Returns:
        Computed `StageOutcome` result.
    """
    services = services or default_reactor_batch_services()
    services.prepare_gas_mechanism()
    equilibrium: dict[str, Any] | None = None
    if services.check_equilibrium is not None:
        equilibrium = services.check_equilibrium()
        if not equilibrium.get('within_tolerance', False):
            logger.warning(
                'Equilibrium check outside tolerance '
                f"(worst_abs_error={equilibrium.get('worst_abs_error')})"
            )
    results: list[ReactorResult] = []
    if candidates is not None:
        named_rows: list[tuple[str, Any]] = []
        for index, row in candidates.iterrows():
            name = (
                judge_catalyst
                if judge_catalyst and str(row.get('candidate_id')) == judge_catalyst
                else f'cat_{index}'
            )
            named_rows.append((name, row))
        if (
            judge_catalyst
            and services.load_reference_candidate is not None
            and {'PFR', 'PBR', 'Fluidized'}.intersection(reactor_types)
            and not any(name == judge_catalyst for name, _ in named_rows)
        ):
            reference = services.load_reference_candidate(judge_catalyst)
            if reference is not None:
                named_rows.append((judge_catalyst, reference))
        for name, row in named_rows:
            candidate = services.simulate_candidate(
                row,
                name,
                temperatures,
                reactor_types,
                forbid_mock=not allow_mock_inputs,
                pathway_mode=pathway_mode,
                multiphysics_results_dir=multiphysics_results_dir,
            )
            results.extend(cast(list[ReactorResult], candidate['sweep']))
    else:
        if not allow_mock_inputs:
            raise RuntimeError('mock catalyst fallback is disabled')
        material_class = 'MoltenMetal' if pathway_mode == 'mmbcr' else 'SolidCatalyst'
        # Test-only C-H activation barrier fixtures (eV): molten Ni-Bi (0.85), FeC (0.65), Cu-Sn (1.10)
        for name, barrier in (
            ('NiBi_10', 0.85),
            ('FeC_supported', 0.65),
            ('CuSn_20', 1.1),
        ):
            mechanism = services.write_mock_mechanism(name, E_act_CH4=barrier)
            results.extend(
                services.run_mock_sweep(
                    name,
                    str(mechanism),
                    temperatures=list(temperatures),
                    reactor_types=list(reactor_types),
                    pathway_mode=pathway_mode,
                    material_class=material_class,
                    multiphysics_results_dir=multiphysics_results_dir,
                )
            )
    state: ReactorBatchState = {'n_simulations': len(results)}
    if equilibrium is not None:
        state['equilibrium_check'] = {
            'within_tolerance': bool(equilibrium.get('within_tolerance')),
            'worst_abs_error': float(equilibrium.get('worst_abs_error') or 0.0),
        }
    if services.build_scorecard is not None:
        from pipeline.economics.tea import (
            estimate_hydrogen_cost,  # pyright: ignore[reportUnknownVariableType]
        )
        from pipeline.validation.viability import evaluate_turquoise

        # Solids are judged on single-pass X of a named (or best non-H-parked)
        # solids catalyst. MMBCR X_eq is reported separately and never ranks.
        scorecard = services.build_scorecard(
            results,
            judge_catalyst=judge_catalyst,
            headline_t_min=headline_t_min,
            headline_t_max=headline_t_max,
        )
        if services.persist_scorecard is not None:
            services.persist_scorecard(scorecard)
        judge_name = (
            scorecard.get('judge_catalyst')
            or scorecard.get('headline_catalyst')
            or 'best_non_h_parked'
        )
        viability_counts: dict[str, int] = {'pass': 0, 'fail': 0, 'unknown': 0}
        for row in results:
            v_status = str(evaluate_turquoise(row).get('status', 'unknown'))
            viability_counts[v_status] = viability_counts.get(v_status, 0) + 1
        headline_conv = scorecard.get('headline_solids_conversion')
        champ_conv = scorecard.get('candidate_champion_conversion')
        state.update(
            {
                'solids_scorecard': scorecard,
                'best_conversion': headline_conv,
                'best_conversion_scope': f'solids_single_pass_judge_{judge_name}',
                'candidate_champion_conversion': champ_conv,
                'candidate_champion_catalyst': scorecard.get(
                    'candidate_champion_catalyst'
                ),
                'mmbcr_max_conversion': scorecard.get('mmbcr_max_conversion'),
                'industrial_viability': viability_counts,
            }
        )
        tea_conv = champ_conv if champ_conv is not None else headline_conv
        if tea_conv is not None and float(tea_conv) > 0.0:
            tea_res = cast(
                dict[str, Any], estimate_hydrogen_cost(float(tea_conv), 'base')
            )
            state['best_h2_cost_usd_kg'] = float(tea_res['h2_cost_usd_kg'])
    elif results:
        from pipeline.reactors.eligibility import rankable_results

        rankable = rankable_results(results)
        state['best_conversion'] = max(
            (float(row.get('CH4_conversion') or 0.0) for row in rankable), default=None
        )
    products: ReactorBatchProducts = {'reactor_results': results}
    return StageOutcome(state=state, products=products)
