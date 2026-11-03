"""Independent branch-and-bound catalyst discovery stage."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from typing import Any, cast

from pipeline.data_models.discovery import AdmissibilitySummary
from pipeline.data_models.stages import DiscoveryProducts, DiscoveryState, TableLike
from pipeline.stages.contracts import StageOutcome


@dataclass(frozen=True)
class DiscoveryServices:
    """Bundle replaceable discovery-stage operations.

    Attributes:
        estimate_space: Configured estimate space value.
        build_config: Configured build config value.
        run_search: Configured run search value.
        annotate_evidence: Configured annotate evidence value.
        select_reactor: Configured select reactor value.
        select_validation: Configured select validation value.
        select_admissible: Optional application-scope filter
            ``(frame) -> (pool, note)`` applied before either slate is drawn
            (ADR 0001: encoded phase must be stable at the pyrolysis T
            band). Selection-only; the coverage denominator is untouched.
    """

    estimate_space: Callable[[], dict[str, int]]
    build_config: Callable[..., Any]
    run_search: Callable[[Any], tuple[Sequence[object], Any]]
    annotate_evidence: Callable[[Any, str], Any]
    select_reactor: Callable[..., TableLike]
    select_validation: Callable[..., TableLike]
    select_admissible: Callable[[Any], tuple[Any, Any]] | None = None


def default_discovery_services() -> DiscoveryServices:
    """Construct production discovery-stage dependencies.

    Returns:
        A `DiscoveryServices` containing the default discovery services result.
    """
    from pipeline.screening.genetic_optimizer import (
        BranchDiscoveryConfig,  # pyright: ignore[reportUnknownVariableType]
        run_branch_discovery,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.screening.stage_selection import (
        annotate_evidence,  # pyright: ignore[reportUnknownVariableType]
        select_for_pyrolysis_reactor,  # pyright: ignore[reportUnknownVariableType]
        select_for_validation,  # pyright: ignore[reportUnknownVariableType]
    )
    from pipeline.search.design_space import estimate_design_space_size
    from pipeline.search.scope import (
        scope_pyrolysis_pool,  # pyright: ignore[reportUnknownVariableType]
    )

    return DiscoveryServices(
        estimate_design_space_size,
        cast(Callable[..., Any], BranchDiscoveryConfig),
        cast(Callable[[Any], tuple[Sequence[object], Any]], run_branch_discovery),
        cast(Callable[[Any, str], Any], annotate_evidence),
        cast(Callable[..., TableLike], select_for_pyrolysis_reactor),
        cast(Callable[..., TableLike], select_for_validation),
        select_admissible=cast(
            Callable[[Any], tuple[Any, Any]], scope_pyrolysis_pool
        ),
    )


def run_discovery_stage(
    *,
    initial_samples: int,
    leaf_size: int,
    max_leaves: int | None,
    top_k_reactor: int,
    top_k_dft: int,
    services: DiscoveryServices | None = None,
) -> StageOutcome[DiscoveryState, DiscoveryProducts]:
    """Search, annotate, and route candidates through explicit dependencies.

    Args:
        initial_samples: Initial samples used by this operation.
        leaf_size: Number of leaf size to use.
        max_leaves: Bound controlling max leaves.
        top_k_reactor: Bound controlling top k reactor.
        top_k_dft: Bound controlling top k dft.
        services: Services used by this operation.

    Returns:
        Computed `StageOutcome` result.
    """
    services = services or default_discovery_services()
    sizes = services.estimate_space()
    config = services.build_config(
        initial_fairchem_samples=initial_samples,
        branch_leaf_size=leaf_size,
        branch_max_leaves=max_leaves,
        expected_space_size=sizes['TOTAL'],
    )
    pareto, database = services.run_search(config)
    from pipeline.screening.small_data_ranker import (
        screening_metric_eligibility,  # pyright: ignore[reportUnknownVariableType]
    )

    eligibility_columns = {'valid', 'E_act', 'E_act_censored', 'pyrolysis_viable'}
    if eligibility_columns.issubset(database.columns):
        eligible_metric = cast(
            Any, screening_metric_eligibility(database, 'turquoise_hydrogen')
        )
        valid = database.loc[eligible_metric].copy()
    else:
        # Injectable unit-stage tables predate the screening contract. The
        # production selector still fails closed when these gates are absent.
        valid = database[database['valid'] == True].copy()
    evidence = services.annotate_evidence(database, 'E_act')
    # Admissibility first (ADR 0001), then both slates from the same pool.
    # Validity is applied per route by the selectors (the validation route
    # may rescue invalid rows), so the pool is drawn from the full table.
    admissibility: AdmissibilitySummary
    if services.select_admissible is not None:
        pool, raw_admissibility = services.select_admissible(database)
        admissibility = cast(AdmissibilitySummary, raw_admissibility)
    else:
        pool = database
        admissibility = {'filter': None}
    # Ensure at least 1 representative per material class in reactor and DFT slates
    reactor = services.select_reactor(pool, top_k_reactor, 'E_act', min_per_class=1)
    validation = services.select_validation(pool, top_k_dft, 'E_act', min_per_class=1)
    raw_dispositions = evidence['candidate_disposition'].value_counts().to_dict()
    dispositions: dict[str, int] = {
        str(k): int(v) for k, v in cast(dict[Any, Any], raw_dispositions).items()
    }
    state: DiscoveryState = {
        'pareto_size': len(pareto),
        'total_evaluated': len(database),
        'valid_count': len(valid),
        'top_catalysts_count': len(reactor),
        'dft_resolution_count': len(validation),
        'admissibility': admissibility,
        'candidate_dispositions': dispositions,
    }
    if len(valid) > 0:
        state['best_E_act'] = float(valid['E_act'].min())
        if 'coking_index' in valid.columns:
            coking = valid['coking_index'].dropna()
            if len(coking) > 0:
                state['best_coking'] = float(coking.max())
    products: DiscoveryProducts = {
        'design_space_sizes': sizes,
        'pareto_genomes': pareto,
        'screening_database': cast(TableLike, database),
        'top_catalysts': reactor,
        'dft_candidates': validation,
    }
    return StageOutcome(state=state, products=products)
