"""Independent branch-and-bound catalyst discovery stage."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Optional

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
    estimate_space: Callable
    build_config: Callable
    run_search: Callable
    annotate_evidence: Callable
    select_reactor: Callable
    select_validation: Callable
    select_admissible: Optional[Callable] = None


def default_discovery_services() -> DiscoveryServices:
    """Construct production discovery-stage dependencies.

    Returns:
        A `DiscoveryServices` containing the default discovery services result.
    """
    from pipeline.search.scope import scope_pyrolysis_pool
    from pipeline.search.design_space import estimate_design_space_size
    from pipeline.screening.genetic_optimizer import (
        BranchDiscoveryConfig, run_branch_discovery)
    from pipeline.screening.stage_selection import (
        annotate_evidence, select_for_reactor, select_for_validation)
    return DiscoveryServices(
        estimate_design_space_size, BranchDiscoveryConfig,
        run_branch_discovery, annotate_evidence,
        select_for_reactor, select_for_validation,
        select_admissible=scope_pyrolysis_pool)


def run_discovery_stage(*, initial_samples: int, leaf_size: int,
                        max_leaves: int | None, top_k_reactor: int,
                        top_k_dft: int,
                        services: DiscoveryServices | None = None) -> StageOutcome:
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
        branch_leaf_size=leaf_size, branch_max_leaves=max_leaves,
        expected_space_size=sizes['TOTAL'])
    pareto, database = services.run_search(config)
    valid = database[database['valid'] == True].copy()
    evidence = services.annotate_evidence(database, 'E_act')
    # Admissibility first (ADR 0001), then both slates from the same pool.
    # Validity is applied per route by the selectors (the validation route
    # may rescue invalid rows), so the pool is drawn from the full table.
    if services.select_admissible is not None:
        pool, admissibility = services.select_admissible(database)
    else:
        pool, admissibility = database, {'filter': None}
    reactor = services.select_reactor(
        pool, top_k_reactor, 'E_act', min_per_class=1)
    validation = services.select_validation(
        pool, top_k_dft, 'E_act', min_per_class=1)
    state = {
        'pareto_size': len(pareto), 'total_evaluated': len(database),
        'valid_count': len(valid), 'top_catalysts_count': len(reactor),
        'dft_resolution_count': len(validation),
        'admissibility': admissibility,
        'candidate_dispositions': evidence[
            'candidate_disposition'].value_counts().to_dict(),
    }
    if len(valid) > 0 and 'E_act' in valid.columns:
        state['best_E_act'] = float(valid['E_act'].min())
        state['best_coking'] = float(valid['coking_index'].max())
    return StageOutcome(state=state, products={
        'design_space_sizes': sizes, 'pareto_genomes': pareto,
        'screening_database': database, 'top_catalysts': reactor,
        'dft_candidates': validation})
