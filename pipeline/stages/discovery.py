"""Independent branch-and-bound catalyst discovery stage."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from pipeline.stages.contracts import StageOutcome


@dataclass(frozen=True)
class DiscoveryServices:
    estimate_space: Callable
    build_config: Callable
    run_search: Callable
    annotate_evidence: Callable
    select_reactor: Callable
    select_validation: Callable


def default_discovery_services() -> DiscoveryServices:
    from pipeline.common.catalyst_spaces import estimate_design_space_size
    from pipeline.screening.genetic_optimizer import (
        BranchDiscoveryConfig, run_branch_discovery)
    from pipeline.screening.stage_selection import (
        annotate_evidence, select_for_reactor, select_for_validation)
    return DiscoveryServices(
        estimate_design_space_size, BranchDiscoveryConfig,
        run_branch_discovery, annotate_evidence,
        select_for_reactor, select_for_validation)


def run_discovery_stage(*, initial_samples: int, leaf_size: int,
                        max_leaves: int | None, top_k_reactor: int,
                        top_k_dft: int,
                        services: DiscoveryServices | None = None) -> StageOutcome:
    """Search, annotate, and route candidates through explicit dependencies."""
    services = services or default_discovery_services()
    sizes = services.estimate_space()
    config = services.build_config(
        initial_fairchem_samples=initial_samples,
        branch_leaf_size=leaf_size, branch_max_leaves=max_leaves,
        expected_space_size=sizes['TOTAL'])
    pareto, database = services.run_search(config)
    valid = database[database['valid'] == True].copy()
    evidence = services.annotate_evidence(database, 'E_act')
    reactor = services.select_reactor(
        database, top_k_reactor, 'E_act', min_per_class=1)
    validation = services.select_validation(
        database, top_k_dft, 'E_act', min_per_class=1)
    state = {
        'pareto_size': len(pareto), 'total_evaluated': len(database),
        'valid_count': len(valid), 'top_catalysts_count': len(reactor),
        'dft_resolution_count': len(validation),
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
