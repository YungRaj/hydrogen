#!/usr/bin/env python3
"""Turquoise-hydrogen catalyst-guided branch discovery."""

from dataclasses import dataclass
import argparse
from typing import Optional

import pandas as pd

from pipeline.search.design_space import (
    extract_elements as _extract_elements_from_genome,
)
from pipeline.search.adaptive_validation import AdvancementCriterion
from pipeline.screening.branch_discovery import (
    BranchDiscoveryConfig as _BranchDiscoveryConfig,
    BranchDiscoveryWorkflow,
    PilotBudget,
    run_guided_branch_discovery,
)
from pipeline.utils import SCREENING_DIR, setup_logger

logger = setup_logger('genetic_optimizer', 'screening/genetic_optimizer.log')
ADVANCEMENT_BARRIER_EV = 0.8


@dataclass
class BranchDiscoveryConfig(_BranchDiscoveryConfig):
    """Configure turquoise-hydrogen catalyst-guided branch discovery."""

    exhaustive_db: str = str(SCREENING_DIR / 'indexed_scan.sqlite')


def run_branch_discovery(
    config: BranchDiscoveryConfig = BranchDiscoveryConfig(),
    existing_db: Optional[pd.DataFrame] = None,
) -> tuple[list[tuple], pd.DataFrame]:
    """Run the standard catalyst-guided turquoise-hydrogen branch search.

    Args:
        config: Branch-search resource and coverage settings.
        existing_db: Optional completed evidence to use for calibration.

    Returns:
        Pareto champion genomes and annotated screening evidence.
    """
    from pipeline.screening.small_data_ranker import turquoise_tree_objectives
    from pipeline.screening.surface_screener import (
        SCREENING_PROTOCOL_ID,
        run_screening,
    )

    workflow = BranchDiscoveryWorkflow(
        application='turquoise_hydrogen',
        outcome='E_act',
        protocol_id=SCREENING_PROTOCOL_ID,
        calibration_db='branch_calibration.csv',
        refill_prefix='branch_calibration_refill',
        evidence_db='branch_ranker_evidence.csv',
        champions_db='branch_champions.csv',
        database_subdir='screening',
        certificate_path=(
            SCREENING_DIR / 'turquoise_hydrogen_coverage_certificate.json'
        ),
        validation_fidelity=ADVANCEMENT_BARRIER_EV,
        diagnostics_label='Pyrolysis',
        screener=run_screening,
        objective_function=turquoise_tree_objectives,
        advancement_criteria=(
            AdvancementCriterion('valid', 1.0, 'max'),
            AdvancementCriterion('pyrolysis_viable', 1.0, 'max'),
            AdvancementCriterion('E_act', ADVANCEMENT_BARRIER_EV, 'min'),
            AdvancementCriterion('model_confidence', 0.5, 'max'),
        ),
        portfolio_path=SCREENING_DIR / 'advancement_portfolio.json',
        logger=logger,
    )
    return run_guided_branch_discovery(config, workflow, existing_db)


def main() -> None:
    """Run a bounded file-backed turquoise throughput probe."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe-seconds', type=float, default=1800.0)
    parser.add_argument('--calibration-candidates', type=int, default=500)
    args = parser.parse_args()
    candidates = int(args.calibration_candidates)
    budget = PilotBudget(
        calibration_candidates=candidates,
        validation_candidates=min(20, candidates),
        portfolio_candidates=10,
        archive_candidates=100,
        branch_leaves=1,
        runtime_s=args.probe_seconds,
        scan_workers=1,
        workers_per_gpu=1,
    )
    config = BranchDiscoveryConfig(
        initial_fairchem_samples=candidates,
        fairchem_eval_top_k=budget.validation_candidates,
        experimental_slate_size=budget.portfolio_candidates,
        htvs_pool_size=budget.archive_candidates,
        branch_max_leaves=budget.branch_leaves,
        max_runtime_s=budget.runtime_s,
        scan_workers=budget.scan_workers,
        workers_per_gpu=budget.workers_per_gpu,
        branch_exploration_interval=1,
        pilot_budget=budget,
    )
    champions, evidence = run_branch_discovery(config)
    print({'champions': len(champions), 'evidence_rows': len(evidence)})


if __name__ == '__main__':
    main()
