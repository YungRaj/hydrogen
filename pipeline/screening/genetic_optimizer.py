#!/usr/bin/env python3
"""Turquoise-hydrogen catalyst-guided branch discovery."""

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from pipeline.search.design_space import (
    extract_elements as _extract_elements_from_genome,
)
from pipeline.search.discovery import fast_non_dominated_sort
from pipeline.search.adaptive_validation import AdvancementCriterion
from pipeline.screening.branch_discovery import (
    BranchDiscoveryConfig as _BranchDiscoveryConfig,
    BranchDiscoveryWorkflow,
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
        logger=logger,
    )
    return run_guided_branch_discovery(config, workflow, existing_db)
