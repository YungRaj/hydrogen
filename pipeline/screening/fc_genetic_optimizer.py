#!/usr/bin/env python3
"""Fuel-cell ORR objectives and catalyst-guided branch discovery."""

from dataclasses import dataclass
from typing import Optional

import pandas as pd

from pipeline.search.adaptive_validation import AdvancementCriterion

from pipeline.search.design_space import (
    extract_elements as _extract_elements_from_genome,
)
from pipeline.screening.branch_discovery import (
    BranchDiscoveryConfig as _BranchDiscoveryConfig,
    BranchDiscoveryWorkflow,
    run_guided_branch_discovery,
)
from pipeline.screening.fc_screener import FENTON_RISK
from pipeline.utils import DFT_DIR, FUEL_CELL_DIR, setup_logger

logger = setup_logger('fc_genetic_optimizer', 'fuel_cell/fc_genetic_optimizer.log')
# CHE Sabatier volcano advancement ceiling (0.40 V ~ Pt(111) benchmark overpotential).
ADVANCEMENT_OVERPOTENTIAL_V = 0.4


@dataclass
class FCBranchDiscoveryConfig(_BranchDiscoveryConfig):
    """Configure fuel-cell ORR catalyst-guided branch discovery."""

    exhaustive_db: str = str(FUEL_CELL_DIR / 'indexed_scan.sqlite')
    calibration_artifact: Optional[str] = str(
        DFT_DIR / 'active_learning_calibration.json'
    )


def _fenton_from_genome(genome: tuple) -> float:
    """Compute Fenton stability score from genome elements."""
    elements = _extract_elements_from_genome(genome)
    risk = sum(FENTON_RISK.get(element, 0) for element in elements)
    # 10 = completely inert to H2O2 -> HO* Fenton radical attack (risk=0, e.g. Pt/Pd);
    # each Fenton-active element subtracts its FENTON_RISK severity down to a 0 floor.
    return float(max(0, 10 - risk))


def _cost_from_genome(genome: tuple) -> float:
    """Compute the minimized abundance-cost objective for a genome."""
    from pipeline.utils import abundance_cost_penalty

    # Negate abundance_cost_penalty (which lies in [-2, 0]) so [0, 2] is minimized.
    return -abundance_cost_penalty(_extract_elements_from_genome(genome))


def run_fc_branch_discovery(
    config: FCBranchDiscoveryConfig,
    existing_db: Optional[pd.DataFrame] = None,
) -> tuple[list[tuple], pd.DataFrame]:
    """Run the standard catalyst-guided ORR branch search.

    Args:
        config: Branch-search resource and coverage settings.
        existing_db: Optional completed evidence to use for calibration.

    Returns:
        Pareto champion genomes and annotated screening evidence.
    """
    from pipeline.screening.fc_screener import (
        SCREENING_PROTOCOL_ID,
        run_orr_screening,
    )
    from pipeline.screening.small_data_ranker import orr_tree_objectives

    workflow = BranchDiscoveryWorkflow(
        application='fuel_cell_orr',
        outcome='orr_overpotential_V',
        protocol_id=SCREENING_PROTOCOL_ID,
        calibration_db='fc_branch_calibration.csv',
        refill_prefix='fc_branch_calibration_refill',
        evidence_db='fc_branch_ranker_evidence.csv',
        champions_db='fc_branch_champions.csv',
        database_subdir='fuel_cell',
        certificate_path=FUEL_CELL_DIR / 'coverage_certificate.json',
        validation_fidelity=ADVANCEMENT_OVERPOTENTIAL_V,
        diagnostics_label='ORR',
        screener=run_orr_screening,
        objective_function=orr_tree_objectives,
        advancement_criteria=(
            AdvancementCriterion('valid', 1.0, 'max'),
            AdvancementCriterion('fc_viable', 1.0, 'max'),
            AdvancementCriterion(
                'orr_overpotential_V', ADVANCEMENT_OVERPOTENTIAL_V, 'min'
            ),
            # fenton_stability >= 7.0 excludes unmitigated Fe-rich Fenton radical generators;
            # model_confidence >= 0.5 requires at least moderate OC20/OC25 OOD coverage.
            AdvancementCriterion('fenton_stability', 7.0, 'max'),
            AdvancementCriterion('model_confidence', 0.5, 'max'),
        ),
        portfolio_path=FUEL_CELL_DIR / 'advancement_portfolio.json',
        logger=logger,
    )
    return run_guided_branch_discovery(config, workflow, existing_db)
