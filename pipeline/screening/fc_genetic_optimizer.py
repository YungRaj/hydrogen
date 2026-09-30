#!/usr/bin/env python3
# Fuel-cell objectives and branch-search orchestration.
"""Fuel-cell ORR objectives and catalyst-guided branch discovery.

The Pareto helpers remain part of objective and archive selection. Candidate
generation is handled exclusively by catalyst-guided branch discovery.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Optional

from pipeline.utils import setup_logger, FUEL_CELL_DIR
from pipeline.search.discovery import (
    add_discovery_metadata,
    candidate_id,
)
logger = setup_logger('fc_genetic_optimizer', 'fuel_cell/fc_genetic_optimizer.log')



@dataclass
class FCBranchDiscoveryConfig:
    """Configure the ORR branch-search adapter.

    Attributes:
        initial_fairchem_samples: Configured initial fairchem samples value.
        fairchem_eval_top_k: Configured fairchem eval top k value.
        n_models: Configured n models value.
        htvs_pool_size: Configured htvs pool size value.
        device: Configured device value.
        exhaustive_batch_size: Configured exhaustive batch size value.
        exhaustive_db: Configured exhaustive db value.
        branch_leaf_size: Configured branch leaf size value.
        branch_probe_count: Configured branch probe count value.
        branch_max_leaves: Configured branch max leaves value.
        expected_space_size: Configured expected space size value.
        max_runtime_s: Configured max runtime s value.
        prior_art_db: Configured prior art db value.
        min_validation_per_class: Configured min validation per class value.
        min_resolved_leaves_per_class: Configured min resolved leaves per class value.
        branch_exploration_interval: Configured branch exploration interval value.
        refresh_pending_priorities: Configured refresh pending priorities value.
        scan_workers: Configured scan workers value.
    """

    initial_fairchem_samples: int = 500
    fairchem_eval_top_k: int = 500
    n_models: int = 3
    htvs_pool_size: int = 20000
    device: str = 'cuda'
    exhaustive_batch_size: int = 65536
    exhaustive_db: str = str(FUEL_CELL_DIR / 'indexed_scan.sqlite')
    branch_leaf_size: int = 1_000_000
    branch_probe_count: int = 9
    branch_max_leaves: Optional[int] = None
    expected_space_size: Optional[int] = None
    max_runtime_s: Optional[float] = None
    prior_art_db: Optional[str] = None
    min_validation_per_class: int = 1
    min_resolved_leaves_per_class: int = 1
    branch_exploration_interval: int = 4
    refresh_pending_priorities: int = 10_000
    scan_workers: int = 8


# ═══════════════════════════════════════════════════════════════════════════════
# ORR-SPECIFIC OBJECTIVES
# ═══════════════════════════════════════════════════════════════════════════════


def _fenton_from_genome(genome: tuple) -> float:
    """Compute Fenton stability score from genome elements."""
    FENTON_RISK = {'Fe': 3, 'Cu': 2, 'Co': 1, 'Mn': 1, 'Cr': 1, 'V': 1}
    elements = _extract_elements_from_genome(genome)
    fenton_risk = sum(FENTON_RISK.get(e, 0) for e in elements)
    return float(max(0, 10 - fenton_risk))




def _cost_from_genome(genome: tuple) -> float:
    """Compute cost penalty from genome elements.

    abundance_cost_penalty() returns [-2, 0] where 0 = abundant, -2 = rare.
    Since archive objectives are minimized, negate the penalty so that:
      abundant → 0 (good)    rare → +2 (bad, penalized)
    """
    from pipeline.utils import abundance_cost_penalty

    elements = _extract_elements_from_genome(genome)
    return -abundance_cost_penalty(elements)


def _extract_elements_from_genome(genome: tuple) -> List[str]:
    """Extract metallic elements from a genome for cost scoring."""
    mat_class = genome[0]
    elements = []
    if mat_class == 'MoltenMetal':
        elements.append(genome[1])
        if genome[2] != 'None':
            elements.append(genome[2])
    elif mat_class == 'SolidCatalyst':
        elements.append(genome[1])
        for d in genome[5]:
            elements.append(d)
    elif mat_class == 'SAC':
        elements.append(genome[1])
    elif mat_class == 'DAC':
        elements.extend([genome[1], genome[2]])
    elif mat_class in ('MOF', 'COF'):
        if genome[1] != 'None':
            elements.append(genome[1])
    elif mat_class == 'Perovskite':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None':
            elements.append(genome[3])
    elif mat_class == 'MetalHydride':
        elements.append(genome[1])
        if genome[3] != 'None':
            elements.append(genome[3])
    elif mat_class == 'MAXPhase':
        elements.extend([genome[1], genome[2]])
        if genome[5] != 'None':
            elements.append(genome[5])
    elif mat_class == 'HEA':
        elements.extend(list(genome[1]))
    elif mat_class == 'Spinel':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None':
            elements.append(genome[3])
    elif mat_class == 'MXene':
        elements.append(genome[1])
        if genome[5] != 'None':
            elements.append(genome[5])
    elif mat_class == 'SAA':
        elements.extend([genome[1], genome[2]])
    elif mat_class == 'MetalFreeCarbon':
        pass  # no metals — zero cost
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# PARETO ARCHIVE SELECTION
# ═══════════════════════════════════════════════════════════════════════════════


def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """Fast non-dominated sorting using vectorized Pareto-front extraction.

    Args:
        objectives: Objectives used by this operation.

    Returns:
        List of computed or validated records.
    """
    n = len(objectives)
    remaining_indices = np.arange(n)
    fronts = []

    while len(remaining_indices) > 0:
        sub_objs = objectives[remaining_indices]
        is_efficient = np.ones(len(sub_objs), dtype=bool)
        for i in range(len(sub_objs)):
            if is_efficient[i]:
                dominated = np.all(sub_objs[i] <= sub_objs, axis=1) & np.any(
                    sub_objs[i] < sub_objs, axis=1
                )
                is_efficient[dominated] = False

        front_sub_idx = np.where(is_efficient)[0]
        front_global_idx = remaining_indices[front_sub_idx].tolist()
        fronts.append(front_global_idx)

        remaining_indices = np.delete(remaining_indices, front_sub_idx)

    return fronts




# ═══════════════════════════════════════════════════════════════════════════════
# BRANCH DISCOVERY
# ═══════════════════════════════════════════════════════════════════════════════


def run_fc_branch_discovery(config: FCBranchDiscoveryConfig, existing_db=None):
    """Run the standard catalyst-guided ORR branch search.

    Args:
        config: Configuration controlling this operation.
        existing_db: Existing db used by this operation.

    Returns:
        Computed result described above.
    """
    from pipeline.search.indexed_space import deterministic_tree_probes
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
    from pipeline.search.exhaustive_search import load_archive_genomes
    from pipeline.screening.fc_screener import run_orr_screening, SCREENING_PROTOCOL_ID

    if existing_db is not None and len(existing_db) > 50:
        evidence = existing_db
    else:
        probes = deterministic_tree_probes(config.initial_fairchem_samples)
        evidence = run_orr_screening(
            probes, db_filename='fc_branch_calibration.csv', workers_per_gpu=2
        )
    from pipeline.screening.small_data_ranker import (
        MIN_TRAINING_ROWS,
        fit_tree_ranker,
        merge_compatible_evidence,
        orr_tree_objectives,
        valid_training_row_count,
    )
    from pipeline.utils import load_screening_db, save_screening_db

    prior_evidence = load_screening_db(
        'fc_branch_ranker_evidence.csv', subdir='fuel_cell'
    )
    evidence = merge_compatible_evidence(
        evidence, prior_evidence, SCREENING_PROTOCOL_ID
    )
    attempted = {str(value) for value in evidence.get('genome', [])}
    refill_limit = max(
        config.initial_fairchem_samples * 3,
        config.initial_fairchem_samples + MIN_TRAINING_ROWS,
    )
    probe_pool = deterministic_tree_probes(refill_limit)
    refill_round = 0
    while valid_training_row_count(evidence, 'fuel_cell_orr') < MIN_TRAINING_ROWS:
        refill = [genome for genome in probe_pool if repr(genome) not in attempted][
            :MIN_TRAINING_ROWS
        ]
        if not refill:
            valid = valid_training_row_count(evidence, 'fuel_cell_orr')
            raise RuntimeError(
                f'calibration exhausted after {len(attempted)} distinct probes; '
                f'only {valid}/{MIN_TRAINING_ROWS} valid ORR rows'
            )
        refill_round += 1
        attempted.update(repr(genome) for genome in refill)
        extra = run_orr_screening(
            refill,
            db_filename=f'fc_branch_calibration_refill_{refill_round}.csv',
            workers_per_gpu=2,
        )
        evidence = merge_compatible_evidence(extra, evidence, SCREENING_PROTOCOL_ID)
    save_screening_db(evidence, 'fc_branch_calibration.csv', subdir='fuel_cell')
    model = fit_tree_ranker(evidence, 'fuel_cell_orr')
    logger.info(
        'ORR ranker diagnostics: %s; acquisition_mode=%s',
        model.diagnostics,
        model.diagnostics.acquisition_mode,
    )
    score_population = lambda pop: orr_tree_objectives(pop, model)
    summary = run_branch_and_bound(
        BranchConfig(
            application='fuel_cell_orr',
            database=config.exhaustive_db,
            leaf_size=config.branch_leaf_size,
            probe_count=config.branch_probe_count,
            scan_batch_size=config.exhaustive_batch_size,
            max_leaves=config.branch_max_leaves,
            expected_population=config.expected_space_size,
            certificate_path=str(FUEL_CELL_DIR / 'coverage_certificate.json'),
            max_runtime_s=config.max_runtime_s,
            min_resolved_leaves_per_class=config.min_resolved_leaves_per_class,
            exploration_interval=config.branch_exploration_interval,
            refresh_pending_priorities=config.refresh_pending_priorities,
            scan_workers=config.scan_workers,
        ),
        score_population,
    )
    logger.info(f"ORR branch discovery: {summary}")
    archive = load_archive_genomes(
        config.exhaustive_db, 'fuel_cell_orr', config.htvs_pool_size
    )
    if not archive:
        return [], add_discovery_metadata(evidence)
    objectives = score_population(archive)
    fronts = fast_non_dominated_sort(objectives)
    champions = [archive[i] for i in fronts[0]]
    from pipeline.search.adaptive_validation import (
        allocate_validation_batch,
        experimental_slate,
        persist_experimental_slate,
        record_screening_frame,
    )

    _, uncertainty = model.predict(archive)
    validate_idx = allocate_validation_batch(
        archive,
        objectives,
        min(config.fairchem_eval_top_k, len(archive)),
        config.exhaustive_db,
        'fuel_cell_orr',
        min_per_class=config.min_validation_per_class,
        uncertainties=uncertainty,
    )
    validated = run_orr_screening(
        [archive[i] for i in validate_idx],
        db_filename='fc_branch_champions.csv',
        workers_per_gpu=2,
    )
    predicted_mean, _ = model.predict([archive[i] for i in validate_idx])
    predictions = {
        candidate_id(archive[i]): float(predicted_mean[position])
        for position, i in enumerate(validate_idx)
    }
    record_screening_frame(
        config.exhaustive_db,
        'fuel_cell_orr',
        predictions,
        validated,
        'orr_overpotential_V',
        'fairchem',
        0.40,
    )
    evidence = pd.concat([evidence, validated], ignore_index=True)
    evidence = merge_compatible_evidence(evidence, None, SCREENING_PROTOCOL_ID)
    save_screening_db(evidence, 'fc_branch_ranker_evidence.csv', subdir='fuel_cell')
    if config.prior_art_db:
        from pipeline.evidence.prior_art import annotate_prior_art

        evidence = annotate_prior_art(evidence, config.prior_art_db)
    slate_idx = experimental_slate(
        archive, objectives, min(config.fairchem_eval_top_k, len(archive))
    )
    persist_experimental_slate(
        config.exhaustive_db, 'fuel_cell_orr', archive, objectives, slate_idx
    )
    evidence.attrs['experimental_slate'] = [archive[i] for i in slate_idx]
    return champions, add_discovery_metadata(evidence)
