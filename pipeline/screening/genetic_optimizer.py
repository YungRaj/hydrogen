#!/usr/bin/env python3
# Turquoise-hydrogen objectives and branch-search orchestration.
"""Turquoise-hydrogen objectives and catalyst-guided branch discovery.

The Pareto helper remains part of archive selection. Candidate generation is
handled exclusively by catalyst-guided branch discovery.
"""

import numpy as np
import pandas as pd
from typing import List, Optional
from dataclasses import dataclass

from pipeline.utils import (
    setup_logger,
    save_screening_db,
    load_screening_db,
    SCREENING_DIR,
)
from pipeline.search.discovery import (
    add_discovery_metadata,
    candidate_id,
)

logger = setup_logger('genetic_optimizer', 'screening/genetic_optimizer.log')


# ═══════════════════════════════════════════════════════════════════════════════
# PARETO ARCHIVE SELECTION
# ═══════════════════════════════════════════════════════════════════════════════




def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """
        Fast non-dominated sort using vectorized Pareto-front extraction.

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
# OBJECTIVE COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════════


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
        pass  # no metals
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# BRANCH DISCOVERY
# ═══════════════════════════════════════════════════════════════════════════════



@dataclass
class BranchDiscoveryConfig:
    """Configure the production branch-search adapter.

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
    device: str = 'cuda:0'
    exhaustive_batch_size: int = 65536
    exhaustive_db: str = str(SCREENING_DIR / 'indexed_scan.sqlite')
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



def run_branch_discovery(
    config: BranchDiscoveryConfig = BranchDiscoveryConfig(),
    existing_db: Optional[pd.DataFrame] = None,
):
    """Run the standard catalyst-guided turquoise-hydrogen branch search.

    Args:
        config: Configuration controlling this operation.
        existing_db: Existing db used by this operation.

    Returns:
        Computed result described above.
    """
    from pipeline.search.indexed_space import deterministic_tree_probes
    from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
    from pipeline.search.exhaustive_search import load_archive_genomes
    from pipeline.screening.surface_screener import run_screening, SCREENING_PROTOCOL_ID

    if existing_db is not None and len(existing_db) > 50:
        evidence = existing_db
    else:
        probes = deterministic_tree_probes(config.initial_fairchem_samples)
        evidence = run_screening(
            probes, db_filename='branch_calibration.csv', workers_per_gpu=2
        )
    from pipeline.screening.small_data_ranker import (
        MIN_TRAINING_ROWS,
        fit_tree_ranker,
        merge_compatible_evidence,
        turquoise_tree_objectives,
        valid_training_row_count,
    )

    prior_evidence = load_screening_db('branch_ranker_evidence.csv')
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
    while valid_training_row_count(evidence, 'turquoise_hydrogen') < MIN_TRAINING_ROWS:
        refill = [genome for genome in probe_pool if repr(genome) not in attempted][
            :MIN_TRAINING_ROWS
        ]
        if not refill:
            valid = valid_training_row_count(evidence, 'turquoise_hydrogen')
            raise RuntimeError(
                f'calibration exhausted after {len(attempted)} distinct probes; '
                f'only {valid}/{MIN_TRAINING_ROWS} valid turquoise-hydrogen rows'
            )
        refill_round += 1
        attempted.update(repr(genome) for genome in refill)
        extra = run_screening(
            refill,
            db_filename=f'branch_calibration_refill_{refill_round}.csv',
            workers_per_gpu=2,
        )
        evidence = merge_compatible_evidence(extra, evidence, SCREENING_PROTOCOL_ID)
    save_screening_db(evidence, 'branch_calibration.csv')
    model = fit_tree_ranker(evidence, 'turquoise_hydrogen')
    logger.info(
        'Pyrolysis ranker diagnostics: %s; acquisition_mode=%s',
        model.diagnostics,
        model.diagnostics.acquisition_mode,
    )
    score_population = lambda pop: turquoise_tree_objectives(pop, model)
    summary = run_branch_and_bound(
        BranchConfig(
            application='turquoise_hydrogen',
            database=config.exhaustive_db,
            leaf_size=config.branch_leaf_size,
            probe_count=config.branch_probe_count,
            scan_batch_size=config.exhaustive_batch_size,
            max_leaves=config.branch_max_leaves,
            expected_population=config.expected_space_size,
            certificate_path=str(
                SCREENING_DIR / 'turquoise_hydrogen_coverage_certificate.json'
            ),
            max_runtime_s=config.max_runtime_s,
            min_resolved_leaves_per_class=config.min_resolved_leaves_per_class,
            exploration_interval=config.branch_exploration_interval,
            refresh_pending_priorities=config.refresh_pending_priorities,
            scan_workers=config.scan_workers,
        ),
        score_population,
    )
    logger.info(f"Branch discovery: {summary}")
    archive = load_archive_genomes(
        config.exhaustive_db, 'turquoise_hydrogen', config.htvs_pool_size
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
        'turquoise_hydrogen',
        min_per_class=config.min_validation_per_class,
        uncertainties=uncertainty,
    )
    validated = run_screening(
        [archive[i] for i in validate_idx],
        db_filename='branch_champions.csv',
        workers_per_gpu=2,
    )
    predicted_mean, _ = model.predict([archive[i] for i in validate_idx])
    predictions = {
        candidate_id(archive[i]): float(predicted_mean[position])
        for position, i in enumerate(validate_idx)
    }
    record_screening_frame(
        config.exhaustive_db,
        'turquoise_hydrogen',
        predictions,
        validated,
        'E_act',
        'fairchem',
        0.8,
    )
    evidence = pd.concat([evidence, validated], ignore_index=True)
    evidence = merge_compatible_evidence(evidence, None, SCREENING_PROTOCOL_ID)
    save_screening_db(evidence, 'branch_ranker_evidence.csv')
    if config.prior_art_db:
        from pipeline.evidence.prior_art import annotate_prior_art

        evidence = annotate_prior_art(evidence, config.prior_art_db)
    slate_idx = experimental_slate(
        archive, objectives, min(config.fairchem_eval_top_k, len(archive))
    )
    persist_experimental_slate(
        config.exhaustive_db, 'turquoise_hydrogen', archive, objectives, slate_idx
    )
    evidence.attrs['experimental_slate'] = [archive[i] for i in slate_idx]
    return champions, add_discovery_metadata(evidence)
