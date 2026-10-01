"""Shared orchestration for catalyst-guided production branch searches."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Optional, Sequence, TypeAlias

import numpy as np
import pandas as pd

from pipeline.search.adaptive_validation import (
    allocate_validation_batch,
    experimental_slate,
    persist_experimental_slate,
    record_screening_frame,
)
from pipeline.search.branch_search import BranchConfig, run_branch_and_bound
from pipeline.search.discovery import (
    add_discovery_metadata,
    candidate_id,
    fast_non_dominated_sort,
)
from pipeline.search.exhaustive_search import load_archive_genomes
from pipeline.search.indexed_space import deterministic_tree_probes
from pipeline.screening.small_data_ranker import (
    MIN_TRAINING_ROWS,
    TreeRanker,
    fit_tree_ranker,
    guide_branch_objectives,
    merge_compatible_evidence,
    valid_training_row_count,
)
from pipeline.utils import load_screening_db, save_screening_db

Genome: TypeAlias = tuple
Population: TypeAlias = Sequence[Genome]
Screener: TypeAlias = Callable[..., pd.DataFrame]
ObjectiveFunction: TypeAlias = Callable[[Population, TreeRanker], np.ndarray]


@dataclass
class BranchDiscoveryConfig:
    """Configure the common catalyst-guided branch workflow."""

    initial_fairchem_samples: int = 500
    fairchem_eval_top_k: int = 500
    htvs_pool_size: int = 20000
    exhaustive_batch_size: int = 65536
    exhaustive_db: str = ""
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


@dataclass(frozen=True)
class BranchDiscoveryWorkflow:
    """Declare the application-specific parts of branch discovery."""

    application: str
    outcome: str
    protocol_id: str
    calibration_db: str
    refill_prefix: str
    evidence_db: str
    champions_db: str
    database_subdir: str
    certificate_path: Path
    validation_fidelity: float
    diagnostics_label: str
    screener: Screener
    objective_function: ObjectiveFunction
    logger: logging.Logger


def run_guided_branch_discovery(
    config: BranchDiscoveryConfig,
    workflow: BranchDiscoveryWorkflow,
    existing_db: Optional[pd.DataFrame] = None,
) -> tuple[list[Genome], pd.DataFrame]:
    """Run calibration, branch scheduling, validation, and evidence persistence.

    Args:
        config: Shared branch-search resource and coverage settings.
        workflow: Application-specific scientific operations and artifact paths.
        existing_db: Optional completed evidence to use instead of fresh probes.

    Returns:
        Pareto champion genomes and the annotated screening evidence.
    """
    if existing_db is not None and len(existing_db) > 50:
        evidence = existing_db
    else:
        evidence = workflow.screener(
            deterministic_tree_probes(config.initial_fairchem_samples),
            db_filename=workflow.calibration_db,
            workers_per_gpu=2,
        )

    evidence = merge_compatible_evidence(
        evidence,
        load_screening_db(workflow.evidence_db, subdir=workflow.database_subdir),
        workflow.protocol_id,
    )
    attempted = {str(value) for value in evidence.get('genome', [])}
    refill_limit = max(
        config.initial_fairchem_samples * 3,
        config.initial_fairchem_samples + MIN_TRAINING_ROWS,
    )
    probe_pool = deterministic_tree_probes(refill_limit)
    refill_round = 0
    while valid_training_row_count(evidence, workflow.application) < MIN_TRAINING_ROWS:
        refill = [genome for genome in probe_pool if repr(genome) not in attempted][
            :MIN_TRAINING_ROWS
        ]
        if not refill:
            valid = valid_training_row_count(evidence, workflow.application)
            raise RuntimeError(
                f'calibration exhausted after {len(attempted)} distinct probes; '
                f'only {valid}/{MIN_TRAINING_ROWS} valid {workflow.application} rows'
            )
        refill_round += 1
        attempted.update(repr(genome) for genome in refill)
        extra = workflow.screener(
            refill,
            db_filename=f'{workflow.refill_prefix}_{refill_round}.csv',
            workers_per_gpu=2,
        )
        evidence = merge_compatible_evidence(extra, evidence, workflow.protocol_id)

    save_screening_db(
        evidence, workflow.calibration_db, subdir=workflow.database_subdir
    )
    model = fit_tree_ranker(evidence, workflow.application)
    workflow.logger.info(
        '%s ranker diagnostics: %s; acquisition_mode=%s',
        workflow.diagnostics_label,
        model.diagnostics,
        model.diagnostics.acquisition_mode,
    )

    def score_population(population: Population) -> np.ndarray:
        return workflow.objective_function(population, model)

    def score_branches(population: Population) -> np.ndarray:
        return guide_branch_objectives(
            score_population(population), evidence, population, workflow.application
        )

    summary = run_branch_and_bound(
        BranchConfig(
            application=workflow.application,
            database=config.exhaustive_db,
            leaf_size=config.branch_leaf_size,
            probe_count=config.branch_probe_count,
            scan_batch_size=config.exhaustive_batch_size,
            max_leaves=config.branch_max_leaves,
            expected_population=config.expected_space_size,
            certificate_path=str(workflow.certificate_path),
            max_runtime_s=config.max_runtime_s,
            min_resolved_leaves_per_class=config.min_resolved_leaves_per_class,
            exploration_interval=config.branch_exploration_interval,
            refresh_pending_priorities=config.refresh_pending_priorities,
            scan_workers=config.scan_workers,
        ),
        score_branches,
    )
    workflow.logger.info('%s branch discovery: %s', workflow.diagnostics_label, summary)
    archive = load_archive_genomes(
        config.exhaustive_db, workflow.application, config.htvs_pool_size
    )
    if not archive:
        return [], add_discovery_metadata(evidence)

    objectives = score_population(archive)
    champions = [archive[index] for index in fast_non_dominated_sort(objectives)[0]]
    validation_count = min(config.fairchem_eval_top_k, len(archive))
    _, uncertainty = model.predict(archive)
    validation_indices = allocate_validation_batch(
        archive,
        objectives,
        validation_count,
        config.exhaustive_db,
        workflow.application,
        min_per_class=config.min_validation_per_class,
        uncertainties=uncertainty,
    )
    validation_genomes = [archive[index] for index in validation_indices]
    validated = workflow.screener(
        validation_genomes,
        db_filename=workflow.champions_db,
        workers_per_gpu=2,
    )
    predicted_mean, _ = model.predict(validation_genomes)
    predictions = {
        candidate_id(genome): float(predicted_mean[position])
        for position, genome in enumerate(validation_genomes)
    }
    record_screening_frame(
        config.exhaustive_db,
        workflow.application,
        predictions,
        validated,
        workflow.outcome,
        'fairchem',
        workflow.validation_fidelity,
    )
    evidence = merge_compatible_evidence(
        pd.concat([evidence, validated], ignore_index=True),
        None,
        workflow.protocol_id,
    )
    save_screening_db(evidence, workflow.evidence_db, subdir=workflow.database_subdir)
    if config.prior_art_db:
        from pipeline.evidence.prior_art import annotate_prior_art

        evidence = annotate_prior_art(evidence, config.prior_art_db)
    slate_indices = experimental_slate(archive, objectives, validation_count)
    persist_experimental_slate(
        config.exhaustive_db,
        workflow.application,
        archive,
        objectives,
        slate_indices,
    )
    evidence.attrs['experimental_slate'] = [archive[index] for index in slate_indices]
    return champions, add_discovery_metadata(evidence)
