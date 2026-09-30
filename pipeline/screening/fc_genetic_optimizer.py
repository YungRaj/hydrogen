#!/usr/bin/env python3
# Fuel-cell objectives and branch-search orchestration.
"""Fuel-cell ORR objectives and catalyst-guided branch discovery.

The Pareto helpers remain part of objective and archive selection. Genetic and
random candidate generation are retired; the compatibility entry point below
fails explicitly instead of retaining an unreachable implementation.
"""

import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Optional

from pipeline.utils import setup_logger, FUEL_CELL_DIR
from pipeline.search.design_space import (
    encode_population,
    FEATURE_DIM,
)
from pipeline.search.discovery import (
    add_discovery_metadata,
    candidate_id,
)
import torch

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
# ORR-SPECIFIC OBJECTIVES & SURROGATE DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════════


class ORRCatalystSurrogate(torch.nn.Module):
    """
    Custom surrogate neural network for Fuel Cell ORR catalyst property prediction.
    Shared backbone → validity, ORR overpotential, and binding stability heads.
    """

    def __init__(
        self, input_dim: int = FEATURE_DIM, hidden_dims: tuple = (512, 256, 128)
    ):
        super().__init__()
        import torch.nn as nn

        layers = []
        prev_dim = input_dim
        for h_dim in hidden_dims:
            layers.extend(
                [
                    nn.Linear(prev_dim, h_dim),
                    nn.BatchNorm1d(h_dim),
                    nn.GELU(),
                    nn.Dropout(0.1),
                ]
            )
            prev_dim = h_dim
        self.backbone = nn.Sequential(*layers)
        self.head_valid = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_orr_eta = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )
        self.head_binding = nn.Sequential(
            nn.Linear(prev_dim, 32), nn.GELU(), nn.Linear(32, 1)
        )

    def forward(self, x: torch.Tensor):
        """Evaluate the neural-network forward pass.

        Args:
            x: Feature matrix or tensor consumed by the fitted model.

        Returns:
            The network output tensor for the supplied batch.
        """
        features = self.backbone(x)
        valid_logit = self.head_valid(features)
        orr_eta = self.head_orr_eta(features)
        binding = self.head_binding(features)
        return valid_logit, orr_eta, binding


class ORRSurrogateEnsemble(torch.nn.Module):
    """Ensemble of ORRCatalystSurrogate models for epistemic uncertainty estimation."""

    def __init__(self, n_models: int = 3, input_dim: int = FEATURE_DIM):
        super().__init__()
        self.models = torch.nn.ModuleList(
            [ORRCatalystSurrogate(input_dim=input_dim) for _ in range(n_models)]
        )


def _fenton_from_genome(genome: tuple) -> float:
    """Compute Fenton stability score from genome elements."""
    FENTON_RISK = {'Fe': 3, 'Cu': 2, 'Co': 1, 'Mn': 1, 'Cr': 1, 'V': 1}
    elements = _extract_elements_from_genome(genome)
    fenton_risk = sum(FENTON_RISK.get(e, 0) for e in elements)
    return float(max(0, 10 - fenton_risk))


def compute_orr_objectives_surrogate(
    population: List[tuple], model, device: str
) -> np.ndarray:
    """
        Compute 4 ORR objectives for a population using the ORR surrogate NN or Ensemble.

    Args:
        population: Ordered values supplying population.
        model: Fitted model used for inference.
        device: CPU or GPU device requested for execution.

    Returns:
        Computed `np.ndarray` result.
    """
    from pipeline.screening.ood import compute_model_confidence, confidence_penalty

    features = encode_population(population)
    import torch

    X = torch.FloatTensor(features).to(device)

    if isinstance(model, ORRSurrogateEnsemble):
        preds_eta_list = []
        preds_bind_list = []
        p_valid_list = []
        for single_model in model.models:
            single_model.eval()
            with torch.no_grad():
                valid_logit, pred_eta, pred_binding = single_model(X)
            p_valid_list.append(torch.sigmoid(valid_logit).cpu().numpy().flatten())
            preds_eta_list.append(pred_eta.cpu().numpy().flatten())
            preds_bind_list.append(pred_binding.cpu().numpy().flatten())

        # Aggregate with acquisition UCB/LCB (kappa = 1.0)
        p_valid = np.column_stack(p_valid_list).mean(axis=1)
        eta_arr = np.column_stack(preds_eta_list)
        # Minimize overpotential → LCB = mean - std
        pred_eta = eta_arr.mean(axis=1) - 1.0 * eta_arr.std(axis=1)
        bind_arr = np.column_stack(preds_bind_list)
        # Maximize binding strength → negate for minimization → UCB = mean + std
        pred_binding = bind_arr.mean(axis=1) + 1.0 * bind_arr.std(axis=1)
    else:
        model.eval()
        with torch.no_grad():
            valid_logit, pred_eta, pred_binding = model(X)
        p_valid = torch.sigmoid(valid_logit).cpu().numpy().flatten()
        pred_eta = pred_eta.cpu().numpy().flatten()
        pred_binding = pred_binding.cpu().numpy().flatten()

    n = len(population)
    objectives = np.zeros((n, 4))

    for i in range(n):
        from pipeline.search.scope import pemfc_cathode_scope

        in_scope = pemfc_cathode_scope(population[i])['status'] == 'candidate'
        if p_valid[i] > 0.3 and in_scope:
            elements = _extract_elements_from_genome(population[i])
            conf = compute_model_confidence(population[i], elements)
            penalty = confidence_penalty(conf)

            objectives[i, 0] = pred_eta[i] + penalty  # additive OOD penalty
            objectives[i, 1] = -_fenton_from_genome(population[i])
            objectives[i, 2] = _cost_from_genome(population[i])
            objectives[i, 3] = -pred_binding[i]
        else:
            objectives[i, :] = [5.0, 0.0, 100.0, 0.0]  # penalty

    return objectives


def _cost_from_genome(genome: tuple) -> float:
    """Compute cost penalty from genome elements.

    abundance_cost_penalty() returns [-2, 0] where 0 = abundant, -2 = rare.
    Since NSGA-II minimizes all objectives, we negate so that:
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
# NSGA-II (reused from methane GA — same algorithm)
# ═══════════════════════════════════════════════════════════════════════════════


def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """NSGA-II fast non-dominated sorting using vectorized Pareto front extraction.

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


def crowding_distance(objectives: np.ndarray, front: List[int]) -> np.ndarray:
    """Compute crowding distances for a Pareto front.

    Args:
        objectives: Objectives used by this operation.
        front: Ordered values supplying front.

    Returns:
        Computed `np.ndarray` result.
    """
    n = len(front)
    if n <= 2:
        return np.full(n, np.inf)

    distances = np.zeros(n)
    m = objectives.shape[1]

    for obj_idx in range(m):
        sorted_indices = np.argsort(objectives[front, obj_idx])
        distances[sorted_indices[0]] = np.inf
        distances[sorted_indices[-1]] = np.inf

        obj_range = (
            objectives[front[sorted_indices[-1]], obj_idx]
            - objectives[front[sorted_indices[0]], obj_idx]
        )
        if obj_range < 1e-10:
            continue

        for i in range(1, n - 1):
            distances[sorted_indices[i]] += (
                objectives[front[sorted_indices[i + 1]], obj_idx]
                - objectives[front[sorted_indices[i - 1]], obj_idx]
            ) / obj_range

    return distances


def nsga2_select(population, objectives, n_select):
    """NSGA-II selection with non-dominated sorting + crowding distance.

    Args:
        population: Population used by this operation.
        objectives: Objectives used by this operation.
        n_select: Number of select to use.

    Returns:
        Computed result described above.
    """
    fronts = fast_non_dominated_sort(objectives)
    selected = []

    for front in fronts:
        if len(selected) + len(front) <= n_select:
            selected.extend(front)
        else:
            remaining = n_select - len(selected)
            if remaining > 0:
                cd = crowding_distance(objectives, front)
                top_cd = np.argsort(-cd)[:remaining]
                selected.extend([front[i] for i in top_cd])
            break

    return selected


# ═══════════════════════════════════════════════════════════════════════════════
# ORR SURROGATE TRAINING
# ═══════════════════════════════════════════════════════════════════════════════






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


def run_fc_genetic_algorithm(*_args: object, **_kwargs: object) -> None:
    """Reject the retired fuel-cell genetic-search entry point.

    Catalyst-guided branch discovery is the only supported candidate search.
    """
    raise RuntimeError(
        "Genetic/random candidate search was retired; use run_fc_branch_discovery()"
    )
