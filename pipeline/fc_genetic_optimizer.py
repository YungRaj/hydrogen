#!/usr/bin/env python3
"""
NSGA-II Genetic Algorithm for Fuel Cell ORR Cathode Catalysts.

4-objective Pareto optimization:
  1. Minimize ORR overpotential (η → 0 = ideal)
  2. Maximize Fenton stability (radical resistance)
  3. Minimize cost (crustal abundance penalty)
  4. Maximize binding strength (dissolution resistance)

Uses the same 25.3B design space as the methane pyrolysis GA,
but trains a separate surrogate NN and evaluates ORR descriptors.
"""

import os
import sys
import ast
import time
import random
import logging
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import List, Optional
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from pipeline.utils import setup_logger, save_json, BASE_DIR
from pipeline.catalyst_spaces import (
    generate_population, crossover, mutate, encode_genome, encode_population,
    ALL_MATERIAL_CLASSES, FEATURE_DIM,
)
import torch

logger = setup_logger('fc_genetic_optimizer', 'fuel_cell/fc_genetic_optimizer.log')


@dataclass
class FCGAConfig:
    """Configuration for fuel cell ORR genetic algorithm."""
    pop_size: int = 1000
    n_generations: int = 3000
    initial_mace_samples: int = 500
    mace_eval_interval: int = 5
    mace_eval_top_k: int = 500
    surrogate_retrain_interval: int = 5
    mutation_rate: float = 0.35
    crossover_rate: float = 0.7
    explore_interval: int = 3           # Run exploration shots every N MACE intervals
    explore_per_class: int = 5          # Random GNN evaluations per class during exploration
    device: str = 'cuda'
    seed: int = 42


# ═══════════════════════════════════════════════════════════════════════════════
# ORR-SPECIFIC OBJECTIVES & SURROGATE DEFINITIONS
# ═══════════════════════════════════════════════════════════════════════════════

class ORRCatalystSurrogate(torch.nn.Module):
    """
    Custom surrogate neural network for Fuel Cell ORR catalyst property prediction.
    Shared backbone → validity, ORR overpotential, and binding stability heads.
    """
    def __init__(self, input_dim: int = FEATURE_DIM, hidden_dims: tuple = (512, 256, 128)):
        super().__init__()
        import torch.nn as nn
        layers = []
        prev_dim = input_dim
        for h_dim in hidden_dims:
            layers.extend([
                nn.Linear(prev_dim, h_dim),
                nn.BatchNorm1d(h_dim),
                nn.GELU(),
                nn.Dropout(0.1),
            ])
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
        features = self.backbone(x)
        valid_logit = self.head_valid(features)
        orr_eta = self.head_orr_eta(features)
        binding = self.head_binding(features)
        return valid_logit, orr_eta, binding


def _fenton_from_genome(genome: tuple) -> float:
    """Compute Fenton stability score from genome elements.
    
    Matches FENTON_RISK table from fc_screener.py.
    Higher score = more stable (less radical degradation risk).
    """
    # Canonical Fenton risk weights (from literature on radical generation in acid)
    FENTON_RISK = {'Fe': 3, 'Cu': 2, 'Co': 1, 'Mn': 1, 'Cr': 1, 'V': 1}
    elements = _extract_elements_from_genome(genome)
    fenton_risk = sum(FENTON_RISK.get(e, 0) for e in elements)
    return float(max(0, 10 - fenton_risk))


def compute_orr_objectives_surrogate(population: List[tuple], model, device: str) -> np.ndarray:
    """
    Compute 4 ORR objectives for a population using the ORR surrogate NN.

    Objectives (all minimized):
      0: ORR overpotential (η) × confidence_penalty — lower is better
      1: -Fenton stability — lower (more negative) = more stable
      2: Cost penalty — lower is cheaper
      3: -Binding strength — lower (more negative) = more durable

    The overpotential is scaled by an OOD confidence penalty so that
    predictions from material classes outside the eSen-SM training
    distribution (MOFs, MetalFreeCarbon, etc.) are discounted.
    """
    from pipeline.ood_detector import compute_model_confidence, confidence_penalty

    features = encode_population(population)
    import torch
    X = torch.FloatTensor(features).to(device)

    model.eval()
    with torch.no_grad():
        valid_logit, pred_eta, pred_binding = model(X)

    n = len(population)
    objectives = np.zeros((n, 4))

    p_valid = torch.sigmoid(valid_logit).cpu().numpy().flatten()
    pred_eta = pred_eta.cpu().numpy().flatten()
    pred_binding = pred_binding.cpu().numpy().flatten()

    for i in range(n):
        if p_valid[i] > 0.3:
            elements = _extract_elements_from_genome(population[i])
            conf = compute_model_confidence(population[i], elements)
            penalty = confidence_penalty(conf)

            objectives[i, 0] = pred_eta[i] + penalty   # additive OOD penalty
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
        if genome[2] != 'None': elements.append(genome[2])
    elif mat_class == 'SolidCatalyst':
        elements.append(genome[1])
        for d in genome[5]: elements.append(d)
    elif mat_class == 'SAC':
        elements.append(genome[1])
    elif mat_class == 'DAC':
        elements.extend([genome[1], genome[2]])
    elif mat_class in ('MOF', 'COF'):
        if genome[1] != 'None': elements.append(genome[1])
    elif mat_class == 'Perovskite':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MetalHydride':
        elements.append(genome[1])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MAXPhase':
        elements.extend([genome[1], genome[2]])
        if genome[5] != 'None': elements.append(genome[5])
    elif mat_class == 'HEA':
        elements.extend(list(genome[1]))
    elif mat_class == 'Spinel':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MXene':
        elements.append(genome[1])
        if genome[5] != 'None': elements.append(genome[5])
    elif mat_class == 'SAA':
        elements.extend([genome[1], genome[2]])
    elif mat_class == 'MetalFreeCarbon':
        pass  # no metals — zero cost
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# NSGA-II (reused from methane GA — same algorithm)
# ═══════════════════════════════════════════════════════════════════════════════

def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """NSGA-II fast non-dominated sorting. All objectives are minimized."""
    n = len(objectives)
    domination_count = np.zeros(n, dtype=int)
    dominated_set = [[] for _ in range(n)]
    fronts = [[]]

    for p in range(n):
        for q in range(p + 1, n):
            p_dom_q = np.all(objectives[p] <= objectives[q]) and np.any(objectives[p] < objectives[q])
            q_dom_p = np.all(objectives[q] <= objectives[p]) and np.any(objectives[q] < objectives[p])
            if p_dom_q:
                dominated_set[p].append(q)
                domination_count[q] += 1
            elif q_dom_p:
                dominated_set[q].append(p)
                domination_count[p] += 1

    for i in range(n):
        if domination_count[i] == 0:
            fronts[0].append(i)

    k = 0
    while fronts[k]:
        next_front = []
        for p in fronts[k]:
            for q in dominated_set[p]:
                domination_count[q] -= 1
                if domination_count[q] == 0:
                    next_front.append(q)
        k += 1
        fronts.append(next_front)

    return [f for f in fronts if f]


def crowding_distance(objectives: np.ndarray, front: List[int]) -> np.ndarray:
    """Compute crowding distances for a Pareto front."""
    n = len(front)
    if n <= 2:
        return np.full(n, np.inf)

    distances = np.zeros(n)
    m = objectives.shape[1]

    for obj_idx in range(m):
        sorted_indices = np.argsort(objectives[front, obj_idx])
        distances[sorted_indices[0]] = np.inf
        distances[sorted_indices[-1]] = np.inf

        obj_range = (objectives[front[sorted_indices[-1]], obj_idx] -
                     objectives[front[sorted_indices[0]], obj_idx])
        if obj_range < 1e-10:
            continue

        for i in range(1, n - 1):
            distances[sorted_indices[i]] += (
                objectives[front[sorted_indices[i + 1]], obj_idx] -
                objectives[front[sorted_indices[i - 1]], obj_idx]
            ) / obj_range

    return distances


def nsga2_select(population, objectives, n_select):
    """NSGA-II selection with non-dominated sorting + crowding distance."""
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

def _train_orr_surrogate(db: pd.DataFrame, device: str):
    """Train ORR surrogate NN from ORR MACE screening data."""
    import torch
    import torch.nn as nn
    from torch.utils.data import DataLoader, TensorDataset

    valid_mask = db['valid'] == True
    valid_db = db[valid_mask].copy()

    if len(valid_db) < 20:
        logger.warning(f"Only {len(valid_db)} valid ORR samples, too few for surrogate")
        return None

    # Parse genomes and encode
    genomes = []
    for _, row in db.iterrows():
        try:
            g = ast.literal_eval(row['genome'])
            genomes.append(g)
        except Exception:
            genomes.append(None)

    features = []
    targets = []  # [valid, orr_eta, binding]

    for i, g in enumerate(genomes):
        if g is None:
            continue
        try:
            feat = encode_genome(g)
            row = db.iloc[i]
            valid = 1.0 if row.get('valid', False) else 0.0
            eta = float(row.get('orr_overpotential_V', 5.0)) if valid else 5.0
            binding = float(row.get('binding_strength', 0.0)) if valid else 0.0
            features.append(feat)
            targets.append([valid, eta, binding])
        except Exception:
            continue

    if len(features) < 20:
        return None

    X = np.array(features)
    Y = np.array(targets)

    # Train ORR surrogate
    model = ORRCatalystSurrogate(input_dim=FEATURE_DIM).to(device)

    X_t = torch.tensor(X, dtype=torch.float32).to(device)
    y_val_t = torch.tensor(Y[:, 0], dtype=torch.float32).unsqueeze(1).to(device)
    y_eta_t = torch.tensor(Y[:, 1], dtype=torch.float32).unsqueeze(1).to(device)
    y_bind_t = torch.tensor(Y[:, 2], dtype=torch.float32).unsqueeze(1).to(device)

    dataset = TensorDataset(X_t, y_val_t, y_eta_t, y_bind_t)
    loader = DataLoader(dataset, batch_size=256, shuffle=True, drop_last=True)

    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    bce_loss = nn.BCEWithLogitsLoss()
    mse_loss = nn.MSELoss()

    model.train()
    for epoch in range(100):
        for batch in loader:
            bX, bV, bE, bB = batch
            optimizer.zero_grad()
            valid_logit, pred_eta, pred_binding = model(bX)

            loss_v = bce_loss(valid_logit, bV)
            mask = (bV > 0.5).squeeze()
            if mask.sum() > 0:
                # Limit slicing indexing issues by checking dimension
                loss_e = mse_loss(pred_eta.squeeze(-1)[mask], bE.squeeze(-1)[mask])
                loss_b = mse_loss(pred_binding.squeeze(-1)[mask], bB.squeeze(-1)[mask])
                loss = loss_v + 2.0 * (loss_e + loss_b)
            else:
                loss = loss_v

            loss.backward()
            optimizer.step()

    model.eval()
    logger.info(f"ORR surrogate trained on {len(features)} samples")
    return model


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN GA LOOP
# ═══════════════════════════════════════════════════════════════════════════════

def run_fc_genetic_algorithm(config: FCGAConfig, existing_db=None):
    """
    Run NSGA-II genetic algorithm for ORR fuel cell cathode catalyst discovery.

    Returns: (pareto_genomes, screening_dataframe)
    """
    random.seed(config.seed)
    np.random.seed(config.seed)

    logger.info(f"ORR FC-GA: Population={config.pop_size}, Gens={config.n_generations}")

    # ── Phase A: Initial MACE ORR Screening ─────────────────────────────────
    if existing_db is not None and len(existing_db) > 50:
        logger.info(f"Using existing ORR database: {len(existing_db)} entries")
        all_mace_results = existing_db
    else:
        logger.info(f"Generating {config.initial_mace_samples} initial ORR MACE samples...")
        initial_pop = generate_population(config.initial_mace_samples)
        from pipeline.fc_screener import run_orr_screening
        all_mace_results = run_orr_screening(
            initial_pop, db_filename="fc_initial_screening.csv", workers_per_gpu=2
        )

    # ── Phase B: Train ORR Surrogate ────────────────────────────────────────
    model = _train_orr_surrogate(all_mace_results, config.device)

    # ── Phase C: Evolutionary Loop ──────────────────────────────────────────
    population = generate_population(config.pop_size)
    fronts = [[]]  # Initialize for logging before first sort
    mace_round = 0  # tracks MACE rounds for exploration scheduling

    for gen in range(1, config.n_generations + 1):
        t_gen = time.time()

        # Evaluate with surrogate
        if model is not None:
            objectives = compute_orr_objectives_surrogate(population, model, config.device)
        else:
            objectives = np.random.rand(len(population), 4)

        # NSGA-II selection
        selected_idx = nsga2_select(population, objectives, config.pop_size // 2)
        parents = [population[i] for i in selected_idx]

        # Generate offspring
        offspring = []
        while len(offspring) < config.pop_size:
            i1 = random.randint(0, len(parents) - 1)
            i2 = random.randint(0, len(parents) - 1)
            if random.random() < config.crossover_rate:
                child = crossover(parents[i1], parents[i2])
            else:
                child = parents[i1]
            child = mutate(child, rate=config.mutation_rate)
            offspring.append(child)

        # ── Class-Diversity Enforcement ─────────────────────────────────
        combined = parents + offspring
        min_per_class = max(2, config.pop_size // 20)
        class_counts = {}
        for g in combined:
            cls = g[0]
            class_counts[cls] = class_counts.get(cls, 0) + 1

        diversity_injections = []
        for cls in ALL_MATERIAL_CLASSES:
            deficit = min_per_class - class_counts.get(cls, 0)
            if deficit > 0:
                fresh = generate_population(deficit, material_class=cls)
                diversity_injections.extend(fresh)

        if diversity_injections:
            combined = combined[:len(combined) - len(diversity_injections)] + diversity_injections

        combined_obj = compute_orr_objectives_surrogate(combined, model, config.device) if model else np.random.rand(len(combined), 4)
        final_idx = nsga2_select(combined, combined_obj, config.pop_size)
        population = [combined[i] for i in final_idx]
        final_obj = combined_obj[final_idx]

        # ── Periodic MACE ORR Validation ────────────────────────────────
        if gen % config.mace_eval_interval == 0 or gen == 1:
            mace_round += 1
            logger.info(f"  Gen {gen}: Running MACE ORR validation on top {config.mace_eval_top_k}...")

            fronts = fast_non_dominated_sort(final_obj)
            top_indices = []
            for front in fronts:
                top_indices.extend(front)
                if len(top_indices) >= config.mace_eval_top_k:
                    break
            top_genomes = [population[i] for i in top_indices[:config.mace_eval_top_k]]

            from pipeline.fc_screener import run_orr_screening
            mace_df = run_orr_screening(
                top_genomes, db_filename=f"fc_mace_gen{gen}.csv", workers_per_gpu=2
            )
            all_mace_results = pd.concat([all_mace_results, mace_df], ignore_index=True)

            # ── Exploration Shots: probe EVERY class with real GNN ───────
            if mace_round % config.explore_interval == 0:
                explore_genomes = []
                for cls in ALL_MATERIAL_CLASSES:
                    explore_genomes.extend(
                        generate_population(config.explore_per_class, material_class=cls)
                    )
                n_explore = len(explore_genomes)
                logger.info(
                    f"  Gen {gen}: EXPLORATION — evaluating {n_explore} random "
                    f"candidates across all 14 classes with real GNN (ORR)..."
                )
                explore_df = run_orr_screening(
                    explore_genomes,
                    db_filename=f"fc_explore_gen{gen}.csv",
                    workers_per_gpu=2
                )
                all_mace_results = pd.concat([all_mace_results, explore_df], ignore_index=True)

                # Inject promising exploration candidates into population
                if 'orr_overpotential' in explore_df.columns:
                    good_explores = explore_df[
                        (explore_df['valid'] == True) &
                        (explore_df['orr_overpotential'] < explore_df['orr_overpotential'].quantile(0.3))
                    ]
                    if len(good_explores) > 0:
                        logger.info(
                            f"    Found {len(good_explores)} promising ORR exploration "
                            f"candidates — injecting into population"
                        )
                        for _, row in good_explores.iterrows():
                            try:
                                g = ast.literal_eval(row['genome'])
                                replace_idx = random.randint(0, len(population) - 1)
                                population[replace_idx] = g
                            except Exception:
                                pass

        # ── Retrain Surrogate ───────────────────────────────────────────
        if gen % config.surrogate_retrain_interval == 0:
            logger.info(f"  Gen {gen}: Retraining ORR surrogate on {len(all_mace_results)} samples...")
            model = _train_orr_surrogate(all_mace_results, config.device)

        # ── Logging ─────────────────────────────────────────────────────
        if gen % 10 == 0 or gen == 1:
            best_eta = final_obj[:, 0].min()
            pareto = len(fronts[0]) if fronts else 0
            n_unique = len(set(str(g) for g in population))
            elapsed = time.time() - t_gen
            logger.info(
                f"  Gen {gen}/{config.n_generations}: "
                f"best_η={best_eta:.4f} V, "
                f"pareto_size={pareto}, "
                f"pop_diversity={n_unique}, "
                f"({elapsed:.1f}s)"
            )

    # ── Return Results ──────────────────────────────────────────────────────
    final_objectives = compute_orr_objectives_surrogate(population, model, config.device) if model else np.random.rand(len(population), 4)
    fronts = fast_non_dominated_sort(final_objectives)
    pareto_genomes = [population[i] for i in fronts[0]]

    return pareto_genomes, all_mace_results
