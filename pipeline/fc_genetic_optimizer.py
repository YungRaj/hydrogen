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
from pipeline.surrogate_model import CatalystSurrogate, train_surrogate, predict_batch

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
    device: str = 'cuda'
    seed: int = 42


# ═══════════════════════════════════════════════════════════════════════════════
# ORR-SPECIFIC OBJECTIVES
# ═══════════════════════════════════════════════════════════════════════════════

def compute_orr_objectives_surrogate(population: List[tuple], model, device: str) -> np.ndarray:
    """
    Compute 4 ORR objectives for a population using the surrogate NN.

    Objectives (all minimized):
      0: ORR overpotential (η) — lower is better
      1: -Fenton stability — lower (more negative) = more stable
      2: Cost penalty — lower is cheaper
      3: -Binding strength — lower (more negative) = more durable

    For candidates the surrogate predicts as invalid, assign penalty values.
    """
    features = encode_population(population)
    import torch
    X = torch.FloatTensor(features).to(device)

    model.eval()
    with torch.no_grad():
        preds = model(X)

    n = len(population)
    objectives = np.zeros((n, 4))

    # The surrogate outputs: [valid_prob, dG_OH, orr_eta, fenton_stability]
    p_valid = torch.sigmoid(preds[:, 0]).cpu().numpy()
    pred_eta = preds[:, 1].cpu().numpy()
    pred_fenton = preds[:, 2].cpu().numpy()
    pred_binding = preds[:, 3].cpu().numpy()

    for i in range(n):
        if p_valid[i] > 0.3:
            objectives[i, 0] = pred_eta[i]         # minimize overpotential
            objectives[i, 1] = -pred_fenton[i]      # minimize -fenton (maximize stability)
            objectives[i, 2] = _cost_from_genome(population[i])
            objectives[i, 3] = -pred_binding[i]     # minimize -binding (maximize durability)
        else:
            objectives[i, :] = [5.0, 0.0, 100.0, 0.0]  # penalty

    return objectives


def _cost_from_genome(genome: tuple) -> float:
    """Compute cost penalty from genome elements."""
    from pipeline.utils import abundance_cost_penalty
    elements = _extract_elements_from_genome(genome)
    return abundance_cost_penalty(elements)


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
    """Train surrogate NN from ORR MACE screening data."""
    import torch

    valid_mask = db['valid'] == True
    valid_db = db[valid_mask].copy()

    if len(valid_db) < 20:
        logger.warning(f"Only {len(valid_db)} valid ORR samples, too few for surrogate")
        return None

    # Parse genomes and encode
    genomes = []
    for _, row in db.iterrows():
        try:
            g = eval(row['genome'])
            genomes.append(g)
        except Exception:
            genomes.append(None)

    features = []
    targets = []  # [valid, orr_eta, fenton, binding]

    for i, g in enumerate(genomes):
        if g is None:
            continue
        try:
            feat = encode_genome(g)
            row = db.iloc[i]
            valid = 1.0 if row.get('valid', False) else 0.0
            eta = float(row.get('orr_overpotential_V', 5.0)) if valid else 5.0
            fenton = float(row.get('fenton_stability', 5.0)) if valid else 0.0
            binding = float(row.get('binding_strength', 0.0)) if valid else 0.0
            features.append(feat)
            targets.append([valid, eta, fenton, binding])
        except Exception:
            continue

    if len(features) < 20:
        return None

    X = torch.FloatTensor(np.array(features)).to(device)
    Y = torch.FloatTensor(np.array(targets)).to(device)

    model = CatalystSurrogate(input_dim=FEATURE_DIM, output_dim=4).to(device)
    model = train_surrogate(model, X, Y, epochs=200, lr=1e-3)
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
        from pipeline.fc_mace_screener import run_orr_screening
        all_mace_results = run_orr_screening(
            initial_pop, db_filename="fc_initial_screening.csv", workers_per_gpu=2
        )

    # ── Phase B: Train ORR Surrogate ────────────────────────────────────────
    model = _train_orr_surrogate(all_mace_results, config.device)

    # ── Phase C: Evolutionary Loop ──────────────────────────────────────────
    population = generate_population(config.pop_size)

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
            logger.info(f"  Gen {gen}: Running MACE ORR validation on top {config.mace_eval_top_k}...")

            fronts = fast_non_dominated_sort(final_obj)
            top_indices = []
            for front in fronts:
                top_indices.extend(front)
                if len(top_indices) >= config.mace_eval_top_k:
                    break
            top_genomes = [population[i] for i in top_indices[:config.mace_eval_top_k]]

            from pipeline.fc_mace_screener import run_orr_screening
            mace_df = run_orr_screening(
                top_genomes, db_filename=f"fc_mace_gen{gen}.csv", workers_per_gpu=2
            )
            all_mace_results = pd.concat([all_mace_results, mace_df], ignore_index=True)

        # ── Retrain Surrogate ───────────────────────────────────────────
        if gen % config.surrogate_retrain_interval == 0:
            logger.info(f"  Gen {gen}: Retraining ORR surrogate on {len(all_mace_results)} samples...")
            model = _train_orr_surrogate(all_mace_results, config.device)

        # ── Logging ─────────────────────────────────────────────────────
        if gen % 10 == 0 or gen == 1:
            best_eta = final_obj[:, 0].min()
            pareto = len(fronts[0]) if 'fronts' in dir() else 0
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
