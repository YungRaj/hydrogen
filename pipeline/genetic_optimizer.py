#!/usr/bin/env python3
"""
NSGA-II Multi-Objective Genetic Algorithm for Catalyst Discovery.

Optimizes 4 objectives simultaneously on the Pareto front:
  1. Minimize activation barrier (E_act)
  2. Maximize coking resistance index
  3. Maximize binding/segregation stability (negative = stable)
  4. Minimize material cost

Uses a surrogate model for rapid population evaluation, with periodic
retraining on MACE-validated data. Top Pareto-front candidates are
passed back to full MACE screening for ground-truth validation.

Implements NSGA-II selection with crowding distance for diverse Pareto fronts.
"""

import os
import time
import random
import numpy as np
import pandas as pd
import logging
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass

from pipeline.utils import (
    setup_logger, print_banner, save_screening_db, load_screening_db,
    abundance_cost_penalty, SCREENING_DIR,
)
from pipeline.catalyst_spaces import (
    generate_population, crossover, mutate, encode_genome, encode_population,
    ALL_MATERIAL_CLASSES, FEATURE_DIM,
)
from pipeline.surrogate_model import CatalystSurrogate, train_surrogate, predict_batch

logger = setup_logger('genetic_optimizer', 'screening/genetic_optimizer.log')


# ═══════════════════════════════════════════════════════════════════════════════
# NSGA-II IMPLEMENTATION
# ═══════════════════════════════════════════════════════════════════════════════

def dominates(obj_a: np.ndarray, obj_b: np.ndarray) -> bool:
    """Return True if solution a dominates solution b (all ≤, at least one <)."""
    return np.all(obj_a <= obj_b) and np.any(obj_a < obj_b)


def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """
    NSGA-II fast non-dominated sort.
    
    Args:
        objectives: (N, M) array where each row is an M-objective vector.
                    All objectives are MINIMIZED.
    Returns:
        List of fronts, where each front is a list of indices.
    """
    N = len(objectives)
    domination_count = np.zeros(N, dtype=int)  # how many solutions dominate i
    dominated_set = [[] for _ in range(N)]      # set of solutions i dominates

    for i in range(N):
        for j in range(i + 1, N):
            if dominates(objectives[i], objectives[j]):
                dominated_set[i].append(j)
                domination_count[j] += 1
            elif dominates(objectives[j], objectives[i]):
                dominated_set[j].append(i)
                domination_count[i] += 1

    fronts = []
    current_front = [i for i in range(N) if domination_count[i] == 0]
    while current_front:
        fronts.append(current_front)
        next_front = []
        for i in current_front:
            for j in dominated_set[i]:
                domination_count[j] -= 1
                if domination_count[j] == 0:
                    next_front.append(j)
        current_front = next_front

    return fronts


def crowding_distance(objectives: np.ndarray, front: List[int]) -> np.ndarray:
    """
    Compute crowding distance for individuals in a front.
    Used to maintain diversity on the Pareto front.
    """
    n = len(front)
    if n <= 2:
        return np.full(n, np.inf)

    distances = np.zeros(n)
    m = objectives.shape[1]

    for k in range(m):
        obj_vals = objectives[front, k]
        sorted_idx = np.argsort(obj_vals)
        distances[sorted_idx[0]] = np.inf
        distances[sorted_idx[-1]] = np.inf

        obj_range = obj_vals[sorted_idx[-1]] - obj_vals[sorted_idx[0]]
        if obj_range < 1e-12:
            continue

        for i in range(1, n - 1):
            distances[sorted_idx[i]] += (
                (obj_vals[sorted_idx[i + 1]] - obj_vals[sorted_idx[i - 1]])
                / obj_range
            )

    return distances


def nsga2_select(population: List[tuple], objectives: np.ndarray,
                 n_select: int) -> List[int]:
    """
    NSGA-II selection: prefer lower rank, then higher crowding distance.
    Returns indices of selected individuals.
    """
    fronts = fast_non_dominated_sort(objectives)
    selected = []

    for front in fronts:
        if len(selected) + len(front) <= n_select:
            selected.extend(front)
        else:
            # Need partial front: select by crowding distance
            cd = crowding_distance(objectives, front)
            sorted_by_cd = sorted(range(len(front)), key=lambda i: -cd[i])
            remaining = n_select - len(selected)
            for i in sorted_by_cd[:remaining]:
                selected.append(front[i])
            break

    return selected


# ═══════════════════════════════════════════════════════════════════════════════
# OBJECTIVE COMPUTATION
# ═══════════════════════════════════════════════════════════════════════════════

def compute_objectives_surrogate(population: List[tuple],
                                  model: CatalystSurrogate,
                                  device: str = 'cuda:0') -> np.ndarray:
    """
    Compute 4 objectives using the surrogate model.
    All objectives are MINIMIZED (negate what should be maximized).
    
    Returns: (N, 4) array
    """
    X = encode_population(population)
    preds = predict_batch(model, X, device=device)

    # Objective 1: E_act (minimize)
    obj1 = preds['E_act']

    # Objective 2: Coking resistance (maximize → negate for minimization)
    obj2 = -preds['coking_index']

    # Objective 3: Stability (minimize segregation energy — more negative = more stable)
    obj3 = preds['segregation_energy']  # already: negative = good

    # Objective 4: Material cost (minimize)
    cost_penalties = np.array([
        abundance_cost_penalty(_extract_elements_from_genome(g))
        for g in population
    ])
    obj4 = -cost_penalties  # abundance_cost_penalty returns [-2, 0]; negate so abundant → small

    objectives = np.column_stack([obj1, obj2, obj3, obj4])

    # Penalty for invalid candidates (push them to worst-case objectives)
    invalid_mask = preds['valid_prob'] < 0.5
    objectives[invalid_mask] = [5.0, 0.0, 1.0, 3.0]  # worst-case values

    return objectives


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
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN GENETIC ALGORITHM
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class GAConfig:
    """Configuration for the genetic algorithm."""
    pop_size: int = 500
    n_generations: int = 200
    mace_eval_interval: int = 50        # Full MACE eval every N generations
    mace_eval_top_k: int = 100          # Top-K from Pareto front for MACE
    surrogate_retrain_interval: int = 50
    mutation_rate: float = 0.3
    crossover_rate: float = 0.7
    tournament_size: int = 5
    initial_mace_samples: int = 200     # Initial MACE samples for surrogate training
    device: str = 'cuda:0'
    seed: int = 42


def tournament_select(population: List[tuple], objectives: np.ndarray,
                      tournament_size: int = 5) -> int:
    """Tournament selection: pick the best from a random subset."""
    candidates = random.sample(range(len(population)), min(tournament_size, len(population)))
    # Prefer lower Pareto rank → select by first-objective dominance
    best = candidates[0]
    for c in candidates[1:]:
        if dominates(objectives[c], objectives[best]):
            best = c
    return best


def run_genetic_algorithm(config: GAConfig = GAConfig(),
                          existing_db: Optional[pd.DataFrame] = None) -> Tuple[List[tuple], pd.DataFrame]:
    """
    Execute the NSGA-II genetic algorithm for catalyst discovery.
    
    Args:
        config: GA configuration
        existing_db: Existing MACE screening results to seed surrogate
        
    Returns:
        (pareto_front_genomes, full_screening_database)
    """
    random.seed(config.seed)
    np.random.seed(config.seed)

    print_banner("NSGA-II MULTI-OBJECTIVE GENETIC ALGORITHM")
    logger.info(f"Population: {config.pop_size}, Generations: {config.n_generations}")

    # ── Phase A: Initial MACE Screening for Surrogate Training ──────────────
    if existing_db is not None and len(existing_db) > 50:
        logger.info(f"Using existing database with {len(existing_db)} entries for surrogate seed")
        all_mace_results = existing_db
    else:
        logger.info(f"Generating {config.initial_mace_samples} initial MACE samples...")
        initial_pop = generate_population(config.initial_mace_samples)
        from pipeline.mace_screener import run_screening
        all_mace_results = run_screening(initial_pop, db_filename="ga_initial_screening.csv",
                                         workers_per_gpu=2)

    # ── Phase B: Train Initial Surrogate ────────────────────────────────────
    model = _train_surrogate_from_db(all_mace_results, config.device)

    # ── Phase C: Evolutionary Loop ──────────────────────────────────────────
    population = generate_population(config.pop_size)
    best_e_act_history = []
    pareto_front_history = []

    for gen in range(config.n_generations):
        t_gen = time.time()

        # Evaluate population with surrogate
        objectives = compute_objectives_surrogate(population, model, config.device)

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

        # ── Class-Diversity Enforcement ─────────────────────────────────────
        # Prevent any single class from dominating the population.
        # Guarantee at least 5% of population from each of the 10 classes.
        combined = parents + offspring
        min_per_class = max(2, config.pop_size // 20)  # 5% of pop per class
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
            # Replace the weakest members (tail of combined) with fresh diverse candidates
            combined = combined[:len(combined) - len(diversity_injections)] + diversity_injections

        combined_obj = compute_objectives_surrogate(combined, model, config.device)
        final_idx = nsga2_select(combined, combined_obj, config.pop_size)
        population = [combined[i] for i in final_idx]
        final_obj = combined_obj[final_idx]

        # Track best activation barrier in population
        best_e_act = final_obj[:, 0].min()
        best_e_act_history.append(best_e_act)

        # Pareto front size
        fronts = fast_non_dominated_sort(final_obj)
        pareto_size = len(fronts[0]) if fronts else 0
        pareto_front_history.append(pareto_size)

        # ── Periodic MACE Validation ────────────────────────────────────────
        if (gen + 1) % config.mace_eval_interval == 0:
            logger.info(f"  Gen {gen+1}: Running MACE validation on top {config.mace_eval_top_k}...")
            pareto_genomes = [population[i] for i in fronts[0][:config.mace_eval_top_k]]
            from pipeline.mace_screener import run_screening
            mace_df = run_screening(
                pareto_genomes,
                db_filename=f"ga_mace_gen{gen+1}.csv",
                workers_per_gpu=2
            )
            all_mace_results = pd.concat([all_mace_results, mace_df], ignore_index=True)

        # ── Periodic Surrogate Retraining ───────────────────────────────────
        if (gen + 1) % config.surrogate_retrain_interval == 0:
            logger.info(f"  Gen {gen+1}: Retraining surrogate on {len(all_mace_results)} samples...")
            model = _train_surrogate_from_db(all_mace_results, config.device)

        # Logging
        if (gen + 1) % 10 == 0 or gen == 0:
            elapsed = time.time() - t_gen
            logger.info(
                f"  Gen {gen+1}/{config.n_generations}: "
                f"best_E_act={best_e_act:.4f} eV, "
                f"pareto_size={pareto_size}, "
                f"pop_diversity={len(set(str(g) for g in population))}, "
                f"({elapsed:.1f}s)"
            )

    # ── Final Pareto Front ──────────────────────────────────────────────────
    final_objectives = compute_objectives_surrogate(population, model, config.device)
    fronts = fast_non_dominated_sort(final_objectives)
    pareto_genomes = [population[i] for i in fronts[0]]

    logger.info(f"\n  GA Complete. Final Pareto front: {len(pareto_genomes)} candidates")
    logger.info(f"  Total MACE evaluations: {len(all_mace_results)}")

    # Save final results
    save_screening_db(all_mace_results, "ga_full_database.csv")

    return pareto_genomes, all_mace_results


def _train_surrogate_from_db(df: pd.DataFrame, device: str) -> CatalystSurrogate:
    """Train surrogate from a MACE screening database DataFrame."""
    valid_df = df.dropna(subset=['E_act', 'coking_index', 'segregation_energy', 'dE_split'])

    if len(valid_df) < 10:
        logger.warning(f"Only {len(valid_df)} valid samples. Surrogate quality may be low.")
        # Fill with defaults for training
        if len(valid_df) == 0:
            # Return untrained model
            model = CatalystSurrogate().to(device)
            return model

    # Parse genomes from string representation
    genomes = []
    for _, row in df.iterrows():
        try:
            g = eval(row['genome'])
            genomes.append(g)
        except Exception:
            genomes.append(('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 0, 0))

    X = encode_population(genomes)

    y_valid = df['valid'].astype(float).values
    y_de_split = df.get('dE_split', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_coking = df.get('coking_index', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_seg = df.get('segregation_energy', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_e_act = df.get('E_act', pd.Series(np.ones(len(df)))).fillna(1.0).values

    model = train_surrogate(
        X, y_valid, y_de_split, y_coking, y_seg, y_e_act,
        epochs=30, batch_size=min(2048, len(X)),
        device=device
    )
    return model


if __name__ == '__main__':
    config = GAConfig(
        pop_size=100,
        n_generations=50,
        initial_mace_samples=50,
        mace_eval_interval=25,
        surrogate_retrain_interval=25,
    )
    pareto, db = run_genetic_algorithm(config)
    print(f"\nPareto front: {len(pareto)} candidates")
    for g in pareto[:5]:
        print(f"  {g}")
