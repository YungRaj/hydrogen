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
import ast
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
    ALL_MATERIAL_CLASSES, FEATURE_DIM, generate_hierarchical_htvs_pool,
)
from pipeline.surrogate_model import CatalystSurrogate, train_surrogate, predict_batch, SurrogateEnsemble, train_ensemble, predict_ensemble

logger = setup_logger('genetic_optimizer', 'screening/genetic_optimizer.log')


# ═══════════════════════════════════════════════════════════════════════════════
# NSGA-II IMPLEMENTATION
# ═══════════════════════════════════════════════════════════════════════════════

def dominates(obj_a: np.ndarray, obj_b: np.ndarray) -> bool:
    """Return True if solution a dominates solution b (all ≤, at least one <)."""
    return np.all(obj_a <= obj_b) and np.any(obj_a < obj_b)


def fast_non_dominated_sort(objectives: np.ndarray) -> List[List[int]]:
    """
    NSGA-II fast non-dominated sort using vectorized Pareto front extraction.
    """
    n = len(objectives)
    remaining_indices = np.arange(n)
    fronts = []

    while len(remaining_indices) > 0:
        sub_objs = objectives[remaining_indices]
        is_efficient = np.ones(len(sub_objs), dtype=bool)
        for i in range(len(sub_objs)):
            if is_efficient[i]:
                dominated = np.all(sub_objs[i] <= sub_objs, axis=1) & np.any(sub_objs[i] < sub_objs, axis=1)
                is_efficient[dominated] = False

        front_sub_idx = np.where(is_efficient)[0]
        front_global_idx = remaining_indices[front_sub_idx].tolist()
        fronts.append(front_global_idx)

        remaining_indices = np.delete(remaining_indices, front_sub_idx)

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
                                  model: object,
                                  device: str = 'cuda:0') -> np.ndarray:
    """
    Compute 4 objectives using the surrogate model or surrogate ensemble.
    All objectives are MINIMIZED (negate what should be maximized).
    
    If model is a SurrogateEnsemble, uses LCB (for E_act and segregation energy)
    and UCB (for coking resistance index) with kappa = 1.0 to drive active learning.
    
    E_act is scaled by OOD confidence penalty to discount predictions
    from material classes outside the eSen-SM training distribution.
    
    Returns: (N, 4) array
    """
    from pipeline.ood_detector import compute_model_confidence, confidence_penalty

    X = encode_population(population)
    
    if isinstance(model, SurrogateEnsemble):
        preds = predict_ensemble(model, X, device=device)
        kappa = 1.0
        e_act_pred = preds['E_act'] - kappa * preds['E_act_std']
        coking_pred = preds['coking_index'] + kappa * preds['coking_index_std']
        seg_pred = preds['segregation_energy'] - kappa * preds['segregation_energy_std']
        valid_prob = preds['valid_prob']
    else:
        preds = predict_batch(model, X, device=device)
        e_act_pred = preds['E_act']
        coking_pred = preds['coking_index']
        seg_pred = preds['segregation_energy']
        valid_prob = preds['valid_prob']

    # Objective 1: E_act (minimize) × confidence penalty
    obj1 = e_act_pred.copy()

    # Objective 2: Coking resistance (maximize → negate for minimization)
    coking = coking_pred.copy()
    py_mode = os.environ.get('PYROLYSIS_MODE', 'ntec')
    if py_mode == 'ntec':
        for i, g in enumerate(population):
            if any(e in {'Ga', 'In', 'Sn', 'Bi'} for e in _extract_elements_from_genome(g)):
                coking[i] += 3.0
    obj2 = -coking

    # Objective 3: Stability (minimize segregation energy — more negative = more stable)
    obj3 = seg_pred.copy()  # already: negative = good

    # Objective 4: Material cost (minimize)
    cost_penalties = np.array([
        abundance_cost_penalty(_extract_elements_from_genome(g))
        for g in population
    ])
    obj4 = -cost_penalties  # abundance_cost_penalty returns [-2, 0]; negate so abundant → small

    # Apply OOD confidence penalty to E_act
    for i, g in enumerate(population):
        elements = _extract_elements_from_genome(g)
        conf = compute_model_confidence(g, elements)
        obj1[i] += confidence_penalty(conf)  # additive OOD penalty

    objectives = np.column_stack([obj1, obj2, obj3, obj4])

    # Penalty for invalid candidates (push them to worst-case objectives)
    invalid_mask = valid_prob < 0.5
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
    elif mat_class == 'Spinel':
        elements.extend([genome[1], genome[2]])
        if genome[3] != 'None': elements.append(genome[3])
    elif mat_class == 'MXene':
        elements.append(genome[1])
        if genome[5] != 'None': elements.append(genome[5])
    elif mat_class == 'SAA':
        elements.extend([genome[1], genome[2]])
    elif mat_class == 'MetalFreeCarbon':
        pass  # no metals
    return [e for e in elements if e != 'None']


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN GENETIC ALGORITHM
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class GAConfig:
    """Configuration for the genetic algorithm."""
    pop_size: int = 500
    n_generations: int = 200
    fairchem_eval_interval: int = 50        # Full Fairchem eval every N generations
    fairchem_eval_top_k: int = 100          # Top-K from Pareto front for Fairchem
    surrogate_retrain_interval: int = 50
    mutation_rate: float = 0.3
    crossover_rate: float = 0.7
    tournament_size: int = 5
    initial_fairchem_samples: int = 200     # Initial Fairchem samples for surrogate training
    explore_interval: int = 3               # Run exploration shots every N Fairchem intervals
    explore_per_class: int = 5              # Random GNN evaluations per class during exploration
    n_models: int = 3                       # Number of models in surrogate ensemble
    htvs_pool_size: int = 20000             # Initial high-throughput virtual screening pool size
    reinjection_interval: int = 20          # Periodically reinject top candidates from new pool
    device: str = 'cuda:0'
    seed: int = 42





def tournament_select(population: List[tuple], objectives: np.ndarray,
                      tournament_size: int = 5) -> int:
    """Tournament selection: pick the best from a random subset."""
    candidates = random.sample(range(len(population)), min(tournament_size, len(population)))
    best = candidates[0]
    for c in candidates[1:]:
        if dominates(objectives[c], objectives[best]):
            best = c
    return best


def run_genetic_algorithm(config: GAConfig = GAConfig(),
                          existing_db: Optional[pd.DataFrame] = None) -> Tuple[List[tuple], pd.DataFrame]:
    """
    Execute the NSGA-II genetic algorithm for catalyst discovery.
    """
    random.seed(config.seed)
    np.random.seed(config.seed)

    print_banner("NSGA-II MULTI-OBJECTIVE GENETIC ALGORITHM")
    logger.info(f"Population: {config.pop_size}, Generations: {config.n_generations}")

    # ── Phase A: Initial Fairchem Screening for Surrogate Training ──────────
    if existing_db is not None and len(existing_db) > 50:
        logger.info(f"Using existing database with {len(existing_db)} entries for surrogate seed")
        all_fairchem_results = existing_db
    else:
        logger.info(f"Generating {config.initial_fairchem_samples} initial Fairchem samples...")
        initial_pop = generate_population(config.initial_fairchem_samples)
        from pipeline.surface_screener import run_screening
        all_fairchem_results = run_screening(initial_pop, db_filename="ga_initial_screening.csv",
                                             workers_per_gpu=2)

    # ── Phase B: Train Initial Surrogate Ensemble ───────────────────────────
    model = _train_ensemble_from_db(all_fairchem_results, config.device, n_models=config.n_models)

    # ── Phase C: Evolutionary Loop ──────────────────────────────────────────
    logger.info(f"Generating HTVS pool of {config.htvs_pool_size} candidates for initial seeding...")
    htvs_pool = generate_hierarchical_htvs_pool(
        config.htvs_pool_size,
        scorer=lambda pop: compute_objectives_surrogate(pop, model, config.device)[:, 0]
    )
    htvs_obj = compute_objectives_surrogate(htvs_pool, model, config.device)

    logger.info(f"Selecting top {config.pop_size} Pareto-optimal seeds using acquisition LCB/UCB values...")
    seed_idx = nsga2_select(htvs_pool, htvs_obj, config.pop_size)
    population = [htvs_pool[i] for i in seed_idx]

    best_e_act_history = []
    pareto_front_history = []
    fairchem_round = 0

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

        # ── Periodic Fairchem Validation ────────────────────────────────────
        if (gen + 1) % config.fairchem_eval_interval == 0:
            fairchem_round += 1
            logger.info(f"  Gen {gen+1}: Running Fairchem validation on top {config.fairchem_eval_top_k}...")
            pareto_genomes = [population[i] for i in fronts[0][:config.fairchem_eval_top_k]]
            from pipeline.surface_screener import run_screening
            fairchem_df = run_screening(
                pareto_genomes,
                db_filename=f"ga_fairchem_gen{gen+1}.csv",
                workers_per_gpu=2
            )
            all_fairchem_results = pd.concat([all_fairchem_results, fairchem_df], ignore_index=True)

            # ── Exploration Shots: probe EVERY class with real GNN ───────────
            if fairchem_round % config.explore_interval == 0:
                explore_genomes = []
                for cls in ALL_MATERIAL_CLASSES:
                    explore_genomes.extend(
                        generate_population(config.explore_per_class, material_class=cls)
                    )
                n_explore = len(explore_genomes)
                logger.info(
                    f"  Gen {gen+1}: EXPLORATION — evaluating {n_explore} random "
                    f"candidates across all 14 classes with real GNN..."
                )
                explore_df = run_screening(
                    explore_genomes,
                    db_filename=f"ga_explore_gen{gen+1}.csv",
                    workers_per_gpu=2
                )
                all_fairchem_results = pd.concat([all_fairchem_results, explore_df], ignore_index=True)

                # Inject any surprisingly good exploration candidates into population
                if 'E_act' in explore_df.columns:
                    good_explores = explore_df[
                        (explore_df['valid'] == True) &
                        (explore_df['E_act'] < explore_df['E_act'].quantile(0.3))
                    ]
                    if len(good_explores) > 0:
                        logger.info(
                            f"    Found {len(good_explores)} promising exploration "
                            f"candidates — injecting into population"
                        )
                        for _, row in good_explores.iterrows():
                            try:
                                g = ast.literal_eval(row['genome'])
                                replace_idx = random.randint(
                                    len(fronts[0]), len(population) - 1
                                ) if len(fronts[0]) < len(population) else random.randint(
                                    0, len(population) - 1
                                )
                                population[replace_idx] = g
                            except Exception:
                                pass

        # ── Periodic HTVS Global Reinjection ────────────────────────────────
        if (gen + 1) % config.reinjection_interval == 0:
            logger.info(f"  Gen {gen+1}: Global HTVS — screening 10,000 fresh candidates...")
            reinject_pool = generate_hierarchical_htvs_pool(
                10000,
                scorer=lambda pop: compute_objectives_surrogate(pop, model, config.device)[:, 0]
            )
            reinject_obj = compute_objectives_surrogate(reinject_pool, model, config.device)

            n_inject = max(10, config.pop_size // 10)
            inject_idx = nsga2_select(reinject_pool, reinject_obj, n_inject)
            inject_genomes = [reinject_pool[i] for i in inject_idx]

            # Merge with current population and select the top pop_size
            combined_pop = population + inject_genomes
            combined_obj = compute_objectives_surrogate(combined_pop, model, config.device)
            keep_idx = nsga2_select(combined_pop, combined_obj, config.pop_size)
            population = [combined_pop[i] for i in keep_idx]
            final_obj = combined_obj[keep_idx]

        # ── Periodic Surrogate Retraining ───────────────────────────────────
        if (gen + 1) % config.surrogate_retrain_interval == 0:
            logger.info(f"  Gen {gen+1}: Retraining surrogate ensemble on {len(all_fairchem_results)} samples...")
            model = _train_ensemble_from_db(all_fairchem_results, config.device, n_models=config.n_models)

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
    logger.info(f"  Total Fairchem evaluations: {len(all_fairchem_results)}")

    # Save final results
    save_screening_db(all_fairchem_results, "ga_full_database.csv")

    return pareto_genomes, all_fairchem_results


def _train_ensemble_from_db(df: pd.DataFrame, device: str, n_models: int = 3) -> SurrogateEnsemble:
    """Train surrogate ensemble from a Fairchem screening database DataFrame."""
    valid_df = df.dropna(subset=['E_act', 'coking_index', 'segregation_energy', 'dE_split'])

    if len(valid_df) < 10:
        logger.warning(f"Only {len(valid_df)} valid samples. Ensemble quality may be low.")
        # Return untrained ensemble
        ensemble = SurrogateEnsemble(n_models=n_models).to(device)
        return ensemble

    # Parse genomes from string representation
    genomes = []
    for _, row in df.iterrows():
        try:
            g = ast.literal_eval(row['genome'])
            genomes.append(g)
        except Exception:
            genomes.append(('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 0, 0))

    X = encode_population(genomes)

    y_valid = df['valid'].astype(float).values
    y_de_split = df.get('dE_split', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_coking = df.get('coking_index', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_seg = df.get('segregation_energy', pd.Series(np.zeros(len(df)))).fillna(0.0).values
    y_e_act = df.get('E_act', pd.Series(np.ones(len(df)))).fillna(1.0).values

    ensemble = train_ensemble(
        X, y_valid, y_de_split, y_coking, y_seg, y_e_act,
        n_models=n_models, epochs=30, batch_size=min(2048, len(X)),
        device=device
    )
    return ensemble


if __name__ == '__main__':
    config = GAConfig(
        pop_size=100,
        n_generations=50,
        initial_fairchem_samples=50,
        fairchem_eval_interval=25,
        surrogate_retrain_interval=25,
    )
    pareto, db = run_genetic_algorithm(config)
    print(f"\nPareto front: {len(pareto)} candidates")
    for g in pareto[:5]:
        print(f"  {g}")
