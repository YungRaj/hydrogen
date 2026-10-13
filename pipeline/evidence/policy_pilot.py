"""Small, reproducible policy comparison on stored screening evidence."""

from __future__ import annotations

import argparse
import ast
import json
import random
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from pipeline.screening.genetic_optimizer import _thermal_full_simulator
from pipeline.search.design_space import (
    generate_random_genome,
    generate_supported_nanoparticle_genome,
)
from pipeline.screening.surface_screener import run_screening


def compare_policies(
    source: Path, *, cohort_size: int = 10, draws: int = 1000
) -> dict[str, Any]:
    """Compare policies on one stored pool.

    Args:
        source: Screening CSV containing candidate evidence.
        cohort_size: Number of candidates selected by each policy.
        draws: Number of random cohorts used for the baseline distribution.

    Returns:
        Aggregate policy means and random-baseline interval.
    """
    frame = pd.read_csv(source).drop_duplicates('genome')
    genomes = [ast.literal_eval(value) for value in frame['genome']]
    thermal = _thermal_full_simulator(genomes, frame)
    frame = frame.join(thermal.set_index('genome'), on='genome')
    frame = frame[frame['full_simulation_status'].eq('complete')].copy()
    size = min(int(cohort_size), len(frame))
    if size <= 0:
        raise ValueError('no complete thermal results available')
    random_scores = []
    rng = np.random.default_rng(20261012)
    for seed in rng.integers(0, 2**31, int(draws)):
        random_scores.append(
            frame.sample(size, random_state=int(seed))['full_simulation_score'].mean()
        )
    random_scores = np.asarray(random_scores, dtype=float)
    cohorts = {
        'random': float(random_scores.mean()),
        'eSen': float(frame.nsmallest(size, 'E_act')['full_simulation_score'].mean()),
        'thermal': float(
            frame.nsmallest(size, 'full_simulation_score')['full_simulation_score'].mean()
        ),
        'uncertainty': float(
            frame.nsmallest(size, 'model_confidence')['full_simulation_score'].mean()
        ),
        'validity': float(
            frame.sort_values(['valid', 'pyrolysis_viable'], ascending=False)
            .head(size)['full_simulation_score'].mean()
        ),
    }
    return {
        'status': 'retrospective_pilot',
        'source': str(source),
        'pool_size': int(len(frame)),
        'cohort_size': size,
        'random_interval_95': [
            float(np.quantile(random_scores, 0.025)),
            float(np.quantile(random_scores, 0.975)),
        ],
        'mean_score_lower_is_better': cohorts,
        'thermal_lift_vs_random': cohorts['random'] - cohorts['thermal'],
        'esen_lift_vs_random': cohorts['random'] - cohorts['eSen'],
    }


def create_fresh_pool(
    output: Path, count: int, *, seed: int = 20261012,
    material_class: str | None = None,
    eligible_metals: tuple[str, ...] = (),
) -> Path:
    """Screen a deterministic random pool into a separate pilot artifact.

    Args:
        output: CSV path for the screened pool.
        count: Number of unique candidates to generate.
        seed: Initial deterministic generator seed.
        material_class: Optional material class restriction.
        eligible_metals: Optional supported-nanoparticle metal allowlist.

    Returns:
        The written output path.
    """
    if count <= 0:
        raise ValueError('fresh pool count must be positive')
    pool: list[tuple[Any, ...]] = []
    seen: set[str] = set()
    next_seed = int(seed)
    while len(pool) < count:
        random_state = random.getstate()
        random.seed(next_seed)
        try:
            candidate = (
                generate_supported_nanoparticle_genome()
                if material_class == 'SupportedNanoparticle'
                else generate_random_genome(material_class)
            )
        finally:
            random.setstate(random_state)
        if material_class == 'SupportedNanoparticle' and eligible_metals:
            candidate = (candidate[0], eligible_metals[len(pool) % len(eligible_metals)], *candidate[2:])
        next_seed += 1
        key = repr(candidate)
        if key not in seen:
            seen.add(key)
            pool.append(candidate)
    frame = run_screening(pool, db_filename=str(output), workers_per_gpu=1)
    output.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output, index=False)
    return output


def main() -> None:
    """Parse pilot arguments, run screening and write the comparison report."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cohort-size', type=int, default=10)
    parser.add_argument('--draws', type=int, default=1000)
    parser.add_argument('--fresh-count', type=int, default=0)
    parser.add_argument('--pool-seed', type=int, default=20261012)
    parser.add_argument('--material-class', default=None)
    parser.add_argument('--eligible-metals', nargs='*', default=())
    args = parser.parse_args()
    source = args.source
    if args.fresh_count:
        source = create_fresh_pool(
            source, args.fresh_count, seed=args.pool_seed,
            material_class=args.material_class,
            eligible_metals=tuple(args.eligible_metals),
        )
    result = compare_policies(source, cohort_size=args.cohort_size, draws=args.draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
