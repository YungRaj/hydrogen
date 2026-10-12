"""Small, reproducible policy comparison on stored screening evidence."""

from __future__ import annotations

import argparse
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.screening.genetic_optimizer import _thermal_full_simulator


def compare_policies(source: Path, *, cohort_size: int = 10, draws: int = 1000) -> dict:
    """Compare random, eSen, and thermal policies on one stored pool."""
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


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cohort-size', type=int, default=10)
    parser.add_argument('--draws', type=int, default=1000)
    args = parser.parse_args()
    result = compare_policies(args.source, cohort_size=args.cohort_size, draws=args.draws)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
