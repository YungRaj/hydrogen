#!/usr/bin/env python3
"""Opt-in real multi-GPU affinity smoke test for eSen workers."""

import os
import subprocess
import time

from pipeline.search.indexed_space import deterministic_tree_probes
from pipeline.screening.surface_screener import run_screening


def main():
    engine = os.environ.get('HYDROGEN_SCREENING_ENGINE', 'batched')
    workers_per_gpu = int(os.environ.get('HYDROGEN_WORKERS_PER_GPU', '1'))
    application = os.environ.get('HYDROGEN_SCREENING_APPLICATION', 'pyrolysis')
    expected = int(subprocess.check_output(
        ['nvidia-smi', '--query-gpu=index', '--format=csv,noheader'],
        text=True).count('\n'))
    if expected < 1:
        raise RuntimeError('no NVIDIA GPUs detected')
    started = time.perf_counter()
    probes = deterministic_tree_probes(max(2 * expected, 14))
    if application == 'pyrolysis':
        frame = run_screening(
            probes, db_filename=f'gpu_affinity_contract_{engine}.csv',
            workers_per_gpu=workers_per_gpu, engine=engine)
    elif application == 'orr':
        from pipeline.screening.fc_screener import run_orr_screening
        frame = run_orr_screening(
            probes, db_filename=f'gpu_affinity_contract_{engine}_orr.csv',
            workers_per_gpu=workers_per_gpu, engine=engine)
    else:
        raise ValueError("HYDROGEN_SCREENING_APPLICATION must be 'pyrolysis' or 'orr'")
    used = {int(value) for value in frame['gpu_id'].dropna().tolist()}
    assert used == set(range(expected)), (used, expected)
    assert frame['screening_protocol'].notna().all()
    print({'application': application, 'engine': engine,
           'workers_per_gpu': workers_per_gpu,
           'visible_gpus': expected, 'used_gpu_ids': sorted(used),
           'candidates': len(frame),
           'elapsed_s': time.perf_counter() - started})


if __name__ == '__main__':
    main()
