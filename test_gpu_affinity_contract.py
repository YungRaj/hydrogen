#!/usr/bin/env python3
"""Opt-in real multi-GPU affinity smoke test for eSen workers."""

import os
import subprocess
import time

from pipeline.search.indexed_space import deterministic_tree_probes
from pipeline.screening.surface_screener import run_screening


def main():
    expected = int(subprocess.check_output(
        ['nvidia-smi', '--query-gpu=index', '--format=csv,noheader'],
        text=True).count('\n'))
    if expected < 1:
        raise RuntimeError('no NVIDIA GPUs detected')
    started = time.perf_counter()
    frame = run_screening(
        deterministic_tree_probes(max(2 * expected, 14)),
        db_filename='gpu_affinity_contract.csv', workers_per_gpu=1)
    used = {int(value) for value in frame['gpu_id'].dropna().tolist()}
    assert used == set(range(expected)), (used, expected)
    assert frame['screening_protocol'].notna().all()
    print({'visible_gpus': expected, 'used_gpu_ids': sorted(used),
           'candidates': len(frame),
           'elapsed_s': time.perf_counter() - started})


if __name__ == '__main__':
    main()
