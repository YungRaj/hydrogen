#!/usr/bin/env python3
"""Benchmark eSEN screening throughput versus candidate threads per GPU.

Each setting screens the same deterministic genomes into an isolated results
directory (no cache hits), while a sampler records GPU utilization and memory.

Usage:
    python scripts/benchmark_screening.py --candidates 64 --workers 1 4 8 16
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import tempfile
import threading
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]


def _sample_gpu(stop: threading.Event, samples: list) -> None:
    while not stop.is_set():
        output = subprocess.run(
            ['nvidia-smi', '--query-gpu=utilization.gpu,memory.used',
             '--format=csv,noheader,nounits'],
            capture_output=True, text=True,
        ).stdout.strip().splitlines()
        for line in output:
            util, memory = (float(x) for x in line.split(','))
            samples.append((util, memory))
        stop.wait(1.0)


def _run_setting(workers: int, candidates: int, seed: int) -> dict:
    """Screen in a fresh interpreter so HYDROGEN_RESULTS_DIR applies at import."""
    with tempfile.TemporaryDirectory(prefix='bench_') as results:
        code = f"""
import json, random, time
from pipeline.search.design_space import generate_random_genome
from pipeline.screening.surface_screener import run_screening
random.seed({seed})
genomes = [generate_random_genome() for _ in range({candidates})]
start = time.time()
frame = run_screening(genomes, db_filename='bench.csv', workers_per_gpu={workers})
print('BENCH', json.dumps({{'elapsed_s': time.time() - start,
                            'rows': len(frame),
                            'valid': int(frame['valid'].sum())}}))
"""
        stop, samples = threading.Event(), []
        sampler = threading.Thread(target=_sample_gpu, args=(stop, samples), daemon=True)
        sampler.start()
        proc = subprocess.run(
            [sys.executable, '-c', code], cwd=REPO, capture_output=True, text=True,
            env={**os.environ, 'HYDROGEN_RESULTS_DIR': results},
        )
        stop.set()
        sampler.join()
    line = next((x for x in proc.stdout.splitlines() if x.startswith('BENCH ')), None)
    if line is None:
        return {'workers_per_gpu': workers, 'error': proc.stderr[-1500:]}
    result = json.loads(line[6:])
    busy = [u for u, _ in samples]
    return {
        'workers_per_gpu': workers,
        **result,
        'candidates_per_hour': 3600.0 * result['rows'] / result['elapsed_s'],
        'mean_gpu_util_pct': sum(busy) / len(busy) if busy else None,
        'peak_vram_MiB': max((m for _, m in samples), default=None),
    }


def main() -> None:
    """Run the benchmark grid and print one JSON line per setting."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidates', type=int, default=64)
    parser.add_argument('--workers', type=int, nargs='+', default=[1, 4, 8, 16])
    parser.add_argument('--seed', type=int, default=20261014)
    parser.add_argument('--output', type=Path, default=None)
    args = parser.parse_args()
    rows = []
    for workers in args.workers:
        row = _run_setting(workers, args.candidates, args.seed)
        rows.append(row)
        print(json.dumps(row), flush=True)
    if args.output:
        args.output.write_text(json.dumps(rows, indent=2) + '\n')


if __name__ == '__main__':
    main()
