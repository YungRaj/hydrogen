"""Shared multi-GPU execution shell for atomistic screening applications."""

from __future__ import annotations

import multiprocessing as mp
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from pipeline.common.utils import print_banner, save_screening_db
from pipeline.screening.worker_supervisor import collect_results, emit, start_heartbeat


@dataclass(frozen=True)
class ScreeningRunSpec:
    """Application-specific labels around the common GPU execution protocol."""

    banner: str
    application: str
    manifest_path: Path
    output_subdir: str
    start_message: str
    completion_label: str
    progress_noun: str = 'candidates'


def execution_layout(device_count: int, workers_per_gpu: int,
                     engine: str) -> tuple[int, int]:
    """Return (process count, candidate threads/process) with validation."""
    if device_count < 1:
        raise RuntimeError('screening requires at least one visible CUDA GPU')
    if workers_per_gpu < 1:
        raise ValueError('workers_per_gpu must be positive')
    if engine not in ('batched', 'legacy'):
        raise ValueError("engine must be 'batched' or 'legacy'")
    processes_per_gpu = 1 if engine == 'batched' else workers_per_gpu
    candidate_threads = workers_per_gpu if engine == 'batched' else 1
    return device_count * processes_per_gpu, candidate_threads


def run_worker_loop(worker_id, task_queue, status_queue, stop_event,
                    candidate_threads: int, batched: bool, calculator,
                    evaluator: Callable, error_record: Callable,
                    batch_service=None, result_context: dict | None = None) -> None:
    """Run the shared leased-task and heartbeat loop inside one GPU process."""
    import queue
    import threading
    from concurrent.futures import ThreadPoolExecutor

    heartbeat_stop = threading.Event()
    heartbeat = start_heartbeat(status_queue, worker_id, heartbeat_stop)
    emit(status_queue, 'ready', worker_id)

    def consume():
        thread_calculator = (batch_service.calculator_proxy()
                             if batch_service is not None else calculator)
        while True:
            try:
                item = task_queue.get(timeout=1.0)
            except queue.Empty:
                if stop_event.is_set():
                    break
                continue
            index, genome = item
            emit(status_queue, 'started', worker_id, index)
            try:
                result = evaluator(genome, thread_calculator)
                result['worker_id'] = worker_id
                result.update(result_context or {})
                emit(status_queue, 'result', worker_id, (index, result))
            except Exception as exc:
                emit(status_queue, 'result', worker_id,
                     (index, error_record(genome, exc)))

    try:
        if batched:
            with ThreadPoolExecutor(max_workers=candidate_threads) as executor:
                futures = [executor.submit(consume)
                           for _ in range(candidate_threads)]
                for future in futures:
                    future.result()
        else:
            consume()
    finally:
        if batch_service is not None:
            batch_service.close()
        heartbeat_stop.set()
        heartbeat.join(timeout=2)


def run_gpu_screening(
    genomes: Sequence[tuple],
    db_filename: str,
    workers_per_gpu: int,
    engine: str,
    worker_target: Callable,
    logger,
    spec: ScreeningRunSpec,
):
    """Execute a screening evaluator with deterministic, supervised GPU work.

    The worker callable retains ownership of model initialization and scientific
    evaluation. This function owns only invariant execution behavior: topology,
    queues, leases, health supervision, ordering, and persistence.
    """
    import pandas as pd
    import torch

    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ.setdefault(name, '1')
    mp.set_start_method('spawn', force=True)

    print_banner(spec.banner)
    logger.info(spec.start_message.format(count=len(genomes)))
    device_count = torch.cuda.device_count()
    if device_count < 1:
        raise RuntimeError(
            f'{spec.completion_label} requires at least one visible CUDA GPU')
    num_workers, candidate_threads = execution_layout(
        device_count, workers_per_gpu, engine)
    gpu_uuids = [f"GPU-{torch.cuda.get_device_properties(i).uuid}"
                 for i in range(device_count)]
    logger.info(
        f"Using {device_count} GPU(s), engine={engine}, "
        f"{num_workers} model process(es), "
        f"{candidate_threads} candidate thread(s)/process")

    task_queue, status_queue, stop_event = mp.Queue(), mp.Queue(), mp.Event()
    for index, genome in enumerate(genomes):
        task_queue.put((index, genome))

    def spawn_worker(worker_id):
        gpu_id = worker_id % device_count
        process = mp.Process(
            target=worker_target,
            args=(worker_id, gpu_id, gpu_uuids[gpu_id], task_queue,
                  status_queue, stop_event, candidate_threads,
                  engine == 'batched'))
        process.start()
        return process

    workers = {worker_id: spawn_worker(worker_id)
               for worker_id in range(num_workers)}
    started = time.time()

    def report_progress(completed, results):
        if completed % 50 == 0 or completed == len(genomes):
            elapsed = time.time() - started
            rate = completed / max(elapsed, 1e-12)
            valid = sum(1 for result in results if result.get('valid', False))
            logger.info(
                f"Progress: {completed}/{len(genomes)} "
                f"({rate:.1f} {spec.progress_noun}/sec, {valid} valid, "
                f"{elapsed:.0f}s elapsed)")

    results = collect_results(
        status_queue, task_queue, stop_event, workers, spawn_worker, genomes,
        spec.application, spec.manifest_path, progress=report_progress)
    for process in workers.values():
        process.join(timeout=30)

    frame = pd.DataFrame(results)
    path = save_screening_db(
        frame, db_filename, subdir=spec.output_subdir)
    logger.info(
        f"{spec.completion_label} complete. {len(frame)} results saved to {path}")
    return frame
