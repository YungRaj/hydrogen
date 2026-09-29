"""Shared multi-GPU execution shell for atomistic screening applications."""

from __future__ import annotations

import multiprocessing as mp
import os
import time
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Sequence

from pipeline.utils import print_banner, save_screening_db
from pipeline.screening.worker_supervisor import collect_results, emit, start_heartbeat


def _record_cache_health(path: Path, requested: int, hits: int,
                         computed: int) -> None:
    """Add cache accounting without discarding worker supervision evidence."""
    if not path.exists():
        return
    payload = json.loads(path.read_text())
    payload['completed'] = requested
    payload['expected'] = requested
    payload['cache_hits'] = hits
    payload['computed'] = computed
    payload['deduplicated_reuses'] = requested - hits - computed
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(payload, indent=2) + '\n')
    temporary.replace(path)


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
    protocol_id: str | None = None
    cache_path: Path | None = None


def execution_layout(device_count: int, workers_per_gpu: int,
                     engine: str) -> tuple[int, int]:
    """Return (process count, candidate threads/process) with validation.

    Args:
        device_count: Number of device count to use.
        workers_per_gpu: Workers per gpu used by this operation.
        engine: Engine used by this operation.

    Returns:
        Ordered tuple of computed values.
    """
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
    """Run the shared leased-task and heartbeat loop inside one GPU process.

    Args:
        worker_id: Worker id used by this operation.
        task_queue: Task queue used by this operation.
        status_queue: Status queue used by this operation.
        stop_event: Stop event used by this operation.
        candidate_threads: Candidate threads used by this operation.
        batched: Whether to enable batched.
        calculator: Atomic calculator used for the evaluation.
        evaluator: Injected callable used to perform evaluator.
        error_record: Injected callable used to perform error record.
        batch_service: Batch service used by this operation.
        result_context: Result context used by this operation.
    """
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

    Args:
        genomes: Sequence of encoded catalyst candidates.
        db_filename: Db filename used by this operation.
        workers_per_gpu: Workers per gpu used by this operation.
        engine: Engine used by this operation.
        worker_target: Injected callable used to perform worker target.
        logger: Logger used by this operation.
        spec: Spec used by this operation.

    Returns:
        Computed result described above.
    """
    import pandas as pd

    requested_genomes = list(genomes)
    cached: dict[int, dict] = {}
    implementation = None
    if spec.cache_path is not None and spec.protocol_id is not None:
        from pipeline.screening.result_cache import (
            implementation_digest, load_cached_results)
        implementation = implementation_digest(worker_target)
        cached = load_cached_results(
            spec.cache_path, spec.application, spec.protocol_id,
            implementation, requested_genomes)

    # Evaluate each uncached canonical candidate once, even if the caller
    # supplied it more than once. Results are expanded back to input order.
    from pipeline.search.discovery import candidate_id
    misses, miss_ids = [], set()
    for index, genome in enumerate(requested_genomes):
        if index in cached:
            continue
        cid = candidate_id(genome)
        if cid not in miss_ids:
            misses.append(genome)
            miss_ids.add(cid)

    if not requested_genomes:
        frame = pd.DataFrame()
        path = save_screening_db(frame, db_filename, subdir=spec.output_subdir)
        logger.info(f'{spec.completion_label}: no candidates; saved empty result to {path}')
        return frame

    if not misses:
        results = []
        for index in range(len(requested_genomes)):
            result = dict(cached[index])
            result['genome'] = str(requested_genomes[index])
            result['material_class'] = requested_genomes[index][0]
            result['cache_hit'] = True
            results.append(result)
        frame = pd.DataFrame(results)
        path = save_screening_db(frame, db_filename, subdir=spec.output_subdir)
        from pipeline.screening.worker_supervisor import write_health_manifest
        write_health_manifest(
            spec.manifest_path, spec.application, 'complete', {}, len(results),
            len(results), [{'kind': 'cache_complete', 'worker_id': -1,
                            'elapsed_s': 0.0,
                            'detail': {'cache_hits': len(results)}}])
        _record_cache_health(spec.manifest_path, len(results), len(results), 0)
        logger.info(
            f'{spec.completion_label} cache hit: {len(frame)} results loaded from '
            f'{spec.cache_path}; saved to {path}')
        return frame

    import torch

    for name in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
        os.environ.setdefault(name, '1')
    mp.set_start_method('spawn', force=True)

    print_banner(spec.banner)
    logger.info(spec.start_message.format(count=len(requested_genomes)))
    if cached:
        logger.info(f'Reusing {len(cached)} content-addressed result(s); '
                    f'evaluating {len(misses)} cache miss(es)')
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
    for index, genome in enumerate(misses):
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
        if completed % 50 == 0 or completed == len(misses):
            elapsed = time.time() - started
            rate = completed / max(elapsed, 1e-12)
            valid = sum(1 for result in results if result.get('valid', False))
            logger.info(
                f"Progress: {completed}/{len(misses)} "
                f"({rate:.1f} {spec.progress_noun}/sec, {valid} valid, "
                f"{elapsed:.0f}s elapsed)")

    results = collect_results(
        status_queue, task_queue, stop_event, workers, spawn_worker, misses,
        spec.application, spec.manifest_path, progress=report_progress)
    for process in workers.values():
        process.join(timeout=30)
    _record_cache_health(
        spec.manifest_path, len(requested_genomes), len(cached), len(misses))

    if spec.cache_path is not None and spec.protocol_id is not None and implementation:
        from pipeline.screening.result_cache import store_completed_results
        store_completed_results(
            spec.cache_path, spec.application, spec.protocol_id,
            implementation, misses, results)

    by_id = {candidate_id(genome): result
             for genome, result in zip(misses, results)}
    ordered = []
    for index, genome in enumerate(requested_genomes):
        if index in cached:
            result = dict(cached[index])
            result['cache_hit'] = True
        else:
            result = dict(by_id[candidate_id(genome)])
            result['cache_hit'] = False
        result['genome'] = str(genome)
        result['material_class'] = genome[0]
        ordered.append(result)
    frame = pd.DataFrame(ordered)
    path = save_screening_db(
        frame, db_filename, subdir=spec.output_subdir)
    logger.info(
        f"{spec.completion_label} complete. {len(frame)} results saved to {path}")
    return frame
