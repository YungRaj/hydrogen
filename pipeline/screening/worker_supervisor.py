"""Health protocol and fail-safe collection for long-running GPU workers."""

from __future__ import annotations

import json
import queue
import threading
import time
from datetime import datetime, timezone
from pathlib import Path


def emit(status_queue, kind: str, worker_id: int, payload=None):
    status_queue.put((kind, worker_id, time.monotonic(), payload))


def start_heartbeat(status_queue, worker_id: int, stop_event,
                    interval_s: float = 10.0) -> threading.Thread:
    def beat():
        while not stop_event.wait(interval_s):
            emit(status_queue, 'heartbeat', worker_id)
    thread = threading.Thread(target=beat, name=f'heartbeat-{worker_id}', daemon=True)
    thread.start()
    return thread


def write_health_manifest(path, application: str, status: str, workers: dict,
                          completed: int, expected: int, events: list,
                          error: str | None = None):
    payload = {
        'schema_version': 1,
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'application': application,
        'status': status,
        'completed': completed,
        'expected': expected,
        'error': error,
        'workers': {
            str(worker_id): {
                'pid': process.pid,
                'exitcode': process.exitcode,
                'alive': process.is_alive(),
            } for worker_id, process in workers.items()
        },
        'events': events[-500:],
    }
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(payload, indent=2) + '\n')
    temporary.replace(target)


def collect_results(status_queue, task_queue, stop_event, workers: dict,
                    spawn_replacement, genomes, application: str,
                    manifest_path, startup_timeout_s: float = 300.0,
                    heartbeat_timeout_s: float = 120.0,
                    result_timeout_s: float = 1800.0,
                    max_restarts_per_worker: int = 1,
                    progress=None):
    """Collect results with startup/heartbeat checks and bounded task recovery."""
    expected = len(genomes)
    pending = set(range(expected))
    results = {}
    active = {worker_id: set() for worker_id in workers}
    launched = {worker_id: time.monotonic() for worker_id in workers}
    last_seen = dict(launched)
    ready = set()
    restarts = {worker_id: 0 for worker_id in workers}
    events = []
    last_result = time.monotonic()

    def record(kind, worker_id, detail=None):
        events.append({'kind': kind, 'worker_id': worker_id,
                       'elapsed_s': round(time.monotonic() - launched.get(worker_id, 0), 3),
                       'detail': detail})

    def recover(worker_id, reason):
        nonlocal workers
        process = workers[worker_id]
        if process.is_alive():
            process.terminate()
            process.join(timeout=10)
        leased = sorted(active.get(worker_id, set()) & pending)
        for index in leased:
            task_queue.put((index, genomes[index]))
        active[worker_id] = set()
        record('worker_failure', worker_id, {'reason': reason, 'requeued': leased})
        if restarts[worker_id] >= max_restarts_per_worker:
            raise RuntimeError(
                f'worker {worker_id} failed after {restarts[worker_id]} restart(s): {reason}')
        restarts[worker_id] += 1
        replacement = spawn_replacement(worker_id)
        workers[worker_id] = replacement
        launched[worker_id] = last_seen[worker_id] = time.monotonic()
        ready.discard(worker_id)
        record('worker_restarted', worker_id, {'restart': restarts[worker_id]})

    try:
        while pending:
            try:
                kind, worker_id, timestamp, payload = status_queue.get(timeout=5.0)
                last_seen[worker_id] = time.monotonic()
                if kind == 'ready':
                    ready.add(worker_id)
                    record('ready', worker_id, payload)
                elif kind == 'heartbeat':
                    pass
                elif kind == 'started':
                    if payload in pending:
                        active.setdefault(worker_id, set()).add(payload)
                elif kind == 'result':
                    index, result = payload
                    active.setdefault(worker_id, set()).discard(index)
                    if index in pending:
                        results[index] = result
                        pending.remove(index)
                        last_result = time.monotonic()
                        if progress:
                            progress(len(results), list(results.values()))
                elif kind == 'fatal':
                    recover(worker_id, str(payload))
            except queue.Empty:
                pass

            now = time.monotonic()
            for worker_id in list(workers):
                process = workers[worker_id]
                if not process.is_alive():
                    recover(worker_id, f'exited with code {process.exitcode}')
                elif worker_id not in ready and now - launched[worker_id] > startup_timeout_s:
                    recover(worker_id, 'startup acknowledgement timeout')
                elif worker_id in ready and now - last_seen[worker_id] > heartbeat_timeout_s:
                    recover(worker_id, 'heartbeat timeout')
            if now - last_result > result_timeout_s:
                raise RuntimeError(
                    f'no completed candidate for {result_timeout_s:.0f}s; {len(pending)} pending')

        stop_event.set()
        write_health_manifest(manifest_path, application, 'complete', workers,
                              len(results), expected, events)
        return [results[index] for index in range(expected)]
    except Exception as exc:
        stop_event.set()
        for process in workers.values():
            if process.is_alive():
                process.terminate()
        write_health_manifest(manifest_path, application, 'failed', workers,
                              len(results), expected, events, str(exc))
        raise
