"""Persistent candidate-keyed queue for high-fidelity validation tasks."""

from __future__ import annotations

import json
import sqlite3
import time
from pathlib import Path

from pipeline.search.discovery import candidate_id


class ValidationTaskQueue:
    """Manage durable, candidate-keyed validation work leases.
    """
    def __init__(self, database: str | Path):
        self.database = str(database)
        Path(self.database).parent.mkdir(parents=True, exist_ok=True)
        with self._connect() as conn:
            conn.execute("""CREATE TABLE IF NOT EXISTS validation_tasks (
                application TEXT NOT NULL, candidate_id TEXT NOT NULL,
                task_type TEXT NOT NULL, protocol_id TEXT NOT NULL,
                genome TEXT NOT NULL, status TEXT NOT NULL,
                attempts INTEGER NOT NULL DEFAULT 0,
                result_path TEXT, error TEXT,
                created_at REAL NOT NULL, updated_at REAL NOT NULL,
                PRIMARY KEY(application, candidate_id, task_type, protocol_id)
            )""")
            conn.execute("CREATE INDEX IF NOT EXISTS validation_task_status_idx "
                         "ON validation_tasks(application, task_type, status)")

    def _connect(self):
        return sqlite3.connect(self.database, timeout=60)

    def enqueue(self, application: str, genome: tuple, task_type: str,
                protocol_id: str) -> str:
        """Insert a candidate-specific validation task if it is not already queued.

        Args:
            application: Scientific application or objective family.
            genome: Encoded catalyst composition and structural configuration.
            task_type: Validation calculation category.
            protocol_id: Scientific protocol identifier.

        Returns:
            The durable identifier of the queued validation task.
        """
        cid = candidate_id(genome)
        now = time.time()
        with self._connect() as conn:
            conn.execute(
                "INSERT OR IGNORE INTO validation_tasks VALUES "
                "(?,?,?,?,?,'queued',0,NULL,NULL,?,?)",
                (application, cid, task_type, protocol_id,
                 json.dumps(genome, separators=(',', ':')), now, now))
        return cid

    def recover_stale(self, stale_after_s: float = 86400) -> int:
        """Return expired task leases to the pending queue.

        Args:
            stale_after_s: Lease age after which an unfinished task may be recovered.

        Returns:
            The computed recover stale numeric value.
        """
        cutoff = time.time() - stale_after_s
        with self._connect() as conn:
            cursor = conn.execute(
                "UPDATE validation_tasks SET status='queued', "
                "error='recovered_stale_running_task', updated_at=? "
                "WHERE status='running' AND updated_at<?",
                (time.time(), cutoff))
            return int(cursor.rowcount)

    def claim(self, application: str, candidate: str, task_type: str,
              protocol_id: str) -> bool:
        """Atomically lease the next pending validation task.

        Args:
            application: Scientific application or objective family.
            candidate: Candidate record being evaluated.
            task_type: Validation calculation category.
            protocol_id: Scientific protocol identifier.

        Returns:
            The leased task record, or `None` when no task is available.
        """
        with self._connect() as conn:
            cursor = conn.execute(
                "UPDATE validation_tasks SET status='running', attempts=attempts+1, "
                "error=NULL, updated_at=? WHERE application=? AND candidate_id=? "
                "AND task_type=? AND protocol_id=? AND status IN ('queued','failed')",
                (time.time(), application, candidate, task_type, protocol_id))
            return cursor.rowcount == 1

    def finish(self, application: str, candidate: str, task_type: str,
               protocol_id: str, converged: bool,
               result_path: str | None = None, error: str | None = None) -> None:
        """Record the terminal status and result for a leased task.

        Args:
            application: Scientific application or objective family.
            candidate: Candidate record being evaluated.
            task_type: Validation calculation category.
            protocol_id: Scientific protocol identifier.
            converged: Whether the external calculation met its convergence contract.
            result_path: Path to the validated calculation result.
            error: Failure detail, if execution did not succeed.
        """
        with self._connect() as conn:
            conn.execute(
                "UPDATE validation_tasks SET status=?, result_path=?, error=?, "
                "updated_at=? WHERE application=? AND candidate_id=? "
                "AND task_type=? AND protocol_id=?",
                ('converged' if converged else 'failed', result_path,
                 None if converged else str(error or 'not_converged')[:500],
                 time.time(), application, candidate, task_type, protocol_id))

    def summary(self, application: str, task_type: str) -> dict:
        """Count validation tasks by status.

        Args:
            application: Scientific application or objective family.
            task_type: Validation calculation category.

        Returns:
            Counts of validation tasks grouped by status.
        """
        with self._connect() as conn:
            rows = conn.execute(
                "SELECT status, COUNT(*) FROM validation_tasks "
                "WHERE application=? AND task_type=? GROUP BY status",
                (application, task_type)).fetchall()
        return {status: int(count) for status, count in rows}
