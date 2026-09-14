#!/usr/bin/env python3
"""Run repository contracts in their owning Conda environments.

The default profile contains only portable suites. Quantum and real-GPU
contracts are explicit opt-ins because they require specialized software and
may consume substantial compute.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import platform
import shutil
import signal
import subprocess
import sys
import time
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REPORT = REPO_ROOT / "results" / "test_summary.json"


@dataclass(frozen=True)
class Suite:
    name: str
    path: str
    environment: str
    timeout_s: int
    category: str = "portable"


@dataclass
class SuiteResult:
    name: str
    path: str
    environment: str
    category: str
    status: str
    returncode: int | None
    duration_s: float
    command: list[str]
    stdout: str = ""
    stderr: str = ""
    reason: str | None = None


FAIRCHEM_ENV = os.environ.get("HYDROGEN_TEST_FAIRCHEM_ENV", "fairchem-env")
QUANTUM_ENV = os.environ.get("HYDROGEN_TEST_QUANTUM_ENV", "quantum-env")

SUITES = (
    Suite("pipeline", "tests/test_pipeline.py", FAIRCHEM_ENV, 180),
    Suite("scientific", "tests/test_scientific_contracts.py", FAIRCHEM_ENV, 300),
    Suite("modular-multiphysics", "tests/test_modular_multiphysics.py", FAIRCHEM_ENV, 180),
    Suite("component-replacement", "tests/test_component_replacement_contracts.py", FAIRCHEM_ENV, 120),
    Suite("architecture", "tests/test_architecture_unit_contracts.py", FAIRCHEM_ENV, 120),
    Suite("coupling", "tests/test_coupling_contracts.py", FAIRCHEM_ENV, 120),
    Suite("experimental-data", "tests/test_experimental_data_contract.py", FAIRCHEM_ENV, 120),
    Suite("test-runner", "tests/test_test_runner.py", FAIRCHEM_ENV, 60),
    Suite("repository-audit", "audit_pipeline.py", FAIRCHEM_ENV, 180),
    Suite("resolution", "tests/test_resolution_contracts.py", QUANTUM_ENV, 180, "resolution"),
    Suite("vqe", "tests/test_vqe_solver_contract.py", QUANTUM_ENV, 900, "vqe"),
    Suite("gpu-affinity", "tests/test_gpu_affinity_contract.py", FAIRCHEM_ENV, 900, "gpu"),
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--suite", action="append", choices=[suite.name for suite in SUITES],
        help="run only the named suite; may be supplied more than once")
    parser.add_argument(
        "--include-resolution", action="store_true",
        help="include candidate-Hamiltonian tests in the quantum environment")
    parser.add_argument(
        "--include-vqe", action="store_true",
        help="include the potentially long real CUDA-Q VQE contract")
    parser.add_argument(
        "--include-gpu", action="store_true",
        help="include the real multi-GPU affinity workload")
    parser.add_argument(
        "--all", action="store_true",
        help="include portable, resolution, VQE, and GPU suites")
    parser.add_argument(
        "--timeout", type=int,
        help="override every selected suite timeout in seconds")
    parser.add_argument(
        "--output", type=Path, default=DEFAULT_REPORT,
        help=f"JSON result path (default: {DEFAULT_REPORT.relative_to(REPO_ROOT)})")
    parser.add_argument(
        "--verbose", action="store_true",
        help="print captured output for successful suites as well as failures")
    parser.add_argument(
        "--dry-run", action="store_true",
        help="print the execution plan without launching a test")
    return parser.parse_args(argv)


def select_suites(args: argparse.Namespace) -> list[Suite]:
    if args.suite:
        requested = set(args.suite)
        return [suite for suite in SUITES if suite.name in requested]
    categories = {"portable"}
    if args.all or args.include_resolution:
        categories.add("resolution")
    if args.all or args.include_vqe:
        categories.add("vqe")
    if args.all or args.include_gpu:
        categories.add("gpu")
    return [suite for suite in SUITES if suite.category in categories]


def environment_command(suite: Suite, conda: str | None) -> tuple[list[str], str | None]:
    active = os.environ.get("CONDA_DEFAULT_ENV")
    if conda:
        return [conda, "run", "--no-capture-output", "-n", suite.environment,
                "python", suite.path], None
    if active == suite.environment:
        return [sys.executable, suite.path], None
    return [], (
        f"conda was not found and active environment {active!r} is not "
        f"the required {suite.environment!r}")


def terminate_process(process: subprocess.Popen[str]) -> None:
    if process.poll() is not None:
        return
    if os.name == "posix":
        os.killpg(process.pid, signal.SIGTERM)
    else:
        process.terminate()
    try:
        process.wait(timeout=5)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            os.killpg(process.pid, signal.SIGKILL)
        else:
            process.kill()
        process.wait()


def run_suite(suite: Suite, conda: str | None, timeout_s: int,
              verbose: bool) -> SuiteResult:
    command, preflight_error = environment_command(suite, conda)
    if preflight_error:
        return SuiteResult(
            suite.name, suite.path, suite.environment, suite.category,
            "failed", None, 0.0, command, reason=preflight_error)

    print(f"RUN  {suite.name:<24} env={suite.environment} timeout={timeout_s}s",
          flush=True)
    started = time.monotonic()
    process = subprocess.Popen(
        command, cwd=REPO_ROOT, text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, start_new_session=(os.name == "posix"))
    try:
        stdout, stderr = process.communicate(timeout=timeout_s)
        duration = time.monotonic() - started
        status = "passed" if process.returncode == 0 else "failed"
        reason = None if status == "passed" else f"exit code {process.returncode}"
    except subprocess.TimeoutExpired:
        terminate_process(process)
        stdout, stderr = process.communicate()
        duration = time.monotonic() - started
        status, reason = "timed_out", f"exceeded {timeout_s}s"

    marker = "PASS" if status == "passed" else "FAIL"
    print(f"{marker} {suite.name:<24} {duration:.1f}s", flush=True)
    if verbose or status != "passed":
        if stdout.strip():
            print(stdout.rstrip())
        if stderr.strip():
            print(stderr.rstrip(), file=sys.stderr)
    return SuiteResult(
        suite.name, suite.path, suite.environment, suite.category, status,
        process.returncode, round(duration, 3), command, stdout, stderr, reason)


def git_revision() -> str | None:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=REPO_ROOT, text=True,
            stderr=subprocess.DEVNULL, timeout=5).strip()
    except (OSError, subprocess.SubprocessError):
        return None


def write_report(path: Path, report: dict) -> None:
    target = path if path.is_absolute() else REPO_ROOT / path
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n")
    temporary.replace(target)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    suites = select_suites(args)
    conda = shutil.which("conda")

    if args.dry_run:
        for suite in suites:
            command, error = environment_command(suite, conda)
            timeout_s = args.timeout or suite.timeout_s
            print(json.dumps({
                "name": suite.name, "category": suite.category,
                "timeout_s": timeout_s, "command": command, "error": error}))
        return 0

    started_at = datetime.now(timezone.utc)
    results = [
        run_suite(suite, conda, args.timeout or suite.timeout_s, args.verbose)
        for suite in suites
    ]
    finished_at = datetime.now(timezone.utc)
    counts = {
        status: sum(result.status == status for result in results)
        for status in ("passed", "failed", "timed_out")
    }
    report = {
        "schema_version": 1,
        "repository": str(REPO_ROOT),
        "git_revision": git_revision(),
        "host": platform.node(),
        "platform": platform.platform(),
        "started_at_utc": started_at.isoformat(),
        "finished_at_utc": finished_at.isoformat(),
        "duration_s": round((finished_at - started_at).total_seconds(), 3),
        "counts": counts,
        "success": counts["failed"] == 0 and counts["timed_out"] == 0,
        "results": [asdict(result) for result in results],
    }
    write_report(args.output, report)
    print(f"\nSummary: {counts['passed']} passed, {counts['failed']} failed, "
          f"{counts['timed_out']} timed out")
    print(f"JSON: {args.output if args.output.is_absolute() else REPO_ROOT / args.output}")
    return 0 if report["success"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
