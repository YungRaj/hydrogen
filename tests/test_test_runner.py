#!/usr/bin/env python3
"""Unit contracts for the multi-environment repository test runner."""

import json
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import run_tests


def test_default_selection_is_portable_only():
    args = run_tests.parse_args([])
    selected = run_tests.select_suites(args)
    assert selected
    assert {suite.category for suite in selected} == {"portable"}
    assert "gpu-affinity" not in {suite.name for suite in selected}
    assert "vqe" not in {suite.name for suite in selected}


def test_opt_ins_and_explicit_selection_are_deterministic():
    args = run_tests.parse_args(["--include-resolution", "--include-vqe"])
    categories = {suite.category for suite in run_tests.select_suites(args)}
    assert categories == {"portable", "resolution", "vqe-smoke"}

    args = run_tests.parse_args(["--profile", "merge"])
    assert [suite.name for suite in run_tests.select_suites(args)] == [
        "scientific", "modular-multiphysics", "coupling", "reactor-merge"]

    args = run_tests.parse_args(["--suite", "coupling", "--suite", "pipeline"])
    assert [suite.name for suite in run_tests.select_suites(args)] == [
        "pipeline", "coupling"]


def test_environment_resolution_is_portable_and_fail_closed():
    suite = run_tests.Suite("sample", "tests/sample.py", "sample-env", 10)
    command, error = run_tests.environment_command(suite, "/portable/conda")
    assert error is None
    assert command == [
        "/portable/conda", "run", "--no-capture-output", "-n", "sample-env",
        "python", "tests/sample.py"]

    with patch.dict("os.environ", {"CONDA_DEFAULT_ENV": "sample-env"}):
        command, error = run_tests.environment_command(suite, None)
    assert error is None
    assert command == [sys.executable, "tests/sample.py"]

    with patch.dict("os.environ", {"CONDA_DEFAULT_ENV": "other-env"}):
        command, error = run_tests.environment_command(suite, None)
    assert command == []
    assert "conda was not found" in error


def test_json_report_is_atomic_and_complete():
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / "nested" / "summary.json"
        report = {"schema_version": 1, "success": True, "results": []}
        run_tests.write_report(target, report)
        assert json.loads(target.read_text()) == report
        assert not target.with_suffix(".json.tmp").exists()


TESTS = (
    test_default_selection_is_portable_only,
    test_opt_ins_and_explicit_selection_are_deterministic,
    test_environment_resolution_is_portable_and_fail_closed,
    test_json_report_is_atomic_and_complete,
)


def main():
    for test in TESTS:
        test()
        print("PASS", test.__name__)
    print(f"{len(TESTS)}/{len(TESTS)} test-runner contracts passed")


if __name__ == "__main__":
    main()
