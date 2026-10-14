#!/usr/bin/env python3
"""Run the pytest-style contract files as one ``run_tests.py`` suite."""

from pathlib import Path
import sys

TESTS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(TESTS_DIR.parent))

PYTEST_CONTRACTS = (
    "test_atomistic_validation_contracts.py",
    "test_attribution_contracts.py",
    "test_electrochemical_adapter_status.py",
    "test_native_electrochemistry_contracts.py",
    "test_native_model_contracts.py",
    "test_objectives.py",
    "test_orr_geometry_contracts.py",
    "test_physics_fixes.py",
    "test_policy_pilot.py",
    "test_prospective_contracts.py",
    "test_qe_gpu_setup_contract.py",
    "test_review_fixes.py",
    "test_screening_gates.py",
    "test_screening_multiphysics.py",
    "test_sensitivity_contracts.py",
    "test_supported_nanoparticle_genome.py",
    "test_volcano_crosscheck.py",
)


if __name__ == "__main__":
    import pytest

    sys.exit(pytest.main([
        "-q", "-p", "no:cacheprovider", "--rootdir", str(TESTS_DIR.parent),
        *(str(TESTS_DIR / name) for name in PYTEST_CONTRACTS),
    ]))
