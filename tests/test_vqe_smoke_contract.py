#!/usr/bin/env python3
"""Bounded CUDA-Q smoke contract without a production convergence claim."""

import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.validation.vqe_transition_state import run_vqe


def test_cudaq_vqe_execution_is_bounded_and_variational():
    # A one-qubit Hamiltonian validates target selection, kernel construction,
    # optimization, and result parsing without paying for chemical accuracy.
    terms = [(-1.0, "Z")]
    result = run_vqe(
        terms, n_qubits=1, n_layers=1, max_iter=24, target="qpp-cpu")
    assert result.get("mock") is not True
    assert math.isfinite(result["energy_Ha"])
    assert result["variational_bound_valid"] is True
    assert result["energy_Ha"] >= -1.0 - 1e-8
    assert result["n_qubits"] == 1
    assert result["n_layers"] == 1
    assert result["max_iter"] == 24
    assert result["target"] == "qpp-cpu"
    # Smoke success is backend evidence only, never candidate evidence.
    assert result["evidence_level"] == "toy_hamiltonian"
    assert result["catalyst_specific_hamiltonian"] is False


if __name__ == "__main__":
    test_cudaq_vqe_execution_is_bounded_and_variational()
    print("PASS test_cudaq_vqe_execution_is_bounded_and_variational")
