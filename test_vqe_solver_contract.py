#!/usr/bin/env python3
"""CUDA-Q environment contract: compare VQE with exact diagonalization."""

import numpy as np

from pipeline.validation.vqe_transition_state import (
    build_ch_splitting_hamiltonian, run_vqe)


PAULI = {
    'I': np.eye(2),
    'X': np.array([[0, 1], [1, 0]], complex),
    'Y': np.array([[0, -1j], [1j, 0]], complex),
    'Z': np.diag([1, -1]),
}


def matrix(pauli_string):
    result = np.array([[1.0]], complex)
    for symbol in pauli_string:
        result = np.kron(result, PAULI[symbol])
    return result


def main():
    terms = build_ch_splitting_hamiltonian()
    hamiltonian = sum(float(coefficient) * matrix(pauli)
                      for coefficient, pauli in terms)
    exact = float(np.linalg.eigvalsh(hamiltonian)[0])
    result = run_vqe(
        terms, target='qpp-cpu', n_layers=3, max_iter=3000)
    assert not result.get('mock', False)
    assert result['energy_Ha'] >= exact - 1e-8
    # Four qubits are cheap enough to demand a useful ansatz/optimizer result.
    assert result['energy_Ha'] - exact < 1e-3
    assert result['benchmarked'] is True
    print({'exact_Ha': exact, 'vqe_Ha': result['energy_Ha'],
           'gap_Ha': result['energy_Ha'] - exact})


if __name__ == '__main__':
    main()
