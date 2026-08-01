#!/usr/bin/env python3
# Optional quantum transition-state validation.
"""
GPU-Accelerated VQE Transition-State Quantum Chemistry (CUDA-Q).

For the single best catalyst from each material class:
  1. Build molecular Hamiltonian for the C-H / O-O bond cleavage transition state
  2. Map to qubit Hamiltonian via Jordan-Wigner transform
  3. UCCSD ansatz VQE on CUDA-Q nvidia target (multi-GPU)
  4. Report correlation energy and comparison to classical reference

This script is designed to run in the quantum-env (CUDA-Q).
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, Optional

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.common.utils import (
    VQE_DIR, Ha_to_eV, setup_logger, print_banner, save_json,
)

logger = setup_logger('vqe_ts', 'vqe/vqe_transition_state.log')


# ═══════════════════════════════════════════════════════════════════════════════
# HAMILTONIAN CONSTRUCTION
# ═══════════════════════════════════════════════════════════════════════════════

def build_ch_splitting_hamiltonian(n_qubits: int = 4) -> list:
    """
    Build a model Hamiltonian for the C-H bond splitting transition state.
    
    Active space: σ(C-H) bonding and σ*(C-H) antibonding orbitals
    plus metal d-orbitals involved in activation.
    
    Uses a Jordan-Wigner-mapped second-quantized Hamiltonian with
    one- and two-body integrals derived from DFT orbital energies.
    
    Returns a list of (coefficient, Pauli_string) tuples.
    """
    # Model Hamiltonian terms (from DFT orbital energies in CH₄/catalyst system)
    # These are representative values for a C-H activation transition state
    hamiltonian_terms = [
        (-39.50, 'IIII'),     # Nuclear repulsion + core energy
        (-0.22, 'ZZII'),      # σ(C-H) bonding orbital
        (-0.18, 'IIZZ'),      # Metal d-orbital
        (0.12, 'ZIZI'),       # σ-d hybridization
        (-0.04, 'XXYY'),      # Exchange coupling
        (-0.04, 'YYXX'),      # Exchange coupling (Hermitian partner)
        (0.17, 'ZZIZ'),       # σ*(C-H) antibonding
        (-0.05, 'IZIZ'),      # d-σ* interaction
        (0.08, 'ZZZZ'),       # Two-body Coulomb
        (-0.02, 'XXII'),      # Single excitation
        (-0.02, 'IIXX'),      # Single excitation
        (0.03, 'ZYZY'),       # Spin-orbit coupling
    ]
    return hamiltonian_terms


def build_orr_hamiltonian(n_qubits: int = 4) -> list:
    """
    Build a model Hamiltonian for the O-O bond cleavage in OOH* intermediate.
    
    Active space: σ(O-O) and π*(O-O) orbitals with metal d participation.
    """
    hamiltonian_terms = [
        (-148.00, 'IIII'),    # Core energy
        (-0.30, 'ZZII'),      # σ(O-O) bonding
        (-0.25, 'IIZZ'),      # π*(O-O) antibonding
        (0.15, 'ZIZI'),       # Metal-O hybridization
        (-0.06, 'XXYY'),      # Exchange
        (-0.06, 'YYXX'),      # Exchange (Hermitian)
        (0.20, 'ZZIZ'),       # σ*(O-O)
        (-0.08, 'IZIZ'),      # Spin coupling
        (0.10, 'ZZZZ'),       # Two-electron Coulomb
        (-0.03, 'XXII'),      # Single excitation
        (-0.03, 'IIXX'),      # Single excitation
    ]
    return hamiltonian_terms


def exact_ground_energy(hamiltonian_terms: list, n_qubits: int) -> float:
    """Classically diagonalize a small Pauli Hamiltonian for VQE validation."""
    if n_qubits > 12:
        raise ValueError('exact VQE benchmark is limited to at most 12 qubits')
    pauli = {
        'I': np.eye(2),
        'X': np.array([[0, 1], [1, 0]], complex),
        'Y': np.array([[0, -1j], [1j, 0]], complex),
        'Z': np.diag([1, -1]),
    }
    hamiltonian = np.zeros((2 ** n_qubits, 2 ** n_qubits), complex)
    for coefficient, word in hamiltonian_terms:
        if len(word) != n_qubits or any(symbol not in pauli for symbol in word):
            raise ValueError('Pauli word does not match the declared qubit count')
        term = np.array([[1.0]], complex)
        for symbol in word:
            term = np.kron(term, pauli[symbol])
        hamiltonian += float(coefficient) * term
    if not np.allclose(hamiltonian, hamiltonian.conj().T):
        raise ValueError('VQE Hamiltonian is not Hermitian')
    return float(np.linalg.eigvalsh(hamiltonian)[0])


# ═══════════════════════════════════════════════════════════════════════════════
# VQE SOLVER
# ═══════════════════════════════════════════════════════════════════════════════

def run_vqe(hamiltonian_terms: list, n_qubits: int = 4,
            n_layers: int = 3, max_iter: int = 3000,
            initial_theta: Optional[list] = None,
            target: str = 'nvidia') -> Dict:
    """
    Run VQE using CUDA-Q with a hardware-efficient ansatz.
    
    Args:
        hamiltonian_terms: List of (coeff, pauli_string) tuples
        n_qubits: Number of qubits
        n_layers: Number of ansatz layers
        max_iter: Maximum COBYLA iterations
        initial_theta: Initial variational parameters
        target: CUDA-Q target ('nvidia' for GPU, 'default' for CPU)
        
    Returns: Dict with optimized energy, parameters, etc.
    """
    try:
        import cudaq
        from cudaq import spin
        HAS_CUDAQ = True
    except ImportError:
        HAS_CUDAQ = False

    if not HAS_CUDAQ:
        return _mock_vqe_result(hamiltonian_terms)

    # Set target
    # CUDA-Q 0.12 names its local state-vector CPU target qpp-cpu.
    resolved_target = 'qpp-cpu' if target == 'default' else target
    cudaq.set_target(resolved_target)

    # Build spin operator
    H = 0.0 * spin.i(0)  # initialize
    for coeff, pauli_str in hamiltonian_terms:
        term = coeff
        for i, p in enumerate(pauli_str):
            if p == 'I':
                term = term * spin.i(i)
            elif p == 'X':
                term = term * spin.x(i)
            elif p == 'Y':
                term = term * spin.y(i)
            elif p == 'Z':
                term = term * spin.z(i)
        H += term

    # Number of variational parameters
    n_params = n_qubits * n_layers * 2  # Ry + Rz per qubit per layer + entangling

    # Define ansatz kernel
    @cudaq.kernel
    def ansatz(thetas: list[float]):
        q = cudaq.qvector(n_qubits)

        # Initial state: Hartree-Fock reference (half-filled)
        for i in range(n_qubits // 2):
            x(q[i])

        # Parameterized layers
        param_idx = 0
        for layer in range(n_layers):
            for i in range(n_qubits):
                ry(thetas[param_idx], q[i])
                param_idx += 1
            for i in range(n_qubits):
                rz(thetas[param_idx], q[i])
                param_idx += 1
            # Entangling gates
            for i in range(n_qubits - 1):
                cx(q[i], q[i + 1])
            if n_qubits > 1:
                cx(q[n_qubits - 1], q[0])

    # Initial parameters
    if initial_theta is None:
        initial_theta = [0.01] * n_params

    # Run VQE optimization using the CUDA-Q 0.12 optimizer contract.
    logger.info(f"  Running CUDA-Q VQE: {n_qubits} qubits, {n_params} parameters, {n_layers} layers")
    optimizer = cudaq.optimizers.COBYLA()
    optimizer.max_iterations = int(max_iter)
    optimizer.initial_parameters = list(initial_theta)
    result = cudaq.vqe(ansatz, H, optimizer, n_params)

    optimal_energy = result.energy if hasattr(result, 'energy') else result[0]
    optimal_params = (result.optimal_parameters
                      if hasattr(result, 'optimal_parameters') else result[1])
    exact_energy = exact_ground_energy(hamiltonian_terms, n_qubits)
    variational_gap = float(optimal_energy) - exact_energy
    variational_valid = variational_gap >= -1e-8
    if not variational_valid:
        raise RuntimeError(
            f'VQE energy violates the variational bound by {-variational_gap:.3e} Ha')
    chemical_accuracy_Ha = 1.6e-3
    benchmark_passed = (
        variational_valid and variational_gap <= chemical_accuracy_Ha)

    logger.info(f"  VQE converged: E = {optimal_energy:.6f} Ha ({optimal_energy * Ha_to_eV:.4f} eV)")

    return {
        'energy_Ha': float(optimal_energy),
        'energy_eV': float(optimal_energy * Ha_to_eV),
        'optimal_params': [float(p) for p in optimal_params] if optimal_params else [],
        'n_qubits': n_qubits,
        'n_layers': n_layers,
        'n_params': n_params,
        'max_iter': max_iter,
        'target': resolved_target,
        'evidence_level': 'toy_hamiltonian',
        'catalyst_specific_hamiltonian': False,
        'benchmarked': benchmark_passed,
        'exact_ground_energy_Ha': exact_energy,
        'variational_gap_Ha': variational_gap,
        'variational_bound_valid': variational_valid,
        'benchmark_tolerance_Ha': chemical_accuracy_Ha,
    }


def _mock_vqe_result(hamiltonian_terms: list) -> Dict:
    """Generate mock VQE results when CUDA-Q is not available."""
    logger.warning("CUDA-Q not available. Generating mock VQE results.")
    # Extract the constant (identity) term as the base energy
    base_energy = sum(c for c, p in hamiltonian_terms if p == 'IIII')
    # Add approximate correlation correction
    correlation = -0.2  # typical correlation energy
    energy = base_energy + correlation

    return {
        'energy_Ha': float(energy),
        'energy_eV': float(energy * Ha_to_eV),
        'optimal_params': [0.01, 0.005],
        'n_qubits': 4,
        'n_layers': 2,
        'mock': True,
        'evidence_level': 'mock',
        'catalyst_specific_hamiltonian': False,
        'benchmarked': False,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED WORKFLOW
# ═══════════════════════════════════════════════════════════════════════════════

def validate_transition_state(catalyst_name: str, reaction_type: str = 'CH_split',
                               target: str = 'nvidia') -> Dict:
    """
    Full VQE transition-state validation for a champion catalyst.
    
    Args:
        catalyst_name: Identifier for the catalyst
        reaction_type: 'CH_split' (methane pyrolysis) or 'ORR' (fuel cell)
        target: CUDA-Q target
        
    Returns: Dict with VQE results
    """
    print_banner(f"CUDA-Q VQE: {catalyst_name} ({reaction_type})")

    if reaction_type == 'CH_split':
        H_terms = build_ch_splitting_hamiltonian()
    elif reaction_type == 'ORR':
        H_terms = build_orr_hamiltonian()
    else:
        raise ValueError(f"Unknown reaction type: {reaction_type}")

    result = run_vqe(
        H_terms, n_qubits=4, n_layers=3, max_iter=3000, target=target)
    result['catalyst_name'] = catalyst_name
    result['reaction_type'] = reaction_type

    save_json(result, f"vqe_{catalyst_name}_{reaction_type}.json", subdir="vqe")
    return result


if __name__ == '__main__':
    # Test with mock
    r1 = validate_transition_state("NiBi_champion", "CH_split", target="default")
    r2 = validate_transition_state("FeN4_champion", "ORR", target="default")
    print(json.dumps(r1, indent=2))
    print(json.dumps(r2, indent=2))
