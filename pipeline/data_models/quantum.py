"""Named data shapes produced by classical and quantum validation workflows."""

from __future__ import annotations

from typing import Any, Literal, TypedDict


class SolverExecution(TypedDict):
    """Resource allocation recorded for a native solver invocation."""

    mpi_ranks: int
    omp_threads: int
    kpoint_pools: int
    image_groups: int


class SolverRunResult(TypedDict):
    """Process-level outcome of a Quantum ESPRESSO invocation."""

    converged: bool
    returncode: int
    output: str
    execution: SolverExecution
    timed_out: bool


class _DFTResultOptional(TypedDict, total=False):
    dft_energy_Ry: float | None
    dft_energy_eV: float | None
    max_force_Ry_bohr: float | None
    converged: bool
    returncode: int
    execution: SolverExecution
    resumed_existing_output: bool
    error: str


class DFTResult(_DFTResultOptional):
    """Candidate-specific DFT screening result with explicit evidence identity."""

    catalyst_name: str
    candidate_id: str
    material_class: str
    genome: str
    evidence_level: str


class NEBResult(TypedDict):
    """Parsed candidate-specific NEB convergence and barrier evidence."""

    converged: bool
    forward_barrier_eV: float | None
    reverse_barrier_eV: float | None
    path_metric: float | None
    candidate_specific: bool


class TransitionStateFrequencyResult(TypedDict):
    """Partial-Hessian frequencies used to validate a transition state."""

    frequencies_cm1: list[float]
    imaginary_count: int
    valid_transition_state: bool
    mode_vectors: list[list[float]]


class _VQEResultOptional(TypedDict, total=False):
    optimal_params: list[float]
    n_params: int
    max_iter: int
    target: str
    mock: bool
    exact_ground_energy_Ha: float
    variational_gap_Ha: float
    variational_bound_valid: bool
    benchmark_tolerance_Ha: float
    catalyst_name: str
    reaction_type: str
    hamiltonian: dict[str, Any]


class VQEResult(_VQEResultOptional):
    """VQE result that distinguishes toy, mock, and candidate-specific evidence."""

    energy_Ha: float
    energy_eV: float
    n_qubits: int
    n_layers: int
    evidence_level: str
    catalyst_specific_hamiltonian: bool
    benchmarked: bool


class ORREnsembleResult(TypedDict, total=False):
    """Multi-site ORR result with completeness and correction provenance."""

    complete: bool
    evidence_level: Literal['corrected_DFT_ensemble', 'incomplete']
    cases: list[dict[str, Any]]
    case_count: int
    expected_cases: int
    invalid_or_incomplete_rows: int
    orr_overpotential_V: float | None
    limiting_step: str
    best_case: dict[str, Any]
    uncertainty_V: float
