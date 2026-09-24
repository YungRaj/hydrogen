"""Static shapes for values exchanged by the six orchestration stages.

These contracts intentionally retain dictionaries and pandas objects at
runtime, so existing artifacts and scientific implementations do not change.
They remove the *type erasure* at stage boundaries and provide an incremental
path toward richer domain objects.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, TypedDict

from pipeline.data_models.quantum import DFTResult, VQEResult
from pipeline.data_models.reactors import ReactorResult


class TableLike(Protocol):
    """Minimal tabular surface passed between screening stages."""

    def __len__(self) -> int:
        """Return the number of candidate records."""
        ...

    def head(self, n: int) -> "TableLike":
        """Return at most the first ``n`` candidate records.

        Args:
            n: Maximum number of records to retain.

        Returns:
            A table-like view containing the leading records.
        """
        ...

    def iterrows(self):
        """Yield index and candidate-record pairs."""
        ...


class _DiscoveryStateOptional(TypedDict, total=False):
    best_E_act: float
    best_coking: float


class DiscoveryState(_DiscoveryStateOptional):
    """Persistable summary emitted by branch-and-bound discovery."""
    pareto_size: int
    total_evaluated: int
    valid_count: int
    top_catalysts_count: int
    dft_resolution_count: int
    admissibility: dict[str, Any]
    candidate_dispositions: dict[str, int]


class DiscoveryProducts(TypedDict):
    """In-process values routed from discovery to later stages."""
    design_space_sizes: dict[str, int]
    pareto_genomes: Any
    screening_database: Any
    top_catalysts: Any
    dft_candidates: Any


class _ReactorBatchStateOptional(TypedDict, total=False):
    equilibrium_check: dict[str, Any]
    solids_scorecard: dict[str, Any]
    best_conversion: float | None
    best_conversion_scope: str
    mmbcr_max_conversion: float | None


class ReactorBatchState(_ReactorBatchStateOptional):
    """Persistable summary emitted by a reactor batch."""
    n_simulations: int


class ReactorBatchProducts(TypedDict):
    """Detailed reactor results retained for downstream consumers."""
    reactor_results: list[ReactorResult]


class DFTState(TypedDict):
    """Persistable DFT execution counts."""
    n_validated: int
    n_converged: int
    n_failed: int


class DFTProducts(TypedDict):
    """Candidate-specific DFT results and isolated failures."""
    dft_results: list[DFTResult]
    failures: list[dict[str, str]]


class VQEState(TypedDict):
    """Persistable VQE execution counts."""
    n_vqe_runs: int


class VQEProducts(TypedDict):
    """Candidate-specific VQE results."""
    vqe_results: list[VQEResult]


class _FuelCellStateOptional(TypedDict, total=False):
    best_power_W_cm2: float
    best_efficiency: float
    min_overpotential_V: float


class FuelCellState(_FuelCellStateOptional):
    """Persistable cathode and PEMFC summary."""
    n_cathodes_screened: int
    n_valid: int
    n_pemfc_simulations: int


class FuelCellProducts(TypedDict):
    """Detailed fuel-cell screening, cell, and stack products."""
    cathode_database: Any
    valid_cathodes: Any
    pemfc_results: list[dict[str, Any]]
    stack_result: dict[str, Any]


class ReportState(TypedDict):
    """Persistable report-stage summary."""
    report_path: str


class ReportProducts(TypedDict):
    """Generated report artifact path."""
    report_path: Path


class SelectedCandidates(TypedDict):
    """Typed restart handoff reconstructed from a screening table."""
    screening_database: Any
    admissibility: dict[str, Any]
    top_catalysts: Any
    dft_candidates: Any
