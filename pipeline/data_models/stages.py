"""Static shapes for values exchanged by the six orchestration stages.

These contracts intentionally retain dictionaries and pandas objects at
runtime, so existing artifacts and scientific implementations do not change.
They remove the *type erasure* at stage boundaries and provide an incremental
path toward richer domain objects.
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from typing import Any, Protocol, Sequence, TypedDict

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


CandidateBatch = TableLike | Sequence[object]
"""Candidate records accepted by validation stages.

Screening normally supplies a pandas table, while focused validation and
tests may supply an ordered sequence of genome records.  Naming this union at
the boundary avoids erasing the handoff to ``Any``.
"""


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
    pareto_genomes: Sequence[object]
    screening_database: TableLike
    top_catalysts: TableLike
    dft_candidates: CandidateBatch


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


class _SelectedCandidatesOptional(TypedDict, total=False):
    admissibility: dict[str, Any]


class SelectedCandidates(_SelectedCandidatesOptional):
    """Typed restart handoff reconstructed from a screening table."""
    screening_database: TableLike
    top_catalysts: TableLike
    dft_candidates: CandidateBatch


@dataclass(frozen=True, slots=True)
class CandidateSelection:
    """Explicit in-memory handoff from discovery to reactor and DFT stages.

    The same structure is constructed from a live discovery result and from a
    persisted screening table.  Downstream phases therefore do not depend on
    which execution route populated their inputs.
    """

    screening_database: TableLike
    top_catalysts: TableLike
    dft_candidates: CandidateBatch
    admissibility: dict[str, Any] | None = None

    @classmethod
    def from_discovery(cls, products: DiscoveryProducts) -> "CandidateSelection":
        """Construct the handoff from live discovery products.

        Args:
            products: Typed products returned by the discovery stage.

        Returns:
            A shared candidate selection for downstream stages.
        """
        return cls(
            screening_database=products['screening_database'],
            top_catalysts=products['top_catalysts'],
            dft_candidates=products['dft_candidates'])

    @classmethod
    def from_restart(cls, selected: SelectedCandidates) -> "CandidateSelection":
        """Construct the same handoff from a persisted screening artifact.

        Args:
            selected: Candidate routes reconstructed by the persistence adapter.

        Returns:
            A shared candidate selection for downstream stages.
        """
        return cls(
            screening_database=selected['screening_database'],
            top_catalysts=selected['top_catalysts'],
            dft_candidates=selected['dft_candidates'],
            admissibility=selected.get('admissibility'))


@dataclass(slots=True)
class PipelineRunContext:
    """Typed transient products available during one orchestrator process.

    Persisted evidence belongs in ``PipelineStateDocument``; large tabular
    products stay here and are reconstructed through the candidate loader on
    phase-only restarts.
    """

    candidates: CandidateSelection | None = None
