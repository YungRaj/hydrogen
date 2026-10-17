"""Static shapes for values exchanged by the six orchestration stages.

These contracts intentionally retain dictionaries and pandas objects at
runtime, so existing artifacts and scientific implementations do not change.
They remove the *type erasure* at stage boundaries and provide an incremental
path toward richer domain objects.
"""

from __future__ import annotations

from pathlib import Path
from dataclasses import dataclass
from typing import Protocol, Sequence, TypedDict

from pipeline.data_models.quantum import DFTResult, VQEResult
from pipeline.data_models.reactors import ReactorResult
from pipeline.data_models.fuel_cells import EmptyStackResult, PEMFCResult, StackResult
from pipeline.data_models.discovery import AdmissibilitySummary, EquilibriumSummary
from pipeline.data_models.reactors import SolidsScorecard


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
    best_coking: float | str
    metric_eligible_count: int
    search_strategy: str
    admissibility: AdmissibilitySummary


class DiscoveryState(_DiscoveryStateOptional):
    """Persistable summary emitted by branch-and-bound discovery."""

    pareto_size: int
    total_evaluated: int
    valid_count: int
    top_catalysts_count: int
    dft_resolution_count: int
    candidate_dispositions: dict[str, int]


class DiscoveryProducts(TypedDict):
    """In-process values routed from discovery to later stages."""

    design_space_sizes: dict[str, int]
    pareto_genomes: Sequence[object]
    screening_database: TableLike
    top_catalysts: TableLike
    dft_candidates: CandidateBatch


class _ReactorBatchStateOptional(TypedDict, total=False):
    n_simulations: int
    catalysts_attempted: int
    catalysts_simulated: int
    complete_sweeps: int
    partial_sweeps: int
    failed_sweeps: int
    equilibrium_check: EquilibriumSummary
    solids_scorecard: SolidsScorecard
    best_conversion: float | None
    best_conversion_scope: str
    candidate_champion_conversion: float | None
    candidate_champion_catalyst: str | None
    best_h2_cost_usd_kg: float
    mmbcr_max_conversion: float | None
    industrial_viability: dict[str, int]
    skipped: bool
    reason: str


class ReactorBatchState(_ReactorBatchStateOptional):
    """Persistable summary emitted by a reactor batch."""


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


class _VQEStateOptional(TypedDict, total=False):
    candidates: list[str]
    n_missing_hamiltonians: int


class VQEState(_VQEStateOptional):
    """Persistable VQE execution counts."""

    n_vqe_runs: int


class _VQEProductsOptional(TypedDict, total=False):
    missing_hamiltonians: list[str]


class VQEProducts(_VQEProductsOptional):
    """Candidate-specific VQE results."""

    vqe_results: list[VQEResult]


class _FuelCellStateOptional(TypedDict, total=False):
    n_pemfc_simulations: int
    best_power_W_cm2: float
    best_efficiency: float
    min_overpotential_V: float
    best_overpotential_V: float
    best_overpotential_gc_V: float
    pourbaix_acid_stable_count: int
    best_catalyst: str
    best_membrane: str
    stack_net_kW: float
    stack_efficiency: float


class FuelCellState(_FuelCellStateOptional):
    """Persistable cathode and PEMFC summary."""

    n_cathodes_screened: int
    n_valid: int



class FuelCellProducts(TypedDict):
    """Detailed fuel-cell screening, cell, and stack products."""

    cathode_database: TableLike
    valid_cathodes: TableLike
    pemfc_results: list[PEMFCResult]
    stack_result: StackResult | EmptyStackResult


class ReportState(TypedDict):
    """Persistable report-stage summary."""

    report_path: str


class ReportProducts(TypedDict):
    """Generated report artifact path."""

    report_path: Path


class _SelectedCandidatesOptional(TypedDict, total=False):
    admissibility: AdmissibilitySummary


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
    admissibility: AdmissibilitySummary | None = None

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
            dft_candidates=products['dft_candidates'],
        )

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
            admissibility=selected.get('admissibility'),
        )


@dataclass(slots=True)
class PipelineRunContext:
    """Typed transient products available during one orchestrator process.

    Persisted evidence belongs in ``PipelineStateDocument``; large tabular
    products stay here and are reconstructed through the candidate loader on
    phase-only restarts.
    """

    candidates: CandidateSelection | None = None
