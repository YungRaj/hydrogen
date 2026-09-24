"""Composable runtime services for pipeline orchestration.

Scientific stages should receive data and return evidence. This module isolates
the process-level effects needed by the legacy CLI so orchestration can be
tested with an in-memory state store and deterministic clock.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
import os
import time
from pathlib import Path
from typing import Any, Callable, Mapping, Protocol, Sequence

from pipeline.data_models.stages import (
    DFTProducts, DFTState, DiscoveryProducts, DiscoveryState,
    FuelCellProducts, FuelCellState, ReactorBatchProducts, ReactorBatchState,
    ReportProducts, ReportState, SelectedCandidates, VQEProducts, VQEState,
)
from pipeline.data_models.campaigns import (
    PipelineStateDocument, validate_pipeline_state)
from pipeline.stages.contracts import StageOutcome


PipelineState = PipelineStateDocument


class DiscoveryStage(Protocol):
    """Callable contract for a replaceable discovery implementation."""

    def __call__(self, *, initial_samples: int, leaf_size: int,
                 max_leaves: int | None, top_k_reactor: int,
                 top_k_dft: int
                 ) -> StageOutcome[DiscoveryState, DiscoveryProducts]:
        """Run discovery and return its typed state and products.

        Args:
            initial_samples: Initial high-fidelity sampling allocation.
            leaf_size: Maximum indexed search leaf size.
            max_leaves: Optional bounded number of leaves to process.
            top_k_reactor: Number of reactor candidates to route.
            top_k_dft: Number of DFT candidates to route.

        Returns:
            Typed discovery state and in-process products.
        """
        ...


class ReactorBatchStage(Protocol):
    """Callable contract for a replaceable reactor-batch implementation."""

    def __call__(
            self, candidates: Any, *, temperatures: Sequence[float],
            reactor_types: Sequence[str], pathway_mode: str,
            multiphysics_results_dir: str, allow_mock_inputs: bool,
            judge_catalyst: str | None = None,
            headline_t_min: float = 1200.0,
            headline_t_max: float | None = None
            ) -> StageOutcome[ReactorBatchState, ReactorBatchProducts]:
        """Run candidate reactor cases and return typed batch evidence.

        Args:
            candidates: Candidate table or compatible batch.
            temperatures: Reactor temperatures in kelvin.
            reactor_types: Physical reactor implementations to execute.
            pathway_mode: Methane-conversion pathway identifier.
            multiphysics_results_dir: Validated full-physics artifact root.
            allow_mock_inputs: Explicit test-only mock permission.
            judge_catalyst: Optional named reference catalyst.
            headline_t_min: Minimum headline temperature in kelvin.
            headline_t_max: Optional maximum headline temperature in kelvin.

        Returns:
            Typed reactor-batch state and detailed results.
        """
        ...


class DFTStage(Protocol):
    """Callable contract for a replaceable DFT implementation."""

    def __call__(self, candidates: Any, *, top_k: int, execute_dft: bool,
                 name_prefix: str = 'dft_cat',
                 error_sink: Callable[[str], None] | None = None
                 ) -> StageOutcome[DFTState, DFTProducts]:
        """Run candidate validation and return typed DFT evidence.

        Args:
            candidates: Candidate table or sequence to validate.
            top_k: Maximum number of candidates to process.
            execute_dft: Whether to launch the production solver.
            name_prefix: Stable prefix for calculation identifiers.
            error_sink: Optional failure-reporting callback.

        Returns:
            Typed DFT execution state, results, and failures.
        """
        ...


class VQEStage(Protocol):
    """Callable contract for a replaceable VQE implementation."""

    def __call__(self, *, top_k: int, execute_quantum: bool
                 ) -> StageOutcome[VQEState, VQEProducts]:
        """Run transition-state validation and return typed VQE evidence.

        Args:
            top_k: Maximum number of candidates to process.
            execute_quantum: Whether to use the configured quantum backend.

        Returns:
            Typed VQE execution state and candidate results.
        """
        ...


class FuelCellStage(Protocol):
    """Callable contract for a replaceable fuel-cell implementation."""

    def __call__(self, *, top_k_pemfc: int, stack_cells: int
                 ) -> StageOutcome[FuelCellState, FuelCellProducts]:
        """Run fuel-cell evaluation and return typed products.

        Args:
            top_k_pemfc: Maximum cathode candidates entering PEMFC modeling.
            stack_cells: Number of cells in the modeled stack.

        Returns:
            Typed fuel-cell state and detailed products.
        """
        ...


class ReportStage(Protocol):
    """Callable contract for a replaceable report implementation."""

    def __call__(self, pipeline_state: Mapping[str, Any]
                 ) -> StageOutcome[ReportState, ReportProducts]:
        """Render the persisted pipeline state into a report artifact.

        Args:
            pipeline_state: Persisted state from completed stages.

        Returns:
            Typed report state and artifact path.
        """
        ...


class CandidateLoader(Protocol):
    """Callable contract for restoring discovery products from disk."""

    def __call__(self, path: str | Path, *, top_k_reactor: int,
                 top_k_dft: int) -> SelectedCandidates | None:
        """Load and route saved candidates, or return ``None`` if absent.

        Args:
            path: Screening-table artifact path.
            top_k_reactor: Number of reactor candidates to restore.
            top_k_dft: Number of validation candidates to restore.

        Returns:
            Typed restored candidates, or ``None`` when the file is absent.
        """
        ...


@dataclass(frozen=True)
class PipelineRuntime:
    """Side-effect boundary used by the master pipeline coordinator."""

    clock: Callable[[], float]
    banner: Callable[[str], None]
    load_state: Callable[[], PipelineState]
    save_state: Callable[[Mapping[str, Any]], None]
    select_pathway_mode: Callable[[str], None]


@dataclass(frozen=True)
class PipelineComponents:
    """Replaceable scientific stage graph used by the coordinator."""

    discovery: DiscoveryStage
    reactor_batch: ReactorBatchStage
    dft: DFTStage
    vqe: VQEStage
    fuel_cell: FuelCellStage
    report: ReportStage
    load_candidates: CandidateLoader


def default_pipeline_runtime() -> PipelineRuntime:
    """Build the filesystem/environment-backed runtime used by the CLI.

    Returns:
        Computed `PipelineRuntime` result.
    """
    from pipeline.utils import load_json, print_banner, save_json

    def load() -> PipelineState:
        return validate_pipeline_state(load_json('pipeline_state.json') or {})

    def save(value: Mapping[str, Any]) -> None:
        save_json(dict(value), 'pipeline_state.json')

    def select_mode(mode: str) -> None:
        # Compatibility boundary for modules that still consume this setting.
        os.environ['PYROLYSIS_MODE'] = mode

    return PipelineRuntime(
        clock=time.time, banner=print_banner, load_state=load,
        save_state=save, select_pathway_mode=select_mode)


def default_pipeline_components() -> PipelineComponents:
    """Bind the production stage graph lazily.

    Returns:
        Computed `PipelineComponents` result.
    """
    from pipeline.stages.candidate_io import load_selected_candidates
    from pipeline.stages.dft import run_dft_stage
    from pipeline.stages.discovery import run_discovery_stage
    from pipeline.stages.fuel_cell import run_fuel_cell_stage
    from pipeline.stages.reactor_batch import run_reactor_batch_stage
    from pipeline.stages.report import run_report_stage
    from pipeline.stages.vqe import run_vqe_stage
    return PipelineComponents(
        discovery=run_discovery_stage, reactor_batch=run_reactor_batch_stage,
        dft=run_dft_stage, vqe=run_vqe_stage,
        fuel_cell=run_fuel_cell_stage, report=run_report_stage,
        load_candidates=load_selected_candidates)


class MemoryStateStore:
    """Small test/runtime adapter that never touches the filesystem."""

    def __init__(self, initial: Mapping[str, Any] | None = None):
        self.state = deepcopy(dict(initial or {}))
        self.history: list[dict] = []

    def load(self) -> PipelineState:
        """Load and validate a serialized surrogate model.

        Returns:
            A validated model reconstructed from disk.
        """
        return deepcopy(self.state)

    def save(self, value: Mapping[str, Any]) -> None:
        """Atomically save the surrogate model.

        Args:
            value: Pipeline state value to copy or persist.
        """
        self.state = deepcopy(dict(value))
        self.history.append(deepcopy(dict(value)))
