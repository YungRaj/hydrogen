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
from typing import Callable, Mapping


@dataclass(frozen=True)
class PipelineRuntime:
    """Side-effect boundary used by the master pipeline coordinator."""

    clock: Callable[[], float]
    banner: Callable[[str], None]
    load_state: Callable[[], dict]
    save_state: Callable[[Mapping], None]
    select_pathway_mode: Callable[[str], None]


@dataclass(frozen=True)
class PipelineComponents:
    """Replaceable scientific stage graph used by the coordinator."""

    discovery: Callable
    reactor_batch: Callable
    dft: Callable
    vqe: Callable
    fuel_cell: Callable
    report: Callable
    load_candidates: Callable


def default_pipeline_runtime() -> PipelineRuntime:
    """Build the filesystem/environment-backed runtime used by the CLI.

    Returns:
        Computed `PipelineRuntime` result.
    """
    from pipeline.common.utils import load_json, print_banner, save_json

    def load() -> dict:
        return load_json('pipeline_state.json') or {}

    def save(value: Mapping) -> None:
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

    def __init__(self, initial: Mapping | None = None):
        self.state = deepcopy(dict(initial or {}))
        self.history: list[dict] = []

    def load(self) -> dict:
        """Load and validate a serialized surrogate model.

        Returns:
            A validated model reconstructed from disk.
        """
        return deepcopy(self.state)

    def save(self, value: Mapping) -> None:
        """Atomically save the surrogate model.

        Args:
            value: Pipeline state value to copy or persist.
        """
        self.state = deepcopy(dict(value))
        self.history.append(deepcopy(dict(value)))
