"""Interfaces for portable external-solver discovery and execution."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable


@dataclass(frozen=True)
class SolverExecutionServices:
    """External process boundary consumed by the multiphysics coordinator."""

    preflight: Callable[[str], dict]
    execute: Callable
    openfoam_version: Callable[[str], str]
    fenics_command: Callable
    cantera_version: Callable[[], str]
