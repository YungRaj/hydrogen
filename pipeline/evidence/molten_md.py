"""Deterministic MLIP molecular-dynamics snapshot generation for liquid alloys."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
from ase import Atoms
from ase import units
from ase.md.langevin import Langevin
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution


@dataclass(frozen=True)
class LiquidMDProtocol:
    """Controls equilibration and decorrelated snapshot collection."""

    temperature_K: float
    timestep_fs: float = 1.0
    equilibration_steps: int = 2_000
    stride_steps: int = 200
    snapshots: int = 10
    friction_per_fs: float = 0.01
    seed: int = 0

    def validate(self) -> None:
        """Reject physically or numerically invalid MD settings."""
        if self.temperature_K <= 0.0:
            raise ValueError('temperature_K must be positive')
        if self.timestep_fs <= 0.0 or self.friction_per_fs <= 0.0:
            raise ValueError('timestep and friction must be positive')
        if self.equilibration_steps < 0 or self.stride_steps <= 0 or self.snapshots <= 0:
            raise ValueError('MD step counts must be non-negative and snapshots positive')


def collect_liquid_snapshots(
    atoms: Atoms,
    protocol: LiquidMDProtocol,
) -> tuple[list[Atoms], dict[str, Any]]:
    """Equilibrate and collect independent liquid snapshots with one calculator.

    Args:
        atoms: Initial liquid structure with an attached calculator.
        protocol: Temperature, integration, and sampling settings.

    Returns:
        Collected snapshots and the resolved protocol metadata.

    The input is copied and never mutated. Each returned snapshot retains the
    calculator so downstream relaxation/NEB code can consume it directly.
    """
    protocol.validate()
    if atoms.calc is None:
        raise ValueError('liquid MD requires an attached calculator')
    if len(atoms) == 0:
        raise ValueError('liquid MD requires at least one atom')

    system = atoms.copy()
    calculator = atoms.calc
    system.calc = calculator
    rng = np.random.RandomState(protocol.seed)
    MaxwellBoltzmannDistribution(
        system,
        temperature_K=protocol.temperature_K,
        rng=rng,
    )
    dynamics = Langevin(
        system,
        timestep=protocol.timestep_fs * units.fs,
        temperature_K=protocol.temperature_K,
        friction=protocol.friction_per_fs / units.fs,
        rng=rng,
    )
    dynamics.run(protocol.equilibration_steps)

    collected: list[Atoms] = []
    for _ in range(protocol.snapshots):
        dynamics.run(protocol.stride_steps)
        snapshot = system.copy()
        snapshot.calc = calculator
        collected.append(snapshot)
    return collected, {
        'temperature_K': protocol.temperature_K,
        'timestep_fs': protocol.timestep_fs,
        'equilibration_steps': protocol.equilibration_steps,
        'stride_steps': protocol.stride_steps,
        'snapshots': protocol.snapshots,
        'seed': protocol.seed,
    }
