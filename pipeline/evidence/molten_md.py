"""Deterministic MLIP molecular-dynamics protocol contract for liquid alloys."""

from __future__ import annotations

from dataclasses import dataclass


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
        if (
            self.equilibration_steps < 0
            or self.stride_steps <= 0
            or self.snapshots <= 0
        ):
            raise ValueError(
                'MD step counts must be non-negative and snapshots positive'
            )
