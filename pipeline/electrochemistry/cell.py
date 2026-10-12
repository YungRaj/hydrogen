"""Analytic first-tier electrochemical cell balance."""

from __future__ import annotations

from dataclasses import dataclass

from pipeline.electrochemistry.transport import concentration_overpotential_V


@dataclass(frozen=True)
class CellOperatingPoint:
    """Inputs for one diagnostic current/voltage point."""

    reversible_voltage_V: float
    anode_overpotential_V: float
    cathode_overpotential_V: float
    current_A: float
    limiting_current_A: float
    resistance_ohm: float
    temperature_K: float
    electrons: int = 2


def cell_voltage(point: CellOperatingPoint) -> float:
    """Return cell voltage including ohmic and concentration losses."""
    if point.current_A < 0.0 or point.resistance_ohm < 0.0:
        raise ValueError('current and resistance must be nonnegative')
    concentration = concentration_overpotential_V(
        point.current_A, point.limiting_current_A, point.temperature_K,
        point.electrons,
    )
    return (
        point.reversible_voltage_V
        + point.anode_overpotential_V
        + abs(point.cathode_overpotential_V)
        + point.current_A * point.resistance_ohm
        + concentration
    )
