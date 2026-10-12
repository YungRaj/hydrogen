"""Analytic first-tier electrolyte transport relations."""

from __future__ import annotations

import math


def ohmic_resistance_ohm(
    gap_m: float, conductivity_S_m: float, area_m2: float
) -> float:
    """Return electrolyte resistance for a uniform gap."""
    if min(gap_m, conductivity_S_m, area_m2) <= 0.0:
        raise ValueError('gap, conductivity, and area must be positive')
    return gap_m / (conductivity_S_m * area_m2)


def limiting_current_A(
    electrons: int, diffusivity_m2_s: float, concentration_mol_m3: float,
    area_m2: float, boundary_layer_m: float,
) -> float:
    """Return the diffusion-limited current from a planar layer."""
    if electrons <= 0 or min(diffusivity_m2_s, concentration_mol_m3,
                              area_m2, boundary_layer_m) <= 0.0:
        raise ValueError('transport inputs must be positive')
    return electrons * 96485.33212 * diffusivity_m2_s * concentration_mol_m3 \
        * area_m2 / boundary_layer_m


def concentration_overpotential_V(
    current_A: float, limiting_current_A_value: float, temperature_K: float,
    electrons: int = 2,
) -> float:
    """Return the concentration overpotential before the limiting current."""
    if limiting_current_A_value <= 0.0 or temperature_K <= 0.0 or electrons <= 0:
        raise ValueError('limiting current, temperature, and electrons are invalid')
    fraction = current_A / limiting_current_A_value
    if not 0.0 <= fraction < 1.0:
        raise ValueError('current must satisfy 0 <= I < I_lim')
    return 8.314462618 * temperature_K / (electrons * 96485.33212) * math.log(
        1.0 / (1.0 - fraction)
    )
