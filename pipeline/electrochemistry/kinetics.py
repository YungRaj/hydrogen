"""Small, explicit electrochemical kinetics equations for screening diagnostics."""

from __future__ import annotations

import math

K_B_EV_K = 8.617333262145e-5


def exchange_current_density(
    descriptor_eV: float,
    optimum_eV: float,
    reference_A_cm2: float,
    temperature_K: float,
) -> float:
    """Return a provisional Sabatier exchange-current estimate."""
    if reference_A_cm2 <= 0.0 or temperature_K <= 0.0:
        raise ValueError('reference current and temperature must be positive')
    if not all(math.isfinite(value) for value in (descriptor_eV, optimum_eV)):
        raise ValueError('descriptor values must be finite')
    return reference_A_cm2 * math.exp(
        -abs(descriptor_eV - optimum_eV) / (K_B_EV_K * temperature_K)
    )


def butler_volmer_current_density(
    exchange_current_A_cm2: float,
    overpotential_V: float,
    temperature_K: float,
    symmetry_factor: float = 0.5,
    electrons: int = 2,
) -> float:
    """Return signed Butler–Volmer current density in A/cm²."""
    if exchange_current_A_cm2 < 0.0 or temperature_K <= 0.0 or electrons <= 0:
        raise ValueError('current, temperature, and electron count are invalid')
    if not 0.0 < symmetry_factor < 1.0:
        raise ValueError('symmetry_factor must lie strictly between zero and one')
    factor = electrons * 96485.33212 / (8.314462618 * temperature_K)
    return exchange_current_A_cm2 * (
        math.exp((1.0 - symmetry_factor) * factor * overpotential_V)
        - math.exp(-symmetry_factor * factor * overpotential_V)
    )
