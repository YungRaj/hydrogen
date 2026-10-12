"""Dimensioned, non-pruning NTEC mechanics estimates."""

from __future__ import annotations

import math


def collision_events_per_volume_s(
    shear_rate_s: float, particle_diameter_m: float, number_density_m3: float
) -> float:
    """Estimate homogeneous-suspension collision events per volume per second."""
    if min(shear_rate_s, particle_diameter_m, number_density_m3) <= 0.0:
        raise ValueError('shear rate, diameter, and number density must be positive')
    return (4.0 / 3.0) * shear_rate_s * particle_diameter_m**3 * number_density_m3**2


def impact_energy_J(
    effective_mass_kg: float, shear_rate_s: float, particle_diameter_m: float
) -> float:
    """Estimate impact energy from the shear velocity scale ``u = G d``."""
    if min(effective_mass_kg, shear_rate_s, particle_diameter_m) <= 0.0:
        raise ValueError('mass, shear rate, and diameter must be positive')
    return 0.5 * effective_mass_kg * (shear_rate_s * particle_diameter_m) ** 2


def detachment_probability(
    adhesion_energy_J_m2: float,
    impact_energy_J_value: float,
    contact_area_m2: float,
) -> float:
    """Return detachment from work of adhesion and impact energy density."""
    if adhesion_energy_J_m2 < 0.0 or impact_energy_J_value <= 0.0 or contact_area_m2 <= 0.0:
        raise ValueError('adhesion, impact energy, and contact area are invalid')
    impact_energy_density = impact_energy_J_value / contact_area_m2
    return math.exp(-adhesion_energy_J_m2 / impact_energy_density)
