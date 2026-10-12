"""Dimensioned, non-pruning NTEC mechanics estimates."""

from __future__ import annotations

import math


def collision_events_per_volume_s(
    shear_rate_s: float, particle_diameter_m: float, number_density_m3: float
) -> float:
    """Estimate homogeneous-suspension collision events per volume per second.

    Args:
        shear_rate_s: Shear rate in inverse seconds.
        particle_diameter_m: Particle diameter in metres.
        number_density_m3: Particle number density per cubic metre.

    Returns:
        Collision events per cubic metre per second.
    """
    if min(shear_rate_s, particle_diameter_m, number_density_m3) <= 0.0:
        raise ValueError('shear rate, diameter, and number density must be positive')
    return (4.0 / 3.0) * shear_rate_s * particle_diameter_m**3 * number_density_m3**2


def impact_energy_J(
    effective_mass_kg: float, shear_rate_s: float, particle_diameter_m: float
) -> float:
    """Estimate impact energy from the shear velocity scale ``u = G d``.

    Args:
        effective_mass_kg: Effective colliding mass.
        shear_rate_s: Shear rate in inverse seconds.
        particle_diameter_m: Particle diameter in metres.

    Returns:
        Impact energy in joules.
    """
    if min(effective_mass_kg, shear_rate_s, particle_diameter_m) <= 0.0:
        raise ValueError('mass, shear rate, and diameter must be positive')
    return 0.5 * effective_mass_kg * (shear_rate_s * particle_diameter_m) ** 2


def detachment_probability(
    adhesion_energy_J_m2: float,
    impact_energy_J_value: float,
    contact_area_m2: float,
) -> float:
    """Return detachment from work of adhesion and impact energy density.

    Args:
        adhesion_energy_J_m2: Work of adhesion per contact area.
        impact_energy_J_value: Impact energy.
        contact_area_m2: Effective contact area.

    Returns:
        Detachment probability between zero and one.
    """
    if adhesion_energy_J_m2 < 0.0 or impact_energy_J_value <= 0.0 or contact_area_m2 <= 0.0:
        raise ValueError('adhesion, impact energy, and contact area are invalid')
    impact_energy_density = impact_energy_J_value / contact_area_m2
    return math.exp(-adhesion_energy_J_m2 / impact_energy_density)
