"""Descriptor-only ORR volcano cross-check.

This is intentionally independent of the screening model.  It provides the
standard one-dimensional OH binding-energy volcano and a CSV export that can
be consumed by CatMAP when that optional dependency is available.
"""

from __future__ import annotations

import math
from typing import Iterable, Mapping


def volcano_rate_log10(dG_OH_eV: float, *, temperature_K: float = 298.15) -> float:
    """Return a dimensionless log10 relative rate from a simple OH volcano.

    The optimum is 0.86 eV (the conventional ORR OH descriptor); deviations
    are penalized symmetrically as a descriptor sanity check, not as a fitted
    screening score.

    Args:
        dG_OH_eV: Hydroxyl adsorption free energy.
        temperature_K: Absolute temperature.

    Returns:
        Base-10 logarithm of the relative volcano rate.
    """
    if temperature_K <= 0 or not math.isfinite(dG_OH_eV):
        raise ValueError('finite descriptor and positive temperature required')
    # Universal ORR Sabatier volcano apex at dG_OH* = (4.92 - 3.20) / 2 = 0.86 eV,
    # where the OOH*-to-OH* scaling offset (3.20 eV) splits the 4.92 eV 4e- reaction
    # energy equally between step 1 (O2 -> OOH*) and step 4 (OH* -> H2O)
    # (Norskov et al., J. Phys. Chem. B 2004; Viswanathan et al., ACS Catal. 2012).
    activation_eV = abs(float(dG_OH_eV) - 0.86)
    # Base-10 relative Arrhenius rate: log10(k / k_max) = -|dG_OH* - 0.86| / (k_B * T * ln(10)),
    # where k_B = 8.617333262e-5 eV/K (CODATA 2018).
    return -activation_eV / (8.617333262e-5 * temperature_K * math.log(10.0))


def rank_descriptors(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    """Annotate and rank rows without changing their source values.

    Returns:
        New rows ordered by descending volcano rate.
    """
    ranked = []
    for row in rows:
        value = float(row['dG_OH_eV'])
        ranked.append({**row, 'volcano_log10_rate': volcano_rate_log10(value)})
    return sorted(ranked, key=lambda row: float(row['volcano_log10_rate']), reverse=True)


