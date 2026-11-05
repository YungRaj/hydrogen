"""DFT validation helpers and CHE ZPE/entropy corrections for fuel-cell ORR catalysts."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from pipeline.validation.dft_validator import (
    parse_convergence,
    parse_total_energy,
)


def converged_energy(output_file: Path | str) -> Optional[float]:
    """Return the final QE energy only for a cleanly converged calculation.

    Args:
        output_file: Filesystem location used for output file.

    Returns:
        Computed `Optional[float]` result.
    """
    path = str(output_file)
    if not parse_convergence(path, require_ionic=True):
        return None
    return parse_total_energy(path)


# Free energy corrections (ZPE + entropy at 298 K, 1 bar)
# From Nørskov et al., J. Phys. Chem. B 108, 17886 (2004)
ZPE_CORRECTIONS_eV: dict[str, float] = {
    'OH*': 0.35,
    'O*': 0.05,
    'OOH*': 0.40,
    'H2O_g': 0.56,
    'H2_g': 0.27,
    'O2_g': 0.10,
}

TS_CORRECTIONS_eV: dict[str, float] = {
    'OH*': 0.07,
    'O*': 0.00,
    'OOH*': 0.10,
    'H2O_g': 0.67,
    'H2_g': 0.40,
    'O2_g': 0.63,
}
