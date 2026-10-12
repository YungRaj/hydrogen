"""Optional pMuTT cross-check boundary for locked thermochemistry results."""

from __future__ import annotations

from typing import Any


def pmutt_status() -> dict[str, Any]:
    """Return availability metadata without changing the production mechanism."""
    try:
        import pmutt
    except ImportError:
        return {'available': False, 'package': 'pymutt', 'version': None}
    return {
        'available': True,
        'package': 'pymutt',
        'version': str(getattr(pmutt, '__version__', 'unknown')),
    }


def compare_gas_reference(
    *,
    current_free_energy_eV: float,
    pmutt_free_energy_eV: float,
    tolerance_eV: float = 0.02,
) -> dict[str, Any]:
    """Compare one pMuTT value with the locked mechanism without overwriting it."""
    if tolerance_eV < 0.0:
        raise ValueError('tolerance_eV must be non-negative')
    delta = float(pmutt_free_energy_eV) - float(current_free_energy_eV)
    return {
        'current_free_energy_eV': float(current_free_energy_eV),
        'pmutt_free_energy_eV': float(pmutt_free_energy_eV),
        'difference_eV': delta,
        'within_tolerance': abs(delta) <= tolerance_eV,
        'tolerance_eV': float(tolerance_eV),
    }
