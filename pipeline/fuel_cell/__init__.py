"""PEM fuel-cell polarization and system-level stack models."""

from __future__ import annotations

from collections.abc import Mapping


def composite_score(result: Mapping[str, object]) -> float:
    """Rank PEMFC results by efficiency and power per volt of overpotential.

    Args:
        result: One PEMFC result record.

    Returns:
        Higher-is-better composite performance score.
    """
    raw_eff = result.get('efficiency_at_peak', 0.0)
    raw_power = result.get('peak_power_W_cm2', 0.0)
    raw_eta = result.get('orr_overpotential_V', 0.4)
    efficiency = float(raw_eff) if isinstance(raw_eff, (int, float)) else 0.0
    power = float(raw_power) if isinstance(raw_power, (int, float)) else 0.0
    eta_val = float(raw_eta) if isinstance(raw_eta, (int, float)) else 0.4
    overpotential = max(eta_val, 0.01)
    return efficiency * power / overpotential

