"""PEM fuel-cell polarization and system-level stack models."""

from typing import Any, Mapping


def composite_score(result: Mapping[str, Any]) -> float:
    """Rank PEMFC results by efficiency and power per volt of overpotential.

    Args:
        result: One PEMFC result record.

    Returns:
        Higher-is-better composite performance score.
    """
    efficiency = float(result.get('efficiency_at_peak', 0.0))
    power = float(result.get('peak_power_W_cm2', 0.0))
    overpotential = max(float(result.get('orr_overpotential_V', 0.4)), 0.01)
    return efficiency * power / overpotential
