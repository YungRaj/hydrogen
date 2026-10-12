"""Small deterministic low/nominal/high sensitivity helper."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import Any


def sensitivity_band(
    evaluate: Callable[[Mapping[str, float]], float],
    cases: Mapping[str, Mapping[str, float]],
) -> dict[str, Any]:
    """Evaluate named low/nominal/high cases without fitting parameters.

    Args:
        evaluate: Callable that maps one parameter mapping to a scalar.
        cases: Named parameter mappings, typically low/nominal/high.

    Returns:
        Named values and their minimum, maximum, and nominal entries.
    """
    if not cases:
        raise ValueError('sensitivity cases cannot be empty')
    values = {name: float(evaluate(parameters)) for name, parameters in cases.items()}
    return {
        'values': values,
        'minimum': min(values.values()),
        'maximum': max(values.values()),
        'nominal': values.get('nominal'),
    }
