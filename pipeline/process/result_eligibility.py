"""Shared baseline for a reactor result that may be ranked or costed.

Scorecard, joint-band search, best-condition selection, the report solids
table, and hydrogen-cost estimates used to apply different subsets of the
same checks. An X>X_eq row could lose the scorecard and still become
``best_condition`` and a $/kg estimate.

Consumers may add filters on top (solids surface, headline T, named
judge, a finite conversion for ranking). They share this baseline:
successful execution, not mock, no equilibrium overshoot, finite
conversion when one was reported, and a passing carbon balance when one
was reported. Missing ``CH4_conversion`` (joint-band yield/lifetime
rows) or ``carbon_balance_ok`` (MMBCR, NTEC, legacy JSON) is not a
failure.
"""

from __future__ import annotations

import math
from typing import Any, Optional


def finite_number(value: Any) -> Optional[float]:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def is_usable_result(record: dict) -> bool:
    """True if *record* may enter a rank, table sort, or TEA estimate."""
    if not isinstance(record, dict):
        return False
    if record.get('status', 'complete') != 'complete':
        return False
    if record.get('mock'):
        return False
    if bool(record.get('exceeds_equilibrium')):
        return False
    if 'CH4_conversion' in record and finite_number(record.get('CH4_conversion')) is None:
        return False
    if 'carbon_balance_ok' in record and record.get('carbon_balance_ok') is False:
        return False
    return True


def is_rankable_result(record: dict) -> bool:
    """Usable baseline plus a finite conversion — used to pick best_condition."""
    return is_usable_result(record) and finite_number(record.get('CH4_conversion')) is not None


def usable_results(records) -> list:
    return [record for record in records if is_usable_result(record)]


def rankable_results(records) -> list:
    return [record for record in records if is_rankable_result(record)]
