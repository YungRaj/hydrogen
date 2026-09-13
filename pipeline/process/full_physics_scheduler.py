"""Coverage-preserving allocation of expensive full-physics calculations."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable


@dataclass(frozen=True)
class FullPhysicsRequest:
    """One proposed calculation with dimensionless, calibrated priorities."""

    case_id: str
    region: str
    expected_improvement: float = 0.0
    uncertainty: float = 0.0
    calibration_error: float = 0.0
    disagreement: float = 0.0
    unproductive_history: float = 0.0

    def validate(self) -> None:
        if not self.case_id or not self.region:
            raise ValueError('case and region identities are required')
        for name in ('expected_improvement', 'uncertainty',
                     'calibration_error', 'disagreement',
                     'unproductive_history'):
            value = float(getattr(self, name))
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f'{name} must be finite and within [0, 1]')

    def priority(self) -> float:
        # Repeatedly unproductive regions lose discretionary priority only.
        benefit = (0.30 * self.expected_improvement + 0.25 * self.uncertainty +
                   0.25 * self.disagreement + 0.20 * self.calibration_error)
        return benefit * (1.0 - 0.75 * self.unproductive_history)


def schedule_full_physics_cases(
        requests: Iterable[FullPhysicsRequest], *, total_budget: int,
        minimum_per_region: int = 1) -> list[dict]:
    """Reserve regional coverage, then allocate remaining budget by priority."""
    values = list(requests)
    if not isinstance(total_budget, int) or total_budget < 1:
        raise ValueError('total_budget must be a positive integer')
    if not isinstance(minimum_per_region, int) or minimum_per_region < 1:
        raise ValueError('minimum_per_region must be a positive integer')
    by_region: dict[str, list[FullPhysicsRequest]] = {}
    seen = set()
    for request in values:
        if not isinstance(request, FullPhysicsRequest):
            raise TypeError('requests must be FullPhysicsRequest values')
        request.validate()
        if request.case_id in seen:
            raise ValueError('full-physics case identities must be unique')
        seen.add(request.case_id)
        by_region.setdefault(request.region, []).append(request)
    required = sum(min(minimum_per_region, len(rows))
                   for rows in by_region.values())
    if total_budget < required:
        raise ValueError('budget cannot satisfy fixed regional coverage')
    ordered = lambda rows: sorted(rows, key=lambda row: (-row.priority(), row.case_id))
    chosen: list[tuple[FullPhysicsRequest, str]] = []
    selected_ids = set()
    for region in sorted(by_region):
        for request in ordered(by_region[region])[:minimum_per_region]:
            chosen.append((request, 'fixed_regional_coverage'))
            selected_ids.add(request.case_id)
    remaining = [row for row in values if row.case_id not in selected_ids]
    for request in ordered(remaining)[:max(0, total_budget - len(chosen))]:
        chosen.append((request, 'adaptive_priority'))
    return [{
        'case_id': request.case_id,
        'region': request.region,
        'priority': request.priority(),
        'selection_reason': reason,
        'candidate_exclusion_authorized': False,
    } for request, reason in chosen]
