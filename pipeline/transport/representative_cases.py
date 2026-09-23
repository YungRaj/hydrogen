"""Deterministic representative-case designs for full multiphysics solvers."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from typing import Mapping, Sequence

import numpy as np


@dataclass(frozen=True)
class ParameterRange:
    """One explicit unit-bearing physical range."""

    minimum: float
    maximum: float
    scale: str = 'linear'

    def validate(self, name: str) -> None:
        """Reject invalid configuration before it reaches scientific execution.

        Args:
            name: Human-readable identifier used in diagnostics and output.
        """
        if (not math.isfinite(self.minimum) or not math.isfinite(self.maximum)
                or self.maximum <= self.minimum):
            raise ValueError(f'invalid parameter range: {name}')
        if self.scale not in {'linear', 'log'}:
            raise ValueError(f'unsupported parameter scale: {name}')
        if self.scale == 'log' and self.minimum <= 0:
            raise ValueError(f'log parameter must be positive: {name}')

    def interpolate(self, fraction: float) -> float:
        """Map a unit-interval coordinate into this parameter range.

        Args:
            fraction: Unit-interval coordinate to map into the parameter range.

        Returns:
            The physical parameter value corresponding to the unit coordinate.
        """
        if self.scale == 'log':
            return float(math.exp(
                math.log(self.minimum) + fraction *
                (math.log(self.maximum) - math.log(self.minimum))))
        return float(self.minimum + fraction * (self.maximum - self.minimum))


def _case_id(mode: str, reactor: str, features: Mapping[str, float]) -> str:
    payload = json.dumps({
        'mode': mode, 'reactor': reactor, 'features': dict(features)},
        sort_keys=True, separators=(',', ':')).encode()
    return hashlib.sha256(payload).hexdigest()[:20]


def design_representative_cases(
        *, pathway_mode: str, reactor_type: str,
        ranges: Mapping[str, ParameterRange], sample_count: int,
        anchors: Sequence[Mapping[str, float]] = (), random_seed: int = 0
        ) -> list[dict]:
    """Generate a Latin-hypercube design plus explicit regime anchors.

        Every sampled dimension occupies every one of ``sample_count`` strata once.
        Anchors are appended after validation and deduplicated by physical identity.

    Args:
        pathway_mode: Configured methane-conversion pathway.
        reactor_type: Physical reactor implementation identifier.
        ranges: Mapping supplying ranges.
        sample_count: Number of sample count to use.
        anchors: Mapping supplying anchors.
        random_seed: Seed controlling deterministic sampling or fitting.

    Returns:
        List of computed or validated records.
    """
    if not pathway_mode or not reactor_type:
        raise ValueError('pathway mode and reactor type are required')
    if not isinstance(sample_count, int) or sample_count < 2:
        raise ValueError('sample_count must be an integer of at least two')
    if not ranges:
        raise ValueError('at least one physical parameter range is required')
    names = tuple(sorted(ranges))
    for name in names:
        if not isinstance(name, str) or not name:
            raise ValueError('parameter names must be non-empty strings')
        ranges[name].validate(name)
    rng = np.random.default_rng(random_seed)
    columns = {}
    for name in names:
        strata = (np.arange(sample_count) + rng.random(sample_count)) / sample_count
        columns[name] = strata[rng.permutation(sample_count)]
    feature_rows = [
        {name: ranges[name].interpolate(columns[name][row]) for name in names}
        for row in range(sample_count)]
    for anchor in anchors:
        if set(anchor) != set(names):
            raise ValueError('every anchor must match the exact parameter schema')
        values = {}
        for name in names:
            value = float(anchor[name])
            if (not math.isfinite(value) or value < ranges[name].minimum or
                    value > ranges[name].maximum):
                raise ValueError(f'anchor lies outside range: {name}')
            values[name] = value
        feature_rows.append(values)
    unique = {}
    for features in feature_rows:
        identity = _case_id(pathway_mode, reactor_type, features)
        unique.setdefault(identity, {
            'case_id': identity, 'pathway_mode': pathway_mode,
            'reactor_type': reactor_type, 'features': features,
            'design_role': ('regime_anchor' if features in anchors
                            else 'stratified_coverage')})
    return list(unique.values())


def assign_case_partitions(cases: Sequence[Mapping], *, validation_count: int,
                           partition_seed: str = 'hydrogen-v1') -> list[dict]:
    """Preassign an exact blind holdout without inspecting solver outcomes.

    Args:
        cases: Mapping supplying cases.
        validation_count: Number of validation count to use.
        partition_seed: Partition seed used by this operation.

    Returns:
        List of computed or validated records.
    """
    values = [dict(case) for case in cases]
    if (not isinstance(validation_count, int) or validation_count < 1 or
            validation_count >= len(values)):
        raise ValueError('validation_count must leave nonempty train and holdout sets')
    identities = [case.get('case_id') for case in values]
    if any(not isinstance(value, str) or not value for value in identities):
        raise ValueError('every designed case needs a stable case_id')
    if len(set(identities)) != len(identities):
        raise ValueError('designed case identities must be unique')
    ranked = sorted(
        identities,
        key=lambda value: hashlib.sha256(
            f'{partition_seed}:{value}'.encode()).hexdigest())
    holdout = set(ranked[:validation_count])
    for case in values:
        case['partition'] = ('validation' if case['case_id'] in holdout
                             else 'training')
    return values
