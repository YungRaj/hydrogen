"""Fail-closed configuration contract for electrochemical methane conversion.

Aqueous and molten electrolytes share one top-level pathway.  Their different
transport, thermodynamic, and kinetic models are selected by ``electrolyte_phase``.
Until measured or validated pathway inputs exist, this module reports missing
evidence and never fabricates conversion, selectivity, or electrical efficiency.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
import math
import os


@dataclass(frozen=True)
class ElectrochemicalConditions:
    electrolyte_phase: str | None = None
    electrolyte_identity: str | None = None
    applied_potential_V: float | None = None
    current_density_A_cm2: float | None = None
    faradaic_efficiency_H2: float | None = None
    methane_conversion: float | None = None
    temperature_K: float | None = None
    pressure_Pa: float | None = None
    measurement_source: str | None = None
    paired_control_source: str | None = None


def conditions_from_environment() -> ElectrochemicalConditions:
    try:
        raw = json.loads(os.environ.get('ELECTROCHEMICAL_CONDITIONS_JSON', '{}'))
        allowed = set(ElectrochemicalConditions.__dataclass_fields__)
        return ElectrochemicalConditions(
            **{key: value for key, value in raw.items() if key in allowed})
    except (TypeError, ValueError, json.JSONDecodeError):
        return ElectrochemicalConditions()


def electrochemical_evidence(conditions: ElectrochemicalConditions) -> dict:
    values = asdict(conditions)
    required = (
        'electrolyte_phase', 'electrolyte_identity', 'applied_potential_V',
        'current_density_A_cm2', 'faradaic_efficiency_H2',
        'methane_conversion', 'temperature_K', 'pressure_Pa',
        'measurement_source', 'paired_control_source',
    )
    missing = [key for key in required if values[key] is None]
    phase = str(conditions.electrolyte_phase or '').lower()
    if phase and phase not in {'aqueous', 'molten'}:
        return {'status': 'invalid', 'missing': [], 'conditions': values,
                'reason': 'electrolyte_phase must be aqueous or molten'}
    numeric = (
        'applied_potential_V', 'current_density_A_cm2',
        'faradaic_efficiency_H2', 'methane_conversion', 'temperature_K',
        'pressure_Pa',
    )
    if any(values[key] is not None and not math.isfinite(float(values[key]))
           for key in numeric):
        return {'status': 'invalid', 'missing': [], 'conditions': values,
                'reason': 'all numerical conditions must be finite'}
    if missing:
        return {'status': 'unknown', 'missing': missing, 'conditions': values}
    if not 0 <= float(conditions.faradaic_efficiency_H2) <= 1 or not 0 <= float(
            conditions.methane_conversion) <= 1:
        return {'status': 'invalid', 'missing': [], 'conditions': values,
                'reason': 'efficiency and conversion must lie in [0, 1]'}
    if float(conditions.current_density_A_cm2) < 0 or \
            float(conditions.temperature_K) <= 0 or float(conditions.pressure_Pa) <= 0:
        return {'status': 'invalid', 'missing': [], 'conditions': values,
                'reason': 'current density must be nonnegative and T/P positive'}
    return {
        'status': 'measured_paired_control', 'missing': [],
        'conditions': values, 'evidence_level': 'measured_operating_point',
    }

