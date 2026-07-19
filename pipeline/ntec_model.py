"""Parameterized NTEC assistance model.

This is deliberately a bounded hypothesis model, not a substitute for measured
triboelectric fields or deactivation data.  With no supplied operating evidence
it returns zero assistance, preventing a material from winning merely because it
belongs to a liquid-metal class.
"""

from dataclasses import dataclass, asdict
import math
import json
import os


@dataclass(frozen=True)
class NTECConditions:
    shear_rate_s: float | None = None
    interfacial_field_V_m: float | None = None
    mechanical_power_W_kg: float | None = None
    carbon_detachment_fraction: float | None = None
    field_measurement_source: str | None = None


def conditions_from_environment() -> NTECConditions:
    """Load NTEC_CONDITIONS_JSON. Bad or absent input deliberately means unknown."""
    try:
        raw = json.loads(os.environ.get('NTEC_CONDITIONS_JSON', '{}'))
        allowed = set(NTECConditions.__dataclass_fields__)
        return NTECConditions(**{k: v for k, v in raw.items() if k in allowed})
    except (TypeError, ValueError, json.JSONDecodeError):
        return NTECConditions()


def ntec_assistance(conditions: NTECConditions) -> dict:
    values = asdict(conditions)
    required = ('shear_rate_s', 'interfacial_field_V_m',
                'mechanical_power_W_kg', 'carbon_detachment_fraction')
    missing = [k for k in required if values[k] is None]
    if missing:
        return {'status': 'unknown', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': missing,
                'conditions': values}
    if not conditions.field_measurement_source:
        missing.append('field_measurement_source')
        return {'status': 'unknown', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': missing,
                'conditions': values}
    if any(not math.isfinite(float(values[k])) or float(values[k]) < 0 for k in required):
        return {'status': 'invalid', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': [], 'conditions': values}

    shear = min(float(conditions.shear_rate_s) / 1e4, 1.0)
    field = min(float(conditions.interfacial_field_V_m) / 1e8, 1.0)
    power = min(float(conditions.mechanical_power_W_kg) / 1e3, 1.0)
    detach = min(max(float(conditions.carbon_detachment_fraction), 0.0), 1.0)
    support = min(shear, field, power)
    return {
        'status': 'modeled',
        # Conservative caps must be recalibrated from paired NTEC/control data.
        'barrier_reduction_eV': 0.25 * support,
        'coking_bonus': 3.0 * support * detach,
        'missing': [], 'conditions': values,
        'calibration_required': True,
    }
