"""Parameterized NTEC assistance model.

This is deliberately a bounded transfer model, not a substitute for measured
triboelectric fields or deactivation data. It returns zero assistance unless the
operating conditions and the effect measured against a paired thermocatalytic
control are supplied. This prevents a material from winning merely because it
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
    paired_control_source: str | None = None
    paired_control_count: int | None = None
    measured_barrier_reduction_eV: float | None = None
    measured_coking_delta_eV: float | None = None


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
    operating = ('shear_rate_s', 'interfacial_field_V_m',
                 'mechanical_power_W_kg', 'carbon_detachment_fraction')
    calibration = ('paired_control_count', 'measured_barrier_reduction_eV',
                   'measured_coking_delta_eV')
    required = operating + calibration
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
    if not conditions.paired_control_source:
        missing.append('paired_control_source')
        return {'status': 'unknown', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': missing,
                'conditions': values}
    if any(not math.isfinite(float(values[k])) or float(values[k]) < 0 for k in required):
        return {'status': 'invalid', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': [], 'conditions': values}
    if int(conditions.paired_control_count) < 1:
        return {'status': 'invalid', 'barrier_reduction_eV': 0.0,
                'coking_bonus': 0.0, 'missing': [], 'conditions': values}

    shear = min(float(conditions.shear_rate_s) / 1e4, 1.0)
    field = min(float(conditions.interfacial_field_V_m) / 1e8, 1.0)
    power = min(float(conditions.mechanical_power_W_kg) / 1e3, 1.0)
    detach = min(max(float(conditions.carbon_detachment_fraction), 0.0), 1.0)
    support = min(shear, field, power)
    return {
        'status': 'paired_control_calibrated',
        # Transfer is bounded by both measured effect and operating support.
        # It remains modeled evidence for a new catalyst, not validation of it.
        'barrier_reduction_eV': min(
            float(conditions.measured_barrier_reduction_eV), 0.25) * support,
        'coking_bonus': min(
            float(conditions.measured_coking_delta_eV), 3.0) * support * detach,
        'missing': [], 'conditions': values,
        'calibration_required': False,
        'evidence_level': 'paired_control_transfer_model',
    }
