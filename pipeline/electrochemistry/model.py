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
from typing import cast

# IAPWS-95 critical temperature of pure water (373.946 C, P_c = 22.064 MPa),
# above which no distinct liquid aqueous electrolyte phase exists.
WATER_CRITICAL_TEMPERATURE_K = 647.096


@dataclass(frozen=True)
class ElectrochemicalConditions:
    """Hold measured electrochemical operating conditions and provenance.

    Attributes:
        electrolyte_phase: Configured electrolyte phase value.
        electrolyte_identity: Configured electrolyte identity value.
        applied_potential_V: Configured applied potential V value.
        current_density_A_cm2: Configured current density A cm2 value.
        faradaic_efficiency_H2: Configured faradaic efficiency H2 value.
        methane_conversion: Configured methane conversion value.
        temperature_K: Configured temperature K value.
        pressure_Pa: Configured pressure Pa value.
        measurement_source: Configured measurement source value.
        paired_control_source: Configured paired control source value.
    """

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
    """Load pathway conditions from the corresponding JSON environment variable.

    Returns:
        A populated conditions object; malformed input yields empty conditions.
    """
    try:
        raw = json.loads(os.environ.get('ELECTROCHEMICAL_CONDITIONS_JSON', '{}'))
        if not isinstance(raw, dict):
            return ElectrochemicalConditions()
        raw_dict = cast(dict[str, object], raw)
        def _opt_str(key: str) -> str | None:
            v = raw_dict.get(key)
            return str(v) if isinstance(v, str) else None

        def _opt_float(key: str) -> float | None:
            v = raw_dict.get(key)
            return float(v) if isinstance(v, (int, float)) else None

        return ElectrochemicalConditions(
            electrolyte_phase=_opt_str('electrolyte_phase'),
            electrolyte_identity=_opt_str('electrolyte_identity'),
            applied_potential_V=_opt_float('applied_potential_V'),
            current_density_A_cm2=_opt_float('current_density_A_cm2'),
            faradaic_efficiency_H2=_opt_float('faradaic_efficiency_H2'),
            methane_conversion=_opt_float('methane_conversion'),
            temperature_K=_opt_float('temperature_K'),
            pressure_Pa=_opt_float('pressure_Pa'),
            measurement_source=_opt_str('measurement_source'),
            paired_control_source=_opt_str('paired_control_source'),
        )
    except (TypeError, ValueError, json.JSONDecodeError):
        return ElectrochemicalConditions()


def electrochemical_evidence(conditions: ElectrochemicalConditions) -> dict[str, object]:
    """Validate electrochemical conditions and classify their evidence status.

    Args:
        conditions: Measured pathway conditions to validate.

    Returns:
        A status record containing validated conditions and missing-evidence details.
    """
    values = asdict(conditions)
    required = (
        'electrolyte_phase',
        'electrolyte_identity',
        'applied_potential_V',
        'current_density_A_cm2',
        'faradaic_efficiency_H2',
        'methane_conversion',
        'temperature_K',
        'pressure_Pa',
        'measurement_source',
        'paired_control_source',
    )
    missing = [key for key in required if values[key] is None]
    phase = str(conditions.electrolyte_phase or '').lower()
    if phase and phase not in {'aqueous', 'molten'}:
        return {
            'status': 'invalid',
            'missing': [],
            'conditions': values,
            'reason': 'electrolyte_phase must be aqueous or molten',
        }
    numeric = (
        'applied_potential_V',
        'current_density_A_cm2',
        'faradaic_efficiency_H2',
        'methane_conversion',
        'temperature_K',
        'pressure_Pa',
    )
    if any(
        values[key] is not None and not math.isfinite(float(values[key]))
        for key in numeric
    ):
        return {
            'status': 'invalid',
            'missing': [],
            'conditions': values,
            'reason': 'all numerical conditions must be finite',
        }
    if missing:
        return {'status': 'unknown', 'missing': missing, 'conditions': values}
    fe_h2 = float(values['faradaic_efficiency_H2'])
    ch4_conv = float(values['methane_conversion'])
    j_a_cm2 = float(values['current_density_A_cm2'])
    temp_k = float(values['temperature_K'])
    press_pa = float(values['pressure_Pa'])
    # Faradaic efficiency FE_H2 and single-pass CH4 conversion X_CH4 are bounded fractions in [0, 1]
    if not 0 <= fe_h2 <= 1 or not 0 <= ch4_conv <= 1:
        return {
            'status': 'invalid',
            'missing': [],
            'conditions': values,
            'reason': 'efficiency and conversion must lie in [0, 1]',
        }
    if j_a_cm2 < 0 or temp_k <= 0 or press_pa <= 0:
        return {
            'status': 'invalid',
            'missing': [],
            'conditions': values,
            'reason': 'current density must be nonnegative and T/P positive',
        }
    if phase == 'aqueous' and temp_k >= WATER_CRITICAL_TEMPERATURE_K:
        return {
            'status': 'invalid',
            'missing': [],
            'conditions': values,
            'reason': 'aqueous electrolyte temperature_K must be below 647.096 K',
        }
    return {
        'status': 'measured_paired_control',
        'missing': [],
        'conditions': values,
        'evidence_level': 'measured_operating_point',
    }

