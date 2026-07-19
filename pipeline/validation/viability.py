"""Fail-closed application-level industrial viability gates.

These are configurable engineering screening defaults, not universal standards.
Missing evidence is UNKNOWN and may not be treated as a pass or used to prune a
branch. Values use fractions rather than percentages unless noted.
"""

from dataclasses import dataclass, asdict
from typing import Mapping
import math


@dataclass(frozen=True)
class TurquoiseHydrogenBounds:
    min_temperature_K: float = 700.0
    max_temperature_K: float = 1300.0
    min_h2_selectivity: float = 0.95
    min_ch4_conversion: float = 0.70
    max_deactivation_fraction_per_h: float = 0.01
    max_coke_fraction: float = 0.05
    max_net_energy_kWh_kg_h2: float = 15.0


@dataclass(frozen=True)
class FuelCellBounds:
    max_orr_overpotential_V: float = 0.40
    min_peak_power_W_cm2: float = 1.00
    min_system_efficiency: float = 0.40
    max_voltage_degradation_uV_h: float = 10.0
    min_measured_hours: float = 100.0


def _number(record: Mapping, names):
    for name in names:
        value = record.get(name)
        if value is not None:
            try:
                value = float(value)
                if math.isfinite(value):
                    return value
            except (TypeError, ValueError):
                pass
    return None


def evaluate_turquoise(record: Mapping,
                       bounds: TurquoiseHydrogenBounds = TurquoiseHydrogenBounds()) -> dict:
    checks = {
        'temperature': (_number(record, ('temperature_K', 'T_K', 'temperature')),
                        bounds.min_temperature_K, bounds.max_temperature_K),
        'h2_selectivity': (_number(record, ('H2_selectivity', 'h2_selectivity')),
                           bounds.min_h2_selectivity, None),
        'ch4_conversion': (_number(record, ('CH4_conversion', 'ch4_conversion')),
                           bounds.min_ch4_conversion, None),
        'deactivation': (_number(record, ('deactivation_fraction_per_h', 'deactivation_rate_per_h')),
                         None, bounds.max_deactivation_fraction_per_h),
        'coke': (_number(record, ('coke_fraction', 'coke_deposition_fraction')),
                 None, bounds.max_coke_fraction),
        'net_energy': (_number(record, ('net_energy_kWh_kg_h2',)),
                       None, bounds.max_net_energy_kWh_kg_h2),
        'experimental_reactor': (_number(record, ('measured_reactor',)), 1.0, None),
    }
    return _evaluate(checks, asdict(bounds))


def evaluate_fuel_cell(record: Mapping,
                       bounds: FuelCellBounds = FuelCellBounds()) -> dict:
    checks = {
        'orr_overpotential': (_number(record, ('orr_overpotential_V', 'orr_overpotential')),
                              None, bounds.max_orr_overpotential_V),
        'peak_power': (_number(record, ('peak_power_W_cm2', 'peak_power_density_W_cm2')),
                       bounds.min_peak_power_W_cm2, None),
        'system_efficiency': (_number(record, ('system_efficiency', 'efficiency')),
                              bounds.min_system_efficiency, None),
        'voltage_degradation': (_number(record, ('voltage_degradation_uV_h',)),
                                None, bounds.max_voltage_degradation_uV_h),
        'measured_hours': (_number(record, ('measured_hours',)), bounds.min_measured_hours, None),
        'experimental_mea': (_number(record, ('measured_mea',)), 1.0, None),
    }
    return _evaluate(checks, asdict(bounds))


def _evaluate(checks, bounds):
    failed, missing, observed = [], [], {}
    for name, (value, minimum, maximum) in checks.items():
        if value is None:
            missing.append(name)
            continue
        observed[name] = value
        if minimum is not None and value < minimum:
            failed.append(f'{name}<{minimum}')
        if maximum is not None and value > maximum:
            failed.append(f'{name}>{maximum}')
    status = 'fail' if failed else ('unknown' if missing else 'pass')
    return {'status': status, 'failed': failed, 'missing': missing,
            'observed': observed, 'bounds': bounds}
