"""Fail-closed application-level industrial viability gates.

These are configurable engineering screening defaults, not universal standards.
Missing evidence is UNKNOWN and may not be treated as a pass or used to prune a
branch. Values use fractions rather than percentages unless noted.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import asdict, dataclass
import math
from typing import Any


@dataclass(frozen=True)
class TurquoiseHydrogenBounds:
    """Define measured viability limits for turquoise-hydrogen operation.

    Attributes:
        min_temperature_K: Configured min temperature K value.
        max_temperature_K: Configured max temperature K value.
        min_h2_selectivity: Configured min h2 selectivity value.
        min_ch4_conversion: Configured min ch4 conversion value.
        max_deactivation_fraction_per_h: Configured max deactivation fraction per h value.
        max_coke_fraction: Configured max coke fraction value.
        max_net_energy_kWh_kg_h2: Configured max net energy kWh kg h2 value.
    """

    min_temperature_K: float = 700.0  # Lower kinetic floor for methane pyrolysis [K]
    max_temperature_K: float = 1300.0  # Upper reactor alloy/refractory wall ceiling [K]
    # Legacy field name: gated quantity is an H-atom balance / reported metric,
    # not true branching H2 selectivity from competing C2 pathways.
    min_h2_selectivity: float = 0.95  # Minimum H-atom balance to gas-phase H2 [-]
    enforce_h2_selectivity_gate: bool = False
    min_ch4_conversion: float = 0.70  # Target industrial single-pass CH4 conversion [-]
    max_deactivation_fraction_per_h: float = 0.01  # Maximum acceptable activity loss rate [h^-1]
    max_coke_fraction: float = 0.05  # Maximum non-harvestable fouling carbon fraction [-]
    max_net_energy_kWh_kg_h2: float = 15.0  # Maximum net thermal/electrical input [kWh/kg H2]
    reject_co2_permitted_runs: bool = True


@dataclass(frozen=True)
class FuelCellBounds:
    """Define measured viability limits for fuel-cell operation.

    Attributes:
        max_orr_overpotential_V: Configured max orr overpotential V value.
        min_peak_power_W_cm2: Configured min peak power W cm2 value.
        min_system_efficiency: Configured min system efficiency value.
        max_voltage_degradation_uV_h: Configured max voltage degradation uV h value.
        min_measured_hours: Configured min measured hours value.
    """

    max_orr_overpotential_V: float = 0.40  # Maximum cathode ORR kinetic overpotential [V]
    min_peak_power_W_cm2: float = 1.00  # Minimum MEA peak power density [W/cm^2] (US DOE target)
    min_system_efficiency: float = 0.40  # Minimum net stack/system LHV efficiency [-]
    max_voltage_degradation_uV_h: float = 10.0  # Maximum voltage decay rate [uV/h]
    min_measured_hours: float = 100.0  # Minimum accelerated durability test duration [h]


def _number(record: Mapping[str, Any], names: Sequence[str]) -> float | None:
    for name in names:
        value = record.get(name)
        if value is not None:
            try:
                num = float(value)
                if math.isfinite(num):
                    return num
            except (TypeError, ValueError):
                pass
    return None


def evaluate_turquoise(
    record: Mapping[str, Any],
    bounds: TurquoiseHydrogenBounds = TurquoiseHydrogenBounds(),
) -> dict[str, Any]:
    """Evaluate measured turquoise-hydrogen performance against viability bounds.

    Args:
        record: Measured performance record to evaluate.
        bounds: Application-specific viability limits.

    Returns:
        A dictionary containing evaluate turquoise outputs, status, and supporting metadata.
    """
    checks: dict[str, tuple[float | None, float | None, float | None]] = {
        'temperature': (
            _number(record, ('temperature_K', 'T_K', 'temperature')),
            bounds.min_temperature_K,
            bounds.max_temperature_K,
        ),
        'ch4_conversion': (
            _number(
                record,
                (
                    'CH4_conversion',
                    'ch4_conversion',
                    'single_pass_CH4_conversion',
                    'cycle_averaged_CH4_conversion',
                ),
            ),
            bounds.min_ch4_conversion,
            None,
        ),
        'deactivation': (
            _number(record, ('deactivation_fraction_per_h', 'deactivation_rate_per_h')),
            None,
            bounds.max_deactivation_fraction_per_h,
        ),
        'coke': (
            _number(record, ('coke_fraction', 'coke_deposition_fraction', 'coke_selectivity')),
            None,
            bounds.max_coke_fraction,
        ),
        'net_energy': (
            _number(record, ('net_energy_kWh_kg_h2', 'specific_energy_kWh_kg_H2')),
            None,
            bounds.max_net_energy_kWh_kg_h2,
        ),
        'experimental_reactor': (_number(record, ('measured_reactor',)), 1.0, None),
    }
    if bounds.enforce_h2_selectivity_gate:
        checks['h2_selectivity'] = (
            _number(record, ('H2_selectivity', 'h2_selectivity', 'H2_atom_balance')),
            bounds.min_h2_selectivity,
            None,
        )
    result = _evaluate(checks, asdict(bounds))
    if bounds.reject_co2_permitted_runs and record.get('co2_permitted') is True:
        result['failed'] = list(result.get('failed') or []) + ['co2_permitted=True']
        result['status'] = 'fail'
    return result


def evaluate_fuel_cell(
    record: Mapping[str, Any],
    bounds: FuelCellBounds = FuelCellBounds(),
) -> dict[str, Any]:
    """Evaluate measured fuel-cell performance against viability bounds.

    Args:
        record: Measured performance record to evaluate.
        bounds: Application-specific viability limits.

    Returns:
        A dictionary containing evaluate fuel cell outputs, status, and supporting metadata.
    """
    checks: dict[str, tuple[float | None, float | None, float | None]] = {
        'orr_overpotential': (
            _number(record, ('orr_overpotential_V', 'orr_overpotential')),
            None,
            bounds.max_orr_overpotential_V,
        ),
        'peak_power': (
            _number(record, ('peak_power_W_cm2', 'peak_power_density_W_cm2')),
            bounds.min_peak_power_W_cm2,
            None,
        ),
        'system_efficiency': (
            _number(record, ('system_efficiency', 'efficiency')),
            bounds.min_system_efficiency,
            None,
        ),
        'voltage_degradation': (
            _number(record, ('voltage_degradation_uV_h',)),
            None,
            bounds.max_voltage_degradation_uV_h,
        ),
        'measured_hours': (
            _number(record, ('measured_hours', 'durability_hours')),
            bounds.min_measured_hours,
            None,
        ),
        'experimental_mea': (_number(record, ('measured_mea',)), 1.0, None),
    }
    return _evaluate(checks, asdict(bounds))


def _evaluate(
    checks: Mapping[str, tuple[float | None, float | None, float | None]],
    bounds: dict[str, Any],
) -> dict[str, Any]:
    failed: list[str] = []
    missing: list[str] = []
    observed: dict[str, float] = {}
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
    return {
        'status': status,
        'failed': failed,
        'missing': missing,
        'observed': observed,
        'bounds': bounds,
    }
