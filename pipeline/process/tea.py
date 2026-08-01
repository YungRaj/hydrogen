"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class TEAAssumptions:
    natural_gas_usd_mmbtu: float
    electricity_usd_kwh: float
    base_energy_kwh_kg_h2: float
    capex_usd_kg_h2: float
    carbon_value_usd_kg_h2: float
    source_id: str


SCENARIOS = {
    'optimistic': TEAAssumptions(2.50, 0.035, 7.0, 0.30, 1.00,
                                 'screening-tea-v2:optimistic'),
    'base': TEAAssumptions(3.50, 0.060, 8.5, 0.50, 0.80,
                           'screening-tea-v2:base'),
    'conservative': TEAAssumptions(6.00, 0.120, 12.0, 1.20, 0.20,
                                   'screening-tea-v2:conservative'),
}


def estimate_hydrogen_cost(conversion: float,
                           scenario: str = 'base') -> dict:
    """Return an assumption-labelled screening estimate, never a measurement."""
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    assumptions = SCENARIOS[scenario]
    conversion = float(conversion)
    energy = assumptions.base_energy_kwh_kg_h2 / conversion
    cost = (
        assumptions.natural_gas_usd_mmbtu * 0.05 / conversion +
        energy * assumptions.electricity_usd_kwh +
        assumptions.capex_usd_kg_h2 -
        assumptions.carbon_value_usd_kg_h2
    )
    return {
        'scenario': scenario,
        'h2_cost_usd_kg': round(cost, 3),
        'conversion': conversion,
        'energy_input_kwh_kg_h2': energy,
        'assumptions': asdict(assumptions),
        'evidence_level': 'screening_scenario_not_measured_tea',
    }


def estimate_scenario_range(conversion: float) -> dict:
    estimates = {
        name: estimate_hydrogen_cost(conversion, name)
        for name in SCENARIOS
    }
    return {
        'estimates': estimates,
        'min_usd_kg': min(x['h2_cost_usd_kg'] for x in estimates.values()),
        'max_usd_kg': max(x['h2_cost_usd_kg'] for x in estimates.values()),
        'evidence_level': 'screening_sensitivity_range',
    }
