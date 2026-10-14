"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import asdict, dataclass


@dataclass(frozen=True)
class TEAAssumptions:
    """Hold a traceable techno-economic scenario.

    Attributes:
        natural_gas_usd_mmbtu: Configured natural gas usd mmbtu value.
        electricity_usd_kwh: Configured electricity usd kwh value.
        base_energy_kwh_kg_h2: Configured base energy kwh kg h2 value.
        capex_usd_kg_h2: Configured capex usd kg h2 value.
        carbon_value_usd_kg_h2: Configured carbon value usd kg h2 value.
        source_id: Configured source id value.
    """

    natural_gas_usd_mmbtu: float
    electricity_usd_kwh: float
    base_energy_kwh_kg_h2: float
    capex_usd_kg_h2: float
    carbon_value_usd_kg_h2: float
    source_id: str


SCENARIOS = {
    'optimistic': TEAAssumptions(
        2.50, 0.035, 7.0, 0.30, 1.00, 'screening-tea-v3:optimistic'
    ),
    'base': TEAAssumptions(3.50, 0.060, 8.5, 0.50, 0.80, 'screening-tea-v3:base'),
    'conservative': TEAAssumptions(
        6.00, 0.120, 12.0, 1.20, 0.20, 'screening-tea-v3:conservative'
    ),
}


# CH4 -> C + 2 H2: 16.04 g CH4 per 4.032 g H2, i.e. 3.98 kg CH4 per kg H2.
# Priced on the HHV basis used for natural-gas $/MMBtu quotes
# (55.5 MJ/kg CH4; 1 MMBtu = 1055.06 MJ): ~0.209 MMBtu per kg H2 at X = 1.
# (v2 used 0.05 MMBtu/kg, about 4x too little feedstock.)
KG_CH4_PER_KG_H2 = 16.043 / (2 * 2.016)
MMBTU_CH4_PER_KG_H2 = KG_CH4_PER_KG_H2 * 55.5 / 1055.06

def estimate_hydrogen_cost(conversion: float, scenario: str = 'base') -> dict:
    """Return an assumption-labelled screening estimate, never a measurement.

    Args:
        conversion: Conversion used by this operation.
        scenario: Scenario used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    assumptions = SCENARIOS[scenario]
    conversion = float(conversion)
    energy = assumptions.base_energy_kwh_kg_h2 / conversion
    cost = (
        assumptions.natural_gas_usd_mmbtu * MMBTU_CH4_PER_KG_H2 / conversion
        + energy * assumptions.electricity_usd_kwh
        + assumptions.capex_usd_kg_h2
        - assumptions.carbon_value_usd_kg_h2
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
    """Evaluate hydrogen cost across named techno-economic scenarios.

    Args:
        conversion: Methane conversion fraction used by the cost model.

    Returns:
        A dictionary containing estimate scenario range outputs, status, and supporting metadata.
    """
    estimates = {name: estimate_hydrogen_cost(conversion, name) for name in SCENARIOS}
    return {
        'estimates': estimates,
        'min_usd_kg': min(x['h2_cost_usd_kg'] for x in estimates.values()),
        'max_usd_kg': max(x['h2_cost_usd_kg'] for x in estimates.values()),
        'evidence_level': 'screening_sensitivity_range',
    }
