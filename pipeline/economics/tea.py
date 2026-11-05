"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import dataclass

from pipeline.data_models.evidence import (
    MultipathwayChampionScorecard,
    TEAAssumptionsDict,
    TEAScenarioEstimate,
    TEAScenarioRangeResult,
)


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


SCENARIOS: dict[str, TEAAssumptions] = {
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
KG_CH4_PER_KG_H2: float = 16.043 / (2 * 2.016)
MMBTU_CH4_PER_KG_H2: float = KG_CH4_PER_KG_H2 * 55.5 / 1055.06
PSA_CH4_RECOVERY: float = 0.96  # Industrial pressure-swing adsorption (PSA) unreacted CH4 recovery fraction
RECYCLE_ENERGY_FRACTION: float = 0.25  # Fractional compression + reheat energy penalty per unit recycle ratio R = (1-X)/X
RECYCLE_CAPEX_SCALING: float = 0.12  # Fractional PSA + recycle compressor CAPEX increment per unit recycle ratio R


def estimate_hydrogen_cost(
    conversion: float,
    scenario: str = 'base',
    *,
    psa_ch4_recovery: float = PSA_CH4_RECOVERY,
    temperature_K: float | None = None,
    solid_c_selectivity: float = 1.0,
    specific_energy_kwh_kg_h2: float | None = None,
) -> TEAScenarioEstimate:
    """Return an assumption-labelled screening estimate with PSA CH4 recycle, never a measurement.

    Models an industrial PSA H2 purification and unreacted CH4 recycle loop:
    - Recycle ratio ``R = (1 - X) / X`` moles recycled per mole converted.
    - Net methane feedstock utilization ``U_CH4 = X / (1 - r_PSA * (1 - X))``
      (equals 1.0 at ``X = 1.0``, and reflects only PSA tail-gas purge loss at
      partial single-pass conversion rather than venting all unreacted CH4).
    - Recycle compression and sensible-reheat energy penalty scales with ``R``
      as ``base_energy * (1 + RECYCLE_ENERGY_FRACTION * R)``.
    - Recycle loop and PSA vessel sizing scales CAPEX as
      ``capex * (1 + RECYCLE_CAPEX_SCALING * R)``.

    Args:
        conversion: Single-pass methane conversion fraction in (0, 1].
        scenario: Scenario used by this operation.
        psa_ch4_recovery: Fraction of unreacted CH4 recovered and recycled by PSA.
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    assumptions = SCENARIOS[scenario]
    conversion = float(conversion)
    r_psa = min(max(float(psa_ch4_recovery), 0.0), 0.999)
    s_c = min(max(float(solid_c_selectivity), 0.0), 1.0)
    recycle_ratio = (1.0 - conversion) / conversion
    net_ch4_utilization = conversion / max(1.0 - r_psa * (1.0 - conversion), 1e-6)
    base_energy = assumptions.base_energy_kwh_kg_h2
    if specific_energy_kwh_kg_h2 is not None and float(specific_energy_kwh_kg_h2) > 0.0:
        base_energy = float(specific_energy_kwh_kg_h2)
    elif temperature_K is not None and float(temperature_K) > 298.15:
        # Scale baseline thermal energy relative to a 1073.15 K (800 °C) reference reactor:
        # 65% endothermic reaction enthalpy ΔH_rxn(T) + 35% sensible feed heating (T - 298.15 K).
        try:
            from pipeline.reactors.carbon_accounting import (
                pyrolysis_reaction_enthalpy_kJ_mol,
            )

            dh_t = pyrolysis_reaction_enthalpy_kJ_mol(float(temperature_K))
            dh_ref = pyrolysis_reaction_enthalpy_kJ_mol(1073.15)
        except ImportError:
            # Linear heat-capacity approximation ΔCp ≈ 0.0155 kJ/(mol·K) above ΔH°_298 = 74.85 kJ/mol
            dh_t = 74.85 + 0.0155 * (float(temperature_K) - 298.15)
            dh_ref = 74.85 + 0.0155 * (1073.15 - 298.15)
        sensible_ratio = max(float(temperature_K) - 298.15, 0.0) / (1073.15 - 298.15)
        thermal_scale = 0.65 * (dh_t / max(dh_ref, 1.0)) + 0.35 * sensible_ratio
        base_energy = base_energy * max(thermal_scale, 0.5)
    energy = base_energy * (
        1.0 + RECYCLE_ENERGY_FRACTION * recycle_ratio
    )
    effective_capex = assumptions.capex_usd_kg_h2 * (
        1.0 + RECYCLE_CAPEX_SCALING * recycle_ratio
    )
    carbon_credit = assumptions.carbon_value_usd_kg_h2 * s_c
    cost = (
        assumptions.natural_gas_usd_mmbtu * MMBTU_CH4_PER_KG_H2 / net_ch4_utilization
        + energy * assumptions.electricity_usd_kwh
        + effective_capex
        - carbon_credit
    )
    assumptions_dict: TEAAssumptionsDict = {
        'natural_gas_usd_mmbtu': assumptions.natural_gas_usd_mmbtu,
        'electricity_usd_kwh': assumptions.electricity_usd_kwh,
        'base_energy_kwh_kg_h2': assumptions.base_energy_kwh_kg_h2,
        'capex_usd_kg_h2': assumptions.capex_usd_kg_h2,
        'carbon_value_usd_kg_h2': assumptions.carbon_value_usd_kg_h2,
        'source_id': assumptions.source_id,
    }
    return {
        'scenario': scenario,
        'h2_cost_usd_kg': round(cost, 3),
        'conversion': conversion,
        'recycle_ratio': round(recycle_ratio, 4),
        'net_ch4_utilization': round(net_ch4_utilization, 4),
        'psa_ch4_recovery': r_psa,
        'solid_c_selectivity': round(s_c, 4),
        'carbon_credit_usd_kg_h2': round(carbon_credit, 4),
        'energy_input_kwh_kg_h2': energy,
        'effective_capex_usd_kg_h2': round(effective_capex, 4),
        'assumptions': assumptions_dict,
        'evidence_level': 'screening_scenario_not_measured_tea',
    }


def estimate_scenario_range(
    conversion: float,
    *,
    temperature_K: float | None = None,
    solid_c_selectivity: float = 1.0,
    specific_energy_kwh_kg_h2: float | None = None,
) -> TEAScenarioRangeResult:
    """Evaluate hydrogen cost across named techno-economic scenarios.

    Args:
        conversion: Methane conversion fraction used by the cost model.
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.

    Returns:
        A dictionary containing estimate scenario range outputs, status, and supporting metadata.
    """
    estimates: dict[str, TEAScenarioEstimate] = {
        name: estimate_hydrogen_cost(
            conversion,
            name,
            temperature_K=temperature_K,
            solid_c_selectivity=solid_c_selectivity,
            specific_energy_kwh_kg_h2=specific_energy_kwh_kg_h2,
        )
        for name in SCENARIOS
    }
    return {
        'estimates': estimates,
        'min_usd_kg': min(x['h2_cost_usd_kg'] for x in estimates.values()),
        'max_usd_kg': max(x['h2_cost_usd_kg'] for x in estimates.values()),
        'evidence_level': 'screening_sensitivity_range',
    }


# Lower heating value of H2 (33.33 kWh/kg H2) for coupled H2 -> PEMFC LCOE.
LHV_H2_KWH_PER_KG: float = 33.33


def build_multipathway_champion_scorecard(
    results_root: str = 'results',
    *,
    stack_cells: int = 400,
    top_k_catalysts: int = 11,
) -> MultipathwayChampionScorecard:
    """Delegate multi-pathway scorecard aggregation to pipeline.evidence.multipathway_scorecard.

    Args:
        results_root: Root directory containing stage outputs.
        stack_cells: Number of series cells in the automotive PEMFC stack.
        top_k_catalysts: Number of top catalysts to evaluate across pathways.

    Returns:
        Multi-pathway scorecard dictionary written to disk.
    """
    from pipeline.evidence.multipathway_scorecard import (
        build_multipathway_champion_scorecard as _build,
    )

    return _build(
        results_root,
        stack_cells=stack_cells,
        top_k_catalysts=top_k_catalysts,
    )

