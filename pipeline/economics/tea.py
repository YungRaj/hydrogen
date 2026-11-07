"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import dataclass

from pipeline.data_models.evidence import (
    MultipathwayChampionScorecard,
    TEAAssumptionsDict,
    TEAScenarioEstimate,
    TEAScenarioRangeResult,
)
from pipeline.data_models.reactors import FeedstockName


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


# Default methane (CH4) stoichiometric reference constants:
# CH4 -> C + 2 H2: 16.04 g CH4 per 4.032 g H2, i.e. 3.98 kg CH4 per kg H2.
# Priced on the HHV basis used for natural-gas $/MMBtu quotes
# (55.5 MJ/kg CH4; 1 MMBtu = 1055.06 MJ): ~0.209 MMBtu per kg H2 at X = 1.
# When a non-methane feedstock is selected in `estimate_hydrogen_cost(feedstock=...)`,
# the feedstock's exact `spec.kg_feed_per_kg_h2` and `spec.hhv_MJ_kg` are used instead.
KG_CH4_PER_KG_H2: float = 16.043 / (2 * 2.016)
MMBTU_CH4_PER_KG_H2: float = KG_CH4_PER_KG_H2 * 55.5 / 1055.06
PSA_FEEDSTOCK_RECOVERY: float = 0.96  # Industrial pressure-swing adsorption (PSA) unreacted hydrocarbon recovery fraction
RECYCLE_ENERGY_FRACTION: float = 0.25  # Fractional compression + reheat energy penalty per unit recycle ratio R = (1-X)/X
RECYCLE_CAPEX_SCALING: float = 0.12  # Fractional PSA + recycle compressor CAPEX increment per unit recycle ratio R


def estimate_hydrogen_cost(
    conversion: float,
    scenario: str = 'base',
    *,
    psa_feedstock_recovery: float = PSA_FEEDSTOCK_RECOVERY,
    temperature_K: float | None = None,
    solid_c_selectivity: float = 1.0,
    specific_energy_kwh_kg_h2: float | None = None,
    feedstock: FeedstockName | str = 'CH4',
) -> TEAScenarioEstimate:
    """Return an assumption-labelled screening estimate with PSA hydrocarbon recycle, never a measurement.

    Models an industrial PSA H2 purification and unreacted feedstock recycle loop:
    - Recycle ratio ``R = (1 - X) / X`` moles recycled per mole converted.
    - Net feedstock utilization ``U_feed = X / (1 - r_PSA * (1 - X))``
      (equals 1.0 at ``X = 1.0``, and reflects only PSA tail-gas purge loss at
      partial single-pass conversion rather than venting all unreacted feed).
    - Recycle compression and sensible-reheat energy penalty scales with ``R``
      as ``base_energy * (1 + RECYCLE_ENERGY_FRACTION * R)``.
    - Recycle loop and PSA vessel sizing scales CAPEX as
      ``capex * (1 + RECYCLE_CAPEX_SCALING * R)``.
    - Defaults strictly to ``feedstock='CH4'``. When another hydrocarbon or biomass
      feedstock is specified, accounts for its exact stoichiometric mass ratio
      (``kg_feed_per_kg_h2``), endothermic reaction enthalpy per mole of H2,
      and solid-carbon co-product yield (``kg_c_per_kg_h2``).

    Args:
        conversion: Single-pass hydrocarbon/biomass conversion fraction in (0, 1].
        scenario: Scenario used by this operation.
        psa_feedstock_recovery: PSA unreacted feedstock recovery fraction in [0, 1).
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.
        feedstock: Hydrocarbon or biomass feedstock formula or alias (strict default ``'CH4'``).

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    from pipeline.reactors.feedstocks import (
        DEFAULT_FEEDSTOCK,
        HYDROCARBON_FEEDSTOCKS,
        HydrocarbonFeedstockSpec,
        feedstock_reaction_enthalpy_kJ_mol,
        resolve_feedstock,
    )

    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    spec: HydrocarbonFeedstockSpec = resolve_feedstock(feedstock)
    ch4_ref_spec: HydrocarbonFeedstockSpec = HYDROCARBON_FEEDSTOCKS[DEFAULT_FEEDSTOCK]
    assumptions: TEAAssumptions = SCENARIOS[scenario]
    conv_val: float = float(conversion)
    r_psa: float = min(max(float(psa_feedstock_recovery), 0.0), 0.999)
    s_c: float = min(max(float(solid_c_selectivity), 0.0), 1.0)
    recycle_ratio: float = (1.0 - conv_val) / conv_val
    net_feedstock_utilization: float = conv_val / max(
        1.0 - r_psa * (1.0 - conv_val), 1e-6
    )
    base_energy: float = assumptions.base_energy_kwh_kg_h2
    dh_per_h2_ratio: float
    if specific_energy_kwh_kg_h2 is not None and float(specific_energy_kwh_kg_h2) > 0.0:
        base_energy = float(specific_energy_kwh_kg_h2)
    elif temperature_K is not None and float(temperature_K) > 298.15:
        # Scale baseline thermal energy relative to a 1073.15 K (800 °C) CH4 reference reactor:
        # 65% endothermic reaction enthalpy per mole of H2 + 35% sensible feed heating (T - 298.15 K).
        dh_ref_ch4: float
        dh_t_feed: float
        try:
            from pipeline.reactors.carbon_accounting import (
                pyrolysis_reaction_enthalpy_kJ_mol,
            )

            dh_ref_ch4 = pyrolysis_reaction_enthalpy_kJ_mol(1073.15)
            dh_t_feed = (
                pyrolysis_reaction_enthalpy_kJ_mol(float(temperature_K))
                if spec.formula == DEFAULT_FEEDSTOCK
                else feedstock_reaction_enthalpy_kJ_mol(float(temperature_K), spec)
            )
        except ImportError:
            dh_ref_ch4 = 74.85 + 0.0155 * (1073.15 - 298.15)
            dh_t_feed = feedstock_reaction_enthalpy_kJ_mol(float(temperature_K), spec)
        dh_per_h2_ratio = (dh_t_feed / spec.nu_H2) / max(
            dh_ref_ch4 / ch4_ref_spec.nu_H2, 1.0
        )
        sensible_ratio: float = max(float(temperature_K) - 298.15, 0.0) / (
            1073.15 - 298.15
        )
        thermal_scale: float = 0.65 * dh_per_h2_ratio + 0.35 * sensible_ratio
        base_energy = base_energy * max(thermal_scale, 0.5)
    elif spec.formula != DEFAULT_FEEDSTOCK:
        # Even without explicit temperature_K, scale reaction-enthalpy share by ΔH°_298 per mol H2
        dh_per_h2_ratio = (
            spec.dh_rxn_298_kJ_mol_h2 / ch4_ref_spec.dh_rxn_298_kJ_mol_h2
        )
        base_energy = base_energy * max(0.65 * dh_per_h2_ratio + 0.35, 0.5)

    energy: float = base_energy * (
        1.0 + RECYCLE_ENERGY_FRACTION * recycle_ratio
    )
    effective_capex: float = assumptions.capex_usd_kg_h2 * (
        1.0 + RECYCLE_CAPEX_SCALING * recycle_ratio
    )
    mmbtu_feed_per_kg_h2: float
    price_mult: float
    c_mass_ratio: float
    if spec.formula == DEFAULT_FEEDSTOCK:
        mmbtu_feed_per_kg_h2 = MMBTU_CH4_PER_KG_H2
        price_mult = 1.0
        c_mass_ratio = 1.0
    else:
        mmbtu_feed_per_kg_h2 = spec.kg_feed_per_kg_h2 * spec.hhv_MJ_kg / 1055.06
        price_mult = spec.base_price_usd_mmbtu_multiplier
        c_mass_ratio = spec.kg_c_per_kg_h2 / ch4_ref_spec.kg_c_per_kg_h2

    carbon_credit: float = assumptions.carbon_value_usd_kg_h2 * c_mass_ratio * s_c
    cost: float = (
        (assumptions.natural_gas_usd_mmbtu * price_mult)
        * mmbtu_feed_per_kg_h2
        / net_feedstock_utilization
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
        'feedstock': spec.formula,
        'feedstock_name': spec.name,
        'kg_feed_per_kg_h2': round(spec.kg_feed_per_kg_h2, 4),
        'kg_c_per_kg_h2': round(spec.kg_c_per_kg_h2, 4),
        'h2_cost_usd_kg': round(cost, 3),
        'conversion': conv_val,
        'recycle_ratio': round(recycle_ratio, 4),
        'net_feedstock_utilization': round(net_feedstock_utilization, 4),
        'psa_feedstock_recovery': r_psa,
        'net_ch4_utilization': round(net_feedstock_utilization, 4),
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
    feedstock: FeedstockName | str = 'CH4',
) -> TEAScenarioRangeResult:
    """Evaluate hydrogen cost across named techno-economic scenarios.

    Args:
        conversion: Hydrocarbon conversion fraction used by the cost model.
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.
        feedstock: Hydrocarbon feedstock formula or alias (strict default ``'CH4'``).

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
            feedstock=feedstock,
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

