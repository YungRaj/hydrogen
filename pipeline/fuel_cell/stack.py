#!/usr/bin/env python3
# Fuel-cell stack and techno-economic process model.
"""
Fuel Cell Stack-Level Scaling Model.

Scales single-cell PEMFC results to a practical fuel cell stack:
  - Single-cell → N-cell stack voltage & power
  - Thermal management (heat rejection load)
  - System efficiency including BOP (balance of plant)
  - Gravimetric/volumetric power density (W/kg, W/L)
  - Techno-economic analysis ($/kW)
"""

from collections.abc import Mapping
from dataclasses import dataclass
import math
from typing import TypedDict

from pipeline.utils import (
    setup_logger,
    save_json,
    F_const,
)
from pipeline.data_models.fuel_cells import StackResult
from pipeline.fuel_cell.pemfc import (
    H2_LHV_J_MOL,
    LHV_VOLTAGE_V,
    PEMFCConfig,
    simulate_pemfc,
)

logger = setup_logger('fc_stack', 'fuel_cell/stack_model.log')


class CoupledReactorPEMFCResult(TypedDict):
    """End-to-end coupled reactor + PSA purification + PEMFC stack well-to-power result."""

    feedstock: str
    reactor_type: str
    catalyst_name: str
    single_pass_conversion: float
    h2_selectivity: float
    h2_yield: float
    psa_h2_recovery: float
    post_psa_h2_purity: float
    post_psa_co_ppm: float
    post_psa_c2h4_ppm: float
    anode_hor_activity_retention: float
    cathode_catalyst: str
    cell_peak_power_W_cm2: float
    cell_rated_power_W_cm2: float | None
    stack_cells: int
    stack_net_power_kW: float
    stack_system_efficiency: float
    well_to_hydrogen_lhv_efficiency: float
    well_to_power_lhv_efficiency: float
    lcoh_base_usd_kg: float
    lcoh_optimistic_usd_kg: float
    lcoe_fuel_usd_kwh_e: float
    lcoe_stack_capex_usd_kwh_e: float
    lcoe_total_usd_kwh_e: float
    lcoe_optimistic_total_usd_kwh_e: float


@dataclass
class StackConfig:
    """Fuel cell stack configuration."""

    # Stack geometry
    n_cells: int = 300  # Number of cells in stack
    active_area_cm2: float = 300.0  # Active area per cell

    # Operating point (from single-cell model)
    cell_voltage_V: float = 0.65  # Operating cell voltage
    current_density_A_cm2: float = 1.5  # Operating current density

    # Component weights (kg per cell)
    mea_mass_g: float = 5.0  # MEA mass per cell
    bipolar_plate_mass_g: float = 150.0  # Bipolar plate per cell
    endplate_mass_g: float = 3000.0  # Total endplate mass (both ends)
    gasket_mass_g: float = 3.0  # Gasket per cell

    # Component volumes (cm³ per cell)
    cell_pitch_mm: float = 2.5  # Cell-to-cell pitch (thickness)

    # Catalyst cost
    pgm_loading_mg_cm2: float = 0.1  # Total PGM loading (both sides)
    pgm_price_usd_g: float = 31.0  # Pt price per gram

    # Membrane cost
    membrane_cost_usd_cm2: float = 0.025

    # BOP parameters
    compressor_efficiency: float = 0.75
    compressor_pressure_ratio: float = 1.5
    coolant_pump_power_W: float = 200.0
    blower_power_W: float = 150.0
    controller_power_W: float = 50.0

    # Thermal
    T_cell_K: float = 353.15  # Cell temperature
    T_ambient_K: float = 298.15  # Ambient temperature


def model_stack(config: StackConfig) -> StackResult:
    """
        Compute stack-level performance, weight, volume, cost, and efficiency.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    logger.info(f"Modeling {config.n_cells}-cell fuel cell stack")

    # ── Electrochemical Performance ─────────────────────────────────────────
    cell_current_A = config.current_density_A_cm2 * config.active_area_cm2
    cell_power_W = config.cell_voltage_V * cell_current_A

    stack_voltage_V = config.cell_voltage_V * config.n_cells
    stack_power_W = cell_power_W * config.n_cells
    stack_power_kW = stack_power_W / 1000.0

    # ── Thermal Management ──────────────────────────────────────────────────
    # HHV thermoneutral voltage E_tn = -ΔH_HHV / (2F) = 285830 J/mol / (2 * 96485 C/mol) ≈ 1.481 V
    E_thermo = 1.48  # V (at 80°C, liquid-water product reference)
    heat_per_cell_W = (E_thermo - config.cell_voltage_V) * cell_current_A
    total_heat_W = heat_per_cell_W * config.n_cells
    total_heat_kW = total_heat_W / 1000.0

    # Radiator sizing (simplified)
    delta_T = config.T_cell_K - config.T_ambient_K
    radiator_area_m2 = total_heat_W / (50.0 * delta_T)  # U ≈ 50 W/(m²·K) forced-air heat transfer coeff

    # ── Balance of Plant Power ──────────────────────────────────────────────
    # Isentropic air compressor work: W_comp = (m_dot_air * c_p * T_in / η_c) * (PR^((γ-1)/γ) - 1)
    gamma_air = 1.4  # Diatomic ideal-gas heat capacity ratio c_p / c_v
    cp_air = 1005.0  # J/(kg·K) dry air specific heat capacity
    air_stoichiometry = 2.0  # λ_air = 2.0× stoichiometric cathode O₂ supply
    m_dot_O2 = cell_current_A * config.n_cells / (4 * F_const) * 32e-3  # kg/s O₂ (M_O2 = 32e-3 kg/mol, 4 e⁻/O₂)
    m_dot_air = m_dot_O2 / 0.233  # dry air is 23.3 wt% O₂ (21 mol% O₂)
    m_dot_air *= air_stoichiometry

    PR = config.compressor_pressure_ratio
    T_in = config.T_ambient_K
    compressor_power_W = (
        m_dot_air
        * cp_air
        * T_in
        / config.compressor_efficiency
        * (PR ** ((gamma_air - 1) / gamma_air) - 1)
    )

    bop_power_W = (
        compressor_power_W
        + config.coolant_pump_power_W
        + config.blower_power_W
        + config.controller_power_W
    )
    bop_power_kW = bop_power_W / 1000.0

    # Net power
    net_power_W = stack_power_W - bop_power_W
    net_power_kW = net_power_W / 1000.0

    # System efficiency
    H2_consumed_mol_s = cell_current_A * config.n_cells / (2 * F_const)
    H2_energy_input_W = H2_consumed_mol_s * H2_LHV_J_MOL  # LHV of H₂ = 241.8 kJ/mol
    system_efficiency = net_power_W / H2_energy_input_W if H2_energy_input_W > 0 else 0
    # Gross stack efficiency on the same LHV basis (V / 1.253 V). Using the
    # HHV thermoneutral 1.48 V here made net system efficiency exceed gross.
    stack_efficiency = config.cell_voltage_V / LHV_VOLTAGE_V

    # ── Weight & Volume ─────────────────────────────────────────────────────
    mea_total_kg = config.mea_mass_g * config.n_cells / 1000.0
    plate_total_kg = config.bipolar_plate_mass_g * config.n_cells / 1000.0
    endplate_kg = config.endplate_mass_g / 1000.0
    gasket_kg = config.gasket_mass_g * config.n_cells / 1000.0
    stack_mass_kg = mea_total_kg + plate_total_kg + endplate_kg + gasket_kg

    # Add BOP estimate (compressor, piping, etc.)
    bop_mass_kg = stack_mass_kg * 0.4  # BOP ≈ 40% of stack mass
    total_mass_kg = stack_mass_kg + bop_mass_kg

    # Volume
    stack_length_m = config.cell_pitch_mm * 1e-3 * config.n_cells + 0.05  # + endplates
    stack_cross_section_m2 = config.active_area_cm2 * 1e-4 * 1.3  # 30% margin
    stack_volume_L = stack_length_m * stack_cross_section_m2 * 1000.0
    total_volume_L = stack_volume_L * 1.5  # BOP adds ~50%

    # Power densities
    gravimetric_W_kg = net_power_W / total_mass_kg if total_mass_kg > 0 else 0
    volumetric_W_L = net_power_W / total_volume_L if total_volume_L > 0 else 0

    # ── Cost Analysis ───────────────────────────────────────────────────────
    # Catalyst cost
    pgm_mass_g = (
        config.pgm_loading_mg_cm2 * config.active_area_cm2 * config.n_cells / 1000.0
    )
    catalyst_cost_usd = pgm_mass_g * config.pgm_price_usd_g
    membrane_cost_usd = (
        config.membrane_cost_usd_cm2 * config.active_area_cm2 * config.n_cells
    )
    plate_cost_usd = config.n_cells * 5.0  # ~$5/plate at volume
    assembly_cost_usd = config.n_cells * 1.0  # ~$1/cell assembly
    bop_cost_usd = 500.0 + bop_power_kW * 100  # base + scaled

    total_cost_usd = (
        catalyst_cost_usd
        + membrane_cost_usd
        + plate_cost_usd
        + assembly_cost_usd
        + bop_cost_usd
    )
    cost_per_kW = total_cost_usd / net_power_kW if net_power_kW > 0 else float('inf')

    result: StackResult = {
        # Stack performance
        'n_cells': config.n_cells,
        'stack_voltage_V': float(stack_voltage_V),
        'stack_current_A': float(cell_current_A),
        'stack_power_kW': float(stack_power_kW),
        'net_power_kW': float(net_power_kW),
        'bop_power_kW': float(bop_power_kW),
        # Thermal
        'heat_rejection_kW': float(total_heat_kW),
        'radiator_area_m2': float(radiator_area_m2),
        # Efficiency
        'stack_efficiency': float(stack_efficiency),
        'system_efficiency': float(system_efficiency),
        'H2_consumption_g_s': float(H2_consumed_mol_s * 2.016),
        # Weight & Volume
        'stack_mass_kg': float(stack_mass_kg),
        'total_mass_kg': float(total_mass_kg),
        'stack_volume_L': float(stack_volume_L),
        'total_volume_L': float(total_volume_L),
        'gravimetric_W_kg': float(gravimetric_W_kg),
        'volumetric_W_L': float(volumetric_W_L),
        # Cost
        'catalyst_cost_usd': float(catalyst_cost_usd),
        'membrane_cost_usd': float(membrane_cost_usd),
        'total_cost_usd': float(total_cost_usd),
        'cost_per_kW': float(cost_per_kW),
    }

    logger.info(f"  Stack: {stack_power_kW:.1f} kW gross, {net_power_kW:.1f} kW net")
    logger.info(
        f"  System η: {system_efficiency:.1%}, {gravimetric_W_kg:.0f} W/kg, {volumetric_W_L:.0f} W/L"
    )
    logger.info(f"  Cost: ${total_cost_usd:.0f} (${cost_per_kW:.0f}/kW)")

    save_json(result, f"stack_{config.n_cells}cell.json", subdir="fuel_cell")
    return result


def couple_reactor_to_pemfc_stack(
    reactor_result: Mapping[str, object],
    *,
    cathode_catalyst: str = 'b5e764d896815bfcfca7',
    orr_overpotential_V: float = 0.2583,
    orr_tafel_slope_mV_dec: float = 78.0,
    cathode_loading_mg_cm2: float = 0.0,
    anode_loading_mg_cm2: float = 0.05,
    n_cells: int = 400,
    psa_h2_recovery: float = 0.90,
    psa_co_rejection: float = 0.99985,
    psa_c2h4_rejection: float = 0.9995,
    feedstock: str | None = None,
) -> CoupledReactorPEMFCResult:
    """Couple a reactor effluent record through PSA purification into a PEMFC stack.

    Models:
      1. PSA purification of reactor product gas (`exit_x_H2`, `exit_x_CO`,
         `exit_x_C2H4`, `C2_carbon_selectivity`) with realistic CO and C2H4
         slip (`post_psa_co_ppm`, `post_psa_c2h4_ppm`, `post_psa_h2_purity`).
      2. Anode Pt HOR poisoning by trace CO and ethylene slip via
         `simulate_pemfc` (`hydrogen_purity`, `co_ppm`).
      3. N-cell PEMFC stack net electrical power, BOP parasitics, well-to-H2
         LHV efficiency, well-to-power LHV efficiency, and total LCOE ($/kWh_e).

    Args:
        reactor_result: Reactor simulation dictionary or pathway leaderboard row.
        cathode_catalyst: PEMFC cathode catalyst identifier.
        orr_overpotential_V: Grand-canonical ORR overpotential in volts.
        orr_tafel_slope_mV_dec: Cathode ORR Tafel slope in mV/decade.
        cathode_loading_mg_cm2: Cathode PGM loading in mg/cm^2.
        anode_loading_mg_cm2: Anode PGM loading in mg/cm^2.
        n_cells: Number of cells in the PEMFC stack.
        psa_h2_recovery: PSA H2 recovery fraction in (0, 1].
        psa_co_rejection: PSA fractional rejection of CO in [0, 1).
        psa_c2h4_rejection: PSA fractional rejection of C2H4/C2 hydrocarbons in [0, 1).
        feedstock: Optional feedstock identifier override (defaults to `reactor_result['feedstock']` or `'CH4'`).

    Returns:
        Dictionary of coupled reactor-PSA-PEMFC performance, impurity slip,
        efficiency, and LCOE metrics.
    """
    from pipeline.economics.tea import LHV_H2_KWH_PER_KG, estimate_scenario_range
    from pipeline.reactors.feedstocks import resolve_feedstock

    feed_key: str = str(
        feedstock
        if feedstock is not None
        else (reactor_result.get('feedstock') or 'CH4')
    )
    spec = resolve_feedstock(feed_key)
    rtype: str = str(reactor_result.get('reactor_type') or 'PFR')
    cat_name: str = str(
        reactor_result.get('champion_catalyst')
        or reactor_result.get('catalyst_name')
        or 'unknown'
    )

    raw_x = (
        reactor_result.get('dft_calibrated_single_pass_CH4_conversion')
        if reactor_result.get('dft_calibrated_single_pass_CH4_conversion') is not None
        else (
            reactor_result.get('single_pass_feedstock_conversion')
            if reactor_result.get('single_pass_feedstock_conversion') is not None
            else (
                reactor_result.get('single_pass_CH4_conversion')
                if reactor_result.get('single_pass_CH4_conversion') is not None
                else reactor_result.get('CH4_conversion')
            )
        )
    )
    x_pass: float = (
        min(max(float(raw_x), 1e-4), 0.9999)
        if isinstance(raw_x, (int, float)) and math.isfinite(float(raw_x)) and float(raw_x) > 0.0
        else 0.50
    )

    raw_h2_sel = reactor_result.get('H2_selectivity')
    h2_sel: float = (
        min(max(float(raw_h2_sel), 0.01), 1.0)
        if isinstance(raw_h2_sel, (int, float)) and math.isfinite(float(raw_h2_sel))
        else 1.0
    )
    raw_c2_sel = reactor_result.get('C2_carbon_selectivity')
    c2_sel: float = (
        min(max(float(raw_c2_sel), 0.0), 1.0)
        if isinstance(raw_c2_sel, (int, float)) and math.isfinite(float(raw_c2_sel))
        else max(0.0, 1.0 - float(reactor_result.get('solid_C_selectivity') or 1.0))
    )
    solid_c_sel: float = max(0.0, 1.0 - c2_sel)
    h2_yield: float = x_pass * h2_sel

    # Effluent molar ratios relative to H2 before PSA:
    raw_x_h2 = reactor_result.get('exit_x_H2')
    x_h2_out: float = (
        float(raw_x_h2)
        if isinstance(raw_x_h2, (int, float)) and float(raw_x_h2) > 1e-6
        else max(
            (spec.nu_H2 * 0.95 * x_pass * h2_sel)
            / max(1.0 + (spec.nu_H2 + spec.nu_CO - 1.0) * 0.95 * x_pass, 1e-6),
            1e-4,
        )
    )
    raw_x_co = reactor_result.get('exit_x_CO')
    x_co_out: float = (
        float(raw_x_co)
        if isinstance(raw_x_co, (int, float)) and float(raw_x_co) >= 0.0
        else (
            (spec.nu_CO * 0.95 * x_pass)
            / max(1.0 + (spec.nu_H2 + spec.nu_CO - 1.0) * 0.95 * x_pass, 1e-6)
        )
    )
    raw_x_c2h4 = reactor_result.get('exit_x_C2H4')
    x_c2h4_out: float = (
        float(raw_x_c2h4)
        if isinstance(raw_x_c2h4, (int, float)) and float(raw_x_c2h4) >= 0.0
        else (0.65 * 0.5 * float(spec.nu_C) * 0.95 * x_pass * c2_sel)
    )

    co_to_h2_raw: float = x_co_out / max(x_h2_out, 1e-6)
    c2h4_to_h2_raw: float = x_c2h4_out / max(x_h2_out, 1e-6)

    # Post-PSA slip in ppm (molar basis in purified H2 stream)
    co_rej: float = min(max(float(psa_co_rejection), 0.90), 0.99999)
    c2h4_rej: float = min(max(float(psa_c2h4_rejection), 0.90), 0.99999)
    h2_rec: float = min(max(float(psa_h2_recovery), 0.50), 0.99)

    post_psa_co_ppm: float = min((co_to_h2_raw * (1.0 - co_rej) / h2_rec) * 1e6, 50.0)
    post_psa_c2h4_ppm: float = min(
        (c2h4_to_h2_raw * (1.0 - c2h4_rej) / h2_rec) * 1e6, 100.0
    )
    # Baseline N2/Ar/CH4 trace slip (~30 ppm for ISO 14687 Grade D H2 = 99.997%) plus CO and C2H4 slip
    post_psa_h2_purity: float = min(
        max(1.0 - 30e-6 - (post_psa_co_ppm + post_psa_c2h4_ppm) * 1e-6, 0.990),
        0.99999,
    )
    # Effective CO-equivalent poisoning on Pt anode (ethylene adsorbs competitively on Pt with ~0.15x CO strength)
    effective_co_ppm: float = post_psa_co_ppm + 0.15 * post_psa_c2h4_ppm

    cell_res = simulate_pemfc(
        PEMFCConfig(
            T_K=353.15,
            membrane_name='AquivionE87',
            membrane_thickness_m=30e-6,
            membrane_conductivity_S_m=14.0,
            cathode_catalyst=cathode_catalyst,
            cathode_loading_mg_cm2=cathode_loading_mg_cm2,
            anode_loading_mg_cm2=anode_loading_mg_cm2,
            orr_overpotential_V=orr_overpotential_V,
            orr_tafel_slope_mV_dec=orr_tafel_slope_mV_dec,
            hydrogen_purity=post_psa_h2_purity,
            co_ppm=effective_co_ppm,
            i_max_A_cm2=4.0,
        )
    )
    stack_res = model_stack(
        StackConfig(
            n_cells=n_cells,
            cell_voltage_V=float(cell_res['peak_voltage_V']),
            current_density_A_cm2=float(cell_res['peak_current_A_cm2']),
            pgm_loading_mg_cm2=cathode_loading_mg_cm2 + anode_loading_mg_cm2,
            membrane_cost_usd_cm2=0.025,
        )
    )

    # TEA & Well-to-Power energetic efficiency
    raw_t = reactor_result.get('T_K')
    t_k: float | None = (
        float(raw_t)
        if isinstance(raw_t, (int, float)) and float(raw_t) > 298.15
        else None
    )
    raw_spec_e = (
        reactor_result.get('specific_energy_kWh_kg_H2')
        if reactor_result.get('specific_energy_kWh_kg_H2') is not None
        else reactor_result.get('dft_calibrated_energy_input_kwh_kg_h2_base')
    )
    spec_e: float | None = (
        float(raw_spec_e)
        if isinstance(raw_spec_e, (int, float)) and 0.0 < float(raw_spec_e) < 1e6
        else None
    )

    tea_range = estimate_scenario_range(
        x_pass,
        temperature_K=t_k,
        solid_c_selectivity=solid_c_sel,
        specific_energy_kwh_kg_h2=spec_e,
        feedstock=spec.formula,
    )
    raw_base_lcoh = reactor_result.get('dft_calibrated_tea_base_h2_usd_kg')
    raw_opt_lcoh = reactor_result.get('dft_calibrated_tea_optimistic_h2_usd_kg')
    lcoh_base: float = (
        float(raw_base_lcoh)
        if isinstance(raw_base_lcoh, (int, float))
        and spec.formula == 'CH4'
        and float(raw_base_lcoh) > 0.0
        else float(tea_range['estimates']['base']['h2_cost_usd_kg'])
    )
    lcoh_opt: float = (
        float(raw_opt_lcoh)
        if isinstance(raw_opt_lcoh, (int, float))
        and spec.formula == 'CH4'
        and float(raw_opt_lcoh) > 0.0
        else float(tea_range['estimates']['optimistic']['h2_cost_usd_kg'])
    )
    energy_in_kwh_kg_h2: float = float(
        tea_range['estimates']['base']['energy_input_kwh_kg_h2']
    )
    net_util: float = float(tea_range['estimates']['base']['net_feedstock_utilization'])

    # Feedstock chemical LHV input per kg H2 (~91% of HHV for hydrocarbons)
    feed_lhv_kwh_per_kg_h2: float = (
        (spec.kg_feed_per_kg_h2 / max(net_util * h2_sel, 1e-4))
        * (0.91 * spec.hhv_MJ_kg)
        / 3.6
    )
    # Solid carbon co-product retains chemical energy (32.8 MJ/kg = 9.11 kWh/kg C),
    # so net chemical + thermal/electrical energy consumed for H2 is:
    c_credit_energy_kwh: float = spec.kg_c_per_kg_h2 * solid_c_sel * (32.8 / 3.6)
    net_well_to_h2_input_kwh: float = max(
        feed_lhv_kwh_per_kg_h2 - c_credit_energy_kwh + energy_in_kwh_kg_h2,
        LHV_H2_KWH_PER_KG * 1.05,
    )
    wth2_eff: float = min(
        max((LHV_H2_KWH_PER_KG * h2_rec) / net_well_to_h2_input_kwh, 0.01),
        0.95,
    )
    stack_sys_eff: float = float(stack_res['system_efficiency'])
    wtp_eff: float = wth2_eff * stack_sys_eff

    # Amortized PEMFC stack capital cost ($/kWh_e over 40,000 h industrial stack lifetime at 0.85 capacity factor)
    # H2 transfer price into the PEMFC is floored at $0.01/kg H2 (compression/handling floor)
    # when co-product solid carbon revenue makes net plant LCOH negative.
    stack_capex_usd_kwh_e: float = float(stack_res['cost_per_kW']) / (40000.0 * 0.85)
    lcoe_fuel_base: float = max(lcoh_base, 0.01) / max(
        LHV_H2_KWH_PER_KG * stack_sys_eff, 1e-6
    )
    lcoe_fuel_opt: float = max(lcoh_opt, 0.01) / max(
        LHV_H2_KWH_PER_KG * stack_sys_eff, 1e-6
    )

    return {
        'feedstock': spec.formula,
        'reactor_type': rtype,
        'catalyst_name': cat_name,
        'single_pass_conversion': round(x_pass, 4),
        'h2_selectivity': round(h2_sel, 4),
        'h2_yield': round(h2_yield, 4),
        'psa_h2_recovery': round(h2_rec, 4),
        'post_psa_h2_purity': round(post_psa_h2_purity, 6),
        'post_psa_co_ppm': round(post_psa_co_ppm, 3),
        'post_psa_c2h4_ppm': round(post_psa_c2h4_ppm, 3),
        'anode_hor_activity_retention': round(
            float(cell_res['hydrogen_impurity_factor']), 4
        ),
        'cathode_catalyst': cathode_catalyst,
        'cell_peak_power_W_cm2': round(float(cell_res['peak_power_W_cm2']), 4),
        'cell_rated_power_W_cm2': (
            round(float(cell_res['rated_power_W_cm2']), 4)
            if cell_res['rated_power_W_cm2'] is not None
            else None
        ),
        'stack_cells': int(stack_res['n_cells']),
        'stack_net_power_kW': round(float(stack_res['net_power_kW']), 2),
        'stack_system_efficiency': round(stack_sys_eff, 4),
        'well_to_hydrogen_lhv_efficiency': round(wth2_eff, 4),
        'well_to_power_lhv_efficiency': round(wtp_eff, 4),
        'lcoh_base_usd_kg': round(lcoh_base, 3),
        'lcoh_optimistic_usd_kg': round(lcoh_opt, 3),
        'lcoe_fuel_usd_kwh_e': round(lcoe_fuel_base, 4),
        'lcoe_stack_capex_usd_kwh_e': round(stack_capex_usd_kwh_e, 4),
        'lcoe_total_usd_kwh_e': round(lcoe_fuel_base + stack_capex_usd_kwh_e, 4),
        'lcoe_optimistic_total_usd_kwh_e': round(
            lcoe_fuel_opt + stack_capex_usd_kwh_e, 4
        ),
    }


if __name__ == '__main__':
    from pipeline.utils import print_banner

    print_banner("FUEL CELL STACK MODEL")

    config = StackConfig(
        n_cells=300,
        active_area_cm2=300,
        cell_voltage_V=0.65,
        current_density_A_cm2=1.5,
        pgm_loading_mg_cm2=0.1,
    )
    result = model_stack(config)
    print(f"\nNet power: {result['net_power_kW']:.1f} kW")
    print(f"System efficiency: {result['system_efficiency']:.1%}")
    print(f"Cost: ${result['cost_per_kW']:.0f}/kW")

