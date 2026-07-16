#!/usr/bin/env python3
"""
Fuel Cell Stack-Level Scaling Model.

Scales single-cell PEMFC results to a practical fuel cell stack:
  - Single-cell → N-cell stack voltage & power
  - Thermal management (heat rejection load)
  - System efficiency including BOP (balance of plant)
  - Gravimetric/volumetric power density (W/kg, W/L)
  - Techno-economic analysis ($/kW)
"""

import numpy as np
from typing import Dict, List, Optional
from dataclasses import dataclass

from pipeline.utils import (
    setup_logger, save_json, F_const, R_gas, FUEL_CELL_DIR,
    METAL_PRICE_USD_KG,
)

logger = setup_logger('fc_stack', 'fuel_cell/stack_model.log')


@dataclass
class StackConfig:
    """Fuel cell stack configuration."""
    # Stack geometry
    n_cells: int = 300               # Number of cells in stack
    active_area_cm2: float = 300.0   # Active area per cell

    # Operating point (from single-cell model)
    cell_voltage_V: float = 0.65     # Operating cell voltage
    current_density_A_cm2: float = 1.5  # Operating current density

    # Component weights (kg per cell)
    mea_mass_g: float = 5.0          # MEA mass per cell
    bipolar_plate_mass_g: float = 150.0  # Bipolar plate per cell
    endplate_mass_g: float = 3000.0  # Total endplate mass (both ends)
    gasket_mass_g: float = 3.0       # Gasket per cell

    # Component volumes (cm³ per cell)
    cell_pitch_mm: float = 2.5       # Cell-to-cell pitch (thickness)

    # Catalyst cost
    pgm_loading_mg_cm2: float = 0.1  # Total PGM loading (both sides)
    pgm_price_usd_g: float = 31.0    # Pt price per gram

    # Membrane cost
    membrane_cost_usd_cm2: float = 0.025

    # BOP parameters
    compressor_efficiency: float = 0.75
    compressor_pressure_ratio: float = 1.5
    coolant_pump_power_W: float = 200.0
    blower_power_W: float = 150.0
    controller_power_W: float = 50.0

    # Thermal
    T_cell_K: float = 353.15         # Cell temperature
    T_ambient_K: float = 298.15      # Ambient temperature


def model_stack(config: StackConfig) -> Dict:
    """
    Compute stack-level performance, weight, volume, cost, and efficiency.
    """
    logger.info(f"Modeling {config.n_cells}-cell fuel cell stack")

    # ── Electrochemical Performance ─────────────────────────────────────────
    cell_current_A = config.current_density_A_cm2 * config.active_area_cm2
    cell_power_W = config.cell_voltage_V * cell_current_A

    stack_voltage_V = config.cell_voltage_V * config.n_cells
    stack_power_W = cell_power_W * config.n_cells
    stack_power_kW = stack_power_W / 1000.0

    # ── Thermal Management ──────────────────────────────────────────────────
    # Thermoneutral voltage
    E_thermo = 1.48  # V (at 80°C)
    heat_per_cell_W = (E_thermo - config.cell_voltage_V) * cell_current_A
    total_heat_W = heat_per_cell_W * config.n_cells
    total_heat_kW = total_heat_W / 1000.0

    # Radiator sizing (simplified)
    delta_T = config.T_cell_K - config.T_ambient_K
    radiator_area_m2 = total_heat_W / (50.0 * delta_T)  # U ≈ 50 W/(m²·K)

    # ── Balance of Plant Power ──────────────────────────────────────────────
    # Air compressor power
    gamma_air = 1.4
    cp_air = 1005.0  # J/(kg·K)
    air_stoichiometry = 2.0  # 2× stoichiometric air
    m_dot_O2 = cell_current_A * config.n_cells / (4 * F_const) * 32e-3  # kg/s O₂
    m_dot_air = m_dot_O2 / 0.233  # air is 23.3% O₂ by mass
    m_dot_air *= air_stoichiometry

    PR = config.compressor_pressure_ratio
    T_in = config.T_ambient_K
    compressor_power_W = (m_dot_air * cp_air * T_in / config.compressor_efficiency *
                          (PR ** ((gamma_air - 1) / gamma_air) - 1))

    bop_power_W = (compressor_power_W + config.coolant_pump_power_W +
                   config.blower_power_W + config.controller_power_W)
    bop_power_kW = bop_power_W / 1000.0

    # Net power
    net_power_W = stack_power_W - bop_power_W
    net_power_kW = net_power_W / 1000.0

    # System efficiency
    H2_consumed_mol_s = cell_current_A * config.n_cells / (2 * F_const)
    H2_energy_input_W = H2_consumed_mol_s * 241800.0  # LHV of H₂ = 241.8 kJ/mol
    system_efficiency = net_power_W / H2_energy_input_W if H2_energy_input_W > 0 else 0
    stack_efficiency = config.cell_voltage_V / E_thermo

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
    pgm_mass_g = config.pgm_loading_mg_cm2 * config.active_area_cm2 * config.n_cells / 1000.0
    catalyst_cost_usd = pgm_mass_g * config.pgm_price_usd_g
    membrane_cost_usd = config.membrane_cost_usd_cm2 * config.active_area_cm2 * config.n_cells
    plate_cost_usd = config.n_cells * 5.0  # ~$5/plate at volume
    assembly_cost_usd = config.n_cells * 1.0  # ~$1/cell assembly
    bop_cost_usd = 500.0 + bop_power_kW * 100  # base + scaled

    total_cost_usd = catalyst_cost_usd + membrane_cost_usd + plate_cost_usd + assembly_cost_usd + bop_cost_usd
    cost_per_kW = total_cost_usd / net_power_kW if net_power_kW > 0 else float('inf')

    result = {
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
    logger.info(f"  System η: {system_efficiency:.1%}, {gravimetric_W_kg:.0f} W/kg, {volumetric_W_L:.0f} W/L")
    logger.info(f"  Cost: ${total_cost_usd:.0f} (${cost_per_kW:.0f}/kW)")

    save_json(result, f"stack_{config.n_cells}cell.json", subdir="fuel_cell")
    return result


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
