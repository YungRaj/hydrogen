#!/usr/bin/env python3
"""
Cantera Reactor-Scale Simulation Models for Methane Pyrolysis.

Three reactor archetypes:
  A. Molten Metal Bubble Column Reactor (MMBCR) — CSTR cascade model
  B. Packed-Bed Catalytic Reactor (PFR) — FlowReactor with surface chemistry
  C. Fluidized Bed Reactor — Two-phase bubble/emulsion model

Each model takes a Cantera mechanism file and operating conditions,
and returns conversion, selectivity, and performance metrics.

This script is designed to run in the cp2k-env (Cantera 3.2).
"""

import os
import sys
import json
import numpy as np
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from dataclasses import dataclass, field

# ─── Dynamic import: Cantera may not be in this env ─────────────────────────
try:
    import cantera as ct
    HAS_CANTERA = True
except ImportError:
    HAS_CANTERA = False

# Local imports (handle case where this is run as standalone)
sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    RESULTS_DIR, MECHANISMS_DIR, REACTOR_DIR,
    setup_logger, print_banner, R_gas, save_json,
)

logger = setup_logger('reactor_models', 'reactor/reactor_simulation.log')


# ═══════════════════════════════════════════════════════════════════════════════
# REACTOR CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════

@dataclass
class ReactorConfig:
    """Configuration for reactor simulation."""
    # Operating conditions
    T_inlet_K: float = 1000.0        # Inlet temperature
    P_inlet_Pa: float = 101325.0     # Inlet pressure (1 atm)
    inlet_composition: str = 'CH4:0.95, Ar:0.05'  # Feed composition

    # Bubble column specific
    column_height_m: float = 1.5     # Molten metal column height
    bubble_diameter_mm: float = 5.0  # Average bubble diameter
    gas_velocity_m_s: float = 0.05   # Superficial gas velocity
    n_cstr_stages: int = 20          # Number of CSTR stages for cascade model

    # Packed bed specific
    bed_length_m: float = 0.5        # Catalyst bed length
    bed_diameter_m: float = 0.05     # Bed diameter (lab-scale tube reactor)
    catalyst_particle_mm: float = 2.0  # Catalyst particle diameter
    bed_porosity: float = 0.4        # Void fraction

    # Fluidized bed specific
    u_mf_m_s: float = 0.02          # Minimum fluidization velocity
    bed_height_m: float = 0.8       # Static bed height
    catalyst_density_kg_m3: float = 2500.0  # Catalyst particle density

    # General
    reactor_type: str = 'MMBCR'     # 'MMBCR', 'PFR', 'Fluidized'
    mechanism_file: str = ''
    catalyst_name: str = 'test'
    max_residence_time_s: float = 60.0


# ═══════════════════════════════════════════════════════════════════════════════
# A. MOLTEN METAL BUBBLE COLUMN REACTOR (MMBCR)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_mmbcr(config: ReactorConfig) -> Dict:
    """
    Simulate a molten metal bubble column reactor as a CSTR cascade.
    
    The methane gas enters as bubbles at the bottom of a column of
    molten metal. As bubbles rise, CH₄ decomposes on the gas-liquid
    interface. The CSTR cascade approximates the axial plug-flow
    behavior of the rising bubbles.
    
    Returns dict with conversion profiles, selectivities, etc.
    """
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'MMBCR')

    logger.info(f"Simulating MMBCR: {config.catalyst_name} at {config.T_inlet_K} K")

    # Load mechanism
    gas = ct.Solution(config.mechanism_file, 'gas')
    surf_name = f'{config.catalyst_name}_surface'
    try:
        surf = ct.Interface(config.mechanism_file, surf_name, [gas])
    except Exception:
        # If no surface phase, use gas-phase only
        surf = None

    # Set initial gas state
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    # Compute residence time per stage
    tau_total = config.column_height_m / config.gas_velocity_m_s
    tau_stage = tau_total / config.n_cstr_stages

    # Bubble geometry → surface-to-volume ratio
    d_b = config.bubble_diameter_mm * 1e-3  # m
    sv_ratio = 6.0 / d_b  # sphere S/V = 6/d (m⁻¹)

    # Track axial profiles
    z_positions = np.linspace(0, config.column_height_m, config.n_cstr_stages + 1)
    conversion_profile = [0.0]
    temperature_profile = [config.T_inlet_K]
    species_profiles = {sp: [gas.X[gas.species_index(sp)] if sp in gas.species_names else 0.0]
                        for sp in ['CH4', 'H2', 'C2H2', 'C2H4', 'C2H6']}

    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # CSTR cascade
    for stage in range(config.n_cstr_stages):
        reactor = ct.IdealGasReactor(gas)
        reactor.volume = 1.0  # normalized volume

        if surf is not None:
            rsurf = ct.ReactorSurface(surf, reactor, A=sv_ratio)

        inlet_res = ct.Reservoir(gas)
        outlet_res = ct.Reservoir(gas)

        mdot = gas.density * config.gas_velocity_m_s * np.pi * (d_b / 2) ** 2
        mfc = ct.MassFlowController(inlet_res, reactor, mdot=max(mdot, 1e-8))
        valve = ct.PressureController(reactor, outlet_res, primary=mfc, K=1e-5)

        net = ct.ReactorNet([reactor])
        net.advance(tau_stage)

        # Update gas state for next stage
        gas.TPX = reactor.thermo.T, reactor.thermo.P, reactor.thermo.X

        # Record profiles
        x_ch4 = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 0.0
        conv = 1.0 - x_ch4 / ch4_initial if ch4_initial > 0 else 0.0
        conversion_profile.append(conv)
        temperature_profile.append(gas.T)

        for sp in species_profiles:
            idx = gas.species_index(sp) if sp in gas.species_names else -1
            species_profiles[sp].append(gas.X[idx] if idx >= 0 else 0.0)

    # Compute final metrics
    final_conv = conversion_profile[-1]
    x_h2 = gas.X[gas.species_index('H2')] if 'H2' in gas.species_names else 0.0
    x_c2h2 = gas.X[gas.species_index('C2H2')] if 'C2H2' in gas.species_names else 0.0
    x_c2h4 = gas.X[gas.species_index('C2H4')] if 'C2H4' in gas.species_names else 0.0
    x_c2h6 = gas.X[gas.species_index('C2H6')] if 'C2H6' in gas.species_names else 0.0

    # H₂ selectivity: fraction of H atoms ending up as H₂
    h_in_ch4 = 4.0 * ch4_initial
    h_in_h2 = 2.0 * x_h2
    h2_selectivity = h_in_h2 / max(h_in_ch4 * final_conv, 1e-10) if final_conv > 0.01 else 0.0

    # Carbon selectivity: fraction of C not forming C₂+ species
    c_in_c2_species = 2.0 * (x_c2h2 + x_c2h4 + x_c2h6)
    c_to_solid = final_conv * ch4_initial - c_in_c2_species
    solid_c_selectivity = c_to_solid / max(final_conv * ch4_initial, 1e-10) if final_conv > 0.01 else 0.0

    result = {
        'reactor_type': 'MMBCR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'P_Pa': config.P_inlet_Pa,
        'column_height_m': config.column_height_m,
        'gas_velocity_m_s': config.gas_velocity_m_s,
        'bubble_diameter_mm': config.bubble_diameter_mm,
        'residence_time_s': tau_total,
        'CH4_conversion': float(final_conv),
        'H2_selectivity': float(np.clip(h2_selectivity, 0, 1)),
        'solid_C_selectivity': float(np.clip(solid_c_selectivity, 0, 1)),
        'exit_x_H2': float(x_h2),
        'exit_x_CH4': float(gas.X[gas.species_index('CH4')]) if 'CH4' in gas.species_names else 0.0,
        'exit_x_C2H2': float(x_c2h2),
        'exit_x_C2H4': float(x_c2h4),
        'exit_x_C2H6': float(x_c2h6),
        'exit_T_K': float(gas.T),
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
        'temperature_profile': temperature_profile,
    }

    logger.info(
        f"  MMBCR result: conversion={final_conv:.2%}, "
        f"H2_selectivity={h2_selectivity:.2%}, τ={tau_total:.1f}s"
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# B. PACKED-BED CATALYTIC REACTOR (PFR)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_pfr(config: ReactorConfig) -> Dict:
    """
    Simulate a packed-bed catalytic reactor using Cantera FlowReactor.
    
    Methane flows through a tube filled with catalyst pellets.
    Surface reactions occur on the catalyst surface area.
    """
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'PFR')

    logger.info(f"Simulating PFR: {config.catalyst_name} at {config.T_inlet_K} K")

    gas = ct.Solution(config.mechanism_file, 'gas')
    surf_name = f'{config.catalyst_name}_surface'
    try:
        surf = ct.Interface(config.mechanism_file, surf_name, [gas])
    except Exception:
        surf = None

    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    # Packed-bed surface area per unit volume
    d_p = config.catalyst_particle_mm * 1e-3  # particle diameter in m
    eps = config.bed_porosity
    sv_ratio = 6.0 * (1.0 - eps) / d_p  # m²/m³

    # Use FlowReactor if available (Cantera 3.x), else CSTR cascade
    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # CSTR cascade approximation for PFR
    n_stages = 50
    bed_cross_area = np.pi * (config.bed_diameter_m / 2) ** 2  # m²
    stage_length = config.bed_length_m / n_stages  # m
    stage_volume = bed_cross_area * stage_length * eps  # void volume

    # Superficial velocity
    u_sup = config.gas_velocity_m_s if config.gas_velocity_m_s > 0 else 0.1
    tau_total = config.bed_length_m * eps / u_sup

    z_positions = np.linspace(0, config.bed_length_m, n_stages + 1)
    conversion_profile = [0.0]

    for stage in range(n_stages):
        tau_stage = tau_total / n_stages

        reactor = ct.IdealGasReactor(gas)
        reactor.volume = stage_volume

        if surf is not None:
            rsurf = ct.ReactorSurface(surf, reactor, A=sv_ratio * stage_volume)

        net = ct.ReactorNet([reactor])
        net.advance(tau_stage)
        gas.TPX = reactor.thermo.T, reactor.thermo.P, reactor.thermo.X

        x_ch4 = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 0.0
        conv = 1.0 - x_ch4 / ch4_initial
        conversion_profile.append(conv)

    final_conv = conversion_profile[-1]
    x_h2 = gas.X[gas.species_index('H2')] if 'H2' in gas.species_names else 0.0

    result = {
        'reactor_type': 'PFR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'bed_length_m': config.bed_length_m,
        'bed_diameter_m': config.bed_diameter_m,
        'catalyst_particle_mm': config.catalyst_particle_mm,
        'residence_time_s': tau_total,
        'WHSV_h-1': 3600.0 / tau_total if tau_total > 0 else 0,
        'CH4_conversion': float(final_conv),
        'exit_x_H2': float(x_h2),
        'z_positions': z_positions.tolist(),
        'conversion_profile': conversion_profile,
    }

    logger.info(
        f"  PFR result: conversion={final_conv:.2%}, τ={tau_total:.1f}s"
    )
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# C. FLUIDIZED BED REACTOR (Simplified Two-Phase Model)
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_fluidized_bed(config: ReactorConfig) -> Dict:
    """
    Simplified two-phase (bubble + emulsion) fluidized bed model.
    
    Bubble phase: Plug-flow, gas exchange with emulsion
    Emulsion phase: Well-mixed CSTR at minimum fluidization
    """
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'Fluidized')

    logger.info(f"Simulating Fluidized Bed: {config.catalyst_name} at {config.T_inlet_K} K")

    gas = ct.Solution(config.mechanism_file, 'gas')
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # Two-phase model parameters
    u0 = max(config.gas_velocity_m_s, 0.05)  # operating velocity
    umf = config.u_mf_m_s
    delta = min(0.5, max(0.01, (u0 - umf) / u0))  # bubble fraction of bed

    # Emulsion phase (CSTR)
    tau_emulsion = config.bed_height_m * (1 - delta) / umf

    reactor_em = ct.IdealGasReactor(gas)
    reactor_em.volume = 1.0

    surf_name = f'{config.catalyst_name}_surface'
    try:
        surf = ct.Interface(config.mechanism_file, surf_name, [gas])
        d_p = config.catalyst_particle_mm * 1e-3
        sv_ratio = 6.0 * 0.55 / d_p  # (1-ε_mf)/d_p
        rsurf = ct.ReactorSurface(surf, reactor_em, A=sv_ratio)
    except Exception:
        pass

    net = ct.ReactorNet([reactor_em])
    net.advance(tau_emulsion)

    gas.TPX = reactor_em.thermo.T, reactor_em.thermo.P, reactor_em.thermo.X

    final_x_ch4 = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 0.0
    final_conv = 1.0 - final_x_ch4 / ch4_initial
    x_h2 = gas.X[gas.species_index('H2')] if 'H2' in gas.species_names else 0.0

    result = {
        'reactor_type': 'Fluidized',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'bed_height_m': config.bed_height_m,
        'u0_m_s': u0,
        'umf_m_s': umf,
        'bubble_fraction': delta,
        'residence_time_s': tau_emulsion,
        'CH4_conversion': float(final_conv),
        'exit_x_H2': float(x_h2),
    }

    logger.info(f"  Fluidized result: conversion={final_conv:.2%}")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# MOCK RESULTS (for testing without Cantera)
# ═══════════════════════════════════════════════════════════════════════════════

def _mock_reactor_result(config: ReactorConfig, reactor_type: str) -> Dict:
    """Generate realistic mock results when Cantera is not available."""
    logger.warning(f"Cantera not available. Generating mock {reactor_type} results.")

    # Physics-based estimate using Arrhenius kinetics
    # Assume E_act ~ 0.8 eV for a decent catalyst
    E_act = 0.8  # eV
    k = 1e13 * np.exp(-E_act / (8.617e-5 * config.T_inlet_K))

    if reactor_type == 'MMBCR':
        tau = config.column_height_m / config.gas_velocity_m_s
    elif reactor_type == 'PFR':
        tau = config.bed_length_m * config.bed_porosity / max(config.gas_velocity_m_s, 0.1)
    else:
        tau = config.bed_height_m / max(config.gas_velocity_m_s, 0.05)

    conversion = 1.0 - np.exp(-k * tau * 1e-12)  # scale k appropriately
    conversion = float(np.clip(conversion, 0.01, 0.99))

    return {
        'reactor_type': reactor_type,
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'residence_time_s': tau,
        'CH4_conversion': conversion,
        'H2_selectivity': 0.95,
        'solid_C_selectivity': 0.90,
        'exit_x_H2': conversion * 0.95 * 2.0 / (1.0 + conversion * 0.95),
        'mock': True,
    }


# ═══════════════════════════════════════════════════════════════════════════════
# UNIFIED SIMULATION INTERFACE
# ═══════════════════════════════════════════════════════════════════════════════

def simulate_reactor(config: ReactorConfig) -> Dict:
    """Run the appropriate reactor simulation based on config.reactor_type."""
    simulators = {
        'MMBCR': simulate_mmbcr,
        'PFR': simulate_pfr,
        'Fluidized': simulate_fluidized_bed,
    }
    sim = simulators.get(config.reactor_type, simulate_mmbcr)
    result = sim(config)

    # Save result
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    fname = f"{config.reactor_type}_{config.catalyst_name}_{int(config.T_inlet_K)}K.json"
    save_json(result, fname, subdir='reactor')

    return result


def run_reactor_sweep(catalyst_name: str, mechanism_file: str,
                      temperatures: List[float] = None,
                      reactor_types: List[str] = None) -> List[Dict]:
    """
    Sweep operating conditions for a catalyst across temperatures and reactor types.
    """
    if temperatures is None:
        import os
        mode = os.environ.get('PYROLYSIS_MODE', 'ntec')
        if mode == 'ntec':
            temperatures = [773.15, 800.0, 900.0, 1000.0]
        else:
            temperatures = [1000.0, 1100.0, 1200.0, 1300.0]
    if reactor_types is None:
        reactor_types = ['MMBCR', 'PFR', 'Fluidized']

    results = []
    for rt in reactor_types:
        for T in temperatures:
            config = ReactorConfig(
                T_inlet_K=T,
                reactor_type=rt,
                mechanism_file=str(mechanism_file),
                catalyst_name=catalyst_name,
            )
            result = simulate_reactor(config)
            results.append(result)

    return results


if __name__ == '__main__':
    print_banner("REACTOR SIMULATION TEST")

    # Test with mock data (no Cantera needed)
    config = ReactorConfig(
        T_inlet_K=1000.0,
        reactor_type='MMBCR',
        catalyst_name='NiBi_test',
    )
    result = simulate_reactor(config)
    print(json.dumps(result, indent=2))
