#!/usr/bin/env python3
# Turquoise-hydrogen reactor process models.
"""
Cantera Reactor-Scale Simulation Models for Methane Pyrolysis.

Three reactor archetypes:
  A. Molten Metal Bubble Column Reactor (MMBCR) — CSTR cascade model
  B. Packed-Bed Catalytic Reactor (PFR) — staged Lagrangian surface model
  C. Fluidized Bed Reactor — emulsion reaction plus bubble bypass model

NTEC and electrochemical pathways consume validated external multiphysics
artifacts. Cantera supplies chemistry within those coupled calculations, but
does not itself solve their mechanical/electrical/charge-transport physics.

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
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.common.utils import (
    RESULTS_DIR, MECHANISMS_DIR, REACTOR_DIR,
    setup_logger, print_banner, R_gas, save_json,
)
from pipeline.process.pathway_modes import (
    DEFAULT_MODE, REACTOR_MODELS, reactor_applicability,
    reactor_types_for_mode, validate_mode_reactors)

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
    column_diameter_m: float = 0.10  # Internal column diameter
    bubble_diameter_mm: float = 5.0  # Average bubble diameter
    gas_velocity_m_s: float = 0.05   # Superficial gas velocity
    gas_holdup_fraction: float = 0.10  # Dispersed gas volume / bed volume
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
    fluidized_bubble_fraction: float | None = None

    # General
    reactor_type: str = 'PFR'
    mechanism_file: str = ''
    catalyst_name: str = 'test'
    material_class: str | None = None
    pathway_mode: str = DEFAULT_MODE
    candidate_id: str = 'unknown'
    multiphysics_results_dir: str | None = None
    multiphysics_artifact: dict | None = None
    max_residence_time_s: float = 60.0
    catalyst_E_act_eV: float = 0.8   # Catalyst activation barrier (used by mock when Cantera unavailable)


def _load_candidate_phases(config: ReactorConfig):
    """Load the gas and required candidate-specific surface phase."""
    gas = ct.Solution(config.mechanism_file, 'gas')
    surf_name = f'{config.catalyst_name}_surface'
    try:
        surface = ct.Interface(config.mechanism_file, surf_name, [gas])
    except Exception as exc:
        raise RuntimeError(
            f'candidate surface phase failed to load: {exc}') from exc
    return gas, surface


def _mechanism_metadata(config: ReactorConfig) -> dict:
    path = Path(config.mechanism_file).with_suffix('.kinetics.json')
    if not path.exists():
        return {'inputs': {'quantitative_status': 'missing_provenance'},
                'carbon_phase_model': 'unknown'}
    try:
        return json.loads(path.read_text())
    except Exception:
        return {'inputs': {'quantitative_status': 'invalid_provenance'},
                'carbon_phase_model': 'unknown'}


def _kinetics_metadata(config: ReactorConfig) -> dict:
    return _mechanism_metadata(config).get('inputs', {})


def _kinetics_evidence(config: ReactorConfig) -> dict:
    """Declare whether reactor output may make a candidate-level decision."""
    metadata = _mechanism_metadata(config)
    status = metadata.get('inputs', {}).get(
        'quantitative_status', 'missing_provenance')
    carbon_model = metadata.get('carbon_phase_model', 'unknown')
    limitations = []
    if status != 'candidate_specific':
        limitations.append('incomplete_candidate_kinetics')
    if carbon_model == 'legacy_gas_tracer':
        limitations.append('legacy_gas_carbon_tracer')
    elif carbon_model == 'unknown':
        limitations.append('unknown_carbon_phase_model')
    complete = not limitations
    return {
        'kinetics_status': status,
        'carbon_phase_model': carbon_model,
        'reactor_evidence_tier': (
            'candidate_specific_kinetics' if complete else
            'diagnostic_screening_template'),
        # Incomplete template kinetics can guide sensitivity/validation but may
        # never eliminate a candidate or count as reactor validation.
        'can_exclude_candidate': bool(complete),
        'reactor_evidence_limitations': limitations,
    }


def _validate_reactor_config(config: ReactorConfig) -> None:
    """Reject geometries that do not represent the selected physical bed."""
    if config.T_inlet_K <= 0 or config.P_inlet_Pa <= 0:
        raise ValueError('reactor temperature and pressure must be positive')
    if config.reactor_type == 'PFR':
        if config.bed_length_m <= 0 or config.bed_diameter_m <= 0 or \
                config.catalyst_particle_mm <= 0 or config.gas_velocity_m_s <= 0:
            raise ValueError('packed-bed dimensions, particle size, and flow must be positive')
        if not 0 < config.bed_porosity < 1:
            raise ValueError('packed-bed porosity must lie strictly between zero and one')
    elif config.reactor_type == 'Fluidized':
        if config.bed_height_m <= 0 or config.catalyst_particle_mm <= 0 or \
                config.u_mf_m_s <= 0:
            raise ValueError('fluidized-bed dimensions, particle size, and umf must be positive')
        if config.gas_velocity_m_s <= config.u_mf_m_s:
            raise ValueError('gas velocity must exceed minimum fluidization velocity')
    elif config.reactor_type == 'MMBCR':
        if config.column_height_m <= 0 or config.column_diameter_m <= 0 or \
                config.bubble_diameter_mm <= 0 or config.gas_velocity_m_s <= 0 or \
                config.n_cstr_stages < 1:
            raise ValueError('MMBCR dimensions, flow, bubbles, and stage count must be positive')
        if not 0 < config.gas_holdup_fraction < 1:
            raise ValueError('MMBCR gas holdup must lie strictly between zero and one')


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
    gas, surf = _load_candidate_phases(config)

    # Set initial gas state
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    # Compute residence time per stage
    tau_total = (config.column_height_m * config.gas_holdup_fraction /
                 config.gas_velocity_m_s)

    # Bubble geometry → surface-to-volume ratio
    d_b = config.bubble_diameter_mm * 1e-3  # m
    interfacial_area_density = (
        6.0 * config.gas_holdup_fraction / d_b)  # m² interface / m³ bed
    column_area = np.pi * (config.column_diameter_m / 2) ** 2
    stage_bed_volume = (column_area * config.column_height_m /
                        config.n_cstr_stages)
    stage_gas_volume = stage_bed_volume * config.gas_holdup_fraction
    stage_interface_area = interfacial_area_density * stage_bed_volume

    # Track axial profiles
    z_positions = np.linspace(0, config.column_height_m, config.n_cstr_stages + 1)
    conversion_profile = [0.0]
    temperature_profile = [config.T_inlet_K]
    species_profiles = {sp: [gas.X[gas.species_index(sp)] if sp in gas.species_names else 0.0]
                        for sp in ['CH4', 'H2', 'C2H2', 'C2H4', 'C2H6']}

    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # CSTR cascade
    for stage in range(config.n_cstr_stages):
        reactor = ct.IdealGasReactor(gas, energy='off')
        reactor.volume = stage_gas_volume

        rsurf = ct.ReactorSurface(surf, reactor, A=stage_interface_area)

        inlet_res = ct.Reservoir(gas)
        outlet_res = ct.Reservoir(gas)

        mdot = gas.density * config.gas_velocity_m_s * column_area
        mfc = ct.MassFlowController(inlet_res, reactor, mdot=max(mdot, 1e-8))
        valve = ct.PressureController(reactor, outlet_res, primary=mfc, K=1e-5)

        net = ct.ReactorNet([reactor])
        net.advance_to_steady_state()

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

    c_in_c2_species = 2.0 * (x_c2h2 + x_c2h4 + x_c2h6)
    c_to_solid = final_conv * ch4_initial - c_in_c2_species
    solid_c_selectivity = c_to_solid / max(final_conv * ch4_initial, 1e-10) if final_conv > 0.01 else 0.0

    result = {
        'reactor_type': 'MMBCR',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'P_Pa': config.P_inlet_Pa,
        'column_height_m': config.column_height_m,
        'column_diameter_m': config.column_diameter_m,
        'gas_velocity_m_s': config.gas_velocity_m_s,
        'gas_holdup_fraction': config.gas_holdup_fraction,
        'bubble_diameter_mm': config.bubble_diameter_mm,
        'residence_time_s': tau_total,
        'CH4_conversion': float(final_conv),
        'H2_selectivity': float(np.clip(h2_selectivity, 0, 1)),
        'solid_C_selectivity': float(np.clip(solid_c_selectivity, 0, 1)),
        **_kinetics_evidence(config),
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
    Simulate a packed bed with a staged Lagrangian PFR approximation.
    
    Methane flows through a tube filled with catalyst pellets.
    Surface reactions occur on the catalyst surface area.
    """
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'PFR')

    logger.info(f"Simulating PFR: {config.catalyst_name} at {config.T_inlet_K} K")

    gas, surf = _load_candidate_phases(config)

    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition

    # Packed-bed surface area per unit volume
    d_p = config.catalyst_particle_mm * 1e-3  # particle diameter in m
    eps = config.bed_porosity
    sv_ratio = 6.0 * (1.0 - eps) / d_p  # m²/m³

    # Follow a gas parcel through sequential surface-reacting control volumes.
    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # CSTR cascade approximation for PFR
    n_stages = 50
    bed_cross_area = np.pi * (config.bed_diameter_m / 2) ** 2  # m²
    stage_length = config.bed_length_m / n_stages  # m
    stage_bed_volume = bed_cross_area * stage_length
    stage_gas_volume = stage_bed_volume * eps

    # Superficial velocity
    u_sup = config.gas_velocity_m_s if config.gas_velocity_m_s > 0 else 0.1
    tau_total = config.bed_length_m * eps / u_sup

    z_positions = np.linspace(0, config.bed_length_m, n_stages + 1)
    conversion_profile = [0.0]

    for stage in range(n_stages):
        tau_stage = tau_total / n_stages

        reactor = ct.IdealGasReactor(gas, energy='off')
        reactor.volume = stage_gas_volume

        rsurf = ct.ReactorSurface(surf, reactor, A=sv_ratio * stage_bed_volume)

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
        **_kinetics_evidence(config),
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
    
    Emulsion phase reacts at minimum-fluidization residence time. The bubble
    fraction is represented as conservative bypass and mixed with the emulsion
    outlet. Interphase mass transfer is not yet resolved, so this remains a
    screening approximation rather than a predictive two-phase CFD model.
    """
    if not HAS_CANTERA:
        return _mock_reactor_result(config, 'Fluidized')

    logger.info(f"Simulating Fluidized Bed: {config.catalyst_name} at {config.T_inlet_K} K")

    gas, surf = _load_candidate_phases(config)
    gas.TPX = config.T_inlet_K, config.P_inlet_Pa, config.inlet_composition
    inlet_x = gas.X.copy()

    ch4_initial = gas.X[gas.species_index('CH4')] if 'CH4' in gas.species_names else 1.0

    # Two-phase model parameters
    u0 = config.gas_velocity_m_s  # operating velocity (> umf is validated)
    umf = config.u_mf_m_s
    delta = (float(config.fluidized_bubble_fraction)
             if config.fluidized_bubble_fraction is not None
             else min(0.5, max(0.01, (u0 - umf) / u0)))
    if not 0 <= delta < 1:
        raise ValueError('fluidized-bed bubble fraction must lie in [0, 1)')

    # Emulsion phase (CSTR)
    tau_emulsion = config.bed_height_m * (1 - delta) / umf

    reactor_em = ct.IdealGasReactor(gas, energy='off')
    reactor_em.volume = 1.0

    d_p = config.catalyst_particle_mm * 1e-3
    sv_ratio = 6.0 * 0.55 / d_p  # (1-ε_mf)/d_p
    rsurf = ct.ReactorSurface(surf, reactor_em, A=sv_ratio)

    net = ct.ReactorNet([reactor_em])
    net.advance(tau_emulsion)

    gas.TPX = reactor_em.thermo.T, reactor_em.thermo.P, reactor_em.thermo.X

    # Conservative two-phase closure: bubbles bypass reaction and mix with the
    # reacted emulsion outlet. This prevents the emulsion proxy from being
    # mislabeled as conversion of the entire gas feed.
    mixed_x = delta * inlet_x + (1.0 - delta) * gas.X
    gas.TPX = gas.T, gas.P, mixed_x

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
        **_kinetics_evidence(config),
        'exit_x_H2': float(x_h2),
    }

    logger.info(f"  Fluidized result: conversion={final_conv:.2%}")
    return result


# ═══════════════════════════════════════════════════════════════════════════════
# MOCK RESULTS (for testing without Cantera)
# ═══════════════════════════════════════════════════════════════════════════════

def _mock_reactor_result(config: ReactorConfig, reactor_type: str) -> Dict:
    """Generate realistic mock results when Cantera is not available.

    Uses the catalyst-specific E_act from config (not a hardcoded value)
    so that mock results still differentiate between catalysts.
    """
    logger.warning(f"Cantera not available. Generating mock {reactor_type} results.")

    # Physics-based estimate using Arrhenius kinetics with actual catalyst E_act
    E_act = config.catalyst_E_act_eV
    k_B_eV = 8.617e-5  # Boltzmann constant in eV/K
    k = 1e13 * np.exp(-E_act / (k_B_eV * config.T_inlet_K))

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
        'catalyst_E_act_eV': E_act,
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

def simulate_ntec_pathway(config: ReactorConfig) -> Dict:
    """Describe NTEC readiness without substituting an unrelated reactor model."""
    from pipeline.process.ntec_model import (
        conditions_from_environment, ntec_assistance)

    assistance = ntec_assistance(conditions_from_environment())
    if config.multiphysics_artifact:
        evidence = config.multiphysics_artifact
        artifact = evidence['artifact']
        outputs = artifact['outputs']
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'NTEC', 'pathway_mode': 'ntec',
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id,
            'T_K': config.T_inlet_K,
            'CH4_conversion': float(outputs['CH4_conversion']),
            'H2_selectivity': float(outputs['H2_selectivity']),
            'solid_C_selectivity': float(outputs['solid_C_selectivity']),
            'specific_energy_kWh_kg_H2': float(
                outputs['specific_energy_kWh_kg_H2']),
            'ntec_assistance': assistance,
            'multiphysics_evidence': evidence,
            'reactor_evidence_tier': 'calibrated_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required',
        'valid': False,
        'reactor_type': 'NTEC',
        'pathway_mode': 'ntec',
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'ntec_assistance': assistance,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_ntec_pathway_kinetics_required',
            'liquid_solid_hydrodynamic_model_required',
        ],
    }


def simulate_electrochemical_pathway(config: ReactorConfig) -> Dict:
    """Report electrochemical evidence without inventing a Cantera conversion."""
    from pipeline.process.electrochemical_model import (
        conditions_from_environment, electrochemical_evidence)

    evidence = electrochemical_evidence(conditions_from_environment())
    phase = evidence['conditions'].get('electrolyte_phase')
    if config.multiphysics_artifact:
        solver_evidence = config.multiphysics_artifact
        artifact = solver_evidence['artifact']
        outputs = artifact['outputs']
        phase = artifact.get('electrolyte_phase', phase)
        return {
            'status': 'complete', 'valid': True,
            'reactor_type': 'Electrochemical',
            'pathway_mode': 'electrochemical',
            'electrolyte_phase': phase,
            'catalyst_name': config.catalyst_name,
            'candidate_id': config.candidate_id,
            'T_K': config.T_inlet_K,
            **{name: float(outputs[name]) for name in (
                'CH4_conversion', 'H2_selectivity',
                'faradaic_efficiency_H2', 'current_density_A_cm2',
                'cell_voltage_V', 'electrical_power_density_W_cm2')},
            'electrochemical_evidence': evidence,
            'multiphysics_evidence': solver_evidence,
            'reactor_evidence_tier': 'mechanistic_multiphysics_screening',
            'can_exclude_candidate': False,
        }
    return {
        'status': 'validation_required',
        'valid': False,
        'reactor_type': 'Electrochemical',
        'pathway_mode': 'electrochemical',
        'electrolyte_phase': phase,
        'catalyst_name': config.catalyst_name,
        'T_K': config.T_inlet_K,
        'electrochemical_evidence': evidence,
        'reactor_evidence_tier': 'pathway_model_pending',
        'can_exclude_candidate': False,
        'limitations': [
            'candidate_specific_electrochemical_kinetics_required',
            f'{phase or "unspecified"}_electrolyte_transport_model_required',
        ],
    }


def simulate_reactor(config: ReactorConfig) -> Dict:
    """Run the appropriate reactor simulation based on config.reactor_type."""
    simulators = {
        'MMBCR': simulate_mmbcr,
        'PFR': simulate_pfr,
        'Fluidized': simulate_fluidized_bed,
        'NTEC': simulate_ntec_pathway,
        'Electrochemical': simulate_electrochemical_pathway,
    }
    if config.reactor_type not in simulators:
        raise ValueError(f'unsupported reactor type: {config.reactor_type}')
    selected_mode = reactor_types_for_mode(config.pathway_mode)
    if config.reactor_type not in selected_mode:
        raise ValueError(
            f'reactor {config.reactor_type!r} is not valid for pathway mode '
            f'{config.pathway_mode!r}; expected {selected_mode}')
    applicable, reason = reactor_applicability(
        config.reactor_type, config.material_class)
    spec = REACTOR_MODELS[config.reactor_type]
    if not applicable:
        result = {
            'status': 'not_applicable', 'valid': False,
            'reactor_type': config.reactor_type,
            'pathway_mode': config.pathway_mode,
            'catalyst_name': config.catalyst_name,
            'material_class': config.material_class,
            'reason': reason, 'can_exclude_candidate': False,
            'reactor_evidence_tier': 'not_applicable_non_excluding',
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
        }
    else:
        from pipeline.process.multiphysics_contract import (
            EXTERNAL_SOLVERS, load_validated_artifact)
        if config.reactor_type in EXTERNAL_SOLVERS:
            loaded = load_validated_artifact(
                config.multiphysics_results_dir, config.candidate_id,
                config.pathway_mode, config.reactor_type, config.T_inlet_K)
            if loaded['valid'] and config.reactor_type == 'Electrochemical':
                from pipeline.process.electrochemical_model import (
                    conditions_from_environment)
                requested = conditions_from_environment()
                artifact_phase = loaded['artifact'].get('electrolyte_phase')
                if (requested.electrolyte_phase is not None and
                        requested.electrolyte_phase.lower() != artifact_phase):
                    loaded = {
                        'valid': False,
                        'reason': 'electrolyte_phase_mismatch',
                        'requested': requested.electrolyte_phase.lower(),
                        'artifact': artifact_phase,
                        'path': loaded['path'],
                    }
            if not loaded['valid']:
                result = {
                    'status': 'validation_required', 'valid': False,
                    'reactor_type': config.reactor_type,
                    'pathway_mode': config.pathway_mode,
                    'catalyst_name': config.catalyst_name,
                    'candidate_id': config.candidate_id,
                    'material_class': config.material_class,
                    'multiphysics_evidence': loaded,
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'required_solver_evidence_missing',
                    'bed_or_interface': spec.bed_or_interface,
                    'cantera_reactor_model': spec.cantera_model,
                    'reaction_domain': spec.reaction_domain,
                    'reactor_model_fidelity': spec.fidelity,
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = (f"{config.reactor_type}_{config.catalyst_name}_"
                         f"{int(config.T_inlet_K)}K.json")
                save_json(result, fname, subdir='reactor')
                return result
            config.multiphysics_artifact = loaded
            outputs = loaded['artifact']['outputs']
            if config.reactor_type == 'Fluidized':
                config.gas_velocity_m_s = float(outputs['gas_velocity_m_s'])
                config.u_mf_m_s = float(outputs['u_mf_m_s'])
                config.fluidized_bubble_fraction = float(
                    outputs['bubble_fraction'])
            elif config.reactor_type == 'MMBCR':
                config.gas_velocity_m_s = float(outputs['gas_velocity_m_s'])
                config.gas_holdup_fraction = float(
                    outputs['gas_holdup_fraction'])
                config.bubble_diameter_mm = float(outputs['bubble_diameter_mm'])
        _validate_reactor_config(config)
        result = simulators[config.reactor_type](config)
        result.setdefault('status', 'complete')
        result.update({
            'pathway_mode': config.pathway_mode,
            'material_class': config.material_class,
            'bed_or_interface': spec.bed_or_interface,
            'cantera_reactor_model': spec.cantera_model,
            'reaction_domain': spec.reaction_domain,
            'reactor_model_fidelity': spec.fidelity,
            'multiphysics_evidence': config.multiphysics_artifact,
        })
        if spec.fidelity != 'validated_predictive':
            limitations = list(result.get('reactor_evidence_limitations', []))
            limitation = f'{config.reactor_type.lower()}_screening_model_not_validated'
            if limitation not in limitations:
                limitations.append(limitation)
            result['reactor_evidence_limitations'] = limitations
            result['can_exclude_candidate'] = False

    # Save result
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    fname = f"{config.reactor_type}_{config.catalyst_name}_{int(config.T_inlet_K)}K.json"
    save_json(result, fname, subdir='reactor')

    return result


def run_reactor_sweep(catalyst_name: str, mechanism_file: str,
                      temperatures: List[float] = None,
                      reactor_types: List[str] = None,
                      catalyst_E_act_eV: float = 0.8,
                      pathway_mode: str = DEFAULT_MODE,
                      material_class: str | None = None,
                      candidate_id: str = 'unknown',
                      multiphysics_results_dir: str | None = None) -> List[Dict]:
    """
    Sweep operating conditions for a catalyst across temperatures and reactor types.

    Args:
        catalyst_E_act_eV: Activation barrier in eV. Passed through to ReactorConfig
            so that mock results (when Cantera is unavailable) still differentiate
            between catalysts.
    """
    if temperatures is None:
        temperatures = [773.15, 900.0, 1100.0, 1300.0]
    if reactor_types is None:
        # The default family is conventional thermocatalysis. MMBCR, NTEC, and
        # electrochemical physics must be selected explicitly.
        reactor_types = list(reactor_types_for_mode(pathway_mode))
    validate_mode_reactors(pathway_mode, reactor_types)

    results = []
    for rt in reactor_types:
        for T in temperatures:
            config = ReactorConfig(
                T_inlet_K=T,
                reactor_type=rt,
                mechanism_file=str(mechanism_file),
                catalyst_name=catalyst_name,
                catalyst_E_act_eV=catalyst_E_act_eV,
                pathway_mode=pathway_mode,
                material_class=material_class,
                candidate_id=candidate_id,
                multiphysics_results_dir=multiphysics_results_dir,
            )
            try:
                result = simulate_reactor(config)
            except Exception as exc:
                # A stiff condition must not discard the other temperatures or
                # reactor types. Preserve it as non-excluding failed evidence.
                result = {
                    'status': 'failed',
                    'valid': False,
                    'reactor_type': rt,
                    'pathway_mode': pathway_mode,
                    'material_class': material_class,
                    'temperature_K': float(T),
                    'catalyst': catalyst_name,
                    'error_type': type(exc).__name__,
                    'error': str(exc),
                    'can_exclude_candidate': False,
                    'reactor_evidence_tier': 'failed_simulation',
                }
                REACTOR_DIR.mkdir(parents=True, exist_ok=True)
                fname = f"{rt}_{catalyst_name}_{int(T)}K.json"
                save_json(result, fname, subdir='reactor')
                logger.error(
                    'Reactor condition failed without excluding candidate: '
                    f'{rt} {catalyst_name} {T} K: {exc}')
            results.append(result)

    return results


if __name__ == '__main__':
    print_banner("REACTOR SIMULATION TEST")

    # Test with mock data (no Cantera needed)
    config = ReactorConfig(
        T_inlet_K=1000.0,
        reactor_type='MMBCR',
        catalyst_name='NiBi_test',
        material_class='MoltenMetal',
        pathway_mode='mmbcr',
    )
    result = simulate_reactor(config)
    print(json.dumps(result, indent=2))
