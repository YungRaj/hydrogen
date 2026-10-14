#!/usr/bin/env python3
# PEM fuel-cell process model.
"""
1D Through-MEA PEMFC Electrochemical Model.

Models a single-cell PEM fuel cell from anode to cathode:
  - Anode: Butler-Volmer HOR kinetics
  - Membrane: Proton conductivity (temperature + hydration dependent)
  - Cathode: Tafel ORR kinetics with mass-transport limiting current
  - Transport: O₂ diffusion through GDL + ionomer film
  - Water: Electro-osmotic drag + back-diffusion balance

Outputs: Polarization curve (V-I), peak power density, efficiency
"""

import numpy as np
from typing import Dict, List
from dataclasses import dataclass

from pipeline.utils import (
    R_gas,
    F_const,
    setup_logger,
    save_json,
)
from pipeline.data_models.fuel_cells import PEMFCResult

logger = setup_logger('pemfc_model', 'fuel_cell/pemfc_simulation.log')


# ═══════════════════════════════════════════════════════════════════════════════
# CLASS-SPECIFIC TAFEL SLOPES (mV/decade)
# ═══════════════════════════════════════════════════════════════════════════════
# The Tafel slope reflects the ORR mechanism, which differs fundamentally
# between catalyst classes. Using a single value distorts power rankings.
# Sources: Nørskov et al. (2004), Jaouen et al. (2025), DOE AMR 2026.

TAFEL_SLOPE_BY_CLASS = {
    # PGM and PGM-derivative classes
    'SolidCatalyst': 65.0,  # Pt-group nanoparticles, direct 4e⁻
    'HEA': 70.0,  # Multi-metallic, mixed mechanism
    'CoreShell': 63.0,  # Strained Pt shell, enhanced 4e⁻
    'Intermetallic': 62.0,  # Ordered alloy, optimized d-band
    # PGM-free single/dual atom
    'SAC': 78.0,  # M-N₄ sites, 2+2e⁻ at low density → 4e⁻ at high
    'DAC': 72.0,  # Dual sites facilitate 4e⁻ pathway
    'SAA': 66.0,  # Single-atom alloy, near-PGM mechanism
    # Oxide/framework classes
    'Perovskite': 110.0,  # Bulk oxygen diffusion mechanism
    'Spinel': 105.0,  # AB₂O₄, peroxide-mediated
    'MOF': 85.0,  # Metal-organic, variable mechanism
    'COF': 90.0,  # Covalent-organic, limited active sites
    # Carbon-based
    'MetalFreeCarbon': 130.0,  # Intrinsic N-doped carbon, surface-mediated 2e⁻
    # Not typically used as ORR cathodes but included for completeness
    'MoltenMetal': 70.0,  # N/A for ORR, default
    'MetalHydride': 70.0,  # N/A for ORR, default
    'MAXPhase': 85.0,  # Carbide surface sites
    'MXene': 80.0,  # 2D carbide, functionalized surface
}


# ═══════════════════════════════════════════════════════════════════════════════
# PEMFC CELL CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════════════


@dataclass
class PEMFCConfig:
    """Single-cell PEMFC operating configuration."""

    # Operating conditions
    T_K: float = 353.15  # Cell temperature (80°C typical)
    P_anode_Pa: float = 150000.0  # Anode pressure (1.5 atm)
    P_cathode_Pa: float = 150000.0  # Cathode pressure (1.5 atm)
    RH_anode: float = 1.0  # Relative humidity anode
    RH_cathode: float = 1.0  # Relative humidity cathode

    # Membrane properties
    membrane_name: str = 'Nafion_212'
    membrane_thickness_m: float = 50e-6  # 50 μm
    membrane_conductivity_S_m: float = 10.0  # Proton conductivity

    # Catalyst layer
    cathode_catalyst: str = 'Pt/C'
    cathode_loading_mg_cm2: float = 0.1  # Pt loading
    cathode_roughness_factor: float = 100.0  # ECSA multiplier
    catalyst_layer_thickness_m: float = 10e-6
    catalyst_layer_oxygen_resistance_s_m: float = 20.0
    ionomer_fraction: float = 0.30
    orr_overpotential_V: float = 0.35  # From DFT screening
    orr_tafel_slope_mV_dec: float = 70.0  # Tafel slope

    anode_catalyst: str = 'Pt/C'
    anode_loading_mg_cm2: float = 0.05
    hor_exchange_current_A_cm2: float = 0.5  # HOR is fast on Pt
    hydrogen_purity: float = 0.99997
    co_ppm: float = 0.0
    sulfur_ppb: float = 0.0
    voltage_degradation_uV_h: float | None = None

    # GDL properties
    gdl_thickness_m: float = 200e-6  # 200 μm
    gdl_porosity: float = 0.7
    gdl_tortuosity: float = 1.5

    # Current range for polarization curve
    i_max_A_cm2: float = 3.0
    n_points: int = 100


# ═══════════════════════════════════════════════════════════════════════════════
# EFFICIENCY BASIS AND RATED POINT
# ═══════════════════════════════════════════════════════════════════════════════

# Lower heating value of H₂ (J/mol) and its equivalent cell voltage
# ΔH_LHV / (2F) ≈ 1.253 V. Cell and stack efficiencies both use this basis.
H2_LHV_J_MOL = 241800.0
LHV_VOLTAGE_V = H2_LHV_J_MOL / (2.0 * F_const)
RATED_VOLTAGE_V = 0.6


def _rated_point(j_array, V_array, voltage: float):
    """Interpolate the current at ``voltage`` on a falling polarization curve.

    Returns ``(current, power, efficiency)``, or ``(None, None, None)`` when the
    swept curve never crosses ``voltage`` (it stays above it up to i_max, or is
    already below it at the first point). Picking the nearest point in that
    case reports a "rated" current at which the cell is not at ``voltage``.
    """
    below = np.nonzero(np.asarray(V_array) <= voltage)[0]
    if len(below) == 0 or (below[0] == 0 and V_array[0] < voltage):
        return None, None, None
    index = int(below[0])
    if index == 0:
        current = float(j_array[0])
    else:
        v_hi, v_lo = float(V_array[index - 1]), float(V_array[index])
        j_lo, j_hi = float(j_array[index - 1]), float(j_array[index])
        fraction = (v_hi - voltage) / (v_hi - v_lo) if v_hi != v_lo else 0.0
        current = j_lo + fraction * (j_hi - j_lo)
    return current, current * voltage, voltage / LHV_VOLTAGE_V


# ═══════════════════════════════════════════════════════════════════════════════
# PEMFC ELECTROCHEMICAL MODEL
# ═══════════════════════════════════════════════════════════════════════════════


def nernst_voltage(T_K: float, P_H2: float, P_O2: float) -> float:
    """
        Compute the Nernst open-circuit voltage.
        E = E₀ + (RT/2F) ln(P_H2 * P_O2^0.5 / P_H2O)

    Args:
        T_K: Absolute temperature in kelvin.
        P_H2: Hydrogen partial pressure in the configured pressure units.
        P_O2: Oxygen partial pressure in the configured pressure units.

    Returns:
        Computed `float` value in the units documented above.
    """
    E0 = 1.229 - 0.85e-3 * (T_K - 298.15)  # temperature correction
    E = E0 + (R_gas * T_K / (2 * F_const)) * np.log(P_H2 * np.sqrt(P_O2))
    return E


def cathode_activation_loss(j: float, j0_cathode: float, tafel_slope_V: float) -> float:
    """
        Cathode activation overpotential using Tafel equation.
        η_act = (b / ln10) * ln(j / j0)   for j > j0
        where b is the Tafel slope in V/decade.

    Args:
        j: J used by this operation.
        j0_cathode: J0 cathode used by this operation.
        tafel_slope_V: Tafel slope v used by this operation.

    Returns:
        Computed `float` value in the units documented above.
    """
    if j <= 0 or j0_cathode <= 0:
        return 0.0
    b_natural = tafel_slope_V / np.log(10.0)  # convert V/decade → V for use with ln
    return b_natural * np.log(max(j / j0_cathode, 1.0))


def anode_activation_loss(j: float, j0_anode: float, T_K: float) -> float:
    """
        Anode activation overpotential (Butler-Volmer, symmetric).
        For HOR on Pt, this is very small.
        η_act = (RT/αF) * arcsinh(j / 2j0)

    Args:
        j: J used by this operation.
        j0_anode: J0 anode used by this operation.
        T_K: Absolute temperature in kelvin.

    Returns:
        Computed `float` value in the units documented above.
    """
    alpha = 0.5
    return (R_gas * T_K / (alpha * F_const)) * np.arcsinh(j / (2.0 * j0_anode))


def ohmic_loss(j: float, R_ohmic_ohm_cm2: float) -> float:
    """Ohmic loss = j * R_total.

    Args:
        j: J used by this operation.
        R_ohmic_ohm_cm2: R ohmic ohm cm2 used by this operation.

    Returns:
        Computed `float` value in the units documented above.
    """
    return j * R_ohmic_ohm_cm2


def mass_transport_loss(j: float, j_L: float) -> float:
    """
        Concentration (mass-transport) overpotential.
        η_conc = c * ln(1 - j/j_L)
        where c ≈ RT/(4F) and j_L is the limiting current density.

    Args:
        j: J used by this operation.
        j_L: J l used by this operation.

    Returns:
        Computed `float` value in the units documented above.
    """
    if j >= j_L * 0.99:
        return 2.0  # effectively infinite loss
    if j <= 0:
        return 0.0
    c = 0.03  # empirical constant (V)
    return -c * np.log(max(1.0 - j / j_L, 0.01))


def oxygen_transport_resistance_s_m(config: PEMFCConfig) -> float:
    """Return the series O₂ transport resistance from channel to catalyst site.

        R_O2 = δ_GDL / D_eff(T, P) + R_CL

        D_O2 in air is 2.1e-5 m²/s at 293.15 K and 1 atm, scaled as T^1.75/P
        (binary gas diffusivity is inversely proportional to pressure), then
        reduced by porosity/tortuosity. ``catalyst_layer_oxygen_resistance_s_m``
        is already a transport resistance (s/m) and adds in series.

    Args:
        config: Configuration controlling this operation.

    Returns:
        Total oxygen transport resistance in s/m.
    """
    D_O2_bulk = (
        2.1e-5
        * (config.T_K / 293.15) ** 1.75
        * (101325.0 / config.P_cathode_Pa)
    )  # m²/s
    D_O2_eff = D_O2_bulk * config.gdl_porosity / config.gdl_tortuosity
    return (
        config.gdl_thickness_m / D_O2_eff
        + config.catalyst_layer_oxygen_resistance_s_m
    )


def compute_limiting_current(config: PEMFCConfig) -> float:
    """
        Compute the mass-transport limiting current density.
        j_L = 4F * c_O2 / R_O2,   R_O2 = δ_GDL / D_eff(T, P) + R_CL

    Args:
        config: Configuration controlling this operation.

    Returns:
        Computed `float` value in the units documented above.
    """
    # O₂ concentration at cathode inlet
    x_O2 = 0.21  # mole fraction in air
    c_O2 = x_O2 * config.P_cathode_Pa / (R_gas * config.T_K)  # mol/m³

    # Limiting current (A/m² → A/cm²)
    j_L = 4 * F_const * c_O2 / oxygen_transport_resistance_s_m(config)
    j_L_cm2 = j_L * 1e-4  # convert m² to cm²
    return j_L_cm2


def simulate_pemfc(config: PEMFCConfig) -> PEMFCResult:
    """
        Compute the full polarization curve for a PEMFC cell.

        Returns dict with:
          - current_density: array (A/cm²)
          - voltage: array (V)
          - power_density: array (W/cm²)
          - peak_power_density: W/cm²
          - efficiency_at_peak: fraction
          - OCV: open circuit voltage (V)

    Args:
        config: Configuration controlling this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    logger.info(f"Simulating PEMFC: {config.cathode_catalyst} at {config.T_K:.0f} K")

    # Open circuit voltage
    P_H2 = config.P_anode_Pa / 101325.0  # in atm
    P_O2 = 0.21 * config.P_cathode_Pa / 101325.0  # partial pressure O₂
    OCV = nernst_voltage(config.T_K, P_H2, P_O2)

    # Cathode exchange current density (from ORR overpotential)
    # j0_cathode = j_ref * exp(-η_ref / b_natural)
    # where j_ref = 1e-3 A/cm² is the reference current density at which the overpotential was evaluated.
    tafel_slope = config.orr_tafel_slope_mV_dec * 1e-3  # V/decade → V
    b_natural = tafel_slope / np.log(10.0)
    if (
        not 0 < config.hydrogen_purity <= 1
        or config.co_ppm < 0
        or config.sulfur_ppb < 0
    ):
        raise ValueError('invalid hydrogen impurity specification')
    j0_cathode = (
        1e-3
        * np.exp(-config.orr_overpotential_V / b_natural)
        * config.cathode_roughness_factor
    )

    # Anode exchange current density
    # Reversible empirical poisoning factor; prospective claims still require
    # the impurity tests enforced by campaign readiness.
    impurity_factor = config.hydrogen_purity * np.exp(
        -0.08 * config.co_ppm - 0.002 * config.sulfur_ppb
    )
    j0_anode = config.hor_exchange_current_A_cm2 * impurity_factor

    # Ohmic resistance
    R_membrane = (
        config.membrane_thickness_m / config.membrane_conductivity_S_m * 1e4
    )  # Ω·cm²
    R_electronic = 0.005  # Ω·cm² (bipolar plates, GDL)
    R_contact = 0.01  # Ω·cm² (contact resistance)
    R_total = R_membrane + R_electronic + R_contact
    # The catalyst-layer O₂ resistance (s/m) is a transport resistance; it is
    # applied in series with the GDL inside compute_limiting_current, not as an
    # area-specific ohmic term.
    R_O2 = oxygen_transport_resistance_s_m(config)

    # Limiting current
    j_L = compute_limiting_current(config)

    # Current density sweep
    j_upper = min(config.i_max_A_cm2, j_L * 0.98)
    j_array = np.linspace(0.001, j_upper, config.n_points)
    V_array = np.zeros_like(j_array)
    P_array = np.zeros_like(j_array)

    for i, j in enumerate(j_array):
        eta_cathode = cathode_activation_loss(j, j0_cathode, tafel_slope)
        eta_anode = anode_activation_loss(j, j0_anode, config.T_K)
        eta_ohmic = ohmic_loss(j, R_total)
        eta_mt = mass_transport_loss(j, j_L)

        V = OCV - eta_cathode - eta_anode - eta_ohmic - eta_mt
        V = max(V, 0.0)
        V_array[i] = V
        P_array[i] = j * V

    # Find peak power
    peak_idx = int(np.argmax(P_array))
    peak_power = P_array[peak_idx]
    peak_current = j_array[peak_idx]
    peak_voltage = V_array[peak_idx]
    # When the maximum sits on the last sweep point and the sweep was cut by
    # i_max (not by j_L), the reported "peak" is only a lower bound set by the
    # sweep range, not a true maximum of the polarization curve.
    peak_limited_by_sweep = bool(
        peak_idx == len(P_array) - 1 and config.i_max_A_cm2 < j_L * 0.98
    )
    if peak_limited_by_sweep:
        logger.warning(
            f"  Peak power for {config.cathode_catalyst} lies at the sweep limit "
            f"i_max={config.i_max_A_cm2} A/cm² (j_L={j_L:.2f} A/cm²); "
            'it is a lower bound, not a curve maximum'
        )

    # Voltage efficiency on the H₂ LHV basis: E_LHV = ΔH_LHV / (2F) ≈ 1.253 V.
    # The stack's system efficiency uses the same LHV, so cell and system
    # efficiencies are directly comparable.
    efficiency_peak = peak_voltage / LHV_VOLTAGE_V

    # Rated point at 0.6 V, only when the curve actually crosses 0.6 V.
    rated_current, rated_power, efficiency_rated = _rated_point(
        j_array, V_array, RATED_VOLTAGE_V)

    result: PEMFCResult = {
        'cathode_catalyst': config.cathode_catalyst,
        'membrane': config.membrane_name,
        'T_K': config.T_K,
        'OCV_V': float(OCV),
        'peak_power_W_cm2': float(peak_power),
        'peak_current_A_cm2': float(peak_current),
        'peak_voltage_V': float(peak_voltage),
        'peak_power_limited_by_sweep': peak_limited_by_sweep,
        'efficiency_at_peak': float(efficiency_peak),
        'efficiency_basis': 'LHV',
        'hydrogen_impurity_factor': float(impurity_factor),
        'oxygen_transport_resistance_s_m': float(R_O2),
        'voltage_degradation_uV_h': config.voltage_degradation_uV_h,
        'evidence_level': 'modeled',
        'requires_mea_validation': True,
        'rated_power_W_cm2': rated_power,
        'rated_current_A_cm2': rated_current,
        'efficiency_at_rated': efficiency_rated,
        'limiting_current_A_cm2': float(j_L),
        'R_ohmic_ohm_cm2': float(R_total),
        'orr_overpotential_V': config.orr_overpotential_V,
        'tafel_slope_mV_dec': config.orr_tafel_slope_mV_dec,
        'current_density': j_array.tolist(),
        'voltage': V_array.tolist(),
        'power_density': P_array.tolist(),
    }

    logger.info(
        f"  Peak power: {peak_power:.4f} W/cm² at {peak_current:.2f} A/cm²  "
        f"η_peak={efficiency_peak:.1%}  OCV={OCV:.3f} V"
    )

    save_json(
        result,
        f"pemfc_{config.cathode_catalyst}_{config.membrane_name}.json",
        subdir="fuel_cell",
    )
    return result


def sweep_membranes(
    cathode_name: str,
    orr_eta: float,
    membranes: List[Dict] | None = None,
    material_class: str | None = None,
) -> list[PEMFCResult]:
    """Sweep membrane types for a given cathode catalyst.

        If material_class is provided, uses the class-specific Tafel slope
        from TAFEL_SLOPE_BY_CLASS instead of the default 70 mV/dec.

    Args:
        cathode_name: Cathode name used by this operation.
        orr_eta: Orr eta used by this operation.
        membranes: Ordered values supplying membranes.
        material_class: Canonical catalyst material-class name.

    Returns:
        List of computed or validated records.
    """
    if membranes is None:
        from pipeline.screening.fc_cathode_screener import MEMBRANE_TYPES

        membranes = MEMBRANE_TYPES

    tafel = TAFEL_SLOPE_BY_CLASS.get(material_class, 70.0) if material_class else 70.0

    results = []
    for mem in membranes:
        config = PEMFCConfig(
            cathode_catalyst=cathode_name,
            membrane_name=mem['name'],
            membrane_thickness_m=mem['thickness_um'] * 1e-6,
            membrane_conductivity_S_m=mem['conductivity_S_cm'] * 100.0,
            orr_overpotential_V=orr_eta,
            orr_tafel_slope_mV_dec=tafel,
        )
        result = simulate_pemfc(config)
        result['membrane_cost_usd_cm2'] = mem['cost_usd_cm2']
        result['power_per_dollar'] = result['peak_power_W_cm2'] / mem['cost_usd_cm2']
        result['material_class'] = material_class or 'unknown'
        results.append(result)

    return results


if __name__ == '__main__':
    from pipeline.utils import print_banner

    print_banner("PEMFC SINGLE-CELL SIMULATION")

    # Test: simulate a Pt/C cathode with Nafion membrane
    config = PEMFCConfig(
        cathode_catalyst='Pt_C',
        orr_overpotential_V=0.35,
        membrane_name='Nafion_212',
    )
    result = simulate_pemfc(config)
    print(f"\nPeak power: {result['peak_power_W_cm2']:.4f} W/cm²")
    print(f"Efficiency at peak: {result['efficiency_at_peak']:.1%}")
