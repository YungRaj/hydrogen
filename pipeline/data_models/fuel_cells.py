"""Explicit records crossing cathode, PEMFC, stack, and reporting boundaries."""

from __future__ import annotations

from typing import Literal, TypedDict


class _PEMFCResultOptional(TypedDict, total=False):
    membrane_cost_usd_cm2: float
    power_per_dollar: float
    material_class: str


class PEMFCResult(_PEMFCResultOptional):
    """One modeled membrane/cathode operating result with unit-bearing names."""

    cathode_catalyst: str
    membrane: str
    T_K: float
    OCV_V: float
    peak_power_W_cm2: float
    peak_current_A_cm2: float
    peak_voltage_V: float
    efficiency_at_peak: float
    hydrogen_impurity_factor: float
    catalyst_layer_resistance_ohm_cm2: float
    voltage_degradation_uV_h: float | None
    evidence_level: Literal['modeled']
    requires_mea_validation: Literal[True]
    rated_power_W_cm2: float
    rated_current_A_cm2: float
    efficiency_at_rated: float
    limiting_current_A_cm2: float
    R_ohmic_ohm_cm2: float
    orr_overpotential_V: float
    tafel_slope_mV_dec: float
    current_density: list[float]
    voltage: list[float]
    power_density: list[float]


class StackResult(TypedDict):
    """Fuel-cell stack performance, thermal, mass, volume, and cost outputs."""

    n_cells: int
    stack_voltage_V: float
    stack_current_A: float
    stack_power_kW: float
    net_power_kW: float
    bop_power_kW: float
    heat_rejection_kW: float
    radiator_area_m2: float
    stack_efficiency: float
    system_efficiency: float
    H2_consumption_g_s: float
    stack_mass_kg: float
    total_mass_kg: float
    stack_volume_L: float
    total_volume_L: float
    gravimetric_W_kg: float
    volumetric_W_L: float
    catalyst_cost_usd: float
    membrane_cost_usd: float
    total_cost_usd: float
    cost_per_kW: float


class EmptyStackResult(TypedDict):
    """No stack was modeled because the PEMFC stage produced no candidates."""
