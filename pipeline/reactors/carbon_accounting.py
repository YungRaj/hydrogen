"""Baker-Helveg-Ermakova-Takenaka off-site carbon accounting and Cantera surface state helpers."""

from __future__ import annotations

from typing import Callable, Mapping, Optional, Protocol, cast

import numpy as np
import numpy.typing as npt

from pipeline.data_models.reactors import (
    CanteraGasPhase,
    CanteraSurfacePhase,
    H2C2SelectivityRateMetrics,
    LoadStatusFields,
)

try:
    import cantera as ct
except ImportError:
    ct = None  # type: ignore[assignment]


class _CarbonConfigLike(Protocol):
    """Minimal reactor configuration contract required by carbon accounting helpers."""

    mechanism_file: str
    catalyst_name: str
    gas_only: bool
    T_inlet_K: float
    P_inlet_Pa: float
    metal_dispersion: float
    metal_loading: float


METAL_MOLAR_MASS_G_MOL: dict[str, float] = {
    'Ni': 58.6934,
    'Fe': 55.845,
    'Co': 58.9332,
}
# Ermakova 2000 / Takenaka: 40-384 gC/gNi over 4-50 h on high-Ni/SiO2 -> ~8-10 gC/(gNi h).
NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H: tuple[float, float] = (8.0, 10.0)
# Same TOS literature: the catalyst dies (encapsulated) after 4-50 h.
NI_TOS_LIFETIME_BAND_H: tuple[float, float] = (4.0, 50.0)
# Baker (1972) / Helveg (2004) / Takenaka (2003) encapsulation-onset markers:
# ENCAP_ONSET_THETA = 0.10 (10% irreversible C_δ site blockage triggers onset warning),
# ENCAP_ONSET_RATIO = 10.0 (filament-to-encapsulating selectivity C_γ/C_δ < 10 signals competitive deactivation),
# ENCAP_DEAD_THETA = 0.50 (50% C_δ site blockage defines end-of-life time-on-stream τ_TOS).
ENCAP_ONSET_THETA: float = 0.1
ENCAP_ONSET_RATIO: float = 10.0
ENCAP_DEAD_THETA: float = 0.5
PFR_PILOT_TAU_BAND_S: tuple[float, float] = (2.0, 12.0)  # Pilot-scale packed-bed gas residence-time window [s]
CH4_PYROLYSIS_DELTA_H_KJ_MOL: float = 74.85  # Standard reaction enthalpy ΔH°_298 for CH4(g) -> C(gr) + 2 H2(g) [kJ/mol]
SOLID_CARBON_DENSITY_KG_M3: float = 2160.0  # Bulk crystalline/turbostratic graphite density [kg/m³] (2.16 g/cm³)

# NASA-7 polynomial coefficients (a1..a6) for Low (200-1000 K) and High (1000-5000/6000 K) temperature intervals:
# H°(T)/(R T) = a1 + a2*T/2 + a3*T²/3 + a4*T³/4 + a5*T⁴/5 + a6/T
_NASA7_CH4_LOW = (
    5.14987613,
    -0.0136709788,
    4.91800599e-05,
    -4.84743026e-08,
    1.66693956e-11,
    -10246.6,
)
_NASA7_CH4_HIGH = (
    0.074851495,
    0.0133909467,
    -5.73285809e-06,
    1.22292535e-09,
    -1.01815230e-13,
    -9468.34459,
)
_NASA7_H2_LOW = (
    2.34433112,
    7.98052075e-03,
    -1.9478151e-05,
    2.01572094e-08,
    -7.37611761e-12,
    -917.935173,
)
_NASA7_H2_HIGH = (
    2.93286575,
    8.26608026e-04,
    -1.46402364e-07,
    1.54100414e-11,
    -6.888048e-16,
    -813.065581,
)
_NASA7_CGR_LOW = (
    -0.310872072,
    4.40353686e-03,
    1.90394118e-06,
    -6.38546966e-09,
    2.98964248e-12,
    -108.650794,
)
_NASA7_CGR_HIGH = (
    1.45571829,
    1.71702216e-03,
    -6.97562786e-07,
    1.35277032e-10,
    -9.67590652e-15,
    -695.138814,
)


def _nasa7_enthalpy_j_mol(
    T_K: float, coeffs_low: tuple[float, ...], coeffs_high: tuple[float, ...]
) -> float:
    from pipeline.utils import R_J_MOL_K

    t = min(max(float(T_K), 200.0), 5000.0)
    a1, a2, a3, a4, a5, a6 = coeffs_low if t <= 1000.0 else coeffs_high
    return R_J_MOL_K * (
        a1 * t
        + 0.5 * a2 * (t ** 2)
        + (a3 / 3.0) * (t ** 3)
        + 0.25 * a4 * (t ** 4)
        + 0.2 * a5 * (t ** 5)
        + a6
    )


def pyrolysis_reaction_enthalpy_kJ_mol(T_K: float = 298.15) -> float:
    """Compute Kirchhoff reaction enthalpy ΔH_rxn(T) [kJ/mol CH4] for CH4 -> C(gr) + 2 H2.

    Evaluates the NASA-7 polynomials used in the Cantera mechanism for CH4(g),
    H2(g), and C(gr) at ``T_K``.

    Args:
        T_K: Absolute temperature in kelvin (default 298.15 K).

    Returns:
        Endothermic reaction enthalpy in kJ/mol CH4.
    """
    h_ch4 = _nasa7_enthalpy_j_mol(T_K, _NASA7_CH4_LOW, _NASA7_CH4_HIGH)
    h_h2 = _nasa7_enthalpy_j_mol(T_K, _NASA7_H2_LOW, _NASA7_H2_HIGH)
    h_cgr = _nasa7_enthalpy_j_mol(T_K, _NASA7_CGR_LOW, _NASA7_CGR_HIGH)
    return float((h_cgr + 2.0 * h_h2 - h_ch4) / 1000.0)


def load_gas_and_surface(
    config: _CarbonConfigLike,
) -> tuple[CanteraGasPhase, object | None, CanteraSurfacePhase | None]:
    """Load gas, optional graphite, and catalytic surface phases from Cantera YAML.

    Args:
        config: Reactor configuration specifying mechanism file and catalyst name.

    Returns:
        Tuple of ``(gas_phase, graphite_phase, surface_phase)``.
    """
    if ct is None:
        raise RuntimeError('Cantera is not installed; cannot load gas and surface')
    ct_mod: object = ct
    gas: CanteraGasPhase = cast(
        CanteraGasPhase, getattr(ct_mod, 'Solution')(config.mechanism_file, 'gas')
    )
    if 'C_graphite' in gas.species_names:
        raise RuntimeError(
            'mechanism still contains gas-phase C_graphite; regenerate YAML'
        )
    graphite: object | None = None
    try:
        graphite = getattr(ct_mod, 'Solution')(config.mechanism_file, 'graphite')
    except Exception as exc:
        if not config.gas_only:
            raise RuntimeError(
                f'graphite phase missing from {config.mechanism_file}; '
                'set gas_only=True only for an explicit gas-only run'
            ) from exc
    surf: CanteraSurfacePhase | None = None
    surf_name = f'{config.catalyst_name}_surface'
    try:
        adjacent: list[object] = [gas]
        if graphite is not None:
            adjacent.append(graphite)
        try:
            surf = cast(
                CanteraSurfacePhase,
                getattr(ct_mod, 'Interface')(
                    config.mechanism_file, surf_name, adjacent
                ),
            )
        except Exception:
            surf = cast(
                CanteraSurfacePhase,
                getattr(ct_mod, 'Interface')(
                    config.mechanism_file, surf_name, [gas]
                ),
            )
    except Exception as exc:
        if not config.gas_only:
            raise RuntimeError(
                f'surface {surf_name!r} missing from {config.mechanism_file} '
                f'(catalyst_name={config.catalyst_name!r} must match the '
                'YAML written by write_full_mechanism); set gas_only=True '
                'only for an explicit gas-only run'
            ) from exc
    return gas, graphite, surf


_load_gas_and_surface = load_gas_and_surface


def load_status_fields(
    config: _CarbonConfigLike, graphite: object | None, surf: object | None
) -> LoadStatusFields:
    """Summarize phase-loading status for reactor result provenance.

    Args:
        config: Reactor configuration for the evaluated run.
        graphite: Loaded graphite phase or ``None``.
        surf: Loaded surface phase or ``None``.

    Returns:
        Dictionary of phase-loading flags and surface name.
    """
    return {
        'gas_only': bool(config.gas_only),
        'surface_loaded': surf is not None,
        'graphite_loaded': graphite is not None,
        'surface_name': f'{config.catalyst_name}_surface',
    }


_load_status_fields = load_status_fields


def species_x(gas: CanteraGasPhase, name: str) -> float:
    """Return mole fraction of ``name`` in ``gas``, or 0.0 when absent.

    Args:
        gas: Cantera gas phase.
        name: Species identifier.

    Returns:
        Species mole fraction in ``[0, 1]``.
    """
    if name not in gas.species_names:
        return 0.0
    return float(gas.X[gas.species_index(name)])


_species_x = species_x


def surface_coverage(surf: Optional[CanteraSurfacePhase], name: str) -> float:
    """Return surface site coverage of ``name`` on ``surf``, or 0.0 when absent.

    Args:
        surf: Cantera surface phase or ``None``.
        name: Surface species identifier.

    Returns:
        Fractional surface site coverage in ``[0, 1]``.
    """
    if surf is None or name not in surf.species_names:
        return 0.0
    return float(surf.coverages[surf.species_index(name)])


_coverage = surface_coverage


def apply_continuous_carbon_removal(
    surf: Optional[CanteraSurfacePhase], rate_1_s: float, dt: float
) -> float:
    """Transport-style continuous removal of C_s → free sites.

    Args:
        surf: Cantera catalytic surface phase or ``None``.
        rate_1_s: First-order carbon removal rate constant in 1/s.
        dt: Integration time step in seconds.

    Returns:
        Approximate coverage fraction of C_s removed during ``dt``.
    """
    if surf is None or rate_1_s <= 0 or dt <= 0:
        return 0.0
    if 'C_s' not in surf.species_names or 'site' not in surf.species_names:
        return 0.0
    cov = np.array(surf.coverages, dtype=float)
    i_c = surf.species_index('C_s')
    i_site = surf.species_index('site')
    c_before = float(cov[i_c])
    removed = float(c_before * (1.0 - np.exp(-rate_1_s * dt)))
    cov[i_c] = c_before - removed
    cov[i_site] += removed
    site_mask = np.array([surf.species(n).size > 0 for n in range(surf.n_species)])
    s = float(cov[site_mask].sum())
    if s > 0:
        cov[site_mask] /= s
    surf.coverages = cov
    return removed


_apply_continuous_carbon_removal = apply_continuous_carbon_removal


def reset_surface_carbon(surf: Optional[CanteraSurfacePhase]) -> None:
    """Mechanical / consumable regen: clear adsorbates and restore free sites.

    Encapsulating carbon (``C_encap_s``, Baker-Helveg Cδ) is time-on-stream
    death and is never cleared.

    Args:
        surf: Cantera catalytic surface phase or ``None``.
    """
    if surf is None or 'C_s' not in surf.species_names:
        return
    encap = surface_coverage(surf, 'C_encap_s') if 'C_encap_s' in surf.species_names else 0.0
    cov = np.zeros(surf.n_species, dtype=float)
    if 'C_encap_s' in surf.species_names:
        cov[surf.species_index('C_encap_s')] = encap
    if 'site' in surf.species_names:
        cov[surf.species_index('site')] = max(0.0, 1.0 - encap)
    surf.coverages = cov


_reset_surface_carbon = reset_surface_carbon


def _surface_carbon_coverage(
    surf: CanteraSurfacePhase, cov: npt.NDArray[np.float64]
) -> float:
    """Total θ of carbon-bearing site-occupying species (C_s, CH_x_s, C_encap_s)."""
    total = 0.0
    for k in range(surf.n_species):
        sp = surf.species(k)
        if sp.size > 0 and sp.composition.get('C', 0) > 0:
            total += float(cov[k])
    return total


def tabulated_x_eq(T_K: float) -> float:
    """Return tabulated 1-bar methane pyrolysis equilibrium conversion at ``T_K``.

    Args:
        T_K: Temperature in kelvin.

    Returns:
        Equilibrium methane conversion fraction in ``[0, 1]``.
    """
    from pipeline.reactors.equilibrium import interpolate_x_eq

    return interpolate_x_eq(T_K)


_tabulated_x_eq = tabulated_x_eq


def local_tp_x_eq(T_K: float, P_Pa: float, y_ch4_0: float) -> float:
    """Return local (T, P) methane pyrolysis equilibrium conversion.

    Args:
        T_K: Temperature in kelvin.
        P_Pa: Pressure in pascals.
        y_ch4_0: Initial methane mole fraction.

    Returns:
        Equilibrium methane conversion fraction in ``[0, 1]``.
    """
    from pipeline.reactors.equilibrium import equilibrium_ch4_conversion_tp

    return equilibrium_ch4_conversion_tp(
        T_K, P_Pa, min(max(float(y_ch4_0), 1e-6), 1.0)
    )


_local_tp_x_eq = local_tp_x_eq


def compute_off_site_carbon_metrics(
    config: _CarbonConfigLike,
    gas: CanteraGasPhase,
    surf: Optional[CanteraSurfacePhase],
    *,
    feed_initial: float,
    ar_initial: float,
    pass_conversion: float,
    n_sites_mol: float,
    n_feed_fed_mol: float,
    n_parcel_in_mol: float,
    pass_time_s: float,
    cov_start: Optional[npt.NDArray[np.float64]],
    removed_surface_carbon_mol: float = 0.0,
    basis: str = 'Gamma*A_stage / (c_CH4*eps*V_stage), PFR parcel',
    mechanism_metadata_fn: Callable[[_CarbonConfigLike], Mapping[str, object]]
    | None = None,
) -> dict[str, object]:
    """Per-pass carbon accounting for the Alves-Xu / Baker-Helveg closure criterion.

    Args:
        config: Reactor configuration for the evaluated bed.
        gas: Cantera gas phase at bed exit.
        surf: Cantera catalytic surface phase or ``None``.
        feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        ar_initial: Inlet argon tracer mole fraction.
        pass_conversion: Single-pass feedstock conversion.
        n_sites_mol: Total catalytic surface sites in the stage/bed in moles.
        n_feed_fed_mol: Hydrocarbon feedstock fed per pass in moles.
        n_parcel_in_mol: Total gas moles in the inlet parcel.
        pass_time_s: Single-pass gas residence time in seconds.
        cov_start: Surface coverage vector at the start of the pass.
        removed_surface_carbon_mol: Surface carbon removed during the pass in moles.
        basis: Text label describing the site-to-gas normalization basis.
        mechanism_metadata_fn: Optional callback returning mechanism sidecar metadata.

    Returns:
        Dictionary of per-pass carbon balance and turnover metrics.
    """
    x_bound = n_sites_mol / n_feed_fed_mol if n_feed_fed_mol > 0 else None
    x_eq = tabulated_x_eq(config.T_inlet_K)
    x_eq_tp = local_tp_x_eq(config.T_inlet_K, config.P_inlet_Pa, feed_initial)
    out: dict[str, object] = {
        'site_inventory_bound_X': x_bound,
        'site_inventory_bound_basis': basis,
        'X_eq_table': x_eq,
        'X_eq_local_tp': x_eq_tp,
        'exceeds_equilibrium': bool(pass_conversion > x_eq + 1e-6),
    }
    if surf is None:
        return out
    cov_end = np.array(surf.coverages, dtype=float)
    if cov_start is None:
        cov_start = np.zeros_like(cov_end)
    x_ar = species_x(gas, 'Ar') or ar_initial
    n_ar = ar_initial * n_parcel_in_mol
    n_c2 = sum(
        2.0 * (species_x(gas, sp) or 0.0) / max(x_ar, 1e-30) * n_ar
        for sp in ('C2H2', 'C2H4', 'C2H6')
    )
    n_solid = max(0.0, pass_conversion * n_feed_fed_mol - n_c2)
    d_surface_c = (
        _surface_carbon_coverage(surf, cov_end)
        - _surface_carbon_coverage(surf, cov_start)
    ) * n_sites_mol
    n_removed_raw = max(0.0, float(removed_surface_carbon_mol))
    available_to_export = max(0.0, n_solid - d_surface_c)
    n_removed = min(n_removed_raw, available_to_export)
    balance_residual = n_solid - d_surface_c - n_removed_raw
    n_gamma = max(0.0, n_solid - d_surface_c - n_removed)
    has_encap = 'C_encap_s' in surf.species_names
    n_delta = 0.0
    if has_encap:
        i = surf.species_index('C_encap_s')
        n_delta = max(0.0, float(cov_end[i] - cov_start[i])) * n_sites_mol
    turnovers = n_solid / n_sites_mol if n_sites_mol > 0 else None
    ratio = (n_gamma / n_delta) if (has_encap and n_delta > 0) else None
    theta_encap = (
        float(cov_end[surf.species_index('C_encap_s')]) if has_encap else None
    )
    lifetime_h: float | None = None
    if has_encap and n_sites_mol > 0 and pass_time_s > 0:
        d_theta_dt = n_delta / n_sites_mol / pass_time_s
        lifetime_h = (
            ENCAP_DEAD_THETA / d_theta_dt / 3600.0 if d_theta_dt > 0 else None
        )
    out.update(
        {
            'carbon_turnovers_per_site': turnovers,
            'turnover_factor_vs_bound': (
                pass_conversion / x_bound if x_bound else None
            ),
            'solid_carbon_mol_per_pass': n_solid,
            'outfeed_carbon_mol_per_pass': n_removed,
            'outfeed_carbon_mol_raw': n_removed_raw,
            'carbon_balance_residual_mol': balance_residual,
            'carbon_balance_ok': bool(
                balance_residual >= -1e-3 * max(n_sites_mol, n_solid, 1e-12)
            ),
            'c_gamma_mol_per_pass': n_gamma if has_encap else None,
            'c_delta_mol_per_pass': n_delta if has_encap else None,
            'c_gamma_to_c_delta_ratio': ratio,
            'exit_theta_C_encap': theta_encap,
            'off_site_carbon_active': has_encap,
        }
    )
    if has_encap and theta_encap is not None:
        out.update(
            {
                'encapsulation_onset': bool(theta_encap >= ENCAP_ONSET_THETA),
                'encapsulation_onset_basis': f'exit theta_encap >= {ENCAP_ONSET_THETA}',
                'c_delta_competitive': (
                    None if ratio is None else bool(ratio < ENCAP_ONSET_RATIO)
                ),
                'c_delta_competitive_basis': f'C_gamma/C_delta < {ENCAP_ONSET_RATIO}',
                'encapsulation_lifetime_h': lifetime_h,
                'encapsulation_lifetime_basis': (
                    f'theta_encap -> {ENCAP_DEAD_THETA} at the pass-averaged '
                    'C_delta rate; linear extrapolation'
                ),
                'tos_lifetime_literature_band_h': list(NI_TOS_LIFETIME_BAND_H),
                'lifetime_within_tos_band': (
                    None
                    if lifetime_h is None
                    else bool(
                        NI_TOS_LIFETIME_BAND_H[0]
                        <= lifetime_h
                        <= NI_TOS_LIFETIME_BAND_H[1]
                    )
                ),
            }
        )
    if has_encap:
        metadata: Mapping[str, object] = (
            mechanism_metadata_fn(config) if mechanism_metadata_fn is not None else {}
        )
        raw_inputs = metadata.get('inputs')
        inputs_map: Mapping[str, object] = (
            cast(Mapping[str, object], raw_inputs)
            if isinstance(raw_inputs, Mapping)
            else {}
        )
        genome = inputs_map.get('genome')
        metal: str | None = None
        try:
            import pipeline.reactors.mechanisms as _mech_mod

            parsed_g = _mech_mod.parse_catalyst_genome(genome)
            raw_m_cls = inputs_map.get('material_class')
            nanoparticle_fn = getattr(_mech_mod, '_nanoparticle_metals')
            metals: set[str] = (
                nanoparticle_fn(
                    parsed_g,
                    str(raw_m_cls) if raw_m_cls is not None else None,
                )
                & _mech_mod.OFF_SITE_CARBON_METALS
            )
            if len(metals) == 1:
                metal = next(iter(metals))
            elif len(metals) > 1:
                if (
                    parsed_g
                    and len(parsed_g) > 1
                    and str(parsed_g[1]) in _mech_mod.OFF_SITE_CARBON_METALS
                ):
                    metal = str(parsed_g[1])
                else:
                    metal = sorted(metals)[0]
        except Exception:
            metal = None
        dispersion = float(config.metal_dispersion)
        yield_rate: float | None = None
        if metal and dispersion > 0 and pass_time_s > 0 and n_sites_mol > 0:
            n_metal = n_sites_mol / dispersion
            yield_rate = (
                (n_gamma / pass_time_s)
                * 3600.0
                * 12.011
                / (n_metal * METAL_MOLAR_MASS_G_MOL[metal])
            )
        eps_parcel = (
            0.45
            if getattr(config, 'reactor_type', '') == 'Fluidized'
            else float(getattr(config, 'bed_porosity', 0.4))
        )
        rho_cat = float(getattr(config, 'catalyst_density_kg_m3', 2500.0))
        bulk_metal_mass_g_per_m3_bed = (
            (1.0 - eps_parcel) * rho_cat * float(config.metal_loading) * 1000.0
        )
        c_total_in = float(config.P_inlet_Pa) / (8.314462618 * float(config.T_inlet_K))
        v_bed_parcel_m3 = (
            n_parcel_in_mol / max(c_total_in * eps_parcel, 1e-12)
        )
        bulk_metal_mass_parcel_g = bulk_metal_mass_g_per_m3_bed * v_bed_parcel_m3
        bulk_yield_rate = (
            (n_gamma / pass_time_s) * 3600.0 * 12.011 / bulk_metal_mass_parcel_g
            if (metal and pass_time_s > 0 and bulk_metal_mass_parcel_g > 0)
            else None
        )
        out.update(
            {
                'filament_metal': metal,
                'filament_yield_gC_per_gMetal_h': yield_rate,
                'filament_yield_basis': (
                    'pass-averaged C_gamma rate; metal moles = sites / dispersion'
                ),
                'bulk_metal_mass_g_per_m3_bed': bulk_metal_mass_g_per_m3_bed,
                'bulk_filament_yield_gC_per_gMetal_bulk_h': bulk_yield_rate,
                'filament_yield_literature_band_gC_per_gNi_h': list(
                    NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H
                ),
                'filament_yield_within_band': (
                    None
                    if yield_rate is None
                    else bool(
                        NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[0]
                        <= yield_rate
                        <= NI_FILAMENT_YIELD_BAND_G_C_PER_G_NI_H[1]
                    )
                ),
            }
        )
    return out


_off_site_carbon_metrics = compute_off_site_carbon_metrics


def h2_atom_balance_metric(
    feed_initial: float,
    final_conv: float,
    x_h2: float,
    *,
    x_ar_feed: Optional[float] = None,
    x_ar_out: Optional[float] = None,
) -> float:
    """Expansion-corrected H-atom balance ratio (2 * n_H2_out) / (4 * n_feed_converted).

    Accounts for the molar gas expansion (1 + y_feed,0 * X, or exact Ar-tracer
    scale y_Ar,0 / y_Ar,out when provided) accompanying the baseline Cantera
    CH4 -> C(s) + 2 H2 pass (where CH4 is the default hydrocarbon feedstock;
    non-methane feedstocks are subsequently mapped by `adapt_reactor_result_to_feedstock`).

    Args:
        feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        final_conv: Final single-pass conversion fraction.
        x_h2: Outlet H2 mole fraction.
        x_ar_feed: Optional inlet argon mole fraction for exact tracer scaling.
        x_ar_out: Optional outlet argon mole fraction for exact tracer scaling.

    Returns:
        Clipped H-atom balance ratio in ``[0, 1]``.
    """
    if final_conv <= 0.01:
        return 0.0
    expansion: float
    if (
        x_ar_feed is not None
        and x_ar_out is not None
        and x_ar_feed > 0.0
        and x_ar_out > 0.0
    ):
        expansion = float(x_ar_feed) / float(x_ar_out)
    else:
        expansion = 1.0 + max(0.0, float(feed_initial)) * max(0.0, float(final_conv))
    h_in_feed_converted: float = 4.0 * float(feed_initial) * float(final_conv)
    h_in_h2_out: float = 2.0 * float(x_h2) * expansion
    return float(np.clip(h_in_h2_out / max(h_in_feed_converted, 1e-10), 0.0, 1.0))


_h2_atom_balance_metric = h2_atom_balance_metric


def compute_h2_c2_selectivity_and_rates(
    gas: Optional[CanteraGasPhase],
    *,
    feed_initial: float = 0.95,
    ar_initial: float = 0.05,
    final_conv: float = 0.0,
    tau_s: float = 1.25,
    P_inlet_Pa: float = 101325.0,
    T_inlet_K: float = 1000.0,
    active_sv_1_m: Optional[float] = None,
    site_density_mol_cm2: float = 2.5e-9,
    bed_porosity_or_holdup: float = 0.4,
    x_h2_override: Optional[float] = None,
    x_feed_override: Optional[float] = None,
    x_c2h6_override: Optional[float] = None,
    x_c2h4_override: Optional[float] = None,
    x_c2h2_override: Optional[float] = None,
    x_ar_override: Optional[float] = None,
    feed_species: str = 'CH4',
) -> H2C2SelectivityRateMetrics:
    """Compute true H2 vs. C2 branching selectivity, H2 yield, and molar H2 production rates.

    Uses the Ar internal standard (``expansion = y_Ar,0 / y_Ar,out``) when ``gas``
    carries ``C2H6``, ``C2H4``, ``C2H2``, ``H2``, and ``Ar`` mole fractions, and
    computes volumetric and site-normalized H2 production rates (TOF_H2).
    Note: The base Cantera gas phase represents the hydrocarbon reactant using
    species ``'CH4'`` (the default feedstock). When a non-methane feedstock is
    configured on ``ReactorConfig.feedstock``, ``adapt_reactor_result_to_feedstock``
    updates both the ``feedstock_*`` and ``CH4_*`` fields for the selected
    hydrocarbon or biomass feedstock.

    Args:
        gas: Outlet Cantera gas phase or ``None``.
        feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        ar_initial: Inlet argon mole fraction.
        final_conv: Single-pass hydrocarbon conversion fraction in ``[0, 1]``.
        tau_s: Reactor gas residence time in seconds.
        P_inlet_Pa: Inlet pressure in pascals.
        T_inlet_K: Inlet temperature in kelvin.
        active_sv_1_m: Active catalytic surface-to-volume ratio in 1/m.
        site_density_mol_cm2: Monolayer site density in mol/cm².
        bed_porosity_or_holdup: Gas void fraction or bubble holdup in ``(0, 1)``.
        x_h2_override: Optional outlet H2 mole fraction when ``gas`` is ``None``.
        x_feed_override: Optional outlet unreacted feedstock mole fraction when ``gas`` is ``None``.
        x_c2h6_override: Optional outlet C2H6 mole fraction when ``gas`` is ``None``.
        x_c2h4_override: Optional outlet C2H4 mole fraction when ``gas`` is ``None``.
        x_c2h2_override: Optional outlet C2H2 mole fraction when ``gas`` is ``None``.
        x_ar_override: Optional outlet Ar mole fraction when ``gas`` is ``None``.
        feed_species: Primary hydrocarbon species name in ``gas`` (default ``'CH4'``).

    Returns:
        Dictionary of selectivity, yield, C2 mole fractions, and H2 rate metrics.
    """
    conv: float = float(np.clip(final_conv, 0.0, 1.0))
    y_feed_0: float = max(float(feed_initial), 1e-12)
    y_ar_0: float = max(float(ar_initial), 0.0)
    x_h2: float
    x_feed_out: float
    x_ar: float
    x_c2h6: float
    x_c2h4: float
    x_c2h2: float
    if gas is not None:
        x_h2 = float(species_x(gas, 'H2') or 0.0)
        x_feed_out = float(species_x(gas, feed_species) or 0.0)
        x_ar = float(species_x(gas, 'Ar') or y_ar_0)
        x_c2h6 = (
            float(species_x(gas, 'C2H6') or 0.0)
            if feed_species != 'C2H6'
            else 0.0
        )
        x_c2h4 = float(species_x(gas, 'C2H4') or 0.0)
        x_c2h2 = float(species_x(gas, 'C2H2') or 0.0)
        if (x_c2h6 + x_c2h4 + x_c2h2) <= 1e-12 and conv > 1e-6 and float(T_inlet_K) >= 850.0:
            eps_void: float = min(max(float(bed_porosity_or_holdup), 0.01), 0.99)
            tau_v: float = max(float(tau_s), 1e-3)
            k_het_eff_s: float = -float(np.log(max(1.0 - min(conv, 0.995), 1e-4))) / tau_v
            k_hom_c2_s: float = 4.5e13 * float(
                np.exp(-356000.0 / (8.314462618 * max(float(T_inlet_K), 300.0)))
            )
            s_c2_void: float = float(
                np.clip(
                    0.08 * (eps_void * k_hom_c2_s) / (k_het_eff_s + eps_void * k_hom_c2_s),
                    0.0,
                    0.08,
                )
            )
            exp_cur: float = (
                y_ar_0 / x_ar
                if (y_ar_0 > 0.0 and x_ar > 1e-12)
                else (1.0 + y_feed_0 * conv)
            )
            n_c2_tot: float = 0.5 * y_feed_0 * conv * s_c2_void
            x_c2h6 = (0.15 * n_c2_tot) / exp_cur
            x_c2h4 = (0.75 * n_c2_tot) / exp_cur
            x_c2h2 = (0.10 * n_c2_tot) / exp_cur
    else:
        expansion_nom: float = 1.0 + y_feed_0 * conv
        x_h2 = (
            float(x_h2_override)
            if x_h2_override is not None
            else (2.0 * y_feed_0 * conv) / expansion_nom
        )
        x_feed_out = (
            float(x_feed_override)
            if x_feed_override is not None
            else (y_feed_0 * (1.0 - conv)) / expansion_nom
        )
        x_ar = (
            float(x_ar_override)
            if x_ar_override is not None
            else y_ar_0 / expansion_nom
        )
        x_c2h6 = float(x_c2h6_override) if x_c2h6_override is not None else 0.0
        x_c2h4 = float(x_c2h4_override) if x_c2h4_override is not None else 0.0
        x_c2h2 = float(x_c2h2_override) if x_c2h2_override is not None else 0.0

    expansion: float = (
        y_ar_0 / x_ar
        if (y_ar_0 > 0.0 and x_ar > 1e-12)
        else (1.0 + y_feed_0 * conv)
    )
    dn_feed_converted: float = max(y_feed_0 * conv, 1e-12)
    n_c_in_c2: float = 2.0 * (x_c2h6 + x_c2h4 + x_c2h2) * expansion
    n_h_in_c2: float = (6.0 * x_c2h6 + 4.0 * x_c2h4 + 2.0 * x_c2h2) * expansion

    c2_c_sel: float
    solid_c_sel: float
    h2_sel: float
    h2_bal: float
    if conv <= 1e-6:
        c2_c_sel = 0.0
        solid_c_sel = 1.0
        h2_sel = 1.0
        h2_bal = 0.0
    else:
        c2_c_sel = float(np.clip(n_c_in_c2 / dn_feed_converted, 0.0, 1.0))
        solid_c_sel = float(np.clip(1.0 - c2_c_sel, 0.0, 1.0))
        h_frac_c2: float = float(np.clip(n_h_in_c2 / (4.0 * dn_feed_converted), 0.0, 1.0))
        h2_sel = float(np.clip(1.0 - h_frac_c2, 0.0, 1.0))
        h2_bal = h2_atom_balance_metric(
            y_feed_0, conv, x_h2, x_ar_feed=y_ar_0, x_ar_out=x_ar
        )

    h2_yield: float = float(conv * h2_sel)

    # Volumetric and site-normalized feedstock conversion & H2 production rates
    c_total_mol_m3: float = float(P_inlet_Pa) / (8.314462618 * max(float(T_inlet_K), 1.0))
    eps_eff: float = min(max(float(bed_porosity_or_holdup), 0.01), 0.99)
    tau_eff: float = max(float(tau_s), 1e-6)
    # Molar feedstock conversion rate and H2 production rate per m³ reactor volume [mol/(m³·s)]
    r_feed_mol_m3_s: float = float(eps_eff * c_total_mol_m3 * y_feed_0 * conv / tau_eff)
    r_h2_mol_m3_s: float = float(2.0 * r_feed_mol_m3_s * h2_sel)
    r_h2_kg_m3_h: float = float(r_h2_mol_m3_s * 0.00201588 * 3600.0)

    tof_h2_1_s: Optional[float] = None
    if active_sv_1_m is not None and float(active_sv_1_m) > 0.0:
        gamma_mol_m2: float = float(site_density_mol_cm2) * 1e4
        sites_mol_m3: float = gamma_mol_m2 * float(active_sv_1_m)
        if sites_mol_m3 > 0.0:
            tof_h2_1_s = float(r_h2_mol_m3_s / sites_mol_m3)

    return {
        'H2_selectivity': h2_sel,
        'H2_atom_balance': h2_bal,
        'H2_yield': h2_yield,
        'solid_C_selectivity': solid_c_sel,
        'C2_carbon_selectivity': c2_c_sel,
        'exit_x_feedstock': x_feed_out,
        'exit_x_CH4': x_feed_out,
        'exit_x_C2H6': x_c2h6,
        'exit_x_C2H4': x_c2h4,
        'exit_x_C2H2': x_c2h2,
        'feedstock_conversion_rate_mol_m3_s': r_feed_mol_m3_s,
        'ch4_conversion_rate_mol_m3_s': r_feed_mol_m3_s,
        'h2_production_rate_mol_m3_s': r_h2_mol_m3_s,
        'h2_production_rate_kg_m3_h': r_h2_kg_m3_h,
        'h2_turnover_frequency_1_s': tof_h2_1_s,
    }


def primary_feed_species(gas: CanteraGasPhase) -> str:
    """Identify the primary inlet hydrocarbon reactant species in a Cantera gas phase.

    Defaults strictly to `'CH4'` when methane is present, and falls back to
    `'C3H8'`, `'C2H6'`, `'C3H6'`, or `'C2H4'` when a multi-carbon feed is loaded.

    Args:
        gas: Cantera gas phase initialized at inlet composition.

    Returns:
        Cantera species name of the primary hydrocarbon feed.
    """
    if species_x(gas, 'CH4') > 1e-6:
        return 'CH4'
    for candidate_sp in ('C3H8', 'C2H6', 'C3H6', 'C2H4'):
        if species_x(gas, candidate_sp) > 1e-6:
            return candidate_sp
    return 'CH4'


def feedstock_extent(
    gas: CanteraGasPhase,
    x_feed_initial: float,
    x_ar_feed: float,
    feed_species: str = 'CH4',
) -> float:
    """Compute argon-tracer-normalized hydrocarbon feedstock conversion extent.

    Evaluates conversion of ``feed_species`` (default ``'CH4'``, or ``'C2H6'`` /
    ``'C3H8'`` / ``'C3H6'`` when an extended multi-carbon mechanism is active).

    Args:
        gas: Cantera gas phase.
        x_feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        x_ar_feed: Inlet argon mole fraction.
        feed_species: Primary hydrocarbon species name in ``gas`` (default ``'CH4'``).

    Returns:
        Feedstock conversion fraction in ``[0, 1]``.
    """
    from pipeline.reactors.equilibrium import ch4_conversion_from_argon_tracer

    return ch4_conversion_from_argon_tracer(
        species_x(gas, feed_species), species_x(gas, 'Ar'), x_feed_initial, x_ar_feed
    )


def clamp_gas_to_equilibrium(
    gas: CanteraGasPhase,
    x_feed_initial: float,
    x_ar_feed: float,
    feed_species: str = 'CH4',
) -> None:
    """Prevent irreversible C_s => C(gr) + site sink from pulling X past X_eq(T, P).

    Args:
        gas: Cantera gas phase modified in place.
        x_feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        x_ar_feed: Inlet argon mole fraction.
        feed_species: Primary hydrocarbon species name in ``gas`` (default ``'CH4'``).
    """
    if x_feed_initial <= 0.0 or x_ar_feed <= 0.0:
        return
    try:
        if feed_species == 'CH4':
            x_eq_1bar = tabulated_x_eq(float(gas.T))
            p_pa = float(gas.P)
            if p_pa > 105000.0:
                x_eq = min(
                    x_eq_1bar,
                    local_tp_x_eq(float(gas.T), p_pa, x_feed_initial),
                )
            else:
                x_eq = x_eq_1bar
        else:
            from pipeline.reactors.feedstocks import feedstock_equilibrium_conversion

            x_eq = feedstock_equilibrium_conversion(
                float(gas.T),
                float(gas.P),
                feed_species,
                y_feed_inlet=x_feed_initial,
            )
    except ValueError:
        return
    conv = feedstock_extent(gas, x_feed_initial, x_ar_feed, feed_species=feed_species)
    if conv <= x_eq:
        return
    x_ar = species_x(gas, 'Ar')
    if x_ar <= 0.0 or feed_species not in gas.species_names:
        return
    scale_to_feed = float(x_ar_feed) / x_ar
    n_vec = np.array(gas.X, dtype=float) * scale_to_feed
    i_feed = gas.species_index(feed_species)
    i_h2 = gas.species_index('H2')
    target_n_feed = float(x_feed_initial) * (1.0 - x_eq)
    delta_feed = max(0.0, target_n_feed - float(n_vec[i_feed]))
    n_vec[i_feed] = target_n_feed
    n_vec[i_h2] = max(0.0, float(n_vec[i_h2]) - 2.0 * delta_feed)
    total = float(n_vec.sum())
    if total > 0.0:
        gas.TPX = float(gas.T), float(gas.P), n_vec / total


_clamp_gas_to_equilibrium = clamp_gas_to_equilibrium


def require_ar_tracer(x_feed_initial: float, x_ar_feed: float) -> None:
    """Validate that positive hydrocarbon feedstock and Ar mole fractions are present at the inlet.

    Args:
        x_feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        x_ar_feed: Inlet argon mole fraction.
    """
    if x_ar_feed <= 0 or x_feed_initial <= 0:
        raise ValueError(
            'solids CH4 conversion uses an Ar mole-fraction tracer; '
            'inlet_composition must include Ar (default CH4:0.95, Ar:0.05)'
        )


_require_ar_tracer = require_ar_tracer


def _disable_reactor_energy(reactor: object) -> str:
    """Hold solids stages at inlet T so they match isothermal MMBCR."""
    if hasattr(reactor, 'energy_enabled'):
        setattr(reactor, 'energy_enabled', False)
        return 'isothermal_energy_disabled'
    raise RuntimeError(
        'Cantera reactor has no energy_enabled; cannot force isothermal solids'
    )


def advance_reactor_stage(
    gas: CanteraGasPhase,
    surf: Optional[CanteraSurfacePhase],
    volume_m3: float,
    surface_area_m2: float,
    dt: float,
) -> None:
    """Advance an isothermal gas+surface stage, falling back to substepped tight tolerances on stiff CVODES failure.

    Args:
        gas: Cantera gas phase modified in place.
        surf: Optional Cantera surface phase modified in place.
        volume_m3: Stage gas volume in m³.
        surface_area_m2: Stage catalytic surface area in m².
        dt: Stage residence time in seconds.
    """
    if ct is None:
        raise RuntimeError('Cantera is not installed; cannot advance reactor stage')
    ct_mod: object = ct
    gas_t, gas_p = float(gas.T), float(gas.P)
    gas_x = np.array(gas.X, dtype=float)
    surf_cov = np.array(surf.coverages, dtype=float) if surf is not None else None

    def _step_stage(step_dt: float, *, tight: bool = False) -> None:
        ideal_gas_reactor_cls = getattr(ct_mod, 'IdealGasReactor')
        reactor_surface_cls = getattr(ct_mod, 'ReactorSurface')
        reactor_net_cls = getattr(ct_mod, 'ReactorNet')
        try:
            reactor = ideal_gas_reactor_cls(gas, clone=False)
        except TypeError:
            reactor = ideal_gas_reactor_cls(gas)
        reactor.volume = volume_m3
        _disable_reactor_energy(reactor)
        if surf is not None:
            try:
                reactor_surface_cls(surf, reactor, A=surface_area_m2, clone=False)
            except TypeError:
                reactor_surface_cls(surf, reactor, A=surface_area_m2)
        net = reactor_net_cls([reactor])
        if tight:
            net.rtol = 1e-11
            net.atol = 1e-18
            net.max_time_step = max(step_dt / 10.0, 1e-7)
        net.advance(step_dt)
        phase = getattr(reactor, 'phase', None) or reactor.thermo
        gas.TPX = phase.T, phase.P, phase.X

    try:
        _step_stage(dt, tight=False)
        return
    except Exception:
        pass

    gas.TPX = gas_t, gas_p, gas_x
    if surf is not None and surf_cov is not None:
        surf.coverages = surf_cov
    n_sub = 25
    sub_dt = dt / n_sub
    for _ in range(n_sub):
        _step_stage(sub_dt, tight=True)


_advance_reactor_stage = advance_reactor_stage


def mole_fraction_drop(
    gas: CanteraGasPhase, x_feed_initial: float, feed_species: str = 'CH4'
) -> float:
    """Compute uncalibrated mole-fraction drop ``1 - y_feed / y_feed,0``.

    Args:
        gas: Cantera gas phase.
        x_feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        feed_species: Primary hydrocarbon species name in ``gas`` (default ``'CH4'``).

    Returns:
        Fractional mole-fraction drop in ``[0, 1]``.
    """
    if x_feed_initial <= 0:
        return 0.0
    return float(max(0.0, 1.0 - species_x(gas, feed_species) / x_feed_initial))


_mole_fraction_drop = mole_fraction_drop


def set_gas_from_feedstock_conversion(
    gas: CanteraGasPhase,
    T_K: float,
    P_Pa: float,
    x_feed_initial: float,
    x_ar_feed: float,
    conversion: float,
) -> None:
    """Set baseline Cantera gas composition from hydrocarbon → C(s) + 2 H2 conversion where solid C exits by flotation.

    Args:
        gas: Cantera gas phase modified in place.
        T_K: Temperature in kelvin.
        P_Pa: Pressure in pascals.
        x_feed_initial: Inlet hydrocarbon feedstock (default methane) mole fraction.
        x_ar_feed: Inlet argon mole fraction.
        conversion: Feedstock conversion fraction in ``[0, 1]``.
    """
    x_clamped = float(np.clip(conversion, 0.0, 1.0))
    n_feed = x_feed_initial * (1.0 - x_clamped)
    n_h2 = 2.0 * x_feed_initial * x_clamped
    n_ar = x_ar_feed
    n_tot = n_feed + n_h2 + n_ar
    if n_tot <= 0:
        return
    gas.TPX = (
        T_K,
        P_Pa,
        f'CH4:{n_feed / n_tot}, H2:{n_h2 / n_tot}, Ar:{n_ar / n_tot}',
    )


def mix_bubble_bypass(
    gas: CanteraGasPhase,
    inlet_x: npt.NDArray[np.float64],
    delta: float,
    x_ar_feed: float,
) -> None:
    """Mix a bypass stream (fraction δ of feed, unreacted) with the emulsion outlet on a molar basis.

    Args:
        gas: Emulsion outlet Cantera gas phase modified in place.
        inlet_x: Inlet gas mole-fraction array.
        delta: Effective bubble bypass fraction in ``[0, 1)``.
        x_ar_feed: Inlet argon mole fraction.
    """
    x_em = np.array(gas.X, dtype=float)
    x_ar_em = species_x(gas, 'Ar')
    if x_ar_em <= 0:
        raise RuntimeError(
            'Ar tracer vanished from emulsion outlet; cannot mix bypass'
        )
    expansion = x_ar_feed / x_ar_em
    moles = (
        delta * np.asarray(inlet_x, dtype=float)
        + (1.0 - delta) * expansion * x_em
    )
    total = float(moles.sum())
    if total <= 0:
        raise RuntimeError('bubble bypass mixing produced no moles')
    gas.TPX = gas.T, gas.P, moles / total


_mix_bubble_bypass = mix_bubble_bypass

