#!/usr/bin/env python3
"""Catalyst-independent CH4 ⇌ C(gr) + 2 H2 equilibrium check vs tabulated conversion."""

from __future__ import annotations

from typing import Dict, List, Optional

from pipeline.common.utils import MECHANISMS_DIR, setup_logger
from pipeline.process.reactor_mechanisms import write_gas_only_mechanism

logger = setup_logger('equilibrium_check', 'reactor/equilibrium_check.log')

# From TURQUOISE_HYDROGEN.md — X_CH4 at 1 bar for CH4 ⇌ C(s) + 2 H2
TABULATED_X_CH4_1BAR = {
    773.15: 0.385,
    900.00: 0.684,
    1000.00: 0.842,
    1100.00: 0.927,
    1200.00: 0.967,
    1300.00: 0.985,
    1400.00: 0.993,
}


def ch4_conversion_from_mole_fractions(x_ch4: float, x_h2: float) -> float:
    """Extent of CH4 → C(s) + 2 H2 from a CH4/H2/inert gas.

    Exact only when H2 comes from that stoichiometry. The solids mechanism
    has an active C2 chain; do not use this for PFR/fluidized. Those use
    ``ch4_conversion_from_argon_tracer``.
    """
    denom = float(x_ch4) + 0.5 * float(x_h2)
    if denom <= 0:
        return 0.0
    return float(min(1.0, max(0.0, 1.0 - float(x_ch4) / denom)))


def ch4_conversion_from_argon_tracer(x_ch4: float, x_ar: float,
                                     x_ch4_feed: float, x_ar_feed: float) -> float:
    """CH4 conversion from an inert Ar tracer. Exact for any gas-C fate.

    X = 1 - (x_CH4/x_Ar) / (x_CH4,0/x_Ar,0). Carbon to C(s) or C2 does not
    change Ar moles, so mole expansion and the C2 chain cancel.
    """
    if x_ar_feed <= 0 or x_ar <= 0 or x_ch4_feed <= 0:
        raise ValueError(
            'Ar tracer conversion requires Ar in the feed and the current gas; '
            'solids reactors must keep an inert mole-fraction tracer')
    ratio0 = float(x_ch4_feed) / float(x_ar_feed)
    ratio = float(x_ch4) / float(x_ar)
    return float(min(1.0, max(0.0, 1.0 - ratio / ratio0)))


def _ch4_conversion_from_mix(gas, ch4_initial_moles: float, graphite_moles: float) -> float:
    """Equilibrium conversion from remaining gas-phase carbon vs initial CH4."""
    if ch4_initial_moles <= 0:
        return 0.0
    # Prefer graphite moles formed when available; fall back to CH4 depletion.
    if graphite_moles > 0:
        return float(min(1.0, max(0.0, graphite_moles / ch4_initial_moles)))
    if 'CH4' not in gas.species_names:
        return 0.0
    x_ch4 = float(gas.X[gas.species_index('CH4')])
    x_h2 = float(gas.X[gas.species_index('H2')]) if 'H2' in gas.species_names else 0.0
    return ch4_conversion_from_mole_fractions(x_ch4, x_h2)


def run_equilibrium_sweep(
    mechanism_file: Optional[str] = None,
    temperatures: Optional[List[float]] = None,
    pressure_Pa: float = 101325.0,
    tabulated: Optional[Dict[float, float]] = None,
    abs_tol: float = 0.08,
) -> Dict:
    """
    Equilibrate gas + condensed graphite at fixed T,P and compare to tables.

    Uses Cantera Mixture (Cantera 3.0) of ideal-gas + fixed-stoichiometry graphite.
    """
    try:
        import cantera as ct
    except ImportError as exc:
        raise RuntimeError('Cantera is required for equilibrium_check') from exc

    if mechanism_file is None:
        mech_path = write_gas_only_mechanism()
        mechanism_file = str(mech_path)
    if temperatures is None:
        temperatures = sorted(TABULATED_X_CH4_1BAR.keys())
    if tabulated is None:
        tabulated = TABULATED_X_CH4_1BAR

    gas = ct.Solution(mechanism_file, 'gas')
    graphite = ct.Solution(mechanism_file, 'graphite')

    rows = []
    worst_abs_err = 0.0
    for T in temperatures:
        gas.TPX = T, pressure_Pa, 'CH4:1.0'
        graphite.TP = T, pressure_Pa
        # Pure CH4 feed: 1 mol gas, 0 mol graphite initially.
        mix = ct.Mixture([(gas, 1.0), (graphite, 0.0)])
        mix.T = T
        mix.P = pressure_Pa
        mix.equilibrate('TP')

        gr_moles = float(mix.phase_moles(1))
        x_calc = _ch4_conversion_from_mix(gas, 1.0, gr_moles)
        # Prefer graphite inventory as conversion for C(s) product.
        x_from_graphite = float(min(1.0, max(0.0, gr_moles)))
        x_ch4_eq = x_from_graphite if gr_moles > 1e-12 else x_calc

        x_tab = tabulated.get(float(T))
        if x_tab is None:
            # nearest tabulated key
            nearest = min(tabulated.keys(), key=lambda t: abs(t - T))
            if abs(nearest - T) < 1.0:
                x_tab = tabulated[nearest]
        abs_err = abs(x_ch4_eq - x_tab) if x_tab is not None else None
        if abs_err is not None:
            worst_abs_err = max(worst_abs_err, abs_err)

        rows.append({
            'T_K': float(T),
            'P_Pa': float(pressure_Pa),
            'X_CH4_calc': x_ch4_eq,
            'X_CH4_table': x_tab,
            'abs_error': abs_err,
            'graphite_moles': gr_moles,
            'gas_X': {sp: float(gas.X[gas.species_index(sp)]) for sp in gas.species_names},
        })
        logger.info(
            f"  T={T:.2f} K  X_calc={x_ch4_eq:.3f}  X_table={x_tab}  "
            f"|err|={abs_err if abs_err is not None else 'n/a'}"
        )

    within_tol = all(
        r['abs_error'] is None or r['abs_error'] <= abs_tol for r in rows
    )
    return {
        'mechanism_file': mechanism_file,
        'abs_tol': abs_tol,
        'worst_abs_error': worst_abs_err,
        'within_tolerance': within_tol,
        'points': rows,
    }


if __name__ == '__main__':
    from pipeline.common.utils import print_banner
    print_banner('EQUILIBRIUM CHECK: CH4 = C(gr) + 2 H2')
    result = run_equilibrium_sweep()
    print(f"within_tolerance={result['within_tolerance']} "
          f"worst_abs_error={result['worst_abs_error']:.4f}")
    for row in result['points']:
        print(f"  {row['T_K']:7.2f} K  calc={row['X_CH4_calc']:.3f}  "
              f"table={row['X_CH4_table']}  err={row['abs_error']}")
