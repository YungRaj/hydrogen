#!/usr/bin/env python3
"""
Catalyst-independent Phase 2 diagnostic: CH4 conversion vs E_act.

Holds reactor geometry, T, P, and mechanism structure fixed; only the
lumped CH4 activation barrier changes. A discriminating reactor should
show a sigmoid: near-equilibrium at low E_act, collapse at high E_act.
A flat curve far below equilibrium means residence time / area / site
density / prefactor dominate — Phase 2 cannot rank catalysts.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from pipeline.common.utils import REACTOR_DIR, setup_logger, save_json
from pipeline.process.reactor_mechanisms import write_full_mechanism
from pipeline.process.reactor_models import (
    DIAGNOSTIC_MATERIAL_CLASS, SINGLE_REACTOR_MODE, ReactorConfig,
    simulate_reactor)
from pipeline.process.equilibrium_check import TABULATED_X_CH4_1BAR

logger = setup_logger('eact_sensitivity', 'reactor/eact_sensitivity.log')


def run_eact_sweep(
    temperatures: Optional[List[float]] = None,
    e_acts_eV: Optional[List[float]] = None,
    reactor_types: Optional[List[str]] = None,
    reactor_config_kwargs: Optional[Dict] = None,
    catalyst_stub: str = 'eact_sweep',
    material_class: Optional[str] = None,
    kinetics=None,
) -> Dict:
    """Sweep E_act on a template mechanism (or a fixed ``kinetics`` record with
    only the CH4 barrier varied) across reactors.

    ``material_class`` applies to every reactor; None uses the per-reactor
    diagnostic class (SolidCatalyst for beds, MoltenMetal for MMBCR). A
    real candidate class makes incompatible reactors ``not_applicable``.

    Args:
        temperatures: Input controlling temperatures.
        e_acts_eV: Input controlling e acts eV.
        reactor_types: Input controlling reactor types.
        reactor_config_kwargs: Input controlling reactor config kwargs.
        catalyst_stub: Input controlling catalyst stub.

    Returns:
        Validated Dict output for this operation.
    """
    if temperatures is None:
        temperatures = [1300.0]
    if e_acts_eV is None:
        e_acts_eV = list(np.round(np.linspace(0.1, 2.0, 20), 3))
    if reactor_types is None:
        reactor_types = ['MMBCR', 'PFR', 'Fluidized']
    extra = dict(reactor_config_kwargs or {})
    # Default: no discrete regen; continuous removal at current production default.
    extra.setdefault('co2_permitted', False)
    extra.setdefault('max_regen_cycles', 3)
    extra.setdefault('fluidized_mode', 'circulating')
    extra.setdefault('regen_mechanism', 'mechanical')

    rows = []
    for E in e_acts_eV:
        cat = f'{catalyst_stub}_{E:.3f}'.replace('.', 'p')
        if kinetics is not None:
            from dataclasses import replace as _replace
            varied = _replace(kinetics, methane_activation_eV=float(E))
            mech = write_full_mechanism(cat, kinetics=varied, T_ref=temperatures[0])
        else:
            mech = write_full_mechanism(cat, E_act_CH4=float(E), T_ref=temperatures[0])
        for rt in reactor_types:
            rt_class = material_class or DIAGNOSTIC_MATERIAL_CLASS[rt]
            for T in temperatures:
                cfg = ReactorConfig(
                    T_inlet_K=float(T),
                    reactor_type=rt,
                    pathway_mode=SINGLE_REACTOR_MODE[rt],
                    material_class=rt_class,
                    mechanism_file=str(mech),
                    catalyst_name=cat,
                    candidate_id=cat,
                    catalyst_E_act_eV=float(E),
                    **extra,
                )
                result = simulate_reactor(cfg)
                x_eq = TABULATED_X_CH4_1BAR.get(float(T))
                if x_eq is None:
                    nearest = min(TABULATED_X_CH4_1BAR, key=lambda t: abs(t - T))
                    x_eq = TABULATED_X_CH4_1BAR[nearest]
                row = {
                    'E_act_eV': float(E),
                    'T_K': float(T),
                    'reactor_type': rt,
                    'material_class': rt_class,
                    'status': result.get('status'),
                    'reason': result.get('reason'),
                    'CH4_conversion': float(result.get('CH4_conversion', 0.0) or 0.0),
                    'X_eq_table': float(x_eq),
                    'exit_theta_C': result.get('exit_theta_C'),
                    'max_theta_C': result.get('max_theta_C'),
                    'carbon_removed_coverage_proxy': result.get(
                        'carbon_removed_coverage_proxy'),
                    'co2_permitted': result.get('co2_permitted'),
                    'mmbcr_carbon_removal_rate_1_s': extra.get(
                        'mmbcr_carbon_removal_rate_1_s'),
                }
                rows.append(row)
                logger.info(
                    f"E={E:.3f} eV  {rt}  T={T}  X={row['CH4_conversion']:.4f}  "
                    f"X_eq≈{x_eq:.3f}"
                )

    # Summarize discriminating power per reactor at each T
    summaries = []
    for rt in reactor_types:
        for T in temperatures:
            sub = [r for r in rows if r['reactor_type'] == rt and r['T_K'] == T]
            xs = [r['CH4_conversion'] for r in sub]
            es = [r['E_act_eV'] for r in sub]
            x_eq = sub[0]['X_eq_table'] if sub else None
            x_low = xs[0] if xs else None   # lowest E_act
            x_high = xs[-1] if xs else None
            span = (max(xs) - min(xs)) if xs else 0.0
            # Monotonic decreasing in E_act?
            mono = all(xs[i] >= xs[i + 1] - 1e-9 for i in range(len(xs) - 1)) if len(xs) > 1 else True
            frac_of_eq = (x_low / x_eq) if (x_eq and x_low is not None and x_eq > 0) else None
            discriminating = bool(
                span > 0.05 and frac_of_eq is not None and frac_of_eq > 0.5 and mono
            )
            summaries.append({
                'reactor_type': rt,
                'T_K': T,
                'X_at_Emin': x_low,
                'X_at_Emax': x_high,
                'span': span,
                'X_eq_table': x_eq,
                'frac_of_eq_at_Emin': frac_of_eq,
                'monotonic_in_E_act': mono,
                'discriminating': discriminating,
                'verdict': (
                    'discriminating_sigmoid_like' if discriminating else
                    ('flat_or_weak' if span <= 0.05 else
                     'responds_but_far_below_equilibrium')
                ),
            })

    out = {
        'diagnostic': 'eact_sensitivity',
        'points': rows,
        'summaries': summaries,
        'reactor_config_kwargs': extra,
        'any_discriminating': any(s['discriminating'] for s in summaries),
    }
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    save_json(out, 'eact_sensitivity_sweep.json', subdir='reactor')
    return out


# PFR clears C_s only by discrete regen. It does not read a continuous
# detachment rate; a flat PFR pair is "not wired", not "detachment = 0".
ABLATION_RESPONSE = {
    'MMBCR': True,
    'Fluidized': True,
    'PFR': False,
}


def run_detachment_ablation(T_K: float = 1300.0, E_act_eV: float = 0.1) -> Dict:
    """
    Ablate continuous carbon takeoff where that knob is actually wired.

    MMBCR: ``mmbcr_carbon_removal_rate_1_s`` is flotation frequency.
    0 fouls the interface; finite keeps it open.
    Fluidized circulating: ``circulating_carbon_removal_rate_1_s`` is
    applied during integrate (substeps), not after ``net.advance()``.
    PFR: discrete regen only. Rows are kept with
    ``responds_to_ablated_variable=false``.

    Args:
        T_K: Input controlling T K.
        E_act_eV: Input controlling E act eV.

    Returns:
        Validated Dict output for this operation.
    """
    cat = 'ablation_cat'
    mech = write_full_mechanism(cat, E_act_CH4=E_act_eV, T_ref=T_K)
    cases = {
        'blocked_no_detach': {'mmbcr_carbon_removal_rate_1_s': 0.0,
                              'circulating_carbon_removal_rate_1_s': 0.0},
        'detach_transport_lump': {'mmbcr_carbon_removal_rate_1_s': 1e3,
                                  'circulating_carbon_removal_rate_1_s': 1e3},
    }
    rows = []
    for name, kwargs in cases.items():
        for rt in ('MMBCR', 'PFR', 'Fluidized'):
            cfg = ReactorConfig(
                T_inlet_K=T_K,
                reactor_type=rt,
                pathway_mode=SINGLE_REACTOR_MODE[rt],
                material_class=DIAGNOSTIC_MATERIAL_CLASS[rt],
                mechanism_file=str(mech),
                catalyst_name=cat,
                candidate_id=cat,
                catalyst_E_act_eV=E_act_eV,
                co2_permitted=False,
                max_regen_cycles=0,
                fluidized_mode='circulating',
                regen_mechanism='mechanical',
                **kwargs,
            )
            result = simulate_reactor(cfg)
            rows.append({
                'case': name,
                'reactor_type': rt,
                'E_act_eV': E_act_eV,
                'T_K': T_K,
                'CH4_conversion': float(result.get('CH4_conversion', 0) or 0),
                'exit_theta_C': result.get('exit_theta_C'),
                'max_theta_C': result.get('max_theta_C'),
                'responds_to_ablated_variable': ABLATION_RESPONSE[rt],
                'ablation_note': (
                    None if ABLATION_RESPONSE[rt] else
                    'PFR clears C_s only by discrete regen; '
                    'max_regen_cycles=0 so both cases are the same produce pass'
                ),
            })
            logger.info(
                f"ablation {name} {rt}: X={rows[-1]['CH4_conversion']:.4f} "
                f"responds={ABLATION_RESPONSE[rt]}"
            )
    by_rt = {}
    for row in rows:
        by_rt.setdefault(row['reactor_type'], {})[row['case']] = row
    interpretations = {}
    for rt, cases in by_rt.items():
        blocked = cases.get('blocked_no_detach', {}).get('CH4_conversion', 0.0)
        detach = cases.get('detach_transport_lump', {}).get('CH4_conversion', 0.0)
        delta = float(detach) - float(blocked)
        if not ABLATION_RESPONSE[rt]:
            meaning = 'not_wired'
        elif abs(delta) > 1e-4:
            meaning = 'wired_and_limiting'
        else:
            meaning = 'wired_not_limiting_at_this_point'
        interpretations[rt] = {
            'delta_X': delta,
            'interpretation': meaning,
        }
        for row in rows:
            if row['reactor_type'] == rt:
                row['ablation_delta_X'] = delta
                row['ablation_interpretation'] = meaning
    out = {
        'diagnostic': 'detachment_ablation',
        'points': rows,
        'wired_reactors': [rt for rt, ok in ABLATION_RESPONSE.items() if ok],
        'unwired_reactors': [rt for rt, ok in ABLATION_RESPONSE.items() if not ok],
        'interpretation': interpretations,
    }
    save_json(out, 'detachment_ablation.json', subdir='reactor')
    return out


if __name__ == '__main__':
    from pipeline.common.utils import print_banner
    print_banner('E_ACT SENSITIVITY SWEEP')
    sweep = run_eact_sweep()
    print(json.dumps(sweep['summaries'], indent=2))
    print_banner('DETACHMENT ABLATION')
    abl = run_detachment_ablation()
    print(json.dumps(abl['points'], indent=2))
