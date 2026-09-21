"""Phase 2 solids scorecard.

Reports single-pass X, active a, WHSV, and Ergun ΔP. H-parked 0.01 eV
cats and MMBCR X_eq are not ranks. A named judge catalyst is a campaign
argument, not a module constant.
"""

from __future__ import annotations

from pipeline.process.result_eligibility import is_usable_result

SOLIDS_TYPES = frozenset({'PFR', 'Fluidized'})
H_PARKED_E_ACT_MAX = 0.05
H_PARKED_ABS_DEH_MIN = 2.0
# High-T end of the standard 4-point pyrolysis sweep (773, 900, 1100, 1300).
# 1300 K is the ADR 0001 band ceiling, not the B5/B6-7 Ni headline. The Ni
# judge uses the 650–700 °C filament ROI (headline_t_min/max = 923.15/973.15).
DEFAULT_HEADLINE_T_MIN = 1200.0
# Alves 2021 / Xu 2021: Ni filaments win at 650–700 °C. Inclusive window.
NI_JUDGE_HEADLINE_T_MIN = 923.15
NI_JUDGE_HEADLINE_T_MAX = 973.15
NI_JUDGE_CATALYST = 'ni_np_lit'


def is_h_parked(record: dict) -> bool:
    e_act = record.get('catalyst_E_act_eV')
    dE_H = record.get('catalyst_dE_H_eV')
    if e_act is None or dE_H is None:
        return False
    return float(e_act) <= H_PARKED_E_ACT_MAX and abs(float(dE_H)) >= H_PARKED_ABS_DEH_MIN


def is_production_reactor_record(record: dict) -> bool:
    """A completed single-condition run. Upstream's non-run records
    (``not_applicable``, ``validation_required``, ``failed``) carry a
    ``status`` and no conversion; they are evidence about the workflow, not
    the catalyst, and never enter the scorecard. Legacy JSON has no
    ``status`` and is accepted on the conversion key."""
    return (
        isinstance(record, dict)
        and record.get('reactor_type') in {'PFR', 'Fluidized', 'MMBCR'}
        and bool(record.get('catalyst_name'))
        and record.get('status', 'complete') == 'complete'
        and 'CH4_conversion' in record
        and 'records' not in record
    )


def is_solids_run(record: dict) -> bool:
    """PFR/fluidized with a loaded Langmuir surface. Missing surface_loaded
    (legacy JSON) or gas_only/mock is not a solids result."""
    return (
        is_production_reactor_record(record)
        and record.get('reactor_type') in SOLIDS_TYPES
        and record.get('surface_loaded') is True
        and not record.get('gas_only', False)
        and not record.get('mock', False)
    )


def is_scoreable_solids(record: dict) -> bool:
    """Solids run that may enter a headline or max-X rank.

    Shared usable baseline (complete, not mock, no X>X_eq, carbon
    balance when reported) plus the solids surface filters. Overshoot
    and carbon-fail rows stay in ``n_solids_records``; they never rank.
    """
    return is_solids_run(record) and is_usable_result(record)


def single_pass_x(record: dict) -> float:
    value = record.get('single_pass_CH4_conversion', record.get('CH4_conversion'))
    return float(value or 0.0)


def _t_k(record: dict) -> float:
    try:
        return float(record.get('T_K') or 0.0)
    except (TypeError, ValueError):
        return 0.0


def in_headline_band(record: dict, t_min: float,
                     t_max: float | None = None) -> bool:
    """True if T is in [t_min, t_max]. ``t_max is None`` means no upper bound."""
    temperature = _t_k(record)
    if temperature < float(t_min) - 1e-9:
        return False
    if t_max is not None and temperature > float(t_max) + 1e-6:
        return False
    return True


def metric_row(record: dict) -> dict:
    return {
        'catalyst_name': record.get('catalyst_name'),
        'reactor_type': record.get('reactor_type'),
        'T_K': record.get('T_K'),
        'single_pass_CH4_conversion': single_pass_x(record),
        'active_sv_1_m': record.get('active_sv_1_m'),
        'WHSV_h-1': record.get('WHSV_h-1'),
        'ergun_delta_p_Pa': record.get('ergun_delta_p_Pa'),
        'ergun_delta_p_bar': record.get('ergun_delta_p_bar'),
        'ergun_ok': record.get('ergun_ok'),
        'catalyst_E_act_eV': record.get('catalyst_E_act_eV'),
        'catalyst_dE_H_eV': record.get('catalyst_dE_H_eV'),
        'h_parked': is_h_parked(record),
        'catalyst_particle_mm': record.get('catalyst_particle_mm'),
        'metal_loading': record.get('metal_loading'),
        'metal_dispersion': record.get('metal_dispersion'),
        'exceeds_equilibrium': bool(record.get('exceeds_equilibrium')),
        'carbon_balance_ok': record.get('carbon_balance_ok'),
    }


def build_solids_scorecard(results, *,
                           judge_catalyst: str | None = None,
                           headline_t_min: float = DEFAULT_HEADLINE_T_MIN,
                           headline_t_max: float | None = None) -> dict:
    solids_src = [r for r in results if is_solids_run(r)]
    solids = [metric_row(r) for r in solids_src]
    mmbcr = [
        r for r in results
        if is_production_reactor_record(r)
        and r.get('reactor_type') == 'MMBCR'
        and is_usable_result(r)
    ]
    eligible = [metric_row(r) for r in solids_src
                if is_scoreable_solids(r) and not is_h_parked(r)]
    judge_src = (
        [r for r in solids_src if r.get('catalyst_name') == judge_catalyst]
        if judge_catalyst else []
    )
    named_present = bool(judge_src)
    if named_present:
        pool = [metric_row(r) for r in judge_src if is_scoreable_solids(r)]
        rank_key = _t_k
        judge_reason = (
            'named judge catalyst; not 0.01 eV H-parked; '
            'MMBCR X_eq is not a solids rank')
    else:
        pool = eligible
        rank_key = lambda r: r['single_pass_CH4_conversion']
        if judge_catalyst:
            judge_reason = (
                f'{judge_catalyst} missing from solids results; '
                'headline is best non-H-parked; MMBCR X_eq is not a rank')
        else:
            judge_reason = (
                'no named judge; headline is best non-H-parked; '
                'MMBCR X_eq is not a solids rank')

    headline = {}
    for reactor_type in ('PFR', 'Fluidized'):
        candidates = [
            r for r in pool
            if r['reactor_type'] == reactor_type
            and in_headline_band(r, headline_t_min, headline_t_max)
        ]
        if candidates:
            headline[reactor_type] = max(candidates, key=rank_key)

    solids_max = (
        max(eligible, key=lambda r: r['single_pass_CH4_conversion'])
        if eligible else None
    )
    h_parked = [r for r in solids if r['h_parked'] and _t_k(r) >= headline_t_min]
    mmbcr_max = max(
        (float(r.get('CH4_conversion') or 0.0) for r in mmbcr),
        default=None,
    )
    judge_x = None
    if headline.get('PFR'):
        judge_x = headline['PFR']['single_pass_CH4_conversion']
    elif headline.get('Fluidized'):
        judge_x = headline['Fluidized']['single_pass_CH4_conversion']

    headline_catalyst = None
    if headline.get('PFR'):
        headline_catalyst = headline['PFR']['catalyst_name']
    elif headline.get('Fluidized'):
        headline_catalyst = headline['Fluidized']['catalyst_name']

    return {
        'judge_catalyst': judge_catalyst if named_present else None,
        'headline_catalyst': headline_catalyst,
        'judge_catalyst_requested': judge_catalyst,
        'judge_reason': judge_reason,
        'headline_t_min': float(headline_t_min),
        'headline_t_max': None if headline_t_max is None else float(headline_t_max),
        'headline_band_note': (
            '1300 K is the ADR 0001 ceiling (773–1300 K), used as the '
            'legacy 4-point headline when t_min=1200 and t_max is open. '
            'The Ni judge headline is the 650–700 °C filament ROI '
            f'({NI_JUDGE_HEADLINE_T_MIN:g}–{NI_JUDGE_HEADLINE_T_MAX:g} K); '
            'X>X_eq and carbon-balance failures are never the headline.'),
        'headline': headline,
        'headline_solids_conversion': judge_x,
        'solids_max_excluding_h_parked': solids_max,
        'h_parked_excluded': h_parked,
        'mmbcr_max_conversion': mmbcr_max,
        'mmbcr_note': 'X_eq by construction for large Da; not a catalyst rank (B5)',
        'n_solids_records': len(solids),
        'n_mmbcr_records': len(mmbcr),
    }


def log_solids_scorecard(scorecard: dict, logger) -> None:
    logger.info(
        f"B1-5 solids scorecard: judge={scorecard.get('judge_catalyst')} "
        f"({scorecard.get('judge_reason')})"
    )
    for reactor_type, row in (scorecard.get('headline') or {}).items():
        x = row.get('single_pass_CH4_conversion') or 0.0
        a = row.get('active_sv_1_m')
        whsv = row.get('WHSV_h-1')
        dp = row.get('ergun_delta_p_bar')
        logger.info(
            f"  {reactor_type} {row.get('catalyst_name')} "
            f"T={row.get('T_K')}  X={x:.2%}  a={a}  "
            f"WHSV={whsv} h-1  dP={dp} bar"
        )
    mmbcr = scorecard.get('mmbcr_max_conversion')
    if mmbcr is not None:
        logger.info(f"  MMBCR max X={mmbcr:.2%} ({scorecard.get('mmbcr_note')})")
