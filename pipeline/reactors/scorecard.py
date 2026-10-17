"""Phase 2 solids scorecard.

Reports single-pass X, active a, WHSV, and Ergun ΔP. H-parked 0.01 eV
cats and MMBCR X_eq are not ranks. A named judge catalyst is a campaign
argument, not a module constant.
"""

from __future__ import annotations

from typing import Iterable, Mapping, Protocol, cast

from pipeline.reactors.eligibility import finite_number, is_usable_result
from pipeline.reactors.modes import ReactorTypeName
from pipeline.data_models.reactors import (
    B5ClosureGateResult,
    SolidsMetricRow,
    SolidsScorecard,
)


class InfoLogger(Protocol):
    """Minimal logging sink accepted by log_solids_scorecard."""

    def info(self, msg: str) -> object:
        """Write an informational message.

        Args:
            msg: Formatted message string to write to the log sink.

        Returns:
            Optional logger return value.
        """
        ...


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
# B5 closure: solids conversion has no authority to rank catalysts until the
# supported-Ni reference passes at 650–700 °C with off-site carbon on and the
# B6-6 A_gamma / theta* sweeps are accepted. Flip only with that evidence.
B5_ACCEPTED = False


def is_h_parked(record: Mapping[str, object]) -> bool:
    """Return whether excessive hydrogen coverage invalidates a solids rank.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    e_act = record.get('catalyst_E_act_eV')
    dE_H = record.get('catalyst_dE_H_eV')
    if not isinstance(e_act, (int, float)) or not isinstance(dE_H, (int, float)):
        return False
    return (
        float(e_act) <= H_PARKED_E_ACT_MAX and abs(float(dE_H)) >= H_PARKED_ABS_DEH_MIN
    )


def is_production_reactor_record(record: Mapping[str, object]) -> bool:
    """A completed single-condition run. Upstream's non-run records
    (``not_applicable``, ``validation_required``, ``failed``) carry a
    ``status`` and no conversion; they are evidence about the workflow, not
    the catalyst, and never enter the scorecard. Legacy JSON has no
    ``status`` and is accepted on the conversion key.

    Returns:
        Validated bool output for this operation.
    """
    return (
        isinstance(record, Mapping)  # pyright: ignore[reportUnnecessaryIsInstance]
        and record.get('reactor_type')
        in {'PFR', 'Fluidized', 'MMBCR', 'NTEC', 'Electrochemical'}
        and bool(record.get('catalyst_name'))
        and record.get('status', 'complete') == 'complete'
        and 'CH4_conversion' in record
        and 'records' not in record
    )


def is_solids_run(record: Mapping[str, object]) -> bool:
    """PFR/fluidized with a loaded Langmuir surface. Missing surface_loaded
    (legacy JSON) or gas_only/mock is not a solids result.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    return (
        is_production_reactor_record(record)
        and record.get('reactor_type') in SOLIDS_TYPES
        and record.get('surface_loaded') is True
        and not record.get('gas_only', False)
        and not record.get('mock', False)
    )


def is_scoreable_solids(record: Mapping[str, object]) -> bool:
    """Solids run that may enter a headline or max-X rank.

    Shared usable baseline (complete, not mock, no X>X_eq, carbon
    balance when reported) plus the solids surface filters. Overshoot
    and carbon-fail rows stay in ``n_solids_records``; they never rank.

    Returns:
        Validated bool output for this operation.
    """
    return is_solids_run(record) and is_usable_result(record)


def single_pass_x(record: Mapping[str, object]) -> float:
    """Read the canonical single-pass methane conversion from a record.

    Args:
        record: Input controlling record.

    Returns:
        Validated float output for this operation.
    """
    value = record.get('single_pass_CH4_conversion', record.get('CH4_conversion'))
    return float(value) if isinstance(value, (int, float)) else 0.0


def _t_k(record: Mapping[str, object]) -> float:
    value = record.get('T_K')
    if isinstance(value, (int, float)):
        return float(value)
    try:
        return float(str(value or 0.0))
    except (TypeError, ValueError):
        return 0.0


def in_headline_band(
    record: Mapping[str, object], t_min: float, t_max: float | None = None
) -> bool:
    """True if T is in [t_min, t_max]. ``t_max is None`` means no upper bound.

    Args:
        record: Input controlling record.

    Returns:
        Validated bool output for this operation.
    """
    temperature = _t_k(record)
    if temperature < float(t_min) - 1e-9:
        return False
    if t_max is not None and temperature > float(t_max) + 1e-6:
        return False
    return True


def _opt_str(value: object) -> str | None:
    return str(value) if isinstance(value, str) else None


def _opt_bool(value: object) -> bool | None:
    return value if isinstance(value, bool) else None


def metric_row(record: Mapping[str, object]) -> SolidsMetricRow:
    """Project a reactor result into the solids-scorecard metric schema.

    Args:
        record: Input controlling record.

    Returns:
        Validated dict output for this operation.
    """
    rt = _opt_str(record.get('reactor_type'))
    return {
        'catalyst_name': _opt_str(record.get('catalyst_name')),
        'reactor_type': cast(ReactorTypeName | None, rt),
        'T_K': finite_number(record.get('T_K')),
        'single_pass_CH4_conversion': single_pass_x(record),
        'active_sv_1_m': finite_number(record.get('active_sv_1_m')),
        'WHSV_h-1': finite_number(record.get('WHSV_h-1')),
        'ergun_delta_p_Pa': finite_number(record.get('ergun_delta_p_Pa')),
        'ergun_delta_p_bar': finite_number(record.get('ergun_delta_p_bar')),
        'ergun_ok': _opt_bool(record.get('ergun_ok')),
        'catalyst_E_act_eV': finite_number(record.get('catalyst_E_act_eV')),
        'catalyst_dE_H_eV': finite_number(record.get('catalyst_dE_H_eV')),
        'h_parked': is_h_parked(record),
        'catalyst_particle_mm': finite_number(record.get('catalyst_particle_mm')),
        'metal_loading': finite_number(record.get('metal_loading')),
        'metal_dispersion': finite_number(record.get('metal_dispersion')),
        'exceeds_equilibrium': bool(record.get('exceeds_equilibrium')),
        'carbon_balance_ok': _opt_bool(record.get('carbon_balance_ok')),
    }


def build_solids_scorecard(
    results: Iterable[Mapping[str, object]],
    *,
    judge_catalyst: str | None = None,
    headline_t_min: float = DEFAULT_HEADLINE_T_MIN,
    headline_t_max: float | None = None,
) -> SolidsScorecard:
    """Build fail-closed PFR, fluidized, and MMBCR comparison summaries.

    Args:
        results: Input controlling results.
        judge_catalyst: Input controlling judge catalyst.
        headline_t_min: Input controlling headline t min.
        headline_t_max: Input controlling headline t max.

    Returns:
        Validated dict output for this operation.
    """
    solids_src = [r for r in results if is_solids_run(r)]
    solids = [metric_row(r) for r in solids_src]
    mmbcr = [
        r
        for r in results
        if is_production_reactor_record(r)
        and r.get('reactor_type') == 'MMBCR'
        and is_usable_result(r)
    ]
    ntec = [
        r
        for r in results
        if is_production_reactor_record(r)
        and r.get('reactor_type') == 'NTEC'
        and is_usable_result(r)
    ]
    electrochemical = [
        r
        for r in results
        if is_production_reactor_record(r)
        and r.get('reactor_type') == 'Electrochemical'
        and is_usable_result(r)
    ]
    judge_src = (
        [r for r in solids_src if r.get('catalyst_name') == judge_catalyst]
        if judge_catalyst
        else []
    )
    named_present = bool(judge_src)
    if named_present:
        pool: list[SolidsMetricRow] = [
            metric_row(r) for r in judge_src if is_scoreable_solids(r)
        ]
        judge_reason = (
            'named judge catalyst; not 0.01 eV H-parked; '
            'MMBCR X_eq is not a solids rank'
        )
    else:
        # B5/B6 remain open: Phase 2 conversion describes a specified judge,
        # but is not authorized to choose among catalysts.
        pool = []
        if judge_catalyst:
            judge_reason = (
                f'{judge_catalyst} missing from solids results; '
                'no catalyst fallback while B5/B6 are open'
            )
        else:
            judge_reason = (
                'no named judge; Phase 2 conversion has no ranking authority '
                'while B5/B6 are open'
            )

    headline: dict[str, SolidsMetricRow] = {}
    for reactor_type in ('PFR', 'Fluidized'):
        candidates = [
            r
            for r in pool
            if r.get('reactor_type') == reactor_type
            and in_headline_band(r, headline_t_min, headline_t_max)
        ]
        if candidates:
            headline[reactor_type] = max(candidates, key=_t_k)

    h_parked = [
        r for r in solids if r.get('h_parked') and _t_k(r) >= headline_t_min
    ]
    mmbcr_max = max(
        (finite_number(r.get('CH4_conversion')) or 0.0 for r in mmbcr),
        default=None,
    )
    ntec_max = max(
        (finite_number(r.get('CH4_conversion')) or 0.0 for r in ntec),
        default=None,
    )
    electrochemical_max = max(
        (finite_number(r.get('CH4_conversion')) or 0.0 for r in electrochemical),
        default=None,
    )
    judge_x: float | None = None
    if 'PFR' in headline:
        judge_x = headline['PFR'].get('single_pass_CH4_conversion')
    elif 'Fluidized' in headline:
        judge_x = headline['Fluidized'].get('single_pass_CH4_conversion')

    headline_catalyst: str | None = None
    if 'PFR' in headline:
        headline_catalyst = headline['PFR'].get('catalyst_name')
    elif 'Fluidized' in headline:
        headline_catalyst = headline['Fluidized'].get('catalyst_name')

    scoreable_non_h_parked = [
        r for r in solids_src if is_scoreable_solids(r) and not is_h_parked(r)
    ]
    candidate_solids = [
        r for r in scoreable_non_h_parked
        if r.get('catalyst_name') != judge_catalyst
    ]
    best_candidate_solids = max(
        candidate_solids,
        key=lambda r: (single_pass_x(r), -_t_k(r)),
        default=None,
    )
    solids_max_ex_h = (
        max((single_pass_x(r) for r in scoreable_non_h_parked), default=None)
        if B5_ACCEPTED
        else None
    )

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
            'X>X_eq and carbon-balance failures are never the headline.'
        ),
        'headline': headline,
        'headline_solids_conversion': judge_x,
        'ranking_authority': B5_ACCEPTED,
        'solids_max_excluding_h_parked': solids_max_ex_h,
        'candidate_champion_conversion': (
            float(single_pass_x(best_candidate_solids))
            if best_candidate_solids is not None
            else None
        ),
        'candidate_champion_catalyst': (
            str(best_candidate_solids.get('catalyst_name'))
            if best_candidate_solids is not None
            else None
        ),
        'candidate_champion_reactor': (
            str(best_candidate_solids.get('reactor_type'))
            if best_candidate_solids is not None
            else None
        ),
        'candidate_champion_T_K': (
            float(_t_k(best_candidate_solids))
            if best_candidate_solids is not None
            else None
        ),
        'candidate_champion_kinetics_status': (
            str(best_candidate_solids.get('kinetics_status') or '') or None
            if best_candidate_solids is not None
            else None
        ),
        'h_parked_excluded': h_parked,
        'mmbcr_max_conversion': mmbcr_max,
        'mmbcr_note': 'X_eq by construction for large Da; not a catalyst rank (B5)',
        'ntec_max_conversion': ntec_max,
        'electrochemical_max_conversion': electrochemical_max,
        'n_solids_records': len(solids),
        'n_mmbcr_records': len(mmbcr),
        'n_ntec_records': len(ntec),
        'n_electrochemical_records': len(electrochemical),
    }


def log_solids_scorecard(scorecard: SolidsScorecard, logger: InfoLogger) -> None:
    """Write the principal scorecard decisions to a supplied logger.

    Args:
        scorecard: Input controlling scorecard.
        logger: Input controlling logger.
    """
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
    cand_x = scorecard.get('candidate_champion_conversion')
    if cand_x is not None:
        logger.info(
            f"  Candidate solids champion: {scorecard.get('candidate_champion_catalyst')} "
            f"({scorecard.get('candidate_champion_reactor')} @ "
            f"{scorecard.get('candidate_champion_T_K')} K) X={cand_x:.2%}"
        )
    mmbcr = scorecard.get('mmbcr_max_conversion')
    if mmbcr is not None:
        logger.info(f"  MMBCR max X={mmbcr:.2%} ({scorecard.get('mmbcr_note')})")
    ntec = scorecard.get('ntec_max_conversion')
    if ntec is not None:
        logger.info(f"  NTEC max X={ntec:.2%} (coupled MOOSE+Cantera screening)")
    echem = scorecard.get('electrochemical_max_conversion')
    if echem is not None:
        logger.info(f"  Electrochemical max X={echem:.2%} (coupled MOOSE+Cantera screening)")


def evaluate_b5_closure_gate(
    scorecard: SolidsScorecard,
    b66_summary: Mapping[str, object] | None = None,
    joint_summary: Mapping[str, object] | None = None,
    *,
    min_roi_conversion: float = 0.50,
) -> B5ClosureGateResult:
    """Evaluate whether B5 solids closure criteria are met by reactor evidence.

    Args:
        scorecard: Solids scorecard dictionary containing headline results.
        b66_summary: Optional B6-6 carbon-criteria summary dictionary.
        joint_summary: Optional B6-7 joint-band summary dictionary.
        min_roi_conversion: Minimum single-pass CH4 conversion required in the ROI.

    Returns:
        Dictionary of B5 gate pass/fail flags and ROI conversion values.
    """
    headline = scorecard.get('headline') or {}
    pfr = headline.get('PFR') or {}
    fluid = headline.get('Fluidized') or {}
    pfr_x = float(pfr.get('single_pass_CH4_conversion') or 0.0)
    fluid_x = float(fluid.get('single_pass_CH4_conversion') or 0.0)
    judge_ok = scorecard.get('judge_catalyst') == NI_JUDGE_CATALYST
    pfr_ok = judge_ok and pfr_x >= min_roi_conversion and not pfr.get('exceeds_equilibrium', False)
    fluid_ok = judge_ok and fluid_x >= min_roi_conversion and not fluid.get('exceeds_equilibrium', False)
    raw_criteria = (b66_summary or {}).get('criteria')
    empty_map: Mapping[str, object] = {}
    criteria: Mapping[str, object] = (
        cast(Mapping[str, object], raw_criteria)
        if isinstance(raw_criteria, Mapping)
        else empty_map
    )
    b66_ok = bool(criteria)
    for key in ('flat', 'linearity', 'encapsulation_onset', 'lifetime', 'yield'):
        entry = criteria.get(key)
        entry_map: Mapping[str, object] = (
            cast(Mapping[str, object], entry)
            if isinstance(entry, Mapping)
            else empty_map
        )
        if not bool(entry_map.get('pass')):
            b66_ok = False
    joint_ok = (
        isinstance(joint_summary, Mapping)
        and joint_summary.get('declaration') == 'simultaneous_hit'
        and int(str(joint_summary.get('n_both_roi') or 0)) > 0
    )
    passed = pfr_ok and fluid_ok and b66_ok and joint_ok
    return {
        'b5_accepted': passed,
        'judge_ok': judge_ok,
        'pfr_roi_conversion': pfr_x,
        'fluidized_roi_conversion': fluid_x,
        'pfr_conversion_ok': pfr_ok,
        'fluidized_conversion_ok': fluid_ok,
        'b66_criteria_ok': b66_ok,
        'b67_joint_roi_ok': joint_ok,
    }



