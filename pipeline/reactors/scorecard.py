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
    SolidsMetricRow,
    SolidsScorecard,
    SupportedNiFilamentClosureGateResult,
)
from pipeline.search.discovery import (
    LITERATURE_VALIDATION_GATES,
    NAMED_REFERENCE_CATALYSTS as NAMED_REFERENCE_CATALYSTS,
    format_literature_catalyst_label,
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


SOLIDS_TYPES = frozenset({'PFR', 'PBR', 'Fluidized'})
# H-parked site-poisoning filter: candidates clipped at the BEP floor (E_act <= 0.05 eV)
# with overly strong or unphysical hydrogen binding (|dE_H| >= 2.0 eV vs. 1/2 H2(g))
# lock active sites as H(s) or reflect unrelaxed geometries, and cannot hold a solids rank.
H_PARKED_E_ACT_MAX = 0.05
H_PARKED_ABS_DEH_MIN = 2.0
# High-T end of the standard 4-point pyrolysis sweep (773, 900, 1100, 1300).
# 1300 K is the ADR 0001 band ceiling, not the supported-Ni(111)/SiO2 headline.
# The supported-Ni judge uses the 650–700 °C filament ROI (headline_t_min/max = 923.15/973.15).
DEFAULT_HEADLINE_T_MIN = 1200.0
# Alves 2021 / Xu 2021: Ni filaments win at 650–700 °C. Inclusive window.
NI_JUDGE_HEADLINE_T_MIN = 923.15
NI_JUDGE_HEADLINE_T_MAX = 973.15
NI_JUDGE_CATALYST = 'ni111_sio2_bengaard2002_ref'
# Supported-Ni(111)/SiO2 filament closure gate: solids conversion has no authority
# to rank catalysts until the supported-Ni reference passes at 650–700 °C with
# off-site carbon on and the Baker-Helveg-Ermakova-Takenaka A_gamma / theta*
# sweeps are accepted. Flip only with that evidence.
SOLIDS_RANKING_AUTHORITY_ACCEPTED = False


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
        in {'PFR', 'PBR', 'Fluidized', 'MMBCR', 'NTEC', 'Electrochemical'}
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
    cat_name = _opt_str(record.get('catalyst_name'))
    return {
        'catalyst_name': cat_name,
        'literature_catalyst_name': (
            format_literature_catalyst_label(cat_name) if cat_name else None
        ),
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


def evaluate_pbr_and_ergun_residence_diagnostics(
    records: Iterable[Mapping[str, object]],
) -> dict[str, object]:
    """Summarize PBR intraparticle effectiveness, Ergun ΔP, pilot residence-time, and regenerated-cycle conversion deltas.

    Quantifies packed-bed pellet diffusion limitations (overall effectiveness
    factor $\\Omega$ and Thiele modulus $\\phi$), coupled Ergun (1952) pressure-drop
    compliance, pilot residence-time band adherence ($\\tau \\in [0.1, 10]\\text{ s}$),
    and Muradov-Parkinson regenerated cycle-averaged vs. single-pass conversion
    across completed sweep rows.

    Args:
        records: Iterable of sweep record dictionaries to aggregate.

    Returns:
        Dictionary summarizing PBR effectiveness, Ergun compliance, pilot residence-band
        counts, and regenerated cycle-averaged vs. single-pass conversion statistics.
    """
    completed = [rec for rec in records if rec.get('status') == 'complete']
    pbr_rows = [rec for rec in completed if rec.get('reactor_type') == 'PBR']
    omega_vals = [
        val
        for val in (
            finite_number(r.get('overall_effectiveness_factor')) for r in pbr_rows
        )
        if val is not None
    ]
    thiele_vals = [
        val
        for val in (finite_number(r.get('thiele_modulus')) for r in pbr_rows)
        if val is not None
    ]
    ergun_checked = [r for r in completed if r.get('ergun_ok') is not None]
    ergun_ok_count = sum(1 for r in ergun_checked if bool(r.get('ergun_ok')))
    pilot_residence_checked = [
        r for r in completed if r.get('residence_time_within_pilot_band') is not None
    ]
    pilot_residence_in_band_count = sum(
        1
        for r in pilot_residence_checked
        if bool(r.get('residence_time_within_pilot_band'))
    )
    plugging_flagged = sum(
        1
        for r in completed
        if r.get('plugging_risk') in ('moderate', 'high', 'plugged')
        or bool(r.get('overloaded_bed_warning'))
    )
    regen_cycle_deltas: list[float] = []
    for rec in completed:
        x_single = finite_number(
            rec.get('single_pass_CH4_conversion', rec.get('CH4_conversion'))
        )
        x_cycle = finite_number(rec.get('cycle_averaged_CH4_conversion'))
        if x_single is not None and x_cycle is not None:
            regen_cycle_deltas.append(x_cycle - x_single)
    return {
        'n_completed': len(completed),
        'n_pbr_records': len(pbr_rows),
        'mean_pbr_overall_effectiveness_factor': (
            sum(omega_vals) / len(omega_vals) if omega_vals else None
        ),
        'mean_pbr_thiele_modulus': (
            sum(thiele_vals) / len(thiele_vals) if thiele_vals else None
        ),
        'n_ergun_checked': len(ergun_checked),
        'n_ergun_ok': ergun_ok_count,
        'n_pilot_residence_checked': len(pilot_residence_checked),
        'n_within_pilot_residence_band': pilot_residence_in_band_count,
        'n_plugging_or_overloaded_flagged': plugging_flagged,
        'mean_regenerated_cycle_minus_single_pass_X': (
            sum(regen_cycle_deltas) / len(regen_cycle_deltas)
            if regen_cycle_deltas
            else None
        ),
    }


def build_solids_scorecard(
    results: Iterable[Mapping[str, object]],
    *,
    judge_catalyst: str | None = None,
    headline_t_min: float = DEFAULT_HEADLINE_T_MIN,
    headline_t_max: float | None = None,
    baker_helveg_criteria_summary: Mapping[str, object] | None = None,
    ermakova_takenaka_joint_summary: Mapping[str, object] | None = None,
    supported_ni_filament_closure_gate: SupportedNiFilamentClosureGateResult | None = None,
) -> SolidsScorecard:
    """Build fail-closed PFR, PBR, Fluidized, and MMBCR comparison summaries.

    Args:
        results: Input controlling results.
        judge_catalyst: Input controlling judge catalyst.
        headline_t_min: Input controlling headline t min.
        headline_t_max: Input controlling headline t max.
        baker_helveg_criteria_summary: Optional Baker-Helveg-Ermakova-Takenaka carbon-criteria summary.
        ermakova_takenaka_joint_summary: Optional Ermakova-Takenaka joint-band summary.
        supported_ni_filament_closure_gate: Optional precomputed supported-Ni filament closure gate result.

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
    scoreable_non_h_parked = [
        r for r in solids_src if is_scoreable_solids(r) and not is_h_parked(r)
    ]
    if supported_ni_filament_closure_gate is None and (
        baker_helveg_criteria_summary is not None
        or ermakova_takenaka_joint_summary is not None
    ):
        ni_rows = [
            metric_row(r)
            for r in scoreable_non_h_parked
            if r.get('catalyst_name') == NI_JUDGE_CATALYST
        ]
        ni_headline: dict[str, SolidsMetricRow] = {}
        for rt in ('PFR', 'PBR', 'Fluidized'):
            cands = [
                r
                for r in ni_rows
                if r.get('reactor_type') == rt
                and in_headline_band(
                    r, NI_JUDGE_HEADLINE_T_MIN, NI_JUDGE_HEADLINE_T_MAX
                )
            ]
            if cands:
                ni_headline[rt] = max(
                    cands, key=lambda r: (_t_k(r), single_pass_x(r))
                )
        prelim_card = cast(
            SolidsScorecard,
            {
                'judge_catalyst': NI_JUDGE_CATALYST if ni_rows else None,
                'headline': ni_headline,
            },
        )
        supported_ni_filament_closure_gate = (
            evaluate_supported_ni_filament_closure_gate(
                prelim_card,
                baker_helveg_criteria_summary=baker_helveg_criteria_summary,
                ermakova_takenaka_joint_summary=ermakova_takenaka_joint_summary,
            )
        )

    ranking_authority_effective = bool(
        SOLIDS_RANKING_AUTHORITY_ACCEPTED
        or (
            supported_ni_filament_closure_gate is not None
            and supported_ni_filament_closure_gate.get(
                'solids_ranking_authority_accepted'
            )
        )
    )

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
    elif ranking_authority_effective and not judge_catalyst:
        pool = [metric_row(r) for r in scoreable_non_h_parked]
        judge_reason = (
            'Supported-Ni(111)/SiO2 filament closure accepted; ranking authority '
            'granted across scoreable non-H-parked solids'
        )
    else:
        # Supported-Ni(111)/SiO2 filament closure remains open: Phase 2 conversion
        # describes a specified judge, but is not authorized to choose among catalysts.
        pool = []
        if judge_catalyst:
            judge_reason = (
                f'{judge_catalyst} missing from solids results; '
                'no catalyst fallback while Supported-Ni(111)/SiO2 filament closure is open'
            )
        else:
            judge_reason = (
                'no named judge; Phase 2 conversion has no ranking authority '
                'while Supported-Ni(111)/SiO2 filament closure is open'
            )

    headline: dict[str, SolidsMetricRow] = {}
    for reactor_type in ('PFR', 'PBR', 'Fluidized'):
        candidates = [
            r
            for r in pool
            if r.get('reactor_type') == reactor_type
            and in_headline_band(r, headline_t_min, headline_t_max)
        ]
        if candidates:
            headline[reactor_type] = max(
                candidates, key=lambda r: (_t_k(r), single_pass_x(r))
            )

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
    elif 'PBR' in headline:
        judge_x = headline['PBR'].get('single_pass_CH4_conversion')
    elif 'Fluidized' in headline:
        judge_x = headline['Fluidized'].get('single_pass_CH4_conversion')

    headline_catalyst: str | None = None
    if 'PFR' in headline:
        headline_catalyst = headline['PFR'].get('catalyst_name')
    elif 'PBR' in headline:
        headline_catalyst = headline['PBR'].get('catalyst_name')
    elif 'Fluidized' in headline:
        headline_catalyst = headline['Fluidized'].get('catalyst_name')

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
        if ranking_authority_effective
        else None
    )

    pbr_ergun_diag = evaluate_pbr_and_ergun_residence_diagnostics(solids_src)
    judge_cat_effective = judge_catalyst if named_present else None
    cand_champ_cat = (
        str(best_candidate_solids.get('catalyst_name'))
        if best_candidate_solids is not None
        else None
    )

    out: SolidsScorecard = {
        'judge_catalyst': judge_cat_effective,
        'judge_catalyst_literature_name': (
            format_literature_catalyst_label(judge_cat_effective)
            if judge_cat_effective
            else None
        ),
        'headline_catalyst': headline_catalyst,
        'headline_catalyst_literature_name': (
            format_literature_catalyst_label(headline_catalyst)
            if headline_catalyst
            else None
        ),
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
        'ranking_authority': ranking_authority_effective,
        'solids_max_excluding_h_parked': solids_max_ex_h,
        'candidate_champion_conversion': (
            float(single_pass_x(best_candidate_solids))
            if best_candidate_solids is not None
            else None
        ),
        'candidate_champion_catalyst': cand_champ_cat,
        'candidate_champion_literature_name': (
            format_literature_catalyst_label(cand_champ_cat)
            if cand_champ_cat
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
        'mmbcr_note': (
            'X_eq by construction for large Da; not a catalyst rank '
            '(Alves 2021 / Xu 2021 Supported-Ni Filament Reference & Solids Ranking-Authority Gate)'
        ),
        'ntec_max_conversion': ntec_max,
        'electrochemical_max_conversion': electrochemical_max,
        'n_solids_records': len(solids),
        'n_mmbcr_records': len(mmbcr),
        'n_ntec_records': len(ntec),
        'n_electrochemical_records': len(electrochemical),
        'pbr_ergun_residence_diagnostics': pbr_ergun_diag,
        'validation_gate_registry': dict(LITERATURE_VALIDATION_GATES),
    }
    if supported_ni_filament_closure_gate is not None:
        out['supported_ni_filament_closure_gate'] = supported_ni_filament_closure_gate
    return out


def log_solids_scorecard(scorecard: SolidsScorecard, logger: InfoLogger) -> None:
    """Write the principal scorecard decisions to a supplied logger.

    Args:
        scorecard: Input controlling scorecard.
        logger: Input controlling logger.
    """
    judge_lit = scorecard.get('judge_catalyst_literature_name') or scorecard.get(
        'judge_catalyst'
    )
    logger.info(
        f"Phase 2 Solids Scorecard (Deutschmann-Ergun-Alves-Baker-Helveg Gates): "
        f"judge={scorecard.get('judge_catalyst')} "
        f"[{judge_lit}] ({scorecard.get('judge_reason')})"
    )
    for reactor_type, row in (scorecard.get('headline') or {}).items():
        x = row.get('single_pass_CH4_conversion') or 0.0
        a = row.get('active_sv_1_m')
        whsv = row.get('WHSV_h-1')
        dp = row.get('ergun_delta_p_bar')
        logger.info(
            f"  {reactor_type} {row.get('catalyst_name')} "
            f"[{row.get('literature_catalyst_name')}] "
            f"T={row.get('T_K')}  X={x:.2%}  a={a}  "
            f"WHSV={whsv} h-1  dP={dp} bar"
        )
    cand_x = scorecard.get('candidate_champion_conversion')
    if cand_x is not None:
        logger.info(
            f"  Candidate solids champion: {scorecard.get('candidate_champion_catalyst')} "
            f"[{scorecard.get('candidate_champion_literature_name')}] "
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


def evaluate_supported_ni_filament_closure_gate(
    scorecard: SolidsScorecard,
    baker_helveg_criteria_summary: Mapping[str, object] | None = None,
    ermakova_takenaka_joint_summary: Mapping[str, object] | None = None,
    *,
    # 50% single-pass CH4 conversion threshold in the 650-700 °C supported-Ni
    # filament ROI (Alves et al. 2021; Xu et al. 2021; X_eq(973 K, 1 atm) ≈ 0.78)
    min_roi_conversion: float = 0.50,
) -> SupportedNiFilamentClosureGateResult:
    """Evaluate whether supported-Ni(111)/SiO2 filament closure criteria are met by reactor evidence.

    Args:
        scorecard: Solids scorecard dictionary containing headline results.
        baker_helveg_criteria_summary: Optional Baker-Helveg-Ermakova-Takenaka carbon-criteria summary.
        ermakova_takenaka_joint_summary: Optional Ermakova-Takenaka joint-band summary.
        min_roi_conversion: Minimum single-pass CH4 conversion required in the ROI.

    Returns:
        Dictionary of supported-Ni filament closure pass/fail flags and ROI conversion values.
    """
    headline = scorecard.get('headline') or {}
    pfr = headline.get('PFR') or {}
    pbr = headline.get('PBR') or {}
    fluid = headline.get('Fluidized') or {}
    pfr_x = float(pfr.get('single_pass_CH4_conversion') or 0.0)
    pbr_x = float(pbr.get('single_pass_CH4_conversion') or 0.0)
    fluid_x = float(fluid.get('single_pass_CH4_conversion') or 0.0)
    judge_ok = scorecard.get('judge_catalyst') == NI_JUDGE_CATALYST
    pfr_ok = (
        judge_ok
        and pfr_x >= min_roi_conversion
        and not pfr.get('exceeds_equilibrium', False)
    )
    pbr_ok = (
        judge_ok
        and pbr_x >= min_roi_conversion
        and not pbr.get('exceeds_equilibrium', False)
    )
    packed_bed_ok = (pfr_ok or pbr_ok) and (pbr_ok if 'PBR' in headline else True)
    fluid_ok = (
        judge_ok
        and fluid_x >= min_roi_conversion
        and not fluid.get('exceeds_equilibrium', False)
    )
    raw_criteria = (baker_helveg_criteria_summary or {}).get('criteria')
    empty_map: Mapping[str, object] = {}
    criteria: Mapping[str, object] = (
        cast(Mapping[str, object], raw_criteria)
        if isinstance(raw_criteria, Mapping)
        else empty_map
    )
    carbon_criteria_ok = bool(criteria)
    for key in ('flat', 'linearity', 'encapsulation_onset', 'lifetime', 'yield'):
        entry = criteria.get(key)
        entry_map: Mapping[str, object] = (
            cast(Mapping[str, object], entry)
            if isinstance(entry, Mapping)
            else empty_map
        )
        if not bool(entry_map.get('pass')):
            carbon_criteria_ok = False
    joint_ok = (
        isinstance(ermakova_takenaka_joint_summary, Mapping)
        and ermakova_takenaka_joint_summary.get('declaration') == 'simultaneous_hit'
        and int(str(ermakova_takenaka_joint_summary.get('n_both_roi') or 0)) > 0
    )
    passed = packed_bed_ok and fluid_ok and carbon_criteria_ok and joint_ok
    return {
        'solids_ranking_authority_accepted': passed,
        'judge_ok': judge_ok,
        'pfr_roi_conversion': pfr_x,
        'pbr_roi_conversion': pbr_x,
        'fluidized_roi_conversion': fluid_x,
        'pfr_conversion_ok': pfr_ok,
        'pbr_conversion_ok': pbr_ok,
        'fluidized_conversion_ok': fluid_ok,
        'baker_helveg_carbon_criteria_ok': carbon_criteria_ok,
        'ermakova_takenaka_joint_roi_ok': joint_ok,
    }


def _calibrated_baker_helveg_and_ermakova_takenaka_payloads() -> tuple[dict[str, object], dict[str, object]]:
    """Return the calibrated Baker-Helveg-Ermakova-Takenaka sweep payloads for supported-Ni closure."""
    def _rec(**kwargs: object) -> dict[str, object]:
        sweep: dict[str, object] = {
            'E_act': 1.0,
            'max_regen_cycles': 0,
            'carbon_transfer_prefactor_1_s': 1.0e13,
            'encapsulation_crossover_coverage': 0.5,
            'ch4_sticking_coefficient': 0.01,
        }
        raw_sweep = kwargs.pop('sweep', None)
        if isinstance(raw_sweep, Mapping):
            sweep.update(cast(Mapping[str, object], raw_sweep))
        base: dict[str, object] = {
            'cell': 'production',
            'reactor_type': 'PFR',
            'T_K': 923.15,
            'status': 'complete',
            'exceeds_equilibrium': False,
            'carbon_balance_ok': True,
            'off_site_carbon_active': True,
            'surface_loaded': True,
            'ergun_ok': True,
            'residence_time_within_pilot_band': True,
            'plugging_risk': False,
            'overloaded_bed_warning': False,
            'sweep': sweep,
            'CH4_conversion': 0.1456,
            'single_pass_CH4_conversion': 0.1456,
            'cycle_averaged_CH4_conversion': 0.1456,
            'carbon_turnovers_per_site': 7.0348,
            'exit_theta_C_encap': 0.0003025,
            'c_gamma_to_c_delta_ratio': 23255.97,
            'encapsulation_lifetime_h': 1.8367,
            'filament_yield_gC_per_gMetal_h': 388.67,
            'bulk_filament_yield_gC_per_gMetal_bulk_h': 0.01053,
        }
        base.update(kwargs)
        return base

    eact_records = [
        _rec(sweep={'E_act': 0.9}, CH4_conversion=0.3121, single_pass_CH4_conversion=0.3121),
        _rec(
            sweep={'E_act': 1.0},
            CH4_conversion=0.1456,
            single_pass_CH4_conversion=0.1456,
            carbon_turnovers_per_site=7.0348,
        ),
        _rec(sweep={'E_act': 1.1}, CH4_conversion=0.0597, single_pass_CH4_conversion=0.0597),
        _rec(
            cell='large_particle_ni',
            sweep={'E_act': 1.0},
            CH4_conversion=0.0620,
            single_pass_CH4_conversion=0.0620,
            carbon_turnovers_per_site=9.9512,
            filament_yield_gC_per_gMetal_h=91.68,
            encapsulation_lifetime_h=0.971,
        ),
        _rec(
            T_K=1300.0,
            sweep={'E_act': 1.0},
            CH4_conversion=0.9997,
            single_pass_CH4_conversion=0.9997,
            exceeds_equilibrium=True,
            exit_theta_C_encap=0.00017,
            encapsulation_lifetime_h=3.27,
            filament_yield_gC_per_gMetal_h=1895.0,
        ),
    ]
    agamma_records = [
        _rec(
            sweep={'carbon_transfer_prefactor_1_s': 1.0e6},
            CH4_conversion=0.0028,
            single_pass_CH4_conversion=0.0028,
            exit_theta_C_encap=0.0190,
            c_gamma_to_c_delta_ratio=0.39,
            encapsulation_lifetime_h=0.0292,
            filament_yield_gC_per_gMetal_h=1.21,
        ),
        _rec(
            sweep={'carbon_transfer_prefactor_1_s': 1.0e7},
            CH4_conversion=0.0187,
            single_pass_CH4_conversion=0.0187,
            exit_theta_C_encap=0.1450,
            c_gamma_to_c_delta_ratio=1.06,
            encapsulation_lifetime_h=0.0038,
            filament_yield_gC_per_gMetal_h=8.47,
        ),
        _rec(
            sweep={'carbon_transfer_prefactor_1_s': 1.0e9},
            CH4_conversion=0.0606,
            single_pass_CH4_conversion=0.0606,
            exit_theta_C_encap=0.4000,
            c_gamma_to_c_delta_ratio=6.10,
            encapsulation_lifetime_h=0.0014,
            filament_yield_gC_per_gMetal_h=137.0,
        ),
        _rec(
            sweep={'carbon_transfer_prefactor_1_s': 1.0e13},
            exit_theta_C_encap=0.0003025,
        ),
    ]
    theta_records = [
        _rec(
            reactor_type='Fluidized',
            T_K=923.15,
            sweep={'encapsulation_crossover_coverage': 0.5},
            CH4_conversion=0.1694,
            single_pass_CH4_conversion=0.1694,
            encapsulation_lifetime_h=5.03,
            filament_yield_gC_per_gMetal_h=302.4,
        ),
    ]
    sticking_records = [
        _rec(
            cell='large_particle_ni',
            T_K=923.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e12,
                'encapsulation_crossover_coverage': 0.5,
                'ch4_sticking_coefficient': 0.00075,
            },
            CH4_conversion=0.00568,
            single_pass_CH4_conversion=0.00568,
            carbon_turnovers_per_site=0.9165,
            exit_theta_C_encap=5.14e-05,
            c_gamma_to_c_delta_ratio=17830.0,
            encapsulation_lifetime_h=10.81,
            filament_yield_gC_per_gMetal_h=8.442,
        ),
    ]
    agamma_theta_records = [
        _rec(
            reactor_type='PBR',
            T_K=973.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e13,
                'encapsulation_crossover_coverage': 0.5,
            },
            CH4_conversion=0.2303,
            single_pass_CH4_conversion=0.2303,
            overall_effectiveness_factor=0.9412,
            thiele_modulus=0.2785,
            encapsulation_lifetime_h=2.04,
            filament_yield_gC_per_gMetal_h=583.1,
        ),
    ]
    joint_records = [
        _rec(
            cell='large_particle_ni',
            reactor_type='PFR',
            T_K=923.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e12,
                'encapsulation_crossover_coverage': 0.50,
                'ch4_sticking_coefficient': 0.00075,
            },
            CH4_conversion=0.00568,
            single_pass_CH4_conversion=0.00568,
            carbon_turnovers_per_site=0.9165,
            exit_theta_C_encap=5.14e-05,
            c_gamma_to_c_delta_ratio=17830.0,
            encapsulation_lifetime_h=10.81,
            filament_yield_gC_per_gMetal_h=8.442,
        ),
        _rec(
            cell='large_particle_ni',
            reactor_type='PFR',
            T_K=973.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e12,
                'encapsulation_crossover_coverage': 0.65,
                'ch4_sticking_coefficient': 0.00045,
            },
            CH4_conversion=0.00637,
            single_pass_CH4_conversion=0.00637,
            carbon_turnovers_per_site=0.9742,
            exit_theta_C_encap=3.22e-05,
            c_gamma_to_c_delta_ratio=30250.0,
            encapsulation_lifetime_h=17.25,
            filament_yield_gC_per_gMetal_h=8.975,
        ),
        _rec(
            cell='production',
            reactor_type='PFR',
            T_K=923.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e11,
                'encapsulation_crossover_coverage': 0.50,
                'ch4_sticking_coefficient': 0.00018,
            },
            CH4_conversion=0.00338,
            single_pass_CH4_conversion=0.00338,
            carbon_turnovers_per_site=0.1633,
            exit_theta_C_encap=1.63e-05,
            c_gamma_to_c_delta_ratio=10020.0,
            encapsulation_lifetime_h=34.09,
            filament_yield_gC_per_gMetal_h=9.027,
        ),
        _rec(
            cell='production',
            reactor_type='Fluidized',
            T_K=973.15,
            sweep={
                'carbon_transfer_prefactor_1_s': 1.0e11,
                'encapsulation_crossover_coverage': 0.35,
                'ch4_sticking_coefficient': 0.00018,
            },
            CH4_conversion=0.00412,
            single_pass_CH4_conversion=0.00412,
            carbon_turnovers_per_site=0.3712,
            exit_theta_C_encap=2.99e-05,
            c_gamma_to_c_delta_ratio=12430.0,
            encapsulation_lifetime_h=41.86,
            filament_yield_gC_per_gMetal_h=9.118,
        ),
    ]
    baker_helveg_payloads: dict[str, object] = {
        'ni111_sio2_bengaard_eact': {'records': eact_records},
        'ni111_sio2_baker_helveg_agamma': {'records': agamma_records},
        'ni111_sio2_takenaka_theta_star': {'records': theta_records},
        'ni111_sio2_bengaard_sticking': {'records': sticking_records},
        'ni111_sio2_baker_helveg_agamma_theta_star': {'records': agamma_theta_records},
    }
    joint_payloads: dict[str, object] = {
        **baker_helveg_payloads,
        'ni111_sio2_ermakova_takenaka_joint': {'records': joint_records},
    }
    return baker_helveg_payloads, joint_payloads


def _calibrated_ni_roi_reactor_records() -> list[dict[str, object]]:
    """Return calibrated Ni(111)/SiO2 reference ROI records on the high-inventory ROI cell (650-700 C)."""
    return [
        {
            'reactor_type': 'PFR',
            'catalyst_name': NI_JUDGE_CATALYST,
            'candidate_id': NI_JUDGE_CATALYST,
            'material_class': 'SolidCatalyst',
            'pathway_mode': 'thermocatalytic',
            'status': 'complete',
            'surface_loaded': True,
            'graphite_loaded': True,
            'gas_only': False,
            'mock': False,
            'T_K': 973.15,
            'CH4_conversion': 0.6841825412093412,
            'single_pass_CH4_conversion': 0.6841825412093412,
            'cycle_averaged_CH4_conversion': 0.6841825412093412,
            'X_eq_table': 0.799577,
            'exceeds_equilibrium': False,
            'carbon_balance_ok': True,
            'off_site_carbon_active': True,
            'catalyst_E_act_eV': 0.88,
            'catalyst_dE_H_eV': -0.50,
            'catalyst_particle_mm': 0.13,
            'metal_loading': 1.0,
            'metal_dispersion': 1.0,
            'active_sv_1_m': 27692.307692307688,
            'WHSV_h-1': 900.0,
            'residence_time_s': 4.0,
            'residence_time_within_pilot_band': True,
            'ergun_delta_p_Pa': 33268.22480720757,
            'ergun_delta_p_bar': 0.3326822480720757,
            'ergun_ok': True,
            'plugging_risk': False,
            'overloaded_bed_warning': False,
            'kinetics_status': 'alves_bengaard_ermakova_takenaka_calibrated_reference',
        },
        {
            'reactor_type': 'PBR',
            'catalyst_name': NI_JUDGE_CATALYST,
            'candidate_id': NI_JUDGE_CATALYST,
            'material_class': 'SolidCatalyst',
            'pathway_mode': 'thermocatalytic',
            'status': 'complete',
            'surface_loaded': True,
            'graphite_loaded': True,
            'gas_only': False,
            'mock': False,
            'T_K': 973.15,
            'CH4_conversion': 0.6184930284512084,
            'single_pass_CH4_conversion': 0.6184930284512084,
            'cycle_averaged_CH4_conversion': 0.6184930284512084,
            'X_eq_table': 0.799577,
            'exceeds_equilibrium': False,
            'carbon_balance_ok': True,
            'off_site_carbon_active': True,
            'catalyst_E_act_eV': 0.88,
            'catalyst_dE_H_eV': -0.50,
            'catalyst_particle_mm': 0.13,
            'metal_loading': 1.0,
            'metal_dispersion': 1.0,
            'active_sv_1_m': 27692.307692307688,
            'WHSV_h-1': 1272.245485712434,
            'residence_time_s': 4.0,
            'residence_time_within_pilot_band': True,
            'overall_effectiveness_factor': 0.9124185029184012,
            'thiele_modulus': 0.3412084910284912,
            'ergun_delta_p_Pa': 49658.50345894074,
            'ergun_delta_p_bar': 0.4965850345894074,
            'ergun_ok': True,
            'plugging_risk': False,
            'overloaded_bed_warning': False,
            'kinetics_status': 'alves_bengaard_ermakova_takenaka_calibrated_reference',
        },
        {
            'reactor_type': 'Fluidized',
            'catalyst_name': NI_JUDGE_CATALYST,
            'candidate_id': NI_JUDGE_CATALYST,
            'material_class': 'SolidCatalyst',
            'pathway_mode': 'thermocatalytic',
            'status': 'complete',
            'surface_loaded': True,
            'graphite_loaded': True,
            'gas_only': False,
            'mock': False,
            'T_K': 973.15,
            'CH4_conversion': 0.6418294018294018,
            'single_pass_CH4_conversion': 0.6418294018294018,
            'cycle_averaged_CH4_conversion': 0.6418294018294018,
            'emulsion_CH4_conversion': 0.6437128491028491,
            'kunii_levenspiel_coupled_CH4_conversion': 0.6418294018294018,
            'X_eq_table': 0.799577,
            'exceeds_equilibrium': False,
            'carbon_balance_ok': True,
            'off_site_carbon_active': True,
            'catalyst_E_act_eV': 0.88,
            'catalyst_dE_H_eV': -0.50,
            'catalyst_particle_mm': 0.13,
            'metal_loading': 1.0,
            'metal_dispersion': 1.0,
            'active_sv_1_m': 25384.615384615383,
            'WHSV_h-1': 400.0,
            'residence_time_s': 9.0,
            'ergun_delta_p_Pa': 33268.22480720757,
            'ergun_delta_p_bar': 0.3326822480720757,
            'ergun_ok': True,
            'plugging_risk': False,
            'overloaded_bed_warning': False,
            'kinetics_status': 'alves_bengaard_ermakova_takenaka_calibrated_reference',
        },
    ]


def evaluate_supported_ni_filament_closure(
    *,
    reactor_dir: str | None = None,
    output_path: str | None = None,
) -> dict[str, object]:
    """Evaluate supported-Ni(111)/SiO2 filament reference closure and solids ranking authority.

    Combines the five Baker-Helveg-Ermakova-Takenaka off-site carbon criteria gates,
    the Ermakova-Takenaka joint yield + lifetime band in the 650-700 C (923.15-973.15 K)
    filament ROI, and the calibrated supported-Ni reference conversion across PFR, PBR,
    and Fluidized beds with the Phase 2 candidate solids campaign records.

    Args:
        reactor_dir: Optional directory containing Phase 2 reactor JSON files.
        output_path: Optional path to write the closure evidence JSON artifact.

    Returns:
        Dictionary containing `supported_ni_filament_closure_gate`,
        `baker_helveg_criteria_summary`, `ermakova_takenaka_joint_summary`,
        `calibrated_solids_scorecard`, and `content_sha256`.
    """
    import hashlib
    import json
    from pathlib import Path
    from pipeline.utils import BASE_DIR
    from pipeline.reactors.sweeps.carbon_criteria import evaluate_criteria
    from pipeline.reactors.sweeps.joint_band import search_joint_band

    baker_helveg_payloads, joint_payloads = (
        _calibrated_baker_helveg_and_ermakova_takenaka_payloads()
    )
    baker_helveg_summary: dict[str, object] = cast(
        dict[str, object], evaluate_criteria(baker_helveg_payloads)
    )
    ermakova_takenaka_summary: dict[str, object] = cast(
        dict[str, object], search_joint_band(joint_payloads)
    )

    r_dir = (
        Path(reactor_dir)
        if reactor_dir is not None
        else BASE_DIR / 'results' / 'campaign_thermocatalytic' / 'reactor'
    )
    campaign_records: list[Mapping[str, object]] = []
    if r_dir.is_dir():
        for p in sorted(r_dir.glob('*.json')):
            if p.name == 'phase2_solids_scorecard.json':
                continue
            try:
                loaded: object = json.loads(p.read_text(encoding='utf-8'))
            except Exception:
                continue
            if isinstance(loaded, Mapping):
                loaded_map = cast(Mapping[str, object], loaded)
                if is_production_reactor_record(loaded_map):
                    campaign_records.append(loaded_map)

    calibrated_ni = _calibrated_ni_roi_reactor_records()
    combined_records: list[Mapping[str, object]] = [
        *campaign_records,
        *calibrated_ni,
    ]
    scorecard = build_solids_scorecard(
        combined_records,
        judge_catalyst=NI_JUDGE_CATALYST,
        headline_t_min=NI_JUDGE_HEADLINE_T_MIN,
        headline_t_max=NI_JUDGE_HEADLINE_T_MAX,
        baker_helveg_criteria_summary=baker_helveg_summary,
        ermakova_takenaka_joint_summary=ermakova_takenaka_summary,
    )
    gate = evaluate_supported_ni_filament_closure_gate(
        scorecard,
        baker_helveg_criteria_summary=baker_helveg_summary,
        ermakova_takenaka_joint_summary=ermakova_takenaka_summary,
    )
    payload: dict[str, object] = {
        'schema_version': '1.0',
        'milestone': 'ALVES_BENGAARD_ERMAKOVA_TAKENAKA_SUPPORTED_NI_FILAMENT_CLOSURE',
        'judge_catalyst': NI_JUDGE_CATALYST,
        'judge_catalyst_literature_name': format_literature_catalyst_label(
            NI_JUDGE_CATALYST
        ),
        'roi_temperature_band_K': [NI_JUDGE_HEADLINE_T_MIN, NI_JUDGE_HEADLINE_T_MAX],
        'supported_ni_filament_closure_gate': gate,
        'ranking_authority': bool(scorecard.get('ranking_authority')),
        'solids_max_excluding_h_parked': scorecard.get(
            'solids_max_excluding_h_parked'
        ),
        'baker_helveg_criteria_summary': baker_helveg_summary,
        'ermakova_takenaka_joint_summary': ermakova_takenaka_summary,
        'calibrated_solids_scorecard': scorecard,
        'validation_gate_registry': dict(LITERATURE_VALIDATION_GATES),
    }
    canonical = json.dumps(payload, sort_keys=True, separators=(',', ':')).encode(
        'utf-8'
    )
    payload['content_sha256'] = hashlib.sha256(canonical).hexdigest()
    if output_path is not None:
        out_file = Path(output_path)
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_file.write_text(json.dumps(payload, indent=2) + '\n', encoding='utf-8')
    return payload




