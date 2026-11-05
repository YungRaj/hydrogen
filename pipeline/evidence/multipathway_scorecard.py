"""End-to-end 7-pathway discovery, reactor, DFT, PEMFC, and TEA champion scorecard."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import cast

from pipeline.data_models.evidence import (
    CoupledPathwayTEARow,
    DFTCalibratedPathwayRow,
    DFTCalibratedReactorMetrics,
    DiscoveryChampionRow,
    MultipathwayChampionScorecard,
    PathwayLeaderboardRow,
    PEMFCStackSummaryRow,
)
from pipeline.data_models.reactors import (
    ReactorResult,
    SupportedNiFilamentClosureSummary,
)
from pipeline.economics.tea import LHV_H2_KWH_PER_KG, estimate_scenario_range


def build_multipathway_champion_scorecard(
    results_root: str = 'results',
    *,
    stack_cells: int = 400,
    top_k_catalysts: int = 11,
) -> MultipathwayChampionScorecard:
    """Build an end-to-end 6-pathway champion portfolio and TEA scorecard.

    Aggregates Phase 1 turquoise discovery rows (annotated with validated
    interfacial carbon displacement ``|c_height_above_surface_A|`` and
    ``interfacial_pyrolysis_score_eV``), all six Phase 2 reactor pathway modes
    (``thermocatalytic_pfr``, ``thermocatalytic_fluidized``,
    ``thermocatalytic``, ``mmbcr``, ``ntec``, ``electrochemical``) with
    scenario-range levelized hydrogen cost estimates, and Phase 5 PEMFC
    cathode/stack techno-economics including coupled H2-to-power LCOE.

    Args:
        results_root: Root directory containing campaign artifacts.
        stack_cells: Number of cells in the modeled PEMFC stack.
        top_k_catalysts: Number of top discovery catalysts to retain.

    Returns:
        JSON-serializable dictionary with ``discovery_champions``,
        ``pathway_leaderboard``, ``pemfc_stack_summary``, and
        ``coupled_turquoise_pemfc_tea``.
    """
    import pandas as pd
    from pipeline.evidence.pyrolysis_benchmark import (
        K_B_EV_K,
        calibrated_interfacial_pyrolysis_score_eV,
    )
    from pipeline.fuel_cell import composite_score
    from pipeline.fuel_cell.pemfc import (
        PEMFCConfig,
        TAFEL_SLOPE_BY_CLASS,
        simulate_pemfc,
    )
    from pipeline.fuel_cell.stack import StackConfig, model_stack
    from pipeline.reactors.eligibility import is_usable_result
    from pipeline.reactors.scorecard import (
        NI_JUDGE_CATALYST,
        evaluate_pbr_and_ergun_residence_diagnostics,
        evaluate_supported_ni_filament_closure,
        is_h_parked,
        is_production_reactor_record,
        is_solids_run,
        single_pass_x,
    )
    from pipeline.screening.stage_selection import select_for_pyrolysis_reactor
    from pipeline.search.discovery import (
        LITERATURE_VALIDATION_GATES,
        NAMED_REFERENCE_CATALYSTS,
        format_literature_catalyst_label,
    )
    from pipeline.search.scope import scope_pyrolysis_pool

    root = Path(results_root)
    campaign_dirs = [
        d
        for d in (
            root,
            root / 'campaign_thermocatalytic',
            root / 'campaign_mmbcr',
            root / 'campaign_ntec',
            root / 'campaign_electrochemical',
        )
        if d.is_dir()
    ]

    # Load active-learning QE DFT calibration artifact if available
    al_cal_path = root / 'dft' / 'active_learning_calibration.json'
    al_turq_class_bias: dict[str, float] = {}
    if al_cal_path.is_file():
        try:
            al_doc = json.loads(al_cal_path.read_text(encoding='utf-8'))
            turq_app = (
                al_doc.get('applications', {}).get('turquoise_hydrogen', {})
                if isinstance(al_doc, dict)
                else {}
            )
            raw_cb = (
                turq_app.get('class_bias_delta_eV')
                or turq_app.get('class_bias_eV')
                or (al_doc.get('class_bias_delta_eV') if isinstance(al_doc, dict) else {})
                or {}
            )
            if isinstance(raw_cb, dict):
                al_turq_class_bias = {
                    str(k): float(v)
                    for k, v in raw_cb.items()
                    if isinstance(v, (int, float)) and math.isfinite(float(v))
                }
        except Exception:
            pass

    # 1. Phase 1 Discovery Champions with Interfacial Carbon Score & QE DFT Validation
    screening_frames: list[pd.DataFrame] = []
    for cdir in campaign_dirs:
        ev_path = cdir / 'screening' / 'branch_ranker_evidence.csv'
        if ev_path.is_file():
            try:
                screening_frames.append(pd.read_csv(ev_path))
            except Exception:
                pass
    if not screening_frames:
        prospective_dir = root / 'screening' / 'prospective'
        if prospective_dir.is_dir():
            for p_csv in sorted(prospective_dir.glob('*_pyrolysis.csv')):
                try:
                    screening_frames.append(pd.read_csv(p_csv))
                except Exception:
                    pass
    discovery_champions: list[DiscoveryChampionRow] = []
    candidate_Lookup: dict[str, dict[str, object]] = {}
    if screening_frames:
        merged_df = pd.concat(screening_frames, ignore_index=True)
        if 'genome' in merged_df.columns:
            merged_df = merged_df.drop_duplicates('genome', keep='first').reset_index(
                drop=True
            )
        pool, _ = scope_pyrolysis_pool(merged_df)
        for _, prow in pool.iterrows():
            cid_key = str(prow.get('candidate_id', '')).strip()
            if cid_key:
                candidate_Lookup[cid_key] = dict(prow.to_dict())
        selected = select_for_pyrolysis_reactor(
            pool, top_k_catalysts, primary='E_act', min_per_class=1
        )
        for rank_idx, (_, row) in enumerate(selected.iterrows(), 1):
            cid_str = str(row.get('candidate_id', ''))
            c_disp = row.get('interfacial_carbon_displacement_A')
            i_score = row.get('interfacial_pyrolysis_score_eV')
            cal_e_act = row.get('calibrated_activation_eV')
            cal_score = row.get('calibrated_interfacial_score_eV')
            if (
                (cal_e_act is None or pd.isna(cal_e_act))
                or (cal_score is None or pd.isna(cal_score))
            ):
                try:
                    cal_info = calibrated_interfacial_pyrolysis_score_eV(row.to_dict())
                    cal_e_act = cal_info['calibrated_activation_eV']
                    cal_score = cal_info['calibrated_interfacial_score_eV']
                    if c_disp is None or pd.isna(c_disp):
                        c_disp = cal_info['interfacial_carbon_displacement_A']
                    if i_score is None or pd.isna(i_score):
                        i_score = cal_info['uncalibrated_interfacial_score_eV']
                except Exception:
                    pass
            coking = row.get('coking_index')
            champ_entry: DiscoveryChampionRow = {
                'rank': rank_idx,
                'candidate_id': cid_str,
                'literature_catalyst_name': format_literature_catalyst_label(
                    cid_str,
                    material_class=str(row.get('material_class', '')),
                    genome=str(row.get('genome', '')),
                ),
                'material_class': str(row.get('material_class', '')),
                'genome': str(row.get('genome', '')),
                'E_act_eV': round(
                    float(raw_e_act)
                    if isinstance(raw_e_act := row.get('E_act', 0.0), (int, float, str))
                    else 0.0,
                    4,
                ),
                'calibrated_activation_eV': (
                    None
                    if cal_e_act is None or pd.isna(cal_e_act)
                    else round(float(cal_e_act), 4)
                ),
                'interfacial_carbon_displacement_A': (
                    None
                    if c_disp is None or pd.isna(c_disp)
                    else round(float(c_disp), 4)
                ),
                'interfacial_pyrolysis_score_eV': (
                    None
                    if i_score is None or pd.isna(i_score)
                    else round(float(i_score), 4)
                ),
                'calibrated_interfacial_score_eV': (
                    None
                    if cal_score is None or pd.isna(cal_score)
                    else round(float(cal_score), 4)
                ),
                'coking_index_eV': (
                    None
                    if coking is None or pd.isna(coking)
                    else round(float(coking), 4)
                ),
            }
            dft_path = root / 'dft' / f'campaign_{cid_str}_dft.json'
            if dft_path.is_file():
                try:
                    dft_rec = json.loads(dft_path.read_text(encoding='utf-8'))
                    if isinstance(dft_rec, dict):
                        champ_entry['qe_dft_converged'] = bool(
                            dft_rec.get('converged', False)
                        )
                        e_ry = dft_rec.get('dft_energy_Ry')
                        e_ev = dft_rec.get('dft_energy_eV')
                        f_max = dft_rec.get('max_force_Ry_bohr')
                        champ_entry['qe_dft_energy_Ry'] = (
                            round(float(e_ry), 8)
                            if isinstance(e_ry, (int, float))
                            else None
                        )
                        champ_entry['qe_dft_energy_eV'] = (
                            round(float(e_ev), 4)
                            if isinstance(e_ev, (int, float))
                            else None
                        )
                        champ_entry['qe_max_force_Ry_bohr'] = (
                            round(float(f_max), 6)
                            if isinstance(f_max, (int, float))
                            else None
                        )
                        neb_conv = dft_rec.get('neb_converged')
                        neb_fwd = dft_rec.get('neb_forward_barrier_eV')
                        if neb_conv is not None:
                            champ_entry['qe_neb_converged'] = bool(neb_conv)
                        if isinstance(neb_fwd, (int, float)) and math.isfinite(
                            float(neb_fwd)
                        ):
                            champ_entry['qe_neb_barrier_eV'] = round(float(neb_fwd), 4)
                except Exception:
                    pass
            discovery_champions.append(champ_entry)
            if cid_str:
                candidate_Lookup[cid_str] = {
                    **candidate_Lookup.get(cid_str, {}),
                    **champ_entry,
                }

    # 2. Phase 2 6-Pathway Reactor & TEA Leaderboard (Raw + DFT-Calibrated)
    reactor_records: list[ReactorResult] = []
    for cdir in campaign_dirs:
        rdir = cdir / 'reactor'
        if not rdir.is_dir():
            continue
        for jpath in sorted(rdir.glob('*.json')):
            if jpath.name == 'phase2_solids_scorecard.json':
                continue
            try:
                rec = json.loads(jpath.read_text(encoding='utf-8'))
            except Exception:
                continue
            if (
                isinstance(rec, dict)
                and is_production_reactor_record(rec)
                and not rec.get('mock', False)
                and is_usable_result(rec)
            ):
                reactor_records.append(cast(ReactorResult, rec))

    pathway_filters = {
        'thermocatalytic_pfr': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'PFR'
            and not is_h_parked(r)
            and r.get('catalyst_name') != NI_JUDGE_CATALYST
        ),
        'thermocatalytic_pbr': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'PBR'
            and not is_h_parked(r)
            and r.get('catalyst_name') != NI_JUDGE_CATALYST
        ),
        'thermocatalytic_fluidized': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'Fluidized'
            and not is_h_parked(r)
            and r.get('catalyst_name') != NI_JUDGE_CATALYST
        ),
        'thermocatalytic': lambda r: (
            is_solids_run(r)
            and not is_h_parked(r)
            and r.get('catalyst_name') != NI_JUDGE_CATALYST
        ),
        'mmbcr': lambda r: r.get('reactor_type') == 'MMBCR',
        'ntec': lambda r: r.get('reactor_type') == 'NTEC',
        'electrochemical': lambda r: r.get('reactor_type') == 'Electrochemical',
    }

    def _resolve_candidate_neb_evidence(
        cat_name: str, cat_meta: dict[str, object]
    ) -> tuple[bool, float | None]:
        """Return (neb_converged, neb_forward_barrier_eV) for a catalyst if present."""
        if not cat_name:
            return False, None
        neb_conv = bool(
            cat_meta.get('qe_neb_converged') or cat_meta.get('neb_converged')
        )
        raw_bar = (
            cat_meta.get('qe_neb_barrier_eV')
            if cat_meta.get('qe_neb_barrier_eV') is not None
            else cat_meta.get('neb_forward_barrier_eV')
        )
        if neb_conv and isinstance(raw_bar, (int, float)) and float(raw_bar) > 0.0:
            return True, float(raw_bar)
        for dft_candidate_path in (
            root / 'dft' / f'campaign_{cat_name}_dft.json',
            root / 'dft' / f'{cat_name}_dft.json',
        ):
            if dft_candidate_path.is_file():
                try:
                    d_rec = json.loads(dft_candidate_path.read_text(encoding='utf-8'))
                    if isinstance(d_rec, dict):
                        c_flag = bool(
                            d_rec.get('neb_converged') or d_rec.get('qe_neb_converged')
                        )
                        b_val = (
                            d_rec.get('neb_forward_barrier_eV')
                            if d_rec.get('neb_forward_barrier_eV') is not None
                            else d_rec.get('qe_neb_barrier_eV')
                        )
                        if (
                            c_flag
                            and isinstance(b_val, (int, float))
                            and float(b_val) > 0.0
                        ):
                            return True, float(b_val)
                except Exception:
                    pass
        return False, None

    def _compute_dft_calibrated_reactor_metrics(
        rec: ReactorResult,
    ) -> DFTCalibratedReactorMetrics:
        cat_name = str(rec.get('catalyst_name', '')).strip()
        cat_meta = candidate_Lookup.get(cat_name, {})
        mat_cls = str(
            rec.get('material_class')
            or cat_meta.get('material_class')
            or ('MoltenMetal' if rec.get('reactor_type') == 'MMBCR' else 'HEA')
        )
        t_val = float(rec.get('T_K') or 1100.0)
        t_eff = max(t_val, 300.0)
        neb_conv, neb_barrier_ev = _resolve_candidate_neb_evidence(cat_name, cat_meta)
        if neb_conv and neb_barrier_ev is not None:
            raw_e_act = (
                rec.get('catalyst_E_act_eV')
                or cat_meta.get('E_act_eV')
                or cat_meta.get('E_act')
                or 0.0
            )
            e_act_bep = (
                float(raw_e_act)
                if isinstance(raw_e_act, (int, float))
                and math.isfinite(float(raw_e_act))
                else 0.0
            )
            delta_e_dft = max(0.0, float(neb_barrier_ev) - e_act_bep)
            kinetics_status = 'dft_neb_complete'
        else:
            # Class-level QE active-learning barrier correction (plus interfacial carbon penalty)
            default_bias = 0.228733 if mat_cls == 'MoltenMetal' else 0.243114
            delta_e_dft = float(al_turq_class_bias.get(mat_cls, default_bias))
            kinetics_status = 'qe_dft_calibrated_bep'
        c_disp_val = cat_meta.get('interfacial_carbon_displacement_A')
        c_disp_a = (
            float(c_disp_val)
            if isinstance(c_disp_val, (int, float)) and math.isfinite(float(c_disp_val))
            else 0.0
        )
        # Effective barrier shift includes the QE AL class bias (or direct NEB shift) and 0.15 eV/A interfacial coking drag
        effective_barrier_shift_ev = max(0.0, delta_e_dft + 0.15 * c_disp_a)
        if str(rec.get('reactor_type', '')) == 'NTEC':
            mp_scr = rec.get('multiphysics_screening')
            mp_out = (
                mp_scr.get('outputs', {})
                if isinstance(mp_scr, dict) and isinstance(mp_scr.get('outputs'), dict)
                else {}
            )
            raw_fe = mp_out.get('field_rate_enhancement')
            if isinstance(raw_fe, (int, float)) and float(raw_fe) > 1.0:
                stark_offset_ev = 0.5 * K_B_EV_K * t_eff * math.log(float(raw_fe))
                effective_barrier_shift_ev = max(
                    0.0, effective_barrier_shift_ev - stark_offset_ev
                )
        rate_scale = math.exp(-effective_barrier_shift_ev / (K_B_EV_K * t_eff))
        x_raw = min(max(float(single_pass_x(rec)), 0.0), 0.9999)
        da_raw = -math.log(max(1.0 - x_raw, 1e-6))

        raw_phi_val = rec.get('thiele_modulus')
        raw_eta_val = rec.get('overall_effectiveness_factor')
        if (
            isinstance(raw_phi_val, (int, float))
            and math.isfinite(float(raw_phi_val))
            and float(raw_phi_val) > 0.0
        ):
            phi_raw = float(raw_phi_val)
            eta_raw = (
                float(raw_eta_val)
                if isinstance(raw_eta_val, (int, float))
                and math.isfinite(float(raw_eta_val))
                and float(raw_eta_val) > 0.0
                else math.tanh(phi_raw) / phi_raw
            )
            phi_cal = max(phi_raw * math.sqrt(rate_scale), 1e-6)
            eta_cal = math.tanh(phi_cal) / phi_cal
            da_cal = da_raw * rate_scale * (eta_cal / max(eta_raw, 1e-6))
            cal_phi_out: float | None = round(phi_cal, 4)
            cal_eta_out: float | None = round(eta_cal, 4)
        else:
            da_cal = da_raw * rate_scale
            cal_phi_out = None
            cal_eta_out = None

        x_cal = 1.0 - math.exp(-da_cal)
        return {
            'material_class': mat_cls,
            'dft_barrier_shift_eV': round(effective_barrier_shift_ev, 4),
            'arrhenius_rate_scaling': round(rate_scale, 4),
            'dft_calibrated_single_pass_CH4_conversion': (
                round(x_cal, 4) if x_cal >= 1e-4 else round(x_cal, 6)
            ),
            '_x_cal_exact': x_cal,
            'dft_calibrated_overall_effectiveness_factor': cal_eta_out,
            'dft_calibrated_thiele_modulus': cal_phi_out,
            'kinetics_status': kinetics_status,
        }

    pathway_leaderboard: dict[str, PathwayLeaderboardRow] = {}
    dft_calibrated_pathway_leaderboard: dict[str, DFTCalibratedPathwayRow] = {}
    for mode_name, predicate in pathway_filters.items():
        matching = [r for r in reactor_records if predicate(r)]
        if not matching:
            pathway_leaderboard[mode_name] = {'status': 'no_usable_records'}
            dft_calibrated_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
            }
            continue
        best_rec = max(
            matching,
            key=lambda r: (single_pass_x(r), -float(r.get('T_K') or 0.0)),
        )
        conv = float(single_pass_x(best_rec))
        t_rec = float(best_rec.get('T_K') or 0.0)
        raw_sc = best_rec.get('solid_C_selectivity')
        s_c_rec = (
            float(raw_sc)
            if isinstance(raw_sc, (int, float)) and 0.0 <= float(raw_sc) <= 1.0
            else 1.0
        )
        raw_spec_e = best_rec.get('specific_energy_kWh_kg_H2')
        spec_e = (
            float(raw_spec_e)
            if isinstance(raw_spec_e, (int, float)) and 0.0 < float(raw_spec_e) < 1e6
            else None
        )
        tea_range = (
            estimate_scenario_range(
                conv,
                temperature_K=t_rec if t_rec > 298.15 else None,
                solid_c_selectivity=s_c_rec,
                specific_energy_kwh_kg_h2=spec_e,
            )
            if 0.0 < conv <= 1.0
            else None
        )
        raw_cycle_x = best_rec.get('cycle_averaged_CH4_conversion')
        raw_eta = best_rec.get('overall_effectiveness_factor')
        raw_phi = best_rec.get('thiele_modulus')
        raw_dp = best_rec.get('ergun_delta_p_bar')
        best_cat_id = str(best_rec.get('catalyst_name', ''))
        best_cat_meta = candidate_Lookup.get(best_cat_id, {})
        pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': best_cat_id,
            'champion_catalyst_literature_name': format_literature_catalyst_label(
                best_cat_id,
                material_class=best_rec.get('material_class')
                or best_cat_meta.get('material_class'),
                genome=best_rec.get('genome') or best_cat_meta.get('genome'),
            ),
            'reactor_type': str(best_rec.get('reactor_type', '')),
            'T_K': t_rec,
            'single_pass_CH4_conversion': round(conv, 4),
            'cycle_averaged_CH4_conversion': (
                round(float(raw_cycle_x), 4)
                if isinstance(raw_cycle_x, (int, float))
                and math.isfinite(float(raw_cycle_x))
                else None
            ),
            'overall_effectiveness_factor': (
                round(float(raw_eta), 4)
                if isinstance(raw_eta, (int, float)) and math.isfinite(float(raw_eta))
                else None
            ),
            'thiele_modulus': (
                round(float(raw_phi), 4)
                if isinstance(raw_phi, (int, float)) and math.isfinite(float(raw_phi))
                else None
            ),
            'ergun_delta_p_bar': (
                round(float(raw_dp), 6)
                if isinstance(raw_dp, (int, float)) and math.isfinite(float(raw_dp))
                else None
            ),
            'H2_selectivity': round(float(best_rec.get('H2_selectivity') or 1.0), 4),
            'reactor_evidence_tier': str(
                best_rec.get('reactor_evidence_tier', 'cantera_1d')
            ),
            'kinetics_status': str(best_rec.get('kinetics_status', 'screening_bep')),
            'tea_base_h2_usd_kg': (
                tea_range['estimates']['base']['h2_cost_usd_kg']
                if tea_range
                else None
            ),
            'tea_optimistic_h2_usd_kg': (
                tea_range['estimates']['optimistic']['h2_cost_usd_kg']
                if tea_range
                else None
            ),
            'tea_conservative_h2_usd_kg': (
                tea_range['estimates']['conservative']['h2_cost_usd_kg']
                if tea_range
                else None
            ),
            'energy_input_kwh_kg_h2_base': (
                round(tea_range['estimates']['base']['energy_input_kwh_kg_h2'], 3)
                if tea_range
                else None
            ),
        }

        # Select champion under DFT-calibrated conversion
        scored_cal = [
            (r, _compute_dft_calibrated_reactor_metrics(r)) for r in matching
        ]
        best_cal_rec, best_cal_metrics = max(
            scored_cal,
            key=lambda item: (
                item[1]['_x_cal_exact'],
                -float(item[0].get('T_K') or 0.0),
            ),
        )
        cal_conv = float(
            best_cal_metrics['dft_calibrated_single_pass_CH4_conversion']
        )
        cal_t_rec = float(best_cal_rec.get('T_K') or 0.0)
        cal_raw_sc = best_cal_rec.get('solid_C_selectivity')
        cal_s_c = (
            float(cal_raw_sc)
            if isinstance(cal_raw_sc, (int, float))
            and 0.0 <= float(cal_raw_sc) <= 1.0
            else 1.0
        )
        cal_raw_spec_e = best_cal_rec.get('specific_energy_kWh_kg_H2')
        cal_spec_e = (
            float(cal_raw_spec_e)
            if isinstance(cal_raw_spec_e, (int, float))
            and 0.0 < float(cal_raw_spec_e) < 1e6
            else None
        )
        cal_tea_range = (
            estimate_scenario_range(
                cal_conv,
                temperature_K=cal_t_rec if cal_t_rec > 298.15 else None,
                solid_c_selectivity=cal_s_c,
                specific_energy_kwh_kg_h2=cal_spec_e,
            )
            if 0.0 < cal_conv <= 1.0
            else None
        )
        cal_dp = best_cal_rec.get('ergun_delta_p_bar')
        best_cal_cat_id = str(best_cal_rec.get('catalyst_name', ''))
        best_cal_meta = candidate_Lookup.get(best_cal_cat_id, {})
        dft_calibrated_pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': best_cal_cat_id,
            'champion_catalyst_literature_name': format_literature_catalyst_label(
                best_cal_cat_id,
                material_class=best_cal_metrics['material_class'],
                genome=best_cal_rec.get('genome') or best_cal_meta.get('genome'),
            ),
            'material_class': best_cal_metrics['material_class'],
            'reactor_type': str(best_cal_rec.get('reactor_type', '')),
            'T_K': cal_t_rec,
            'uncalibrated_single_pass_CH4_conversion': round(
                float(single_pass_x(best_cal_rec)), 4
            ),
            'dft_barrier_shift_eV': best_cal_metrics['dft_barrier_shift_eV'],
            'arrhenius_rate_scaling': best_cal_metrics['arrhenius_rate_scaling'],
            'dft_calibrated_single_pass_CH4_conversion': cal_conv,
            'dft_calibrated_overall_effectiveness_factor': best_cal_metrics[
                'dft_calibrated_overall_effectiveness_factor'
            ],
            'dft_calibrated_thiele_modulus': best_cal_metrics[
                'dft_calibrated_thiele_modulus'
            ],
            'ergun_delta_p_bar': (
                round(float(cal_dp), 6)
                if isinstance(cal_dp, (int, float)) and math.isfinite(float(cal_dp))
                else None
            ),
            'H2_selectivity': round(
                float(best_cal_rec.get('H2_selectivity') or 1.0), 4
            ),
            'kinetics_status': best_cal_metrics.get(
                'kinetics_status', 'qe_dft_calibrated_bep'
            ),
            'dft_calibrated_tea_base_h2_usd_kg': (
                cal_tea_range['estimates']['base']['h2_cost_usd_kg']
                if cal_tea_range
                else None
            ),
            'dft_calibrated_tea_optimistic_h2_usd_kg': (
                cal_tea_range['estimates']['optimistic']['h2_cost_usd_kg']
                if cal_tea_range
                else None
            ),
            'dft_calibrated_tea_conservative_h2_usd_kg': (
                cal_tea_range['estimates']['conservative']['h2_cost_usd_kg']
                if cal_tea_range
                else None
            ),
            'dft_calibrated_energy_input_kwh_kg_h2_base': (
                round(
                    cal_tea_range['estimates']['base']['energy_input_kwh_kg_h2'], 3
                )
                if cal_tea_range
                else None
            ),
        }

    # 3. Phase 5 PEMFC Cathode & Stack TEA (Screening + QE GC/Pourbaix Confirmed)
    pemfc_records: list[dict[str, object]] = []
    for cdir in campaign_dirs:
        fcdir = cdir / 'fuel_cell'
        if not fcdir.is_dir():
            continue
        for jpath in sorted(fcdir.glob('pemfc_*.json')):
            try:
                rec = json.loads(jpath.read_text(encoding='utf-8'))
            except Exception:
                continue
            if isinstance(rec, dict) and 'peak_power_W_cm2' in rec:
                pemfc_records.append(rec)

    pgm_free_classes = {
        'SAC',
        'DAC',
        'MOF',
        'COF',
        'MetalFreeCarbon',
        'Perovskite',
        'Spinel',
        'MAXPhase',
        'MXene',
    }
    _low_pgm_classes = {
        'SolidCatalyst',
        'SAA',
        'HEA',
        'MetalHydride',
        'MoltenMetal',
    }
    pemfc_stack_summary: PEMFCStackSummaryRow = {'status': 'no_pemfc_records'}
    coupled_tea: dict[str, CoupledPathwayTEARow] = {}
    if pemfc_records:
        dedup_pemfc = {
            (str(r.get('cathode_catalyst')), str(r.get('membrane'))): r
            for r in pemfc_records
        }
        ranked_pemfc = sorted(
            dedup_pemfc.values(), key=composite_score, reverse=True
        )
        best_cell = ranked_pemfc[0]
        mat_cls = str(best_cell.get('material_class', ''))
        default_pgm = 0.05 if mat_cls in pgm_free_classes else 0.10
        raw_pgm = best_cell.get('pgm_loading_mg_cm2')
        cell_pgm = (
            float(raw_pgm) if isinstance(raw_pgm, (int, float)) else default_pgm
        )
        raw_mem_cost = best_cell.get('membrane_cost_usd_cm2')
        cell_mem_cost = (
            float(raw_mem_cost)
            if isinstance(raw_mem_cost, (int, float))
            else 0.025
        )
        raw_pk_v = best_cell.get('peak_voltage_V')
        raw_pk_i = best_cell.get('peak_current_A_cm2')
        stack_res = model_stack(
            StackConfig(
                n_cells=stack_cells,
                cell_voltage_V=(
                    float(raw_pk_v) if isinstance(raw_pk_v, (int, float)) else 0.65
                ),
                current_density_A_cm2=(
                    float(raw_pk_i) if isinstance(raw_pk_i, (int, float)) else 1.5
                ),
                pgm_loading_mg_cm2=cell_pgm,
                membrane_cost_usd_cm2=cell_mem_cost,
            )
        )
        sys_eff = float(stack_res.get('system_efficiency', 0.50))
        best_cath_id = str(best_cell.get('cathode_catalyst', ''))
        best_cath_meta = candidate_Lookup.get(best_cath_id, {})
        raw_pk_p = best_cell.get('peak_power_W_cm2')
        raw_orr_eta = best_cell.get('orr_overpotential_V')
        raw_eff_pk = best_cell.get('efficiency_at_peak')
        pemfc_stack_summary = {
            'status': 'complete',
            'n_single_cell_simulations': len(dedup_pemfc),
            'champion_cathode': best_cath_id,
            'champion_cathode_literature_name': format_literature_catalyst_label(
                best_cath_id,
                material_class=str(
                    best_cell.get('material_class')
                    or best_cath_meta.get('material_class')
                    or ''
                ),
                genome=str(
                    best_cell.get('genome') or best_cath_meta.get('genome') or ''
                ),
            ),
            'champion_membrane': str(best_cell.get('membrane', '')),
            'peak_power_W_cm2': round(
                float(raw_pk_p) if isinstance(raw_pk_p, (int, float)) else 0.0, 4
            ),
            'orr_overpotential_V': round(
                float(raw_orr_eta)
                if isinstance(raw_orr_eta, (int, float))
                else 0.0,
                4,
            ),
            'cell_efficiency_at_peak': round(
                float(raw_eff_pk) if isinstance(raw_eff_pk, (int, float)) else 0.0,
                4,
            ),
            'stack_cells': int(stack_res.get('n_cells', stack_cells)),
            'stack_net_power_kW': round(float(stack_res.get('net_power_kW', 0.0)), 2),
            'stack_system_efficiency': round(sys_eff, 4),
            'stack_gravimetric_W_kg': round(
                float(stack_res.get('gravimetric_W_kg', 0.0)), 1
            ),
            'stack_volumetric_W_L': round(
                float(stack_res.get('volumetric_W_L', 0.0)), 1
            ),
            'stack_cost_usd_per_kW': round(
                float(stack_res.get('cost_per_kW', 0.0)), 2
            ),
        }
        for mode_name, entry in pathway_leaderboard.items():
            if entry.get('status') != 'complete':
                continue
            base_h2 = entry.get('tea_base_h2_usd_kg')
            opt_h2 = entry.get('tea_optimistic_h2_usd_kg')
            if base_h2 is not None and opt_h2 is not None and sys_eff > 0:
                coupled_tea[mode_name] = {
                    'h2_cost_base_usd_kg': base_h2,
                    'h2_cost_optimistic_usd_kg': opt_h2,
                    'pemfc_fuel_cost_base_usd_kwh_e': round(
                        float(base_h2) / (LHV_H2_KWH_PER_KG * sys_eff), 4
                    ),
                    'pemfc_fuel_cost_optimistic_usd_kwh_e': round(
                        float(opt_h2) / (LHV_H2_KWH_PER_KG * sys_eff), 4
                    ),
                }

    # 4. QE-Confirmed Grand-Canonical + Pourbaix Acid-Stable ORR PEMFC Stack & Coupled TEA
    orr_dft_path = root / 'dft' / 'orr_confirmatory_dft_summary.json'
    dft_calibrated_pemfc_stack_summary: PEMFCStackSummaryRow = {
        'status': 'no_confirmatory_orr_dft'
    }
    dft_calibrated_coupled_tea: dict[str, CoupledPathwayTEARow] = {}
    if orr_dft_path.is_file():
        try:
            orr_doc = json.loads(orr_dft_path.read_text(encoding='utf-8'))
            orr_records = (
                orr_doc.get('records', []) if isinstance(orr_doc, dict) else []
            )
            stable_converged = [
                r
                for r in orr_records
                if isinstance(r, dict)
                and bool(r.get('qe_converged', False))
                and bool(r.get('pourbaix_acid_stable', False))
                and isinstance(r.get('orr_overpotential_gc_V'), (int, float))
            ]
            simulated_gc_cells: list[dict[str, object]] = []
            for r in stable_converged:
                cid = str(r.get('candidate_id', ''))
                m_cls = str(r.get('material_class', 'SolidCatalyst'))
                eta_gc = float(r['orr_overpotential_gc_V'])
                tafel = float(TAFEL_SLOPE_BY_CLASS.get(m_cls, 70.0))
                cath_loading = 0.0 if m_cls in pgm_free_classes else 0.05
                cell_sim: dict[str, object] = dict(
                    simulate_pemfc(
                        PEMFCConfig(
                            T_K=353.15,
                            membrane_name='AquivionE87',
                            membrane_thickness_m=30e-6,
                            membrane_conductivity_S_m=14.0,
                            cathode_catalyst=cid,
                            cathode_loading_mg_cm2=cath_loading,
                            anode_loading_mg_cm2=0.05,
                            orr_overpotential_V=eta_gc,
                            orr_tafel_slope_mV_dec=tafel,
                            i_max_A_cm2=4.0,
                        )
                    )
                )
                cell_sim['material_class'] = m_cls
                cell_sim['genome'] = str(r.get('genome', ''))
                cell_sim['orr_overpotential_V'] = eta_gc
                cell_sim['orr_overpotential_che_V'] = round(
                    float(r.get('predicted_overpotential_V', eta_gc)), 4
                )
                cell_sim['orr_overpotential_gc_V'] = round(eta_gc, 4)
                cell_sim['orr_tafel_slope_mV_dec'] = tafel
                cell_sim['pourbaix_dissolution_margin_V'] = round(
                    float(r.get('pourbaix_dissolution_margin_V', 0.0)), 4
                )
                cell_sim['qe_dft_energy_Ry'] = r.get('qe_dft_energy_Ry')
                cell_sim['membrane_cost_usd_cm2'] = 0.025
                simulated_gc_cells.append(cell_sim)

            if simulated_gc_cells:
                ranked_gc = sorted(
                    simulated_gc_cells, key=composite_score, reverse=True
                )
                best_gc = ranked_gc[0]
                pgm_free_gc = [
                    c
                    for c in ranked_gc
                    if str(c.get('material_class', '')) in pgm_free_classes
                ]
                low_pgm_gc = [
                    c
                    for c in ranked_gc
                    if str(c.get('material_class', '')) in _low_pgm_classes
                    or str(c.get('material_class', '')) not in pgm_free_classes
                ]
                best_pgm_free = pgm_free_gc[0] if pgm_free_gc else None
                best_low_pgm = low_pgm_gc[0] if low_pgm_gc else None
                raw_gc_v = best_gc.get('peak_voltage_V')
                raw_gc_i = best_gc.get('peak_current_A_cm2')
                raw_gc_pgm = best_gc.get('pgm_loading_mg_cm2')
                raw_gc_mem = best_gc.get('membrane_cost_usd_cm2')
                gc_stack = model_stack(
                    StackConfig(
                        n_cells=stack_cells,
                        cell_voltage_V=(
                            float(raw_gc_v)
                            if isinstance(raw_gc_v, (int, float))
                            else 0.65
                        ),
                        current_density_A_cm2=(
                            float(raw_gc_i)
                            if isinstance(raw_gc_i, (int, float))
                            else 1.5
                        ),
                        pgm_loading_mg_cm2=(
                            float(raw_gc_pgm)
                            if isinstance(raw_gc_pgm, (int, float))
                            else 0.10
                        ),
                        membrane_cost_usd_cm2=(
                            float(raw_gc_mem)
                            if isinstance(raw_gc_mem, (int, float))
                            else 0.025
                        ),
                    )
                )
                gc_sys_eff = float(gc_stack.get('system_efficiency', 0.50))
                gc_cath_id = str(best_gc.get('cathode_catalyst', ''))
                raw_gc_pk_p = best_gc.get('peak_power_W_cm2')
                raw_gc_rt_p = best_gc.get('rated_power_W_cm2')
                raw_gc_eff = best_gc.get('efficiency_at_peak')
                raw_gc_ry = best_gc.get('qe_dft_energy_Ry')
                dft_calibrated_pemfc_stack_summary = {
                    'status': 'complete',
                    'n_qe_converged_acid_stable_cathodes': len(simulated_gc_cells),
                    'champion_cathode': gc_cath_id,
                    'champion_cathode_literature_name': format_literature_catalyst_label(
                        gc_cath_id,
                        material_class=str(best_gc.get('material_class', '')),
                        genome=str(best_gc.get('genome', '')),
                    ),
                    'champion_material_class': str(best_gc.get('material_class', '')),
                    'champion_genome': str(best_gc.get('genome', '')),
                    'champion_membrane': str(best_gc.get('membrane', 'AquivionE87')),
                    'orr_overpotential_che_V': float(
                        cast(float, best_gc['orr_overpotential_che_V'])
                    ),
                    'orr_overpotential_gc_V': float(
                        cast(float, best_gc['orr_overpotential_gc_V'])
                    ),
                    'orr_tafel_slope_mV_dec': float(
                        cast(float, best_gc['orr_tafel_slope_mV_dec'])
                    ),
                    'pourbaix_dissolution_margin_V': float(
                        cast(float, best_gc['pourbaix_dissolution_margin_V'])
                    ),
                    'qe_dft_energy_Ry': (
                        float(raw_gc_ry)
                        if isinstance(raw_gc_ry, (int, float))
                        else None
                    ),
                    'peak_power_W_cm2': round(
                        float(raw_gc_pk_p)
                        if isinstance(raw_gc_pk_p, (int, float))
                        else 0.0,
                        4,
                    ),
                    'rated_power_at_0p6V_W_cm2': (
                        round(float(raw_gc_rt_p), 4)
                        if isinstance(raw_gc_rt_p, (int, float))
                        else None
                    ),
                    'cell_efficiency_at_peak': round(
                        float(raw_gc_eff)
                        if isinstance(raw_gc_eff, (int, float))
                        else 0.0,
                        4,
                    ),
                    'stack_cells': int(gc_stack.get('n_cells', stack_cells)),
                    'stack_net_power_kW': round(
                        float(gc_stack.get('net_power_kW', 0.0)), 2
                    ),
                    'stack_system_efficiency': round(gc_sys_eff, 4),
                    'stack_gravimetric_W_kg': round(
                        float(gc_stack.get('gravimetric_W_kg', 0.0)), 1
                    ),
                    'stack_volumetric_W_L': round(
                        float(gc_stack.get('volumetric_W_L', 0.0)), 1
                    ),
                    'stack_cost_usd_per_kW': round(
                        float(gc_stack.get('cost_per_kW', 0.0)), 2
                    ),
                    'pgm_free_champion': (
                        {
                            'cathode_catalyst': str(
                                best_pgm_free.get('cathode_catalyst', '')
                            ),
                            'cathode_catalyst_literature_name': format_literature_catalyst_label(
                                str(best_pgm_free.get('cathode_catalyst', '')),
                                material_class=str(
                                    best_pgm_free.get('material_class', '')
                                ),
                                genome=str(best_pgm_free.get('genome', '')),
                            ),
                            'material_class': str(
                                best_pgm_free.get('material_class', '')
                            ),
                            'genome': str(best_pgm_free.get('genome', '')),
                            'orr_overpotential_gc_V': float(
                                cast(float, best_pgm_free['orr_overpotential_gc_V'])
                            ),
                            'pourbaix_dissolution_margin_V': float(
                                cast(
                                    float,
                                    best_pgm_free['pourbaix_dissolution_margin_V'],
                                )
                            ),
                            'peak_power_W_cm2': round(
                                float(
                                    cast(
                                        float,
                                        best_pgm_free.get('peak_power_W_cm2', 0.0),
                                    )
                                ),
                                4,
                            ),
                            'cell_efficiency_at_peak': round(
                                float(
                                    cast(
                                        float,
                                        best_pgm_free.get(
                                            'efficiency_at_peak', 0.0
                                        ),
                                    )
                                ),
                                4,
                            ),
                        }
                        if best_pgm_free is not None
                        else None
                    ),
                    'low_pgm_champion': (
                        {
                            'cathode_catalyst': str(
                                best_low_pgm.get('cathode_catalyst', '')
                            ),
                            'cathode_catalyst_literature_name': format_literature_catalyst_label(
                                str(best_low_pgm.get('cathode_catalyst', '')),
                                material_class=str(
                                    best_low_pgm.get('material_class', '')
                                ),
                                genome=str(best_low_pgm.get('genome', '')),
                            ),
                            'material_class': str(
                                best_low_pgm.get('material_class', '')
                            ),
                            'genome': str(best_low_pgm.get('genome', '')),
                            'orr_overpotential_gc_V': float(
                                cast(float, best_low_pgm['orr_overpotential_gc_V'])
                            ),
                            'pourbaix_dissolution_margin_V': float(
                                cast(
                                    float,
                                    best_low_pgm['pourbaix_dissolution_margin_V'],
                                )
                            ),
                            'peak_power_W_cm2': round(
                                float(
                                    cast(
                                        float,
                                        best_low_pgm.get('peak_power_W_cm2', 0.0),
                                    )
                                ),
                                4,
                            ),
                            'cell_efficiency_at_peak': round(
                                float(
                                    cast(
                                        float,
                                        best_low_pgm.get(
                                            'efficiency_at_peak', 0.0
                                        ),
                                    )
                                ),
                                4,
                            ),
                        }
                        if best_low_pgm is not None
                        else None
                    ),
                }
        except Exception:
            pass

    eff_for_cal_tea = float(
        dft_calibrated_pemfc_stack_summary.get('stack_system_efficiency')
        or pemfc_stack_summary.get('stack_system_efficiency')
        or 0.0
    )
    if eff_for_cal_tea > 0.0:
        for mode_name, cal_entry in dft_calibrated_pathway_leaderboard.items():
            if cal_entry.get('status') != 'complete':
                continue
            base_h2_cal = cal_entry.get('dft_calibrated_tea_base_h2_usd_kg')
            opt_h2_cal = cal_entry.get('dft_calibrated_tea_optimistic_h2_usd_kg')
            cons_h2_cal = cal_entry.get('dft_calibrated_tea_conservative_h2_usd_kg')
            if base_h2_cal is not None and opt_h2_cal is not None:
                dft_calibrated_coupled_tea[mode_name] = {
                    'dft_calibrated_h2_cost_base_usd_kg': base_h2_cal,
                    'dft_calibrated_h2_cost_optimistic_usd_kg': opt_h2_cal,
                    'dft_calibrated_h2_cost_conservative_usd_kg': cons_h2_cal,
                    'dft_calibrated_pemfc_fuel_cost_base_usd_kwh_e': round(
                        float(base_h2_cal) / (LHV_H2_KWH_PER_KG * eff_for_cal_tea),
                        4,
                    ),
                    'dft_calibrated_pemfc_fuel_cost_optimistic_usd_kwh_e': round(
                        float(opt_h2_cal) / (LHV_H2_KWH_PER_KG * eff_for_cal_tea),
                        4,
                    ),
                }

    supported_ni_closure_summary: SupportedNiFilamentClosureSummary | None = None
    try:
        supported_ni_closure_summary = evaluate_supported_ni_filament_closure()
    except Exception:
        pass

    return {
        'discovery_champions': discovery_champions,
        'pathway_leaderboard': pathway_leaderboard,
        'dft_calibrated_pathway_leaderboard': dft_calibrated_pathway_leaderboard,
        'pbr_ergun_residence_diagnostics': evaluate_pbr_and_ergun_residence_diagnostics(
            reactor_records
        ),
        'supported_ni_filament_closure_summary': supported_ni_closure_summary,
        'validation_gate_registry': dict(LITERATURE_VALIDATION_GATES),
        'literature_validation_gates': dict(LITERATURE_VALIDATION_GATES),
        'named_reference_catalysts': dict(NAMED_REFERENCE_CATALYSTS),
        'pemfc_stack_summary': pemfc_stack_summary,
        'dft_calibrated_pemfc_stack_summary': dft_calibrated_pemfc_stack_summary,
        'coupled_turquoise_pemfc_tea': coupled_tea,
        'dft_calibrated_coupled_turquoise_pemfc_tea': dft_calibrated_coupled_tea,
        'evidence_level': 'screening_multipathway_scorecard',
    }
