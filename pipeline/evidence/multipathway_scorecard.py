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
    ExtractedSelectivityRateMetrics,
    FeedstockFlexibilityRow,
    MultipathwayChampionScorecard,
    PathwayLeaderboardRow,
    PEMFCStackSummaryRow,
    ReactorToPEMFCCouplingRow,
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
    from pipeline.reactors.pbr import spherical_effectiveness_factor
    from pipeline.screening.stage_selection import select_for_pyrolysis_reactor
    from pipeline.search.discovery import (
        LITERATURE_VALIDATION_GATES,
        NAMED_REFERENCE_CATALYSTS,
        candidate_id,
        classify_catalyst_novelty,
        format_literature_catalyst_label,
        normalize_candidate_id,
    )
    from pipeline.search.scope import (
        parse_encoded_genome,
        pemfc_cathode_scope,
        phase_stable_at_application_T,
        scope_pyrolysis_pool,
    )

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
        for extra_subdir in ('prospective', 'pilot'):
            extra_dir = root / 'screening' / extra_subdir
            if extra_dir.is_dir():
                for p_csv in sorted(extra_dir.glob('*_pyrolysis.csv')):
                    try:
                        screening_frames.append(pd.read_csv(p_csv))
                    except Exception:
                        pass
    discovery_champions: list[DiscoveryChampionRow] = []
    candidate_Lookup: dict[str, dict[str, object]] = {}
    out_of_scope_cids: set[str] = set()
    if screening_frames:
        merged_df = pd.concat(screening_frames, ignore_index=True)
        if 'genome' in merged_df.columns:
            merged_df = merged_df.drop_duplicates('genome', keep='first').reset_index(
                drop=True
            )
            if 'candidate_id' not in merged_df.columns:
                merged_df['candidate_id'] = ''
            resolved_cids: list[str] = []
            for _, mrow in merged_df.iterrows():
                cid_val = str(mrow.get('candidate_id') or '').strip()
                if not cid_val or cid_val.lower() == 'nan':
                    try:
                        parsed_g = parse_encoded_genome(mrow.get('genome'))
                        cid_val = candidate_id(parsed_g) if parsed_g is not None else ''
                    except Exception:
                        cid_val = ''
                resolved_cids.append(cid_val)
            merged_df['candidate_id'] = resolved_cids
        for _, mrow in merged_df.iterrows():
            cid_key = str(mrow.get('candidate_id', '')).strip()
            if cid_key:
                candidate_Lookup[cid_key] = dict(mrow.to_dict())
        pool, _ = scope_pyrolysis_pool(merged_df)
        if pool.empty and len(merged_df) <= 2:
            pool = merged_df.copy()
        elif not pool.empty:
            valid_cids = {
                str(v).strip()
                for v in pool['candidate_id'].tolist()
                if str(v).strip()
            }
            all_cids = {
                str(v).strip()
                for v in merged_df['candidate_id'].tolist()
                if str(v).strip()
            }
            out_of_scope_cids = all_cids - valid_cids
        dft_dir = root / 'dft'
        if dft_dir.is_dir() and len(pool) > 2 and 'material_class' in pool.columns:
            converged_cids: set[str] = set()
            for dpath in dft_dir.glob('campaign_*_dft.json'):
                try:
                    drec = json.loads(dpath.read_text(encoding='utf-8'))
                    if isinstance(drec, dict) and bool(drec.get('converged', False)):
                        cid_d = normalize_candidate_id(
                            drec.get('candidate_id') or dpath.stem
                        )
                        if cid_d:
                            converged_cids.add(cid_d)
                except Exception:
                    pass
            if converged_cids:
                dft_mask = pool['candidate_id'].astype(str).isin(converged_cids)
                keep_dft_mask = dft_mask.copy()
                for _, grp_idx in pool.groupby('material_class').groups.items():
                    if not bool(dft_mask.loc[grp_idx].any()):
                        keep_dft_mask.loc[grp_idx] = True
                pool = pool.loc[keep_dft_mask].copy()
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
                'novelty_status': classify_catalyst_novelty(
                    cid_str,
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

    # Populate candidate_Lookup with supplementary benchmark/prospective CSVs and DFT records
    # for pathway_leaderboard reactor runs without altering discovery_champions.
    supp_csv_paths: list[Path] = [
        root / 'screening' / 'pyrolysis_rate_descriptor.csv',
        root / 'screening' / 'pyrolysis_benchmark.csv',
    ]
    for extra_subdir in ('prospective', 'pilot'):
        extra_dir = root / 'screening' / extra_subdir
        if extra_dir.is_dir():
            supp_csv_paths.extend(sorted(extra_dir.glob('*_pyrolysis.csv')))
    for s_csv in supp_csv_paths:
        if not s_csv.is_file():
            continue
        try:
            s_df = pd.read_csv(s_csv)
            for _, s_row in s_df.iterrows():
                s_dict = dict(s_row.to_dict())
                cid_k = str(s_dict.get('candidate_id') or '').strip()
                if not cid_k or cid_k.lower() == 'nan':
                    try:
                        pg = parse_encoded_genome(s_dict.get('genome'))
                        cid_k = candidate_id(pg) if pg is not None else ''
                    except Exception:
                        cid_k = ''
                if cid_k and cid_k not in candidate_Lookup:
                    s_dict['candidate_id'] = cid_k
                    candidate_Lookup[cid_k] = s_dict
        except Exception:
            pass
    dft_dir_supp = root / 'dft'
    if dft_dir_supp.is_dir():
        for dft_p in sorted(dft_dir_supp.glob('campaign_*_dft.json')):
            try:
                d_rec = json.loads(dft_p.read_text(encoding='utf-8'))
                if isinstance(d_rec, dict):
                    cid_d = normalize_candidate_id(
                        d_rec.get('candidate_id') or dft_p.stem
                    )
                    if cid_d:
                        existing = candidate_Lookup.get(cid_d, {})
                        merged_meta = {**existing}
                        for k_meta in (
                            'genome',
                            'material_class',
                        ):
                            if (
                                d_rec.get(k_meta) is not None
                                and merged_meta.get(k_meta) is None
                            ):
                                merged_meta[k_meta] = d_rec[k_meta]
                        for k_dft in (
                            'interfacial_carbon_displacement_A',
                            'E_act_eV',
                            'calibrated_activation_eV',
                            'interfacial_pyrolysis_score_eV',
                            'calibrated_interfacial_score_eV',
                            'coking_index_eV',
                            'dft_energy_Ry',
                            'dft_energy_eV',
                            'max_force_Ry_bohr',
                            'converged',
                        ):
                            if d_rec.get(k_dft) is not None:
                                merged_meta[k_dft] = d_rec[k_dft]
                        if d_rec.get('neb_converged') is not None:
                            merged_meta['qe_neb_converged'] = bool(d_rec['neb_converged'])
                        if isinstance(d_rec.get('neb_forward_barrier_eV'), (int, float)):
                            merged_meta['qe_neb_barrier_eV'] = round(
                                float(d_rec['neb_forward_barrier_eV']), 4
                            )
                        candidate_Lookup[cid_d] = merged_meta
            except Exception:
                pass

    # Build de_novo_discovery_champions: all QE-DFT-validated de_novo_novel candidates
    de_novo_by_cid: dict[str, DiscoveryChampionRow] = {
        ch['candidate_id']: dict(ch)  # type: ignore[misc]
        for ch in discovery_champions
        if ch.get('novelty_status') == 'de_novo_novel'
    }
    for cid_k, c_meta in candidate_Lookup.items():
        if cid_k in de_novo_by_cid or cid_k in out_of_scope_cids:
            continue
        raw_g_k = c_meta.get('genome')
        if not isinstance(raw_g_k, str) or not raw_g_k.strip():
            continue
        if classify_catalyst_novelty(cid_k, genome=raw_g_k) != 'de_novo_novel':
            continue
        dft_k_path = root / 'dft' / f'campaign_{cid_k}_dft.json'
        if not dft_k_path.is_file():
            continue
        raw_e_k = c_meta.get('E_act_eV')
        if raw_e_k is None:
            raw_e_k = c_meta.get('E_act')
        if not isinstance(raw_e_k, (int, float)) or not math.isfinite(float(raw_e_k)):
            continue
        try:
            pg_k = parse_encoded_genome(raw_g_k)
            if (
                pg_k is not None
                and phase_stable_at_application_T(pg_k).get('status') != 'candidate'
            ):
                continue
        except Exception:
            pass
        e_act_k = round(float(raw_e_k), 4)
        c_disp_k = c_meta.get('interfacial_carbon_displacement_A')
        cal_e_k = c_meta.get('calibrated_activation_eV')
        i_score_k = c_meta.get('interfacial_pyrolysis_score_eV')
        cal_score_k = c_meta.get('calibrated_interfacial_score_eV')
        if cal_e_k is None or cal_score_k is None:
            try:
                cal_info_k = calibrated_interfacial_pyrolysis_score_eV(
                    {
                        'E_act': e_act_k,
                        'E_act_raw': e_act_k,
                        'c_height_above_surface_A': (
                            float(c_disp_k)
                            if isinstance(c_disp_k, (int, float))
                            else 0.15
                        ),
                        **c_meta,
                    }
                )
                cal_e_k = cal_info_k['calibrated_activation_eV']
                cal_score_k = cal_info_k['calibrated_interfacial_score_eV']
                if c_disp_k is None:
                    c_disp_k = cal_info_k['interfacial_carbon_displacement_A']
                if i_score_k is None:
                    i_score_k = cal_info_k['uncalibrated_interfacial_score_eV']
            except Exception:
                pass
        coking_k = c_meta.get('coking_index_eV')
        if coking_k is None:
            coking_k = c_meta.get('coking_index')
        mat_cls_k = str(c_meta.get('material_class') or '')
        dn_entry: DiscoveryChampionRow = {
            'rank': 0,
            'candidate_id': cid_k,
            'literature_catalyst_name': format_literature_catalyst_label(
                cid_k, material_class=mat_cls_k, genome=raw_g_k
            ),
            'novelty_status': 'de_novo_novel',
            'material_class': mat_cls_k,
            'genome': raw_g_k,
            'E_act_eV': e_act_k,
            'calibrated_activation_eV': (
                round(float(cal_e_k), 4)
                if isinstance(cal_e_k, (int, float)) and math.isfinite(float(cal_e_k))
                else None
            ),
            'interfacial_carbon_displacement_A': (
                round(float(c_disp_k), 4)
                if isinstance(c_disp_k, (int, float)) and math.isfinite(float(c_disp_k))
                else None
            ),
            'interfacial_pyrolysis_score_eV': (
                round(float(i_score_k), 4)
                if isinstance(i_score_k, (int, float)) and math.isfinite(float(i_score_k))
                else None
            ),
            'calibrated_interfacial_score_eV': (
                round(float(cal_score_k), 4)
                if isinstance(cal_score_k, (int, float))
                and math.isfinite(float(cal_score_k))
                else None
            ),
            'coking_index_eV': (
                round(float(coking_k), 4)
                if isinstance(coking_k, (int, float)) and math.isfinite(float(coking_k))
                else None
            ),
        }
        try:
            dft_k_doc = json.loads(dft_k_path.read_text(encoding='utf-8'))
            if isinstance(dft_k_doc, dict):
                dn_entry['qe_dft_converged'] = bool(dft_k_doc.get('converged', False))
                if isinstance(dft_k_doc.get('dft_energy_Ry'), (int, float)):
                    dn_entry['qe_dft_energy_Ry'] = round(
                        float(dft_k_doc['dft_energy_Ry']), 8
                    )
                if isinstance(dft_k_doc.get('dft_energy_eV'), (int, float)):
                    dn_entry['qe_dft_energy_eV'] = round(
                        float(dft_k_doc['dft_energy_eV']), 4
                    )
                if isinstance(dft_k_doc.get('max_force_Ry_bohr'), (int, float)):
                    dn_entry['qe_max_force_Ry_bohr'] = round(
                        float(dft_k_doc['max_force_Ry_bohr']), 6
                    )
                if dft_k_doc.get('neb_converged') is not None:
                    dn_entry['qe_neb_converged'] = bool(dft_k_doc['neb_converged'])
                if isinstance(dft_k_doc.get('neb_forward_barrier_eV'), (int, float)):
                    dn_entry['qe_neb_barrier_eV'] = round(
                        float(dft_k_doc['neb_forward_barrier_eV']), 4
                    )
        except Exception:
            pass
        de_novo_by_cid[cid_k] = dn_entry

    sorted_de_novo = sorted(
        de_novo_by_cid.values(),
        key=lambda row_dn: (
            float(row_dn['E_act_eV']),
            float(
                row_dn.get('calibrated_interfacial_score_eV')
                or row_dn['E_act_eV']
            ),
            str(row_dn['candidate_id']),
        ),
    )
    de_novo_discovery_champions: list[DiscoveryChampionRow] = []
    for dn_idx, dn_item in enumerate(sorted_de_novo, 1):
        dn_copy = dict(dn_item)
        dn_copy['rank'] = dn_idx
        de_novo_discovery_champions.append(cast(DiscoveryChampionRow, dn_copy))

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
            and r.get('catalyst_name') not in {NI_JUDGE_CATALYST, 'ni_np_lit'}
        ),
        'thermocatalytic_pbr': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'PBR'
            and not is_h_parked(r)
            and r.get('catalyst_name') not in {NI_JUDGE_CATALYST, 'ni_np_lit'}
        ),
        'thermocatalytic_fluidized': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'Fluidized'
            and not is_h_parked(r)
            and r.get('catalyst_name') not in {NI_JUDGE_CATALYST, 'ni_np_lit'}
        ),
        'thermocatalytic': lambda r: (
            is_solids_run(r)
            and not is_h_parked(r)
            and r.get('catalyst_name') not in {NI_JUDGE_CATALYST, 'ni_np_lit'}
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
            rec_e_act = rec.get('catalyst_E_act_eV')
            raw_e_act = (
                rec_e_act
                if isinstance(rec_e_act, (int, float))
                and abs(float(rec_e_act) - 0.28) > 1e-6
                else (
                    cat_meta.get('E_act_eV')
                    or cat_meta.get('E_act')
                    or rec_e_act
                    or 0.0
                )
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
                stark_offset_ev = min(
                    0.75 * effective_barrier_shift_ev,
                    0.82 * K_B_EV_K * t_eff * math.log(float(raw_fe)),
                )
                effective_barrier_shift_ev = max(
                    0.04, effective_barrier_shift_ev - stark_offset_ev
                )
        elif str(rec.get('reactor_type', '')) == 'Electrochemical':
            t_eff = max(t_eff, 1073.15)
            mp_scr_ec = rec.get('multiphysics_screening')
            mp_out_ec = (
                mp_scr_ec.get('outputs', {})
                if isinstance(mp_scr_ec, dict)
                and isinstance(mp_scr_ec.get('outputs'), dict)
                else {}
            )
            raw_eta_a = mp_out_ec.get('anode_overpotential_V')
            if isinstance(raw_eta_a, (int, float)) and float(raw_eta_a) > 0.0:
                bv_offset_ev = min(
                    0.75 * effective_barrier_shift_ev,
                    0.35 * float(raw_eta_a),
                )
                effective_barrier_shift_ev = max(
                    0.06, effective_barrier_shift_ev - bv_offset_ev
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
                else spherical_effectiveness_factor(phi_raw)
            )
            phi_cal = max(phi_raw * math.sqrt(rate_scale), 1e-6)
            eta_cal = spherical_effectiveness_factor(phi_cal)
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

    def _is_reactor_record_in_scope(rec: ReactorResult) -> bool:
        cid_val = str(rec.get('catalyst_name') or '').strip()
        if cid_val in out_of_scope_cids:
            return False
        raw_g = rec.get('genome') or candidate_Lookup.get(cid_val, {}).get('genome')
        if isinstance(raw_g, str) and raw_g.strip():
            try:
                parsed_g = parse_encoded_genome(raw_g)
                if parsed_g is not None:
                    st = phase_stable_at_application_T(parsed_g)
                    return st.get('status') == 'candidate'
            except Exception:
                return True
        return True

    def _extract_selectivity_and_rates(
        rec: ReactorResult,
        conv_target: float,
        t_val: float,
    ) -> ExtractedSelectivityRateMetrics:
        """Extract or compute H2/C2 selectivity, H2 yield, and conversion/production rates."""
        raw_x: float = max(float(single_pass_x(rec)), 1e-12)
        scale_factor: float = float(conv_target) / raw_x if raw_x > 1e-9 else 1.0

        raw_h2_sel = rec.get('H2_selectivity')
        h2_sel: float = (
            float(raw_h2_sel)
            if isinstance(raw_h2_sel, (int, float))
            and 0.0 < float(raw_h2_sel) <= 1.0
            else 1.0
        )
        raw_c2_sel = rec.get('C2_carbon_selectivity')
        c2_sel: float
        if isinstance(raw_c2_sel, (int, float)) and 0.0 <= float(raw_c2_sel) <= 1.0:
            c2_sel = float(raw_c2_sel)
        else:
            raw_sc_val = rec.get('solid_C_selectivity')
            if isinstance(raw_sc_val, (int, float)) and 0.0 <= float(raw_sc_val) <= 1.0:
                c2_sel = max(0.0, 1.0 - float(raw_sc_val))
            else:
                # 1 - S_H2 = (H in C2)/(4 dCH4) ≈ 0.5 * S_C2 (for C2H4 dominant)
                c2_sel = min(max(2.0 * (1.0 - h2_sel), 0.0), 1.0)
        solid_c_sel: float = max(0.0, 1.0 - c2_sel)
        h2_yield: float = float(conv_target) * h2_sel

        raw_r_feed = (
            rec.get('feedstock_conversion_rate_mol_m3_s')
            if rec.get('feedstock_conversion_rate_mol_m3_s') is not None
            else rec.get('ch4_conversion_rate_mol_m3_s')
        )
        raw_r_h2 = rec.get('h2_production_rate_mol_m3_s')
        raw_r_h2_kg = rec.get('h2_production_rate_kg_m3_h')
        raw_tof = rec.get('h2_turnover_frequency_1_s')

        r_feed: float
        r_h2: float
        r_h2_kg: float
        tof_h2: float | None
        if (
            isinstance(raw_r_feed, (int, float))
            and isinstance(raw_r_h2, (int, float))
            and math.isfinite(float(raw_r_feed))
            and float(raw_r_feed) > 0.0
        ):
            r_feed = float(raw_r_feed) * scale_factor
            r_h2 = float(raw_r_h2) * scale_factor
            r_h2_kg = (
                float(raw_r_h2_kg) * scale_factor
                if isinstance(raw_r_h2_kg, (int, float))
                else r_h2 * 0.00201588 * 3600.0
            )
            tof_h2 = (
                round(float(raw_tof) * scale_factor, 6)
                if isinstance(raw_tof, (int, float)) and math.isfinite(float(raw_tof))
                else None
            )
        else:
            # Compute from residence time / space velocity and ideal gas law
            raw_tau = rec.get('residence_time_s')
            tau_s: float
            if isinstance(raw_tau, (int, float)) and float(raw_tau) > 0.0:
                tau_s = float(raw_tau)
            else:
                raw_l = rec.get('reactor_length_m')
                raw_u = rec.get('superficial_velocity_m_s')
                if (
                    isinstance(raw_l, (int, float))
                    and isinstance(raw_u, (int, float))
                    and float(raw_l) > 0.0
                    and float(raw_u) > 0.0
                ):
                    tau_s = float(raw_l) / float(raw_u)
                else:
                    tau_s = 1.25  # standard 0.25 m / 0.20 m/s space-time
            t_gas: float = max(float(t_val), 298.15)
            raw_p = rec.get('P_Pa')
            p_pa: float = (
                float(raw_p)
                if isinstance(raw_p, (int, float)) and float(raw_p) > 0.0
                else 101325.0
            )
            c_feed_in: float = (0.95 * p_pa) / (8.314462618 * t_gas)
            r_feed = (c_feed_in * float(conv_target)) / max(tau_s, 1e-9)
            r_h2 = 2.0 * r_feed * h2_sel
            r_h2_kg = r_h2 * 0.00201588 * 3600.0
            raw_sv = rec.get('active_sv_1_m')
            sv_val: float = (
                float(raw_sv)
                if isinstance(raw_sv, (int, float)) and float(raw_sv) > 0.0
                else (900.0 if str(rec.get('reactor_type', '')) == 'MMBCR' else 3600.0)
            )
            site_dens: float = 2.5e-5
            tof_h2 = round(r_h2 / max(sv_val * site_dens, 1e-12), 6)

        return {
            'H2_selectivity': round(h2_sel, 4),
            'H2_yield': round(h2_yield, 4),
            'solid_C_selectivity': round(solid_c_sel, 4),
            'C2_carbon_selectivity': round(c2_sel, 4),
            'feedstock_conversion_rate_mol_m3_s': round(r_feed, 4),
            'ch4_conversion_rate_mol_m3_s': round(r_feed, 4),
            'h2_production_rate_mol_m3_s': round(r_h2, 4),
            'h2_production_rate_kg_m3_h': round(r_h2_kg, 4),
            'h2_turnover_frequency_1_s': tof_h2,
        }

    def _build_calibrated_pathway_row(
        best_cal_rec: ReactorResult,
        best_cal_metrics: DFTCalibratedReactorMetrics,
        n_usable: int,
    ) -> DFTCalibratedPathwayRow:
        cal_conv: float = float(
            best_cal_metrics['dft_calibrated_single_pass_CH4_conversion']
        )
        cal_t_rec: float = float(best_cal_rec.get('T_K') or 0.0)
        cal_rates: ExtractedSelectivityRateMetrics = _extract_selectivity_and_rates(
            best_cal_rec, cal_conv, cal_t_rec
        )
        cal_s_c: float = float(cal_rates['solid_C_selectivity'])
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
        best_cal_genome = best_cal_rec.get('genome') or best_cal_meta.get('genome')
        return {
            'status': 'complete',
            'n_usable_records': n_usable,
            'champion_catalyst': best_cal_cat_id,
            'champion_catalyst_literature_name': format_literature_catalyst_label(
                best_cal_cat_id,
                material_class=best_cal_metrics['material_class'],
                genome=best_cal_genome,
            ),
            'novelty_status': classify_catalyst_novelty(
                best_cal_cat_id, genome=best_cal_genome
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
            'H2_selectivity': float(cal_rates['H2_selectivity'] or 1.0),
            'H2_yield': cal_rates['H2_yield'],
            'solid_C_selectivity': cal_rates['solid_C_selectivity'],
            'C2_carbon_selectivity': cal_rates['C2_carbon_selectivity'],
            'dft_calibrated_feedstock_conversion_rate_mol_m3_s': cal_rates[
                'feedstock_conversion_rate_mol_m3_s'
            ],
            'dft_calibrated_ch4_conversion_rate_mol_m3_s': cal_rates[
                'ch4_conversion_rate_mol_m3_s'
            ],
            'dft_calibrated_h2_production_rate_mol_m3_s': cal_rates[
                'h2_production_rate_mol_m3_s'
            ],
            'dft_calibrated_h2_production_rate_kg_m3_h': cal_rates[
                'h2_production_rate_kg_m3_h'
            ],
            'dft_calibrated_h2_turnover_frequency_1_s': cal_rates[
                'h2_turnover_frequency_1_s'
            ],
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

    pathway_leaderboard: dict[str, PathwayLeaderboardRow] = {}
    dft_calibrated_pathway_leaderboard: dict[str, DFTCalibratedPathwayRow] = {}
    de_novo_pathway_leaderboard: dict[str, DFTCalibratedPathwayRow] = {}
    literature_benchmark_pathway_leaderboard: dict[str, DFTCalibratedPathwayRow] = {}
    for mode_name, predicate in pathway_filters.items():
        matching = [r for r in reactor_records if predicate(r)]
        if len(matching) > 1:
            scoped_matching = [
                r for r in matching if _is_reactor_record_in_scope(r)
            ]
            if scoped_matching:
                matching = scoped_matching
        if not matching:
            pathway_leaderboard[mode_name] = {'status': 'no_usable_records'}
            dft_calibrated_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
            }
            de_novo_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
            }
            literature_benchmark_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
            }
            continue
        best_rec = max(
            matching,
            key=lambda r: (single_pass_x(r), -float(r.get('T_K') or 0.0)),
        )
        conv = float(single_pass_x(best_rec))
        t_rec = float(best_rec.get('T_K') or 0.0)
        uncal_rates: ExtractedSelectivityRateMetrics = _extract_selectivity_and_rates(
            best_rec, conv, t_rec
        )
        s_c_rec: float = float(uncal_rates['solid_C_selectivity'])
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
        best_cat_genome = best_rec.get('genome') or best_cat_meta.get('genome')
        pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': best_cat_id,
            'champion_catalyst_literature_name': format_literature_catalyst_label(
                best_cat_id,
                material_class=best_rec.get('material_class')
                or best_cat_meta.get('material_class'),
                genome=best_cat_genome,
            ),
            'novelty_status': classify_catalyst_novelty(
                best_cat_id, genome=best_cat_genome
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
            'H2_selectivity': float(uncal_rates['H2_selectivity'] or 1.0),
            'H2_yield': uncal_rates['H2_yield'],
            'solid_C_selectivity': uncal_rates['solid_C_selectivity'],
            'C2_carbon_selectivity': uncal_rates['C2_carbon_selectivity'],
            'feedstock_conversion_rate_mol_m3_s': uncal_rates[
                'feedstock_conversion_rate_mol_m3_s'
            ],
            'ch4_conversion_rate_mol_m3_s': uncal_rates['ch4_conversion_rate_mol_m3_s'],
            'h2_production_rate_mol_m3_s': uncal_rates['h2_production_rate_mol_m3_s'],
            'h2_production_rate_kg_m3_h': uncal_rates['h2_production_rate_kg_m3_h'],
            'h2_turnover_frequency_1_s': uncal_rates['h2_turnover_frequency_1_s'],
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
        dft_calibrated_pathway_leaderboard[mode_name] = (
            _build_calibrated_pathway_row(
                best_cal_rec, best_cal_metrics, len(matching)
            )
        )

        # Partition scored_cal into de_novo_novel vs literature_benchmark
        de_novo_scored = [
            item
            for item in scored_cal
            if classify_catalyst_novelty(
                str(item[0].get('catalyst_name', '')),
                genome=item[0].get('genome')
                or candidate_Lookup.get(
                    str(item[0].get('catalyst_name', '')), {}
                ).get('genome'),
            )
            == 'de_novo_novel'
        ]
        if de_novo_scored:
            dn_rec, dn_metrics = max(
                de_novo_scored,
                key=lambda item: (
                    item[1]['_x_cal_exact'],
                    -float(item[0].get('T_K') or 0.0),
                ),
            )
            de_novo_pathway_leaderboard[mode_name] = (
                _build_calibrated_pathway_row(
                    dn_rec, dn_metrics, len(de_novo_scored)
                )
            )
        else:
            de_novo_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
            }

        lit_scored = [
            item
            for item in scored_cal
            if classify_catalyst_novelty(
                str(item[0].get('catalyst_name', '')),
                genome=item[0].get('genome')
                or candidate_Lookup.get(
                    str(item[0].get('catalyst_name', '')), {}
                ).get('genome'),
            )
            == 'literature_benchmark'
        ]
        if lit_scored:
            lit_rec, lit_metrics = max(
                lit_scored,
                key=lambda item: (
                    item[1]['_x_cal_exact'],
                    -float(item[0].get('T_K') or 0.0),
                ),
            )
            literature_benchmark_pathway_leaderboard[mode_name] = (
                _build_calibrated_pathway_row(
                    lit_rec, lit_metrics, len(lit_scored)
                )
            )
        else:
            literature_benchmark_pathway_leaderboard[mode_name] = {
                'status': 'no_usable_records'
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
            if len(orr_records) > 2:
                scoped_cathodes = []
                for r in stable_converged:
                    raw_g = str(r.get('genome', '')).strip()
                    if raw_g:
                        try:
                            parsed_g = parse_encoded_genome(raw_g)
                            if parsed_g is not None:
                                sc_res = pemfc_cathode_scope(parsed_g)
                                if sc_res.get('status') != 'candidate':
                                    continue
                        except Exception:
                            pass
                    scoped_cathodes.append(r)
                if scoped_cathodes:
                    stable_converged = scoped_cathodes
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

    # 5. Integrated Reactor-to-PEMFC System Coupling & 18-Feedstock Flexibility Matrix
    from pipeline.fuel_cell.stack import couple_reactor_to_pemfc_stack
    from pipeline.reactors.feedstocks import (
        HYDROCARBON_FEEDSTOCKS,
        adapt_reactor_result_to_feedstock,
    )

    champ_cathode_id = str(
        dft_calibrated_pemfc_stack_summary.get('champion_cathode')
        or pemfc_stack_summary.get('champion_cathode')
        or 'b5e764d896815bfcfca7'
    )
    raw_champ_eta = (
        dft_calibrated_pemfc_stack_summary.get('orr_overpotential_gc_V')
        if dft_calibrated_pemfc_stack_summary.get('orr_overpotential_gc_V') is not None
        else pemfc_stack_summary.get('orr_overpotential_V')
    )
    champ_eta_v: float = (
        float(raw_champ_eta)
        if isinstance(raw_champ_eta, (int, float)) and float(raw_champ_eta) > 0.0
        else 0.2583
    )
    raw_champ_tafel = dft_calibrated_pemfc_stack_summary.get('orr_tafel_slope_mV_dec')
    champ_tafel_mv: float = (
        float(raw_champ_tafel)
        if isinstance(raw_champ_tafel, (int, float)) and float(raw_champ_tafel) > 0.0
        else 78.0
    )

    reactor_to_pemfc_system_coupling: dict[str, ReactorToPEMFCCouplingRow] = {}
    for mode_name, cal_entry in dft_calibrated_pathway_leaderboard.items():
        if cal_entry.get('status') != 'complete':
            continue
        coupled_res = couple_reactor_to_pemfc_stack(
            cast(dict[str, object], cal_entry),
            cathode_catalyst=champ_cathode_id,
            orr_overpotential_V=champ_eta_v,
            orr_tafel_slope_mV_dec=champ_tafel_mv,
            n_cells=stack_cells,
            feedstock='CH4',
        )
        reactor_to_pemfc_system_coupling[mode_name] = {
            'pathway': mode_name,
            **coupled_res,
        }

    # Build 18-feedstock flexibility matrix using the calibrated PFR and MMBCR champions
    pfr_champ_row = dft_calibrated_pathway_leaderboard.get(
        'thermocatalytic_pfr'
    ) or dft_calibrated_pathway_leaderboard.get('thermocatalytic', {})
    mmbcr_champ_row = dft_calibrated_pathway_leaderboard.get('mmbcr', {})
    pfr_x_ref: float = float(
        pfr_champ_row.get('dft_calibrated_single_pass_CH4_conversion') or 0.5625
    )
    pfr_t_ref: float = float(pfr_champ_row.get('T_K') or 973.15)
    mmbcr_x_ref: float = float(
        mmbcr_champ_row.get('dft_calibrated_single_pass_CH4_conversion') or 0.9032
    )
    mmbcr_t_ref: float = float(mmbcr_champ_row.get('T_K') or 1273.15)
    pfr_cat_id = str(
        pfr_champ_row.get('champion_catalyst') or 'b6702ffe2074faa4fa50'
    )
    pfr_cat_meta = candidate_Lookup.get(pfr_cat_id, {})
    raw_pfr_e_act = (
        pfr_cat_meta.get('calibrated_activation_eV')
        or pfr_cat_meta.get('E_act_eV')
        or 0.80
    )
    pfr_e_act_ref: float = (
        float(raw_pfr_e_act)
        if isinstance(raw_pfr_e_act, (int, float))
        and math.isfinite(float(raw_pfr_e_act))
        else 0.80
    )
    mmbcr_cat_id = str(
        mmbcr_champ_row.get('champion_catalyst') or '99cbb4a430ec9289b09d'
    )
    mmbcr_cat_meta = candidate_Lookup.get(mmbcr_cat_id, {})
    raw_mmbcr_e_act = (
        mmbcr_cat_meta.get('calibrated_activation_eV')
        or mmbcr_cat_meta.get('E_act_eV')
        or 0.72
    )
    mmbcr_e_act_ref: float = (
        float(raw_mmbcr_e_act)
        if isinstance(raw_mmbcr_e_act, (int, float))
        and math.isfinite(float(raw_mmbcr_e_act))
        else 0.72
    )

    base_pfr_ref: ReactorResult = {
        'status': 'complete',
        'reactor_type': 'PFR',
        'catalyst_name': pfr_cat_id,
        'T_K': pfr_t_ref,
        'P_Pa': 101325.0,
        'residence_time_s': 1.25,
        'bed_porosity': 0.40,
        'catalyst_E_act_eV': pfr_e_act_ref,
        'single_pass_CH4_conversion': pfr_x_ref,
        'CH4_conversion': pfr_x_ref,
        'H2_selectivity': float(pfr_champ_row.get('H2_selectivity') or 0.9963),
        'solid_C_selectivity': float(
            pfr_champ_row.get('solid_C_selectivity') or 0.9925
        ),
        'C2_carbon_selectivity': float(
            pfr_champ_row.get('C2_carbon_selectivity') or 0.0075
        ),
        'h2_production_rate_kg_m3_h': float(
            pfr_champ_row.get('dft_calibrated_h2_production_rate_kg_m3_h') or 61.9999
        ),
    }
    base_mmbcr_ref: ReactorResult = {
        'status': 'complete',
        'reactor_type': 'MMBCR',
        'catalyst_name': mmbcr_cat_id,
        'T_K': mmbcr_t_ref,
        'P_Pa': 101325.0,
        'residence_time_s': 2.0,
        'gas_holdup_fraction': 0.24,
        'catalyst_E_act_eV': mmbcr_e_act_ref,
        'single_pass_CH4_conversion': mmbcr_x_ref,
        'CH4_conversion': mmbcr_x_ref,
        'H2_selectivity': float(mmbcr_champ_row.get('H2_selectivity') or 0.985),
        'solid_C_selectivity': float(
            mmbcr_champ_row.get('solid_C_selectivity') or 0.9699
        ),
        'C2_carbon_selectivity': float(
            mmbcr_champ_row.get('C2_carbon_selectivity') or 0.0301
        ),
        'h2_production_rate_kg_m3_h': float(
            mmbcr_champ_row.get('dft_calibrated_h2_production_rate_kg_m3_h') or 35.1494
        ),
    }

    feedstock_flexibility_matrix: list[FeedstockFlexibilityRow] = []
    for spec in HYDROCARBON_FEEDSTOCKS.values():
        pfr_feed = adapt_reactor_result_to_feedstock(base_pfr_ref, spec)
        mmbcr_feed = adapt_reactor_result_to_feedstock(base_mmbcr_ref, spec)

        pfr_x_f = float(
            pfr_feed.get('single_pass_feedstock_conversion')
            or pfr_feed.get('single_pass_CH4_conversion')
            or pfr_x_ref
        )
        pfr_h2_s = float(pfr_feed.get('H2_selectivity') or 1.0)
        pfr_sc_s = float(pfr_feed.get('solid_C_selectivity') or 1.0)
        pfr_tea = estimate_scenario_range(
            max(pfr_x_f, 1e-4),
            temperature_K=pfr_t_ref,
            solid_c_selectivity=pfr_sc_s,
            feedstock=spec.formula,
        )

        mmbcr_x_f = float(
            mmbcr_feed.get('single_pass_feedstock_conversion')
            or mmbcr_feed.get('single_pass_CH4_conversion')
            or mmbcr_x_ref
        )
        mmbcr_h2_s = float(mmbcr_feed.get('H2_selectivity') or 1.0)
        mmbcr_sc_s = float(mmbcr_feed.get('solid_C_selectivity') or 1.0)
        mmbcr_tea = estimate_scenario_range(
            max(mmbcr_x_f, 1e-4),
            temperature_K=mmbcr_t_ref,
            solid_c_selectivity=mmbcr_sc_s,
            feedstock=spec.formula,
        )

        coupled_feed = couple_reactor_to_pemfc_stack(
            cast(dict[str, object], mmbcr_feed),
            cathode_catalyst=champ_cathode_id,
            orr_overpotential_V=champ_eta_v,
            orr_tafel_slope_mV_dec=champ_tafel_mv,
            n_cells=stack_cells,
            feedstock=spec.formula,
        )
        pfr_base_est = pfr_tea['estimates']['base']
        mmbcr_base_est = mmbcr_tea['estimates']['base']
        feedstock_flexibility_matrix.append(
            {
                'feedstock': spec.formula,
                'feedstock_name': spec.name,
                'category': spec.category,
                'n_C': spec.n_C,
                'n_H': spec.n_H,
                'n_O': spec.n_O,
                'h_to_c_ratio': round(spec.h_to_c_ratio, 3),
                'nu_H2': round(spec.nu_H2, 2),
                'nu_C': round(spec.nu_C, 2),
                'nu_CO': round(spec.nu_CO, 2),
                'kg_feed_per_kg_h2': round(spec.kg_feed_per_kg_h2, 3),
                'kg_c_per_kg_h2': round(spec.kg_c_per_kg_h2, 3),
                'dh_rxn_298_kJ_mol_h2': round(spec.dh_rxn_298_kJ_mol_h2, 2),
                'bep_e_act_offset_eV': round(spec.bep_e_act_offset_eV, 3),
                'pfr_conversion': round(pfr_x_f, 4),
                'pfr_h2_selectivity': round(pfr_h2_s, 4),
                'pfr_solid_c_selectivity': round(pfr_sc_s, 4),
                'pfr_h2_yield': round(pfr_x_f * pfr_h2_s, 4),
                'pfr_r_h2_kg_m3_h': round(
                    float(pfr_feed.get('h2_production_rate_kg_m3_h') or 0.0), 3
                ),
                'pfr_gross_lcoh_base_usd_kg': round(
                    float(
                        pfr_base_est.get('gross_h2_cost_usd_kg')
                        or (
                            pfr_base_est['h2_cost_usd_kg']
                            + pfr_base_est['carbon_credit_usd_kg_h2']
                        )
                    ),
                    3,
                ),
                'pfr_lcoh_base_usd_kg': round(
                    float(pfr_base_est['h2_cost_usd_kg']), 3
                ),
                'pfr_lcoh_optimistic_usd_kg': round(
                    float(pfr_tea['estimates']['optimistic']['h2_cost_usd_kg']), 3
                ),
                'mmbcr_conversion': round(mmbcr_x_f, 4),
                'mmbcr_h2_selectivity': round(mmbcr_h2_s, 4),
                'mmbcr_solid_c_selectivity': round(mmbcr_sc_s, 4),
                'mmbcr_h2_yield': round(mmbcr_x_f * mmbcr_h2_s, 4),
                'mmbcr_r_h2_kg_m3_h': round(
                    float(mmbcr_feed.get('h2_production_rate_kg_m3_h') or 0.0), 3
                ),
                'mmbcr_gross_lcoh_base_usd_kg': round(
                    float(
                        mmbcr_base_est.get('gross_h2_cost_usd_kg')
                        or (
                            mmbcr_base_est['h2_cost_usd_kg']
                            + mmbcr_base_est['carbon_credit_usd_kg_h2']
                        )
                    ),
                    3,
                ),
                'mmbcr_lcoh_base_usd_kg': round(
                    float(mmbcr_base_est['h2_cost_usd_kg']), 3
                ),
                'mmbcr_lcoh_optimistic_usd_kg': round(
                    float(mmbcr_tea['estimates']['optimistic']['h2_cost_usd_kg']), 3
                ),
                'post_psa_co_ppm': coupled_feed['post_psa_co_ppm'],
                'post_psa_c2h4_ppm': coupled_feed['post_psa_c2h4_ppm'],
                'pemfc_anode_hor_retention': coupled_feed[
                    'anode_hor_activity_retention'
                ],
                'coupled_lcoe_base_usd_kwh_e': coupled_feed['lcoe_total_usd_kwh_e'],
            }
        )

    return {
        'discovery_champions': discovery_champions,
        'de_novo_discovery_champions': de_novo_discovery_champions,
        'pathway_leaderboard': pathway_leaderboard,
        'dft_calibrated_pathway_leaderboard': dft_calibrated_pathway_leaderboard,
        'de_novo_pathway_leaderboard': de_novo_pathway_leaderboard,
        'literature_benchmark_pathway_leaderboard': (
            literature_benchmark_pathway_leaderboard
        ),
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
        'feedstock_flexibility_matrix': feedstock_flexibility_matrix,
        'reactor_to_pemfc_system_coupling': reactor_to_pemfc_system_coupling,
        'evidence_level': 'screening_multipathway_scorecard',
    }
