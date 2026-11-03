"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math


@dataclass(frozen=True)
class TEAAssumptions:
    """Hold a traceable techno-economic scenario.

    Attributes:
        natural_gas_usd_mmbtu: Configured natural gas usd mmbtu value.
        electricity_usd_kwh: Configured electricity usd kwh value.
        base_energy_kwh_kg_h2: Configured base energy kwh kg h2 value.
        capex_usd_kg_h2: Configured capex usd kg h2 value.
        carbon_value_usd_kg_h2: Configured carbon value usd kg h2 value.
        source_id: Configured source id value.
    """

    natural_gas_usd_mmbtu: float
    electricity_usd_kwh: float
    base_energy_kwh_kg_h2: float
    capex_usd_kg_h2: float
    carbon_value_usd_kg_h2: float
    source_id: str


SCENARIOS = {
    'optimistic': TEAAssumptions(
        2.50, 0.035, 7.0, 0.30, 1.00, 'screening-tea-v3:optimistic'
    ),
    'base': TEAAssumptions(3.50, 0.060, 8.5, 0.50, 0.80, 'screening-tea-v3:base'),
    'conservative': TEAAssumptions(
        6.00, 0.120, 12.0, 1.20, 0.20, 'screening-tea-v3:conservative'
    ),
}


# CH4 -> C + 2 H2: 16.04 g CH4 per 4.032 g H2, i.e. 3.98 kg CH4 per kg H2.
# Priced on the HHV basis used for natural-gas $/MMBtu quotes
# (55.5 MJ/kg CH4; 1 MMBtu = 1055.06 MJ): ~0.209 MMBtu per kg H2 at X = 1.
# (v2 used 0.05 MMBtu/kg, about 4x too little feedstock.)
KG_CH4_PER_KG_H2 = 16.043 / (2 * 2.016)
MMBTU_CH4_PER_KG_H2 = KG_CH4_PER_KG_H2 * 55.5 / 1055.06
PSA_CH4_RECOVERY: float = 0.96  # Industrial pressure-swing adsorption (PSA) unreacted CH4 recovery fraction
RECYCLE_ENERGY_FRACTION: float = 0.25  # Fractional compression + reheat energy penalty per unit recycle ratio R = (1-X)/X
RECYCLE_CAPEX_SCALING: float = 0.12  # Fractional PSA + recycle compressor CAPEX increment per unit recycle ratio R


def estimate_hydrogen_cost(
    conversion: float,
    scenario: str = 'base',
    *,
    psa_ch4_recovery: float = PSA_CH4_RECOVERY,
    temperature_K: float | None = None,
    solid_c_selectivity: float = 1.0,
    specific_energy_kwh_kg_h2: float | None = None,
) -> dict:
    """Return an assumption-labelled screening estimate with PSA CH4 recycle, never a measurement.

    Models an industrial PSA H2 purification and unreacted CH4 recycle loop:
    - Recycle ratio ``R = (1 - X) / X`` moles recycled per mole converted.
    - Net methane feedstock utilization ``U_CH4 = X / (1 - r_PSA * (1 - X))``
      (equals 1.0 at ``X = 1.0``, and reflects only PSA tail-gas purge loss at
      partial single-pass conversion rather than venting all unreacted CH4).
    - Recycle compression and sensible-reheat energy penalty scales with ``R``
      as ``base_energy * (1 + RECYCLE_ENERGY_FRACTION * R)``.
    - Recycle loop and PSA vessel sizing scales CAPEX as
      ``capex * (1 + RECYCLE_CAPEX_SCALING * R)``.

    Args:
        conversion: Single-pass methane conversion fraction in (0, 1].
        scenario: Scenario used by this operation.
        psa_ch4_recovery: Fraction of unreacted CH4 recovered and recycled by PSA.
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    assumptions = SCENARIOS[scenario]
    conversion = float(conversion)
    r_psa = min(max(float(psa_ch4_recovery), 0.0), 0.999)
    s_c = min(max(float(solid_c_selectivity), 0.0), 1.0)
    recycle_ratio = (1.0 - conversion) / conversion
    net_ch4_utilization = conversion / max(1.0 - r_psa * (1.0 - conversion), 1e-6)
    base_energy = assumptions.base_energy_kwh_kg_h2
    if specific_energy_kwh_kg_h2 is not None and float(specific_energy_kwh_kg_h2) > 0.0:
        base_energy = float(specific_energy_kwh_kg_h2)
    elif temperature_K is not None and float(temperature_K) > 298.15:
        # Scale baseline thermal energy relative to a 1073.15 K (800 °C) reference reactor:
        # 65% endothermic reaction enthalpy ΔH_rxn(T) + 35% sensible feed heating (T - 298.15 K).
        try:
            from pipeline.reactors.carbon_accounting import (
                pyrolysis_reaction_enthalpy_kJ_mol,
            )

            dh_t = pyrolysis_reaction_enthalpy_kJ_mol(float(temperature_K))
            dh_ref = pyrolysis_reaction_enthalpy_kJ_mol(1073.15)
        except ImportError:
            # Linear heat-capacity approximation ΔCp ≈ 0.0155 kJ/(mol·K) above ΔH°_298 = 74.85 kJ/mol
            dh_t = 74.85 + 0.0155 * (float(temperature_K) - 298.15)
            dh_ref = 74.85 + 0.0155 * (1073.15 - 298.15)
        sensible_ratio = max(float(temperature_K) - 298.15, 0.0) / (1073.15 - 298.15)
        thermal_scale = 0.65 * (dh_t / max(dh_ref, 1.0)) + 0.35 * sensible_ratio
        base_energy = base_energy * max(thermal_scale, 0.5)
    energy = base_energy * (
        1.0 + RECYCLE_ENERGY_FRACTION * recycle_ratio
    )
    effective_capex = assumptions.capex_usd_kg_h2 * (
        1.0 + RECYCLE_CAPEX_SCALING * recycle_ratio
    )
    carbon_credit = assumptions.carbon_value_usd_kg_h2 * s_c
    cost = (
        assumptions.natural_gas_usd_mmbtu * MMBTU_CH4_PER_KG_H2 / net_ch4_utilization
        + energy * assumptions.electricity_usd_kwh
        + effective_capex
        - carbon_credit
    )
    return {
        'scenario': scenario,
        'h2_cost_usd_kg': round(cost, 3),
        'conversion': conversion,
        'recycle_ratio': round(recycle_ratio, 4),
        'net_ch4_utilization': round(net_ch4_utilization, 4),
        'psa_ch4_recovery': r_psa,
        'solid_c_selectivity': round(s_c, 4),
        'carbon_credit_usd_kg_h2': round(carbon_credit, 4),
        'energy_input_kwh_kg_h2': energy,
        'effective_capex_usd_kg_h2': round(effective_capex, 4),
        'assumptions': asdict(assumptions),
        'evidence_level': 'screening_scenario_not_measured_tea',
    }


def estimate_scenario_range(
    conversion: float,
    *,
    temperature_K: float | None = None,
    solid_c_selectivity: float = 1.0,
    specific_energy_kwh_kg_h2: float | None = None,
) -> dict:
    """Evaluate hydrogen cost across named techno-economic scenarios.

    Args:
        conversion: Methane conversion fraction used by the cost model.
        temperature_K: Optional reactor operating temperature in kelvin.
        solid_c_selectivity: Fraction of converted carbon captured as solid carbon in [0, 1].
        specific_energy_kwh_kg_h2: Optional pathway-specific energy input in kWh/kg H2.

    Returns:
        A dictionary containing estimate scenario range outputs, status, and supporting metadata.
    """
    estimates = {
        name: estimate_hydrogen_cost(
            conversion,
            name,
            temperature_K=temperature_K,
            solid_c_selectivity=solid_c_selectivity,
            specific_energy_kwh_kg_h2=specific_energy_kwh_kg_h2,
        )
        for name in SCENARIOS
    }
    return {
        'estimates': estimates,
        'min_usd_kg': min(x['h2_cost_usd_kg'] for x in estimates.values()),
        'max_usd_kg': max(x['h2_cost_usd_kg'] for x in estimates.values()),
        'evidence_level': 'screening_sensitivity_range',
    }


# Lower heating value of H2 (33.33 kWh/kg H2) for coupled H2 -> PEMFC LCOE.
LHV_H2_KWH_PER_KG = 33.33


def build_multipathway_champion_scorecard(
    results_root: str = 'results',
    *,
    stack_cells: int = 400,
    top_k_catalysts: int = 10,
) -> dict:
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
    import json
    from pathlib import Path
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
        evaluate_pbr_and_b4_diagnostics,
        is_h_parked,
        is_production_reactor_record,
        is_solids_run,
        single_pass_x,
    )
    from pipeline.screening.stage_selection import select_for_pyrolysis_reactor
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
    discovery_champions: list[dict] = []
    candidate_Lookup: dict[str, dict] = {}
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
                candidate_Lookup[cid_key] = prow.to_dict()
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
            champ_entry: dict = {
                'rank': rank_idx,
                'candidate_id': cid_str,
                'material_class': str(row.get('material_class', '')),
                'genome': str(row.get('genome', '')),
                'E_act_eV': round(float(row['E_act']), 4),
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
                except Exception:
                    pass
            discovery_champions.append(champ_entry)
            if cid_str:
                candidate_Lookup[cid_str] = {
                    **candidate_Lookup.get(cid_str, {}),
                    **champ_entry,
                }

    # 2. Phase 2 6-Pathway Reactor & TEA Leaderboard (Raw + DFT-Calibrated)
    reactor_records: list[dict] = []
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
                reactor_records.append(rec)

    pathway_filters = {
        'thermocatalytic_pfr': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'PFR'
            and not is_h_parked(r)
            and r.get('catalyst_name') != 'ni_np_lit'
        ),
        'thermocatalytic_pbr': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'PBR'
            and not is_h_parked(r)
            and r.get('catalyst_name') != 'ni_np_lit'
        ),
        'thermocatalytic_fluidized': lambda r: (
            is_solids_run(r)
            and r.get('reactor_type') == 'Fluidized'
            and not is_h_parked(r)
            and r.get('catalyst_name') != 'ni_np_lit'
        ),
        'thermocatalytic': lambda r: (
            is_solids_run(r)
            and not is_h_parked(r)
            and r.get('catalyst_name') != 'ni_np_lit'
        ),
        'mmbcr': lambda r: r.get('reactor_type') == 'MMBCR',
        'ntec': lambda r: r.get('reactor_type') == 'NTEC',
        'electrochemical': lambda r: r.get('reactor_type') == 'Electrochemical',
    }

    def _compute_dft_calibrated_reactor_metrics(rec: dict) -> dict:
        cat_name = str(rec.get('catalyst_name', '')).strip()
        cat_meta = candidate_Lookup.get(cat_name, {})
        mat_cls = str(
            rec.get('material_class')
            or cat_meta.get('material_class')
            or ('MoltenMetal' if rec.get('reactor_type') == 'MMBCR' else 'HEA')
        )
        t_val = float(rec.get('T_K') or 1100.0)
        t_eff = max(t_val, 300.0)
        # Class-level QE active-learning barrier correction (plus interfacial carbon penalty)
        default_bias = 0.228733 if mat_cls == 'MoltenMetal' else 0.243114
        delta_e_dft = float(al_turq_class_bias.get(mat_cls, default_bias))
        c_disp_val = cat_meta.get('interfacial_carbon_displacement_A')
        c_disp_a = (
            float(c_disp_val)
            if isinstance(c_disp_val, (int, float)) and math.isfinite(float(c_disp_val))
            else 0.0
        )
        # Effective barrier shift includes the QE AL class bias and 0.15 eV/A interfacial coking drag
        effective_barrier_shift_ev = max(0.0, delta_e_dft + 0.15 * c_disp_a)
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
            'dft_calibrated_single_pass_CH4_conversion': round(x_cal, 4),
            'dft_calibrated_overall_effectiveness_factor': cal_eta_out,
            'dft_calibrated_thiele_modulus': cal_phi_out,
        }

    pathway_leaderboard: dict[str, dict] = {}
    dft_calibrated_pathway_leaderboard: dict[str, dict] = {}
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
        pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': str(best_rec.get('catalyst_name', '')),
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
                item[1]['dft_calibrated_single_pass_CH4_conversion'],
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
        dft_calibrated_pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': str(best_cal_rec.get('catalyst_name', '')),
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
            'kinetics_status': 'qe_dft_calibrated_bep',
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
    pemfc_records: list[dict] = []
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
    pemfc_stack_summary: dict = {'status': 'no_pemfc_records'}
    coupled_tea: dict[str, dict] = {}
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
        cell_pgm = float(best_cell.get('pgm_loading_mg_cm2', default_pgm))
        cell_mem_cost = float(best_cell.get('membrane_cost_usd_cm2', 0.025))
        stack_res = model_stack(
            StackConfig(
                n_cells=stack_cells,
                cell_voltage_V=float(best_cell.get('peak_voltage_V', 0.65)),
                current_density_A_cm2=float(best_cell.get('peak_current_A_cm2', 1.5)),
                pgm_loading_mg_cm2=cell_pgm,
                membrane_cost_usd_cm2=cell_mem_cost,
            )
        )
        sys_eff = float(stack_res.get('system_efficiency', 0.50))
        pemfc_stack_summary = {
            'status': 'complete',
            'n_single_cell_simulations': len(dedup_pemfc),
            'champion_cathode': str(best_cell.get('cathode_catalyst', '')),
            'champion_membrane': str(best_cell.get('membrane', '')),
            'peak_power_W_cm2': round(float(best_cell.get('peak_power_W_cm2', 0.0)), 4),
            'orr_overpotential_V': round(
                float(best_cell.get('orr_overpotential_V', 0.0)), 4
            ),
            'cell_efficiency_at_peak': round(
                float(best_cell.get('efficiency_at_peak', 0.0)), 4
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
    dft_calibrated_pemfc_stack_summary: dict = {'status': 'no_confirmatory_orr_dft'}
    dft_calibrated_coupled_tea: dict[str, dict] = {}
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
            simulated_gc_cells: list[dict] = []
            for r in stable_converged:
                cid = str(r.get('candidate_id', ''))
                m_cls = str(r.get('material_class', 'SolidCatalyst'))
                eta_gc = float(r['orr_overpotential_gc_V'])
                tafel = float(TAFEL_SLOPE_BY_CLASS.get(m_cls, 70.0))
                cath_loading = 0.0 if m_cls in pgm_free_classes else 0.05
                cell_sim = dict(
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
                    if str(c.get('material_class', '')) not in pgm_free_classes
                ]
                best_pgm_free = pgm_free_gc[0] if pgm_free_gc else None
                best_low_pgm = low_pgm_gc[0] if low_pgm_gc else None
                gc_stack = model_stack(
                    StackConfig(
                        n_cells=stack_cells,
                        cell_voltage_V=float(best_gc.get('peak_voltage_V', 0.65)),
                        current_density_A_cm2=float(
                            best_gc.get('peak_current_A_cm2', 1.5)
                        ),
                        pgm_loading_mg_cm2=float(
                            best_gc.get('pgm_loading_mg_cm2', 0.10)
                        ),
                        membrane_cost_usd_cm2=float(
                            best_gc.get('membrane_cost_usd_cm2', 0.025)
                        ),
                    )
                )
                gc_sys_eff = float(gc_stack.get('system_efficiency', 0.50))
                dft_calibrated_pemfc_stack_summary = {
                    'status': 'complete',
                    'n_qe_converged_acid_stable_cathodes': len(simulated_gc_cells),
                    'champion_cathode': str(best_gc.get('cathode_catalyst', '')),
                    'champion_material_class': str(best_gc.get('material_class', '')),
                    'champion_genome': str(best_gc.get('genome', '')),
                    'champion_membrane': str(best_gc.get('membrane', 'AquivionE87')),
                    'orr_overpotential_che_V': best_gc['orr_overpotential_che_V'],
                    'orr_overpotential_gc_V': best_gc['orr_overpotential_gc_V'],
                    'orr_tafel_slope_mV_dec': best_gc['orr_tafel_slope_mV_dec'],
                    'pourbaix_dissolution_margin_V': best_gc[
                        'pourbaix_dissolution_margin_V'
                    ],
                    'qe_dft_energy_Ry': best_gc.get('qe_dft_energy_Ry'),
                    'peak_power_W_cm2': round(
                        float(best_gc.get('peak_power_W_cm2', 0.0)), 4
                    ),
                    'rated_power_at_0p6V_W_cm2': (
                        round(float(best_gc['rated_power_W_cm2']), 4)
                        if best_gc.get('rated_power_W_cm2') is not None
                        else None
                    ),
                    'cell_efficiency_at_peak': round(
                        float(best_gc.get('efficiency_at_peak', 0.0)), 4
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
                            'material_class': str(
                                best_pgm_free.get('material_class', '')
                            ),
                            'genome': str(best_pgm_free.get('genome', '')),
                            'orr_overpotential_gc_V': best_pgm_free[
                                'orr_overpotential_gc_V'
                            ],
                            'pourbaix_dissolution_margin_V': best_pgm_free[
                                'pourbaix_dissolution_margin_V'
                            ],
                            'peak_power_W_cm2': round(
                                float(best_pgm_free.get('peak_power_W_cm2', 0.0)), 4
                            ),
                            'cell_efficiency_at_peak': round(
                                float(best_pgm_free.get('efficiency_at_peak', 0.0)), 4
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
                            'material_class': str(
                                best_low_pgm.get('material_class', '')
                            ),
                            'genome': str(best_low_pgm.get('genome', '')),
                            'orr_overpotential_gc_V': best_low_pgm[
                                'orr_overpotential_gc_V'
                            ],
                            'pourbaix_dissolution_margin_V': best_low_pgm[
                                'pourbaix_dissolution_margin_V'
                            ],
                            'peak_power_W_cm2': round(
                                float(best_low_pgm.get('peak_power_W_cm2', 0.0)), 4
                            ),
                            'cell_efficiency_at_peak': round(
                                float(best_low_pgm.get('efficiency_at_peak', 0.0)), 4
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

    return {
        'discovery_champions': discovery_champions,
        'pathway_leaderboard': pathway_leaderboard,
        'dft_calibrated_pathway_leaderboard': dft_calibrated_pathway_leaderboard,
        'pbr_b4_diagnostics': evaluate_pbr_and_b4_diagnostics(reactor_records),
        'pemfc_stack_summary': pemfc_stack_summary,
        'dft_calibrated_pemfc_stack_summary': dft_calibrated_pemfc_stack_summary,
        'coupled_turquoise_pemfc_tea': coupled_tea,
        'dft_calibrated_coupled_turquoise_pemfc_tea': dft_calibrated_coupled_tea,
        'evidence_level': 'screening_multipathway_scorecard',
    }

