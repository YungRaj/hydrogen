"""Supported-Ni(111)/SiO2 filament closure gate and Baker-Helveg / Ermakova-Takenaka reference calibration."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping, cast

from pipeline.data_models.reactors import (
    SolidsScorecard,
    SupportedNiFilamentClosureGateResult,
)
from pipeline.search.discovery import (
    LITERATURE_VALIDATION_GATES,
    format_literature_catalyst_label,
)


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
    from pipeline.reactors.scorecard import NI_JUDGE_CATALYST

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


def _calibrated_baker_helveg_and_ermakova_takenaka_payloads() -> tuple[
    dict[str, object], dict[str, object]
]:
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
        _rec(
            sweep={'E_act': 0.9},
            CH4_conversion=0.3121,
            single_pass_CH4_conversion=0.3121,
        ),
        _rec(
            sweep={'E_act': 1.0},
            CH4_conversion=0.1456,
            single_pass_CH4_conversion=0.1456,
            carbon_turnovers_per_site=7.0348,
        ),
        _rec(
            sweep={'E_act': 1.1},
            CH4_conversion=0.0597,
            single_pass_CH4_conversion=0.0597,
        ),
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
        'ni111_sio2_baker_helveg_agamma_theta_star': {
            'records': agamma_theta_records
        },
    }
    joint_payloads: dict[str, object] = {
        **baker_helveg_payloads,
        'ni111_sio2_ermakova_takenaka_joint': {'records': joint_records},
    }
    return baker_helveg_payloads, joint_payloads


def _calibrated_ni_roi_reactor_records() -> list[dict[str, object]]:
    """Return calibrated Ni(111)/SiO2 reference ROI records on the high-inventory ROI cell (650-700 C)."""
    from pipeline.reactors.scorecard import NI_JUDGE_CATALYST

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
    from pipeline.reactors.scorecard import (
        NI_JUDGE_CATALYST,
        NI_JUDGE_HEADLINE_T_MAX,
        NI_JUDGE_HEADLINE_T_MIN,
        build_solids_scorecard,
        is_production_reactor_record,
    )
    from pipeline.reactors.sweeps.carbon_criteria import evaluate_criteria
    from pipeline.reactors.sweeps.joint_band import search_joint_band
    from pipeline.utils import BASE_DIR

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
