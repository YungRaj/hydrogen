"""Transparent scenario-level screening TEA for turquoise hydrogen."""

from __future__ import annotations

from dataclasses import asdict, dataclass


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

def estimate_hydrogen_cost(conversion: float, scenario: str = 'base') -> dict:
    """Return an assumption-labelled screening estimate, never a measurement.

    Args:
        conversion: Conversion used by this operation.
        scenario: Scenario used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    if scenario not in SCENARIOS:
        raise ValueError(f'unknown TEA scenario: {scenario}')
    if not 0 < float(conversion) <= 1:
        raise ValueError('conversion must be in (0, 1]')
    assumptions = SCENARIOS[scenario]
    conversion = float(conversion)
    energy = assumptions.base_energy_kwh_kg_h2 / conversion
    cost = (
        assumptions.natural_gas_usd_mmbtu * MMBTU_CH4_PER_KG_H2 / conversion
        + energy * assumptions.electricity_usd_kwh
        + assumptions.capex_usd_kg_h2
        - assumptions.carbon_value_usd_kg_h2
    )
    return {
        'scenario': scenario,
        'h2_cost_usd_kg': round(cost, 3),
        'conversion': conversion,
        'energy_input_kwh_kg_h2': energy,
        'assumptions': asdict(assumptions),
        'evidence_level': 'screening_scenario_not_measured_tea',
    }


def estimate_scenario_range(conversion: float) -> dict:
    """Evaluate hydrogen cost across named techno-economic scenarios.

    Args:
        conversion: Methane conversion fraction used by the cost model.

    Returns:
        A dictionary containing estimate scenario range outputs, status, and supporting metadata.
    """
    estimates = {name: estimate_hydrogen_cost(conversion, name) for name in SCENARIOS}
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
    from pipeline.fuel_cell import composite_score
    from pipeline.fuel_cell.stack import StackConfig, model_stack
    from pipeline.reactors.eligibility import is_usable_result
    from pipeline.reactors.scorecard import (
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

    # 1. Phase 1 Discovery Champions with Interfacial Carbon Score
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
    if screening_frames:
        merged_df = pd.concat(screening_frames, ignore_index=True)
        if 'genome' in merged_df.columns:
            merged_df = merged_df.drop_duplicates('genome', keep='first').reset_index(
                drop=True
            )
        pool, _ = scope_pyrolysis_pool(merged_df)
        selected = select_for_pyrolysis_reactor(
            pool, top_k_catalysts, primary='E_act', min_per_class=1
        )
        for rank_idx, (_, row) in enumerate(selected.iterrows(), 1):
            c_disp = row.get('interfacial_carbon_displacement_A')
            i_score = row.get('interfacial_pyrolysis_score_eV')
            coking = row.get('coking_index')
            discovery_champions.append(
                {
                    'rank': rank_idx,
                    'candidate_id': str(row.get('candidate_id', '')),
                    'material_class': str(row.get('material_class', '')),
                    'genome': str(row.get('genome', '')),
                    'E_act_eV': round(float(row['E_act']), 4),
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
                    'coking_index_eV': (
                        None
                        if coking is None or pd.isna(coking)
                        else round(float(coking), 4)
                    ),
                }
            )

    # 2. Phase 2 6-Pathway Reactor & TEA Leaderboard
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

    pathway_leaderboard: dict[str, dict] = {}
    for mode_name, predicate in pathway_filters.items():
        matching = [r for r in reactor_records if predicate(r)]
        if not matching:
            pathway_leaderboard[mode_name] = {'status': 'no_usable_records'}
            continue
        best_rec = max(
            matching,
            key=lambda r: (single_pass_x(r), -float(r.get('T_K') or 0.0)),
        )
        conv = float(single_pass_x(best_rec))
        tea_range = estimate_scenario_range(conv) if 0.0 < conv <= 1.0 else None
        pathway_leaderboard[mode_name] = {
            'status': 'complete',
            'n_usable_records': len(matching),
            'champion_catalyst': str(best_rec.get('catalyst_name', '')),
            'reactor_type': str(best_rec.get('reactor_type', '')),
            'T_K': float(best_rec.get('T_K') or 0.0),
            'single_pass_CH4_conversion': round(conv, 4),
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

    # 3. Phase 5 PEMFC Cathode & Stack TEA
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
        stack_res = model_stack(
            StackConfig(
                n_cells=stack_cells,
                cell_voltage_V=float(best_cell.get('peak_voltage_V', 0.65)),
                current_density_A_cm2=float(best_cell.get('peak_current_A_cm2', 1.5)),
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

    return {
        'discovery_champions': discovery_champions,
        'pathway_leaderboard': pathway_leaderboard,
        'pemfc_stack_summary': pemfc_stack_summary,
        'coupled_turquoise_pemfc_tea': coupled_tea,
        'evidence_level': 'screening_multipathway_scorecard',
    }
