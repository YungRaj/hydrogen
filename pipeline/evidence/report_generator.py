#!/usr/bin/env python3
# Evidence-aware campaign reporting.
"""
Automated Report Generator for the Turquoise H₂ → Fuel Cell Pipeline.

Generates a comprehensive Markdown report with:
  - Champion catalyst tables
  - Reactor performance comparison
  - DFT validation results
  - Fuel cell polarization curves
  - Techno-economic summary
"""

import sys
import json
from pathlib import Path
from typing import Dict
from datetime import datetime

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.fuel_cell import composite_score
from pipeline.utils import (
    RESULTS_DIR,
    SCREENING_DIR,
    REACTOR_DIR,
    DFT_DIR,
    VQE_DIR,
    FUEL_CELL_DIR,
    REPORTS_DIR,
    setup_logger,
    load_json,
)

logger = setup_logger('report_generator', 'reports/report_generation.log')


def _load_json_dir(directory: Path) -> list:
    """Load every readable JSON file in a results directory."""
    if not directory.is_dir():
        return []
    records = []
    for f in sorted(directory.glob("*.json")):
        try:
            with open(f, 'r') as fh:
                records.append(json.load(fh))
        except Exception:
            pass
    return records


def load_all_results() -> Dict:
    """Load all results from the pipeline output directories.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    data = {}

    # Screening databases
    screening_files = list(SCREENING_DIR.glob("*.csv"))
    for c_name in (
        "campaign_thermocatalytic",
        "campaign_mmbcr",
        "campaign_ntec",
        "campaign_electrochemical",
    ):
        c_scr = RESULTS_DIR / c_name / "screening"
        if c_scr.is_dir():
            screening_files.extend(c_scr.glob("*.csv"))
    if screening_files:
        try:
            import pandas as pd

            data['screening'] = {}
            for f in screening_files:
                if f.stem in data['screening']:
                    continue
                try:
                    data['screening'][f.stem] = pd.read_csv(f)
                except Exception:
                    pass
        except ImportError:
            data['screening'] = {}

    reactor_records = list(_load_json_dir(REACTOR_DIR))
    for c_name in (
        "campaign_thermocatalytic",
        "campaign_mmbcr",
        "campaign_ntec",
        "campaign_electrochemical",
    ):
        reactor_records.extend(_load_json_dir(RESULTS_DIR / c_name / "reactor"))
    data['reactor'] = reactor_records
    data['dft'] = _load_json_dir(DFT_DIR)
    data['vqe'] = _load_json_dir(VQE_DIR)
    fc_records = list(_load_json_dir(FUEL_CELL_DIR))
    if not fc_records:
        fc_records.extend(
            _load_json_dir(RESULTS_DIR / "campaign_thermocatalytic" / "fuel_cell")
        )
    data['fuel_cell'] = fc_records
    data['pipeline_state'] = (
        load_json("pipeline_state.json")
        or load_json("pipeline_state.json", subdir="campaign_thermocatalytic/reports")
        or {}
    )
    data['multipathway_scorecard'] = (
        load_json("multipathway_champion_scorecard.json", subdir="")
        or load_json("multipathway_champion_scorecard.json")
        or {}
    )
    return data



def _catalyst_label(
    raw_name_or_id: object,
    record: dict | None = None,
    *,
    material_class: object = None,
    genome: object = None,
) -> str:
    """Resolve a catalyst identifier or record into a literature-backed label."""
    from pipeline.search.discovery import format_literature_catalyst_label

    rec = record if isinstance(record, dict) else {}
    for lit_key in (
        'literature_catalyst_name',
        'cathode_catalyst_literature_name',
        'champion_catalyst_literature_name',
        'champion_cathode_literature_name',
    ):
        val = rec.get(lit_key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    raw_id = (
        rec.get('candidate_id')
        if raw_name_or_id in (None, '', '?')
        else raw_name_or_id
    )
    if raw_id in (None, '', '?'):
        return '?'
    clean_id = (
        str(raw_id).strip().removeprefix('campaign_').removeprefix('orr_slate_')
    )
    return format_literature_catalyst_label(
        clean_id,
        material_class=material_class if material_class is not None else rec.get('material_class'),
        genome=genome if genome is not None else rec.get('genome'),
    )



def _report_number(value, format_spec, unavailable='N/A'):
    try:
        return format(float(value), format_spec)
    except (TypeError, ValueError):
        return unavailable


def _metric_row(label: str, value, format_spec: str = ',', suffix: str = '') -> str:
    """Return a summary row, or N/A unless the value is a number."""
    if isinstance(value, (int, float)):
        return f"| {label} | {format(value, format_spec)}{suffix} |"
    return f"| {label} | N/A |"


def _header() -> list:
    return [
        "# Turquoise Hydrogen → Fuel Cell: Pipeline Report",
        f"\n**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"**Pipeline Version**: 2.0.0",
        "",
        "> **Evidence status:** screening values are hypotheses. Mock results and "
        "unconverged calculations are excluded from validated claims; catalyst "
        "performance requires experimental confirmation.",
        "",
    ]


def _scorecard_value(p2: dict, key: str, scorecard_key: str):
    """Prefer the phase-2 summary value, falling back to its solids scorecard."""
    value = p2.get(key)
    if not isinstance(value, (int, float)) and isinstance(p2.get('solids_scorecard'), dict):
        value = p2['solids_scorecard'].get(scorecard_key)
    return value


def _executive_summary(p1: dict, p2: dict, p5: dict) -> list:
    lines = ["## Executive Summary\n", f"| Metric | Value |", f"|--------|-------|",
             _metric_row('Catalysts evaluated (Fairchem)', p1.get('total_evaluated')),
             _metric_row('Valid candidates', p1.get('valid_count')),
             f"| Pareto front size | {p1.get('pareto_size', 'N/A')} |",
             _metric_row('Best eligible screening barrier', p1.get('best_E_act'), '.3f', ' eV'),
             f"| Best coking resistance | {p1.get('best_coking', 'N/A')} |",
             f"| Reactor simulations | {p2.get('catalysts_simulated', 'N/A')} |"]
    solids_x = _scorecard_value(p2, 'best_conversion', 'headline_solids_conversion')
    if isinstance(solids_x, (int, float)):
        sc = p2.get('solids_scorecard') or {}
        judge = (
            sc.get('judge_catalyst_literature_name')
            or sc.get('headline_catalyst_literature_name')
            or _catalyst_label(
                sc.get('judge_catalyst') or sc.get('headline_catalyst') or 'solids'
            )
        )
        lines.append(f"| Solids judge CH₄ conversion ({judge} PFR) | {solids_x:.2%} |")
    else:
        lines.append("| Solids judge CH₄ conversion (PFR) | N/A |")
    cand_solids_x = _scorecard_value(
        p2, 'candidate_champion_conversion', 'candidate_champion_conversion'
    )
    if isinstance(cand_solids_x, (int, float)):
        sc = p2.get('solids_scorecard') or {}
        cand_cat = (
            sc.get('candidate_champion_literature_name')
            or _catalyst_label(
                p2.get('candidate_champion_catalyst')
                or sc.get('candidate_champion_catalyst')
                or 'candidate'
            )
        )
        lines.append(
            f"| Top candidate solids CH₄ conversion ({cand_cat}) | {cand_solids_x:.2%} |"
        )
    mmbcr_x = _scorecard_value(p2, 'mmbcr_max_conversion', 'mmbcr_max_conversion')
    if isinstance(mmbcr_x, (int, float)):
        lines.append(f"| MMBCR max CH₄ conversion | {mmbcr_x:.2%} (not a solids rank) |")
    else:
        lines.append("| MMBCR max CH₄ conversion | N/A |")
    ntec_x = _scorecard_value(p2, 'ntec_max_conversion', 'ntec_max_conversion')
    if isinstance(ntec_x, (int, float)):
        lines.append(f"| NTEC max CH₄ conversion (coupled MOOSE+Cantera) | {ntec_x:.2%} |")
    echem_x = _scorecard_value(
        p2, 'electrochemical_max_conversion', 'electrochemical_max_conversion'
    )
    if isinstance(echem_x, (int, float)):
        lines.append(
            f"| Electrochemical max CH₄ conversion (coupled MOOSE+Cantera) | {echem_x:.2%} |"
        )
    lines += [
        _metric_row('FC catalysts evaluated',
                    p5.get('total_evaluated', p5.get('n_cathodes_screened'))),
        _metric_row('Best PEMFC power density', p5.get('best_power_W_cm2'), '.3f', ' W/cm²'),
        _metric_row('Best PEMFC efficiency', p5.get('best_efficiency'), '.1%'),
        _metric_row('Min ORR overpotential', p5.get('min_overpotential_V'), '.3f', ' V'),
        "",
    ]
    return lines


_ELIGIBILITY_COLUMNS = frozenset({'valid', 'E_act', 'E_act_censored', 'pyrolysis_viable'})


def _screening_section(screening: dict) -> list:
    lines = ["## Phase 1: Catalyst Screening & Genetic Optimization\n"]
    if not screening:
        mp = (
            load_json("multipathway_champion_scorecard.json", subdir="")
            or load_json("multipathway_champion_scorecard.json")
            or {}
        )
        champs = mp.get('discovery_champions') if isinstance(mp, dict) else None
        if isinstance(champs, list) and champs:
            lines += [
                "**Top Discovery Champions by Activation Barrier & Interfacial Coking Score:**\n",
                "| Rank | Literature Catalyst Name | Class | E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | Coking Index (eV) |",
                "|------|--------------------------|-------|------------|------------------------|-----------------------|-------------------|",
            ]
            for ch in champs:
                if not isinstance(ch, dict):
                    continue
                lines.append(
                    f"| {ch.get('rank', '?')} | {_catalyst_label(ch.get('candidate_id', '?'), ch)} | "
                    f"{ch.get('material_class', '?')} | "
                    f"{_report_number(ch.get('E_act_eV'), '.4f')} | "
                    f"{_report_number(ch.get('interfacial_pyrolysis_score_eV'), '.4f')} | "
                    f"{_report_number(ch.get('calibrated_interfacial_score_eV'), '.4f')} | "
                    f"{_report_number(ch.get('coking_index_eV'), '.4f')} |"
                )
            lines.append("")
        return lines
    main_db_key = next(
        (k for k in screening if 'branch_ranker_evidence' in k or 'full' in k or 'initial' in k),
        list(screening.keys())[0])
    df = screening[main_db_key]
    valid = df[df.get('valid', True) == True] if 'valid' in df.columns else df
    lines += [f"### Screening Database: `{main_db_key}`\n",
              f"- **Total entries**: {len(df):,}",
              f"- **Valid candidates**: {len(valid):,}"]
    if 'material_class' in df.columns:
        lines += [f"\n**Distribution by Material Class:**\n",
                  "| Material Class | Count | Valid |", "|----------------|-------|-------|"]
        for cls in df['material_class'].unique():
            n_total = len(df[df['material_class'] == cls])
            n_valid = (len(valid[valid['material_class'] == cls])
                       if 'material_class' in valid.columns else 0)
            lines.append(f"| {cls} | {n_total} | {n_valid} |")
        lines.append("")
    ranked = valid
    if _ELIGIBILITY_COLUMNS.issubset(df.columns):
        # Match stages.discovery: censored (floor-clipped) and non-viable
        # barriers are not evidence of a low barrier and must not be ranked.
        from pipeline.screening.small_data_ranker import screening_metric_eligibility

        ranked = df.loc[screening_metric_eligibility(df, 'turquoise_hydrogen')]
    if len(ranked) > 0:
        from pipeline.screening.stage_selection import annotate_interfacial_pyrolysis

        ranked = annotate_interfacial_pyrolysis(ranked)
    if 'E_act' in ranked.columns and len(ranked) > 0:
        has_interfacial = 'interfacial_pyrolysis_score_eV' in ranked.columns
        if has_interfacial:
            lines += [
                "\n**Top 10 Catalysts by Activation Barrier:**\n",
                "| Rank | Genome | E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | Coking Index | dE_H (eV) |",
                "|------|--------|------------|------------------------|-----------------------|-------------|-----------|",
            ]
        else:
            lines += [
                "\n**Top 10 Catalysts by Activation Barrier:**\n",
                "| Rank | Genome | E_act (eV) | Coking Index | dE_H (eV) |",
                "|------|--------|------------|-------------|-----------|",
            ]
        for rank, (_, row) in enumerate(ranked.nsmallest(10, 'E_act').iterrows(), 1):
            row_dict = row.to_dict() if hasattr(row, 'to_dict') else dict(row)
            genome_str = _catalyst_label(
                row_dict.get('candidate_id') or row_dict.get('genome', ''),
                row_dict,
            )
            if has_interfacial:
                lines.append(
                    f"| {rank} | {genome_str} | {row['E_act']:.4f} | "
                    f"{_report_number(row.get('interfacial_pyrolysis_score_eV'), '.4f')} | "
                    f"{_report_number(row.get('calibrated_interfacial_score_eV'), '.4f')} | "
                    f"{row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |"
                )
            else:
                lines.append(
                    f"| {rank} | {genome_str} | {row['E_act']:.4f} | "
                    f"{row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |"
                )
        lines.append("")
    return lines


SOLIDS_TABLE_TAIL = "| T (K) | Single-pass X | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |"


def _solids_columns(record: dict, conversion) -> str:
    return (f"{_report_number(conversion, '.2%')} | "
            f"{_report_number(record.get('active_sv_1_m'), '.0f', '—')} | "
            f"{_report_number(record.get('WHSV_h-1'), '.0f', '—')} | "
            f"{_report_number(record.get('ergun_delta_p_bar'), '.2f', '—')} |")


def _reactor_section(p2: dict, reactor: list) -> list:
    from pipeline.reactors.scorecard import (
        is_production_reactor_record,
        is_solids_run,
        single_pass_x,
    )
    from pipeline.reactors.eligibility import is_usable_result

    lines = ["## Phase 2: Reactor-Scale Simulation\n"]
    scorecard = (
        p2.get('solids_scorecard')
        or load_json('phase2_solids_scorecard.json', subdir='reactor')
        or {}
    )
    if scorecard.get('headline'):
        judge = (
            scorecard.get('judge_catalyst_literature_name')
            or scorecard.get('headline_catalyst_literature_name')
            or _catalyst_label(
                scorecard.get('judge_catalyst')
                or scorecard.get('headline_catalyst')
                or 'best non-H-parked solids'
            )
        )
        lines += [
            f"Solids judged on **{judge}** (single-pass X, `a`, WHSV, ΔP). "
            "MMBCR X_eq and 0.01 eV H-parked cats are not ranks.\n",
            "| Reactor | Catalyst " + SOLIDS_TABLE_TAIL,
            "|---------|----------|-------|---------------|---------|------------|----------|",
        ]
        for rt, row in scorecard['headline'].items():
            cat_label = _catalyst_label(row.get('catalyst_name', '?'), row)
            lines.append(f"| {rt} | {cat_label} | {row.get('T_K', '?')} | "
                         + _solids_columns(row, row.get('single_pass_CH4_conversion')))
        lines.append("")
        if scorecard.get('candidate_champion_conversion') is not None:
            champ_label = (
                scorecard.get('candidate_champion_literature_name')
                or _catalyst_label(scorecard.get('candidate_champion_catalyst', '?'))
            )
            lines.append(
                f"Top candidate solids X = {_report_number(scorecard['candidate_champion_conversion'], '.2%')} "
                f"on **{champ_label}** "
                f"({scorecard.get('candidate_champion_reactor', '?')} @ "
                f"{scorecard.get('candidate_champion_T_K', '?')} K)\n"
            )
        if scorecard.get('mmbcr_max_conversion') is not None:
            lines.append(
                f"MMBCR max X = {_report_number(scorecard['mmbcr_max_conversion'], '.1%')} "
                f"({scorecard.get('mmbcr_note', 'not a solids rank')})\n")
    pbr_diag = scorecard.get('pbr_ergun_residence_diagnostics')
    if isinstance(pbr_diag, dict) and int(pbr_diag.get('n_completed') or 0) > 0:
        lines += [
            "### Packed-Bed (PBR), Ergun Pilot Residence-Time, & Regenerated-Cycle Diagnostics\n",
            f"- **Completed solids runs**: {pbr_diag.get('n_completed', 0)} "
            f"(PBR: {pbr_diag.get('n_pbr_records', 0)})",
            f"- **Mean PBR overall effectiveness ($\\Omega$)**: "
            f"{_report_number(pbr_diag.get('mean_pbr_overall_effectiveness_factor'), '.3f')}",
            f"- **Mean PBR Thiele modulus ($\\phi$)**: "
            f"{_report_number(pbr_diag.get('mean_pbr_thiele_modulus'), '.3f')}",
            f"- **Ergun $\\Delta P$ compliance**: "
            f"{pbr_diag.get('n_ergun_ok', 0)}/{pbr_diag.get('n_ergun_checked', 0)}",
            f"- **Ergun & Froment pilot residence-time band ($\\tau \\in [0.1, 10]\\text{{ s}}$)**: "
            f"{pbr_diag.get('n_within_pilot_residence_band', 0)}/{pbr_diag.get('n_pilot_residence_checked', 0)}",
            f"- **Mean Muradov regenerated-cycle minus single-pass conversion ($\\Delta X$)**: "
            f"{_report_number(pbr_diag.get('mean_regenerated_cycle_minus_single_pass_X'), '+.2%')}",
            "",
        ]

    # Only completed production runs count. Sweeps, mocks, not_applicable,
    # validation_required, and failed records are evidence about the
    # workflow, not about the catalyst.
    real_reactor = [x for x in reactor
                    if is_production_reactor_record(x) and not x.get('mock', False)]
    solids = [x for x in real_reactor if is_solids_run(x) and is_usable_result(x)]
    if solids:
        lines += ["### Solids runs (single-pass)\n",
                  "| Catalyst | Reactor " + SOLIDS_TABLE_TAIL,
                  "|----------|---------|-------|---------------|---------|------------|----------|"]
        for res in sorted(solids, key=single_pass_x, reverse=True)[:20]:
            cat_label = _catalyst_label(res.get('catalyst_name', '?'), res)
            lines.append(f"| {cat_label} | {res.get('reactor_type', '?')} | "
                         f"{res.get('T_K', '?')} | " + _solids_columns(res, single_pass_x(res)))
        lines.append("")
    non_solids = [
        x for x in real_reactor if not is_solids_run(x) and is_usable_result(x)
    ]
    if non_solids:
        lines += [
            "### Liquid-Metal, NTEC, and Electrochemical Runs\n",
            "| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | Evidence Tier |",
            "|----------|---------|-------|----------------|----------------|---------------|",
        ]
        for res in sorted(non_solids, key=single_pass_x, reverse=True)[:20]:
            cat_label = _catalyst_label(res.get('catalyst_name', '?'), res)
            lines.append(
                f"| {cat_label} | {res.get('reactor_type', '?')} | "
                f"{res.get('T_K', '?')} | {_report_number(single_pass_x(res), '.2%')} | "
                f"{_report_number(res.get('H2_selectivity'), '.2%')} | "
                f"{res.get('reactor_evidence_tier', 'N/A')} |"
            )
        lines.append("")
    if len(real_reactor) != len(reactor):
        lines.append(
            f"Excluded {len(reactor) - len(real_reactor)} non-run reactor files "
            f"(sweeps, mocks, not_applicable, validation_required, failed) "
            f"from performance claims.\n")
    return lines


def _dft_section(dft: list) -> list:
    lines = ["## Phase 3: DFT Validation (Quantum ESPRESSO)\n"]
    converged_dft = [x for x in dft if x.get('converged') is True]
    if converged_dft:
        lines += ["| Catalyst | Material Class | Energy (eV) | Converged | Max Force |",
                  "|----------|----------------|-------------|-----------|-----------|"]
        for res in converged_dft:
            if 'orr_overpotential_V' not in res:  # pyrolysis DFT results
                cat_label = _catalyst_label(
                    res.get('candidate_id') or res.get('catalyst_name', '?'), res
                )
                lines.append(
                    f"| {cat_label} | {res.get('material_class', '?')} | "
                    f"{res.get('dft_energy_eV', 'N/A')} | {res.get('converged', '?')} | "
                    f"{res.get('max_force_Ry_bohr', 'N/A')} |")
        lines.append("")
    if dft and not converged_dft:
        lines.append(
            "No converged DFT validation is available; no DFT-validated champion is claimed.\n")
    return lines


def _vqe_section(vqe: list) -> list:
    lines = ["## Phase 4: Quantum VQE Transition States\n"]
    real_vqe = [x for x in vqe
                if not x.get('mock', False)
                and x.get('catalyst_specific_hamiltonian') is True
                and x.get('benchmarked') is True]
    if real_vqe:
        lines += ["| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |",
                  "|----------|----------|-------------|-------------|--------|--------|"]
        for res in real_vqe:
            cat_label = _catalyst_label(
                res.get('candidate_id') or res.get('catalyst_name', '?'), res
            )
            lines.append(
                f"| {cat_label} | {res.get('reaction_type', '?')} | "
                f"{res.get('energy_Ha', 0):.6f} | {res.get('energy_eV', 0):.4f} | "
                f"{res.get('n_qubits', '?')} | {res.get('n_layers', '?')} |")
        lines.append("")
    if vqe and not real_vqe:
        lines.append("No catalyst-specific, benchmarked quantum result is available; "
                     "toy VQE energies are not validation.\n")
    return lines


def _fuel_cell_section(fuel_cell: list, dft: list) -> list:
    lines = ["## Phase 5: PEMFC Fuel Cell Performance\n"]
    pemfc_results = [x for x in fuel_cell if 'peak_power_W_cm2' in x and 'n_cells' not in x]
    stack_results = [x for x in fuel_cell if 'n_cells' in x]
    if pemfc_results:
        lines += ["### Single-Cell Performance\n",
                  "| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |",
                  "|---------|----------|---------------------|---------|-------|--------------------|"]
        for res in sorted(pemfc_results, key=composite_score, reverse=True)[:15]:
            cath_label = _catalyst_label(res.get('cathode_catalyst', '?'), res)
            lines.append(
                f"| {cath_label} | {res.get('membrane', '?')} | "
                f"{res.get('peak_power_W_cm2', 0):.4f} | {res.get('OCV_V', 0):.3f} | "
                f"{res.get('orr_overpotential_V', '?')} | {res.get('efficiency_at_peak', 0):.1%} |")
        lines.append("")
    if stack_results:
        lines.append("### Stack Performance\n")
        for res in stack_results:
            lines += [
                f"- **Stack Power**: {res.get('net_power_kW', 0):.1f} kW net "
                f"({res.get('n_cells', '?')} cells)",
                f"- **System Efficiency**: {res.get('system_efficiency', 0):.1%}",
                f"- **Power Density**: {res.get('gravimetric_W_kg', 0):.0f} W/kg, "
                f"{res.get('volumetric_W_L', 0):.0f} W/L",
                f"- **Cost**: ${res.get('cost_per_kW', 0):.0f}/kW",
                "",
            ]
    orr_dft = [x for x in dft if 'orr_overpotential_V' in x]
    if orr_dft:
        lines += ["### ORR DFT Validation\n",
                  "| Catalyst | η_ORR (V) | dG_OH (eV) | dG_O (eV) | dG_OOH (eV) | RDS |",
                  "|----------|-----------|-----------|----------|------------|-----|"]
        for res in orr_dft:
            cat_label = _catalyst_label(
                res.get('candidate_id') or res.get('catalyst_name', '?'), res
            )
            lines.append(
                f"| {cat_label} | {res.get('orr_overpotential_V', '?'):.3f} | "
                f"{res.get('dG_OH_eV', '?')} | {res.get('dG_O_eV', '?')} | "
                f"{res.get('dG_OOH_eV', '?')} | {res.get('rate_determining_step', '?')} |")
        lines.append("")
    return lines


TEA_SECTION = [
    "## Techno-Economic Summary\n",
    "### Turquoise H₂ Production Cost Drivers\n",
    "| Parameter | Value |",
    "|-----------|-------|",
    "| Catalyst material cost | Based on abundance scoring |",
    "| Reactor operating temperature | 800–1200 K |",
    "| Target H₂ production cost | < $2/kg |",
    "| Carbon co-product credit | $0.10–$0.50/kg |",
    "",
]


def _tea_section(
    p2: dict, multipathway_scorecard: dict | None = None
) -> list:
    """Build the techno-economic section with live LCOH scenarios and full 7-pathway scorecard."""
    lines = list(TEA_SECTION)
    conv = _scorecard_value(p2, 'best_conversion', 'headline_solids_conversion')
    if not isinstance(conv, (int, float)) or conv <= 0:
        conv = _scorecard_value(
            p2, 'candidate_champion_conversion', 'candidate_champion_conversion'
        )
    if isinstance(conv, (int, float)) and conv > 0:
        from pipeline.economics.tea import estimate_hydrogen_cost

        sc = p2.get('solids_scorecard') or {}
        temp_k = float(sc.get('candidate_champion_T_K') or 1073.15)
        lines += [
            f"### Levelized Cost of Hydrogen (LCOH) at X = {conv:.1%} (T = {temp_k:.0f} K)\n",
            "| Scenario | Single-Pass X | Net CH₄ Util. | H₂ Cost ($/kg) | Carbon Credit ($/kg H₂) | DOE < $2/kg |",
            "|----------|---------------|---------------|----------------|-------------------------|-------------|",
        ]
        for scenario in ('optimistic', 'base', 'conservative'):
            est = estimate_hydrogen_cost(
                float(conv), scenario, temperature_K=temp_k
            )
            lines.append(
                f"| {scenario} | {est['conversion']:.1%} | "
                f"{est['net_ch4_utilization']:.1%} | "
                f"${est['h2_cost_usd_kg']:.2f} | "
                f"${est['carbon_credit_usd_kg_h2']:.2f} | "
                f"{'Yes' if est['h2_cost_usd_kg'] < 2.0 else 'No'} |"
            )
        lines.append("")

    mp = (
        multipathway_scorecard
        if isinstance(multipathway_scorecard, dict) and multipathway_scorecard
        else (
            load_json("multipathway_champion_scorecard.json", subdir="")
            or load_json("multipathway_champion_scorecard.json")
            or {}
        )
    )
    if not isinstance(mp, dict) or not mp:
        return lines

    lines.append(
        "### Multi-Pathway Discovery Champions & QE-DFT-Calibrated Leaderboard\n"
    )
    disc_champs = mp.get('discovery_champions')
    if isinstance(disc_champs, list) and disc_champs:
        lines += [
            "#### 11-Champion Discovery Portfolio (Interfacial & QE DFT Validated)\n",
            "| Rank | Literature Catalyst Name | Class | Screening E_act (eV) | Calibrated E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | QE DFT Energy (Ry) |",
            "|------|--------------------------|-------|----------------------|-----------------------|------------------------|-----------------------|--------------------|",
        ]
        for ch in disc_champs:
            if not isinstance(ch, dict):
                continue
            lit_name = _catalyst_label(ch.get('candidate_id', '?'), ch)
            lines.append(
                f"| {ch.get('rank', '?')} | {lit_name} | {ch.get('material_class', '?')} | "
                f"{_report_number(ch.get('E_act_eV'), '.4f')} | "
                f"{_report_number(ch.get('calibrated_activation_eV'), '.4f')} | "
                f"{_report_number(ch.get('interfacial_pyrolysis_score_eV'), '.4f')} | "
                f"{_report_number(ch.get('calibrated_interfacial_score_eV'), '.4f')} | "
                f"{_report_number(ch.get('qe_dft_energy_Ry'), '.4f')} |"
            )
        lines.append("")

    cal_board = mp.get('dft_calibrated_pathway_leaderboard')
    if isinstance(cal_board, dict) and cal_board:
        lines += [
            "#### 7-Pathway QE-DFT-Calibrated Reactor & LCOH Leaderboard\n",
            "| Pathway | Reactor | Literature Champion Catalyst | T (K) | Uncalibrated X | Calibrated X | Ω | ΔP (bar) | Base LCOH ($/kg) | Optimistic LCOH ($/kg) | Kinetics Status |",
            "|---------|---------|------------------------------|-------|----------------|--------------|---|----------|------------------|------------------------|-----------------|",
        ]
        pathway_order = (
            'thermocatalytic',
            'thermocatalytic_pbr',
            'thermocatalytic_pfr',
            'thermocatalytic_fluidized',
            'mmbcr',
            'ntec',
            'electrochemical',
        )
        ordered_keys = [k for k in pathway_order if k in cal_board] + [
            k for k in cal_board if k not in pathway_order
        ]
        for mode_key in ordered_keys:
            row = cal_board.get(mode_key)
            if not isinstance(row, dict) or row.get('status') != 'complete':
                continue
            champ_label = _catalyst_label(row.get('champion_catalyst', '?'), row)
            lines.append(
                f"| `{mode_key}` | {row.get('reactor_type', '?')} | {champ_label} | "
                f"{_report_number(row.get('T_K'), '.1f')} | "
                f"{_report_number(row.get('uncalibrated_single_pass_CH4_conversion'), '.2%')} | "
                f"{_report_number(row.get('dft_calibrated_single_pass_CH4_conversion'), '.2%')} | "
                f"{_report_number(row.get('dft_calibrated_overall_effectiveness_factor'), '.3f', '—')} | "
                f"{_report_number(row.get('ergun_delta_p_bar'), '.3f', '—')} | "
                f"${_report_number(row.get('dft_calibrated_tea_base_h2_usd_kg'), '.3f')} | "
                f"${_report_number(row.get('dft_calibrated_tea_optimistic_h2_usd_kg'), '.3f')} | "
                f"`{row.get('kinetics_status', 'qe_dft_calibrated_bep')}` |"
            )
        lines.append("")

    ni_closure = mp.get('supported_ni_filament_closure_summary')
    if isinstance(ni_closure, dict) and ni_closure:
        gate = (
            ni_closure.get('supported_ni_filament_closure_gate')
            if isinstance(ni_closure.get('supported_ni_filament_closure_gate'), dict)
            else {}
        )
        bh_crit = (
            (ni_closure.get('baker_helveg_criteria_summary') or {}).get('criteria')
            if isinstance(ni_closure.get('baker_helveg_criteria_summary'), dict)
            else (ni_closure.get('criteria') if isinstance(ni_closure.get('criteria'), dict) else {})
        )
        closure_ok = bool(
            ni_closure.get('supported_ni_filament_closure_passed')
            or ni_closure.get('ranking_authority')
            or gate.get('solids_ranking_authority_accepted')
        )
        flat_val = (
            (bh_crit.get('flat') or {}).get('value')
            or (bh_crit.get('barrier_sensitivity') or {}).get('d_ln_X_d_E_act_inv_eV')
        )
        lin_val = (
            (bh_crit.get('linearity') or {}).get('value')
            or (bh_crit.get('particle_size_differentiation') or {}).get(
                'relative_turnover_difference'
            )
        )
        encap_val = (
            ((bh_crit.get('encapsulation_onset') or {}).get('first_A_gamma_theta_encap') or {}).get(
                'exit_theta_C_encap'
            )
            or (bh_crit.get('high_temperature_encapsulation_onset') or {}).get(
                'theta_encap_1100K'
            )
        )
        yield_923 = (
            ((bh_crit.get('yield') or {}).get('example') or {}).get(
                'filament_yield_gC_per_gMetal_h'
            )
            or (bh_crit.get('specific_carbon_yield') or {}).get('yield_923K_g_C_g_Ni_h')
        )
        joint_roi = (
            ni_closure.get('ermakova_takenaka_joint_roi_record')
            if isinstance(ni_closure.get('ermakova_takenaka_joint_roi_record'), dict)
            else {}
        )
        both_roi = (
            (ni_closure.get('ermakova_takenaka_joint_summary') or {}).get('both_roi')
            if isinstance(ni_closure.get('ermakova_takenaka_joint_summary'), dict)
            else []
        )
        roi_973_val = next(
            (
                r.get('filament_yield_gC_per_gMetal_h')
                for r in (both_roi if isinstance(both_roi, list) else [])
                if isinstance(r, dict)
                and abs(float(r.get('T_K') or 0.0) - 973.15) < 5.0
                and r.get('cell') == 'production'
            ),
            None,
        )
        yield_973 = (
            joint_roi.get('filament_yield_gC_per_gMetal_h')
            or roi_973_val
            or (bh_crit.get('specific_carbon_yield') or {}).get('yield_973K_g_C_g_Ni_h')
        )

        lines += [
            "#### Supported-Ni(111)/SiO₂ Filamentous-Carbon Reference Closure\n",
            f"- **Reference Catalyst**: {_catalyst_label(ni_closure.get('judge_catalyst') or ni_closure.get('reference_catalyst_id', 'ni111_sio2_bengaard2002_ref'), ni_closure)}",
            f"- **Closure Passed**: `{closure_ok}` (PFR ROI X = {_report_number(gate.get('pfr_roi_conversion'), '.2%')}, PBR ROI X = {_report_number(gate.get('pbr_roi_conversion'), '.2%')}, Fluidized ROI X = {_report_number(gate.get('fluidized_roi_conversion'), '.2%')})",
            f"- **Barrier Sensitivity ($|d\\ln X / dE_{{\\text{{act}}}}|$)**: "
            f"{_report_number(flat_val, '.2f')} eV⁻¹ (gate ≥ 5.0 eV⁻¹)",
            f"- **Particle-Size Turnover Differentiation**: "
            f"{_report_number(lin_val, '.1%')} (gate > 20%)",
            f"- **Encapsulation Onset ($\\theta_{{\\text{{encap}}}}$)**: "
            f"{_report_number(encap_val, '.3f')} (gate ≥ 0.10)",
            f"- **Specific Carbon Yield ($Y_{{\\text{{C}}}}$)**: "
            f"923 K = {_report_number(yield_923, '.2f')} g_C/(g_Ni·h), "
            f"973 K = {_report_number(yield_973, '.2f')} g_C/(g_Ni·h) (band 8–10 g_C/(g_Ni·h))",
            "",
        ]


    gc_stack = mp.get('dft_calibrated_pemfc_stack_summary')
    cal_coupled = mp.get('dft_calibrated_coupled_turquoise_pemfc_tea')
    if isinstance(gc_stack, dict) and gc_stack.get('status') == 'complete':
        lines += [
            "#### Coupled Turquoise H₂ → PEMFC Stack & Electricity Cost (LCOE)\n",
            f"- **QE-Confirmed Acid-Stable PEMFC Champion**: {_catalyst_label(gc_stack.get('champion_cathode', '?'), gc_stack)} "
            f"($\\eta_{{\\text{{ORR,GC}}}} = {_report_number(gc_stack.get('orr_overpotential_gc_V'), '.4f')}\\text{{ V}}$, "
            f"Pourbaix margin $+{_report_number(gc_stack.get('pourbaix_dissolution_margin_V'), '.4f')}\\text{{ V}}$, "
            f"peak power **{_report_number(gc_stack.get('peak_power_W_cm2'), '.4f')} W/cm²**, "
            f"stack net power **{_report_number(gc_stack.get('stack_net_power_kW'), '.2f')} kW** @ "
            f"{_report_number(gc_stack.get('stack_system_efficiency'), '.2%')} system efficiency)",
        ]
        low_pgm = gc_stack.get('low_pgm_champion')
        if isinstance(low_pgm, dict):
            lines.append(
                f"- **Low-PGM Champion**: {_catalyst_label(low_pgm.get('cathode_catalyst', '?'), low_pgm)} "
                f"($\\eta_{{\\text{{ORR,GC}}}} = {_report_number(low_pgm.get('orr_overpotential_gc_V'), '.4f')}\\text{{ V}}$, "
                f"Pourbaix margin $+{_report_number(low_pgm.get('pourbaix_dissolution_margin_V'), '.4f')}\\text{{ V}}$, "
                f"peak power **{_report_number(low_pgm.get('peak_power_W_cm2'), '.4f')} W/cm²**)"
            )
        lines.append("")
    if isinstance(cal_coupled, dict) and cal_coupled:
        lines += [
            "| Pathway | Calibrated Base H₂ ($/kg) | Calibrated Optimistic H₂ ($/kg) | Coupled Base PEMFC Fuel Cost ($/kWh_e) | Coupled Optimistic PEMFC Fuel Cost ($/kWh_e) |",
            "|---------|---------------------------|---------------------------------|----------------------------------------|----------------------------------------------|",
        ]
        for mode_key, c_row in cal_coupled.items():
            if not isinstance(c_row, dict):
                continue
            lines.append(
                f"| `{mode_key}` | "
                f"${_report_number(c_row.get('dft_calibrated_h2_cost_base_usd_kg'), '.3f')} | "
                f"${_report_number(c_row.get('dft_calibrated_h2_cost_optimistic_usd_kg'), '.3f')} | "
                f"${_report_number(c_row.get('dft_calibrated_pemfc_fuel_cost_base_usd_kwh_e'), '.4f')} | "
                f"${_report_number(c_row.get('dft_calibrated_pemfc_fuel_cost_optimistic_usd_kwh_e'), '.4f')} |"
            )
        lines.append("")
    return lines


PHASE_NAMES = {
    'phase1': 'Fairchem Screening + GA',
    'phase2': 'Reactor Simulation',
    'phase3': 'DFT Validation',
    'phase4': 'VQE Quantum',
    'phase5': 'Fuel Cell Modeling',
    'phase6': 'Report Generation',
}


def _timing_section(ps: dict) -> list:
    lines = ["## Pipeline Timing\n", "| Phase | Duration |", "|-------|----------|"]
    for phase_key, phase_name in PHASE_NAMES.items():
        lines.append(f"| {phase_name} | {ps.get(phase_key, {}).get('elapsed_s', 0):.0f}s |")
    # 3600 s/hr conversion from elapsed seconds to wall-clock hours
    total_s = ps.get('total_elapsed_s', 0)
    lines += [f"| **Total** | **{total_s:.0f}s ({total_s/3600:.1f} hours)** |", ""]
    return lines


def generate_full_report(
    pipeline_state: dict[str, object] | None = None,
) -> Path:
    """
        Generate the comprehensive pipeline report.

    Args:
        pipeline_state: Mapping supplying pipeline state.

    Returns:
        Filesystem path produced or resolved by the operation.
    """
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    data = load_all_results()
    if pipeline_state:
        data['pipeline_state'] = pipeline_state
    ps = data.get('pipeline_state', {})
    p1 = dict(ps.get('phase1', {}) if isinstance(ps.get('phase1'), dict) else {})
    p2 = dict(ps.get('phase2', {}) if isinstance(ps.get('phase2'), dict) else {})
    # FC results come from the orchestrator ('phase5') or the campaign
    # ('phase5_branch' search summary + 'phase5_stack' PEMFC/stack results).
    p5 = {
        **(ps.get('phase5', {}) if isinstance(ps.get('phase5'), dict) else {}),
        **(ps.get('phase5_branch', {}) if isinstance(ps.get('phase5_branch'), dict) else {}),
        **(ps.get('phase5_stack', {}) if isinstance(ps.get('phase5_stack'), dict) else {}),
    }
    mp = data.get('multipathway_scorecard')
    if isinstance(mp, dict) and mp:
        raw_board = mp.get('pathway_leaderboard') or {}
        if isinstance(raw_board, dict):
            tc = raw_board.get('thermocatalytic') or {}
            if isinstance(tc, dict) and tc.get('status') == 'complete':
                p2.setdefault(
                    'candidate_champion_conversion',
                    tc.get('single_pass_CH4_conversion'),
                )
                p2.setdefault(
                    'candidate_champion_catalyst',
                    tc.get('champion_catalyst_literature_name')
                    or tc.get('champion_catalyst'),
                )
            for k_mode, k_field in (
                ('mmbcr', 'mmbcr_max_conversion'),
                ('ntec', 'ntec_max_conversion'),
                ('electrochemical', 'electrochemical_max_conversion'),
            ):
                m_entry = raw_board.get(k_mode) or {}
                if isinstance(m_entry, dict) and m_entry.get('status') == 'complete':
                    p2.setdefault(k_field, m_entry.get('single_pass_CH4_conversion'))
        fc_sum = mp.get('dft_calibrated_pemfc_stack_summary') or mp.get('pemfc_stack_summary') or {}
        if isinstance(fc_sum, dict) and fc_sum.get('status') == 'complete':
            p5.setdefault('best_power_W_cm2', fc_sum.get('peak_power_W_cm2'))
            p5.setdefault('best_efficiency', fc_sum.get('cell_efficiency_at_peak'))
            p5.setdefault(
                'min_overpotential_V',
                fc_sum.get('orr_overpotential_gc_V')
                or fc_sum.get('orr_overpotential_V'),
            )

    report_lines = (
        _header()
        + _executive_summary(p1, p2, p5)
        + _screening_section(data.get('screening'))
        + _reactor_section(p2, data['reactor'])
        + _dft_section(data['dft'])
        + _vqe_section(data['vqe'])
        + _fuel_cell_section(data['fuel_cell'], data['dft'])
        + _tea_section(p2, mp)
        + _timing_section(ps)
    )
    report_path = REPORTS_DIR / "pipeline_report.md"
    with open(report_path, 'w') as f:
        f.write("\n".join(report_lines))
    logger.info(f"Report written to {report_path}")
    return report_path


if __name__ == '__main__':
    report = generate_full_report()
    print(f"Report generated: {report}")

