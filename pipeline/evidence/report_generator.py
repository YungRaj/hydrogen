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
    records = []
    for f in list(directory.glob("*.json")):
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
    if screening_files:
        import pandas as pd

        data['screening'] = {}
        for f in screening_files:
            try:
                data['screening'][f.stem] = pd.read_csv(f)
            except Exception:
                pass

    data['reactor'] = _load_json_dir(REACTOR_DIR)
    data['dft'] = _load_json_dir(DFT_DIR)
    data['vqe'] = _load_json_dir(VQE_DIR)
    data['fuel_cell'] = _load_json_dir(FUEL_CELL_DIR)
    data['pipeline_state'] = load_json("pipeline_state.json") or {}
    return data


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
        judge = sc.get('judge_catalyst') or sc.get('headline_catalyst') or 'solids'
        lines.append(f"| Solids judge CH₄ conversion ({judge} PFR) | {solids_x:.2%} |")
    else:
        lines.append("| Solids judge CH₄ conversion (PFR) | N/A |")
    cand_solids_x = _scorecard_value(
        p2, 'candidate_champion_conversion', 'candidate_champion_conversion'
    )
    if isinstance(cand_solids_x, (int, float)):
        sc = p2.get('solids_scorecard') or {}
        cand_cat = (
            p2.get('candidate_champion_catalyst')
            or sc.get('candidate_champion_catalyst')
            or 'candidate'
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
    if 'E_act' in ranked.columns and len(ranked) > 0:
        lines += [f"\n**Top 10 Catalysts by Activation Barrier:**\n",
                  "| Rank | Genome | E_act (eV) | Coking Index | dE_H (eV) |",
                  "|------|--------|------------|-------------|-----------|"]
        for rank, (_, row) in enumerate(ranked.nsmallest(10, 'E_act').iterrows(), 1):
            genome_str = str(row.get('genome', ''))[:50]
            lines.append(
                f"| {rank} | {genome_str} | {row['E_act']:.4f} | "
                f"{row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |")
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
        judge = (scorecard.get('judge_catalyst') or scorecard.get('headline_catalyst')
                 or 'best non-H-parked solids')
        lines += [
            f"Solids judged on **{judge}** (single-pass X, `a`, WHSV, ΔP). "
            "MMBCR X_eq and 0.01 eV H-parked cats are not ranks.\n",
            "| Reactor | Catalyst " + SOLIDS_TABLE_TAIL,
            "|---------|----------|-------|---------------|---------|------------|----------|",
        ]
        for rt, row in scorecard['headline'].items():
            lines.append(f"| {rt} | {row.get('catalyst_name', '?')} | {row.get('T_K', '?')} | "
                         + _solids_columns(row, row.get('single_pass_CH4_conversion')))
        lines.append("")
        if scorecard.get('candidate_champion_conversion') is not None:
            lines.append(
                f"Top candidate solids X = {_report_number(scorecard['candidate_champion_conversion'], '.2%')} "
                f"on **{scorecard.get('candidate_champion_catalyst', '?')}** "
                f"({scorecard.get('candidate_champion_reactor', '?')} @ "
                f"{scorecard.get('candidate_champion_T_K', '?')} K)\n"
            )
        if scorecard.get('mmbcr_max_conversion') is not None:
            lines.append(
                f"MMBCR max X = {_report_number(scorecard['mmbcr_max_conversion'], '.1%')} "
                f"({scorecard.get('mmbcr_note', 'not a solids rank')})\n")

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
            lines.append(f"| {res.get('catalyst_name', '?')} | {res.get('reactor_type', '?')} | "
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
            lines.append(
                f"| {res.get('catalyst_name', '?')} | {res.get('reactor_type', '?')} | "
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
                lines.append(
                    f"| {res.get('catalyst_name', '?')} | {res.get('material_class', '?')} | "
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
            lines.append(
                f"| {res.get('catalyst_name', '?')} | {res.get('reaction_type', '?')} | "
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
            lines.append(
                f"| {res.get('cathode_catalyst', '?')} | {res.get('membrane', '?')} | "
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
            lines.append(
                f"| {res.get('catalyst_name', '?')} | {res.get('orr_overpotential_V', '?'):.3f} | "
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
    total_s = ps.get('total_elapsed_s', 0)
    lines += [f"| **Total** | **{total_s:.0f}s ({total_s/3600:.1f} hours)** |", ""]
    return lines


def generate_full_report(pipeline_state: Dict = None) -> Path:
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
    p1 = ps.get('phase1', {})
    p2 = ps.get('phase2', {})
    # FC results come from the orchestrator ('phase5') or the campaign
    # ('phase5_branch' search summary + 'phase5_stack' PEMFC/stack results).
    p5 = {**ps.get('phase5', {}), **ps.get('phase5_branch', {}), **ps.get('phase5_stack', {})}

    report_lines = (
        _header()
        + _executive_summary(p1, p2, p5)
        + _screening_section(data.get('screening'))
        + _reactor_section(p2, data['reactor'])
        + _dft_section(data['dft'])
        + _vqe_section(data['vqe'])
        + _fuel_cell_section(data['fuel_cell'], data['dft'])
        + TEA_SECTION
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
