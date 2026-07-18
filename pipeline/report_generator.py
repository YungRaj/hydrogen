#!/usr/bin/env python3
"""
Automated Report Generator for the Turquoise H₂ → Fuel Cell Pipeline.

Generates a comprehensive Markdown report with:
  - Champion catalyst tables
  - Reactor performance comparison
  - DFT validation results
  - Fuel cell polarization curves
  - Techno-economic summary
"""

import os
import sys
import json
import glob
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional
from datetime import datetime

sys.path.insert(0, str(Path(__file__).parent.parent))
from pipeline.utils import (
    RESULTS_DIR, SCREENING_DIR, REACTOR_DIR, DFT_DIR,
    VQE_DIR, FUEL_CELL_DIR, REPORTS_DIR,
    setup_logger, save_json, load_json,
)

logger = setup_logger('report_generator', 'reports/report_generation.log')


def load_all_results() -> Dict:
    """Load all results from the pipeline output directories."""
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

    # Reactor results
    reactor_files = list(REACTOR_DIR.glob("*.json"))
    data['reactor'] = []
    for f in reactor_files:
        try:
            with open(f, 'r') as fh:
                data['reactor'].append(json.load(fh))
        except Exception:
            pass

    # DFT results
    dft_files = list(DFT_DIR.glob("*.json"))
    data['dft'] = []
    for f in dft_files:
        try:
            with open(f, 'r') as fh:
                data['dft'].append(json.load(fh))
        except Exception:
            pass

    # VQE results
    vqe_files = list(VQE_DIR.glob("*.json"))
    data['vqe'] = []
    for f in vqe_files:
        try:
            with open(f, 'r') as fh:
                data['vqe'].append(json.load(fh))
        except Exception:
            pass

    # Fuel cell results
    fc_files = list(FUEL_CELL_DIR.glob("*.json"))
    data['fuel_cell'] = []
    for f in fc_files:
        try:
            with open(f, 'r') as fh:
                data['fuel_cell'].append(json.load(fh))
        except Exception:
            pass

    # Pipeline state
    data['pipeline_state'] = load_json("pipeline_state.json") or {}

    return data


def generate_full_report(pipeline_state: Dict = None) -> Path:
    """
    Generate the comprehensive pipeline report.
    """
    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    data = load_all_results()

    if pipeline_state:
        data['pipeline_state'] = pipeline_state

    report_lines = []
    r = report_lines.append  # shorthand

    # ─── Header ─────────────────────────────────────────────────────────────
    r("# Turquoise Hydrogen → Fuel Cell: Pipeline Report")
    r(f"\n**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    r(f"**Pipeline Version**: 2.0.0")
    r("")

    # ─── Executive Summary ──────────────────────────────────────────────────
    r("## Executive Summary\n")
    ps = data.get('pipeline_state', {})

    p1 = ps.get('phase1', {})
    p2 = ps.get('phase2', {})
    # FC results may come from orchestrator ('phase5') or campaign ('phase5_ga'/'phase5_stack')
    p5 = {**ps.get('phase5', {}), **ps.get('phase5_ga', {}), **ps.get('phase5_stack', {})}

    r(f"| Metric | Value |")
    r(f"|--------|-------|")
    r(f"| Catalysts evaluated (MACE) | {p1.get('total_evaluated', 'N/A'):,} |" if isinstance(p1.get('total_evaluated'), (int, float)) else "| Catalysts evaluated (MACE) | N/A |")
    r(f"| Valid candidates | {p1.get('valid_count', 'N/A'):,} |" if isinstance(p1.get('valid_count'), (int, float)) else "| Valid candidates | N/A |")
    r(f"| Pareto front size | {p1.get('pareto_size', 'N/A')} |")
    r(f"| Best activation barrier | {p1.get('best_E_act', 'N/A')} eV |")
    r(f"| Best coking resistance | {p1.get('best_coking', 'N/A')} |")
    r(f"| Reactor simulations | {p2.get('catalysts_simulated', 'N/A')} |")
    r(f"| Best CH₄ conversion | {p2.get('best_conversion', 'N/A'):.1%} |" if isinstance(p2.get('best_conversion'), (int, float)) else f"| Best CH₄ conversion | N/A |")
    r(f"| FC catalysts evaluated | {p5.get('total_evaluated', p5.get('n_cathodes_screened', 'N/A')):,} |" if isinstance(p5.get('total_evaluated', p5.get('n_cathodes_screened')), (int, float)) else f"| FC catalysts evaluated | N/A |")
    r(f"| Best PEMFC power density | {p5.get('best_power_W_cm2', 'N/A')} W/cm² |" if isinstance(p5.get('best_power_W_cm2'), (int, float)) else f"| Best PEMFC power density | N/A |")
    r(f"| Best PEMFC efficiency | {p5.get('best_efficiency', 'N/A'):.1%} |" if isinstance(p5.get('best_efficiency'), (int, float)) else f"| Best PEMFC efficiency | N/A |")
    r(f"| Min ORR overpotential | {p5.get('min_overpotential_V', 'N/A')} V |" if isinstance(p5.get('min_overpotential_V'), (int, float)) else f"| Min ORR overpotential | N/A |")
    r("")

    # ─── Phase 1: Catalyst Screening ────────────────────────────────────────
    r("## Phase 1: Catalyst Screening & Genetic Optimization\n")

    if 'screening' in data and data['screening']:
        main_db_key = next(
            (k for k in data['screening'] if 'full' in k or 'initial' in k),
            list(data['screening'].keys())[0] if data['screening'] else None
        )
        if main_db_key:
            df = data['screening'][main_db_key]
            valid = df[df.get('valid', True) == True] if 'valid' in df.columns else df

            r(f"### Screening Database: `{main_db_key}`\n")
            r(f"- **Total entries**: {len(df):,}")
            r(f"- **Valid candidates**: {len(valid):,}")

            if 'material_class' in df.columns:
                r(f"\n**Distribution by Material Class:**\n")
                r("| Material Class | Count | Valid |")
                r("|----------------|-------|-------|")
                for cls in df['material_class'].unique():
                    n_total = len(df[df['material_class'] == cls])
                    n_valid = len(valid[valid['material_class'] == cls]) if 'material_class' in valid.columns else 0
                    r(f"| {cls} | {n_total} | {n_valid} |")
                r("")

            if 'E_act' in valid.columns and len(valid) > 0:
                r(f"\n**Top 10 Catalysts by Activation Barrier:**\n")
                r("| Rank | Genome | E_act (eV) | Coking Index | dE_H (eV) |")
                r("|------|--------|------------|-------------|-----------|")
                top = valid.nsmallest(10, 'E_act')
                for rank, (_, row) in enumerate(top.iterrows(), 1):
                    genome_str = str(row.get('genome', ''))[:50]
                    r(f"| {rank} | {genome_str} | {row['E_act']:.4f} | {row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |")
                r("")

    # ─── Phase 2: Reactor Simulation ────────────────────────────────────────
    r("## Phase 2: Reactor-Scale Simulation\n")

    if data['reactor']:
        r("| Catalyst | Reactor | T (K) | CH₄ Conv. | H₂ Select. | τ (s) |")
        r("|----------|---------|-------|-----------|------------|-------|")
        for res in sorted(data['reactor'], key=lambda x: x.get('CH4_conversion', 0), reverse=True)[:20]:
            r(f"| {res.get('catalyst_name', '?')} | {res.get('reactor_type', '?')} | "
              f"{res.get('T_K', '?')} | {res.get('CH4_conversion', 0):.1%} | "
              f"{res.get('H2_selectivity', 'N/A')} | {res.get('residence_time_s', '?'):.1f} |")
        r("")

    # ─── Phase 3: DFT Validation ───────────────────────────────────────────
    r("## Phase 3: DFT Validation (Quantum ESPRESSO)\n")

    if data['dft']:
        r("| Catalyst | Material Class | Energy (eV) | Converged | Max Force |")
        r("|----------|----------------|-------------|-----------|-----------|")
        for res in data['dft']:
            if 'orr_overpotential_V' not in res:  # pyrolysis DFT results
                r(f"| {res.get('catalyst_name', '?')} | {res.get('material_class', '?')} | "
                  f"{res.get('dft_energy_eV', 'N/A')} | {res.get('converged', '?')} | "
                  f"{res.get('max_force_Ry_bohr', 'N/A')} |")
        r("")

    # ─── Phase 4: VQE Results ──────────────────────────────────────────────
    r("## Phase 4: Quantum VQE Transition States\n")

    if data['vqe']:
        r("| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |")
        r("|----------|----------|-------------|-------------|--------|--------|")
        for res in data['vqe']:
            r(f"| {res.get('catalyst_name', '?')} | {res.get('reaction_type', '?')} | "
              f"{res.get('energy_Ha', 0):.6f} | {res.get('energy_eV', 0):.4f} | "
              f"{res.get('n_qubits', '?')} | {res.get('n_layers', '?')} |")
        r("")

    # ─── Phase 5: Fuel Cell Results ────────────────────────────────────────
    r("## Phase 5: PEMFC Fuel Cell Performance\n")

    pemfc_results = [r_ for r_ in data['fuel_cell'] if 'peak_power_W_cm2' in r_ and 'n_cells' not in r_]
    stack_results = [r_ for r_ in data['fuel_cell'] if 'n_cells' in r_]

    if pemfc_results:
        # Sort by the new composite efficiency-overvoltage-power score
        def fc_composite_score(r):
            eff = r.get('efficiency_at_peak', 0.0)
            power = r.get('peak_power_W_cm2', 0.0)
            eta = max(r.get('orr_overpotential_V', 0.4), 0.01)
            return (eff * power) / eta

        r("### Single-Cell Performance\n")
        r("| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |")
        r("|---------|----------|---------------------|---------|-------|--------------------|")
        sorted_pemfc = sorted(pemfc_results, key=fc_composite_score, reverse=True)
        for res in sorted_pemfc[:15]:
            r(f"| {res.get('cathode_catalyst', '?')} | {res.get('membrane', '?')} | "
              f"{res.get('peak_power_W_cm2', 0):.4f} | {res.get('OCV_V', 0):.3f} | "
              f"{res.get('orr_overpotential_V', '?')} | {res.get('efficiency_at_peak', 0):.1%} |")
        r("")

    if stack_results:
        r("### Stack Performance\n")
        for res in stack_results:
            r(f"- **Stack Power**: {res.get('net_power_kW', 0):.1f} kW net ({res.get('n_cells', '?')} cells)")
            r(f"- **System Efficiency**: {res.get('system_efficiency', 0):.1%}")
            r(f"- **Power Density**: {res.get('gravimetric_W_kg', 0):.0f} W/kg, {res.get('volumetric_W_L', 0):.0f} W/L")
            r(f"- **Cost**: ${res.get('cost_per_kW', 0):.0f}/kW")
            r("")

    # ─── ORR DFT Results ───────────────────────────────────────────────────
    orr_dft = [r_ for r_ in data['dft'] if 'orr_overpotential_V' in r_]
    if orr_dft:
        r("### ORR DFT Validation\n")
        r("| Catalyst | η_ORR (V) | dG_OH (eV) | dG_O (eV) | dG_OOH (eV) | RDS |")
        r("|----------|-----------|-----------|----------|------------|-----|")
        for res in orr_dft:
            r(f"| {res.get('catalyst_name', '?')} | {res.get('orr_overpotential_V', '?'):.3f} | "
              f"{res.get('dG_OH_eV', '?')} | {res.get('dG_O_eV', '?')} | "
              f"{res.get('dG_OOH_eV', '?')} | {res.get('rate_determining_step', '?')} |")
        r("")

    # ─── Techno-Economics ──────────────────────────────────────────────────
    r("## Techno-Economic Summary\n")
    r("### Turquoise H₂ Production Cost Drivers\n")
    r("| Parameter | Value |")
    r("|-----------|-------|")
    r("| Catalyst material cost | Based on abundance scoring |")
    r("| Reactor operating temperature | 800–1200 K |")
    r("| Target H₂ production cost | < $2/kg |")
    r("| Carbon co-product credit | $0.10–$0.50/kg |")
    r("")

    # ─── Pipeline Timing ───────────────────────────────────────────────────
    r("## Pipeline Timing\n")
    r("| Phase | Duration |")
    r("|-------|----------|")
    for phase_key in ['phase1', 'phase2', 'phase3', 'phase4', 'phase5', 'phase6']:
        phase_data = ps.get(phase_key, {})
        elapsed = phase_data.get('elapsed_s', 0)
        phase_name = {
            'phase1': 'MACE Screening + GA',
            'phase2': 'Reactor Simulation',
            'phase3': 'DFT Validation',
            'phase4': 'VQE Quantum',
            'phase5': 'Fuel Cell Modeling',
            'phase6': 'Report Generation',
        }.get(phase_key, phase_key)
        r(f"| {phase_name} | {elapsed:.0f}s |")
    total_s = ps.get('total_elapsed_s', 0)
    r(f"| **Total** | **{total_s:.0f}s ({total_s/3600:.1f} hours)** |")
    r("")

    # Write report
    report_content = "\n".join(report_lines)
    report_path = REPORTS_DIR / "pipeline_report.md"
    with open(report_path, 'w') as f:
        f.write(report_content)

    logger.info(f"Report written to {report_path}")
    return report_path


if __name__ == '__main__':
    report = generate_full_report()
    print(f"Report generated: {report}")
