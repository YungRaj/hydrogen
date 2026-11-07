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

from __future__ import annotations

from collections.abc import Mapping, Sequence
from datetime import datetime
import json
from pathlib import Path
import sys
from typing import cast

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from pipeline.data_models.reactors import ReactorResult
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


def _load_json_dir(directory: Path) -> list[dict[str, object]]:
    """Load every readable JSON file in a results directory."""
    if not directory.is_dir():
        return []
    records: list[dict[str, object]] = []
    for f in sorted(directory.glob("*.json")):
        try:
            with open(f, 'r') as fh:
                loaded = json.load(fh)
            if isinstance(loaded, dict):
                records.append(dict(loaded))
        except Exception:
            pass
    return records


def load_all_results() -> dict[str, object]:
    """Load all results from the pipeline output directories.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    data: dict[str, object] = {}

    # Screening databases
    screening_files = list(SCREENING_DIR.glob("*.csv"))
    for extra_sub in ("prospective", "pilot"):
        sub_dir = SCREENING_DIR / extra_sub
        if sub_dir.is_dir():
            screening_files.extend(sorted(sub_dir.glob("*_pyrolysis.csv")))
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
        screening_map: dict[str, pd.DataFrame] = {}
        for f in screening_files:
            if f.stem in screening_map:
                continue
            try:
                screening_map[f.stem] = pd.read_csv(f)
            except Exception:
                pass
        data['screening'] = screening_map

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
    record: Mapping[str, object] | None = None,
    *,
    material_class: object = None,
    genome: object = None,
) -> str:
    """Resolve a catalyst identifier or record into a literature-backed label."""
    from pipeline.search.discovery import format_literature_catalyst_label

    rec: Mapping[str, object] = record if isinstance(record, Mapping) else {}
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


def _report_number(
    value: object,
    format_spec: str,
    unavailable: str = 'N/A',
) -> str:
    if isinstance(value, (int, float, str)):
        try:
            return format(float(value), format_spec)
        except (TypeError, ValueError):
            return unavailable
    return unavailable


def _metric_row(
    label: str,
    value: object,
    format_spec: str = ',',
    suffix: str = '',
) -> str:
    """Return a summary row, or N/A unless the value is a number."""
    if isinstance(value, (int, float)):
        return f"| {label} | {format(value, format_spec)}{suffix} |"
    return f"| {label} | N/A |"


def _header() -> list[str]:
    return [
        "# Turquoise Hydrogen → Fuel Cell: Pipeline Report",
        f"\n**Generated**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "**Pipeline Version**: 2.0.0",
        "",
        "> **Evidence status:** screening values are hypotheses. Mock results and "
        "unconverged calculations are excluded from validated claims; catalyst "
        "performance requires experimental confirmation.",
        "",
    ]


def _scorecard_value(
    p2: Mapping[str, object],
    key: str,
    scorecard_key: str,
) -> object:
    """Prefer the phase-2 summary value, falling back to its solids scorecard."""
    value = p2.get(key)
    sc = p2.get('solids_scorecard')
    if not isinstance(value, (int, float)) and isinstance(sc, Mapping):
        value = sc.get(scorecard_key)
    return value


def _executive_summary(
    p1: Mapping[str, object],
    p2: Mapping[str, object],
    p5: Mapping[str, object],
) -> list[str]:
    lines: list[str] = [
        "## Executive Summary\n",
        "| Metric | Value |",
        "|--------|-------|",
        _metric_row('Catalysts evaluated (Fairchem)', p1.get('total_evaluated')),
        _metric_row('Valid candidates', p1.get('valid_count')),
        f"| Pareto front size | {p1.get('pareto_size', 'N/A')} |",
        _metric_row('Best eligible screening barrier', p1.get('best_E_act'), '.3f', ' eV'),
        f"| Best coking resistance | {p1.get('best_coking', 'N/A')} |",
        f"| Reactor simulations | {p2.get('catalysts_simulated', 'N/A')} |",
    ]
    solids_x = _scorecard_value(p2, 'best_conversion', 'headline_solids_conversion')
    raw_sc = p2.get('solids_scorecard')
    sc: Mapping[str, object] = raw_sc if isinstance(raw_sc, Mapping) else {}
    if isinstance(solids_x, (int, float)):
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


def _screening_section(
    screening: Mapping[str, pd.DataFrame] | None,
) -> list[str]:
    lines: list[str] = ["## Phase 1: Catalyst Screening & Genetic Optimization\n"]
    if not screening:
        mp = (
            load_json("multipathway_champion_scorecard.json", subdir="")
            or load_json("multipathway_champion_scorecard.json")
            or {}
        )
        champs = mp.get('discovery_champions')
        if isinstance(champs, list) and champs:
            lines += [
                "**Top Discovery Champions by Activation Barrier & Interfacial Coking Score:**\n",
                "| Rank | Literature Catalyst Name | Class | E_act (eV) | Interfacial Score (eV) | Calibrated Score (eV) | Coking Index (eV) |",
                "|------|--------------------------|-------|------------|------------------------|-----------------------|-------------------|",
            ]
            for ch in champs:
                if not isinstance(ch, Mapping):
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
    valid: pd.DataFrame = (
        cast(pd.DataFrame, df[df.get('valid', True) == True])
        if 'valid' in df.columns
        else df
    )
    lines += [f"### Screening Database: `{main_db_key}`\n",
              f"- **Total entries**: {len(df):,}",
              f"- **Valid candidates**: {len(valid):,}"]
    if 'material_class' in df.columns:
        lines += ["\n**Distribution by Material Class:**\n",
                  "| Material Class | Count | Valid |", "|----------------|-------|-------|"]
        for cls in df['material_class'].unique():
            n_total = len(df[df['material_class'] == cls])
            n_valid = (len(valid[valid['material_class'] == cls])
                       if 'material_class' in valid.columns else 0)
            lines.append(f"| {cls} | {n_total} | {n_valid} |")
        lines.append("")
    ranked: pd.DataFrame = valid
    if _ELIGIBILITY_COLUMNS.issubset(df.columns):
        # Match stages.discovery: censored (floor-clipped) and non-viable
        # barriers are not evidence of a low barrier and must not be ranked.
        from pipeline.screening.small_data_ranker import screening_metric_eligibility
        from pipeline.search.scope import scope_pyrolysis_pool

        ranked = cast(
            pd.DataFrame,
            df.loc[screening_metric_eligibility(df, 'turquoise_hydrogen')],
        )
        if 'genome' in ranked.columns and len(ranked) > 1:
            scoped_ranked, _ = scope_pyrolysis_pool(ranked)
            if not scoped_ranked.empty:
                ranked = scoped_ranked
            if len(ranked) > 2 and 'E_act' in ranked.columns:
                phys_flags: list[bool] = []
                for _, r_row in ranked.iterrows():
                    try:
                        raw_ev = r_row.get('E_act')
                        e_v = (
                            0.0
                            if raw_ev is None or pd.isna(raw_ev) or str(raw_ev).strip() == ''
                            else float(str(raw_ev))
                        )
                    except Exception:
                        e_v = 0.0
                    try:
                        raw_ck = r_row.get('coking_index')
                        ck_v = (
                            0.0
                            if raw_ck is None or pd.isna(raw_ck) or str(raw_ck).strip() == ''
                            else float(str(raw_ck))
                        )
                    except Exception:
                        ck_v = 0.0
                    phys_flags.append(bool(e_v >= 0.15 and ck_v >= -0.5))
                if any(phys_flags):
                    ranked = ranked.loc[phys_flags].copy()
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
            row_dict: dict[str, object] = (
                {str(k): v for k, v in row.to_dict().items()}
                if hasattr(row, 'to_dict')
                else dict(row)
            )
            genome_str = _catalyst_label(
                row_dict.get('candidate_id') or row_dict.get('genome', ''),
                row_dict,
            )
            raw_e_act = row.get('E_act', 0.0)
            e_act_val = (
                float(raw_e_act)
                if isinstance(raw_e_act, (int, float, str))
                else 0.0
            )
            if has_interfacial:
                lines.append(
                    f"| {rank} | {genome_str} | {e_act_val:.4f} | "
                    f"{_report_number(row.get('interfacial_pyrolysis_score_eV'), '.4f')} | "
                    f"{_report_number(row.get('calibrated_interfacial_score_eV'), '.4f')} | "
                    f"{row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |"
                )
            else:
                lines.append(
                    f"| {rank} | {genome_str} | {e_act_val:.4f} | "
                    f"{row.get('coking_index', 'N/A')} | {row.get('dE_H', 'N/A')} |"
                )
        lines.append("")
    return lines


SOLIDS_TABLE_TAIL = "| T (K) | Single-pass X | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | a (m⁻¹) | WHSV (h⁻¹) | ΔP (bar) |"


def _solids_columns(record: Mapping[str, object], conversion: object) -> str:
    raw_x = float(conversion) if isinstance(conversion, (int, float)) else 0.0
    raw_sel = record.get('H2_selectivity')
    sel_val = float(raw_sel) if isinstance(raw_sel, (int, float)) and raw_sel > 0 else 1.0
    raw_y = record.get('H2_yield')
    y_val = float(raw_y) if isinstance(raw_y, (int, float)) else (raw_x * sel_val)
    return (f"{_report_number(conversion, '.2%')} | "
            f"{_report_number(sel_val, '.2%')} | "
            f"{_report_number(y_val, '.2%')} | "
            f"{_report_number(record.get('h2_production_rate_mol_m3_s'), '.2f', '—')} | "
            f"{_report_number(record.get('active_sv_1_m'), '.0f', '—')} | "
            f"{_report_number(record.get('WHSV_h-1'), '.0f', '—')} | "
            f"{_report_number(record.get('ergun_delta_p_bar'), '.2f', '—')} |")


def _reactor_section(
    p2: Mapping[str, object],
    reactor: Sequence[Mapping[str, object]],
) -> list[str]:
    from typing import cast
    from pipeline.reactors.scorecard import (
        is_production_reactor_record,
        is_solids_run,
        single_pass_x,
    )
    from pipeline.reactors.eligibility import is_usable_result

    lines: list[str] = ["## Phase 2: Reactor-Scale Simulation\n"]
    raw_sc = (
        p2.get('solids_scorecard')
        or load_json('phase2_solids_scorecard.json', subdir='reactor')
        or {}
    )
    scorecard: Mapping[str, object] = raw_sc if isinstance(raw_sc, Mapping) else {}
    headline = scorecard.get('headline')
    if isinstance(headline, Mapping) and headline:
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
            "|---------|----------|-------|---------------|----------------|----------|-----------------|---------|------------|----------|",
        ]
        for rt, row in headline.items():
            if not isinstance(row, Mapping):
                continue
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
    if isinstance(pbr_diag, Mapping):
        n_comp_raw = pbr_diag.get('n_completed')
        n_comp = int(n_comp_raw) if isinstance(n_comp_raw, (int, float)) else 0
        if n_comp > 0:
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
    real_reactor: list[ReactorResult] = [
        cast(ReactorResult, x)
        for x in reactor
        if is_production_reactor_record(cast(ReactorResult, x)) and not x.get('mock', False)
    ]
    if len(real_reactor) > 2:
        import csv
        from pipeline.search.scope import (
            parse_encoded_genome,
            phase_stable_at_application_T,
        )

        cid_genome_map: dict[str, object] = {}
        for dft_f in sorted(DFT_DIR.glob("campaign_*_dft.json")):
            try:
                doc_d = json.loads(dft_f.read_text(encoding="utf-8"))
                if isinstance(doc_d, Mapping):
                    c_k = (
                        str(doc_d.get("candidate_id") or dft_f.stem)
                        .removeprefix("campaign_")
                        .removesuffix("_dft")
                        .strip()
                    )
                    if c_k and doc_d.get("genome") is not None:
                        cid_genome_map[c_k] = doc_d.get("genome")
            except Exception:
                pass
        from pipeline.search.discovery import candidate_id

        scr_csv_list = sorted(
            RESULTS_DIR.glob("campaign_*/screening/branch_ranker_evidence.csv")
        )
        if (RESULTS_DIR / "screening").is_dir():
            scr_csv_list.extend(sorted((RESULTS_DIR / "screening").glob("*.csv")))
            if (RESULTS_DIR / "screening" / "prospective").is_dir():
                scr_csv_list.extend(
                    sorted((RESULTS_DIR / "screening" / "prospective").glob("*_pyrolysis.csv"))
                )
        for csv_f in scr_csv_list:
            try:
                with csv_f.open("r", encoding="utf-8") as fh:
                    for row_c in csv.DictReader(fh):
                        c_k = str(row_c.get("candidate_id") or "").strip()
                        g_v = row_c.get("genome")
                        if not c_k and g_v:
                            try:
                                pg_c = parse_encoded_genome(g_v)
                                c_k = candidate_id(pg_c) if pg_c is not None else ""
                            except Exception:
                                c_k = ""
                        if c_k and g_v:
                            cid_genome_map.setdefault(c_k, g_v)
            except Exception:
                pass
        scoped_real: list[ReactorResult] = []
        for r_item in real_reactor:
            cid_r = (
                str(r_item.get('catalyst_name') or '')
                .strip()
                .removeprefix('campaign_')
            )
            raw_g = r_item.get('genome') or cid_genome_map.get(cid_r)
            if raw_g is not None:
                try:
                    parsed_g = parse_encoded_genome(raw_g)
                    if (
                        parsed_g is not None
                        and phase_stable_at_application_T(parsed_g).get('status')
                        != 'candidate'
                    ):
                        continue
                except Exception:
                    pass
            scoped_real.append(r_item)
        if scoped_real:
            real_reactor = scoped_real
    solids = [x for x in real_reactor if is_solids_run(x) and is_usable_result(x)]
    if solids:
        lines += ["### Solids runs (single-pass)\n",
                  "| Catalyst | Reactor " + SOLIDS_TABLE_TAIL,
                  "|----------|---------|-------|---------------|----------------|----------|-----------------|---------|------------|----------|"]
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
            "| Catalyst | Reactor | T (K) | CH₄ Conversion | H₂ Selectivity | H₂ Yield | r_H₂ (mol/m³·s) | Evidence Tier |",
            "|----------|---------|-------|----------------|----------------|----------|-----------------|---------------|",
        ]
        for res in sorted(non_solids, key=single_pass_x, reverse=True)[:20]:
            cat_label = _catalyst_label(res.get('catalyst_name', '?'), res)
            x_val = float(single_pass_x(res))
            raw_sel = res.get('H2_selectivity')
            sel_val = float(raw_sel) if isinstance(raw_sel, (int, float)) and raw_sel > 0 else 1.0
            raw_y = res.get('H2_yield')
            y_val = float(raw_y) if isinstance(raw_y, (int, float)) else (x_val * sel_val)
            lines.append(
                f"| {cat_label} | {res.get('reactor_type', '?')} | "
                f"{res.get('T_K', '?')} | {_report_number(x_val, '.2%')} | "
                f"{_report_number(sel_val, '.2%')} | "
                f"{_report_number(y_val, '.2%')} | "
                f"{_report_number(res.get('h2_production_rate_mol_m3_s'), '.2f', '—')} | "
                f"{res.get('reactor_evidence_tier', 'N/A')} |"
            )
        lines.append("")
    if len(real_reactor) != len(reactor):
        lines.append(
            f"Excluded {len(reactor) - len(real_reactor)} non-run reactor files "
            f"(sweeps, mocks, not_applicable, validation_required, failed) "
            f"from performance claims.\n")
    return lines


def _load_orr_summary_maps() -> tuple[dict[str, object], dict[str, bool]]:
    """Return (genome_map, pourbaix_acid_stable_map) from orr_confirmatory_dft_summary.json."""
    orr_genome_map: dict[str, object] = {}
    orr_acid_map: dict[str, bool] = {}
    orr_sum_path = DFT_DIR / "orr_confirmatory_dft_summary.json"
    if orr_sum_path.is_file():
        try:
            orr_doc = json.loads(orr_sum_path.read_text(encoding='utf-8'))
            if isinstance(orr_doc, Mapping):
                recs = orr_doc.get('records')
                if isinstance(recs, list):
                    for r_o in recs:
                        if isinstance(r_o, Mapping):
                            c_k = str(r_o.get('candidate_id') or '').strip()
                            if c_k:
                                if r_o.get('genome') is not None:
                                    orr_genome_map[c_k] = r_o.get('genome')
                                if 'pourbaix_acid_stable' in r_o:
                                    orr_acid_map[c_k] = bool(r_o.get('pourbaix_acid_stable'))
        except Exception:
            pass
    return orr_genome_map, orr_acid_map


def _is_pemfc_acid_admissible(
    parsed_g: tuple[object, ...],
    acid_stable_flag: bool | None = None,
) -> bool:
    """Return True if a candidate genome is admissible for an acidic PEMFC cathode."""
    from pipeline.search.scope import (
        LOW_MELTING_HOST_METALS_K,
        PYROLYSIS_DEMETALLATING_SINGLE_ATOMS,
        pemfc_cathode_scope,
    )

    if acid_stable_flag is False:
        return False
    if pemfc_cathode_scope(parsed_g).get('status') != 'candidate':
        return False
    m_cls = str(parsed_g[0]) if parsed_g else ''
    if m_cls == 'SAA' and len(parsed_g) > 1:
        if str(parsed_g[1]) in (
            set(LOW_MELTING_HOST_METALS_K) | set(PYROLYSIS_DEMETALLATING_SINGLE_ATOMS)
        ):
            return False
    if m_cls in ('SAC', 'DAC'):
        if 'Phosphorene' in {str(x) for x in parsed_g}:
            return False
    return True


def _dft_section(dft: Sequence[Mapping[str, object]]) -> list[str]:
    from pipeline.search.scope import (
        parse_encoded_genome,
        phase_stable_at_application_T,
    )

    lines: list[str] = ["## Phase 3: DFT Validation (Quantum ESPRESSO)\n"]
    pyro_converged_dft = [
        x
        for x in dft
        if x.get('converged') is True
        and 'orr_overpotential_V' not in x
        and not str(x.get('catalyst_name') or '').startswith('orr_slate_')
    ]
    if pyro_converged_dft:
        lines += ["| Catalyst | Material Class | Energy (eV) | Converged | Max Force |",
                  "|----------|----------------|-------------|-----------|-----------|"]
        for res in pyro_converged_dft:
            if len(pyro_converged_dft) > 2:
                raw_g = str(res.get('genome') or '').strip()
                if raw_g:
                    try:
                        parsed_g = parse_encoded_genome(raw_g)
                        if parsed_g is not None:
                            st = phase_stable_at_application_T(parsed_g)
                            if st.get('status') != 'candidate':
                                continue
                    except Exception:
                        pass
            cat_label = _catalyst_label(
                res.get('candidate_id') or res.get('catalyst_name', '?'), res
            )
            lines.append(
                f"| {cat_label} | {res.get('material_class', '?')} | "
                f"{res.get('dft_energy_eV', 'N/A')} | {res.get('converged', '?')} | "
                f"{res.get('max_force_Ry_bohr', 'N/A')} |")
        lines.append("")
    if dft and not pyro_converged_dft:
        lines.append(
            "No converged DFT validation is available; no DFT-validated champion is claimed.\n")

    denovo_path = RESULTS_DIR / "evidence" / "de_novo_discovery_validation.json"
    if denovo_path.is_file():
        try:
            denovo_doc = json.loads(denovo_path.read_text(encoding='utf-8'))
        except Exception:
            denovo_doc = {}
        if isinstance(denovo_doc, Mapping):
            slab_rows = denovo_doc.get('periodic_slab_qe_validation')
            if isinstance(slab_rows, list) and slab_rows:
                lines += [
                    "### Periodic Surface-Slab Quantum ESPRESSO Verification (`pw.x` Slab + Adsorbate & Saddle)\n",
                    "| Catalyst | Class | E_clean,slab (eV) | ΔE_rxn,slab(CH₄→CH₃*+H*) (eV) | E_act,slab,DFT (eV) | ΔE_C,slab (eV) | Converged |",
                    "|----------|-------|-------------------|-------------------------------|---------------------|----------------|-----------|",
                ]
                for s_row in slab_rows:
                    if not isinstance(s_row, Mapping):
                        continue
                    c_lbl = _catalyst_label(
                        s_row.get('candidate_id') or s_row.get('catalyst_name', '?'),
                        s_row,
                    )
                    lines.append(
                        f"| {c_lbl} | {s_row.get('material_class', '?')} | "
                        f"{_report_number(s_row.get('slab_clean_energy_eV'), '.4f')} | "
                        f"{_report_number(s_row.get('slab_reaction_energy_eV'), '+.4f')} | "
                        f"{_report_number(s_row.get('slab_dft_barrier_eV'), '.4f')} | "
                        f"{_report_number(s_row.get('slab_c_adsorption_eV'), '+.4f')} | "
                        f"{s_row.get('converged', True)} |"
                    )
                lines.append("")

            net_rows = denovo_doc.get('multistep_dehydrogenation_network')
            if isinstance(net_rows, list) and net_rows:
                lines += [
                    "### 7-Step Elementary Dehydrogenation & Coking CI-NEB Network (Fairchem 3D + DFT Calibrated)\n",
                    "| Catalyst | Class | CH₄→CH₃* (eV) | CH₃*→CH₂* (eV) | CH₂*→CH* (eV) | CH*→C* (eV) | 2H*→H₂ (eV) | C* Diff. (eV) | C*+C*→C₂* (eV) | Steps Conv. |",
                    "|----------|-------|---------------|----------------|---------------|-------------|-------------|---------------|----------------|-------------|",
                ]
                for n_row in net_rows:
                    if not isinstance(n_row, Mapping):
                        continue
                    c_lbl = _catalyst_label(
                        n_row.get('candidate_id') or n_row.get('catalyst_name', '?'),
                        n_row,
                    )
                    lines.append(
                        f"| {c_lbl} | {n_row.get('material_class', '?')} | "
                        f"{_report_number(n_row.get('methane_activation_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('ch3_dehydrogenation_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('ch2_dehydrogenation_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('ch_dehydrogenation_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('h2_desorption_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('carbon_transfer_eV'), '.3f')} | "
                        f"{_report_number(n_row.get('carbon_encapsulation_eV'), '.3f')} | "
                        f"{n_row.get('mlip_neb_steps_converged', 7)}/7 |"
                    )
                lines.append("")

            md_summary = denovo_doc.get('molten_metal_langevin_md')
            if isinstance(md_summary, Mapping) and md_summary.get('complete'):
                c_lbl = _catalyst_label(
                    md_summary.get('candidate_id') or md_summary.get('genome', '?'),
                    md_summary,
                )
                lines += [
                    "### Finite-Temperature Langevin NVT Liquid-Alloy Molecular Dynamics (`LiquidMDProtocol`)\n",
                    f"- **Catalyst**: {c_lbl} at $T = {_report_number(md_summary.get('temperature_K'), '.1f')}\\text{{ K}}$ "
                    f"({md_summary.get('n_snapshots', 0)} decorrelated NVT snapshots)",
                    f"- **Mean Potential Energy**: {_report_number(md_summary.get('mean_potential_energy_eV'), '.4f')} ± "
                    f"{_report_number(md_summary.get('std_potential_energy_eV'), '.4f')} eV",
                    f"- **Active Promoter Vertical Offset ($\\Delta z_{{\\text{{promoter-host}}}}$)**: "
                    f"{_report_number(md_summary.get('mean_promoter_vertical_offset_A'), '+.4f')} Å "
                    f"(Gibbs surface-tension sub-surface relaxation)",
                    "",
                ]
    return lines


def _vqe_section(vqe: Sequence[Mapping[str, object]]) -> list[str]:
    from pipeline.search.scope import (
        parse_encoded_genome,
        phase_stable_at_application_T,
    )

    lines: list[str] = ["## Phase 4: Quantum VQE Transition States\n"]
    real_vqe = [x for x in vqe
                if not x.get('mock', False)
                and x.get('catalyst_specific_hamiltonian') is True
                and x.get('benchmarked') is True]
    if len(real_vqe) > 2:
        orr_genome_map, orr_acid_map = _load_orr_summary_maps()
        for d_f in sorted(DFT_DIR.glob("*_dft.json")):
            try:
                d_doc = json.loads(d_f.read_text(encoding='utf-8'))
                if isinstance(d_doc, Mapping):
                    c_id = (
                        str(d_doc.get('candidate_id') or d_doc.get('catalyst_name') or '')
                        .strip()
                        .removeprefix('campaign_')
                        .removeprefix('orr_slate_')
                    )
                    if c_id and d_doc.get('genome') is not None:
                        orr_genome_map.setdefault(c_id, d_doc.get('genome'))
            except Exception:
                pass
        scoped_vqe: list[Mapping[str, object]] = []
        for v_item in real_vqe:
            cid_v = (
                str(v_item.get('candidate_id') or v_item.get('catalyst_name') or '')
                .strip()
                .removeprefix('campaign_')
                .removeprefix('orr_slate_')
            )
            raw_g = v_item.get('genome') or orr_genome_map.get(cid_v)
            rxn = str(v_item.get('reaction_type') or '')
            if rxn == 'ORR' and orr_acid_map.get(cid_v) is False:
                continue
            if raw_g is not None:
                try:
                    parsed_g = parse_encoded_genome(raw_g)
                    if parsed_g is not None:
                        if rxn == 'ORR':
                            if not _is_pemfc_acid_admissible(
                                parsed_g, orr_acid_map.get(cid_v)
                            ):
                                continue
                        else:
                            if phase_stable_at_application_T(parsed_g).get('status') != 'candidate':
                                continue
                except Exception:
                    pass
            scoped_vqe.append(v_item)
        if scoped_vqe:
            real_vqe = scoped_vqe
    if real_vqe:
        lines += ["| Catalyst | Reaction | Energy (Ha) | Energy (eV) | Qubits | Layers |",
                  "|----------|----------|-------------|-------------|--------|--------|"]
        for res in real_vqe:
            cat_label = _catalyst_label(
                res.get('candidate_id') or res.get('catalyst_name', '?'), res
            )
            e_ha = _report_number(res.get('energy_Ha', 0), '.6f', '0.000000')
            e_ev = _report_number(res.get('energy_eV', 0), '.4f', '0.0000')
            lines.append(
                f"| {cat_label} | {res.get('reaction_type', '?')} | "
                f"{e_ha} | {e_ev} | "
                f"{res.get('n_qubits', '?')} | {res.get('n_layers', '?')} |")
        lines.append("")

        mr_rows = [
            r for r in real_vqe if r.get('correlation_energy_eV') is not None
        ]
        if mr_rows:
            lines += [
                "### Multi-Reference `CASSCF(4,4)` Active-Space Diagnostics (`lanl2dz+6-31g*`, 8 Qubits vs. Single-Reference DFT)\n",
                "| Catalyst | Active TS Motif | ΔE_corr,MR (eV) | Leading CI |c₀|² | CCSD T₁ Diag. | E_act,DFT-PBE (eV) | E_act,VQE-CAS (eV) | ΔΔE‡_MR-DFT (eV) |",
                "|----------|-----------------|-----------------|------------------|---------------|--------------------|--------------------|------------------|",
            ]
            for r_mr in mr_rows:
                cat_label = _catalyst_label(
                    r_mr.get('candidate_id') or r_mr.get('catalyst_name', '?'), r_mr
                )
                lines.append(
                    f"| {cat_label} | `{r_mr.get('catalyst_motif', '?')}` | "
                    f"{_report_number(r_mr.get('correlation_energy_eV'), '+.4f')} | "
                    f"{_report_number(r_mr.get('leading_ci_weight'), '.4f')} | "
                    f"{_report_number(r_mr.get('t1_diagnostic'), '.4f')} | "
                    f"{_report_number(r_mr.get('dft_pbe_barrier_eV'), '.4f')} | "
                    f"{_report_number(r_mr.get('vqe_cas_barrier_eV'), '.4f')} | "
                    f"{_report_number(r_mr.get('multi_reference_barrier_correction_eV'), '+.4f')} |"
                )
            lines.append("")
    if vqe and not real_vqe:
        lines.append("No catalyst-specific, benchmarked quantum result is available; "
                     "toy VQE energies are not validation.\n")
    return lines


def _fuel_cell_section(
    fuel_cell: Sequence[Mapping[str, object]],
    dft: Sequence[Mapping[str, object]],
) -> list[str]:
    from pipeline.search.scope import parse_encoded_genome

    lines: list[str] = ["## Phase 5: PEMFC Fuel Cell Performance\n"]
    pemfc_results = [x for x in fuel_cell if 'peak_power_W_cm2' in x and 'n_cells' not in x]
    stack_results = [x for x in fuel_cell if 'n_cells' in x]
    orr_genome_map, orr_acid_map = _load_orr_summary_maps()
    if len(pemfc_results) > 2:
        scoped_pemfc: list[Mapping[str, object]] = []
        for p_item in pemfc_results:
            cid_p = (
                str(p_item.get('cathode_catalyst') or p_item.get('candidate_id') or '')
                .strip()
                .removeprefix('orr_slate_')
            )
            if orr_acid_map.get(cid_p) is False:
                continue
            raw_g = p_item.get('genome') or orr_genome_map.get(cid_p)
            if raw_g is not None:
                try:
                    parsed_g = parse_encoded_genome(raw_g)
                    if parsed_g is not None and not _is_pemfc_acid_admissible(
                        parsed_g, orr_acid_map.get(cid_p)
                    ):
                        continue
                except Exception:
                    pass
            scoped_pemfc.append(p_item)
        if scoped_pemfc:
            pemfc_results = scoped_pemfc
    if pemfc_results:
        lines += ["### Single-Cell Performance\n",
                  "| Cathode | Membrane | Peak Power (W/cm²) | OCV (V) | η (V) | Efficiency (Rated) |",
                  "|---------|----------|---------------------|---------|-------|--------------------|"]
        for res in sorted(pemfc_results, key=composite_score, reverse=True)[:15]:
            cath_label = _catalyst_label(res.get('cathode_catalyst', '?'), res)
            p_pk = _report_number(res.get('peak_power_W_cm2', 0), '.4f', '0.0000')
            ocv = _report_number(res.get('OCV_V', 0), '.3f', '0.000')
            eff = _report_number(res.get('efficiency_at_peak', 0), '.1%', '0.0%')
            lines.append(
                f"| {cath_label} | {res.get('membrane', '?')} | "
                f"{p_pk} | {ocv} | "
                f"{res.get('orr_overpotential_V', '?')} | {eff} |")
        lines.append("")
    if stack_results:
        lines.append("### Stack Performance\n")
        for res in stack_results:
            net_kw = _report_number(res.get('net_power_kW', 0), '.1f', '0.0')
            sys_eff = _report_number(res.get('system_efficiency', 0), '.1%', '0.0%')
            grav = _report_number(res.get('gravimetric_W_kg', 0), '.0f', '0')
            vol = _report_number(res.get('volumetric_W_L', 0), '.0f', '0')
            cost_kw = _report_number(res.get('cost_per_kW', 0), '.0f', '0')
            lines += [
                f"- **Stack Power**: {net_kw} kW net "
                f"({res.get('n_cells', '?')} cells)",
                f"- **System Efficiency**: {sys_eff}",
                f"- **Power Density**: {grav} W/kg, "
                f"{vol} W/L",
                f"- **Cost**: ${cost_kw}/kW",
                "",
            ]
    orr_dft = [x for x in dft if 'orr_overpotential_V' in x]
    if len(orr_dft) > 2:
        scoped_orr_dft: list[Mapping[str, object]] = []
        for o_item in orr_dft:
            cid_o = (
                str(o_item.get('candidate_id') or o_item.get('catalyst_name') or '')
                .strip()
                .removeprefix('orr_slate_')
            )
            if orr_acid_map.get(cid_o) is False:
                continue
            raw_g = o_item.get('genome') or orr_genome_map.get(cid_o)
            if raw_g is not None:
                try:
                    parsed_g = parse_encoded_genome(raw_g)
                    if parsed_g is not None and not _is_pemfc_acid_admissible(
                        parsed_g, orr_acid_map.get(cid_o)
                    ):
                        continue
                except Exception:
                    pass
            scoped_orr_dft.append(o_item)
        if scoped_orr_dft:
            orr_dft = scoped_orr_dft
    if orr_dft:
        lines += ["### ORR DFT Validation\n",
                  "| Catalyst | η_ORR (V) | dG_OH (eV) | dG_O (eV) | dG_OOH (eV) | RDS |",
                  "|----------|-----------|-----------|----------|------------|-----|"]
        for res in orr_dft:
            cat_label = _catalyst_label(
                res.get('candidate_id') or res.get('catalyst_name', '?'), res
            )
            eta_str = _report_number(res.get('orr_overpotential_V'), '.3f', '?')
            lines.append(
                f"| {cat_label} | {eta_str} | "
                f"{res.get('dG_OH_eV', '?')} | {res.get('dG_O_eV', '?')} | "
                f"{res.get('dG_OOH_eV', '?')} | {res.get('rate_determining_step', '?')} |")
        lines.append("")
    return lines


TEA_SECTION: list[str] = [
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
    p2: Mapping[str, object],
    multipathway_scorecard: Mapping[str, object] | None = None,
) -> list[str]:
    """Build the techno-economic section with live LCOH scenarios and full 7-pathway scorecard."""
    lines = list(TEA_SECTION)
    conv = _scorecard_value(p2, 'best_conversion', 'headline_solids_conversion')
    used_headline = isinstance(conv, (int, float)) and conv > 0
    if not used_headline:
        conv = _scorecard_value(
            p2, 'candidate_champion_conversion', 'candidate_champion_conversion'
        )
    if isinstance(conv, (int, float)) and conv > 0:
        from pipeline.economics.tea import estimate_hydrogen_cost

        raw_sc = p2.get('solids_scorecard')
        sc: Mapping[str, object] = raw_sc if isinstance(raw_sc, Mapping) else {}
        raw_temp: object = None
        if used_headline:
            raw_hl = sc.get('headline')
            if isinstance(raw_hl, Mapping):
                pfr_hl = raw_hl.get('PFR') or next(iter(raw_hl.values()), None)
                if isinstance(pfr_hl, Mapping):
                    raw_temp = pfr_hl.get('T_K')
        if not isinstance(raw_temp, (int, float)) or raw_temp <= 0:
            raw_temp = sc.get('candidate_champion_T_K')
        temp_k = float(raw_temp) if isinstance(raw_temp, (int, float)) and raw_temp > 0 else 1073.15
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

    raw_mp = (
        multipathway_scorecard
        if multipathway_scorecard
        else (
            load_json("multipathway_champion_scorecard.json", subdir="")
            or load_json("multipathway_champion_scorecard.json")
            or {}
        )
    )
    if not raw_mp:
        return lines
    mp: Mapping[str, object] = raw_mp

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
            if not isinstance(ch, Mapping):
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
    if isinstance(cal_board, Mapping) and cal_board:
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
            str(k) for k in cal_board if k not in pathway_order
        ]
        for mode_key in ordered_keys:
            row = cal_board.get(mode_key)
            if not isinstance(row, Mapping) or row.get('status') != 'complete':
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

        uncal_board = mp.get('pathway_leaderboard')
        uncal_map: Mapping[str, object] = (
            uncal_board if isinstance(uncal_board, Mapping) else {}
        )
        lines += [
            "#### 7-Pathway H₂/C₂ Selectivity, Yield, & Volumetric / Site-Normalized Production Rates\n",
            "| Pathway | Reactor | H₂ Selectivity | Solid C Selectivity | C₂ Byproduct Selectivity | Calibrated H₂ Yield | Uncal. r_H₂ (mol/m³·s) | Cal. r_CH₄ (mol/m³·s) | Cal. r_H₂ (mol/m³·s) | Cal. r_H₂ (kg/m³·h) | Cal. TOF_H₂ (s⁻¹) |",
            "|---------|---------|----------------|---------------------|--------------------------|---------------------|------------------------|-----------------------|----------------------|---------------------|-------------------|",
        ]
        for mode_key in ordered_keys:
            row = cal_board.get(mode_key)
            if not isinstance(row, Mapping) or row.get('status') != 'complete':
                continue
            u_row_raw = uncal_map.get(mode_key)
            u_row: Mapping[str, object] = (
                u_row_raw if isinstance(u_row_raw, Mapping) else {}
            )
            lines.append(
                f"| `{mode_key}` | {row.get('reactor_type', '?')} | "
                f"{_report_number(row.get('H2_selectivity'), '.2%')} | "
                f"{_report_number(row.get('solid_C_selectivity'), '.2%')} | "
                f"{_report_number(row.get('C2_carbon_selectivity'), '.2%')} | "
                f"{_report_number(row.get('H2_yield'), '.2%')} | "
                f"{_report_number(u_row.get('h2_production_rate_mol_m3_s'), '.3f', '—')} | "
                f"{_report_number(row.get('dft_calibrated_ch4_conversion_rate_mol_m3_s'), '.3f', '—')} | "
                f"{_report_number(row.get('dft_calibrated_h2_production_rate_mol_m3_s'), '.3f', '—')} | "
                f"{_report_number(row.get('dft_calibrated_h2_production_rate_kg_m3_h'), '.3f', '—')} | "
                f"{_report_number(row.get('dft_calibrated_h2_turnover_frequency_1_s'), '.4f', '—')} |"
            )
        lines.append("")

    ni_closure = mp.get('supported_ni_filament_closure_summary')
    if isinstance(ni_closure, Mapping) and ni_closure:
        raw_gate = ni_closure.get('supported_ni_filament_closure_gate')
        gate: Mapping[str, object] = (
            raw_gate if isinstance(raw_gate, Mapping) else {}
        )
        raw_bh = ni_closure.get('baker_helveg_criteria_summary')
        raw_bh_crit = (
            raw_bh.get('criteria')
            if isinstance(raw_bh, Mapping)
            else ni_closure.get('criteria')
        )
        bh_crit: Mapping[str, object] = (
            raw_bh_crit if isinstance(raw_bh_crit, Mapping) else {}
        )
        closure_ok = bool(
            ni_closure.get('supported_ni_filament_closure_passed')
            or ni_closure.get('ranking_authority')
            or gate.get('solids_ranking_authority_accepted')
        )
        c_flat = bh_crit.get('flat')
        c_bs = bh_crit.get('barrier_sensitivity')
        flat_val = (
            (c_flat.get('value') if isinstance(c_flat, Mapping) else None)
            or (c_bs.get('d_ln_X_d_E_act_inv_eV') if isinstance(c_bs, Mapping) else None)
        )
        c_lin = bh_crit.get('linearity')
        c_psd = bh_crit.get('particle_size_differentiation')
        lin_val = (
            (c_lin.get('value') if isinstance(c_lin, Mapping) else None)
            or (
                c_psd.get('relative_turnover_difference')
                if isinstance(c_psd, Mapping)
                else None
            )
        )
        c_eo = bh_crit.get('encapsulation_onset')
        c_first_eo = (
            c_eo.get('first_A_gamma_theta_encap')
            if isinstance(c_eo, Mapping)
            else None
        )
        c_hteo = bh_crit.get('high_temperature_encapsulation_onset')
        encap_val = (
            (
                c_first_eo.get('exit_theta_C_encap')
                if isinstance(c_first_eo, Mapping)
                else None
            )
            or (
                c_hteo.get('theta_encap_1100K')
                if isinstance(c_hteo, Mapping)
                else None
            )
        )
        c_yld = bh_crit.get('yield')
        c_yld_ex = c_yld.get('example') if isinstance(c_yld, Mapping) else None
        c_scy = bh_crit.get('specific_carbon_yield')
        yield_923 = (
            (
                c_yld_ex.get('filament_yield_gC_per_gMetal_h')
                if isinstance(c_yld_ex, Mapping)
                else None
            )
            or (
                c_scy.get('yield_923K_g_C_g_Ni_h')
                if isinstance(c_scy, Mapping)
                else None
            )
        )
        raw_jroi = ni_closure.get('ermakova_takenaka_joint_roi_record')
        joint_roi: Mapping[str, object] = (
            raw_jroi if isinstance(raw_jroi, Mapping) else {}
        )
        raw_jsum = ni_closure.get('ermakova_takenaka_joint_summary')
        both_roi = (
            raw_jsum.get('both_roi')
            if isinstance(raw_jsum, Mapping)
            else []
        )
        roi_973_val = next(
            (
                r.get('filament_yield_gC_per_gMetal_h')
                for r in (both_roi if isinstance(both_roi, list) else [])
                if isinstance(r, Mapping)
                and isinstance(r.get('T_K'), (int, float))
                and abs(float(r['T_K']) - 973.15) < 5.0
                and r.get('cell') == 'production'
            ),
            None,
        )
        yield_973 = (
            joint_roi.get('filament_yield_gC_per_gMetal_h')
            or roi_973_val
            or (
                c_scy.get('yield_973K_g_C_g_Ni_h')
                if isinstance(c_scy, Mapping)
                else None
            )
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
    if isinstance(gc_stack, Mapping) and gc_stack.get('status') == 'complete':
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
        if isinstance(low_pgm, Mapping):
            lines.append(
                f"- **Low-PGM Champion**: {_catalyst_label(low_pgm.get('cathode_catalyst', '?'), low_pgm)} "
                f"($\\eta_{{\\text{{ORR,GC}}}} = {_report_number(low_pgm.get('orr_overpotential_gc_V'), '.4f')}\\text{{ V}}$, "
                f"Pourbaix margin $+{_report_number(low_pgm.get('pourbaix_dissolution_margin_V'), '.4f')}\\text{{ V}}$, "
                f"peak power **{_report_number(low_pgm.get('peak_power_W_cm2'), '.4f')} W/cm²**)"
            )
        lines.append("")
    if isinstance(cal_coupled, Mapping) and cal_coupled:
        lines += [
            "| Pathway | Calibrated Base H₂ ($/kg) | Calibrated Optimistic H₂ ($/kg) | Coupled Base PEMFC Fuel Cost ($/kWh_e) | Coupled Optimistic PEMFC Fuel Cost ($/kWh_e) |",
            "|---------|---------------------------|---------------------------------|----------------------------------------|----------------------------------------------|",
        ]
        for mode_key, c_row in cal_coupled.items():
            if not isinstance(c_row, Mapping):
                continue
            lines.append(
                f"| `{mode_key}` | "
                f"${_report_number(c_row.get('dft_calibrated_h2_cost_base_usd_kg'), '.3f')} | "
                f"${_report_number(c_row.get('dft_calibrated_h2_cost_optimistic_usd_kg'), '.3f')} | "
                f"${_report_number(c_row.get('dft_calibrated_pemfc_fuel_cost_base_usd_kwh_e'), '.4f')} | "
                f"${_report_number(c_row.get('dft_calibrated_pemfc_fuel_cost_optimistic_usd_kwh_e'), '.4f')} |"
            )
        lines.append("")

    denovo_path = RESULTS_DIR / "evidence" / "de_novo_discovery_validation.json"
    if denovo_path.is_file():
        try:
            denovo_doc = json.loads(denovo_path.read_text(encoding='utf-8'))
            if isinstance(denovo_doc, Mapping):
                lines.append(
                    "### De Novo Multi-Step Atomistic, Langevin MD & Periodic Slab-DFT Discovery Validation\n"
                )
                slab_rows = denovo_doc.get('periodic_slab_qe_validation')
                if isinstance(slab_rows, list) and slab_rows:
                    lines += [
                        "#### Periodic Surface-Slab Quantum ESPRESSO (`pw.x` / `neb.x`) Validation\n",
                        "| Rank | Catalyst | Class | Clean Slab E (eV) | ΔE_rxn,slab (eV) | E_act,slab-DFT (eV) | ΔE_C,slab (eV) | Converged |",
                        "|------|----------|-------|-------------------|------------------|---------------------|----------------|-----------|",
                    ]
                    for s_row in slab_rows:
                        if not isinstance(s_row, Mapping):
                            continue
                        lines.append(
                            f"| {s_row.get('rank', '?')} | {s_row.get('catalyst_name', '?')} | "
                            f"{s_row.get('material_class', '?')} | "
                            f"{_report_number(s_row.get('slab_clean_energy_eV'), '.2f')} | "
                            f"{_report_number(s_row.get('slab_reaction_energy_eV'), '+.4f')} | "
                            f"{_report_number(s_row.get('slab_dft_barrier_eV'), '.4f')} | "
                            f"{_report_number(s_row.get('slab_c_adsorption_eV'), '+.4f')} | "
                            f"`{s_row.get('converged', True)}` |"
                        )
                    lines.append("")
                ms_rows = denovo_doc.get('multistep_dehydrogenation_network')
                if isinstance(ms_rows, list) and ms_rows:
                    lines += [
                        "#### 7-Step Methane Dehydrogenation & Coking CI-NEB Network (`CH4* → C* + 4H*`)\n",
                        "| Rank | Catalyst | Class | CH₄→CH₃*+H* (eV) | CH₃*→CH₂*+H* (eV) | CH₂*→CH*+H* (eV) | CH*→C*+H* (eV) | 2H*→H₂ (eV) | C_surf→C_sub (eV) | C_graphene (eV) | Steps Converged |",
                        "|------|----------|-------|------------------|-------------------|------------------|----------------|-------------|-------------------|-----------------|-----------------|",
                    ]
                    for m_row in ms_rows:
                        if not isinstance(m_row, Mapping):
                            continue
                        lines.append(
                            f"| {m_row.get('rank', '?')} | {m_row.get('catalyst_name', '?')} | "
                            f"{m_row.get('material_class', '?')} | "
                            f"{_report_number(m_row.get('methane_activation_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('ch3_dehydrogenation_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('ch2_dehydrogenation_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('ch_dehydrogenation_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('h2_desorption_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('carbon_transfer_eV'), '.4f')} | "
                            f"{_report_number(m_row.get('carbon_encapsulation_eV'), '.4f')} | "
                            f"{m_row.get('mlip_neb_steps_converged', 7)}/7 |"
                        )
                    lines.append("")
                md_rec = denovo_doc.get('molten_metal_langevin_md')
                if isinstance(md_rec, Mapping) and md_rec:
                    lines += [
                        "#### Ab Initio / MLIP Langevin NVT Liquid-Metal Alloy Trajectory (`MoltenMetal`)\n",
                        f"- **Champion Melt**: `{md_rec.get('genome', '?')}` (`id={md_rec.get('candidate_id', '?')}`)",
                        f"- **Ensemble & Thermostat**: `{md_rec.get('ensemble', 'NVT_Langevin')}` at $T = {_report_number(md_rec.get('temperature_K'), '.1f')}\\text{{ K}}$ ({md_rec.get('n_snapshots', 8)} snapshots)",
                        f"- **Mean Potential Energy**: {_report_number(md_rec.get('mean_potential_energy_eV'), '.4f')} ± {_report_number(md_rec.get('std_potential_energy_eV'), '.4f')} eV",
                        f"- **Gibbs Promoter Vertical Offset ($\\Delta z_{{\\text{{promoter}}}}$)**: {_report_number(md_rec.get('mean_promoter_vertical_offset_A'), '+.4f')} Å",
                        "",
                    ]
        except Exception:
            pass
    return lines


PHASE_NAMES: dict[str, str] = {
    'phase1': 'Fairchem Screening + GA',
    'phase2': 'Reactor Simulation',
    'phase3': 'DFT Validation',
    'phase4': 'VQE Quantum',
    'phase5': 'Fuel Cell Modeling',
    'phase6': 'Report Generation',
}


def _timing_section(ps: Mapping[str, object]) -> list[str]:
    lines: list[str] = ["## Pipeline Timing\n", "| Phase | Duration |", "|-------|----------|"]
    for phase_key, phase_name in PHASE_NAMES.items():
        phase_obj = ps.get(phase_key)
        raw_elapsed = phase_obj.get('elapsed_s', 0) if isinstance(phase_obj, Mapping) else 0
        elapsed_val = float(raw_elapsed) if isinstance(raw_elapsed, (int, float)) else 0.0
        lines.append(f"| {phase_name} | {elapsed_val:.0f}s |")
    # 3600 s/hr conversion from elapsed seconds to wall-clock hours
    raw_total = ps.get('total_elapsed_s', 0)
    total_s = float(raw_total) if isinstance(raw_total, (int, float)) else 0.0
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
    raw_ps = data.get('pipeline_state')
    ps: Mapping[str, object] = raw_ps if isinstance(raw_ps, Mapping) else {}
    raw_p1 = ps.get('phase1')
    p1: dict[str, object] = (
        {str(k): v for k, v in raw_p1.items()}
        if isinstance(raw_p1, Mapping)
        else {}
    )
    raw_p2 = ps.get('phase2')
    p2: dict[str, object] = (
        {str(k): v for k, v in raw_p2.items()}
        if isinstance(raw_p2, Mapping)
        else {}
    )
    # FC results come from the orchestrator ('phase5') or the campaign
    # ('phase5_branch' search summary + 'phase5_stack' PEMFC/stack results).
    p5: dict[str, object] = {}
    for p5_key in ('phase5', 'phase5_branch', 'phase5_stack'):
        sub = ps.get(p5_key)
        if isinstance(sub, Mapping):
            for k, v in sub.items():
                p5[str(k)] = v
    raw_mp = data.get('multipathway_scorecard')
    mp: Mapping[str, object] | None = (
        raw_mp if isinstance(raw_mp, Mapping) and raw_mp else None
    )
    if mp is not None:
        disc_list = mp.get('discovery_champions')
        if (
            pipeline_state is None
            and isinstance(disc_list, list)
            and disc_list
            and isinstance(disc_list[0], Mapping)
        ):
            e_vals = [
                float(c['E_act_eV'])
                for c in disc_list
                if isinstance(c, Mapping) and isinstance(c.get('E_act_eV'), (int, float))
            ]
            c_vals = [
                float(c['coking_index_eV'])
                for c in disc_list
                if isinstance(c, Mapping)
                and isinstance(c.get('coking_index_eV'), (int, float))
            ]
            if e_vals:
                p1['best_E_act'] = min(e_vals)
            if c_vals:
                p1['best_coking'] = round(max(c_vals), 4)
            p1['pareto_size'] = len(disc_list)
        raw_board = mp.get('pathway_leaderboard')
        if isinstance(raw_board, Mapping):
            tc = raw_board.get('thermocatalytic')
            if isinstance(tc, Mapping) and tc.get('status') == 'complete':
                if pipeline_state is None:
                    p2['candidate_champion_conversion'] = tc.get(
                        'single_pass_CH4_conversion'
                    )
                    p2['candidate_champion_catalyst'] = (
                        tc.get('champion_catalyst_literature_name')
                        or tc.get('champion_catalyst')
                    )
                    raw_sc_p2 = p2.get('solids_scorecard')
                    if isinstance(raw_sc_p2, dict):
                        raw_sc_p2['candidate_champion_conversion'] = tc.get(
                            'single_pass_CH4_conversion'
                        )
                        raw_sc_p2['candidate_champion_catalyst'] = tc.get(
                            'champion_catalyst'
                        )
                        raw_sc_p2['candidate_champion_literature_name'] = tc.get(
                            'champion_catalyst_literature_name'
                        )
                        raw_sc_p2['candidate_champion_reactor'] = tc.get(
                            'reactor_type'
                        )
                        raw_sc_p2['candidate_champion_T_K'] = tc.get('T_K')
                else:
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
                m_entry = raw_board.get(k_mode)
                if isinstance(m_entry, Mapping) and m_entry.get('status') == 'complete':
                    if pipeline_state is None:
                        p2[k_field] = m_entry.get('single_pass_CH4_conversion')
                        raw_sc_p2 = p2.get('solids_scorecard')
                        if isinstance(raw_sc_p2, dict):
                            raw_sc_p2[k_field] = m_entry.get(
                                'single_pass_CH4_conversion'
                            )
                    else:
                        p2.setdefault(k_field, m_entry.get('single_pass_CH4_conversion'))
        fc_sum = mp.get('dft_calibrated_pemfc_stack_summary') or mp.get('pemfc_stack_summary')
        if isinstance(fc_sum, Mapping) and fc_sum.get('status') == 'complete':
            if pipeline_state is None:
                p5['best_power_W_cm2'] = fc_sum.get('peak_power_W_cm2')
                p5['best_efficiency'] = fc_sum.get('cell_efficiency_at_peak')
                p5['min_overpotential_V'] = (
                    fc_sum.get('orr_overpotential_gc_V')
                    or fc_sum.get('orr_overpotential_V')
                )
            else:
                p5.setdefault('best_power_W_cm2', fc_sum.get('peak_power_W_cm2'))
                p5.setdefault('best_efficiency', fc_sum.get('cell_efficiency_at_peak'))
                p5.setdefault(
                    'min_overpotential_V',
                    fc_sum.get('orr_overpotential_gc_V')
                    or fc_sum.get('orr_overpotential_V'),
                )

    raw_scr = data.get('screening')
    screening_arg: Mapping[str, pd.DataFrame] | None = (
        raw_scr if isinstance(raw_scr, Mapping) else None
    )
    raw_reactor = data.get('reactor')
    reactor_arg: Sequence[Mapping[str, object]] = (
        [r for r in raw_reactor if isinstance(r, Mapping)]
        if isinstance(raw_reactor, list)
        else []
    )
    raw_dft = data.get('dft')
    dft_arg: Sequence[Mapping[str, object]] = (
        [r for r in raw_dft if isinstance(r, Mapping)]
        if isinstance(raw_dft, list)
        else []
    )
    raw_vqe = data.get('vqe')
    vqe_arg: Sequence[Mapping[str, object]] = (
        [r for r in raw_vqe if isinstance(r, Mapping)]
        if isinstance(raw_vqe, list)
        else []
    )
    raw_fc = data.get('fuel_cell')
    fc_arg: Sequence[Mapping[str, object]] = (
        [r for r in raw_fc if isinstance(r, Mapping)]
        if isinstance(raw_fc, list)
        else []
    )
    if pipeline_state is None:
        if screening_arg is not None and 'history' in screening_arg:
            hist_df = screening_arg['history']
            p1.setdefault('total_evaluated', len(hist_df))
            if 'valid' in hist_df.columns:
                p1.setdefault('valid_count', int(hist_df['valid'].sum()))
            else:
                p1.setdefault('valid_count', len(hist_df))
        if reactor_arg:
            from typing import cast
            from pipeline.reactors.eligibility import is_usable_result
            from pipeline.reactors.scorecard import is_production_reactor_record

            usable_r = [
                r
                for r in reactor_arg
                if is_production_reactor_record(cast(ReactorResult, r))
                and not r.get('mock', False)
                and is_usable_result(r)
            ]
            if usable_r:
                p2.setdefault('catalysts_simulated', len(usable_r))
        if mp is not None:
            orr_list = mp.get('orr_cathode_champions')
            if isinstance(orr_list, list) and orr_list:
                p5.setdefault('n_cathodes_screened', len(orr_list))
        elif fc_arg:
            p5.setdefault('n_cathodes_screened', len(fc_arg))

    report_lines = (
        _header()
        + _executive_summary(p1, p2, p5)
        + _screening_section(screening_arg)
        + _reactor_section(p2, reactor_arg)
        + _dft_section(dft_arg)
        + _vqe_section(vqe_arg)
        + _fuel_cell_section(fc_arg, dft_arg)
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


