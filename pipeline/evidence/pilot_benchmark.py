"""Blinded candidate-held-out evidence benchmark for discovery selection.

It compares a genome-only learned ranker with a frozen chemistry heuristic and
random selection. Outcomes are legacy computational screening values, so this
is a pipeline diagnostic rather than experimental discovery evidence.
"""

from __future__ import annotations

import ast
import hashlib
import json
import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from pipeline.common.catalyst_spaces import encode_population


@dataclass(frozen=True)
class PilotSpec:
    application: str
    paths: tuple[str, ...]
    outcome: str
    folds: int = 5
    selection_fraction: float = 0.20
    random_trials: int = 20_000


def _genome(value) -> tuple:
    genome = value if isinstance(value, tuple) else ast.literal_eval(str(value))
    if not isinstance(genome, tuple) or not genome:
        raise ValueError("invalid genome")
    return genome


def _identity(genome: tuple) -> str:
    return hashlib.sha256(repr(genome).encode()).hexdigest()


def load_legacy_outcomes(spec: PilotSpec) -> pd.DataFrame:
    """Load and candidate-deduplicate valid legacy outcomes."""
    frames = [pd.read_csv(path) for path in spec.paths if Path(path).is_file()]
    if not frames:
        raise FileNotFoundError(f"no pilot inputs found for {spec.application}")
    frame = pd.concat(frames, ignore_index=True)
    missing = {"genome", "valid", spec.outcome} - set(frame.columns)
    if missing:
        raise ValueError(f"missing pilot columns: {sorted(missing)}")
    frame = frame[frame.valid.eq(True)].copy()
    frame[spec.outcome] = pd.to_numeric(frame[spec.outcome], errors="coerce")
    frame = frame[np.isfinite(frame[spec.outcome])]
    parsed_rows = []
    for _, row in frame.iterrows():
        try:
            genome = _genome(row.genome)
        except (ValueError, SyntaxError):
            continue
        record = row.to_dict()
        record.update(parsed_genome=genome, candidate_id=_identity(genome))
        parsed_rows.append(record)
    frame = pd.DataFrame(parsed_rows)
    grouped = []
    for _, rows in frame.groupby("candidate_id", sort=True):
        first = rows.iloc[0].copy()
        first[spec.outcome] = float(rows[spec.outcome].median())
        first["replicates"] = int(len(rows))
        grouped.append(first)
    return pd.DataFrame(grouped).reset_index(drop=True)


def _ridge_predict(train: pd.DataFrame, test: pd.DataFrame, outcome: str,
                   alpha: float = 10.0) -> np.ndarray:
    """Deterministic genome-only ridge model, scaled using training data only."""
    x_train = encode_population(train.parsed_genome.tolist()).astype(float)
    x_test = encode_population(test.parsed_genome.tolist()).astype(float)
    mean, scale = x_train.mean(axis=0), x_train.std(axis=0)
    scale[scale < 1e-12] = 1.0
    x_train = np.column_stack([np.ones(len(x_train)), (x_train - mean) / scale])
    x_test = np.column_stack([np.ones(len(x_test)), (x_test - mean) / scale])
    penalty = np.eye(x_train.shape[1]) * alpha
    penalty[0, 0] = 0.0
    beta = np.linalg.pinv(x_train.T @ x_train + penalty) @ x_train.T @ train[outcome].to_numpy(float)
    return x_test @ beta


def _elements(genome: tuple) -> set[str]:
    known = {"Fe", "Co", "Ni", "Cu", "Pt", "Pd", "Ir", "Ru", "Rh", "Ag",
             "Au", "Mo", "W", "V", "Cr", "Mn", "Sn", "Ga", "In"}
    found = set()
    def visit(value):
        if isinstance(value, (tuple, list)):
            for item in value:
                visit(item)
        elif value in known:
            found.add(value)
    visit(genome[1:])
    return found


def expert_score(genome: tuple, application: str) -> float:
    """Frozen outcome-blind conventional-chemistry baseline; lower is better."""
    cls, elems = genome[0], _elements(genome)
    if application == "turquoise_hydrogen":
        prior = {"MoltenMetal": 0.0, "SolidCatalyst": 0.1, "SAA": 0.2,
                 "HEA": 0.3, "MAXPhase": 0.6, "MXene": 0.7,
                 "SAC": 0.8, "DAC": 0.9}.get(cls, 1.5)
        inactive = 0.0 if elems & {"Ni", "Fe", "Co", "Ru", "Rh", "Pt", "Pd"} else 0.8
        costly = 0.35 if elems & {"Pt", "Pd", "Ir", "Rh", "Ru", "Au"} else 0.0
        return prior + inactive + costly
    prior = {"SAC": 0.0, "DAC": 0.1, "SAA": 0.2, "SolidCatalyst": 0.3,
             "MetalFreeCarbon": 0.4, "Spinel": 0.6, "Perovskite": 0.7}.get(cls, 1.3)
    inactive = 0.0 if elems & {"Pt", "Pd", "Ir", "Fe", "Co", "Ni"} else 0.8
    fenton = 0.25 * len(elems & {"Fe", "Cu", "Co", "Mn", "Cr", "V"})
    return prior + inactive + fenton


def _stats(outcomes: np.ndarray, selected: np.ndarray, cutoff: float) -> dict:
    values = outcomes[selected]
    return {"selected": int(len(selected)), "hits": int(np.sum(values <= cutoff)),
            "hit_rate": float(np.mean(values <= cutoff)),
            "best_outcome": float(np.min(values)), "mean_outcome": float(np.mean(values))}


def run_pilot(spec: PilotSpec) -> dict:
    data = load_legacy_outcomes(spec)
    if len(data) < max(20, spec.folds * 3):
        raise ValueError("too few deduplicated valid outcomes for pilot")
    cutoff = float(data[spec.outcome].quantile(0.20))
    membership = np.array([int(cid[:12], 16) % spec.folds for cid in data.candidate_id])
    folds, learned_rows, expert_rows, random_by_fold = [], [], [], []
    rng = np.random.default_rng(20260719)
    for fold in range(spec.folds):
        test, train = data[membership == fold].copy(), data[membership != fold].copy()
        if len(test) < 2 or len(train) < 10:
            continue
        budget = max(1, math.ceil(len(test) * spec.selection_fraction))
        predictions = _ridge_predict(train, test, spec.outcome)
        learned_idx = np.argsort(predictions, kind="stable")[:budget]
        heuristic = np.array([expert_score(g, spec.application) for g in test.parsed_genome])
        expert_idx = np.lexsort((test.candidate_id.to_numpy(), heuristic))[:budget]
        outcomes = test[spec.outcome].to_numpy(float)
        learned, expert = _stats(outcomes, learned_idx, cutoff), _stats(outcomes, expert_idx, cutoff)
        random = [_stats(outcomes, rng.choice(len(test), budget, replace=False), cutoff)
                  for _ in range(spec.random_trials)]
        folds.append({"fold": fold, "train": len(train), "test": len(test), "budget": budget,
                      "learned": learned, "expert": expert,
                      "random_mean_hit_rate": float(np.mean([r["hit_rate"] for r in random]))})
        learned_rows.append(learned); expert_rows.append(expert); random_by_fold.append(random)

    def aggregate(rows):
        selected, hits = sum(r["selected"] for r in rows), sum(r["hits"] for r in rows)
        return {"selected": selected, "hits": hits, "hit_rate": hits / selected,
                "best_outcome": min(r["best_outcome"] for r in rows),
                "mean_outcome": float(np.average([r["mean_outcome"] for r in rows],
                                                  weights=[r["selected"] for r in rows]))}
    learned, expert = aggregate(learned_rows), aggregate(expert_rows)
    # Each random trial receives exactly the same aggregate selection budget as
    # the learned and expert methods across all folds.
    random_rates = np.array([
        sum(fold_rows[trial]["hits"] for fold_rows in random_by_fold) /
        sum(fold_rows[trial]["selected"] for fold_rows in random_by_fold)
        for trial in range(spec.random_trials)
    ])
    random_rate = float(random_rates.mean())
    return {
        "application": spec.application, "evidence_level": "legacy_computational_screening",
        "valid_deduplicated_candidates": len(data), "folds_completed": len(folds),
        "selection_fraction": spec.selection_fraction, "hit_definition": "best 20% by hidden outcome",
        "hit_cutoff": cutoff, "outcome": spec.outcome, "learned": learned, "expert": expert,
        "random": {"trials": len(random_rates), "mean_hit_rate": random_rate,
                   "hit_rate_95pct": [float(np.quantile(random_rates, .025)),
                                      float(np.quantile(random_rates, .975))]},
        "enrichment_vs_random": learned["hit_rate"] / random_rate if random_rate else None,
        "enrichment_vs_expert": learned["hit_rate"] / expert["hit_rate"] if expert["hit_rate"] else None,
        "beats_random": learned["hit_rate"] > float(np.quantile(random_rates, .975)),
        "beats_expert": learned["hit_rate"] > expert["hit_rate"], "folds": folds,
        "limitations": [
            "Outcomes are legacy computational screening values, not DFT, reactor, or MEA measurements.",
            "This is candidate-held-out but not a prospective publication-time split.",
            "The expert baseline is a frozen heuristic, not a panel of human experts.",
            "Small samples make this a pipeline diagnostic, not a discovery-performance claim.",
        ],
    }


def default_specs() -> tuple[PilotSpec, PilotSpec]:
    return (
        PilotSpec("turquoise_hydrogen", (
            "results/screening/ga_initial_screening.csv", "results/screening/ga_fairchem_gen1.csv",
            "results/screening/ga_fairchem_gen2.csv", "results/screening/ga_mace_gen1.csv",
            "results/screening/ga_mace_gen2.csv"), "E_act"),
        PilotSpec("fuel_cell_orr", (
            "results/fuel_cell/cathode_screening.csv", "results/fuel_cell/fc_initial_screening.csv",
            "results/fuel_cell/fc_fairchem_gen1.csv", "results/fuel_cell/fc_mace_gen1.csv"),
                  "orr_overpotential_V"),
    )


def write_report(results: list[dict], output_dir: str = "results/pilot") -> tuple[Path, Path]:
    root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
    json_path, md_path = root / "selection_benchmark.json", root / "selection_benchmark.md"
    json_path.write_text(json.dumps({"schema_version": 1, "results": results}, indent=2) + "\n")
    lines = ["# Pilot selection benchmark", "",
             "Internal computational enrichment test; not experimental validation.", ""]
    for r in results:
        lines += [f"## {r['application']}", "",
                  f"- Deduplicated candidates: {r['valid_deduplicated_candidates']}",
                  f"- Selection budget: {r['learned']['selected']} across {r['folds_completed']} held-out folds",
                  f"- Learned hit rate: {r['learned']['hit_rate']:.1%}",
                  f"- Expert-heuristic hit rate: {r['expert']['hit_rate']:.1%}",
                  f"- Random mean hit rate: {r['random']['mean_hit_rate']:.1%}",
                  f"- Enrichment vs random: {r['enrichment_vs_random']:.2f}x",
                  f"- Beats random 95% bound: {r['beats_random']}",
                  f"- Beats expert heuristic: {r['beats_expert']}", ""]
    md_path.write_text("\n".join(lines) + "\n")
    return json_path, md_path
