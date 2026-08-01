"""Deterministic screening rankers for small catalyst datasets."""

from __future__ import annotations

import ast
from dataclasses import dataclass

import numpy as np

from pipeline.common.catalyst_spaces import encode_population

TREE_ENSEMBLE_SIZE = 256


def merge_compatible_evidence(current, ledger, protocol_id: str):
    """Merge prior same-protocol observations without inventing compatibility."""
    import pandas as pd

    frames = [current]
    if ledger is not None and len(ledger) and 'screening_protocol' in ledger.columns:
        compatible = ledger[ledger['screening_protocol'] == protocol_id]
        if len(compatible):
            frames.insert(0, compatible)
    merged = pd.concat(frames, ignore_index=True)
    if 'genome' in merged.columns:
        merged = merged.drop_duplicates('genome', keep='last').reset_index(drop=True)
    return merged


@dataclass
class TreeRanker:
    application: str
    model: object
    target_columns: tuple[str, ...]

    def _primary(self, raw) -> np.ndarray:
        raw = np.asarray(raw, float)
        if self.application == "turquoise_hydrogen":
            return raw.reshape(-1)
        d_oh, d_o, d_ooh = raw.T
        # Keep the continuous CHE value: clipping destroys rank information.
        return 1.23 + np.maximum.reduce([
            d_ooh - 4.92, d_o - d_ooh, d_oh - d_o, -d_oh])

    def predict(self, genomes, *, uncertainty: bool = True) -> tuple[np.ndarray, np.ndarray]:
        # sklearn validates and converts X on every individual DecisionTree
        # call.  The encoder already guarantees a finite dense matrix, so do
        # that conversion once instead of hundreds of times per scan batch.
        x = np.asarray(encode_population(genomes), dtype=np.float32, order="C")
        total = np.zeros(len(x), dtype=float)
        total_sq = np.zeros(len(x), dtype=float) if uncertainty else None
        for tree in self.model.estimators_:
            member = self._primary(tree.predict(x, check_input=False))
            total += member
            if total_sq is not None:
                total_sq += member * member
        count = len(self.model.estimators_)
        mean = total / count
        if total_sq is None:
            return mean, np.zeros(len(x), dtype=float)
        # Streaming moments avoid a (batch x trees) allocation in every shard.
        variance = np.maximum(0.0, total_sq / count - mean * mean)
        return mean, np.sqrt(variance)


def fit_tree_ranker(frame, application: str, random_state: int = 20260721) -> TreeRanker:
    """Fit the application-specific form validated by prospective pilots."""
    from sklearn.ensemble import ExtraTreesRegressor
    if application == "turquoise_hydrogen":
        columns = ("E_act",)
    elif application == "fuel_cell_orr":
        columns = ("dG_OH_eV", "dG_O_eV", "dG_OOH_eV")
    else:
        raise ValueError(f"unknown application {application}")
    rows, genomes = [], []
    for _, row in frame.iterrows():
        try:
            values = [float(row[column]) for column in columns]
            genome = ast.literal_eval(row["genome"]) if isinstance(row["genome"], str) else tuple(row["genome"])
            if bool(row.get("valid", True)) and np.all(np.isfinite(values)):
                rows.append(values); genomes.append(genome)
        except (ValueError, TypeError, SyntaxError, KeyError):
            continue
    if len(rows) < 20:
        raise ValueError(f"tree ranker requires at least 20 valid rows; got {len(rows)}")
    y = np.asarray(rows, float)
    if len(columns) == 1:
        y = y[:, 0]
    # With tens to hundreds of calibration rows, 1024 trees add repeated
    # inference work without useful independent evidence.  The fixed 256-tree
    # ensemble retains deterministic uncertainty and ranking while keeping a
    # production-scale indexed scan tractable.
    model = ExtraTreesRegressor(n_estimators=TREE_ENSEMBLE_SIZE, min_samples_leaf=1,
                                max_features=1.0, random_state=random_state, n_jobs=-1)
    model.fit(encode_population(genomes), y)
    return TreeRanker(application, model, columns)


def turquoise_tree_objectives(genomes, ranker: TreeRanker) -> np.ndarray:
    from pipeline.common.utils import abundance_cost_penalty
    from pipeline.screening.genetic_optimizer import _extract_elements_from_genome
    primary, _ = ranker.predict(genomes, uncertainty=False)
    costs = [
        -abundance_cost_penalty(_extract_elements_from_genome(genome))
        for genome in genomes
    ]
    return np.column_stack([primary, np.zeros(len(genomes)), np.zeros(len(genomes)),
                            costs])


def orr_tree_objectives(genomes, ranker: TreeRanker) -> np.ndarray:
    from pipeline.screening.fc_genetic_optimizer import _cost_from_genome, _fenton_from_genome
    from pipeline.common.application_scope import pemfc_cathode_scope
    primary, _ = ranker.predict(genomes, uncertainty=False)
    objectives = np.column_stack([
        primary, [-_fenton_from_genome(g) for g in genomes],
        [_cost_from_genome(g) for g in genomes], np.zeros(len(genomes))])
    for i, genome in enumerate(genomes):
        if pemfc_cathode_scope(genome)["status"] != "candidate":
            objectives[i] = [5.0, 0.0, 100.0, 0.0]
    return objectives
