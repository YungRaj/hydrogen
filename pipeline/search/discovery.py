"""Deterministic, coverage-aware search candidate acquisition.

The evolutionary optimizers are useful local exploiters, but a GA alone cannot
show which parts of a combinatorial space were never visited.  This module
provides canonical candidate identities, chemistry-region coverage accounting,
and a deterministic quality/novelty acquisition rule shared by both campaigns.
"""

from __future__ import annotations

import hashlib
import json
import ast
from collections import Counter
from typing import Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


def canonicalize_genome(genome: tuple) -> tuple:
    """Return a stable representation without changing site semantics."""
    g = list(genome)
    if not g:
        raise ValueError("empty catalyst genome")
    if g[0] == "SolidCatalyst":
        # Dopant order is not represented by the generated slab; permutations
        # therefore describe the same candidate.  Preserve multiplicity.
        g[5] = tuple(sorted(g[5]))
        g[4] = round(float(g[4]), 6)
    elif g[0] == "MoltenMetal":
        # Zero loading and self-promotion both describe the pure host.
        if g[2] == "None" or float(g[3]) == 0.0 or g[2] == g[1]:
            g[2], g[3] = "None", 0.0
    elif g[0] == "Perovskite":
        # Dopant identity has no physical meaning at zero fraction.
        if float(g[4]) == 0.0:
            g[3] = "None"
    elif g[0] == "HEA":
        g[1] = tuple(sorted(g[1]))
    elif g[0] == "MetalHydride" and g[3] == g[1]:
        # A repeated secondary metal is the single-metal composition.
        g[3] = "None"
    elif g[0] == "MAXPhase" and g[5] == g[1]:
        g[5] = "None"
    elif g[0] == "MXene" and g[5] == g[1]:
        g[5] = "None"
    return tuple(g)


def candidate_id(genome: tuple) -> str:
    """Content-addressed ID suitable for resumable, sharded campaigns."""
    payload = json.dumps(canonicalize_genome(genome), separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()[:20]


def discovery_region(genome: tuple) -> Tuple[str, ...]:
    """A coarse chemistry cell used to prevent global-rank mode collapse."""
    g = canonicalize_genome(genome)
    cls = g[0]
    if cls == "SolidCatalyst":
        return cls, str(g[1]), str(g[2]), str(g[3])
    if cls == "MoltenMetal":
        return cls, str(g[1]), str(g[2])
    if cls in ("SAC", "DAC"):
        return tuple([cls] + [str(x) for x in g[1:4]])
    if cls in ("MOF", "COF"):
        return cls, str(g[1]), str(g[2]), str(g[3])
    if cls == "HEA":
        return cls, "-".join(g[1]), str(g[2])
    return tuple([cls] + [str(x) for x in g[1:3]])


def coverage_summary(genomes: Iterable[tuple]) -> dict:
    """Return JSON-serializable class and chemistry-cell coverage."""
    canonical = {candidate_id(g): canonicalize_genome(g) for g in genomes}
    class_counts = Counter(g[0] for g in canonical.values())
    regions = {discovery_region(g) for g in canonical.values()}
    return {
        "unique_candidates": len(canonical),
        "unique_regions": len(regions),
        "class_counts": dict(sorted(class_counts.items())),
    }


def add_discovery_metadata(frame):
    """Add stable identity and coverage-cell columns to a screening DataFrame."""
    if frame is None or "genome" not in frame.columns:
        return frame
    ids, regions = [], []
    for raw in frame["genome"]:
        try:
            genome = ast.literal_eval(raw) if isinstance(raw, str) else tuple(raw)
            ids.append(candidate_id(genome))
            regions.append("|".join(discovery_region(genome)))
        except (ValueError, SyntaxError, TypeError):
            ids.append(None)
            regions.append(None)
    frame = frame.copy()
    frame["candidate_id"] = ids
    frame["discovery_region"] = regions
    return frame


def _quality_score(objectives: np.ndarray) -> np.ndarray:
    """Robustly scale minimization objectives and combine without units."""
    obj = np.asarray(objectives, dtype=float)
    if obj.ndim != 2 or len(obj) == 0:
        return np.zeros(len(obj), dtype=float)
    finite = np.where(np.isfinite(obj), obj, np.nan)
    lo = np.nanpercentile(finite, 5, axis=0)
    hi = np.nanpercentile(finite, 95, axis=0)
    scaled = np.clip((finite - lo) / np.maximum(hi - lo, 1e-12), 0.0, 1.0)
    scaled = np.nan_to_num(scaled, nan=1.0, posinf=1.0, neginf=0.0)
    return 1.0 - scaled.mean(axis=1)


def select_discovery_batch(
    candidates: Sequence[tuple],
    objectives: np.ndarray,
    n_select: int,
    evaluated: Iterable[tuple] = (),
    uncertainties: Optional[Sequence[float]] = None,
    confidence: Optional[Sequence[float]] = None,
) -> List[int]:
    """Select a deterministic batch balancing viability, novelty, and OOD value.

    First take the best candidate from as many unseen chemistry cells as the
    budget permits.  Remaining slots use a combined score.  Low confidence is
    treated as a reason to *validate* a candidate, not as proof it is poor.
    """
    if n_select <= 0 or not candidates:
        return []
    n_select = min(n_select, len(candidates))
    quality = _quality_score(objectives)
    seen_ids = {candidate_id(g) for g in evaluated}
    seen_regions = {discovery_region(g) for g in evaluated}
    uncertainty = np.zeros(len(candidates)) if uncertainties is None else np.asarray(uncertainties, float)
    if len(uncertainty) != len(candidates):
        raise ValueError("uncertainties must match candidates")
    uncertainty_span = np.ptp(uncertainty)
    if uncertainty_span > 0:
        uncertainty = (uncertainty - uncertainty.min()) / uncertainty_span
    conf = np.ones(len(candidates)) if confidence is None else np.asarray(confidence, float)
    if len(conf) != len(candidates):
        raise ValueError("confidence must match candidates")

    ids = [candidate_id(g) for g in candidates]
    regions = [discovery_region(g) for g in candidates]
    eligible = [i for i, cid in enumerate(ids) if cid not in seen_ids]
    by_region = {}
    for i in eligible:
        by_region.setdefault(regions[i], []).append(i)

    # Region champions: quality remains important, with deterministic ID ties.
    champions = [
        max(indices, key=lambda i: (quality[i], uncertainty[i], ids[i]))
        for region, indices in sorted(by_region.items()) if region not in seen_regions
    ]
    champions.sort(key=lambda i: (-quality[i], -uncertainty[i], ids[i]))
    selected = champions[:n_select]
    selected_set = set(selected)

    # OOD value is capped: unfamiliar chemistry receives validation opportunity
    # but cannot swamp a batch solely because the model knows little about it.
    score = 0.60 * quality + 0.20 * uncertainty + 0.20 * np.clip(1.0 - conf, 0.0, 1.0)
    remainder = sorted(
        (i for i in eligible if i not in selected_set),
        key=lambda i: (-score[i], ids[i]),
    )
    selected.extend(remainder[: n_select - len(selected)])
    return selected
