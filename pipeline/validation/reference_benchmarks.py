"""Fail-closed validation of physical solver outputs against literature.

This module deliberately separates numerical verification (for example VQE
against exact diagonalization) from physical validation against experiment or
published electronic-structure calculations.  Only derived observables with
matching identity, units, conditions/protocol, and evidence quality may be
compared.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
from typing import Any

from pipeline.data_models.reference_validation import (
    ComputedObservation,
    ReferenceComparison,
    ReferenceObservation,
)


def _required_text(record: dict[str, Any], key: str) -> str:
    value = record.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"reference {key!r} must be non-empty text")
    return value.strip()


def load_reference_manifest(path: str | Path) -> dict[str, ReferenceObservation]:
    """Load and validate a curated quantum-reference manifest.

    Args:
        path: JSON manifest containing literature reference observations.

    Returns:
        References keyed by their stable reference IDs.
    """
    data = json.loads(Path(path).read_text())
    if data.get("schema_version") != 1:
        raise ValueError("reference manifest schema_version must be 1")
    records = data.get("references")
    if not isinstance(records, list) or not records:
        raise ValueError("reference manifest must contain references")

    references: dict[str, ReferenceObservation] = {}
    for raw in records:
        if not isinstance(raw, dict):
            raise ValueError("each reference must be an object")
        reference_id = _required_text(raw, "reference_id")
        if reference_id in references:
            raise ValueError(f"duplicate reference_id: {reference_id}")
        kind = raw.get("reference_kind")
        if kind not in {"experiment", "published_computation"}:
            raise ValueError(f"invalid reference_kind for {reference_id}")
        solver = raw.get("solver_family")
        if solver not in {"quantum_espresso", "cudaq_vqe"}:
            raise ValueError(f"invalid solver_family for {reference_id}")
        conditions = raw.get("conditions")
        if not isinstance(conditions, dict) or not conditions:
            raise ValueError(f"conditions are required for {reference_id}")
        value = float(raw["value"])
        tolerance = float(raw["tolerance"])
        if not math.isfinite(value) or not math.isfinite(tolerance) or tolerance <= 0:
            raise ValueError(f"finite value and positive tolerance required for {reference_id}")
        references[reference_id] = ReferenceObservation(
            reference_id=reference_id,
            material_id=_required_text(raw, "material_id"),
            observable=_required_text(raw, "observable"),
            value=value,
            unit=_required_text(raw, "unit"),
            tolerance=tolerance,
            reference_kind=kind,
            solver_family=solver,
            citation=_required_text(raw, "citation"),
            source_url=_required_text(raw, "source_url"),
            conditions=conditions,
            protocol=_required_text(raw, "protocol"),
            tolerance_basis=_required_text(raw, "tolerance_basis"),
        )
    return references


def computed_observation(record: dict[str, Any]) -> ComputedObservation:
    """Parse a solver-derived observable, rejecting incomplete provenance.

    Args:
        record: Serialized physical observable and its solver evidence flags.

    Returns:
        A typed observation ready for a like-for-like reference comparison.
    """
    required_text = (
        "reference_id", "material_id", "observable", "unit", "solver_family",
        "protocol", "artifact_sha256",
    )
    values = {key: _required_text(record, key) for key in required_text}
    if values["solver_family"] not in {"quantum_espresso", "cudaq_vqe"}:
        raise ValueError("computed observation has unsupported solver_family")
    digest = values["artifact_sha256"].lower()
    if len(digest) != 64 or any(character not in "0123456789abcdef" for character in digest):
        raise ValueError("artifact_sha256 must be a hexadecimal SHA-256 digest")
    value = float(record["value"])
    if not math.isfinite(value):
        raise ValueError("computed value must be finite")
    conditions = record.get("conditions")
    if not isinstance(conditions, dict) or not conditions:
        raise ValueError("computed observation requires physical conditions")
    calculation_conditions = record.get("calculation_conditions")
    if not isinstance(calculation_conditions, dict) or not calculation_conditions:
        raise ValueError("computed observation requires calculation conditions")
    return ComputedObservation(
        reference_id=values["reference_id"],
        material_id=values["material_id"],
        observable=values["observable"],
        value=value,
        unit=values["unit"],
        solver_family=values["solver_family"],  # type: ignore[arg-type]
        protocol=values["protocol"],
        conditions=conditions,
        calculation_conditions=calculation_conditions,
        converged=record.get("converged") is True,
        mock=record.get("mock") is True,
        candidate_specific=record.get("candidate_specific") is True,
        benchmarked=record.get("benchmarked") is True,
        artifact_sha256=digest,
    )


def compare_with_reference(
    reference: ReferenceObservation,
    observed: ComputedObservation,
) -> ReferenceComparison:
    """Compare a completed physical calculation to a like-for-like reference.

    Args:
        reference: Curated literature value and comparison protocol.
        observed: Completed, provenance-bound physical solver observation.

    Returns:
        Absolute-error comparison evaluated against the declared tolerance.

    Raises:
        ValueError: If the result is not scientifically comparable or lacks the
            evidence needed for the declared solver family.
    """
    identity_fields = ("reference_id", "material_id", "observable", "unit",
                       "solver_family", "protocol")
    mismatched = [field for field in identity_fields
                  if getattr(reference, field) != getattr(observed, field)]
    if mismatched:
        raise ValueError("reference comparison mismatch: " + ", ".join(mismatched))
    if dict(reference.conditions) != dict(observed.conditions):
        raise ValueError("reference comparison mismatch: conditions")
    if not observed.converged:
        raise ValueError("an unconverged calculation cannot validate a reference")
    if observed.mock:
        raise ValueError("mock calculations cannot validate a reference")
    if not observed.candidate_specific:
        raise ValueError("toy or generic models cannot validate a material reference")
    if reference.solver_family == "cudaq_vqe" and not observed.benchmarked:
        raise ValueError("VQE must first pass its exact-solver numerical benchmark")

    error = abs(observed.value - reference.value)
    return ReferenceComparison(
        reference_id=reference.reference_id,
        absolute_error=error,
        tolerance=reference.tolerance,
        passed=error <= reference.tolerance,
    )


def sha256_file(path: str | Path) -> str:
    """Return the SHA-256 identity used to bind an observation to its artifact.

    Args:
        path: Solver artifact whose immutable content identity is required.

    Returns:
        Lowercase hexadecimal SHA-256 digest.
    """
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()
