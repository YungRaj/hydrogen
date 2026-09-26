#!/usr/bin/env python3
"""Opt-in grading of real QE/VQE observables against curated references.

Set HYDROGEN_QUANTUM_REFERENCE_RESULTS to a JSON file with an ``observations``
array.  This suite intentionally consumes completed, checksum-bound artifacts;
it does not hide an expensive production calculation inside a unit test.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.validation.reference_benchmarks import (
    compare_with_reference,
    computed_observation,
    load_reference_manifest,
    sha256_file,
)


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = ROOT / "tests" / "fixtures" / "quantum_reference_manifest.json"


def main() -> None:
    result_path = os.environ.get("HYDROGEN_QUANTUM_REFERENCE_RESULTS")
    if not result_path:
        raise SystemExit(
            "HYDROGEN_QUANTUM_REFERENCE_RESULTS must point to real solver observations")
    manifest_path = Path(result_path).expanduser().resolve()
    payload = json.loads(manifest_path.read_text())
    observations = payload.get("observations")
    if not isinstance(observations, list) or not observations:
        raise AssertionError("result manifest must contain observations")
    references = load_reference_manifest(REFERENCES)
    for item in observations:
        if not isinstance(item, dict):
            raise AssertionError("each solver observation must be an object")
        raw = dict(item)
        artifact_value = raw.pop("artifact_path", None)
        if not isinstance(artifact_value, str) or not artifact_value:
            raise AssertionError("each observation requires artifact_path")
        artifact = Path(artifact_value).expanduser()
        if not artifact.is_absolute():
            artifact = manifest_path.parent / artifact
        artifact = artifact.resolve()
        if not artifact.is_file():
            raise AssertionError(f"solver artifact does not exist: {artifact}")
        if raw.get("artifact_sha256") != sha256_file(artifact):
            raise AssertionError(f"solver artifact checksum mismatch: {artifact}")
        observed = computed_observation(raw)
        if observed.reference_id not in references:
            raise AssertionError(f"unknown reference_id: {observed.reference_id}")
        comparison = compare_with_reference(references[observed.reference_id], observed)
        assert comparison.passed, comparison
        print("PASS", comparison)


if __name__ == "__main__":
    main()
