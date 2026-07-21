"""Auditable raw, canonical, and sampled-admissible design-space accounting."""

from __future__ import annotations

import argparse
import json
from math import comb
from pathlib import Path

from pipeline.common import catalyst_spaces as cs
from pipeline.common.design_space_provenance import validate_provenance
from pipeline.search.discovery import candidate_id
from pipeline.search.indexed_space import CLASS_ORDER, CLASS_SIZES, candidate_at_class, is_physically_admissible


def canonical_class_sizes() -> dict[str, int]:
    """Return exact counts after known representational equivalences collapse.

    These counts do not assert physical feasibility or symmetry uniqueness.
    They only remove equivalences explicitly implemented by ``candidate_id``.
    """
    sizes = dict(CLASS_SIZES)
    dopants = len(cs.SOLID_DOPANTS)
    sizes["SolidCatalyst"] = (
        len(cs.SOLID_ACTIVE_METALS) * len(cs.SOLID_SUPPORTS) *
        len(cs.SOLID_FACETS) * 20 * (dopants * (dopants + 1) // 2) * 4 * 3
    )

    positive_loadings = sum(float(x) > 0 for x in cs.MOLTEN_PROMOTER_AT_PCT)
    molten_per_temperature = 0
    for host in cs.MOLTEN_HOSTS:
        real_promoters = [p for p in cs.MOLTEN_PROMOTERS if p not in ("None", host)]
        molten_per_temperature += 1 + len(real_promoters) * positive_loadings
    sizes["MoltenMetal"] = molten_per_temperature * len(cs.MOLTEN_TEMPERATURES_K)

    positive_fractions = sum(float(x) > 0 for x in cs.PEROVSKITE_DOPANT_FRAC)
    sizes["Perovskite"] = (
        len(cs.PEROVSKITE_A_SITE) * len(cs.PEROVSKITE_B_SITE) *
        len(cs.PEROVSKITE_DEFECTS) *
        (1 + (len(cs.PEROVSKITE_B_SITE) + 1) * positive_fractions)
    )

    # The primary metal is present in HYDRIDE_SECOND_METAL; self and None
    # canonicalize to the same single-metal composition.
    sizes["MetalHydride"] = (
        len(cs.HYDRIDE_METALS) * len(cs.HYDRIDE_TYPES) *
        (len(cs.HYDRIDE_SECOND_METAL) - 1) * len(cs.HYDRIDE_ADDITIVES) * 9
    )
    sizes["MAXPhase"] = (
        len(cs.MAX_M_ELEMENTS) * len(cs.MAX_A_ELEMENTS) * len(cs.MAX_X_ELEMENTS) *
        len(cs.MAX_N_VALUES) * len(cs.MAX_M_ELEMENTS) * 3
    )

    mxene_site_choices = sum(
        len(cs.MXENE_SAC_METALS) - (1 if metal in cs.MXENE_SAC_METALS else 0)
        for metal in cs.MXENE_M_ELEMENTS
    )
    sizes["MXene"] = (
        mxene_site_choices * len(cs.MXENE_X_ELEMENTS) * len(cs.MXENE_N_VALUES) *
        len(cs.MXENE_TERMINATIONS)
    )
    return sizes


def _van_der_corput(value: int) -> float:
    result, denominator = 0.0, 1.0
    while value:
        value, digit = divmod(value, 2)
        denominator *= 2.0
        result += digit / denominator
    return result


def _sample_indices(size: int, count: int) -> list[int]:
    count = min(max(1, count), size)
    indices, seen, sequence = [], set(), 1
    for fixed in (0, size // 2, size - 1):
        if fixed not in seen and len(indices) < count:
            seen.add(fixed)
            indices.append(fixed)
    while len(indices) < count:
        index = min(size - 1, int(_van_der_corput(sequence) * size))
        sequence += 1
        if index not in seen:
            seen.add(index)
            indices.append(index)
    return indices


def audit_design_space(sample_per_class: int = 2048,
                       min_raw_per_class: int = 1000,
                       min_projected_admissible_per_class: int = 1000) -> dict:
    """Audit all classes without materializing the combinatorial population."""
    canonical_sizes = canonical_class_sizes()
    classes, failures = {}, []
    for material_class in CLASS_ORDER:
        raw_size = CLASS_SIZES[material_class]
        indices = _sample_indices(raw_size, sample_per_class)
        genomes = [candidate_at_class(material_class, i) for i in indices]
        admissible = [g for g in genomes if is_physically_admissible(g)[0]]
        unique_ids = {candidate_id(g) for g in genomes}
        fraction = len(admissible) / len(genomes)
        projected = int(round(raw_size * fraction))
        record = {
            "raw_cartesian_count": raw_size,
            "canonical_count": canonical_sizes[material_class],
            "canonical_fraction": canonical_sizes[material_class] / raw_size,
            "sampled_count": len(genomes),
            "sampled_unique_canonical_count": len(unique_ids),
            "sampled_admissible_count": len(admissible),
            "sampled_admissible_fraction": fraction,
            "projected_admissible_count": projected,
        }
        classes[material_class] = record
        if raw_size < min_raw_per_class:
            failures.append(f"raw_class_too_small:{material_class}:{raw_size}")
        if projected < min_projected_admissible_per_class:
            failures.append(f"admissible_class_too_small:{material_class}:{projected}")
        if not admissible:
            failures.append(f"no_sampled_admissible_candidates:{material_class}")

    provenance = validate_provenance(CLASS_ORDER)
    failures.extend(provenance["failures"])
    return {
        "valid": not failures,
        "failures": failures,
        "raw_cartesian_total": sum(CLASS_SIZES.values()),
        "canonical_total": sum(canonical_sizes.values()),
        "classes_represented": len(classes),
        "minimum_raw_per_class": min(x["raw_cartesian_count"] for x in classes.values()),
        "minimum_projected_admissible_per_class": min(
            x["projected_admissible_count"] for x in classes.values()),
        "classes": classes,
        "provenance": provenance,
        "interpretation": (
            "Canonical counts remove encoded equivalences only. Admissible counts "
            "are deterministic sample projections, not exhaustive feasibility proofs."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples-per-class", type=int, default=2048)
    parser.add_argument("--output", default="")
    args = parser.parse_args()
    report = audit_design_space(args.samples_per_class)
    payload = json.dumps(report, indent=2, sort_keys=True)
    if args.output:
        path = Path(args.output)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(payload + "\n")
    print(payload)
    if not report["valid"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
