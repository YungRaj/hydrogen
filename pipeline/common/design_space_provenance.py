"""Machine-readable provenance for the encoded catalyst design grammar.

The references support the material families and descriptor choices. They do
not imply that every Cartesian product member has been reported or synthesized.
"""

from __future__ import annotations


SOURCES = {
    "project_turquoise_review": {
        "kind": "curated_project_review",
        "citation": "TURQUOISE_HYDROGEN.md",
        "scope": "methane-pyrolysis catalyst families, reactors, and operating ranges",
    },
    "project_fuel_cell_review": {
        "kind": "curated_project_review",
        "citation": "FUEL_CELL.md",
        "scope": "ORR catalyst families, coordination motifs, MEAs, and stacks",
    },
    "materials_project": {
        "kind": "external_database_method",
        "citation": "Jain et al., APL Materials 1, 011002 (2013), doi:10.1063/1.4812323",
        "scope": "inorganic compositions and crystal prototypes",
    },
    "catalysis_hub": {
        "kind": "external_database_method",
        "citation": "Winther et al., Scientific Data 6, 75 (2019), doi:10.1038/s41597-019-0081-y",
        "scope": "surface structures and adsorption-energy data organization",
    },
    "open_catalyst_2020": {
        "kind": "external_dataset_method",
        "citation": "Chanussot et al., ACS Catalysis 11, 6059-6072 (2021), doi:10.1021/acscatal.0c04525",
        "scope": "adsorbate-surface chemical diversity and ML screening",
    },
}


def _axes(names, *sources):
    return {
        name: {
            "selection_basis": "project_curated_extrapolation",
            "source_ids": list(sources),
        }
        for name in names
    }


DESIGN_AXIS_PROVENANCE = {
    "MoltenMetal": _axes(
        ("host", "promoter", "promoter_at_pct", "temperature_K"),
        "project_turquoise_review", "materials_project"),
    "SolidCatalyst": _axes(
        ("active_metal", "support", "facet", "strain", "dopant_1",
         "dopant_2", "substitution_count", "vacancy_count"),
        "project_turquoise_review", "catalysis_hub", "open_catalyst_2020"),
    "SAC": _axes(
        ("metal", "coordination", "substrate", "axial_ligand"),
        "project_fuel_cell_review", "catalysis_hub"),
    "DAC": _axes(
        ("metal_1", "metal_2", "coordination", "substrate"),
        "project_fuel_cell_review", "catalysis_hub"),
    "MOF": _axes(
        ("metal_node", "linker", "cavity", "pore_size_A"),
        "materials_project", "project_turquoise_review", "project_fuel_cell_review"),
    "COF": _axes(
        ("metal", "linkage", "cavity", "pore_size_A"),
        "project_turquoise_review", "project_fuel_cell_review"),
    "Perovskite": _axes(
        ("a_site", "b_site", "dopant", "dopant_fraction", "defect"),
        "materials_project", "project_fuel_cell_review"),
    "MetalHydride": _axes(
        ("primary_metal", "hydride_type", "secondary_metal", "additive", "temperature_K"),
        "materials_project", "project_turquoise_review"),
    "MAXPhase": _axes(
        ("m_element", "a_element", "x_element", "n", "dopant", "facet"),
        "materials_project", "project_turquoise_review", "project_fuel_cell_review"),
    "HEA": _axes(
        ("components", "structure", "facet", "temperature_K"),
        "materials_project", "catalysis_hub", "project_turquoise_review"),
    "Spinel": _axes(
        ("a_metal", "b_metal", "dopant", "morphology", "carbon_support"),
        "materials_project", "project_fuel_cell_review"),
    "MXene": _axes(
        ("m_element", "x_element", "n", "termination", "single_atom_metal"),
        "materials_project", "project_fuel_cell_review"),
    "SAA": _axes(
        ("trace_metal", "host_metal", "facet", "loading_ppm"),
        "catalysis_hub", "project_turquoise_review", "project_fuel_cell_review"),
    "MetalFreeCarbon": _axes(
        ("nitrogen_type", "nitrogen_fraction", "defect", "substrate", "co_dopant"),
        "project_turquoise_review", "project_fuel_cell_review"),
}


def validate_provenance(material_classes) -> dict:
    """Fail closed if a class, axis, or referenced source lacks provenance."""
    failures = []
    for material_class in material_classes:
        axes = DESIGN_AXIS_PROVENANCE.get(material_class)
        if not axes:
            failures.append(f"missing_class:{material_class}")
            continue
        for axis, record in axes.items():
            source_ids = record.get("source_ids", ())
            if not source_ids:
                failures.append(f"missing_sources:{material_class}:{axis}")
            for source_id in source_ids:
                if source_id not in SOURCES:
                    failures.append(f"unknown_source:{material_class}:{axis}:{source_id}")
    unknown_classes = sorted(set(DESIGN_AXIS_PROVENANCE) - set(material_classes))
    failures.extend(f"unknown_class:{name}" for name in unknown_classes)
    return {
        "valid": not failures,
        "failures": failures,
        "classes": len(DESIGN_AXIS_PROVENANCE),
        "axes": sum(len(axes) for axes in DESIGN_AXIS_PROVENANCE.values()),
        "sources": len(SOURCES),
        "interpretation": (
            "Family/axis provenance only; Cartesian candidates remain hypotheses "
            "until prior-art, synthesis, and performance validation."
        ),
    }
