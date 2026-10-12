from pipeline.electrochemistry.native_models import (
    electrochemical_screening_point,
    ntec_screening_point,
)


def test_native_electrochemical_point_is_non_pruning() -> None:
    result = electrochemical_screening_point({
        'descriptor_eV': 0.0, 'optimum_descriptor_eV': 0.0,
        'reference_current_A_cm2': 1.0, 'temperature_K': 1000.0,
        'overpotential_V': 0.1, 'reversible_voltage_V': 1.0,
        'area_m2': 1e-4, 'diffusivity_m2_s': 1e-9,
        'concentration_mol_m3': 1000.0, 'boundary_layer_m': 1e-4,
        'gap_m': 1e-3, 'conductivity_S_m': 10.0,
    })
    assert result['reactor_evidence_tier'] == 'native_screening_model'
    assert result['can_exclude_candidate'] is False


def test_native_ntec_point_is_non_pruning() -> None:
    result = ntec_screening_point({
        'shear_rate_s': 100.0, 'particle_diameter_m': 1e-4,
        'number_density_m3': 1e12, 'effective_mass_kg': 1e-9,
        'adhesion_energy_J': 1e-19,
    })
    assert result['detachment_probability'] > 0.0
    assert result['can_exclude_candidate'] is False
