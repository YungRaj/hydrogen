from types import SimpleNamespace

from pipeline.reactors.electrochemical_models import (
    simulate_electrochemical_pathway,
    simulate_ntec_pathway,
)


def _config(**kwargs):
    base = {
        'multiphysics_artifact': None,
        'catalyst_name': 'test', 'candidate_id': 'test', 'T_inlet_K': 1000.0,
        'native_electrochemical_parameters': None,
        'native_ntec_parameters': None,
    }
    base.update(kwargs)
    return SimpleNamespace(**base)


def test_native_electrochemical_status_follows_model_result():
    params = {
        'descriptor_eV': 0.0, 'optimum_descriptor_eV': 0.0,
        'reference_current_A_cm2': 1.0, 'temperature_K': 1000.0,
        'overpotential_V': 0.001, 'reversible_voltage_V': 1.0,
        'area_m2': 1e-4, 'diffusivity_m2_s': 1e-9,
        'concentration_mol_m3': 1000.0, 'boundary_layer_m': 1e-4,
        'gap_m': 1e-3, 'conductivity_S_m': 10.0,
        'ch4_feed_mol_s': 1e6,
    }
    result = simulate_electrochemical_pathway(
        _config(native_electrochemical_parameters=params)
    )
    assert result['status'] == 'complete'
    assert result['valid'] is True


def test_native_ntec_remains_non_validating_diagnostic():
    params = {
        'shear_rate_s': 100.0, 'particle_diameter_m': 1e-4,
        'number_density_m3': 1e12, 'effective_mass_kg': 1e-9,
        'adhesion_energy_J_m2': 1.0, 'contact_area_m2': 1e-10,
        'active_site_density_m2': 1e18,
    }
    result = simulate_ntec_pathway(_config(native_ntec_parameters=params))
    assert result['status'] == 'validation_required'
    assert result['valid'] is False
