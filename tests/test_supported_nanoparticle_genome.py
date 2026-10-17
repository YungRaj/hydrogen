from pipeline.search.design_space import (
    extract_elements,
    generate_supported_nanoparticle_genome,
)
from pipeline.screening.surface_screener import generate_structure
from pipeline.reactors.mechanisms import off_site_carbon_allowed
from pipeline.reactors.modes import reactor_applicability
from pipeline.screening.supported_nanoparticle import SUPPORTED_NANOPARTICLE_PROTOCOL


def test_supported_nanoparticle_genome_exposes_size_and_loading():
    genome = generate_supported_nanoparticle_genome()
    assert genome[0] == 'SupportedNanoparticle'
    assert 1.0 <= genome[4] <= 20.0
    assert 1.0 <= genome[5] <= 30.0
    assert genome[1] in extract_elements(genome)


def test_supported_nanoparticle_uses_validated_surface_builder():
    genome = ('SupportedNanoparticle', 'Ni', 'SiO2', 'fcc111', 5.0, 10.0, (), 1, 0)
    slab, active, material_class = generate_structure(genome)
    assert material_class == 'SupportedNanoparticle'
    assert active
    assert slab.info['particle_diameter_nm'] == 5.0
    assert slab.info['loading_wt_pct'] == 10.0


def test_supported_nanoparticle_enters_b6_gate_for_ni_fe_co():
    assert off_site_carbon_allowed(
        ('SupportedNanoparticle', 'Ni', 'SiO2', 'fcc111', 5.0, 10.0, (), 1, 0)
    )


def test_supported_nanoparticle_is_applicable_to_pfr():
    applicable, reason = reactor_applicability('PFR', 'SupportedNanoparticle')
    assert applicable
    assert reason is None


def test_supported_nanoparticle_protocol_is_versioned():
    assert SUPPORTED_NANOPARTICLE_PROTOCOL == 'supported-nanoparticle-v1'


def test_supported_nanoparticle_varies_slab_morphology_with_size_and_loading():
    small_np = ('SupportedNanoparticle', 'Ni', 'SiO2', 'fcc111', 2.5, 5.0, (), 1, 0)
    large_np = ('SupportedNanoparticle', 'Ni', 'SiO2', 'fcc111', 15.0, 25.0, (), 1, 0)
    small_slab, _, _ = generate_structure(small_np)
    large_slab, _, _ = generate_structure(large_np)
    assert small_slab.info['dispersion'] > large_slab.info['dispersion']
    assert small_slab.info['step_site_fraction'] > large_slab.info['step_site_fraction']
    assert small_slab.info['effective_strain'] != large_slab.info['effective_strain']
    assert small_slab.info['effective_surface_vacancies'] == 1
    assert large_slab.info['effective_surface_vacancies'] == 0
    assert len(small_slab) == len(large_slab) - 1
    assert not (small_slab.cell.lengths() == large_slab.cell.lengths()).all()

