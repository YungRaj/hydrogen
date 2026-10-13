from pipeline.search.design_space import (
    extract_elements,
    generate_supported_nanoparticle_genome,
)
from pipeline.screening.surface_screener import generate_structure
from pipeline.reactors.mechanisms import off_site_carbon_allowed
from pipeline.reactors.modes import reactor_applicability


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
