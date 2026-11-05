"""Portable contracts for optional atomistic validation helpers."""

import pytest

pytest.importorskip('ase')

from ase import Atoms

from pipeline.evidence.molten_md import LiquidMDProtocol
from pipeline.evidence.saddle_refinement import refine_saddle
from pipeline.evidence import molten_neb
from pipeline.evidence.thermochemistry_crosscheck import compare_gas_reference


def test_md_protocol_rejects_invalid_temperature() -> None:
    with pytest.raises(ValueError, match='temperature_K'):
        LiquidMDProtocol(0.0).validate()


def test_saddle_refinement_requires_calculator() -> None:
    with pytest.raises(ValueError, match='attached calculator'):
        refine_saddle(Atoms('H'))


def test_path_saddle_selects_highest_internal_image(monkeypatch) -> None:
    images = [Atoms('H', positions=[(0, 0, float(i))]) for i in range(5)]
    for i, image in enumerate(images):
        image.get_potential_energy = lambda i=i: float(i)
    monkeypatch.setattr(molten_neb, 'refine_saddle', lambda atoms, **kwargs: {
        'converged': True, 'steps': 1, 'energy_eV': 0.0, 'max_force_eV_A': 0.0,
    })
    result = molten_neb.refine_path_saddle(images)
    assert result['image_index'] == 3


def test_bond_scan_rejects_invalid_indices() -> None:
    with pytest.raises(IndexError):
        molten_neb.run_bond_stretch_scan(
            Atoms('CH', positions=((0, 0, 0), (1.1, 0, 0))),
            0, 2, object(), 0.0,
        )


def test_bond_scan_tracks_requested_distances() -> None:
    from ase.calculators.emt import EMT
    state = Atoms('CH', positions=((0, 0, 0), (1.1, 0, 0)))
    result = molten_neb.run_bond_stretch_scan(
        state, 0, 1, EMT(), 0.0, start_A=1.1, stop_A=1.3,
        step_A=0.1, steps=5,
    )
    actual = [row['actual_CH_A'] for row in result['records']]
    assert actual == pytest.approx([1.1, 1.2, 1.3], abs=0.02)


def test_thermochemistry_comparison_is_non_mutating() -> None:
    result = compare_gas_reference(
        current_free_energy_eV=-1.0,
        pmutt_free_energy_eV=-0.98,
        tolerance_eV=0.03,
    )
    assert result['within_tolerance'] is True
    assert result['difference_eV'] == pytest.approx(0.02)


def test_molten_methane_endpoints_preserve_slab_atom_matching(monkeypatch) -> None:
    from ase.build import fcc111
    from ase.calculators.emt import EMT

    slab = fcc111('Cu', size=(3, 3, 3), vacuum=10.0)
    slab.pbc = True
    monkeypatch.setattr(
        molten_neb,
        'relax_with_record',
        lambda atoms, label, fmax, steps: {f'relax_{label}_converged': True, f'relax_{label}_steps': 1},
    )
    cfg = molten_neb.NEBConfig(minimum_product_CH_A=1.2, maximum_product_CH_A=3.5)
    initial, final, meta = molten_neb.methane_endpoints(slab, 0, EMT(), cfg)
    assert meta['endpoint_match']['passed'] is True
    assert meta['endpoint_match']['max_surface_displacement_A'] == pytest.approx(0.0)
    assert len(initial) == len(final) == len(slab) + 5


def test_run_liquid_alloy_md_collects_decorrelated_snapshots() -> None:
    from ase.calculators.emt import EMT
    from pipeline.evidence.molten_md import LiquidMDProtocol, run_liquid_alloy_md

    proto = LiquidMDProtocol(
        temperature_K=1273.15,
        timestep_fs=1.0,
        equilibration_steps=4,
        stride_steps=2,
        snapshots=3,
        friction_per_fs=0.01,
        seed=7,
    )
    res = run_liquid_alloy_md(('MoltenMetal', 'Cu', 'Ni', 22.0, 1273.15), EMT(), proto)
    assert res['complete'] is True
    assert res['n_snapshots'] == 3
    assert isinstance(res['mean_potential_energy_eV'], float)


def test_crystallographic_oxide_and_max_slabs_expose_active_surface_cations() -> None:
    from pipeline.screening.structure_builders import generate_structure

    perov_slab, perov_idx, _ = generate_structure(
        ('Perovskite', 'La', 'Ni', 'Fe', 0.25, 'O_vacancy')
    )
    assert perov_slab[perov_idx[0]].symbol in {'Fe', 'Ni'}
    assert float(perov_slab.positions[perov_idx[0], 2]) >= float(perov_slab.positions[:, 2].max()) - 0.85

    spinel_slab, spinel_idx, _ = generate_structure(
        ('Spinel', 'Mg', 'Co', 'Ce', 'nanoparticle', 'N-graphene')
    )
    assert spinel_slab[spinel_idx[0]].symbol in {'Ce', 'Co'}
    assert float(spinel_slab.positions[spinel_idx[0], 2]) >= float(spinel_slab.positions[:, 2].max()) - 0.85

    max_slab, max_idx, _ = generate_structure(
        ('MAXPhase', 'Co', 'Al', 'C', 2, 'Ni', 'basal_0001')
    )
    assert max_slab[max_idx[0]].symbol in {'Ni', 'Co'}
    assert 'C' in max_slab.get_chemical_symbols()
    assert 'Al' in max_slab.get_chemical_symbols()


def test_catalyst_embedded_casscf_8qubit_hamiltonian_and_diagnostics(tmp_path):
    """Verify 8-qubit CASSCF(4,4) catalyst-embedded TS Hamiltonian generation and diagnostics."""
    import pytest

    pytest.importorskip('pyscf')
    pytest.importorskip('openfermion')
    from pipeline.search.discovery import candidate_id
    from pipeline.validation.candidate_hamiltonian import (
        build_candidate_hamiltonian,
        generate_catalyst_ts_fcidump,
    )

    genome = ('SAA', 'Pt', 'Ni', '110', 0.05)
    cid = candidate_id(genome)
    out_dir = tmp_path / cid
    meta = generate_catalyst_ts_fcidump(
        cid,
        genome,
        reaction_type='CH_split',
        output_dir=out_dir,
        active_electrons=4,
        active_orbitals=4,
    )
    assert meta['active_electrons'] == 4
    assert meta['active_orbitals'] == 4
    assert 'Pt' in str(meta['catalyst_motif'])
    assert float(meta['correlation_energy_Ha']) < 0.0
    assert 0.0 < float(meta['leading_ci_weight']) <= 1.0

    ham = build_candidate_hamiltonian(out_dir / 'FCIDUMP')
    assert ham['n_qubits'] == 8
    assert ham['catalyst_specific_hamiltonian'] is True
    assert ham['catalyst_motif'] == meta['catalyst_motif']
