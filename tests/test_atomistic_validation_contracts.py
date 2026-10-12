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


def test_thermochemistry_comparison_is_non_mutating() -> None:
    result = compare_gas_reference(
        current_free_energy_eV=-1.0,
        pmutt_free_energy_eV=-0.98,
        tolerance_eV=0.03,
    )
    assert result['within_tolerance'] is True
    assert result['difference_eV'] == pytest.approx(0.02)
