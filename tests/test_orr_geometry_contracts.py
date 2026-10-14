from ase import Atoms

from pipeline.screening.fc_screener import _adsorbate_geometry_valid


def test_orr_geometry_rejects_broken_ooh():
    slab = Atoms('NiOOH', positions=[(0, 0, 0), (0, 0, 2), (4, 0, 2), (4, 0, 3)])
    valid, reason = _adsorbate_geometry_valid(slab, 1, 'ooh')
    assert not valid
    assert reason == 'ooh_bond_broken'


def test_orr_geometry_accepts_normal_ooh():
    slab = Atoms('NiOOH', positions=[(0, 0, 0), (0, 0, 2), (1.2, 0, 2), (1.2, 0, 2.97)])
    valid, reason = _adsorbate_geometry_valid(slab, 1, 'ooh')
    assert valid
    assert reason == ''
