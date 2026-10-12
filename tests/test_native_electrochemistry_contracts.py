from math import isclose

import pytest

from pipeline.electrochemistry.cell import CellOperatingPoint, cell_voltage
from pipeline.electrochemistry.kinetics import butler_volmer_current_density
from pipeline.electrochemistry.ntec_mechanics import (
    collision_events_per_volume_s,
    detachment_probability,
    impact_energy_J,
)
from pipeline.electrochemistry.transport import limiting_current_A


def test_butler_volmer_is_zero_at_zero_overpotential() -> None:
    assert butler_volmer_current_density(1.0, 0.0, 1000.0) == pytest.approx(0.0)


def test_transport_and_cell_balance_are_dimensionally_bounded() -> None:
    limit = limiting_current_A(2, 1e-9, 1000.0, 1e-4, 1e-4)
    point = CellOperatingPoint(1.0, 0.2, -0.1, limit * 0.1, limit, 2.0, 1000.0)
    assert cell_voltage(point) > 1.2


def test_ntec_mechanics_monotonicity() -> None:
    low = impact_energy_J(1e-9, 100.0, 1e-4)
    high = impact_energy_J(1e-9, 200.0, 1e-4)
    assert high > low
    assert detachment_probability(1e-19, high) > detachment_probability(1e-19, low)
    assert collision_events_per_volume_s(100.0, 1e-4, 1e12) > 0.0
