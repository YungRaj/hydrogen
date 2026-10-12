import pytest

from pipeline.electrochemistry.sensitivity import sensitivity_band


def test_sensitivity_band_is_deterministic_and_bounded() -> None:
    result = sensitivity_band(lambda values: values['x'] ** 2, {
        'low': {'x': 1.0}, 'nominal': {'x': 2.0}, 'high': {'x': 3.0}
    })
    assert result['nominal'] == pytest.approx(4.0)
    assert result['minimum'] == pytest.approx(1.0)
    assert result['maximum'] == pytest.approx(9.0)
