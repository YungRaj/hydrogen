import pytest

from pipeline.electrochemistry.attribution import compare_regeneration


def test_regeneration_attribution_reports_deltas() -> None:
    baseline = {'CH4_conversion': 0.2, 'H2_selectivity': 0.8,
                'solid_C_selectivity': 0.7, 'specific_energy_kWh_kg_H2': 4.0}
    enabled = {'CH4_conversion': 0.3, 'H2_selectivity': 0.85,
               'solid_C_selectivity': 0.75, 'specific_energy_kWh_kg_H2': 4.2}
    result = compare_regeneration(baseline, enabled)
    assert result['delta']['CH4_conversion'] == pytest.approx(0.1)
    assert result['mechanism'] == 'carbon_regeneration'
