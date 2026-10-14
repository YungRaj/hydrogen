"""Contracts for the uncalibrated NTEC/electrochemical multiphysics screening tier.

No MOOSE or Cantera process runs here: the solver boundary is the
deterministic ``FixtureSolver`` from ``tests/fixtures/reactor_cases.py``.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from pipeline.reactors.mechanisms import write_full_mechanism  # noqa: E402
from pipeline.simulation.external_runner import (  # noqa: E402
    SCREENING_RESULT, SCREENING_TIER, run_screening_case)
from pipeline.simulation.physical_case import load_physical_case  # noqa: E402
from pipeline.simulation.screening_cases import (  # noqa: E402
    build_screening_case, electrolyte_phase_for)
from tests.fixtures.reactor_cases import FixtureSolver  # noqa: E402

CASES = (
    ('NTEC', 'ntec', 473.15),
    ('Electrochemical', 'electrochemical', 353.15),
)


def _mechanism(tmp_path: Path) -> Path:
    return write_full_mechanism('cand', E_act_CH4=0.9, output_dir=tmp_path / 'mech')


@pytest.mark.parametrize('reactor_type,mode,temperature', CASES)
def test_screening_case_is_complete_but_uncalibrated(tmp_path, reactor_type, mode,
                                                     temperature):
    target = build_screening_case(
        tmp_path / 'case', candidate_id='cand', mode=mode,
        reactor_type=reactor_type, temperature_K=temperature,
        mechanism_path=_mechanism(tmp_path))
    case = json.loads(target.read_text())
    assert 'calibration' not in case
    assert case['kinetics']['interface_phase'] == 'cand_surface'
    assert (target.parent / case['kinetics']['mechanism']).is_file()
    assert case['operating']['temperature_K'] == temperature
    assert case['screening']['calibrated'] is False
    loaded = load_physical_case(
        target, candidate_id='cand', mode=mode, reactor_type=reactor_type,
        temperature_K=temperature, require_calibration=False)
    assert loaded['candidate_id'] == 'cand'
    with pytest.raises(ValueError, match='calibration'):
        load_physical_case(
            target, candidate_id='cand', mode=mode, reactor_type=reactor_type,
            temperature_K=temperature)


def test_uncalibrated_loading_still_enforces_physics(tmp_path):
    target = build_screening_case(
        tmp_path / 'case', candidate_id='cand', mode='electrochemical',
        reactor_type='Electrochemical', temperature_K=353.15,
        mechanism_path=_mechanism(tmp_path))
    case = json.loads(target.read_text())
    del case['parameter_sources']['properties.ionic_conductivity_S_m']
    case['operating']['temperature_K'] = 700.0  # aqueous above water's T_c
    target.write_text(json.dumps(case))
    with pytest.raises(ValueError) as error:
        load_physical_case(
            target, candidate_id='cand', mode='electrochemical',
            reactor_type='Electrochemical', temperature_K=700.0,
            require_calibration=False)
    assert 'parameter_sources.properties.ionic_conductivity_S_m' in str(error.value)
    assert 'aqueous_temperature_below_water_critical_point' in str(error.value)


def test_electrolyte_phase_follows_temperature_unless_declared(monkeypatch):
    monkeypatch.delenv('ELECTROCHEMICAL_CONDITIONS_JSON', raising=False)
    assert electrolyte_phase_for(353.15) == 'aqueous'
    assert electrolyte_phase_for(823.15) == 'molten'
    monkeypatch.setenv('ELECTROCHEMICAL_CONDITIONS_JSON',
                       json.dumps({'electrolyte_phase': 'molten'}))
    assert electrolyte_phase_for(353.15) == 'molten'


def test_build_refuses_non_empty_directories(tmp_path):
    case = tmp_path / 'case'
    case.mkdir()
    (case / 'stale.txt').write_text('x')
    with pytest.raises(ValueError, match='not empty'):
        build_screening_case(case, candidate_id='cand', mode='ntec',
                             reactor_type='NTEC', temperature_K=473.15,
                             mechanism_path=_mechanism(tmp_path))


@pytest.mark.parametrize('reactor_type,mode,temperature', CASES)
def test_run_screening_case_is_labelled_and_never_persists_artifacts(
        tmp_path, reactor_type, mode, temperature):
    target = build_screening_case(
        tmp_path / 'case', candidate_id='cand', mode=mode,
        reactor_type=reactor_type, temperature_K=temperature,
        mechanism_path=_mechanism(tmp_path))
    result = run_screening_case(
        mode=mode, reactor_type=reactor_type, candidate_id='cand',
        temperature_K=temperature, case_dir=target.parent,
        execution=FixtureSolver(reactor_type, 'cand', temperature).services())
    assert result['valid'] is True
    assert result['tier'] == SCREENING_TIER == 'uncalibrated_multiphysics_screening'
    assert result['calibrated'] is False
    assert result['candidate_exclusion_authorized'] is False
    assert 'model_validation' not in result and 'calibration' not in result
    assert result['backend_solvers'] == {
        'moose': 'test-fixture-moose', 'cantera': 'test-fixture-cantera'}
    assert result['solver_coupling']['converged'] is True
    assert len(result['input_sha256']) == 64
    assert 0 <= result['outputs']['CH4_conversion'] <= 1
    assert json.loads((target.parent / SCREENING_RESULT).read_text()) == result
    # No production artifact (<candidate>/<mode>/<reactor>_<T>K.json) anywhere.
    assert not list(tmp_path.rglob(f'{reactor_type}_*K.json'))


def test_run_screening_case_rejects_thermal_reactors(tmp_path):
    with pytest.raises(ValueError, match='no uncalibrated screening tier'):
        run_screening_case(mode='mmbcr', reactor_type='MMBCR', candidate_id='cand',
                           temperature_K=1200.0, case_dir=tmp_path)


class _UnphysicalSolver(FixtureSolver):
    """Fixture whose final outputs violate the NTEC physical checks."""

    def execute(self, command, case, timeout, backend='solver'):
        super().execute(command, case, timeout, backend)
        outputs = json.loads((case / 'hydrogen_outputs.json').read_text())
        outputs['specific_energy_kWh_kg_H2'] = -1.0
        (case / 'hydrogen_outputs.json').write_text(json.dumps(outputs))


def test_unphysical_solver_output_raises_runtime_error(tmp_path):
    target = build_screening_case(
        tmp_path / 'case', candidate_id='cand', mode='ntec', reactor_type='NTEC',
        temperature_K=473.15, mechanism_path=_mechanism(tmp_path))
    with pytest.raises(RuntimeError, match='physical checks'):
        run_screening_case(
            mode='ntec', reactor_type='NTEC', candidate_id='cand',
            temperature_K=473.15, case_dir=target.parent,
            execution=_UnphysicalSolver('NTEC', 'cand', 473.15).services())
    assert not (target.parent / SCREENING_RESULT).exists()


def _screening_record(reactor_type: str) -> dict:
    outputs = {
        'NTEC': {'CH4_conversion': 0.03, 'H2_selectivity': 1.0,
                 'solid_C_selectivity': 1.0, 'specific_energy_kWh_kg_H2': 40.0},
        'Electrochemical': {'CH4_conversion': 0.02, 'H2_selectivity': 0.5,
                            'faradaic_efficiency_H2': 1.0,
                            'current_density_A_cm2': 0.001, 'cell_voltage_V': 0.9,
                            'electrical_power_density_W_cm2': 0.0009},
    }[reactor_type]
    record = {'valid': True, 'tier': SCREENING_TIER,
              'candidate_exclusion_authorized': False, 'outputs': outputs}
    if reactor_type == 'Electrochemical':
        record['electrolyte_phase'] = 'molten'
    return record


@pytest.mark.parametrize('reactor_type,mode,temperature', CASES)
def test_adapters_report_uncalibrated_tier_and_calibrated_wins(
        reactor_type, mode, temperature):
    from pipeline.reactors.electrochemical_models import (
        simulate_electrochemical_pathway, simulate_ntec_pathway)
    from pipeline.reactors.reactor_core import ReactorConfig

    simulate = simulate_ntec_pathway if reactor_type == 'NTEC' \
        else simulate_electrochemical_pathway
    screening = _screening_record(reactor_type)
    config = ReactorConfig(reactor_type=reactor_type, pathway_mode=mode,
                           candidate_id='cand', T_inlet_K=temperature,
                           screening_multiphysics=screening)
    result = simulate(config)
    assert result['status'] == 'complete' and result['valid'] is True
    assert result['reactor_evidence_tier'] == 'uncalibrated_multiphysics_screening'
    assert result['can_exclude_candidate'] is False
    assert result['multiphysics_screening'] is screening
    assert result['limitations'] == ['uncalibrated_representative_operating_point']
    for name, value in screening['outputs'].items():
        assert result[name] == value
    if reactor_type == 'Electrochemical':
        assert result['electrolyte_phase'] == 'molten'

    artifact_outputs = dict(screening['outputs'], CH4_conversion=0.25)
    config.multiphysics_artifact = {'artifact': {
        'outputs': artifact_outputs, 'electrolyte_phase': 'aqueous'}}
    calibrated = simulate(config)
    assert calibrated['CH4_conversion'] == 0.25
    assert calibrated['reactor_evidence_tier'] in {
        'calibrated_multiphysics_screening', 'mechanistic_multiphysics_screening'}
    assert 'multiphysics_screening' not in calibrated


def test_reactor_stage_runs_and_reuses_screening_cases(tmp_path):
    from pipeline.reactors.models import ReactorConfig, simulate_reactor

    mechanism = _mechanism(tmp_path)
    root = tmp_path / 'multiphysics'
    services = FixtureSolver('NTEC', 'cand', 473.15).services()

    def config():
        return ReactorConfig(reactor_type='NTEC', pathway_mode='ntec',
                             catalyst_name='cand', candidate_id='cand',
                             material_class='MoltenMetal', T_inlet_K=473.15,
                             mechanism_file=str(mechanism),
                             multiphysics_results_dir=str(root))

    with (
        patch('pipeline.simulation.screening_cases.screening_available',
              return_value=True),
        patch('pipeline.simulation.external_runner.default_solver_execution_services',
              return_value=services) as factory,
    ):
        first = simulate_reactor(config())
        second = simulate_reactor(config())
    assert first['status'] == 'complete'
    assert first['reactor_evidence_tier'] == 'uncalibrated_multiphysics_screening'
    assert first['can_exclude_candidate'] is False
    assert first['multiphysics_screening']['reused'] is False
    assert second['multiphysics_screening']['reused'] is True
    assert factory.call_count == 1
    case_dir = root / 'screening_cases' / 'cand' / 'NTEC' / '473.15K'
    assert (case_dir / SCREENING_RESULT).is_file()
    assert not (root / 'cand').exists()  # no production artifact tree


def test_reactor_stage_records_screening_failure(tmp_path):
    from pipeline.reactors.models import ReactorConfig, simulate_reactor

    with patch('pipeline.simulation.screening_cases.screening_available',
               return_value=True):
        result = simulate_reactor(ReactorConfig(
            reactor_type='Electrochemical', pathway_mode='electrochemical',
            catalyst_name='cand', candidate_id='cand', T_inlet_K=353.15,
            mechanism_file=str(tmp_path / 'missing.yaml'),
            multiphysics_results_dir=str(tmp_path / 'multiphysics')))
    assert result['status'] == 'validation_required'
    assert result['reactor_evidence_tier'] == 'required_solver_evidence_missing'
    assert 'mechanism file is missing' in result['multiphysics_evidence']['screening_error']
    assert not (tmp_path / 'multiphysics').exists()
