"""Reference-campaign configuration and isolated sweep artifacts."""

from dataclasses import replace
import json
from pathlib import Path

import pandas as pd
import pytest

from pipeline.orchestrator import PipelineConfig, normalized_pipeline_config
from pipeline.process import reactor_mechanisms, reactor_models, yaml_sweep
from pipeline.process.phase2_scorecard import build_solids_scorecard
from pipeline.stages.reactor_batch import (
    ReactorBatchServices, default_reactor_batch_services, run_reactor_batch_stage,
)


@pytest.mark.parametrize('quick_mode', [False, True])
def test_normalization_preserves_requested_reference_temperatures(quick_mode):
    config = PipelineConfig(
        reactor_temperatures=(923.15, 973.15), quick_mode=quick_mode)
    effective = normalized_pipeline_config(config)
    assert effective is not config
    assert effective.reactor_temperatures == (923.15, 973.15)
    assert config.reactor_temperatures == (923.15, 973.15)


def test_existing_reference_keeps_its_name_and_does_not_mutate_candidates():
    candidates = pd.DataFrame([
        {'candidate_id': 'ordinary', 'E_act': 0.8},
        {'candidate_id': 'ni_np_lit', 'E_act': 1.0},
    ])
    original = candidates.copy(deep=True)
    names = []

    def simulate(row, name, temperatures, reactors, **kwargs):
        names.append(name)
        return {'sweep': [dict(
            catalyst_name=name, reactor_type='PFR', T_K=temperature,
            CH4_conversion=0.15, single_pass_CH4_conversion=0.15,
            surface_loaded=True, status='complete')
            for temperature in temperatures]}

    def unexpected_reference_load(name):
        raise AssertionError('an existing reference must not be loaded again')

    services = ReactorBatchServices(
        prepare_gas_mechanism=lambda: None, simulate_candidate=simulate,
        write_mock_mechanism=lambda *args, **kwargs: None,
        run_mock_sweep=lambda *args, **kwargs: [],
        build_scorecard=build_solids_scorecard,
        load_reference_candidate=unexpected_reference_load)
    result = run_reactor_batch_stage(
        candidates, temperatures=(923.15, 973.15), reactor_types=('PFR',),
        pathway_mode='thermocatalytic_pfr', multiphysics_results_dir='',
        allow_mock_inputs=False, judge_catalyst='ni_np_lit',
        headline_t_min=923.15, headline_t_max=973.15, services=services)
    assert names == ['cat_0', 'ni_np_lit']
    assert result.state['solids_scorecard']['judge_catalyst'] == 'ni_np_lit'
    assert result.state['best_conversion'] == 0.15
    pd.testing.assert_frame_equal(candidates, original)


def test_default_campaign_loads_ni_reference_and_produces_a_real_headline(
        tmp_path, monkeypatch):
    pytest.importorskip('cantera')
    monkeypatch.setattr(reactor_mechanisms, 'MECHANISMS_DIR', tmp_path)
    monkeypatch.setattr(reactor_models, 'save_json', lambda *args, **kwargs: None)
    config = normalized_pipeline_config(PipelineConfig())
    candidates = pd.DataFrame(columns=['candidate_id', 'E_act'])
    services = replace(
        default_reactor_batch_services(), prepare_gas_mechanism=lambda: None,
        check_equilibrium=None, persist_scorecard=None)
    result = run_reactor_batch_stage(
        candidates, temperatures=config.reactor_temperatures,
        reactor_types=('PFR',), pathway_mode='thermocatalytic_pfr',
        multiphysics_results_dir=str(tmp_path), allow_mock_inputs=False,
        judge_catalyst=config.solids_judge_catalyst,
        headline_t_min=config.solids_headline_t_min,
        headline_t_max=config.solids_headline_t_max, services=services)
    card = result.state['solids_scorecard']
    assert card['judge_catalyst'] == 'ni_np_lit'
    assert card['headline']['PFR']['T_K'] == 973.15
    assert 0 < result.state['best_conversion'] < 1
    rows = result.products['reactor_results']
    assert len(rows) == len(config.reactor_temperatures)
    assert all(row['status'] == 'complete' and not row.get('mock') for row in rows)
    assert candidates.empty


@pytest.mark.parametrize('reactor_type,judge', [
    ('MMBCR', 'ni_np_lit'), ('PFR', None),
])
def test_reference_is_only_added_to_enabled_solids_campaigns(reactor_type, judge):
    def unexpected_reference_load(name):
        raise AssertionError('this campaign must not load the solids reference')

    services = replace(
        default_reactor_batch_services(), prepare_gas_mechanism=lambda: None,
        check_equilibrium=None, persist_scorecard=None,
        load_reference_candidate=unexpected_reference_load)
    result = run_reactor_batch_stage(
        pd.DataFrame(), temperatures=(973.15,), reactor_types=(reactor_type,),
        pathway_mode=reactor_models.SINGLE_REACTOR_MODE[reactor_type],
        multiphysics_results_dir='', allow_mock_inputs=False,
        judge_catalyst=judge, services=services)
    assert result.products['reactor_results'] == []


@pytest.mark.parametrize('swept', [False, True])
def test_different_sweep_jobs_cannot_replace_each_others_mechanisms(
        tmp_path, monkeypatch, swept):
    ct = pytest.importorskip('cantera')
    monkeypatch.setattr(reactor_mechanisms, 'MECHANISMS_DIR', tmp_path / 'shared')
    monkeypatch.setattr(yaml_sweep, 'SWEEPS_DIR', tmp_path / 'sweeps')
    monkeypatch.setattr(yaml_sweep, '_print_table', lambda *args: None)

    def simulate(name, mechanism, temperatures, reactor_types, **kwargs):
        gas = ct.Solution(mechanism, 'gas')
        graphite = ct.Solution(mechanism, 'graphite')
        surface = ct.Interface(mechanism, f'{name}_surface', [gas, graphite])
        assert 'C_encap_s' in surface.species_names
        return [dict(catalyst_name=name, reactor_type='PFR', T_K=temperatures[0],
                     status='complete', CH4_conversion=0.1, surface_loaded=True)]

    monkeypatch.setattr(reactor_models, 'run_reactor_sweep', simulate)
    snapshots = []
    for name, barrier in [('first', 0.7), ('second', 1.3)]:
        activation = '' if swept else f'    E_act: {barrier}\n'
        grid = (f'sweep:\n  kinetics:\n    E_act: [{barrier}]\n'
                '  policy:\n    max_regen_cycles: [0, 3]\n') if swept else ''
        source = tmp_path / f'{name}.yaml'
        source.write_text(
            f'name: {name}\ncatalyst:\n  name: ni_np_lit\n'
            '  material_class: SolidCatalyst\n'
            '  genome: "(\'SolidCatalyst\', \'Ni\', \'SiO2\', \'fcc111\', 0.0, (), 1, 0)"\n'
            '  kinetics:\n' + activation + '    dE_H: -0.5\n' + grid +
            'conditions:\n  temperatures_K: [923.15]\n  reactors: [PFR]\n'
            'cells:\n  - name: production\n    catalyst_particle_mm: 0.13\n'
            '    metal_loading: 0.5\n    metal_dispersion: 0.3\n',
            encoding='utf-8')
        result = yaml_sweep.run_sweep(source)
        assert len(result['mechanism_files']) == 1
        mechanism = Path(result['mechanism_file'])
        assert all(Path(row['mechanism_file']) == mechanism
                   for row in result['records'])
        sidecar = mechanism.with_suffix('.kinetics.json')
        assert json.loads(sidecar.read_text())['inputs']['methane_activation_eV'] == barrier
        snapshots.append((mechanism, mechanism.read_bytes(), sidecar.read_bytes()))

    assert snapshots[0][0] != snapshots[1][0]
    for mechanism, yaml_bytes, metadata_bytes in snapshots:
        assert mechanism.read_bytes() == yaml_bytes
        assert mechanism.with_suffix('.kinetics.json').read_bytes() == metadata_bytes
