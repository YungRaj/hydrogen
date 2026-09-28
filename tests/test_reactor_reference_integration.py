"""Reference-campaign configuration and isolated sweep artifacts."""

from dataclasses import replace
import json
import math
from pathlib import Path
import sys
import tempfile
from unittest.mock import patch

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from pipeline.orchestrator import PipelineConfig, normalized_pipeline_config
from pipeline.reactors import mechanisms as reactor_mechanisms
from pipeline.reactors import models as reactor_models
from pipeline.reactors.sweeps import runner as yaml_sweep
from pipeline.reactors.scorecard import build_solids_scorecard
from pipeline.stages.reactor_batch import (
    ReactorBatchServices, default_reactor_batch_services, run_reactor_batch_stage,
)


def test_normalization_preserves_requested_reference_temperatures():
    """Preserve explicit reference temperatures in normal and quick modes."""
    for quick_mode in (False, True):
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


def test_default_campaign_loads_ni_reference_and_produces_a_real_headline():
    """Reproduce the documented Ni reference through the real Cantera path.

    These are regression targets for the B6-5 literature-parameterized model,
    not claims that the reduced mechanism has been experimentally calibrated.
    The explicit evidence-tier assertions keep that distinction fail-closed.
    """
    import cantera  # noqa: F401 - availability is part of this suite's contract

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        with patch.object(reactor_mechanisms, 'MECHANISMS_DIR', tmp_path), \
                patch.object(reactor_models, 'save_json',
                             lambda *args, **kwargs: None):
            config = normalized_pipeline_config(PipelineConfig())
            candidates = pd.DataFrame(columns=['candidate_id', 'E_act'])
            services = replace(
                default_reactor_batch_services(),
                prepare_gas_mechanism=lambda: None,
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

        # Documented B6-5 production-cell results.  A tight but non-bitwise
        # tolerance catches unit, phase, prefactor, or reactor-wiring changes
        # while allowing harmless integrator/platform roundoff.
        reference = {
            923.15: {
                'CH4_conversion': 0.1456,
                'carbon_turnovers_per_site': 7.03,
                'exit_theta_C': 1.8e-5,
                'exit_theta_C_encap': 3.0e-4,
                'c_gamma_to_c_delta_ratio': 2.3e4,
                'filament_yield_gC_per_gMetal_h': 389.0,
            },
            973.15: {
                'CH4_conversion': 0.2633,
                'carbon_turnovers_per_site': 12.1,
                'exit_theta_C': 1.1e-5,
                'exit_theta_C_encap': 3.5e-4,
                'c_gamma_to_c_delta_ratio': 3.4e4,
                'filament_yield_gC_per_gMetal_h': 667.0,
            },
        }
        by_temperature = {row['T_K']: row for row in rows}
        for temperature, expected in reference.items():
            row = by_temperature[temperature]
            for field, target in expected.items():
                assert math.isclose(
                    row[field], target, rel_tol=0.02, abs_tol=1e-8), (
                        temperature, field, row[field], target)

            # Methane carbon must appear in the two explicitly represented
            # condensed-carbon channels; graphite is not an ideal-gas tracer.
            assert row['carbon_phase_model'] == \
                'condensed_graphite_plus_surface_C_s'
            assert row['graphite_loaded'] is True
            assert row['carbon_balance_ok'] is True
            assert row['outfeed_carbon_mol_per_pass'] == 0.0
            assert 0.0 <= row['c_gamma_mol_per_pass'] <= \
                row['solid_carbon_mol_per_pass']
            assert 0.0 <= row['c_delta_mol_per_pass'] <= \
                row['solid_carbon_mol_per_pass']
            # With no circulating outfeed, the reported balance residual is
            # precisely carbon transferred off-site into the graphite sink.
            assert math.isclose(
                row['carbon_balance_residual_mol'],
                row['c_gamma_mol_per_pass'],
                rel_tol=1e-10, abs_tol=1e-14)

            # Exercise the actual packed-bed/surface path and its axial
            # solution, not merely successful Cantera YAML loading.
            assert row['surface_loaded'] is True
            assert row['surface_name'] == 'ni_np_lit_surface'
            assert row['bed_or_interface'] == 'fixed_packed_bed'
            assert row['cantera_reactor_model'] == \
                'staged_lagrangian_ideal_gas_reactors'
            assert row['conversion_basis'] == 'argon_tracer'
            assert row['co2_permitted'] is False
            assert all(
                later >= earlier
                for earlier, later in zip(
                    row['conversion_profile'], row['conversion_profile'][1:]))

            # Deactivation and off-site carbon remain reported, while this
            # uncalibrated screening mechanism is forbidden from presenting
            # itself as validated predictive reactor evidence.
            assert row['off_site_carbon_active'] is True
            assert row['regen_cycles_completed'] == 0
            assert row['encapsulation_onset'] is False
            assert row['kinetics_status'] == 'screening_template_incomplete'
            assert row['reactor_evidence_tier'] == \
                'diagnostic_screening_template'
            assert 'incomplete_candidate_kinetics' in \
                row['reactor_evidence_limitations']
            assert 'pfr_screening_model_not_validated' in \
                row['reactor_evidence_limitations']

        assert by_temperature[973.15]['CH4_conversion'] > \
            by_temperature[923.15]['CH4_conversion']


def test_reference_is_only_added_to_enabled_solids_campaigns():
    """Do not inject the Ni judge into MMBCR or judge-disabled campaigns."""
    def unexpected_reference_load(name):
        raise AssertionError('this campaign must not load the solids reference')

    for reactor_type, judge in (('MMBCR', 'ni_np_lit'), ('PFR', None)):
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


def test_different_sweep_jobs_cannot_replace_each_others_mechanisms():
    """Keep every static or swept job's mechanism immutable and distinct."""
    import cantera as ct

    with tempfile.TemporaryDirectory() as tmp:
        tmp_path = Path(tmp)
        for swept in (False, True):
            _exercise_isolated_sweep(tmp_path / str(swept), swept, ct)


def _exercise_isolated_sweep(tmp_path: Path, swept: bool, ct) -> None:
    """Execute one isolated mechanism test case without pytest fixtures."""
    tmp_path.mkdir()

    def simulate(name, mechanism, temperatures, reactor_types, **kwargs):
        gas = ct.Solution(mechanism, 'gas')
        graphite = ct.Solution(mechanism, 'graphite')
        surface = ct.Interface(mechanism, f'{name}_surface', [gas, graphite])
        assert 'C_encap_s' in surface.species_names
        return [dict(catalyst_name=name, reactor_type='PFR', T_K=temperatures[0],
                     status='complete', CH4_conversion=0.1, surface_loaded=True)]

    with patch.object(reactor_mechanisms, 'MECHANISMS_DIR', tmp_path / 'shared'), \
            patch.object(yaml_sweep, 'SWEEPS_DIR', tmp_path / 'sweeps'), \
            patch.object(yaml_sweep, '_print_table', lambda *args: None), \
            patch.object(reactor_models, 'run_reactor_sweep', simulate):
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
            assert json.loads(sidecar.read_text())[\
                'inputs']['methane_activation_eV'] == barrier
            snapshots.append((mechanism, mechanism.read_bytes(), sidecar.read_bytes()))

    assert snapshots[0][0] != snapshots[1][0]
    for mechanism, yaml_bytes, metadata_bytes in snapshots:
        assert mechanism.read_bytes() == yaml_bytes
        assert mechanism.with_suffix('.kinetics.json').read_bytes() == metadata_bytes


def main() -> None:
    """Run reference-integration contracts without a pytest dependency."""
    tests = [value for name, value in sorted(globals().items())
             if name.startswith('test_') and callable(value)]
    for test in tests:
        test()
        print('PASS', test.__name__)
    print(f'{len(tests)}/{len(tests)} reactor reference contracts passed')


if __name__ == '__main__':
    main()
