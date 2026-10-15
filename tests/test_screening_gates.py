"""Regression contracts for screening-steering gates (Handoff 28/29 fixes)."""

import math

import pandas as pd
import pytest

from pipeline.screening import genetic_optimizer


def test_required_screening_fields_fail_closed():
    row = pd.Series({'temperature_K': 923.15, 'dE_H': float('nan'), 'E_act_raw': 'x'})
    assert genetic_optimizer._required_finite(row, 'temperature_K') == 923.15
    for field in ('dE_H', 'E_act_raw', 'missing'):
        with pytest.raises(ValueError):
            genetic_optimizer._required_finite(row, field)


@pytest.mark.parametrize(
    'flags',
    [{'exceeds_equilibrium': True}, {'carbon_balance_ok': False}],
)
def test_unphysical_pfr_result_is_never_complete(flags):
    result = {'status': 'complete', 'CH4_conversion': 0.99, **flags}
    outcome = genetic_optimizer._full_simulation_outcome(result)
    assert outcome['full_simulation_status'] == 'validation_required'
    assert math.isnan(outcome['full_simulation_score'])
    assert outcome['full_simulation_exceeds_equilibrium'] is bool(
        flags.get('exceeds_equilibrium', False)
    )


def test_physical_pfr_result_is_scored():
    outcome = genetic_optimizer._full_simulation_outcome({
        'status': 'complete', 'CH4_conversion': 0.25,
        'exceeds_equilibrium': False, 'carbon_balance_ok': True,
        'off_site_carbon_active': True,
    })
    assert outcome['full_simulation_status'] == 'complete'
    assert outcome['full_simulation_score'] == pytest.approx(0.75)
    assert outcome['full_simulation_off_site_carbon_active'] is True


def test_thermal_scores_do_not_steer_before_b5_closure(monkeypatch):
    import pipeline.reactors.scorecard as scorecard
    from pipeline.screening.branch_discovery import _gate_thermal_steering

    simulation = pd.DataFrame({
        'genome': ["('SAC', 0)"],
        'full_simulation_score': [0.1],
        'full_simulation_status': ['complete'],
        'full_simulation_exceeds_equilibrium': [False],
    })
    monkeypatch.setattr(scorecard, 'B5_ACCEPTED', False)
    closed = _gate_thermal_steering(simulation)
    assert not {'full_simulation_score', 'full_simulation_status'} & set(closed.columns)
    assert closed['thermal_diagnostic_score'].tolist() == [0.1]
    assert 'full_simulation_exceeds_equilibrium' in closed.columns
    monkeypatch.setattr(scorecard, 'B5_ACCEPTED', True)
    accepted = _gate_thermal_steering(simulation)
    assert accepted['full_simulation_score'].tolist() == [0.1]


def test_b5_closure_is_not_accepted_without_evidence():
    from pipeline.reactors.scorecard import (
        B5_ACCEPTED,
        NI_JUDGE_CATALYST,
        evaluate_b5_closure_gate,
    )

    assert B5_ACCEPTED is False
    incomplete = evaluate_b5_closure_gate({
        'judge_catalyst': NI_JUDGE_CATALYST,
        'headline': {
            'PFR': {'single_pass_CH4_conversion': 0.1456},
            'Fluidized': {'single_pass_CH4_conversion': 0.1694},
        },
    })
    assert incomplete['b5_accepted'] is False
    complete = evaluate_b5_closure_gate(
        {
            'judge_catalyst': NI_JUDGE_CATALYST,
            'headline': {
                'PFR': {'single_pass_CH4_conversion': 0.62, 'exceeds_equilibrium': False},
                'Fluidized': {'single_pass_CH4_conversion': 0.58, 'exceeds_equilibrium': False},
            },
        },
        b66_summary={
            'criteria': {
                k: {'pass': True}
                for k in ('flat', 'linearity', 'encapsulation_onset', 'lifetime', 'yield')
            }
        },
        joint_summary={'declaration': 'simultaneous_hit', 'n_both_roi': 1},
    )
    assert complete['b5_accepted'] is True


def test_tea_feedstock_uses_methane_stoichiometry():
    from pipeline.economics.tea import (
        MMBTU_CH4_PER_KG_H2, SCENARIOS, estimate_hydrogen_cost,
    )

    # ~4 kg CH4 (HHV 55.5 MJ/kg) per kg H2 -> ~0.21 MMBtu.
    assert MMBTU_CH4_PER_KG_H2 == pytest.approx(0.209, abs=0.002)
    base = SCENARIOS['base']
    full = estimate_hydrogen_cost(1.0, 'base')['h2_cost_usd_kg']
    expected = (
        base.natural_gas_usd_mmbtu * MMBTU_CH4_PER_KG_H2
        + base.base_energy_kwh_kg_h2 * base.electricity_usd_kwh
        + base.capex_usd_kg_h2 - base.carbon_value_usd_kg_h2
    )
    assert full == pytest.approx(expected, abs=1e-3)
    # The optimistic scenario must not report a negative hydrogen cost.
    assert estimate_hydrogen_cost(1.0, 'optimistic')['h2_cost_usd_kg'] > 0


def test_archive_pool_never_ranks_objectives_against_each_other(tmp_path):
    import sqlite3

    import numpy as np

    from pipeline.search.exhaustive_search import (
        ScanConfig, load_archive_genomes, run_streaming_scan,
    )

    def scorer(genomes):
        n = len(genomes)
        primary = np.arange(n, dtype=float) + 1.0     # always > 0
        constant = np.zeros(n)                         # carries no ranking signal
        cost = -1000.0 - np.arange(n, dtype=float)     # different, negative scale
        return np.column_stack([primary, constant, cost])

    database = str(tmp_path / 'scan.sqlite')
    run_streaming_scan(
        ScanConfig('app', database, stop=400, batch_size=400, global_archive_size=30),
        scorer,
    )
    with sqlite3.connect(database) as conn:
        best = conn.execute(
            "SELECT genome FROM global_archive ORDER BY primary_score LIMIT 1"
        ).fetchone()[0]
    pool = load_archive_genomes(database, 'app', limit=6)
    import json

    assert pool[0] == tuple(
        tuple(v) if isinstance(v, list) else v for v in json.loads(best)
    )
    assert len(pool) == 6


def test_campaign_vqe_requires_candidate_hamiltonians(tmp_path, monkeypatch):
    import subprocess
    import sys
    from pathlib import Path

    import pipeline.utils

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import run_production_campaign as campaign

    calls = []
    monkeypatch.setattr(pipeline.utils, 'save_json', lambda *a, **k: None)
    monkeypatch.setattr(subprocess, 'run', lambda *a, **k: calls.append(a))
    champions = pd.DataFrame({'candidate_id': ['cand_a', 'cand_b']})

    state = {}
    campaign._run_phase4_vqe(state, champions, tmp_path)
    assert state['phase4']['skipped'] is True
    assert len(state['phase4']['missing_hamiltonians']) == 2
    assert calls == []  # never runs CUDA-Q on a placeholder/toy model

    fcidump = campaign._candidate_hamiltonian_path(tmp_path, 'cand_b')
    fcidump.parent.mkdir(parents=True)
    fcidump.write_text('')
    seen = []
    monkeypatch.setattr(
        campaign, '_run_vqe', lambda name, root, path: seen.append((name, path)) or {}
    )
    state = {}
    campaign._run_phase4_vqe(state, champions, tmp_path)
    assert seen == [('cand_b', fcidump)]
    assert state['phase4']['candidates'] == ['cand_b']


def test_campaign_reactor_judge_reference_is_scorecard_only():
    import sys
    from pathlib import Path
    from types import SimpleNamespace

    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import run_production_campaign as campaign

    assert 923.15 in campaign.REACTOR_TEMPERATURES_K
    assert 973.15 in campaign.REACTOR_TEMPERATURES_K
    calls = []

    def fake_simulate(row, name, temperatures, **kwargs):
        calls.append(name)
        return {'sweep': [{'catalyst_name': name, 'T_K': 923.15}]}

    args = SimpleNamespace(mode='thermocatalytic_pfr', multiphysics_results_dir=None)
    records = []
    campaign._simulate_judge_reference(args, fake_simulate, records)
    assert calls == ['ni_np_lit']
    assert records == [{'catalyst_name': 'ni_np_lit', 'T_K': 923.15}]
    # Idempotent: an existing judge sweep is not simulated twice.
    campaign._simulate_judge_reference(args, fake_simulate, records)
    assert calls == ['ni_np_lit']


def test_campaign_and_orchestrator_share_temperatures_judge_and_selectors():
    from pathlib import Path
    import run_production_campaign as campaign
    from pipeline.orchestrator import PipelineConfig
    from pipeline.reactors.modes import MODE_CHOICES
    from pipeline.reactors.scorecard import (
        NI_JUDGE_CATALYST,
        NI_JUDGE_HEADLINE_T_MAX,
        NI_JUDGE_HEADLINE_T_MIN,
    )

    for mode in MODE_CHOICES:
        cfg = PipelineConfig(pyrolysis_mode=mode)
        assert tuple(campaign._reactor_temperatures(mode)) == cfg.reactor_temperatures
        assert cfg.solids_judge_catalyst == NI_JUDGE_CATALYST
        assert cfg.solids_headline_t_min == NI_JUDGE_HEADLINE_T_MIN
        assert cfg.solids_headline_t_max == NI_JUDGE_HEADLINE_T_MAX
    assert 'select_for_pyrolysis_reactor' in Path(campaign.__file__).read_text()
    assert 'run_report_stage' in Path(campaign.__file__).read_text()

