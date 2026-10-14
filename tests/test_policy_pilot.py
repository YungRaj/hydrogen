import pandas as pd

from pipeline.evidence import policy_pilot


def test_policy_pilot_reports_all_policy_cohorts(tmp_path, monkeypatch):
    source = tmp_path / 'screening.csv'
    pd.DataFrame(
        {
            'genome': ["('A', 0)", "('B', 0)", "('C', 0)"],
            'E_act': [0.3, 0.2, 0.1],
            'model_confidence': [0.2, 0.5, 0.8],
            'valid': [True, True, True],
            'pyrolysis_viable': [True, True, True],
        }
    ).to_csv(source, index=False)

    def fake_simulator(genomes, frame):
        return pd.DataFrame(
            {
                'genome': frame['genome'],
                'full_simulation_score': [0.3, 0.2, 0.1],
                'full_simulation_status': ['complete'] * len(frame),
            }
        )

    monkeypatch.setattr(policy_pilot, '_thermal_full_simulator', fake_simulator)
    result = policy_pilot.compare_policies(source, cohort_size=2, draws=10)
    assert result['status'] == 'retrospective_pilot'
    assert set(result['mean_score_lower_is_better']) == {
        'random', 'eSen', 'thermal_oracle_ceiling', 'uncertainty', 'validity'
    }
    assert result['cohort_size'] == 2
    assert 'thermal_lift_vs_random' not in result


def test_fresh_pool_novelty_ledger_excludes_previous_genomes(tmp_path, monkeypatch):
    candidates = iter([('old',), ('new-a',), ('new-b',)])
    monkeypatch.setattr(policy_pilot, 'generate_random_genome', lambda _: next(candidates))
    monkeypatch.setattr(
        policy_pilot,
        'run_screening',
        lambda pool, **_: pd.DataFrame({'genome': [repr(value) for value in pool]}),
    )
    ledger = tmp_path / 'seen.txt'
    ledger.write_text("('old',)\n")
    output = policy_pilot.create_fresh_pool(tmp_path / 'pool.csv', 2, novelty_ledger=ledger)
    assert set(pd.read_csv(output)['genome']) == {"('new-a',)", "('new-b',)"}
