"""Dependency-light contracts for prospective scoring and training gates."""

from pathlib import Path
import sys
import tempfile
import unittest

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import run_divide_conquer_pilot as campaign
from pipeline.search.discovery import candidate_id
from pipeline.screening.small_data_ranker import screening_metric_eligibility


class ProspectiveContracts(unittest.TestCase):
    def test_censored_and_nonviable_pyrolysis_are_ineligible(self) -> None:
        frame = pd.DataFrame({
            'valid': [True, True, True, False],
            'E_act': [0.01, 0.02, 0.5, 0.1],
            'E_act_censored': [True, False, False, False],
            'pyrolysis_viable': [True, False, True, True],
        })
        self.assertEqual(
            screening_metric_eligibility(frame, 'turquoise_hydrogen').tolist(),
            [False, False, True, False],
        )

    def test_tied_hit_class_is_bounded_by_top_k(self) -> None:
        genomes = [
            ('MetalFreeCarbon', 'graphitic', 0.01 + index / 1000, 'none', 'graphene', 'N')
            for index in range(10)
        ]
        identities = [candidate_id(genome) for genome in genomes]
        frame = pd.DataFrame({
            'genome': [repr(genome) for genome in genomes],
            'valid': True,
            'E_act': [0.5] * 10,
            'E_act_censored': False,
            'pyrolysis_viable': True,
        })
        record = {
            'application': 'turquoise_hydrogen',
            'eligible_ids': identities,
            'budget': 2,
            'policy_selected_ids': {
                name: identities[:2]
                for name in ('catalyst', 'uncertainty', 'validity')
            },
            'random_seed': 7,
            'random_trials': 20,
        }
        result, _, _ = campaign._score_record(record, frame, 'portable')
        self.assertEqual(result['hit_count'], 2)
        self.assertLessEqual(max(result['policy_hits'].values()), 2)

    def test_missing_training_baseline_fails_clearly(self) -> None:
        original = campaign.BASELINE_OUTCOMES
        with tempfile.TemporaryDirectory() as directory:
            campaign.BASELINE_OUTCOMES = {
                **original,
                'turquoise_hydrogen': Path(directory) / 'missing.csv',
            }
            try:
                with self.assertRaisesRegex(RuntimeError, 'missing tracked'):
                    campaign._training('turquoise_hydrogen')
            finally:
                campaign.BASELINE_OUTCOMES = original

    def test_observation_ledger_is_idempotent_and_conflict_safe(self) -> None:
        original = campaign.ROOT
        row = {
            'batch_id': 'batch',
            'application': 'fuel_cell_orr',
            'candidate_id': 'candidate',
            'valid': True,
            'outcome': 0.4,
        }
        with tempfile.TemporaryDirectory() as directory:
            campaign.ROOT = Path(directory)
            try:
                campaign._commit_observations([row])
                campaign._commit_observations([row])
                self.assertEqual(
                    len((campaign.ROOT / 'outcomes.jsonl').read_text().splitlines()),
                    1,
                )
                with self.assertRaisesRegex(RuntimeError, 'conflicting observation'):
                    campaign._commit_observations([{**row, 'outcome': 0.5}])
            finally:
                campaign.ROOT = original


if __name__ == '__main__':
    unittest.main()
