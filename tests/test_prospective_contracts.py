"""Dependency-light contracts for prospective scoring and training gates."""

from pathlib import Path
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import run_divide_conquer_pilot as campaign
from pipeline.search.discovery import candidate_id
from pipeline.screening.small_data_ranker import screening_metric_eligibility
from pipeline.screening.small_data_ranker import training_metric_eligibility


class ProspectiveContracts(unittest.TestCase):
    def test_pool_deduplicates_canonical_candidate_ids(self) -> None:
        first = ('MoltenMetal', 'In', 'Dy', 0.0, 950)
        duplicate = ('MoltenMetal', 'In', 'In', 0.5, 950)
        self.assertEqual(candidate_id(first), candidate_id(duplicate))
        with patch.object(campaign, 'CLASS_ORDER', ('MoltenMetal',)), \
             patch.object(campaign, 'CLASS_SIZES', {'MoltenMetal': 3}), \
             patch.object(campaign, '_prior', return_value=set()), \
             patch.object(campaign, '_probe_indices', return_value=[0, 1, 2]), \
             patch.object(
                 campaign,
                 'candidate_at_class',
                 side_effect=(first, duplicate, ('MoltenMetal', 'Ga', 'In', 0.5, 950)),
             ):
            pool = campaign._pool(2)
        self.assertEqual(len(pool), 2)
        self.assertEqual(len({candidate_id(genome) for genome in pool}), 2)

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
        self.assertEqual(
            training_metric_eligibility(frame, 'turquoise_hydrogen').tolist(),
            [True, False, True, False],
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

    def test_policy_lock_is_enforced(self) -> None:
        original = campaign.POLICY_LOCK
        with tempfile.TemporaryDirectory() as directory:
            campaign.POLICY_LOCK = Path(directory) / 'lock.json'
            try:
                expected = campaign._policy_source_hash()
                campaign.POLICY_LOCK.write_text(json.dumps({
                    'policy_source_sha256': expected,
                }))
                self.assertEqual(campaign._require_locked_policy(), expected)
                campaign.POLICY_LOCK.write_text(json.dumps({
                    'policy_source_sha256': '0' * 64,
                }))
                with self.assertRaisesRegex(RuntimeError, 'hash mismatch'):
                    campaign._require_locked_policy()
            finally:
                campaign.POLICY_LOCK = original

    def test_policy_hash_is_independent_of_working_directory(self) -> None:
        expected = campaign._policy_source_hash()
        previous = Path.cwd()
        with tempfile.TemporaryDirectory() as directory:
            try:
                __import__('os').chdir(directory)
                self.assertEqual(campaign._policy_source_hash(), expected)
            finally:
                __import__('os').chdir(previous)
        self.assertEqual(campaign._require_locked_policy(), expected)

    def test_confirmatory_progress_excludes_smoke_and_mismatched_hash(self) -> None:
        original_root, original_lock = campaign.ROOT, campaign.POLICY_LOCK
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            campaign.ROOT = root / 'campaign'
            campaign.POLICY_LOCK = root / 'lock.json'
            digest = campaign._policy_source_hash()
            campaign.POLICY_LOCK.write_text(json.dumps({
                'policy_source_sha256': digest,
                'confirmatory_applications': ['fuel_cell_orr'],
                'exploratory_applications': ['turquoise_hydrogen'],
                'planned_batches': 10,
                'random_one_sided_alpha': 0.05,
            }))
            try:
                rows = []
                for batch, role, policy_hash in (
                    ('prereg-01', 'confirmatory', digest),
                    ('smoke-v4-01', 'smoke', digest),
                    ('prereg-old', 'confirmatory', '0' * 64),
                ):
                    folder = campaign.ROOT / 'batches' / batch
                    folder.mkdir(parents=True)
                    manifest = folder / 'manifest.json'
                    manifest.write_text(json.dumps({
                        'batch_id': batch,
                        'evidence_role': role,
                        'policy_source_sha256': policy_hash,
                    }))
                    analysis = folder / 'analysis.json'
                    analysis.write_text(json.dumps({
                        'batch_id': batch,
                        'manifest_sha256': campaign._hash(manifest),
                    }))
                    rows.append({
                        'batch_id': batch,
                        'manifest_sha256': campaign._hash(manifest),
                    })
                campaign.ROOT.mkdir(exist_ok=True)
                (campaign.ROOT / 'batches.jsonl').write_text(
                    ''.join(json.dumps(row) + '\n' for row in rows)
                )
                output = campaign.confirmatory_analysis()
                progress = json.loads(output.read_text())
                self.assertEqual(progress['batch_ids'], ['prereg-01'])
                self.assertEqual(progress['completed_batches'], 1)
                self.assertTrue(progress['inferential_statistics_withheld'])
                self.assertNotIn('p_value', output.read_text())
            finally:
                campaign.ROOT, campaign.POLICY_LOCK = original_root, original_lock

    def test_confirmatory_verdict_waits_for_exactly_ten_orr_batches(self) -> None:
        original_root, original_lock = campaign.ROOT, campaign.POLICY_LOCK
        original_outcome = campaign._outcome
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            campaign.ROOT = root / 'campaign'
            campaign.POLICY_LOCK = root / 'lock.json'
            campaign._outcome = lambda batch, app: root / f'{batch}-{app}.csv'
            digest = campaign._policy_source_hash()
            campaign.POLICY_LOCK.write_text(json.dumps({
                'policy_source_sha256': digest,
                'confirmatory_applications': ['fuel_cell_orr'],
                'exploratory_applications': ['turquoise_hydrogen'],
                'planned_batches': 10,
                'random_one_sided_alpha': 0.05,
            }))
            try:
                ledger = []
                for index in range(10):
                    batch = f'prereg-{index + 1:02d}'
                    folder = campaign.ROOT / 'batches' / batch
                    folder.mkdir(parents=True)
                    manifest = folder / 'manifest.json'
                    manifest.write_text(json.dumps({
                        'batch_id': batch,
                        'evidence_role': 'confirmatory',
                        'policy_source_sha256': digest,
                        'records': [{'application': 'fuel_cell_orr'}],
                    }))
                    outcome_hashes = {}
                    for app in ('fuel_cell_orr',):
                        outcome = campaign._outcome(batch, app)
                        outcome.write_text('value\n1\n')
                        outcome_hashes[app] = campaign._hash(outcome)
                    analysis = folder / 'analysis.json'
                    analysis.write_text(json.dumps({
                        'batch_id': batch,
                        'manifest_sha256': campaign._hash(manifest),
                        'outcome_sha256': outcome_hashes,
                    }))
                    ledger.append({
                        'batch_id': batch,
                        'manifest_sha256': campaign._hash(manifest),
                    })
                campaign.ROOT.mkdir(exist_ok=True)
                (campaign.ROOT / 'batches.jsonl').write_text(
                    ''.join(json.dumps(row) + '\n' for row in ledger)
                )
                fake_result = {
                    'policy_hits': {'catalyst': 2, 'uncertainty': 1, 'validity': 0}
                }
                with patch.object(
                    campaign,
                    '_score_record',
                    return_value=(fake_result, [], __import__('numpy').zeros(50)),
                ):
                    output = campaign.confirmatory_analysis()
                report = json.loads(output.read_text())
                self.assertEqual(report['status'], 'complete')
                self.assertEqual(len(report['batch_ids']), 10)
                self.assertTrue(report['acceptance_passed'])
                self.assertEqual(
                    [row['application'] for row in report['pooled']],
                    ['fuel_cell_orr'],
                )
            finally:
                campaign.ROOT, campaign.POLICY_LOCK = original_root, original_lock
                campaign._outcome = original_outcome


if __name__ == '__main__':
    unittest.main()
