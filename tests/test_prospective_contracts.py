"""Dependency-light contracts for prospective scoring and training gates."""

from pathlib import Path
import hashlib
import json
import sys
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd
from ase import Atoms

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

import run_divide_conquer_pilot as campaign
from pipeline.evidence.orr_dft_slate import verify_slate
from pipeline.evidence.orr_quality_diagnostic import verify_diagnostic
from pipeline.evidence.pyrolysis_benchmark import (
    analyze_screening_results,
    dataset_digest,
    effective_barrier_eV,
    energetic_span_eV,
    h2_chemical_potential_eV,
    rate_free_energy_span_eV,
    roc_auc,
    validate_benchmark,
    verify_analysis_lock,
)
from pipeline.evidence.molten_neb import NEBConfig, _active_site, _path_continuous
from pipeline.search.discovery import candidate_id
from pipeline.screening.small_data_ranker import screening_metric_eligibility
from pipeline.screening.small_data_ranker import training_metric_eligibility


class ProspectiveContracts(unittest.TestCase):
    def test_molten_neb_selects_promoter_and_rejects_path_jumps(self) -> None:
        structure = Atoms(
            'BiCuBi',
            positions=((0.0, 0.0, 0.0), (1.0, 1.0, 0.0), (2.0, 2.0, 0.0)),
            cell=(10.0, 10.0, 10.0),
            pbc=True,
        )
        self.assertEqual(_active_site(structure, [0, 1, 2], 'Cu'), 1)
        close = [structure.copy(), structure.copy()]
        close[1].positions += 0.1
        self.assertTrue(_path_continuous(close, NEBConfig().max_image_displacement_A))
        close[1].positions[0, 0] += 2.0
        self.assertFalse(_path_continuous(close, NEBConfig().max_image_displacement_A))

    def test_failed_molten_neb_evidence_preserves_its_historical_lock(self) -> None:
        lock_path = REPO_ROOT / 'docs/evidence/molten_neb_validation_lock.json'
        lock = json.loads(lock_path.read_text())
        self.assertEqual(lock['model'], 'esen-sm-conserving-all-oc25')
        self.assertEqual(lock['acceptance']['minimum_valid_paths'], 4)
        self.assertIn('Do not read or score', lock['held_out_rule'])
        result = json.loads((
            REPO_ROOT / 'results/evidence/molten_neb_validation.json'
        ).read_text())
        self.assertEqual(
            result['lock_sha256'], hashlib.sha256(lock_path.read_bytes()).hexdigest()
        )
        self.assertEqual(result['attempted_paths'], 5)
        self.assertEqual(result['valid_paths'], 0)
        self.assertFalse(result['acceptance_passed'])
        self.assertFalse(result['held_out_conditions_consumed'])

    def test_ni111_neb_workflow_lock_matches_current_implementation(self) -> None:
        lock = json.loads((
            REPO_ROOT / 'docs/evidence/ni111_mlip_neb_validation_lock.json'
        ).read_text())
        implementation = REPO_ROOT / 'pipeline/evidence/molten_neb.py'
        self.assertEqual(lock['kind'], 'ni111_reference')
        self.assertEqual(lock['reference_barrier_eV'], 1.05)
        self.assertEqual(lock['neb_config']['images'], 11)
        self.assertEqual(lock['neb_config']['plain_path_steps'], 1000)
        self.assertEqual(
            hashlib.sha256(implementation.read_bytes()).hexdigest(),
            lock['implementation_sha256'],
        )
        result = json.loads((
            REPO_ROOT / 'results/evidence/ni111_mlip_neb_validation.json'
        ).read_text())
        self.assertEqual(result['lock_sha256'], hashlib.sha256(
            (REPO_ROOT / 'docs/evidence/ni111_mlip_neb_validation_lock.json').read_bytes()
        ).hexdigest())
        self.assertFalse(result['acceptance_passed'])
        self.assertFalse(result['held_out_conditions_consumed'])
        self.assertTrue(result['record']['valid'])
        self.assertTrue(result['record']['climbing_converged'])
        self.assertGreater(result['record']['absolute_error_eV'], 0.5)

    def test_pyrolysis_benchmark_is_complete_and_locked(self) -> None:
        benchmark = json.loads((
            REPO_ROOT / 'docs/evidence/pyrolysis_benchmark_dataset.json'
        ).read_text())
        self.assertEqual(validate_benchmark(benchmark), [])
        self.assertEqual(len(benchmark['records']), 24)
        self.assertTrue(all(row['exact_encoding'] for row in benchmark['records']))

    def test_complete_pyrolysis_benchmark_contract_can_lock(self) -> None:
        records = []
        for index in range(24):
            label = 'active' if index < 12 else 'poor_or_poisoning'
            genome = (
                'SolidCatalyst', 'Ni', 'Al2O3', 'fcc111',
                -0.10 + index * 0.008, ('B', 'N'), 1, 0,
            )
            records.append({
                'condition_id': f'condition-{index}',
                'catalyst_identity': f'synthetic catalyst {index}',
                'temperature_K': 900.0,
                'reactor_context': 'synthetic contract fixture',
                'source_id': f'doi:10.0000/fixture-{index}',
                'source_family': f'source-{index}',
                'label': label,
                'label_basis': 'synthetic contract fixture only',
                'source_url': f'https://doi.org/10.0000/fixture-{index}',
                'exact_encoding': True,
                'genome': repr(genome),
            })
        document = {'schema_version': 1, 'locked': True, 'records': records}
        document['dataset_sha256'] = dataset_digest(document)
        self.assertEqual(validate_benchmark(document), [])

        duplicate = {**document, 'records': [dict(row) for row in records]}
        duplicate['records'][1]['condition_id'] = duplicate['records'][0]['condition_id']
        duplicate['dataset_sha256'] = dataset_digest(duplicate)
        self.assertTrue(any(
            'duplicates condition_id' in error
            for error in validate_benchmark(duplicate)
        ))

        incomplete = {**document, 'records': [dict(row) for row in records]}
        del incomplete['records'][0]['temperature_K']
        incomplete['dataset_sha256'] = dataset_digest(incomplete)
        self.assertTrue(any(
            'temperature_K' in error for error in validate_benchmark(incomplete)
        ))

        unsupported = {**document, 'records': [dict(row) for row in records]}
        unsupported['records'][0]['genome'] = repr(
            ('MoltenMetal', 'Hg', 'None', 0.0, 973.15)
        )
        unsupported['dataset_sha256'] = dataset_digest(unsupported)
        self.assertTrue(any(
            'outside the design space' in error
            for error in validate_benchmark(unsupported)
        ))

        mismatched_source = {**document, 'records': [dict(row) for row in records]}
        mismatched_source['records'][0]['source_url'] = 'https://doi.org/10.0000/wrong'
        mismatched_source['dataset_sha256'] = dataset_digest(mismatched_source)
        self.assertTrue(any(
            'source_id and source_url disagree' in error
            for error in validate_benchmark(mismatched_source)
        ))

    def test_pyrolysis_energetic_span_penalizes_product_poisoning(self) -> None:
        balanced = {'E_act_raw': 0.8, 'dE_H': -0.2, 'dE_C': -0.3}
        hydrogen_parked = {'E_act_raw': 0.1, 'dE_H': -2.0, 'dE_C': -0.3}
        carbon_parked = {'E_act_raw': 0.1, 'dE_H': -0.2, 'dE_C': -3.0}
        self.assertEqual(energetic_span_eV(balanced), 0.8)
        self.assertEqual(energetic_span_eV(hydrogen_parked), 4.0)
        self.assertEqual(energetic_span_eV(carbon_parked), 3.0)
        with self.assertRaisesRegex(ValueError, 'energetic-span input'):
            energetic_span_eV({'E_act': 0.8, 'dE_H': -0.2, 'dE_C': -0.3})

    def test_rate_average_gives_fast_sites_their_arrhenius_weight(self) -> None:
        effective = effective_barrier_eV([0.4, 1.5, 1.5, 1.5, 1.5], 1273.15)
        self.assertLess(effective, 0.6)
        self.assertGreater(effective, 0.4)
        self.assertAlmostEqual(effective_barrier_eV([0.8] * 5, 1000.0), 0.8)

    def test_h2_chemical_potential_uses_locked_one_bar_standard_state(self) -> None:
        self.assertAlmostEqual(h2_chemical_potential_eV(1000.0), -1.42061, places=5)
        self.assertAlmostEqual(h2_chemical_potential_eV(1300.0), -1.95043, places=5)
        with self.assertRaisesRegex(ValueError, '298 <= T <= 2500'):
            h2_chemical_potential_eV(200.0)

    def test_rate_free_energy_span_rejects_unvalidated_pure_hosts(self) -> None:
        replicas = json.dumps([
            {'E_act_raw': 0.4, 'dE_H': -0.5, 'dE_C': 1.0},
            {'E_act_raw': 1.5, 'dE_H': -0.5, 'dE_C': 1.0},
        ])
        alloy = {
            'genome': "('MoltenMetal', 'Bi', 'Ni', 27.0, 1273.15)",
            'temperature_K': 1273.15,
            'molten_replicas_json': replicas,
        }
        self.assertLess(rate_free_energy_span_eV(alloy), 0.6)
        pure = {**alloy, 'genome': "('MoltenMetal', 'Bi', 'None', 0.0, 1273.15)"}
        with self.assertRaisesRegex(ValueError, 'outside validated scope'):
            rate_free_energy_span_eV(pure)

    def test_roc_auc_handles_ties(self) -> None:
        self.assertEqual(roc_auc([True, False], [1.0, 0.0]), 1.0)
        self.assertEqual(roc_auc([True, False], [0.0, 1.0]), 0.0)
        self.assertEqual(roc_auc([True, False], [1.0, 1.0]), 0.5)

    def test_descriptor_analysis_fails_closed_on_missing_screening(self) -> None:
        benchmark = json.loads((
            REPO_ROOT / 'docs/evidence/pyrolysis_benchmark_dataset.json'
        ).read_text())
        with self.assertRaisesRegex(ValueError, 'missing screening result'):
            analyze_screening_results(benchmark, [])

    def test_pyrolysis_analysis_lock_matches_sources(self) -> None:
        lock = REPO_ROOT / 'docs/evidence/pyrolysis_descriptor_benchmark_lock.json'
        self.assertEqual(verify_analysis_lock(lock, REPO_ROOT), [])

    def test_rate_descriptor_result_preserves_lock_and_screening_input(self) -> None:
        result = json.loads((
            REPO_ROOT / 'docs/evidence/pyrolysis_rate_descriptor_result.json'
        ).read_text())
        paths = {
            'analysis_lock_sha256': 'docs/evidence/pyrolysis_descriptor_benchmark_lock.json',
            'screening_csv_sha256': 'results/screening/pyrolysis_rate_descriptor.csv',
        }
        for key, relative in paths.items():
            actual = hashlib.sha256((REPO_ROOT / relative).read_bytes()).hexdigest()
            self.assertEqual(result[key], actual)
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            lock = REPO_ROOT / 'docs/evidence/pyrolysis_descriptor_benchmark_lock.json'
            altered = json.loads(lock.read_text())
            relative = 'pipeline/evidence/pyrolysis_benchmark.py'
            target = root / relative
            target.parent.mkdir(parents=True)
            target.write_text('changed')
            altered['sha256'] = {relative: altered['sha256'][relative]}
            copied_lock = root / 'lock.json'
            copied_lock.write_text(json.dumps(altered))
            self.assertEqual(
                verify_analysis_lock(copied_lock, root),
                [f'locked file changed: {relative}'],
            )

    def test_orr_quality_diagnostic_is_checksum_locked(self) -> None:
        source = REPO_ROOT / 'results/prospective_search/orr_quality_diagnostic.json'
        self.assertTrue(verify_diagnostic(source))
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / 'diagnostic.json'
            document = json.loads(source.read_text())
            document['totals']['invalid'] += 1
            altered.write_text(json.dumps(document))
            self.assertFalse(verify_diagnostic(altered))

    def test_confirmatory_orr_dft_slate_is_checksum_locked(self) -> None:
        source = REPO_ROOT / 'results/dft/orr_confirmatory_slate.json'
        self.assertTrue(verify_slate(source))
        with tempfile.TemporaryDirectory() as directory:
            altered = Path(directory) / 'slate.json'
            document = json.loads(source.read_text())
            document['records'][0]['predicted_overpotential_V'] += 0.1
            altered.write_text(json.dumps(document))
            self.assertFalse(verify_slate(altered))

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

    def test_pyrolysis_reactor_handoff_fails_closed_but_validation_rescues(self) -> None:
        from pipeline.screening.stage_selection import (
            select_for_pyrolysis_reactor,
            select_for_validation,
        )

        frame = pd.DataFrame({
            'material_class': ['SAC', 'DAC', 'SolidCatalyst'],
            'valid': [True, True, True],
            'E_act': [0.01, 0.20, 0.60],
            'E_act_censored': [True, False, False],
            'pyrolysis_viable': [True, False, True],
            'needs_dft_validation': [True, False, False],
        })
        reactor = select_for_pyrolysis_reactor(frame, 3)
        validation = select_for_validation(frame, 3, 'E_act')
        self.assertEqual(reactor['material_class'].tolist(), ['SolidCatalyst'])
        self.assertEqual(set(validation['material_class']), {'SAC', 'DAC', 'SolidCatalyst'})

        missing_gates = frame.drop(columns=['E_act_censored', 'pyrolysis_viable'])
        self.assertTrue(select_for_pyrolysis_reactor(missing_gates, 3).empty)

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
