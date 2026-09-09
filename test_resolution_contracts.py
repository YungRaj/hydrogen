#!/usr/bin/env python3
"""Contracts for evidence curation and candidate-specific Hamiltonians."""

import hashlib
import json
import tempfile
from pathlib import Path
from unittest.mock import patch


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_curated_prior_art_manifest():
    from pipeline.evidence.prior_art import PriorArtRegistry
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        source = root / 'source.txt'
        source.write_text('immutable source snapshot\n')
        records = []
        for genome, year, split, source_id in (
            (('SAC', 'Fe', 'N4', 'N-graphene', 'OH'), 2020,
             'training', 'doi:training'),
            (('SAC', 'Co', 'N4', 'N-graphene', 'OH'), 2025,
             'holdout', 'doi:holdout')):
            records.append({
                'genome': repr(genome), 'source_type': 'literature',
                'source_id': source_id, 'citation': source_id,
                'evidence_level': 'reported', 'publication_year': year,
                'split': split, 'source_path': source.name,
                'source_sha256': digest(source),
            })
        manifest = root / 'prior_art.json'
        manifest.write_text(json.dumps({
            'schema_version': 1, 'training_cutoff_year': 2022,
            'records': records}))
        registry = PriorArtRegistry(str(root / 'prior.sqlite'))
        report = registry.import_curated_manifest(str(manifest))
        assert report['imported'] == 2
        assert report['training'] == report['holdout'] == 1
        assert registry.count() == 2


def test_candidate_hamiltonian_from_fcidump():
    import numpy as np
    from pyscf.tools import fcidump
    from pipeline.validation.candidate_hamiltonian import (
        build_candidate_hamiltonian)
    from pipeline.validation.vqe_transition_state import exact_ground_energy
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        geometry = root / 'candidate.xyz'
        geometry.write_text('1\ncontract geometry\nH 0 0 0\n')
        source = root / 'FCIDUMP'
        h1 = np.array([[-1.0, 0.05], [0.05, -0.4]])
        h2 = np.zeros((2, 2, 2, 2))
        h2[0, 0, 0, 0] = 0.7
        h2[1, 1, 1, 1] = 0.5
        fcidump.from_integrals(str(source), h1, h2, 2, 2, nuc=0.2, ms=0)
        source.with_suffix('.json').write_text(json.dumps({
            'schema_version': 1, 'candidate_id': 'candidate-1',
            'geometry_id': 'candidate-1/reactant/v1',
            'geometry_path': geometry.name,
            'geometry_sha256': digest(geometry),
            'electronic_structure_protocol': 'contract/CASCI',
            'basis': 'contract-basis', 'active_orbitals': 2,
            'active_electrons': 2, 'charge': 0, 'multiplicity': 1,
            'frozen_orbitals': [0], 'integral_source': 'contract fixture',
            'fermion_to_qubit_mapping': 'Jordan-Wigner',
            'fcidump_sha256': digest(source),
        }))
        result = build_candidate_hamiltonian(source)
        assert result['candidate_id'] == 'candidate-1'
        assert result['n_qubits'] == 4
        assert result['catalyst_specific_hamiltonian'] is True
        assert result['pauli_terms']
        assert isinstance(exact_ground_energy(result['pauli_terms'], 4), float)
        with patch(
                'pipeline.validation.vqe_transition_state.run_vqe',
                return_value={'benchmarked': True, 'mock': False,
                              'energy_Ha': -1.0}), patch(
                'pipeline.validation.vqe_transition_state.save_json'):
            from pipeline.validation.vqe_transition_state import (
                validate_transition_state)
            promoted = validate_transition_state(
                'candidate-1', candidate_hamiltonian=str(source),
                target='qpp-cpu')
        assert promoted['catalyst_specific_hamiltonian'] is True
        assert promoted['evidence_level'] == 'candidate_specific_vqe_benchmarked'


def main():
    tests = (test_curated_prior_art_manifest,
             test_candidate_hamiltonian_from_fcidump)
    for test in tests:
        test()
        print('PASS', test.__name__)


if __name__ == '__main__':
    main()
