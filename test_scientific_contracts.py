#!/usr/bin/env python3
"""Independent scientific and computational contracts for the pipeline.

This suite intentionally tests equations from first principles instead of
reusing the implementation under test.  Expensive eSen checks are enabled with
RUN_ESEN_CONTRACTS=1; CUDA-Q is tested separately in quantum-env when present.
"""

from __future__ import annotations

import math
import os
import tempfile
import queue
import threading
import time
from pathlib import Path
from unittest.mock import patch

import numpy as np


def test_screening_relaxations_record_and_enforce_force_convergence():
    from ase import Atoms
    from ase.calculators.calculator import Calculator, all_changes
    from pipeline.screening.relaxation import require_relaxation

    class Harmonic(Calculator):
        implemented_properties = ['energy', 'forces']
        def calculate(self, atoms=None, properties=('energy',),
                      system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            positions = np.asarray(atoms.positions)
            self.results = {'energy': float(0.5 * (positions ** 2).sum()),
                            'forces': -positions}

    converging = Atoms('H', positions=[[0.2, 0, 0]])
    converging.calc = Harmonic()
    result = {'valid': False}
    assert require_relaxation(result, converging, 'contract', 0.01, 50)
    assert result['relax_contract_converged'] is True
    assert result['relax_contract_max_force_eV_A'] <= 0.01
    assert len(result['relax_contract_geometry_sha256']) == 64

    incomplete = Atoms('H', positions=[[1.0, 0, 0]])
    incomplete.calc = Harmonic()
    result = {'valid': True}
    assert not require_relaxation(
        result, incomplete, 'contract', 0.01, 0, recovery=False)
    assert result['valid'] is False and result['needs_dft_validation'] is True
    assert result['relax_contract_termination'] == 'recovery_exhausted'


def test_worker_supervisor_requeues_leased_task_and_writes_manifest():
    from pipeline.screening.worker_supervisor import collect_results, emit

    class FakeProcess:
        pid = 123
        exitcode = None
        alive = True
        def is_alive(self): return self.alive
        def terminate(self): self.alive = False; self.exitcode = -15
        def join(self, timeout=None): pass

    statuses, tasks, stopped = queue.Queue(), queue.Queue(), threading.Event()
    initial = FakeProcess()
    emit(statuses, 'ready', 0)
    emit(statuses, 'started', 0, 0)
    emit(statuses, 'fatal', 0, 'injected failure')

    def restart(worker_id):
        replacement = FakeProcess()
        def complete():
            index, genome = tasks.get(timeout=1)
            emit(statuses, 'ready', worker_id)
            emit(statuses, 'started', worker_id, index)
            emit(statuses, 'result', worker_id, (index, {
                'genome': repr(genome), 'valid': True}))
        threading.Thread(target=complete, daemon=True).start()
        return replacement

    with tempfile.TemporaryDirectory() as tmp:
        manifest = Path(tmp) / 'health.json'
        results = collect_results(
            statuses, tasks, stopped, {0: initial}, restart, [('SAC',)],
            'contract', manifest, startup_timeout_s=1,
            heartbeat_timeout_s=1, result_timeout_s=2)
        payload = __import__('json').loads(manifest.read_text())
    assert results[0]['valid'] is True
    assert payload['status'] == 'complete' and payload['completed'] == 1
    assert any(event['kind'] == 'worker_restarted' for event in payload['events'])


def test_candidate_cantera_mechanism_records_kinetics_provenance():
    from pipeline.process.reactor_mechanisms import (
        CandidateKinetics, write_full_mechanism)

    row = {
        'E_act': 0.72, 'dE_H': -0.30, 'dE_CH3': -0.55, 'dE_C': -1.10,
        'screening_protocol': 'contract:relax-v3',
    }
    kinetics = CandidateKinetics.from_screening_row(row, candidate_id='abc')
    path = write_full_mechanism('contract_candidate', kinetics=kinetics)
    text = path.read_text()
    assert 'species: [CH4, H2, C2H2, C2H4, C2H6, Ar, C_graphite]' in text
    assert 'adjacent-phases: [gas]' in text
    assert 'carbon_phase_model' not in text
    metadata = __import__('json').loads(
        path.with_suffix('.kinetics.json').read_text())
    inputs = metadata['inputs']
    assert inputs['methane_activation_eV'] == 0.72
    assert inputs['h_adsorption_eV'] == -0.30
    assert inputs['ch3_adsorption_eV'] == -0.55
    assert inputs['c_adsorption_eV'] == -1.10
    assert inputs['quantitative_status'] == 'screening_template_incomplete'
    assert inputs['provenance']['ch2_dehydrogenation_eV'] == 'template_default'
    from pipeline.process.reactor_models import ReactorConfig, _kinetics_evidence
    evidence = _kinetics_evidence(ReactorConfig(
        mechanism_file=str(path), catalyst_name='contract_candidate'))
    assert evidence['reactor_evidence_tier'] == 'diagnostic_screening_template'
    assert evidence['can_exclude_candidate'] is False

    try:
        import cantera as ct
    except ImportError:
        return
    gas = ct.Solution(str(path), 'gas')
    surface = ct.Interface(str(path), 'contract_candidate_surface', [gas])
    assert 'C_graphite' in gas.species_names
    reactions = surface.reactions()
    assert all(reaction.reversible for reaction in reactions[:-1])
    assert reactions[-1].reversible is False


def test_stage_selection_rescues_incomplete_evidence_without_feeding_reactor():
    import pandas as pd
    from pipeline.screening.stage_selection import (
        select_for_reactor, select_for_validation)

    frame = pd.DataFrame([
        {'genome': "('Z', 0)", 'material_class': 'Z', 'valid': False,
         'error': 'Unconverged clean relaxation', 'needs_dft_validation': True},
        {'genome': "('A', 0)", 'material_class': 'A', 'valid': True,
         'E_act': 0.20},
        {'genome': "('B', 0)", 'material_class': 'B', 'valid': True,
         'E_act': 0.01, 'E_act_censored': True,
         'needs_dft_validation': True},
        {'genome': "('C', 0)", 'material_class': 'C', 'valid': False,
         'error': 'Contains toxic/radioactive elements: Hg',
         'candidate_disposition': 'hard_excluded'},
        {'genome': "('D', 0)", 'material_class': 'D', 'valid': True,
         'E_act': 0.50},
    ])
    reactor = select_for_reactor(frame, 10, 'E_act')
    assert set(reactor.material_class) == {'A', 'D'}
    assert reactor.candidate_disposition.eq('quantitative_screening').all()
    validation = select_for_validation(frame, 10, 'E_act')
    assert {'Z', 'B'}.issubset(set(validation.material_class))
    assert 'C' not in set(validation.material_class)


def test_refactored_protocol_executor_and_reactor_stage_contracts():
    from pipeline.screening.gpu_executor import execution_layout
    from pipeline.screening.protocols import ORR_PROTOCOL, PYROLYSIS_PROTOCOL
    from pipeline.stages.reactor import simulate_candidate

    assert execution_layout(3, 2, 'batched') == (3, 2)
    assert execution_layout(3, 2, 'legacy') == (6, 1)
    assert PYROLYSIS_PROTOCOL.clean.fmax_eV_A == 0.08
    assert PYROLYSIS_PROTOCOL.adsorbate.steps == 100
    assert ORR_PROTOCOL.reference.steps == 200
    row = {'candidate_id': 'cid', 'E_act': 0.72, 'dE_H': -0.3,
           'dE_CH3': -0.5, 'dE_C': -1.0,
           'screening_protocol': PYROLYSIS_PROTOCOL.protocol_id}
    fake_sweep = [
        {'CH4_conversion': 0.1, 'mock': False},
        {'CH4_conversion': 0.4, 'mock': False},
    ]
    with patch('pipeline.process.reactor_mechanisms.write_full_mechanism',
               return_value=Path('contract.yaml')) as write, \
         patch('pipeline.process.reactor_models.run_reactor_sweep',
               return_value=fake_sweep) as sweep:
        result = simulate_candidate(
            row, 'contract', [800.0, 900.0], ['PFR'], forbid_mock=True)
    assert result['candidate_id'] == 'cid'
    assert result['best_condition']['CH4_conversion'] == 0.4
    assert write.call_args.kwargs['kinetics'].methane_activation_eV == 0.72
    assert sweep.call_args.kwargs['reactor_types'] == ['PFR']


def test_arrhenius_matches_joule_and_ev_forms():
    from pipeline.common.utils import (
        R_gas, eV_to_J, arrhenius_rate, k_B_eV, tst_prefactor)

    prefactor, barrier_ev, temperature = 2.5e13, 0.83, 973.15
    expected_ev = prefactor * math.exp(-barrier_ev / (k_B_eV * temperature))
    barrier_j_mol = barrier_ev * eV_to_J * 6.02214076e23
    expected_molar = prefactor * math.exp(
        -barrier_j_mol / (R_gas * temperature))
    observed = arrhenius_rate(prefactor, barrier_ev, temperature)
    assert math.isclose(observed, expected_ev, rel_tol=2e-15)
    assert math.isclose(observed, expected_molar, rel_tol=2e-9)
    assert arrhenius_rate(prefactor, 0.0, temperature) == prefactor
    assert arrhenius_rate(prefactor, barrier_ev, temperature + 100) > observed
    assert arrhenius_rate(prefactor, barrier_ev + 0.1, temperature) < observed
    assert math.isclose(tst_prefactor(298.15), 6.212e12, rel_tol=2e-3)


def test_bep_is_linear_inside_bounds_and_explicitly_censored_outside():
    from pipeline.common.utils import bep_activation_energy

    # Default relation: Ea = 0.87 + 0.75*dE, until the declared screen bounds.
    for reaction_energy in (-0.5, 0.0, 1.0):
        expected = 0.87 + 0.75 * reaction_energy
        assert math.isclose(
            bep_activation_energy(reaction_energy), expected, rel_tol=1e-14)
    assert bep_activation_energy(-100.0) == 0.01
    assert bep_activation_energy(100.0) == 5.0


def _independent_orr(dg_oh, dg_o, dg_ooh):
    steps = (dg_ooh - 4.92, dg_o - dg_ooh, dg_oh - dg_o, -dg_oh)
    limiting_potential = -max(steps)
    return max(1.229 - limiting_potential, 0.0), int(np.argmax(steps)), steps


def test_che_stoichiometry_limiting_potential_and_nernst_terms():
    from pipeline.common.utils import orr_overpotential
    from pipeline.validation.orr_workflows import ORRCorrections, apply_orr_corrections

    names = ('step_1_OOH', 'step_2_O', 'step_3_OH', 'step_4_H2O')
    for values in ((1.23, 2.46, 3.69), (0.8, 1.7, 3.4), (1.5, 2.2, 4.1)):
        expected_eta, expected_step, steps = _independent_orr(*values)
        eta, step = orr_overpotential(*values)
        assert math.isclose(sum(steps), -4.92, abs_tol=1e-14)
        assert math.isclose(eta, expected_eta, abs_tol=1e-14)
        assert step == names[expected_step]

    base = ORRCorrections(source_id='test:independent', temperature_K=298.15)
    acid = apply_orr_corrections(1.0, 2.0, 3.0, base)
    shifted = apply_orr_corrections(
        1.0, 2.0, 3.0,
        ORRCorrections(source_id='test:independent', temperature_K=298.15,
                       electrode_potential_V=0.2, pH=1.0))
    nernst = 8.617333262e-5 * 298.15 * math.log(10.0)
    # Adsorption descriptors remain referenced at U=0. Every elementary ORR
    # step transfers one proton/electron pair and receives the same CHE shift.
    assert acid['dG_OH_eV'] == shifted['dG_OH_eV']
    assert acid['dG_O_eV'] == shifted['dG_O_eV']
    for name in acid['orr_steps_at_condition_eV']:
        difference = (shifted['orr_steps_at_condition_eV'][name] -
                      acid['orr_steps_at_condition_eV'][name])
        assert math.isclose(difference, 0.2 + nernst, rel_tol=1e-12)


def test_orr_site_coverage_ensemble_is_complete_only_with_all_cases():
    from ase import Atoms
    from pipeline.validation.orr_workflows import (
        ORRCorrections, build_orr_validation_plan, evaluate_orr_ensemble)

    slab = Atoms('Pt3', positions=[[0, 0, 0], [2.7, 0, 0], [1.35, 2.3, 0]],
                 cell=[8, 8, 12], pbc=[True, True, False])
    plan = build_orr_validation_plan(slab, coverages=(0.25, 0.5))
    assert {task['adsorbate'] for task in plan} == {'OH', 'O', 'OOH'}
    assert {task['coverage_ML'] for task in plan} == {0.25, 0.5}
    rows = [
        {'site_id': 'atop_0000', 'coverage_ML': 0.25,
         'adsorbate': adsorbate, 'dG_eV': value, 'converged': True}
        for adsorbate, value in [('OH', 0.9), ('O', 1.8), ('OOH', 3.5)]
    ]
    corrections = [
        ORRCorrections(source_id='solvation:model-a'),
        ORRCorrections(solvation_OH_eV=-0.2, solvation_OOH_eV=-0.25,
                       source_id='solvation:model-b'),
    ]
    complete = evaluate_orr_ensemble(rows, corrections, expected_cases=2)
    assert complete['complete'] is True
    assert complete['evidence_level'] == 'corrected_DFT_ensemble'
    assert complete['uncertainty']['range_V'] >= 0
    incomplete = evaluate_orr_ensemble(rows[:-1], corrections, expected_cases=2)
    assert incomplete['complete'] is False
    assert incomplete['orr_overpotential_V'] is None


def test_qe_force_parser_and_frequency_campaign_contracts():
    from ase import Atoms
    from pipeline.validation.qe_workflows import parse_atomic_forces
    from pipeline.validation.production_workflow import prepare_frequency_jobs

    output = """
     Forces acting on atoms (cartesian axes, Ry/au):
     atom    1 type  1   force =     0.01000000   -0.02000000    0.03000000
     atom    2 type  1   force =    -0.01000000    0.02000000   -0.03000000
     Total force =     0.074833
     convergence has been achieved
     JOB DONE.
    """
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / 'force.out'
        path.write_text(output)
        forces = parse_atomic_forces(path, expected_atoms=2)
        transition_state = Atoms(
            'H2', positions=[[0, 0, 4], [0, 0, 4.8]],
            cell=[8, 8, 8], pbc=True)
        prepared = prepare_frequency_jobs(
            Path(tmp) / 'frequency', transition_state, 'contract',
            active_indices=[1])
        manifest = __import__('json').loads(
            (Path(tmp) / 'frequency/frequency_forces/manifest.json').read_text())
    conversion = 13.605693122994 / 0.529177210903
    assert forces.shape == (2, 3)
    assert math.isclose(forces[0, 0], 0.01 * conversion, rel_tol=1e-12)
    assert len(manifest['jobs']) == 6
    assert manifest['active_indices'] == [1]
    assert len(manifest['manifest_sha256']) == 64
    assert prepared['status'] == 'force_jobs_pending'


def test_qe_relax_and_neb_preserve_fixed_slab_atoms():
    from ase import Atoms
    from ase.constraints import FixAtoms
    from pipeline.validation.qe_workflows import (
        write_qe_neb_input, write_qe_relax_input)

    atoms = Atoms('Ni2H', positions=[[0, 0, 1], [1, 0, 2], [1, 0, 3]],
                  cell=[6, 6, 10], pbc=True)
    atoms.set_constraint(FixAtoms(indices=[0]))
    with tempfile.TemporaryDirectory() as tmp:
        relax = Path(tmp) / 'relax.in'
        neb = Path(tmp) / 'neb.in'
        write_qe_relax_input(atoms, relax, 'constraint')
        images = [atoms.copy() for _ in range(5)]
        write_qe_neb_input(images, neb, 'constraint')
        relax_text, neb_text = relax.read_text(), neb.read_text()
    assert 'Ni 0.000000000000 0.000000000000 1.000000000000 0 0 0' in relax_text
    assert 'Ni 1.000000000000 0.000000000000 2.000000000000 1 1 1' in relax_text
    assert neb_text.count(
        'Ni 0.000000000000 0.000000000000 1.000000000000 0 0 0') == 5


def test_pyrolysis_manifest_preserves_unresolved_steps_and_validated_identity():
    import json
    from pipeline.process.reactor_mechanisms import CandidateKinetics
    from pipeline.validation.production_workflow import (
        PYROLYSIS_ELEMENTARY_STEPS, pyrolysis_campaign_status)

    with tempfile.TemporaryDirectory() as tmp:
        manifest = Path(tmp) / 'campaign.json'
        manifest.write_text(json.dumps({
            'candidate_id': 'candidate-123',
            'campaign_dir': 'calculations',
            'steps': {},
        }))
        status = pyrolysis_campaign_status(manifest)
    assert status['complete'] is False
    assert set(status['unresolved_kinetics_fields']) == set(
        PYROLYSIS_ELEMENTARY_STEPS.values())

    row = {'E_act': 0.9, 'screening_protocol': 'contract'}
    resolved = {field: 0.5 + index * 0.1 for index, field in enumerate(
        PYROLYSIS_ELEMENTARY_STEPS.values())}
    validation = {
        'candidate_id': 'candidate-123', 'complete': True,
        'evidence_level': 'converged_dft_neb_frequency',
        'resolved_kinetics_eV': resolved,
    }
    kinetics = CandidateKinetics.from_screening_row(
        row, candidate_id='candidate-123', validation=validation)
    values = kinetics.resolved()
    assert values['quantitative_status'] == 'candidate_specific'
    assert values['methane_activation_eV'] == resolved['methane_activation_eV']
    bad = dict(validation, candidate_id='different-candidate')
    try:
        CandidateKinetics.from_screening_row(
            row, candidate_id='candidate-123', validation=bad)
    except ValueError:
        pass
    else:
        raise AssertionError('mismatched validation identity was accepted')


def test_qe_inputs_use_verified_cutoffs_references_and_parallel_contracts():
    from pipeline.validation.dft_validator import (
        generate_molecule_input, generate_slab_scf_input)
    from pipeline.validation.qe_workflows import (
        QEExecutionConfig, build_qe_command, verify_sssp)

    verified = verify_sssp(['H', 'C', 'Fe', 'O'])
    assert verified['valid'] and verified['ecutwfc_Ry'] > 0
    slab = generate_slab_scf_input(
        ['Fe', 'C'], [(0, 0, 0), (1, 1, 1)],
        [[10, 0, 0], [0, 10, 0], [0, 0, 18]], ecutwfc=1,
        kpoints=(2, 2, 1))
    assert f"ecutwfc = {verify_sssp(['Fe', 'C'])['ecutwfc_Ry']}" in slab
    assert "nspin = 2" in slab and "2 2 1  1 1 0" in slab
    molecule = generate_molecule_input(
        ['H', 'H'], [(7.5, 7.5, 7.13), (7.5, 7.5, 7.87)],
        15.0, 'h2_contract', calculation='scf')
    assert "occupations = 'fixed'" in molecule
    assert 'smearing' not in molecule and 'K_POINTS {gamma}' in molecule

    with tempfile.TemporaryDirectory() as tmp:
        input_path = Path(tmp) / 'pw.in'
        input_path.write_text(molecule)
        pw = Path(tmp) / 'pw.x'
        mpi = Path(tmp) / 'mpirun'
        pw.write_text('#!/bin/sh\nexit 0\n')
        mpi.write_text('#!/bin/sh\nexit 0\n')
        pw.chmod(0o755)
        mpi.chmod(0o755)
        with patch.dict(os.environ, {'MPIEXEC': str(mpi)}):
            command = build_qe_command(
                str(pw), str(input_path),
                QEExecutionConfig(mpi_ranks=4, omp_threads=1, kpoint_pools=2))
        assert '-np' in command and command[command.index('-np') + 1] == '4'
        assert command[command.index('-nk') + 1] == '2'


def test_external_executable_resolution_is_machine_portable():
    from pipeline.common.executables import resolve_executable

    with tempfile.TemporaryDirectory() as tmp:
        executable = Path(tmp) / 'portable-tool'
        executable.write_text('#!/bin/sh\nexit 0\n')
        executable.chmod(0o755)
        with patch.dict(os.environ, {'PORTABLE_TOOL': str(executable)}):
            assert resolve_executable(
                'missing-name', env_var='PORTABLE_TOOL') == str(executable.resolve())
        with patch.dict(os.environ, {'PATH': tmp}, clear=False):
            assert resolve_executable('portable-tool') == str(executable.resolve())
        with patch.dict(os.environ, {'PORTABLE_TOOL': str(Path(tmp) / 'absent')}):
            try:
                resolve_executable('portable-tool', env_var='PORTABLE_TOOL')
            except RuntimeError as exc:
                assert 'PORTABLE_TOOL is set' in str(exc)
            else:
                raise AssertionError('invalid executable override was accepted')
        with patch.dict(os.environ, {}, clear=True), \
             patch('pipeline.common.executables.shutil.which', return_value=None):
            try:
                resolve_executable(
                    'absent-tool', env_var='ABSENT_TOOL', conda_env='tool-env')
            except RuntimeError as exc:
                message = str(exc)
                assert 'ABSENT_TOOL' in message and 'README.md' in message
            else:
                raise AssertionError('missing executable did not fail explicitly')


def test_qe_relaxation_cannot_pass_on_electronic_convergence_alone():
    from pipeline.validation.dft_validator import parse_convergence

    with tempfile.TemporaryDirectory() as tmp:
        output = Path(tmp) / 'relax.out'
        output.write_text(
            'convergence has been achieved in 8 iterations\n'
            '! total energy = -10.0 Ry\nJOB DONE\n')
        # A relaxation requires ionic/force convergence, not merely one SCF.
        assert not parse_convergence(str(output), require_ionic=True)


def test_partial_hessian_recovers_one_transition_mode():
    from pipeline.validation.qe_workflows import partial_hessian

    # One atom, diagonal Hessian [-1, 2, 3] eV/A^2. Central force differences
    # obey F(+d)-F(-d)=-2*d*H for each displacement coordinate.
    displacement = 0.01
    hessian = np.diag([-1.0, 2.0, 3.0])
    difference = -2 * displacement * hessian.T
    plus = (difference / 2).reshape(3, 1, 3)
    minus = (-difference / 2).reshape(3, 1, 3)
    result = partial_hessian(plus, minus, displacement, np.array([1.0]))
    assert result['imaginary_count'] == 1
    assert result['valid_transition_state']
    assert sum(x < 0 for x in result['frequencies_cm1']) == 1


def _pauli_matrix(pauli_string):
    matrices = {
        'I': np.eye(2),
        'X': np.array([[0, 1], [1, 0]], complex),
        'Y': np.array([[0, -1j], [1j, 0]], complex),
        'Z': np.diag([1, -1]),
    }
    result = np.array([[1.0]], complex)
    for symbol in pauli_string:
        result = np.kron(result, matrices[symbol])
    return result


def test_vqe_toy_hamiltonians_are_hermitian_and_fail_closed_as_evidence():
    from pipeline.validation.vqe_transition_state import (
        _mock_vqe_result, build_ch_splitting_hamiltonian,
        build_orr_hamiltonian)

    for terms in (build_ch_splitting_hamiltonian(), build_orr_hamiltonian()):
        assert terms and all(len(pauli) == 4 for _, pauli in terms)
        matrix = sum(float(c) * _pauli_matrix(p) for c, p in terms)
        assert np.allclose(matrix, matrix.conj().T)
        exact_ground = float(np.linalg.eigvalsh(matrix)[0])
        mock = _mock_vqe_result(terms)
        assert mock['mock'] and mock['evidence_level'] == 'mock'
        assert not mock['catalyst_specific_hamiltonian'] and not mock['benchmarked']
        # A claimed variational energy may not lie below the exact eigenvalue.
        assert mock['energy_Ha'] >= exact_ground - 1e-10


def test_fairchem_device_affinity_is_forwarded_not_merely_recorded():
    from pipeline.screening.surface_calculator import get_ocp_calculator

    with patch(
        'fairchem.core.calculate.pretrained_mlip.get_predict_unit'
    ) as get_unit, patch('fairchem.core.FAIRChemCalculator'):
        get_unit.return_value = object()
        get_ocp_calculator('esen-sm-conserving-all-oc25', device='cuda:2')
        assert get_unit.call_args.kwargs.get('device') == 'cuda:2'


def test_ranker_throughput_and_uncertainty_do_not_change_mean():
    import pandas as pd
    from pipeline.search.indexed_space import deterministic_tree_probes
    from pipeline.screening.small_data_ranker import fit_tree_ranker

    training = deterministic_tree_probes(28)
    frame = pd.DataFrame({
        'genome': [repr(x) for x in training], 'valid': True,
        'E_act': np.linspace(0.1, 2.0, len(training)),
    })
    ranker = fit_tree_ranker(frame, 'turquoise_hydrogen')
    candidates = deterministic_tree_probes(8192)
    started = time.perf_counter()
    fast, omitted = ranker.predict(candidates, uncertainty=False)
    elapsed = time.perf_counter() - started
    full, uncertainty = ranker.predict(candidates, uncertainty=True)
    assert np.allclose(fast, full, rtol=1e-12, atol=1e-12)
    assert np.all(omitted == 0) and np.all(uncertainty >= 0)
    assert len(candidates) / elapsed > 10_000


def test_esen_energy_force_and_invariance_contracts():
    if os.environ.get('RUN_ESEN_CONTRACTS') != '1':
        return
    from ase import Atoms
    from pipeline.screening.surface_calculator import get_ocp_calculator

    # The contract process is masked to one physical GPU; use its logical name.
    calc = get_ocp_calculator(device='cuda')
    assert calc is not None
    atoms = Atoms('H2', positions=[[0, 0, 0], [0, 0, 0.75]],
                  cell=[10, 10, 10], pbc=True)
    atoms.calc = calc
    energy = float(atoms.get_potential_energy())
    forces = np.asarray(atoms.get_forces())
    translated = atoms.copy()
    translated.positions += [1.2, 2.3, 3.4]
    translated.calc = calc
    assert math.isclose(
        energy, float(translated.get_potential_energy()), rel_tol=1e-5,
        abs_tol=1e-5)
    assert np.allclose(forces, translated.get_forces(), rtol=1e-4, atol=1e-4)
    assert np.linalg.norm(forces.sum(axis=0)) < 1e-3

    step = 1e-3
    plus, minus = atoms.copy(), atoms.copy()
    plus.positions[1, 2] += step
    minus.positions[1, 2] -= step
    plus.calc = calc
    minus.calc = calc
    numerical_force = -(
        plus.get_potential_energy() - minus.get_potential_energy()) / (2 * step)
    assert math.isclose(forces[1, 2], numerical_force, rel_tol=2e-2, abs_tol=2e-2)


def test_esen_dynamic_batch_matches_single_system_inference():
    if os.environ.get('RUN_ESEN_CONTRACTS') != '1':
        return
    from concurrent.futures import ThreadPoolExecutor
    from ase.build import molecule
    from pipeline.screening.batched_calculator import BatchedInferenceService
    from pipeline.screening.surface_calculator import get_ocp_calculator

    base = get_ocp_calculator(device='cuda')
    assert base is not None
    structures = [molecule('H2'), molecule('H2O')]
    for atoms in structures:
        atoms.set_cell([10, 10, 10])
        atoms.center()
        atoms.pbc = True
    expected = []
    for atoms in structures:
        atoms.calc = base
        expected.append((atoms.get_potential_energy(), atoms.get_forces().copy()))

    service = BatchedInferenceService(base, batch_wait_ms=20)
    def evaluate(atoms):
        atoms = atoms.copy()
        atoms.calc = service.calculator_proxy()
        return atoms.get_potential_energy(), atoms.get_forces()
    try:
        with ThreadPoolExecutor(max_workers=2) as executor:
            observed = list(executor.map(evaluate, structures))
    finally:
        service.close()
    for reference, batched in zip(expected, observed):
        assert math.isclose(reference[0], batched[0], abs_tol=2e-5)
        assert np.allclose(reference[1], batched[1], atol=2e-5, rtol=2e-5)


TESTS = [value for name, value in sorted(globals().items())
         if name.startswith('test_') and callable(value)]


if __name__ == '__main__':
    failures = []
    for test in TESTS:
        try:
            test()
            print(f'PASS {test.__name__}')
        except Exception as exc:
            failures.append((test.__name__, exc))
            print(f'FAIL {test.__name__}: {type(exc).__name__}: {exc}')
    print(f'\n{len(TESTS)-len(failures)}/{len(TESTS)} scientific contracts passed')
    raise SystemExit(1 if failures else 0)
