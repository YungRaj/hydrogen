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
import time
from pathlib import Path
from unittest.mock import patch

import numpy as np


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
    assert math.isclose(
        acid['dG_OH_eV'] - shifted['dG_OH_eV'], 0.2 + nernst,
        rel_tol=1e-12)
    assert math.isclose(
        acid['dG_O_eV'] - shifted['dG_O_eV'], 2 * (0.2 + nernst),
        rel_tol=1e-12)


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
        command = build_qe_command(
            '/home/ilhanraja/miniconda3/envs/qe-env/bin/pw.x', str(input_path),
            QEExecutionConfig(mpi_ranks=4, omp_threads=1, kpoint_pools=2))
        assert '-np' in command and command[command.index('-np') + 1] == '4'
        assert command[command.index('-nk') + 1] == '2'


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

    calc = get_ocp_calculator(device='cuda:0')
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
