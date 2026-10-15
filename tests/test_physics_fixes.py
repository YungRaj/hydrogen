"""Focused regression tests for the reactor, fuel-cell, DFT and pipeline fixes.

No GPU, no QE/pw.x, no network. Tests needing Cantera, PySCF or OpenFermion
skip cleanly when those packages are absent.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

np = pytest.importorskip('numpy')


# ═══════════════════════════════════════════════════════════════════════════
# #1  Equilibrium X_eq is interpolated, not snapped to the nearest table point
# ═══════════════════════════════════════════════════════════════════════════


def test_x_eq_interpolates_between_table_points():
    from pipeline.reactors.equilibrium import TABULATED_X_CH4_1BAR, interpolate_x_eq

    for temperature, value in TABULATED_X_CH4_1BAR.items():
        assert interpolate_x_eq(temperature) == pytest.approx(value)
    # Ni judge band: 923.15 K used to snap to 900 K (0.684) and 973.15 K to
    # 1000 K (0.842).
    x900, x1000 = TABULATED_X_CH4_1BAR[900.0], TABULATED_X_CH4_1BAR[1000.0]
    assert interpolate_x_eq(923.15) == pytest.approx(
        x900 + (x1000 - x900) * 0.2315)
    assert interpolate_x_eq(973.15) == pytest.approx(
        x900 + (x1000 - x900) * 0.7315)
    assert x900 < interpolate_x_eq(923.15) < interpolate_x_eq(973.15) < x1000


@pytest.mark.parametrize('temperature', [700.0, 773.0, 1400.5, 2000.0, float('nan')])
def test_x_eq_rejects_temperatures_outside_table(temperature):
    from pipeline.reactors.equilibrium import interpolate_x_eq

    with pytest.raises(ValueError):
        interpolate_x_eq(temperature)


def test_reactor_core_x_eq_and_overshoot_flag_use_interpolation():
    from pipeline.reactors.equilibrium import interpolate_x_eq
    from pipeline.reactors.reactor_core import (
        ReactorConfig, _off_site_carbon_metrics, _tabulated_x_eq)

    assert _tabulated_x_eq(923.15) == pytest.approx(interpolate_x_eq(923.15))

    def flag(conversion, temperature):
        return _off_site_carbon_metrics(
            ReactorConfig(T_inlet_K=temperature), None, None,
            ch4_initial=0.95, ar_initial=0.05, pass_conversion=conversion,
            n_sites_mol=1.0, n_ch4_fed_mol=1.0, n_parcel_in_mol=1.0,
            pass_time_s=1.0, cov_start=None)['exceeds_equilibrium']

    # 0.70 is below X_eq(923.15 K) ≈ 0.721; nearest-point lookup (0.684)
    # wrongly flagged it as an equilibrium overshoot.
    assert flag(0.70, 923.15) is False
    assert flag(0.73, 923.15) is True
    # At 973.15 K the old lookup (0.842) hid real overshoots above ≈ 0.80.
    assert flag(0.82, 973.15) is True


# ═══════════════════════════════════════════════════════════════════════════
# #2  Bulk alloy QE input is a real 4-site fcc conventional cube
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def fake_sssp(monkeypatch):
    import pipeline.validation.qe_workflows as qe

    def verify(elements):
        return {
            'valid': True, 'errors': [], 'ecutwfc_Ry': 0.0, 'ecutrho_Ry': 0.0,
            'records': {element: {'filename': f'{element}.upf'}
                        for element in elements},
        }

    monkeypatch.setattr(qe, 'verify_sssp', verify)
    monkeypatch.setattr(qe, 'SSSP_DIR', '/portable/sssp')


def _parse_bulk_input(text):
    ibrav = int(text.split('ibrav = ')[1].split()[0])
    celldm = float(text.split('celldm(1) = ')[1].split()[0])
    block = text.split('ATOMIC_POSITIONS {crystal}\n')[1].split('K_POINTS')[0]
    atoms = []
    for line in block.strip().splitlines():
        symbol, *coords = line.split()
        atoms.append((symbol, tuple(float(value) for value in coords)))
    return ibrav, celldm, atoms


def _min_cubic_distance(atoms, a):
    best = math.inf
    for i in range(len(atoms)):
        for j in range(i + 1, len(atoms)):
            delta = [x - y for x, y in zip(atoms[i][1], atoms[j][1])]
            delta = [d - round(d) for d in delta]
            best = min(best, a * math.sqrt(sum(d * d for d in delta)))
    return best


def test_solid_catalyst_bulk_cell_has_physical_nearest_neighbour(fake_sssp):
    from pipeline.validation.dft_validator import _generate_alloy_input

    text = _generate_alloy_input(
        ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0), 'ni_bulk')
    ibrav, celldm, atoms = _parse_bulk_input(text)
    assert ibrav == 1  # simple cubic conventional cell, not ibrav=2 primitive
    assert len(atoms) == 4 and f'nat = {len(atoms)}' in text
    bohr_to_angstrom = 0.529177
    nearest = _min_cubic_distance(atoms, celldm) * bohr_to_angstrom
    # fcc Ni nearest neighbour a/sqrt(2) = 2.49 Å (was 1.24 Å with ibrav=2).
    assert nearest == pytest.approx(2.49, abs=0.02)


def test_molten_metal_cell_uses_host_lattice_and_four_sites(fake_sssp):
    from pipeline.validation.dft_validator import (
        _FCC_EQUIVALENT_LATTICE_BOHR, _generate_alloy_input)

    text = _generate_alloy_input(('MoltenMetal', 'Bi', 'Ni', 10.0, 1000), 'nibi')
    ibrav, celldm, atoms = _parse_bulk_input(text)
    assert ibrav == 1
    assert sorted(symbol for symbol, _ in atoms) == ['Bi', 'Bi', 'Bi', 'Ni']
    # The Bi host sets the cell, not the Ni promoter on one site.
    assert celldm == pytest.approx(_FCC_EQUIVALENT_LATTICE_BOHR['Bi'])


def test_fcc_equivalent_lattice_constants_match_atomic_volumes():
    from pipeline.validation.dft_validator import _FCC_EQUIVALENT_LATTICE_BOHR

    # bcc Fe (a=2.8665 Å, 2 atoms/cell) has the same atomic volume as an fcc
    # cell of a=3.61 Å = 6.82 Bohr; the old 5.42 Bohr was the bcc a itself.
    volume = 2.8665 ** 3 / 2
    expected = (4 * volume) ** (1 / 3) / 0.529177
    assert _FCC_EQUIVALENT_LATTICE_BOHR['Fe'] == pytest.approx(expected, abs=0.02)
    assert all(5.0 < value < 11.0 for value in _FCC_EQUIVALENT_LATTICE_BOHR.values())


# ═══════════════════════════════════════════════════════════════════════════
# #3  Two-electron integrals are reordered to OpenFermion's convention
# ═══════════════════════════════════════════════════════════════════════════


def test_chemist_to_openfermion_order_index_mapping():
    from pipeline.validation.candidate_hamiltonian import (
        _chemist_to_openfermion_order)

    rng = np.random.default_rng(7)
    eri = rng.normal(size=(3, 3, 3, 3))
    reordered = _chemist_to_openfermion_order(eri)
    for p, q, r, s in np.ndindex(eri.shape):
        assert reordered[p, q, r, s] == eri[p, s, q, r]
    assert reordered.flags['C_CONTIGUOUS']


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _diagonal_energy(pauli_terms, occupied):
    """Expectation of the diagonal Pauli terms in a JW basis determinant."""
    energy = 0.0
    for coefficient, word in pauli_terms:
        if any(symbol in 'XY' for symbol in word):
            continue
        sign = 1.0
        for qubit, symbol in enumerate(word):
            if symbol == 'Z' and qubit in occupied:
                sign = -sign
        energy += coefficient * sign
    return energy


def test_candidate_hamiltonian_keeps_inter_orbital_coulomb(tmp_path):
    pytest.importorskip('pyscf.tools.fcidump')
    pytest.importorskip('openfermion')
    from pyscf.tools import fcidump
    from pipeline.validation.candidate_hamiltonian import build_candidate_hamiltonian

    geometry = tmp_path / 'candidate.xyz'
    geometry.write_text('1\nfixture\nH 0 0 0\n')
    source = tmp_path / 'FCIDUMP'
    h1 = np.diag([-1.0, -0.4])
    eri = np.zeros((2, 2, 2, 2))
    eri[0, 0, 0, 0], eri[1, 1, 1, 1] = 0.7, 0.5
    eri[0, 0, 1, 1] = eri[1, 1, 0, 0] = 0.3  # J01 = (00|11), no exchange
    nuc = 0.2
    fcidump.from_integrals(str(source), h1, eri, 2, 2, nuc=nuc, ms=0)
    source.with_suffix('.json').write_text(json.dumps({
        'schema_version': 1, 'candidate_id': 'cand', 'geometry_id': 'cand/v1',
        'geometry_path': geometry.name, 'geometry_sha256': _sha(geometry),
        'electronic_structure_protocol': 'fixture', 'basis': 'fixture',
        'active_orbitals': 2, 'active_electrons': 2, 'charge': 0,
        'multiplicity': 1, 'frozen_orbitals': [0], 'integral_source': 'fixture',
        'fermion_to_qubit_mapping': 'Jordan-Wigner', 'fcidump_sha256': _sha(source),
    }))
    terms = build_candidate_hamiltonian(source)['pauli_terms']
    # OpenFermion spin-orbital order: 2p = alpha(p), 2p+1 = beta(p).
    # |0α 1α>: h00 + h11 + J01 = -1.1 (+ nuc). The chemist-order bug dropped
    # J01 here (the (pq|rs) coefficient multiplied a+_p a+_p a_q a_q = 0).
    assert _diagonal_energy(terms, {0, 2}) == pytest.approx(-1.0 - 0.4 + 0.3 + nuc)
    assert _diagonal_energy(terms, {0, 3}) == pytest.approx(-1.0 - 0.4 + 0.3 + nuc)
    # Doubly occupied orbital 0: 2 h00 + J00.
    assert _diagonal_energy(terms, {0, 1}) == pytest.approx(-2.0 + 0.7 + nuc)


# ═══════════════════════════════════════════════════════════════════════════
# #4  Restart reads the evidence table branch discovery actually writes
# ═══════════════════════════════════════════════════════════════════════════


def test_restart_evidence_file_matches_branch_discovery_writer():
    from pipeline.stages.candidate_io import SCREENING_EVIDENCE_DB

    source = (ROOT / 'pipeline' / 'screening' / 'genetic_optimizer.py').read_text()
    assert f"evidence_db='{SCREENING_EVIDENCE_DB}'" in source
    orchestrator = (ROOT / 'pipeline' / 'orchestrator.py').read_text()
    assert 'ga_full_database.csv' not in orchestrator


def test_restored_candidates_carry_content_addressed_identity(tmp_path):
    pd = pytest.importorskip('pandas')
    from pipeline.search.discovery import candidate_id
    from pipeline.stages.candidate_io import load_selected_candidates

    genome = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    path = tmp_path / 'branch_ranker_evidence.csv'
    pd.DataFrame([{
        'genome': str(genome), 'material_class': 'SolidCatalyst', 'E_act': 0.9,
        'E_act_censored': False, 'pyrolysis_viable': True, 'valid': True,
    }]).to_csv(path, index=False)
    restored = load_selected_candidates(path, top_k_reactor=5, top_k_dft=5)
    assert list(restored['screening_database']['candidate_id']) == [
        candidate_id(genome)]
    assert load_selected_candidates(tmp_path / 'missing.csv',
                                    top_k_reactor=1, top_k_dft=1) is None


# ═══════════════════════════════════════════════════════════════════════════
# #5  Fuel-cell stage refuses mock cathode rows unless explicitly allowed
# ═══════════════════════════════════════════════════════════════════════════


def _fuel_cell_services(cathodes):
    from pipeline.stages.fuel_cell import FuelCellServices

    def sweep(name, eta, **kwargs):
        return [{'orr_overpotential_V': eta, 'efficiency_at_peak': 0.4,
                 'peak_power_W_cm2': 1.0, 'peak_voltage_V': 0.7,
                 'peak_current_A_cm2': 1.4}]

    return FuelCellServices(
        screen_cathodes=lambda: cathodes, sweep_membranes=sweep,
        build_stack_config=lambda **kwargs: kwargs,
        model_stack=lambda config: {'net_power_kW': 1.0})


def test_fuel_cell_stage_rejects_mock_screening_rows():
    pd = pytest.importorskip('pandas')
    from pipeline.stages.fuel_cell import run_fuel_cell_stage

    cathodes = pd.DataFrame([
        {'valid': True, 'name': 'random', 'orr_overpotential_V': 0.3, 'mock': True},
    ])
    with pytest.raises(RuntimeError, match='mock'):
        run_fuel_cell_stage(top_k_pemfc=1, stack_cells=10,
                            services=_fuel_cell_services(cathodes))
    outcome = run_fuel_cell_stage(top_k_pemfc=1, stack_cells=10,
                                  services=_fuel_cell_services(cathodes),
                                  allow_mock_inputs=True)
    assert outcome.state['n_pemfc_simulations'] == 1


def test_orchestrator_threads_mock_opt_in_to_fuel_cell_stage():
    from pipeline.orchestrator import PipelineConfig, run_pipeline
    from pipeline.stages.contracts import StageOutcome
    from pipeline.stages.orchestration import (
        MemoryStateStore, PipelineComponents, PipelineRuntime)

    seen = []

    def fuel(**kwargs):
        seen.append(kwargs)
        return StageOutcome({'n_valid': 0}, {
            'cathode_database': [], 'valid_cathodes': [], 'pemfc_results': [],
            'stack_result': {}})

    def never(*args, **kwargs):
        raise AssertionError('unselected stage executed')

    store = MemoryStateStore()
    ticks = iter(float(value) for value in range(100))
    runtime = PipelineRuntime(
        clock=lambda: next(ticks), banner=lambda value: None,
        load_state=store.load, save_state=store.save,
        select_pathway_mode=lambda mode: None)
    components = PipelineComponents(
        discovery=never, reactor_batch=never, dft=never, vqe=never,
        fuel_cell=fuel, report=never, load_candidates=never)
    for allow in (False, True):
        run_pipeline(PipelineConfig(allow_mock_inputs=allow), start_phase=5,
                     end_phase=5, runtime=runtime, components=components)
    assert [call['allow_mock_inputs'] for call in seen] == [False, True]


# ═══════════════════════════════════════════════════════════════════════════
# #6  DFT resume identity, task-queue lookup, and candidate-keyed job names
# ═══════════════════════════════════════════════════════════════════════════

_CONVERGED_OUT = (
    "convergence has been achieved in 8 iterations\n"
    "!    total energy              =     -100.00000000 Ry\n"
    "JOB DONE.\n"
)


@pytest.fixture
def dft_sandbox(tmp_path, monkeypatch):
    import pipeline.validation.dft_validator as dft

    saved = []
    monkeypatch.setattr(dft, 'DFT_DIR', tmp_path)
    monkeypatch.setattr(dft, 'save_json',
                        lambda value, name, subdir=None: saved.append((name, value)))
    monkeypatch.setattr(dft, '_generate_alloy_input',
                        lambda genome, name: f"&CONTROL\n  calculation = 'scf'\n/ {name}\n")
    return dft, tmp_path, saved


def test_dft_does_not_reuse_output_from_a_different_input(dft_sandbox):
    dft, root, _ = dft_sandbox
    name = 'job'
    genome = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    calc = root / name
    calc.mkdir()
    (calc / f'{name}.in').write_text('a different candidate input\n')
    (calc / f'{name}.out').write_text(_CONVERGED_OUT)
    result = dft.validate_catalyst(name, genome, run_dft=False)
    assert 'resumed_existing_output' not in result
    assert result.get('dft_energy_Ry') is None
    assert Path(result['discarded_stale_output']).read_text() == _CONVERGED_OUT
    assert not (calc / f'{name}.out').exists()


def test_dft_resumes_output_for_identical_input(dft_sandbox):
    dft, root, _ = dft_sandbox
    name = 'job'
    genome = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    calc = root / name
    calc.mkdir()
    (calc / f'{name}.in').write_text(dft._generate_alloy_input(genome, name))
    (calc / f'{name}.out').write_text(_CONVERGED_OUT)
    result = dft.validate_catalyst(name, genome, run_dft=False)
    assert result['resumed_existing_output'] is True
    assert result['converged'] is True
    assert result['dft_energy_Ry'] == pytest.approx(-100.0)


def test_dft_structure_failure_still_writes_result(dft_sandbox, monkeypatch):
    dft, _, saved = dft_sandbox
    screener = pytest.importorskip('pipeline.screening.surface_screener')

    def boom(genome):
        raise RuntimeError('no structure')

    monkeypatch.setattr(screener, 'generate_structure', boom)
    result = dft.validate_catalyst('perov', ('Perovskite', 'La', 'Mn', 'Sr', 0.1),
                                   run_dft=False)
    assert result['converged'] is False and 'no structure' in result['error']
    assert saved and saved[-1][0] == 'perov_dft.json'


def test_task_queue_lookup_reports_terminal_state(tmp_path):
    from pipeline.search.discovery import candidate_id
    from pipeline.validation.task_queue import ValidationTaskQueue

    queue = ValidationTaskQueue(tmp_path / 'tasks.sqlite')
    genome = ('SolidCatalyst', 'Ni', 'SiO2', 'fcc111', 0.0, (), 1, 0)
    cid = queue.enqueue('app', genome, 'dft', 'proto')
    assert cid == candidate_id(genome)
    assert queue.task('app', 'unknown', 'dft', 'proto') is None
    assert queue.task('app', cid, 'dft', 'proto')['status'] == 'queued'
    assert queue.claim('app', cid, 'dft', 'proto') is True
    assert queue.recover_stale(stale_after_s=0) == 1
    assert queue.claim('app', cid, 'dft', 'proto') is True
    queue.finish('app', cid, 'dft', 'proto', True, result_path='/r.json')
    state = queue.task('app', cid, 'dft', 'proto')
    assert state == {'status': 'converged', 'attempts': 2,
                     'result_path': '/r.json', 'error': None}
    assert queue.claim('app', cid, 'dft', 'proto') is False


# ═══════════════════════════════════════════════════════════════════════════
# #7 / #8 / #15  PEMFC transport, efficiency basis, rated point; stack basis
# ═══════════════════════════════════════════════════════════════════════════


@pytest.fixture
def no_fuel_cell_io(monkeypatch):
    import pipeline.fuel_cell.pemfc as pemfc
    import pipeline.fuel_cell.stack as stack

    monkeypatch.setattr(pemfc, 'save_json', lambda *args, **kwargs: None)
    monkeypatch.setattr(stack, 'save_json', lambda *args, **kwargs: None)
    return pemfc, stack


def test_limiting_current_pressure_and_catalyst_layer_resistance(no_fuel_cell_io):
    pemfc, _ = no_fuel_cell_io
    base = dict(catalyst_layer_oxygen_resistance_s_m=0.0)
    low = pemfc.compute_limiting_current(pemfc.PEMFCConfig(P_cathode_Pa=1.0e5, **base))
    high = pemfc.compute_limiting_current(pemfc.PEMFCConfig(P_cathode_Pa=3.0e5, **base))
    # GDL-limited j_L ∝ c_O2·D ∝ P·(1/P): pressure-independent once D ∝ 1/P.
    assert high == pytest.approx(low, rel=1e-12)
    loose = pemfc.compute_limiting_current(pemfc.PEMFCConfig())
    tight = pemfc.compute_limiting_current(
        pemfc.PEMFCConfig(catalyst_layer_oxygen_resistance_s_m=200.0))
    assert tight < loose  # the declared CL O2 resistance now matters
    # Default j_L is no longer the ~28 A/cm² pure-GDL value.
    assert 5.0 < loose < 15.0


def test_peak_power_flags_sweep_truncation(no_fuel_cell_io):
    pemfc, _ = no_fuel_cell_io
    truncated = pemfc.simulate_pemfc(pemfc.PEMFCConfig(orr_overpotential_V=0.25))
    assert truncated['peak_power_limited_by_sweep'] is True
    assert truncated['peak_current_A_cm2'] == pytest.approx(3.0)
    full = pemfc.simulate_pemfc(
        pemfc.PEMFCConfig(orr_overpotential_V=0.25, i_max_A_cm2=50.0))
    assert full['peak_power_limited_by_sweep'] is False
    assert full['peak_current_A_cm2'] < full['limiting_current_A_cm2']
    assert full['peak_power_W_cm2'] >= truncated['peak_power_W_cm2']


def test_pemfc_efficiency_is_lhv_and_rated_point_is_interpolated(no_fuel_cell_io):
    pemfc, _ = no_fuel_cell_io
    result = pemfc.simulate_pemfc(pemfc.PEMFCConfig(orr_overpotential_V=0.35))
    lhv_voltage = 241800.0 / (2 * 96485.33212)
    assert result['efficiency_basis'] == 'LHV'
    assert result['efficiency_at_peak'] == pytest.approx(
        result['peak_voltage_V'] / lhv_voltage)
    # Curve crosses 0.6 V: the rated current is where V = 0.6 V.
    j, v = np.array(result['current_density']), np.array(result['voltage'])
    assert np.interp(result['rated_current_A_cm2'], j, v) == pytest.approx(0.6, abs=1e-3)
    assert result['rated_power_W_cm2'] == pytest.approx(
        0.6 * result['rated_current_A_cm2'])


def test_rated_point_is_none_when_curve_never_reaches_0p6_v(no_fuel_cell_io):
    pemfc, _ = no_fuel_cell_io
    strong = pemfc.simulate_pemfc(pemfc.PEMFCConfig(orr_overpotential_V=0.0))
    assert min(strong['voltage']) > 0.6
    assert strong['rated_current_A_cm2'] is None
    assert strong['rated_power_W_cm2'] is None
    assert strong['efficiency_at_rated'] is None
    assert pemfc._rated_point([0.1, 0.2], [0.5, 0.4], 0.6) == (None, None, None)


def test_stack_system_efficiency_never_exceeds_gross_stack_efficiency(no_fuel_cell_io):
    _, stack = no_fuel_cell_io
    for voltage, current in ((0.65, 1.5), (0.8, 0.3), (0.5, 2.5)):
        result = stack.model_stack(stack.StackConfig(
            cell_voltage_V=voltage, current_density_A_cm2=current))
        assert result['stack_efficiency'] == pytest.approx(
            voltage / (241800.0 / (2 * 96485.33212)))
        assert result['system_efficiency'] <= result['stack_efficiency'] + 1e-12


# ═══════════════════════════════════════════════════════════════════════════
# #10  Fluidized emulsion residence time is interstitial (includes ε_mf)
# ═══════════════════════════════════════════════════════════════════════════


def test_fluidized_residence_time_includes_emulsion_voidage(tmp_path):
    pytest.importorskip('cantera')
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import ReactorConfig, simulate_fluidized_bed
    from pipeline.reactors.reactor_core import FLUIDIZED_EMULSION_VOIDAGE

    path = write_full_mechanism('fix_fluid_tau', E_act_CH4=1.0, output_dir=tmp_path)
    config = ReactorConfig(
        reactor_type='Fluidized', mechanism_file=str(path),
        catalyst_name='fix_fluid_tau', T_inlet_K=973.15, max_regen_cycles=0)
    result = simulate_fluidized_bed(config)
    delta = result['bubble_fraction']
    expected = (config.bed_height_m * (1 - delta) * FLUIDIZED_EMULSION_VOIDAGE
                / config.u_mf_m_s)
    assert result['residence_time_s'] == pytest.approx(expected)


# ═══════════════════════════════════════════════════════════════════════════
# VQE stage: candidate identities only, never placeholder champion_<i>
# ═══════════════════════════════════════════════════════════════════════════


def test_vqe_stage_never_uses_placeholder_names(tmp_path):
    from pipeline.stages.vqe import run_vqe_stage

    seen = []
    (tmp_path / 'cand_b').mkdir()
    (tmp_path / 'cand_b' / 'FCIDUMP').write_text('')
    outcome = run_vqe_stage(
        [{'candidate_id': 'cand_a'}, {'candidate_id': 'cand_b'}],
        top_k=5, execute_quantum=False, hamiltonian_root=tmp_path,
        validator=lambda name, reaction, **kwargs: seen.append((name, kwargs)) or {})
    assert [name for name, _ in seen] == ['cand_b']
    assert seen[0][1]['candidate_hamiltonian'] == str(tmp_path / 'cand_b' / 'FCIDUMP')
    assert outcome.state['n_missing_hamiltonians'] == 1
    with pytest.raises(ValueError):
        run_vqe_stage([{'candidate_id': '../escape'}], top_k=1,
                      execute_quantum=False, hamiltonian_root=tmp_path,
                      validator=lambda *a, **k: {})


def test_b2_b3_b4_reactor_contracts(tmp_path):
    pytest.importorskip('cantera')
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import (
        ReactorConfig,
        simulate_fluidized_bed,
        simulate_mmbcr,
        simulate_pfr,
    )
    from pipeline.reactors.reactor_core import _validate_carbon_policy

    # B3: validator enforces 0 < k0 <= 1 m/s
    with pytest.raises(ValueError, match='mmbcr_interfacial_k0_m_s'):
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=0.0))
    with pytest.raises(ValueError, match='mmbcr_interfacial_k0_m_s'):
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=1.5))

    path = write_full_mechanism('b2_b3_b4_cat', E_act_CH4=0.84, output_dir=tmp_path)

    # B3: MMBCR is bounded by X_eq and kinetics-limited at Chen-scale E_a (~0.84 eV)
    mmbcr_res = simulate_mmbcr(
        ReactorConfig(
            reactor_type='MMBCR',
            mechanism_file=str(path),
            catalyst_name='b2_b3_b4_cat',
            T_inlet_K=1300.0,
            catalyst_E_act_eV=0.84,
        )
    )
    assert 0.0 < mmbcr_res['CH4_conversion'] < mmbcr_res['X_eq_table']
    assert mmbcr_res['CH4_conversion'] < 0.5 * mmbcr_res['X_eq_table']

    # B2 & B4: PFR reports both single-pass and cycle-averaged conversion,
    # pilot residence-time band (2-12 s), and endothermic heat duty.
    pfr_res = simulate_pfr(
        ReactorConfig(
            reactor_type='PFR',
            mechanism_file=str(path),
            catalyst_name='b2_b3_b4_cat',
            T_inlet_K=973.15,
            max_regen_cycles=2,
            regen_coverage_threshold=0.0,
        )
    )
    assert 'single_pass_CH4_conversion' in pfr_res
    assert 'cycle_averaged_CH4_conversion' in pfr_res
    assert pfr_res['residence_time_within_pilot_band'] is True
    assert pfr_res['endothermic_heat_duty_kJ_per_mol_CH4_fed'] == pytest.approx(
        74.85 * pfr_res['CH4_conversion']
    )

    fluid_res = simulate_fluidized_bed(
        ReactorConfig(
            reactor_type='Fluidized',
            mechanism_file=str(path),
            catalyst_name='b2_b3_b4_cat',
            T_inlet_K=973.15,
        )
    )
    assert 'single_pass_CH4_conversion' in fluid_res
    assert 'cycle_averaged_CH4_conversion' in fluid_res

