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
            feed_initial=0.95, ar_initial=0.05, pass_conversion=conversion,
            n_sites_mol=1.0, n_feed_fed_mol=1.0, n_parcel_in_mol=1.0,
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


def test_muradov_upham_ergun_reactor_contracts(tmp_path):
    pytest.importorskip('cantera')
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import (
        ReactorConfig,
        simulate_fluidized_bed,
        simulate_mmbcr,
        simulate_pfr,
    )
    from pipeline.reactors.reactor_core import _validate_carbon_policy

    # Upham (2017) / Chen (2023): validator enforces 0 < k0 <= 1 m/s
    with pytest.raises(ValueError, match='mmbcr_interfacial_k0_m_s'):
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=0.0))
    with pytest.raises(ValueError, match='mmbcr_interfacial_k0_m_s'):
        _validate_carbon_policy(ReactorConfig(mmbcr_interfacial_k0_m_s=1.5))

    path = write_full_mechanism('muradov_upham_ergun_cat', E_act_CH4=0.84, output_dir=tmp_path)

    # Upham (2017) / Chen (2023): MMBCR is bounded by X_eq and kinetics-limited at Chen-scale E_a (~0.84 eV)
    mmbcr_res = simulate_mmbcr(
        ReactorConfig(
            reactor_type='MMBCR',
            mechanism_file=str(path),
            catalyst_name='muradov_upham_ergun_cat',
            T_inlet_K=1300.0,
            catalyst_E_act_eV=0.84,
        )
    )
    assert 0.0 < mmbcr_res['CH4_conversion'] < mmbcr_res['X_eq_table']
    assert mmbcr_res['CH4_conversion'] < 0.5 * mmbcr_res['X_eq_table']

    # Muradov (2005) & Ergun (1952): PFR reports both single-pass and cycle-averaged conversion,
    # pilot residence-time band (2-12 s), and endothermic heat duty.
    pfr_res = simulate_pfr(
        ReactorConfig(
            reactor_type='PFR',
            mechanism_file=str(path),
            catalyst_name='muradov_upham_ergun_cat',
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
            catalyst_name='muradov_upham_ergun_cat',
            T_inlet_K=973.15,
        )
    )
    assert 'single_pass_CH4_conversion' in fluid_res
    assert 'cycle_averaged_CH4_conversion' in fluid_res


def test_mlip_neb_upgrade_and_stability_hull_gate(tmp_path):
    pd = pytest.importorskip('pandas')
    from pipeline.reactors.mechanisms import CandidateKinetics
    from pipeline.screening.objectives import operating_stability_and_hull_gate
    from pipeline.screening.surface_screener import upgrade_candidates_with_mlip_neb

    df = pd.DataFrame([{
        'candidate_id': 'ni_cu_test',
        'genome': str(('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Cu',), 1, 0)),
        'material_class': 'SolidCatalyst',
        'E_act': 0.88,
        'dE_H': -0.35,
        'dE_CH3': -0.45,
        'dE_C': -0.50,
        'screening_protocol': 'esen-v1',
        'valid': True,
    }])
    upgraded = upgrade_candidates_with_mlip_neb(
        df, kinetics_validation_dir=tmp_path
    )
    val_path = tmp_path / 'ni_cu_test' / 'pyrolysis_validation.json'
    assert val_path.is_file()
    val_payload = json.loads(val_path.read_text())
    assert val_payload['complete'] is True
    assert val_payload['evidence_level'] == 'converged_mlip_neb'
    kinetics = CandidateKinetics.from_screening_row(
        upgraded.iloc[0], candidate_id='ni_cu_test', validation=val_payload
    )
    resolved = kinetics.resolved()
    assert resolved['quantitative_status'] == 'mlip_neb_complete'

    stable_gate = operating_stability_and_hull_gate(
        ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, ('Cu',), 1, 0),
        'turquoise_hydrogen',
    )
    assert stable_gate['phase_stable'] is True
    assert stable_gate['synthesizable'] is True

    unstable_gate = operating_stability_and_hull_gate(
        ('MOF', 'Zr', 'BDC', 'fcu', 'none', 0),
        'turquoise_hydrogen',
    )
    assert unstable_gate['phase_stable'] is False
    assert unstable_gate['synthesizable'] is False


def test_sota_scorecard_solvation_and_campaign_preset(tmp_path):
    import run_production_campaign as campaign
    from pipeline.reactors.scorecard import (
        NI_JUDGE_CATALYST,
        build_solids_scorecard,
    )
    from pipeline.screening.fc_screener import (
        compute_grand_canonical_orr_corrections,
        compute_pourbaix_stability_margin,
    )

    records = [
        {
            'reactor_type': 'PFR',
            'catalyst_name': NI_JUDGE_CATALYST,
            'T_K': 973.15,
            'CH4_conversion': 0.1694,
            'single_pass_CH4_conversion': 0.1694,
            'catalyst_E_act_eV': 0.90,
            'catalyst_dE_H_eV': -0.52,
            'status': 'complete',
            'surface_loaded': True,
            'exceeds_equilibrium': False,
        },
        {
            'reactor_type': 'PFR',
            'catalyst_name': 'cand_co_ce_01',
            'T_K': 1100.0,
            'CH4_conversion': 0.642,
            'single_pass_CH4_conversion': 0.642,
            'catalyst_E_act_eV': 0.62,
            'catalyst_dE_H_eV': -0.28,
            'status': 'complete',
            'surface_loaded': True,
            'exceeds_equilibrium': False,
            'kinetics_status': 'mlip_neb_complete',
        },
    ]
    card = build_solids_scorecard(
        records,
        judge_catalyst=NI_JUDGE_CATALYST,
        headline_t_min=923.15,
        headline_t_max=973.15,
    )
    assert card['headline_solids_conversion'] == pytest.approx(0.1694)
    assert card['candidate_champion_conversion'] == pytest.approx(0.642)
    assert card['candidate_champion_catalyst'] == 'cand_co_ce_01'
    assert card['solids_max_excluding_h_parked'] is None

    gc_pt = compute_grand_canonical_orr_corrections(
        0.85, 1.75, 3.85,
        genome=('SolidCatalyst', 'Pt', 'C', 'fcc111', 0.0, ('Co',), 1, 0),
    )
    assert gc_pt['pourbaix_acid_stable'] is True
    assert gc_pt['pourbaix_dissolution_margin_V'] > 0.30
    assert 0.0 < gc_pt['orr_overpotential_gc_V'] < 2.0

    pb_mn = compute_pourbaix_stability_margin(
        ('SolidCatalyst', 'Mn', 'C', 'bcc110', 0.0, (), 1, 0),
        U_V_vs_RHE=0.80,
        pH=0.0,
    )
    assert pb_mn['pourbaix_acid_stable'] is False

    parser = campaign._build_parser()
    args = parser.parse_args([
        '--sota-campaign',
        '--results-dir',
        str(tmp_path / 'sota_out'),
    ])
    campaign._resolve_paths(args)
    assert args.no_dft is True
    assert args.no_vqe is True
    assert args.branch_class_floor == 4
    assert args.branch_max_leaves == 56
    assert args.calibration_probes == 1000
    assert args.validation_batch == 500
    assert args.min_validation_per_class == 4


def test_muradov_upham_ergun_reactor_guardrails_and_diagnostics(tmp_path):
    from pipeline.reactors.reactor_core import (
        ReactorConfig,
        _validate_carbon_policy,
        _validate_reactor_config,
    )

    # Ergun (1952) & Froment (2011): PFR residence time bound and max_residence_time_s validation
    with pytest.raises(ValueError, match='max_residence_time_s must be positive'):
        _validate_reactor_config(
            ReactorConfig(reactor_type='PFR', max_residence_time_s=0.0)
        )
    with pytest.raises(ValueError, match='exceeds max_residence_time_s'):
        _validate_reactor_config(
            ReactorConfig(
                reactor_type='PFR',
                bed_length_m=5.0,
                bed_porosity=0.5,
                gas_velocity_m_s=0.01,
                max_residence_time_s=60.0,
            )
        )

    # Upham (2017) & Chen (2023): MMBCR interfacial k0 guardrail
    with pytest.raises(ValueError, match='mmbcr_interfacial_k0_m_s'):
        _validate_carbon_policy(
            ReactorConfig(reactor_type='MMBCR', mmbcr_interfacial_k0_m_s=2.5)
        )

    pytest.importorskip('cantera')
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import simulate_fluidized_bed, simulate_pfr

    mech_path = write_full_mechanism('muradov_ergun_diag_cat', E_act_CH4=0.88, output_dir=tmp_path)

    pfr_res = simulate_pfr(
        ReactorConfig(
            reactor_type='PFR',
            mechanism_file=str(mech_path),
            catalyst_name='muradov_ergun_diag_cat',
            T_inlet_K=973.15,
            bed_length_m=0.5,
            bed_porosity=0.4,
            gas_velocity_m_s=0.05,
        )
    )
    assert 'overloaded_bed_warning' in pfr_res
    assert 'plugging_risk' in pfr_res
    assert 'carbon_outfeed_policy_summary' in pfr_res
    assert 'regen_mechanism=' in pfr_res['carbon_outfeed_policy_summary']

    fb_res = simulate_fluidized_bed(
        ReactorConfig(
            reactor_type='Fluidized',
            mechanism_file=str(mech_path),
            catalyst_name='muradov_ergun_diag_cat',
            T_inlet_K=973.15,
            gas_velocity_m_s=0.08,
            u_mf_m_s=0.02,
        )
    )
    assert fb_res['circulating_outfeed_active'] is True
    assert 'carbon_outfeed_policy_summary' in fb_res
    assert 'circulating_outfeed=True' in fb_res['carbon_outfeed_policy_summary']


def test_assumption_audit_fixes(tmp_path):
    from pipeline.economics.tea import estimate_hydrogen_cost
    from pipeline.electrochemistry.native_models import (
        electrochemical_screening_point,
        turquoise_reversible_voltage_V,
    )
    from pipeline.reactors.equilibrium import (
        equilibrium_ch4_conversion_tp,
        interpolate_x_eq,
    )
    from pipeline.reactors.mechanisms import (
        off_site_carbon_allowed,
        write_full_mechanism,
        write_gas_only_mechanism,
    )
    from pipeline.reactors.pbr import DEFAULT_STICKING_COEFF, pbr_pellet_effectiveness
    from pipeline.reactors.reactor_core import (
        ReactorConfig,
        _mmbcr_hydrodynamics,
        adiabatic_temperature_drop_k,
        class_site_utilization_factor,
        kunii_levenspiel_bubble_exchange,
        wen_yu_umf_m_s,
    )
    from pipeline.screening.mlip_neb import ZPE_CORRECTION_CH_EV
    from pipeline.screening.structure_builders import generate_structure

    # 1. Exact T/P/dilution equilibrium conversion (Le Chatelier scaling)
    x_1bar = equilibrium_ch4_conversion_tp(1000.0, 100000.0, 1.0)
    x_10bar = equilibrium_ch4_conversion_tp(1000.0, 1000000.0, 1.0)
    x_diluted = equilibrium_ch4_conversion_tp(1000.0, 100000.0, 0.50)
    assert x_1bar == pytest.approx(interpolate_x_eq(1000.0), rel=1e-4)
    assert x_10bar < x_1bar
    assert x_diluted > x_1bar

    # 2. Doped/bimetallic SolidCatalyst & SupportedNanoparticle off-site carbon admission
    cu_ni_genome = ('SolidCatalyst', 'Cu', 'Al2O3', 'fcc111', 0.0, ('Ni',), 1, 0)
    cu_zn_genome = ('SolidCatalyst', 'Cu', 'Al2O3', 'fcc111', 0.0, ('Zn',), 1, 0)
    assert off_site_carbon_allowed(cu_ni_genome, 'SolidCatalyst') is True
    assert off_site_carbon_allowed(cu_zn_genome, 'SolidCatalyst') is False

    # 3. Reversible C2 gas pyrolysis and graphite adjacency on regeneration
    gas_yaml = Path(write_gas_only_mechanism()).read_text()
    assert '2 CH4 <=> C2H6 + H2' in gas_yaml
    assert 'C2H6 <=> C2H4 + H2' in gas_yaml
    assert 'C2H4 <=> C2H2 + H2' in gas_yaml
    sac_regen_yaml = Path(
        write_full_mechanism(
            'sac_regen',
            E_act_CH4=0.85,
            output_dir=tmp_path,
            off_site_carbon=False,
            carbon_regeneration_rate_1_s=0.05,
        )
    ).read_text()
    assert 'adjacent-phases: [gas, graphite]' in sac_regen_yaml

    # 4. Class-specific active site utilization and adiabatic temperature drop
    assert class_site_utilization_factor('SolidCatalyst') == pytest.approx(1.0)
    assert class_site_utilization_factor('SAC') < 0.25
    assert adiabatic_temperature_drop_k(0.95, 0.50) > 400.0

    # 5. PBR sticking parity (0.01) and linear Thiele modulus scaling with pellet diameter
    assert DEFAULT_STICKING_COEFF == pytest.approx(0.01)
    p1 = pbr_pellet_effectiveness(
        ReactorConfig(reactor_type='PBR', catalyst_particle_mm=1.0),
        1073.15,
        150000.0,
        0.25,
        0.40,
    )
    p3 = pbr_pellet_effectiveness(
        ReactorConfig(reactor_type='PBR', catalyst_particle_mm=3.0),
        1073.15,
        150000.0,
        0.25,
        0.40,
    )
    assert p3['thiele_modulus'] == pytest.approx(3.0 * p1['thiele_modulus'], rel=1e-4)

    # 6. Fluidized Wen-Yu u_mf and Kunii-Levenspiel interphase exchange
    cfg_fb = ReactorConfig(reactor_type='Fluidized', T_inlet_K=1073.15)
    assert wen_yu_umf_m_s(cfg_fb) > 0.0
    kl = kunii_levenspiel_bubble_exchange(cfg_fb, 0.25, 0.10, 0.02)
    assert 0.0 < kl['effective_bypass_fraction'] < 0.25

    # 7. MMBCR composition-dependent liquid metal density & surface tension
    from pipeline.reactors.mechanisms import CandidateKinetics

    mech_bi = write_full_mechanism(
        'bi_ni_melt',
        kinetics=CandidateKinetics(
            methane_activation_eV=0.85,
            candidate_id='bi_ni_melt',
            material_class='MoltenMetal',
            genome=('MoltenMetal', 'Bi', 'Ni', 20.0, 1273.15),
            screening_protocol='esen-v1',
        ),
        output_dir=tmp_path,
    )
    mech_ga = write_full_mechanism(
        'ga_ni_melt',
        kinetics=CandidateKinetics(
            methane_activation_eV=0.85,
            candidate_id='ga_ni_melt',
            material_class='MoltenMetal',
            genome=('MoltenMetal', 'Ga', 'Ni', 20.0, 1273.15),
            screening_protocol='esen-v1',
        ),
        output_dir=tmp_path,
    )
    hydro_bi = _mmbcr_hydrodynamics(
        ReactorConfig(reactor_type='MMBCR', mechanism_file=str(mech_bi))
    )
    hydro_ga = _mmbcr_hydrodynamics(
        ReactorConfig(reactor_type='MMBCR', mechanism_file=str(mech_ga))
    )
    assert hydro_bi['melt_density_kg_m3'] > 8500.0
    assert hydro_ga['melt_density_kg_m3'] < 6500.0
    assert hydro_ga['melt_surface_tension_N_m'] > hydro_bi['melt_surface_tension_N_m']

    # 8. MoltenMetal liquid corrugation and Gibbs promoter sub-surface offset
    assert ZPE_CORRECTION_CH_EV == pytest.approx(-0.16)
    molten_atoms, molten_surf, _ = generate_structure(
        ('MoltenMetal', 'Bi', 'Ni', 27.0, 1273.15), 0
    )
    ni_z = [
        float(molten_atoms.positions[i, 2])
        for i in molten_surf
        if molten_atoms[i].symbol == 'Ni'
    ]
    bi_z = [
        float(molten_atoms.positions[i, 2])
        for i in molten_surf
        if molten_atoms[i].symbol == 'Bi'
    ]
    assert np.mean(ni_z) < np.mean(bi_z)

    # 9. Temperature-dependent Nernst potential and PSA CH4 recycle TEA
    assert turquoise_reversible_voltage_V(300.0, 4) > turquoise_reversible_voltage_V(
        1000.0, 4
    )
    ec_pt = electrochemical_screening_point({
        'descriptor_eV': 0.0,
        'optimum_descriptor_eV': 0.0,
        'reference_current_A_cm2': 0.1,
        'temperature_K': 900.0,
        'overpotential_V': 0.15,
        'area_m2': 1e-4,
        'diffusivity_m2_s': 1e-9,
        'concentration_mol_m3': 1000.0,
        'boundary_layer_m': 1e-4,
        'gap_m': 1e-3,
        'conductivity_S_m': 10.0,
        'ch4_feed_mol_s': 1e-3,
    })
    assert 'thermodynamic_reversible_voltage_V' in ec_pt

    tea_50 = estimate_hydrogen_cost(0.50, 'base')
    assert tea_50['recycle_ratio'] == pytest.approx(1.0)
    assert tea_50['net_ch4_utilization'] > 0.90


def test_second_pass_assumption_fixes(tmp_path):
    """Verify second-pass physical, chemical, kinetic, atomistic, electrochemistry, PEMFC, and TEA fixes."""
    from pipeline.economics.tea import estimate_hydrogen_cost
    from pipeline.electrochemistry.cell import CellOperatingPoint
    from pipeline.electrochemistry.kinetics import butler_volmer_current_density
    from pipeline.electrochemistry.native_models import electrochemical_screening_point
    from pipeline.fuel_cell.pemfc import (
        PEMFCConfig,
        simulate_pemfc,
        springer_conductivity_multiplier,
        sweep_membranes,
        water_saturation_pressure_Pa,
    )
    from pipeline.reactors.carbon_accounting import (
        SOLID_CARBON_DENSITY_KG_M3,
        pyrolysis_reaction_enthalpy_kJ_mol,
    )
    from pipeline.reactors.pbr import pbr_pellet_effectiveness
    from pipeline.reactors.reactor_core import (
        ReactorConfig,
        _policy_metadata,
        adiabatic_temperature_drop_k,
    )
    from pipeline.reactors.scorecard import evaluate_supported_ni_filament_closure_gate
    from pipeline.reactors.sweeps.activation import ABLATION_RESPONSE
    from pipeline.screening.fc_cathode_screener import (
        SOLVATION_CORRECTION_EV,
        build_orr_adsorbates,
    )
    from pipeline.screening.small_data_ranker import _extract_secondary_targets
    from pipeline.screening.structure_builders import generate_structure
    from pipeline.search.scope import phase_stable_at_application_T

    # 1. Unified solid carbon density (2160 kg/m^3 = 2.16 g/cm^3) and Kirchhoff enthalpy ΔH(T)
    assert SOLID_CARBON_DENSITY_KG_M3 == pytest.approx(2160.0)
    dh_298 = pyrolysis_reaction_enthalpy_kJ_mol(298.15)
    dh_1073 = pyrolysis_reaction_enthalpy_kJ_mol(1073.15)
    assert dh_298 == pytest.approx(74.85, rel=0.02)
    assert dh_1073 > 85.0  # Kirchhoff increase at high T
    assert adiabatic_temperature_drop_k(0.95, 0.50, T_K=1073.15) > adiabatic_temperature_drop_k(
        0.95, 0.50
    )

    # 2. PBR effectiveness fouling with theta_c and supported-Ni filament closure gate PBR support
    cfg_pbr = ReactorConfig(reactor_type='PBR', catalyst_particle_mm=2.0)
    clean_eff = pbr_pellet_effectiveness(cfg_pbr, 1073.15, 150000.0, 0.25, 0.40, theta_c=0.0)
    fouled_eff = pbr_pellet_effectiveness(cfg_pbr, 1073.15, 150000.0, 0.25, 0.40, theta_c=0.75)
    assert fouled_eff['thiele_modulus'] < clean_eff['thiele_modulus']
    assert 'PBR' in ABLATION_RESPONSE
    filament_gate = evaluate_supported_ni_filament_closure_gate({})
    assert 'pbr_roi_conversion' in filament_gate and 'pbr_conversion_ok' in filament_gate

    # 3. MMBCR _policy_metadata preserves composition-resolved melt properties
    from pipeline.reactors.mechanisms import CandidateKinetics, write_full_mechanism

    mech_ga = write_full_mechanism(
        'ga_ni_policy',
        kinetics=CandidateKinetics(
            methane_activation_eV=0.85,
            candidate_id='ga_ni_policy',
            material_class='MoltenMetal',
            genome=('MoltenMetal', 'Ga', 'Ni', 20.0, 1273.15),
            screening_protocol='esen-v1',
        ),
        output_dir=tmp_path,
    )
    pol_ga = _policy_metadata(ReactorConfig(reactor_type='MMBCR', mechanism_file=str(mech_ga)))
    assert pol_ga['melt_surface_tension_N_m'] > 0.60
    assert 'Ga_Ni' in str(pol_ga['melt_property_basis'])

    # 4. Dopant/promoter active site ordering in generate_alloy_slab and low-melting solid host rejection
    saa_atoms, saa_surf, _ = generate_structure(('SAA', 'Ni', 'Cu', 'fcc111', 1), 0)
    assert saa_atoms[saa_surf[0]].symbol == 'Ni'
    assert (
        phase_stable_at_application_T(
            ('SolidCatalyst', 'Ga', 'Al2O3', 'fcc111', 0.0, (), 0, 0), T_K=1073.15
        )['status']
        == 'out_of_scope'
    )
    assert (
        phase_stable_at_application_T(
            ('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 0, 0), T_K=1073.15
        )['status']
        == 'candidate'
    )

    # 5. Small-data ranker signed desorption overbinding penalty (-dg_c)
    sec_weak_c = _extract_secondary_targets(
        {'dE_H': 0.0, 'dE_C': 3.50, 'temperature_K': 1000.0},
        'turquoise_hydrogen',
    )
    sec_overbound_c = _extract_secondary_targets(
        {'dE_H': 0.0, 'dE_C': 2.00, 'temperature_K': 1000.0},
        'turquoise_hydrogen',
    )
    assert sec_weak_c['desorption_overbinding_penalty_eV'] == pytest.approx(0.0)
    assert sec_overbound_c['desorption_overbinding_penalty_eV'] > 0.50

    # 6. FC cathode screener solvation corrections and bent OOH* geometry
    assert SOLVATION_CORRECTION_EV['OH'] == pytest.approx(-0.30)
    assert SOLVATION_CORRECTION_EV['OOH'] == pytest.approx(-0.25)
    ooh_atoms = build_orr_adsorbates(saa_atoms, saa_surf)['OOH']
    # Distal O is displaced laterally (bent O-O-H), not colinear along z
    o1_pos = ooh_atoms.positions[-3]
    o2_pos = ooh_atoms.positions[-2]
    assert abs(float(o2_pos[0] - o1_pos[0])) > 0.5

    # 7. Electrochemistry 4-electron stoichiometry and 1-electron RDS Butler-Volmer default
    assert CellOperatingPoint(1.0, 0.1, -0.1, 0.1, 1.0, 1.0, 1000.0).electrons == 4
    bv_default = butler_volmer_current_density(0.1, 0.1, 1000.0)
    bv_n1 = butler_volmer_current_density(0.1, 0.1, 1000.0, electrons=1)
    assert bv_default == pytest.approx(bv_n1)
    ec_res = electrochemical_screening_point({
        'descriptor_eV': 0.0,
        'optimum_descriptor_eV': 0.0,
        'reference_current_A_cm2': 0.1,
        'temperature_K': 900.0,
        'overpotential_V': 0.15,
        'area_m2': 1e-4,
        'diffusivity_m2_s': 1e-9,
        'concentration_mol_m3': 1000.0,
        'boundary_layer_m': 1e-4,
        'gap_m': 1e-3,
        'conductivity_S_m': 10.0,
        'ch4_feed_mol_s': 1e-3,
    })
    assert ec_res['stoichiometric_electrons'] == 4
    assert ec_res['charge_transfer_electrons'] == 1

    # 8. PEMFC Springer conductivity hydration scaling, wet Nernst OCV, and PGM-free loading
    assert springer_conductivity_multiplier(353.15, 1.0, 1.0) == pytest.approx(1.0)
    assert springer_conductivity_multiplier(353.15, 0.4, 0.4) < 0.5
    assert water_saturation_pressure_Pa(353.15) > 40000.0
    pem_wet = simulate_pemfc(PEMFCConfig(cathode_catalyst='Pt_wet', RH_anode=1.0, RH_cathode=1.0))
    pem_dry = simulate_pemfc(PEMFCConfig(cathode_catalyst='Pt_dry', RH_anode=0.3, RH_cathode=0.3))
    assert pem_dry['R_ohmic_ohm_cm2'] > pem_wet['R_ohmic_ohm_cm2']
    assert pem_wet['OCV_wet_V'] < pem_wet['OCV_V']
    sac_sweeps = sweep_membranes('Fe_N_C', 0.38, material_class='SAC')
    assert sac_sweeps[0]['pgm_loading_mg_cm2'] == pytest.approx(0.05)

    # 9. TEA carbon selectivity credit scaling and temperature-dependent duty
    tea_full_c = estimate_hydrogen_cost(0.80, 'base', solid_c_selectivity=1.0)
    tea_half_c = estimate_hydrogen_cost(0.80, 'base', solid_c_selectivity=0.5)
    assert tea_half_c['h2_cost_usd_kg'] > tea_full_c['h2_cost_usd_kg']
    tea_high_t = estimate_hydrogen_cost(0.80, 'base', temperature_K=1373.15)
    tea_low_t = estimate_hydrogen_cost(0.80, 'base', temperature_K=973.15)
    assert tea_high_t['energy_input_kwh_kg_h2'] > tea_low_t['energy_input_kwh_kg_h2']

    # 10. PBR / Ergun residence / Muradov regenerated-cycle sweep diagnostics and class-anchored DFT pyrolysis calibration
    from pipeline.evidence.pyrolysis_benchmark import (
        calibrated_interfacial_pyrolysis_score_eV,
        class_dft_barrier_offset_eV,
    )
    from pipeline.reactors.scorecard import build_solids_scorecard
    from pipeline.reactors.sweeps.carbon_criteria import (
        evaluate_pbr_and_ergun_residence_diagnostics,
    )
    from pipeline.screening.stage_selection import annotate_interfacial_pyrolysis

    pbr_rec = {
        'status': 'complete',
        'reactor_type': 'PBR',
        'catalyst_name': 'Ni_Al2O3_ref',
        'T_K': 973.15,
        'CH4_conversion': 0.42,
        'surface_loaded': True,
        'exceeds_equilibrium': False,
        'overall_effectiveness_factor': 0.72,
        'thiele_modulus': 1.4,
        'ergun_ok': True,
        'residence_time_within_pilot_band': True,
        'plugging_risk': 'low',
        'single_pass_CH4_conversion': 0.42,
        'cycle_averaged_CH4_conversion': 0.39,
    }
    diag = evaluate_pbr_and_ergun_residence_diagnostics([pbr_rec])
    assert diag['n_pbr_records'] == 1
    assert diag['mean_pbr_overall_effectiveness_factor'] == pytest.approx(0.72)
    assert diag['n_within_pilot_residence_band'] == 1
    assert diag['mean_regenerated_cycle_minus_single_pass_X'] == pytest.approx(-0.03)
    scorecard_with_diag = build_solids_scorecard(
        [pbr_rec], judge_catalyst='Ni_Al2O3_ref', headline_t_min=923.15
    )
    assert scorecard_with_diag['pbr_ergun_residence_diagnostics']['n_pbr_records'] == 1
    assert scorecard_with_diag['pbr_ergun_residence_diagnostics'][
        'mean_pbr_overall_effectiveness_factor'
    ] == pytest.approx(0.72)

    ni_anchor = class_dft_barrier_offset_eV('SolidCatalyst')
    assert ni_anchor['anchor_offset_eV'] == pytest.approx(1.05 - 0.386)
    cal_score = calibrated_interfacial_pyrolysis_score_eV({
        'material_class': 'SolidCatalyst',
        'E_act_raw': 0.386,
        'c_height_above_surface_A': 0.20,
    })
    assert cal_score['calibrated_interfacial_score_eV'] == pytest.approx(1.25)

    pd = pytest.importorskip('pandas')
    annotated_df = annotate_interfacial_pyrolysis(
        pd.DataFrame([
            {
                'material_class': 'SolidCatalyst',
                'E_act_raw': 0.386,
                'c_height_above_surface_A': -0.20,
            },
            {
                'material_class': 'MoltenMetal',
                'E_act_raw': 0.562,
                'c_height_above_surface_A': 0.15,
            },
        ])
    )
    molten_anchor = class_dft_barrier_offset_eV('MoltenMetal')
    assert annotated_df.loc[0, 'calibrated_activation_eV'] == pytest.approx(1.05)
    assert annotated_df.loc[0, 'calibrated_interfacial_score_eV'] == pytest.approx(1.25)
    assert annotated_df.loc[1, 'calibrated_activation_eV'] == pytest.approx(
        0.562 + molten_anchor['anchor_offset_eV']
    )
    assert annotated_df.loc[1, 'calibrated_interfacial_score_eV'] == pytest.approx(
        0.562 + molten_anchor['anchor_offset_eV'] + 0.15
    )

    # 11. Report generator surfaces calibrated interfacial scores, PBR / Ergun residence diagnostics, and live LCOH scenarios
    from pipeline.evidence.report_generator import (
        _reactor_section,
        _screening_section,
        _tea_section,
    )

    scr_lines = '\n'.join(
        _screening_section({
            'branch_ranker_evidence': pd.DataFrame([
                {
                    'genome': "('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 0, 0)",
                    'material_class': 'SolidCatalyst',
                    'E_act': 0.386,
                    'E_act_raw': 0.386,
                    'c_height_above_surface_A': 0.20,
                    'E_act_censored': False,
                    'pyrolysis_viable': True,
                    'valid': True,
                }
            ])
        })
    )
    assert 'Interfacial Score (eV)' in scr_lines
    assert 'Calibrated Score (eV)' in scr_lines
    assert '1.2500' in scr_lines

    rx_lines = '\n'.join(
        _reactor_section({'solids_scorecard': scorecard_with_diag}, [pbr_rec])
    )
    assert (
        'Packed-Bed (PBR), Ergun Pilot Residence-Time, & Regenerated-Cycle Diagnostics'
        in rx_lines
    )
    assert '0.720' in rx_lines

    tea_lines = '\n'.join(
        _tea_section({'best_conversion': 0.65, 'solids_scorecard': scorecard_with_diag})
    )
    assert 'Levelized Cost of Hydrogen (LCOH)' in tea_lines
    assert '| base |' in tea_lines

    # 12. Multipathway champion scorecard includes calibrated interfacial scores and PBR / Ergun residence diagnostics
    import json
    from pipeline.economics.tea import build_multipathway_champion_scorecard

    mp_root = tmp_path / 'mp_results'
    (mp_root / 'screening').mkdir(parents=True, exist_ok=True)
    (mp_root / 'reactor').mkdir(parents=True, exist_ok=True)
    pd.DataFrame([
        {
            'candidate_id': 'ni_al2o3_champ',
            'genome': "('SolidCatalyst', 'Ni', 'Al2O3', 'fcc111', 0.0, (), 0, 0)",
            'material_class': 'SolidCatalyst',
            'E_act': 0.386,
            'E_act_raw': 0.386,
            'c_height_above_surface_A': 0.20,
            'coking_index': 1.85,
            'E_act_censored': False,
            'pyrolysis_viable': True,
            'valid': True,
        }
    ]).to_csv(mp_root / 'screening' / 'branch_ranker_evidence.csv', index=False)
    (mp_root / 'reactor' / 'pbr_champ.json').write_text(
        json.dumps({
            **pbr_rec,
            'ergun_delta_p_bar': 0.0125,
            'H2_selectivity': 0.995,
        }),
        encoding='utf-8',
    )
    mp_scorecard = build_multipathway_champion_scorecard(str(mp_root))
    assert len(mp_scorecard['discovery_champions']) == 1
    champ0 = mp_scorecard['discovery_champions'][0]
    assert champ0['calibrated_activation_eV'] == pytest.approx(1.05)
    assert champ0['calibrated_interfacial_score_eV'] == pytest.approx(1.25)
    pbr_lead = mp_scorecard['pathway_leaderboard']['thermocatalytic_pbr']
    assert pbr_lead['status'] == 'complete'
    assert pbr_lead['cycle_averaged_CH4_conversion'] == pytest.approx(0.39)
    assert pbr_lead['overall_effectiveness_factor'] == pytest.approx(0.72)
    assert pbr_lead['thiele_modulus'] == pytest.approx(1.4)
    assert pbr_lead['ergun_delta_p_bar'] == pytest.approx(0.0125)
    assert mp_scorecard['pbr_ergun_residence_diagnostics']['n_pbr_records'] == 1

    # 13. IAPWS-95 aqueous critical point (647.096 K) enforced in electrochemical_evidence
    from pipeline.electrochemistry.model import (
        ElectrochemicalConditions,
        electrochemical_evidence,
    )

    base_ec_kwargs = {
        'electrolyte_identity': '1 M KOH',
        'applied_potential_V': 1.1,
        'current_density_A_cm2': 0.2,
        'faradaic_efficiency_H2': 0.90,
        'methane_conversion': 0.15,
        'pressure_Pa': 101325.0,
        'measurement_source': 'bench_test',
        'paired_control_source': 'bench_control',
    }
    aq_supercritical = electrochemical_evidence(
        ElectrochemicalConditions(
            electrolyte_phase='aqueous',
            temperature_K=700.0,
            **base_ec_kwargs,
        )
    )
    assert aq_supercritical['status'] == 'invalid'
    assert '647.096' in str(aq_supercritical['reason'])
    molten_high_t = electrochemical_evidence(
        ElectrochemicalConditions(
            electrolyte_phase='molten',
            temperature_K=900.0,
            **base_ec_kwargs,
        )
    )
    assert molten_high_t['status'] == 'measured_paired_control'

    # 14. Industrial viability gates support production reactor record field aliases
    from pipeline.validation.viability import evaluate_fuel_cell, evaluate_turquoise

    prod_rx_viable = evaluate_turquoise({
        'T_K': 1073.15,
        'single_pass_CH4_conversion': 0.78,
        'deactivation_rate_per_h': 0.004,
        'coke_selectivity': 0.02,
        'specific_energy_kWh_kg_H2': 11.8,
        'measured_reactor': 1.0,
    })
    assert prod_rx_viable['status'] == 'pass'
    assert prod_rx_viable['observed']['ch4_conversion'] == pytest.approx(0.78)
    assert prod_rx_viable['observed']['net_energy'] == pytest.approx(11.8)

    fc_viable = evaluate_fuel_cell({
        'orr_overpotential_V': 0.34,
        'peak_power_W_cm2': 1.15,
        'system_efficiency': 0.48,
        'voltage_degradation_uV_h': 6.5,
        'durability_hours': 150.0,
        'measured_mea': 1.0,
    })
    assert fc_viable['status'] == 'pass'
    assert fc_viable['observed']['measured_hours'] == pytest.approx(150.0)

    # 15. ReactorBatchState and FuelCellState populate complete telemetry contracts
    from pipeline.stages.fuel_cell import FuelCellServices, run_fuel_cell_stage
    from pipeline.stages.reactor_batch import (
        ReactorBatchServices,
        run_reactor_batch_stage,
    )

    rx_candidates = pd.DataFrame([
        {'candidate_id': 'ni111_sio2_bengaard2002_ref', 'E_act': 0.90},
        {'candidate_id': 'cand_co_ce_01', 'E_act': 0.62},
    ])
    rx_services = ReactorBatchServices(
        prepare_gas_mechanism=lambda: None,
        simulate_candidate=lambda row, name, temps, rtypes, **kw: {
            'sweep': [
                {
                    'reactor_type': 'PFR',
                    'catalyst_name': name,
                    'T_K': 973.15,
                    'CH4_conversion': 0.26 if name == 'ni111_sio2_bengaard2002_ref' else 0.75,
                    'single_pass_CH4_conversion': 0.26 if name == 'ni111_sio2_bengaard2002_ref' else 0.75,
                    'deactivation_rate_per_h': 0.005,
                    'coke_selectivity': 0.02,
                    'specific_energy_kWh_kg_H2': 11.2,
                    'measured_reactor': 1.0,
                    'status': 'complete',
                    'surface_loaded': True,
                    'exceeds_equilibrium': False,
                }
            ]
        },
        write_mock_mechanism=lambda *a, **kw: None,
        run_mock_sweep=lambda *a, **kw: [],
        build_scorecard=build_solids_scorecard,
    )
    rx_outcome = run_reactor_batch_stage(
        rx_candidates,
        temperatures=(973.15,),
        reactor_types=('PFR',),
        pathway_mode='thermocatalytic_pfr',
        multiphysics_results_dir=str(tmp_path),
        allow_mock_inputs=False,
        judge_catalyst='ni111_sio2_bengaard2002_ref',
        headline_t_min=923.15,
        headline_t_max=973.15,
        services=rx_services,
    )
    assert rx_outcome.state['best_conversion'] == pytest.approx(0.26)
    assert rx_outcome.state['candidate_champion_conversion'] == pytest.approx(0.75)
    assert rx_outcome.state['candidate_champion_catalyst'] == 'cat_1'
    assert rx_outcome.state['best_h2_cost_usd_kg'] > 0.0
    assert rx_outcome.state['industrial_viability']['pass'] == 1
    assert rx_outcome.state['industrial_viability']['fail'] == 1

    fc_cathodes = pd.DataFrame([
        {
            'valid': True,
            'name': 'PtCo_SAA',
            'orr_overpotential_V': 0.31,
            'orr_overpotential_gc_V': 0.29,
            'pourbaix_acid_stable': True,
            'material_class': 'SAA',
        }
    ])
    fc_services = FuelCellServices(
        screen_cathodes=lambda: fc_cathodes,
        sweep_membranes=lambda name, eta, **kw: [
            {
                'cathode_catalyst': name,
                'membrane': 'Nafion_211',
                'orr_overpotential_V': eta,
                'peak_power_W_cm2': 1.28,
                'efficiency_at_peak': 0.51,
                'peak_voltage_V': 0.68,
                'peak_current_A_cm2': 1.5,
            }
        ],
        build_stack_config=lambda **kw: kw,
        model_stack=lambda cfg: {'net_power_kW': 84.5, 'system_efficiency': 0.49},
    )
    fc_outcome = run_fuel_cell_stage(
        top_k_pemfc=1, stack_cells=380, services=fc_services
    )
    assert fc_outcome.state['best_catalyst'] == 'PtCo_SAA'
    assert fc_outcome.state['best_membrane'] == 'Nafion_211'
    assert fc_outcome.state['best_overpotential_V'] == pytest.approx(0.31)
    assert fc_outcome.state['best_overpotential_gc_V'] == pytest.approx(0.29)
    assert fc_outcome.state['pourbaix_acid_stable_count'] == 1
    assert fc_outcome.state['stack_net_kW'] == pytest.approx(84.5)
    assert fc_outcome.state['stack_efficiency'] == pytest.approx(0.49)


def test_literature_report_scorecard_and_champion_neb_calibration(tmp_path):
    """Verify literature-backed catalyst labels in report_generator, 7-pathway scorecard rendering, and QE NEB barrier calibration."""
    import json
    from pipeline.economics.tea import build_multipathway_champion_scorecard
    from pipeline.evidence.report_generator import (
        _catalyst_label,
        _dft_section,
        _tea_section,
    )
    from pipeline.validation.dft_validator import (
        build_champion_neb_slab,
        generate_champion_neb_input,
        parse_champion_neb_barrier,
    )

    pd = pytest.importorskip('pandas')

    # 1. Literature catalyst label resolution handles campaign_ and orr_slate_ prefixes
    lbl = _catalyst_label(
        'campaign_cand_08da71e8eb47',
        material_class='MoltenMetal',
        genome="('MoltenMetal', 'Bi', 'Ni', 25.0, 1073.15)",
    )
    assert 'Liquid Metal' in lbl and 'Bi' in lbl and 'Ni' in lbl

    dft_lines = '\n'.join(
        _dft_section([
            {
                'catalyst_name': 'campaign_cand_08da71e8eb47',
                'material_class': 'MoltenMetal',
                'genome': "('MoltenMetal', 'Bi', 'Ni', 25.0, 1073.15)",
                'dft_energy_Ry': -128.45,
                'dft_energy_eV': -1747.65,
                'converged': True,
                'max_force_Ry_bohr': 0.0004,
            }
        ])
    )
    assert 'Liquid Metal' in dft_lines and 'Bi' in dft_lines

    # 2. Candidate-specific QE NEB slab builder, input writer, and barrier parser
    genome_mm = ('MoltenMetal', 'Bi', 'Ni', 25.0, 1073.15)
    slab, active_idx = build_champion_neb_slab(genome_mm)
    assert len(slab) == 8
    assert 0 <= active_idx < len(slab)
    assert slab[active_idx].symbol == 'Ni'

    neb_meta = generate_champion_neb_input(
        'champ_nibi', genome_mm, n_images=5, output_dir=tmp_path
    )
    neb_in_file = tmp_path / 'champ_nibi.neb.in'
    assert neb_in_file.is_file()
    assert neb_meta['n_images'] == 5
    neb_in_text = neb_in_file.read_text(encoding='utf-8')
    assert 'BEGIN_PATH_INPUT' in neb_in_text
    assert 'num_of_images=5' in neb_in_text

    fake_neb_out = tmp_path / 'champ_nibi.neb.out'
    fake_neb_out.write_text(
        '     activation energy (->) =   0.9850 eV\n'
        '     activation energy (<-) =   0.6120 eV\n'
        '     path length            =   2.4500 bohr\n'
        '     neb: convergence achieved in 18 iterations\n'
        '   JOB DONE.\n',
        encoding='utf-8',
    )
    neb_parsed = parse_champion_neb_barrier(str(fake_neb_out))
    assert neb_parsed['neb_converged'] is True
    assert neb_parsed['neb_forward_barrier_eV'] == pytest.approx(0.985)
    assert neb_parsed['neb_reverse_barrier_eV'] == pytest.approx(0.612)
    assert neb_parsed['kinetics_status'] == 'dft_neb_complete'

    # 3. Multipathway scorecard upgrades kinetics_status to dft_neb_complete when converged NEB barrier exists
    mp_root = tmp_path / 'mp_neb_results'
    (mp_root / 'screening').mkdir(parents=True, exist_ok=True)
    (mp_root / 'reactor').mkdir(parents=True, exist_ok=True)
    (mp_root / 'dft').mkdir(parents=True, exist_ok=True)

    pd.DataFrame([
        {
            'candidate_id': 'cand_neb_01',
            'genome': "('MoltenMetal', 'Bi', 'Ni', 25.0, 1073.15)",
            'material_class': 'MoltenMetal',
            'E_act': 0.75,
            'E_act_raw': 0.75,
            'c_height_above_surface_A': 0.10,
            'coking_index': 1.90,
            'E_act_censored': False,
            'pyrolysis_viable': True,
            'valid': True,
        }
    ]).to_csv(mp_root / 'screening' / 'branch_ranker_evidence.csv', index=False)

    (mp_root / 'dft' / 'campaign_cand_neb_01_dft.json').write_text(
        json.dumps({
            'catalyst_name': 'campaign_cand_neb_01',
            'candidate_id': 'cand_neb_01',
            'material_class': 'MoltenMetal',
            'converged': True,
            'dft_energy_Ry': -210.5,
            'dft_energy_eV': -2864.0,
            'neb_converged': True,
            'neb_forward_barrier_eV': 0.985,
            'kinetics_status': 'dft_neb_complete',
        }),
        encoding='utf-8',
    )

    (mp_root / 'reactor' / 'mmbcr_neb.json').write_text(
        json.dumps({
            'status': 'complete',
            'reactor_type': 'MMBCR',
            'catalyst_name': 'cand_neb_01',
            'material_class': 'MoltenMetal',
            'T_K': 1173.15,
            'CH4_conversion': 0.82,
            'single_pass_CH4_conversion': 0.82,
            'H2_selectivity': 0.998,
            'surface_loaded': True,
            'exceeds_equilibrium': False,
        }),
        encoding='utf-8',
    )

    mp_scorecard = build_multipathway_champion_scorecard(str(mp_root))
    assert mp_scorecard['discovery_champions'][0]['qe_neb_converged'] is True
    assert mp_scorecard['discovery_champions'][0]['qe_neb_barrier_eV'] == pytest.approx(0.985)
    mmbcr_cal = mp_scorecard['dft_calibrated_pathway_leaderboard']['mmbcr']
    assert mmbcr_cal['status'] == 'complete'
    assert mmbcr_cal['kinetics_status'] == 'dft_neb_complete'
    # Barrier shift = (0.985 - 0.75) + 0.15 * 0.10 = 0.250 eV
    assert mmbcr_cal['dft_barrier_shift_eV'] == pytest.approx(0.25, abs=1e-3)

    tea_md = '\n'.join(_tea_section({'best_conversion': 0.82}, mp_scorecard))
    assert 'Multi-Pathway Discovery Champions & QE-DFT-Calibrated Leaderboard' in tea_md
    assert 'dft_neb_complete' in tea_md
    assert '7-Pathway H₂/C₂ Selectivity, Yield, & Volumetric / Site-Normalized Production Rates' in tea_md


def test_h2_c2_selectivity_and_volumetric_production_rates(tmp_path: Path) -> None:
    """Verify expansion-corrected H-atom balance, H2/C2 branching selectivity, H2 yield, and volumetric/site-normalized rates."""
    import pytest
    from pipeline.data_models.reactors import (
        H2C2SelectivityRateMetrics,
        ReactorResult,
    )
    from pipeline.reactors.carbon_accounting import (
        compute_h2_c2_selectivity_and_rates,
        h2_atom_balance_metric,
    )
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import ReactorConfig, simulate_mmbcr

    # 1. Expansion-corrected H2 atom balance fixes the 1 -> 2 mol dilution artifact
    x_conv: float = 0.7383
    n_tot_out: float = 1.0 + 0.95 * x_conv
    x_h2_out: float = (2.0 * 0.95 * x_conv) / n_tot_out
    x_ar_out: float = 0.05 / n_tot_out
    # Without Ar tracer argument (using exact stoichiometry expansion 1 + 0.95*X)
    assert h2_atom_balance_metric(0.95, x_conv, x_h2_out) == pytest.approx(1.0, abs=1e-5)
    # With Ar tracer argument
    assert h2_atom_balance_metric(
        0.95, x_conv, x_h2_out, x_ar_feed=0.05, x_ar_out=x_ar_out
    ) == pytest.approx(1.0, abs=1e-5)

    # 2. Direct selectivity & volumetric/site-normalized rate calculation with C2H4 byproduct
    metrics: H2C2SelectivityRateMetrics = compute_h2_c2_selectivity_and_rates(
        None,
        feed_initial=0.95,
        ar_initial=0.05,
        final_conv=0.50,
        tau_s=1.25,
        P_inlet_Pa=101325.0,
        T_inlet_K=1173.15,
        active_sv_1_m=3600.0,
        site_density_mol_cm2=2.5e-9,
        bed_porosity_or_holdup=0.4,
        x_h2_override=0.625,
        x_feed_override=0.320,
        x_c2h4_override=0.005,
        x_ar_override=0.034,
    )
    assert 0.90 < metrics['H2_selectivity'] < 1.0
    assert metrics['C2_carbon_selectivity'] > 0.0
    assert metrics['solid_C_selectivity'] + metrics['C2_carbon_selectivity'] == pytest.approx(
        1.0, abs=1e-5
    )
    assert metrics['H2_yield'] == pytest.approx(
        0.50 * metrics['H2_selectivity'], abs=1e-4
    )
    assert metrics['ch4_conversion_rate_mol_m3_s'] > 0.0
    assert metrics['h2_production_rate_mol_m3_s'] == pytest.approx(
        2.0 * metrics['ch4_conversion_rate_mol_m3_s'] * metrics['H2_selectivity'],
        rel=1e-3,
    )
    assert metrics['h2_production_rate_kg_m3_h'] > 0.0
    tof_val: float | None = metrics['h2_turnover_frequency_1_s']
    assert tof_val is not None and tof_val > 0.0

    # 3. MMBCR simulation includes competitive bubble-interior C2 pyrolysis and full rate telemetry
    mech_path: str = write_full_mechanism(
        'bi73ni27_test', E_act_CH4=0.45, output_dir=tmp_path
    )
    cfg: ReactorConfig = ReactorConfig(
        T_inlet_K=1250.0,
        P_inlet_Pa=101325.0,
        reactor_type='MMBCR',
        catalyst_name='bi73ni27_test',
        mechanism_file=str(mech_path),
        material_class='MoltenMetal',
        pathway_mode='mmbcr',
        catalyst_E_act_eV=0.45,
    )
    mmbcr_res: ReactorResult = simulate_mmbcr(cfg)
    assert mmbcr_res['CH4_conversion'] > 0.20
    assert 0.85 <= mmbcr_res['H2_selectivity'] <= 1.0
    assert mmbcr_res['H2_atom_balance'] == pytest.approx(
        mmbcr_res['H2_selectivity'], abs=1e-4
    )
    assert mmbcr_res['C2_carbon_selectivity'] > 0.0
    solid_c_sel_mmbcr: float | None = mmbcr_res['solid_C_selectivity']
    assert solid_c_sel_mmbcr is not None
    assert solid_c_sel_mmbcr + mmbcr_res['C2_carbon_selectivity'] == pytest.approx(
        1.0, abs=1e-4
    )
    assert mmbcr_res['H2_yield'] == pytest.approx(
        mmbcr_res['CH4_conversion'] * mmbcr_res['H2_selectivity'],
        abs=1e-3,
    )
    assert mmbcr_res['ch4_conversion_rate_mol_m3_s'] > 0.0
    assert mmbcr_res['h2_production_rate_mol_m3_s'] > 0.0
    assert mmbcr_res['h2_production_rate_kg_m3_h'] > 0.0
    tof_mmbcr: float | None = mmbcr_res['h2_turnover_frequency_1_s']
    assert tof_mmbcr is not None and tof_mmbcr > 0.0
    assert mmbcr_res['c2_tracked'] is True


def test_multi_hydrocarbon_feedstocks_default_to_methane(tmp_path: Path) -> None:
    """Verify multi-hydrocarbon cracking (18 feedstocks including biomass & waste plastic) with CH4 as strict default."""
    import pytest
    from pipeline.data_models.evidence import TEAScenarioEstimate
    from pipeline.data_models.reactors import FeedstockName, ReactorResult
    from pipeline.economics.tea import estimate_hydrogen_cost
    from pipeline.reactors.feedstocks import (
        DEFAULT_FEEDSTOCK,
        adapt_reactor_result_to_feedstock,
        feedstock_equilibrium_conversion,
        resolve_feedstock,
    )
    from pipeline.reactors.mechanisms import write_full_mechanism
    from pipeline.reactors.models import (
        ReactorConfig,
        simulate_mmbcr,
        simulate_multi_feedstock_comparison,
    )

    # 1. Default feedstock is strictly CH4, and all 18 hydrocarbon / waste / biomass feedstocks resolve
    assert DEFAULT_FEEDSTOCK == 'CH4'
    assert resolve_feedstock().formula == 'CH4'
    assert resolve_feedstock('methane').formula == 'CH4'
    assert resolve_feedstock('flare_gas').formula == 'ASSOC_GAS'
    assert resolve_feedstock('propane').formula == 'C3H8'
    assert resolve_feedstock('pentane').formula == 'C5H12'
    assert resolve_feedstock('hexane').formula == 'C6H14'
    assert resolve_feedstock('cyclohexane').formula == 'C6H12'
    assert resolve_feedstock('octane').formula == 'C8H18'
    assert resolve_feedstock('diesel').formula == 'C12H26'
    assert resolve_feedstock('cetane').formula == 'C16H34'
    assert resolve_feedstock('waste_plastic').formula == 'C10H20_PE'
    assert resolve_feedstock('biomethane').formula == 'BIOMETHANE'
    assert resolve_feedstock('biogas').formula == 'BIOGAS'
    assert resolve_feedstock('waste_cooking_oil').formula == 'C18H34O2_LIPID'
    assert resolve_feedstock('glycerol').formula == 'C3H8O3_GLYCEROL'
    assert resolve_feedstock('bio_oil').formula == 'C6H10O2_BIOOIL'
    assert resolve_feedstock('biomass').formula == 'C6H10O4_BIOMASS'
    assert resolve_feedstock('diesel').nu_H2 == pytest.approx(13.0)
    assert resolve_feedstock('diesel').nu_C == pytest.approx(12.0)
    assert resolve_feedstock('biomass').nu_H2 == pytest.approx(5.0)
    assert resolve_feedstock('biomass').nu_C == pytest.approx(2.0)
    assert resolve_feedstock('biomass').nu_CO == pytest.approx(4.0)

    with pytest.raises(ValueError, match='Unsupported hydrocarbon feedstock'):
        resolve_feedstock('unobtainium')

    # 2. Equilibrium conversion at 850 K increases from CH4 -> C2H6 -> C3H8 -> C6H14 -> C12H26
    x_eq_ch4: float = feedstock_equilibrium_conversion(850.0, 101325.0, 'CH4')
    x_eq_c2: float = feedstock_equilibrium_conversion(850.0, 101325.0, 'C2H6')
    x_eq_c3: float = feedstock_equilibrium_conversion(850.0, 101325.0, 'C3H8')
    x_eq_c6: float = feedstock_equilibrium_conversion(850.0, 101325.0, 'C6H14')
    x_eq_c12: float = feedstock_equilibrium_conversion(850.0, 101325.0, 'diesel')
    assert x_eq_ch4 < x_eq_c2 <= x_eq_c3 <= x_eq_c6 <= x_eq_c12

    # 3. Default ReactorConfig uses CH4 and multi-feedstock comparison adapts across all 18 feedstocks
    mech_path: str = write_full_mechanism(
        'bi73ni27_feed_test', E_act_CH4=0.65, output_dir=tmp_path
    )
    cfg_default: ReactorConfig = ReactorConfig(
        T_inlet_K=1050.0,
        P_inlet_Pa=101325.0,
        reactor_type='MMBCR',
        catalyst_name='bi73ni27_feed_test',
        mechanism_file=str(mech_path),
        material_class='MoltenMetal',
        pathway_mode='mmbcr',
        catalyst_E_act_eV=0.65,
    )
    assert cfg_default.feedstock == 'CH4'
    res_ch4: ReactorResult = simulate_mmbcr(cfg_default)
    assert res_ch4['feedstock'] == 'CH4'
    assert res_ch4['feedstock_conversion'] == pytest.approx(res_ch4['CH4_conversion'], abs=1e-6)
    assert res_ch4['single_pass_feedstock_conversion'] == pytest.approx(
        res_ch4['CH4_conversion'], abs=1e-6
    )
    assert res_ch4['exit_x_feedstock'] == pytest.approx(res_ch4['exit_x_CH4'], abs=1e-6)
    assert res_ch4['feedstock_conversion_rate_mol_m3_s'] == pytest.approx(
        res_ch4['ch4_conversion_rate_mol_m3_s'], abs=1e-6
    )
    assert res_ch4['endothermic_heat_duty_kJ_per_mol_feed_fed'] == pytest.approx(
        res_ch4['endothermic_heat_duty_kJ_per_mol_CH4_fed'], abs=1e-6
    )

    comp: dict[FeedstockName, ReactorResult] = simulate_multi_feedstock_comparison(
        cfg_default
    )
    assert len(comp) == 18
    assert comp['CH4']['CH4_conversion'] == pytest.approx(res_ch4['CH4_conversion'], abs=1e-6)
    assert comp['BIOMETHANE']['feedstock_conversion'] == pytest.approx(
        res_ch4['CH4_conversion'], abs=1e-6
    )
    assert comp['C3H8']['feedstock_conversion'] > comp['CH4']['feedstock_conversion']
    assert comp['C3H8']['single_pass_feedstock_conversion'] == pytest.approx(
        comp['C3H8']['feedstock_conversion'], abs=1e-6
    )
    assert comp['C3H8']['exit_x_feedstock'] == pytest.approx(
        comp['C3H8']['exit_x_CH4'], abs=1e-6
    )
    assert comp['C6H14']['feedstock_conversion'] > comp['C3H8']['feedstock_conversion']
    assert comp['C12H26']['h2_production_rate_mol_m3_s'] > comp['CH4']['h2_production_rate_mol_m3_s']
    assert comp['C6H10O4_BIOMASS']['exit_x_CO'] > 0.0
    assert comp['CH4']['exit_x_CO'] == pytest.approx(0.0)
    from dataclasses import replace

    direct_adapted: ReactorResult = adapt_reactor_result_to_feedstock(res_ch4, 'C6H14')
    assert direct_adapted['feedstock_conversion'] == pytest.approx(
        comp['C6H14']['feedstock_conversion'], abs=1e-6
    )
    direct_sim_c6: ReactorResult = simulate_mmbcr(replace(cfg_default, feedstock='C6H14'))
    assert direct_sim_c6['feedstock'] == 'C6H14'
    assert direct_sim_c6['feedstock_conversion'] == pytest.approx(
        comp['C6H14']['feedstock_conversion'], abs=1e-6
    )

    # 4. TEA defaults strictly to CH4 and scales carbon co-product credit & cheap waste/biomass feedstocks
    tea_default: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.55, 'base', temperature_K=1100.0
    )
    tea_ch4: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.55, 'base', temperature_K=1100.0, feedstock='CH4'
    )
    tea_flare: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.55, 'base', temperature_K=1100.0, feedstock='flare_gas'
    )
    tea_plastic: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.85, 'base', temperature_K=1100.0, feedstock='waste_plastic'
    )
    tea_biomass: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.85, 'base', temperature_K=1100.0, feedstock='biomass'
    )
    tea_propane: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.55, 'base', temperature_K=1100.0, feedstock='C3H8', psa_feedstock_recovery=0.95
    )
    tea_diesel: TEAScenarioEstimate = estimate_hydrogen_cost(
        0.55, 'base', temperature_K=1100.0, feedstock='diesel'
    )
    assert tea_default['h2_cost_usd_kg'] == tea_ch4['h2_cost_usd_kg']
    assert tea_default.get('net_feedstock_utilization') == tea_default['net_ch4_utilization']
    assert tea_propane.get('psa_feedstock_recovery') == pytest.approx(0.95)
    assert tea_flare['h2_cost_usd_kg'] < tea_ch4['h2_cost_usd_kg']
    assert tea_plastic['h2_cost_usd_kg'] < tea_diesel['h2_cost_usd_kg']
    assert tea_biomass['h2_cost_usd_kg'] < tea_diesel['h2_cost_usd_kg']
    assert tea_diesel['kg_c_per_kg_h2'] > tea_propane['kg_c_per_kg_h2'] > tea_ch4['kg_c_per_kg_h2']
    assert tea_diesel['carbon_credit_usd_kg_h2'] > tea_ch4['carbon_credit_usd_kg_h2']


def test_ntec_electrochemical_multicarbon_mechanisms_and_pemfc_coupling(
    tmp_path: Path,
) -> None:
    """Verify Option 1 (NTEC & Electrochemical coupled operating points), Option 2 (explicit Cantera C2/C3/biogas mechanisms), Option 3 (Reactor-PSA-PEMFC coupling), and Option 4 (report tables)."""
    import pytest
    from pipeline.electrochemistry.native_models import (
        coupled_electrochemical_operating_point,
        coupled_ntec_operating_point,
    )
    from pipeline.evidence.report_generator import _tea_section
    from pipeline.fuel_cell.stack import couple_reactor_to_pemfc_stack
    from pipeline.reactors.mechanisms import (
        write_feedstock_mechanism,
        write_full_mechanism,
    )

    # 1. Option 1: Coupled NTEC and Electrochemical native operating point solvers
    ntec_pt = coupled_ntec_operating_point(
        {
            'catalyst_E_act_eV': 0.55,
            'catalyst_dE_C_eV': 2.80,
            'T_K': 473.15,
            'P_Pa': 506625.0,
            'interface_area_m2': 2.7,
        }
    )
    assert 0.30 < ntec_pt['CH4_conversion'] < 1.0
    assert ntec_pt['field_rate_enhancement'] > 10.0
    assert ntec_pt['detachment_regime'] == 'detaching'

    ec_pt = coupled_electrochemical_operating_point(
        {
            'catalyst_E_act_eV': 0.55,
            'catalyst_dE_C_eV': 2.80,
            'T_K': 823.15,
            'P_Pa': 506625.0,
            'electrode_area_m2': 1.8,
            'applied_potential_V': 0.65,
        }
    )
    assert 0.30 < ec_pt['CH4_conversion'] < 1.0
    assert ec_pt['specific_energy_kWh_kg_H2'] < 20.0
    assert ec_pt['electrolyte_phase'] == 'molten'

    # 2. Option 2: Explicit Cantera 3.x multi-carbon (C2H6, C3H8, associated_gas) and biogas mechanisms
    ct = pytest.importorskip('cantera')
    ch4_mech = write_feedstock_mechanism(
        'test_ch4_cat', 'CH4', E_act_CH4=0.75, output_dir=tmp_path
    )
    ch4_ref = write_full_mechanism(
        'test_ch4_cat_ref', E_act_CH4=0.75, output_dir=tmp_path
    )
    assert Path(ch4_mech).is_file() and Path(ch4_ref).is_file()

    c3_mech = write_feedstock_mechanism(
        'test_c3_cat', 'C3H8', E_act_CH4=0.75, output_dir=tmp_path
    )
    gas_c3 = ct.Solution(str(c3_mech), 'gas')
    surf_c3 = ct.Interface(str(c3_mech), 'test_c3_cat_surface', [gas_c3])
    assert 'C3H8' in gas_c3.species_names and 'C3H6' in gas_c3.species_names
    assert surf_c3.n_reactions >= 8

    biogas_mech = write_feedstock_mechanism(
        'test_biogas_cat', 'biogas', E_act_CH4=0.75, output_dir=tmp_path
    )
    gas_bg = ct.Solution(str(biogas_mech), 'gas')
    surf_bg = ct.Interface(str(biogas_mech), 'test_biogas_cat_surface', [gas_bg])
    for sp in ('CH4', 'CO2', 'CO', 'H2O', 'H2'):
        assert sp in gas_bg.species_names
    for ssp in ('O_s', 'OH_s', 'CO_s'):
        assert ssp in surf_bg.species_names

    # 3. Option 3: Integrated Reactor-PSA-PEMFC coupling and CO/C2H4 impurity tolerance
    coupled_ch4 = couple_reactor_to_pemfc_stack(
        {
            'reactor_type': 'MMBCR',
            'catalyst_name': '99cbb4a430ec9289b09d',
            'T_K': 1273.15,
            'dft_calibrated_single_pass_CH4_conversion': 0.9032,
            'H2_selectivity': 0.985,
            'C2_carbon_selectivity': 0.0301,
            'dft_calibrated_tea_base_h2_usd_kg': 1.249,
            'dft_calibrated_tea_optimistic_h2_usd_kg': 0.678,
        },
        feedstock='CH4',
    )
    assert coupled_ch4['post_psa_co_ppm'] == pytest.approx(0.0)
    assert coupled_ch4['post_psa_c2h4_ppm'] > 0.0
    assert coupled_ch4['anode_hor_activity_retention'] > 0.90
    assert coupled_ch4['stack_net_power_kW'] > 100.0
    assert 0.25 < coupled_ch4['well_to_power_lhv_efficiency'] < 0.70
    assert 0.03 < coupled_ch4['lcoe_total_usd_kwh_e'] < 0.20

    coupled_biomass = couple_reactor_to_pemfc_stack(
        {
            'reactor_type': 'MMBCR',
            'catalyst_name': '99cbb4a430ec9289b09d',
            'T_K': 1273.15,
            'single_pass_feedstock_conversion': 0.95,
            'H2_selectivity': 0.97,
            'C2_carbon_selectivity': 0.04,
        },
        feedstock='biomass',
    )
    assert coupled_biomass['post_psa_co_ppm'] > coupled_ch4['post_psa_co_ppm']
    assert (
        coupled_biomass['anode_hor_activity_retention']
        < coupled_ch4['anode_hor_activity_retention']
    )

    # 4. Option 4: Report generator renders 18-feedstock flexibility and Reactor-PSA-PEMFC coupling tables
    tea_md = '\n'.join(
        _tea_section(
            {'best_conversion': 0.9032},
            {
                'reactor_to_pemfc_system_coupling': {
                    'mmbcr': {'pathway': 'mmbcr', **coupled_ch4}
                },
                'feedstock_flexibility_matrix': [
                    {
                        'feedstock': 'CH4',
                        'feedstock_name': 'Methane (natural gas)',
                        'category': 'natural_gas',
                        'h_to_c_ratio': 4.0,
                        'kg_feed_per_kg_h2': 3.979,
                        'kg_c_per_kg_h2': 2.979,
                        'dh_rxn_298_kJ_mol_h2': 37.43,
                        'pfr_conversion': 0.5625,
                        'pfr_lcoh_base_usd_kg': 1.246,
                        'mmbcr_conversion': 0.9032,
                        'mmbcr_lcoh_base_usd_kg': 1.249,
                        'mmbcr_lcoh_optimistic_usd_kg': 0.678,
                        'post_psa_co_ppm': 0.0,
                        'coupled_lcoe_base_usd_kwh_e': coupled_ch4[
                            'lcoe_total_usd_kwh_e'
                        ],
                    }
                ],
            },
        )
    )
    assert (
        'Integrated Reactor-PSA-PEMFC Well-to-Power Coupling & CO/C₂H₄ Impurity Tolerance'
        in tea_md
    )
    assert (
        '18-Feedstock Multi-Hydrocarbon, Waste Polyolefin & Biomass Flexibility Matrix'
        in tea_md
    )
