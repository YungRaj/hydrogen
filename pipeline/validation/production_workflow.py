"""Resumable, fail-closed orchestration for production QE validation."""

from __future__ import annotations

import json
import hashlib
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read as ase_read
from ase.mep import NEB

from pipeline.validation.dft_validator import parse_convergence
from pipeline.validation.qe_workflows import (
    QEExecutionConfig,
    _fixed_atom_indices,
    parse_atomic_forces,
    parse_neb_result,
    partial_hessian,
    relaxed_structure,
    run_neb,
    run_pw,
    write_qe_force_input,
    write_qe_neb_input,
    write_qe_relax_input,
)


ORR_STAGES = ('clean', 'OH', 'O', 'OOH', 'h2', 'h2o')
PYROLYSIS_ELEMENTARY_STEPS = {
    'methane_activation': 'methane_activation_eV',
    'methyl_dehydrogenation': 'ch3_dehydrogenation_eV',
    'methylene_dehydrogenation': 'ch2_dehydrogenation_eV',
    'methylidyne_dehydrogenation': 'ch_dehydrogenation_eV',
    'hydrogen_desorption': 'h2_desorption_eV',
    'carbon_transport': 'carbon_transfer_eV',
}


def _manifest_digest(payload: dict) -> str:
    return hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()


def _pw_execution(execution: QEExecutionConfig | None) -> QEExecutionConfig | None:
    """Remove NEB-only image grouping when launching pw.x stages."""
    if execution is None:
        return None
    return QEExecutionConfig(
        mpi_ranks=execution.mpi_ranks,
        omp_threads=execution.omp_threads,
        kpoint_pools=execution.kpoint_pools,
        image_groups=1,
    )


def qe_output_status(path: str | Path) -> str:
    """Classify an output without mistaking a partial energy for evidence."""
    target = Path(path)
    if not target.exists() or target.stat().st_size == 0:
        return 'missing'
    text = target.read_text(errors='replace').lower()
    if 'error in routine' in text or 'convergence not achieved' in text:
        return 'failed'
    if parse_convergence(str(target), require_ionic=True):
        return 'converged'
    return 'incomplete'


def qe_scf_output_status(path: str | Path) -> str:
    """Classify a fixed-geometry force calculation without requiring BFGS."""
    target = Path(path)
    if not target.exists() or target.stat().st_size == 0:
        return 'missing'
    text = target.read_text(errors='replace').lower()
    if 'error in routine' in text or 'convergence not achieved' in text:
        return 'failed'
    return 'converged' if parse_convergence(str(target)) else 'incomplete'


def orr_campaign_status(calc_dir: str | Path, catalyst_name: str) -> dict:
    root = Path(calc_dir)
    stages = {
        stage: qe_output_status(root / f'{catalyst_name}_{stage}.out')
        for stage in ORR_STAGES
    }
    complete = all(value == 'converged' for value in stages.values())
    return {
        'stages': stages,
        'complete': complete,
        'orr_result_allowed': complete,
        'next_stage': next((key for key, value in stages.items()
                            if value != 'converged'), None),
    }


def run_orr_sequence(calc_dir: str | Path, catalyst_name: str,
                     timeout_s: int = 86400, restart_incomplete: bool = False,
                     execution: QEExecutionConfig | None = None) -> dict:
    """Run missing ORR stages in order and resume cleanly completed outputs.

    Nonempty incomplete outputs are left untouched by default because they may
    belong to a calculation currently running in another process.
    """
    root = Path(calc_dir)
    for stage in ORR_STAGES:
        input_path = root / f'{catalyst_name}_{stage}.in'
        output_path = root / f'{catalyst_name}_{stage}.out'
        state = qe_output_status(output_path)
        if state == 'converged':
            continue
        if state == 'incomplete' and not restart_incomplete:
            break
        if not input_path.is_file():
            raise FileNotFoundError(f'missing QE input: {input_path}')
        outcome = run_pw(str(input_path), str(output_path), timeout_s=timeout_s,
                         execution=_pw_execution(execution))
        if not outcome['converged']:
            break
    return orr_campaign_status(root, catalyst_name)


def methane_neb_status(calc_dir: str | Path) -> dict:
    root = Path(calc_dir)
    endpoints = {
        name: qe_output_status(root / f'{name}.relax.out')
        for name in ('initial', 'final')
    }
    neb_path = root / 'candidate.neb.out'
    neb = parse_neb_result(str(neb_path)) if neb_path.exists() else {
        'converged': False, 'forward_barrier_eV': None,
        'reverse_barrier_eV': None, 'candidate_specific': True,
    }
    frequency_path = root / 'transition_state_frequency.json'
    frequency = {'valid_transition_state': False, 'status': 'missing'}
    if frequency_path.exists():
        import json
        frequency = json.loads(frequency_path.read_text())
    endpoints_converged = all(x == 'converged' for x in endpoints.values())
    frequency_valid = bool(
        frequency.get('complete') and
        frequency.get('candidate_specific') and
        frequency.get('method') == 'central_finite_difference_qe_forces' and
        frequency.get('valid_transition_state'))
    return {
        'endpoints': endpoints,
        'endpoints_converged': endpoints_converged,
        'neb': neb,
        'frequency': frequency,
        'complete': bool(endpoints_converged and neb.get('converged') and
                         neb.get('candidate_specific') and frequency_valid),
    }


def prepare_methane_neb(calc_dir: str | Path, prefix: str,
                        initial: Atoms, final: Atoms) -> dict:
    """Prepare both endpoint relaxations from explicit candidate geometries."""
    if initial.get_chemical_symbols() != final.get_chemical_symbols():
        raise ValueError('NEB endpoints must have identical atoms and ordering')
    if len(initial) == 0 or not np.allclose(initial.cell.array, final.cell.array):
        raise ValueError('NEB endpoints require the same nonempty periodic cell')
    if _fixed_atom_indices(initial) != _fixed_atom_indices(final):
        raise ValueError('NEB endpoints require identical fixed-atom constraints')
    root = Path(calc_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_qe_relax_input(initial, str(root / 'initial.relax.in'),
                         f'{prefix}_initial')
    write_qe_relax_input(final, str(root / 'final.relax.in'),
                         f'{prefix}_final')
    metadata = {
        'schema_version': 1,
        'prefix': prefix,
        'atom_count': len(initial),
        'symbols': initial.get_chemical_symbols(),
        'status': 'endpoints_prepared',
    }
    (root / 'step_metadata.json').write_text(
        json.dumps(metadata, indent=2, sort_keys=True))
    return methane_neb_status(root)


def prepare_frequency_jobs(calc_dir: str | Path, transition_state: Atoms,
                           prefix: str, displacement_A: float = 0.01,
                           active_indices: list[int] | None = None) -> dict:
    """Prepare central finite-difference force jobs for a proposed TS."""
    if len(transition_state) == 0 or displacement_A <= 0:
        raise ValueError('transition state and positive displacement are required')
    root = Path(calc_dir)
    force_root = root / 'frequency_forces'
    force_root.mkdir(parents=True, exist_ok=True)
    active = (list(range(len(transition_state))) if active_indices is None
              else [int(index) for index in active_indices])
    if not active or len(set(active)) != len(active) or \
            any(index < 0 or index >= len(transition_state) for index in active):
        raise ValueError('active frequency atom indices are invalid')
    jobs = []
    for dof in range(3 * len(active)):
        active_position, axis = divmod(dof, 3)
        atom_index = active[active_position]
        for sign, suffix in ((1.0, 'plus'), (-1.0, 'minus')):
            displaced = transition_state.copy()
            displaced.positions[atom_index, axis] += sign * displacement_A
            stem = f'dof_{dof:04d}_{suffix}'
            write_qe_force_input(
                displaced, str(force_root / f'{stem}.in'), f'{prefix}_{stem}')
            jobs.append({'dof': dof, 'sign': suffix, 'stem': stem})
    manifest = {
        'schema_version': 1,
        'prefix': prefix,
        'displacement_A': displacement_A,
        'atom_count': len(transition_state),
        'active_indices': active,
        'masses_amu': transition_state.get_masses()[active].tolist(),
        'jobs': jobs,
    }
    manifest['manifest_sha256'] = _manifest_digest(manifest)
    (force_root / 'manifest.json').write_text(
        json.dumps(manifest, indent=2, sort_keys=True))
    return frequency_status(root)


def frequency_status(calc_dir: str | Path) -> dict:
    """Report finite-difference progress and accept only a first-order TS."""
    root = Path(calc_dir)
    manifest_path = root / 'frequency_forces/manifest.json'
    if not manifest_path.is_file():
        return {'status': 'missing', 'valid_transition_state': False,
                'complete': False}
    manifest = json.loads(manifest_path.read_text())
    result_path = root / 'transition_state_frequency.json'
    if result_path.is_file():
        result = json.loads(result_path.read_text())
        if result.get('force_manifest_sha256') == manifest.get('manifest_sha256'):
            result.setdefault('status', 'complete')
            return result
    states = {
        job['stem']: qe_scf_output_status(
            root / 'frequency_forces' / f"{job['stem']}.out")
        for job in manifest['jobs']}
    return {
        'status': 'forces_converged' if all(x == 'converged' for x in states.values())
                  else 'force_jobs_pending',
        'complete': False,
        'valid_transition_state': False,
        'job_counts': {state: list(states.values()).count(state)
                       for state in sorted(set(states.values()))},
        'jobs': states,
    }


def run_frequency_sequence(calc_dir: str | Path, timeout_s: int = 86400,
                           execution: QEExecutionConfig | None = None,
                           restart_incomplete: bool = False) -> dict:
    """Run/resume finite-difference jobs, then construct a validated Hessian."""
    root = Path(calc_dir)
    force_root = root / 'frequency_forces'
    manifest_path = force_root / 'manifest.json'
    if not manifest_path.is_file():
        return frequency_status(root)
    manifest = json.loads(manifest_path.read_text())
    for job in manifest['jobs']:
        input_path = force_root / f"{job['stem']}.in"
        output_path = force_root / f"{job['stem']}.out"
        state = qe_scf_output_status(output_path)
        if state == 'converged':
            continue
        if state == 'incomplete' and not restart_incomplete:
            return frequency_status(root)
        outcome = run_pw(str(input_path), str(output_path), timeout_s=timeout_s,
                         execution=_pw_execution(execution))
        if not outcome['converged']:
            return frequency_status(root)
    n_total = int(manifest['atom_count'])
    active = [int(index) for index in manifest['active_indices']]
    n_atoms = len(active)
    plus = np.empty((3 * n_atoms, n_atoms, 3), dtype=float)
    minus = np.empty_like(plus)
    for job in manifest['jobs']:
        forces = parse_atomic_forces(
            str(force_root / f"{job['stem']}.out"), expected_atoms=n_total)
        (plus if job['sign'] == 'plus' else minus)[int(job['dof'])] = forces[active]
    result = partial_hessian(
        plus, minus, float(manifest['displacement_A']),
        np.asarray(manifest['masses_amu'], dtype=float))
    result.update({
        'status': 'complete', 'complete': True,
        'candidate_specific': True,
        'method': 'central_finite_difference_qe_forces',
        'displacement_A': float(manifest['displacement_A']),
        'active_indices': active,
        'force_manifest_sha256': manifest['manifest_sha256'],
    })
    (root / 'transition_state_frequency.json').write_text(
        json.dumps(result, indent=2, sort_keys=True))
    return result


def run_methane_neb(calc_dir: str | Path, prefix: str,
                    n_images: int = 7, timeout_s: int = 86400,
                    execution: QEExecutionConfig | None = None) -> dict:
    """Start NEB only after both candidate-specific endpoints converge."""
    root = Path(calc_dir)
    status = methane_neb_status(root)
    if not status['endpoints_converged']:
        return status
    if status['neb'].get('converged'):
        return status
    initial = relaxed_structure(str(root / 'initial.relax.out'))
    final = relaxed_structure(str(root / 'final.relax.out'))
    if initial.get_chemical_symbols() != final.get_chemical_symbols():
        raise RuntimeError('relaxed NEB endpoints have different atom ordering')
    images = [initial] + [initial.copy() for _ in range(n_images - 2)] + [final]
    NEB(images, method='improvedtangent').interpolate(method='idpp')
    input_path = root / 'candidate.neb.in'
    output_path = root / 'candidate.neb.out'
    write_qe_neb_input(images, str(input_path), prefix)
    run_neb(str(input_path), str(output_path), timeout_s=timeout_s,
            execution=execution)
    return methane_neb_status(root)


def advance_methane_neb(calc_dir: str | Path, prefix: str,
                        n_images: int = 7, timeout_s: int = 86400,
                        execution: QEExecutionConfig | None = None,
                        restart_incomplete: bool = False) -> dict:
    """Advance endpoints, NEB, and prepared frequency jobs in safe order."""
    root = Path(calc_dir)
    for endpoint in ('initial', 'final'):
        input_path = root / f'{endpoint}.relax.in'
        output_path = root / f'{endpoint}.relax.out'
        state = qe_output_status(output_path)
        if state == 'converged':
            continue
        if not input_path.is_file():
            return methane_neb_status(root)
        if state == 'incomplete' and not restart_incomplete:
            return methane_neb_status(root)
        outcome = run_pw(str(input_path), str(output_path), timeout_s=timeout_s,
                         execution=_pw_execution(execution))
        if not outcome['converged']:
            return methane_neb_status(root)
    status = run_methane_neb(
        root, prefix, n_images=n_images, timeout_s=timeout_s,
        execution=execution)
    if status['neb'].get('converged'):
        run_frequency_sequence(
            root, timeout_s=timeout_s, execution=execution,
            restart_incomplete=restart_incomplete)
    return methane_neb_status(root)


def prepare_pyrolysis_campaign(manifest_path: str | Path) -> dict:
    """Prepare every explicitly supplied methane elementary-step endpoint pair.

    The JSON manifest maps step names to ASE-readable initial/final structures.
    An optional transition_state structure prepares finite-difference jobs.
    No chemically generic endpoint is invented for an arbitrary catalyst.
    """
    source = Path(manifest_path)
    manifest = json.loads(source.read_text())
    root_value = Path(manifest['campaign_dir'])
    root = root_value if root_value.is_absolute() else source.parent / root_value
    root = root.resolve()
    catalyst_id = str(manifest.get('candidate_id', 'candidate'))
    prepared = {}
    for step_name, kinetics_field in PYROLYSIS_ELEMENTARY_STEPS.items():
        spec = manifest.get('steps', {}).get(step_name)
        if not spec:
            prepared[step_name] = {
                'status': 'unresolved', 'kinetics_field': kinetics_field,
                'reason': 'explicit_candidate_geometries_missing'}
            continue
        step_root = root / step_name
        initial = ase_read(str((source.parent / spec['initial']).resolve()))
        final = ase_read(str((source.parent / spec['final']).resolve()))
        prepare_methane_neb(step_root, f'{catalyst_id}_{step_name}', initial, final)
        if spec.get('transition_state'):
            transition_state = ase_read(
                str((source.parent / spec['transition_state']).resolve()))
            prepare_frequency_jobs(
                step_root, transition_state, f'{catalyst_id}_{step_name}_ts',
                displacement_A=float(spec.get('displacement_A', 0.01)),
                active_indices=spec.get('frequency_active_indices'))
        prepared[step_name] = {
            'status': 'prepared', 'kinetics_field': kinetics_field,
            **methane_neb_status(step_root)}
    return {'candidate_id': catalyst_id, 'campaign_dir': str(root),
            'steps': prepared}


def pyrolysis_campaign_status(manifest_path: str | Path) -> dict:
    """Collect barriers only from fully converged, frequency-validated steps."""
    source = Path(manifest_path)
    manifest = json.loads(source.read_text())
    root_value = Path(manifest['campaign_dir'])
    root = root_value if root_value.is_absolute() else source.parent / root_value
    root = root.resolve()
    steps, resolved = {}, {}
    for step_name, kinetics_field in PYROLYSIS_ELEMENTARY_STEPS.items():
        if step_name not in manifest.get('steps', {}):
            steps[step_name] = {'status': 'unresolved',
                                'reason': 'not_supplied'}
            continue
        status = methane_neb_status(root / step_name)
        steps[step_name] = status
        if status['complete']:
            barrier = status['neb'].get('forward_barrier_eV')
            if barrier is not None and np.isfinite(float(barrier)):
                resolved[kinetics_field] = float(barrier)
    return {
        'candidate_id': str(manifest.get('candidate_id', 'candidate')),
        'steps': steps,
        'resolved_kinetics_eV': resolved,
        'unresolved_kinetics_fields': sorted(
            set(PYROLYSIS_ELEMENTARY_STEPS.values()) - set(resolved)),
        'complete': len(resolved) == len(PYROLYSIS_ELEMENTARY_STEPS),
        'evidence_level': 'converged_dft_neb_frequency' if
                          len(resolved) == len(PYROLYSIS_ELEMENTARY_STEPS)
                          else 'incomplete',
    }


def advance_pyrolysis_campaign(manifest_path: str | Path,
                                timeout_s: int = 86400,
                                execution: QEExecutionConfig | None = None,
                                restart_incomplete: bool = False) -> dict:
    """Advance all supplied elementary steps and retain unresolved fields."""
    source = Path(manifest_path)
    manifest = json.loads(source.read_text())
    root_value = Path(manifest['campaign_dir'])
    root = root_value if root_value.is_absolute() else source.parent / root_value
    root = root.resolve()
    catalyst_id = str(manifest.get('candidate_id', 'candidate'))
    for step_name in PYROLYSIS_ELEMENTARY_STEPS:
        if step_name not in manifest.get('steps', {}):
            continue
        advance_methane_neb(
            root / step_name, f'{catalyst_id}_{step_name}',
            n_images=int(manifest['steps'][step_name].get('n_images', 7)),
            timeout_s=timeout_s, execution=execution,
            restart_incomplete=restart_incomplete)
    result = pyrolysis_campaign_status(source)
    output = root / 'pyrolysis_validation.json'
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True))
    return result
