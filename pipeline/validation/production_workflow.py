"""Resumable, fail-closed orchestration for production QE validation."""

from __future__ import annotations

from pathlib import Path

from ase.mep import NEB

from pipeline.validation.dft_validator import parse_convergence
from pipeline.validation.qe_workflows import (
    parse_neb_result,
    relaxed_structure,
    run_neb,
    run_pw,
    write_qe_neb_input,
)


ORR_STAGES = ('clean', 'OH', 'O', 'OOH', 'h2', 'h2o')


def qe_output_status(path: str | Path) -> str:
    """Classify an output without mistaking a partial energy for evidence."""
    target = Path(path)
    if not target.exists() or target.stat().st_size == 0:
        return 'missing'
    text = target.read_text(errors='replace').lower()
    if 'error in routine' in text or 'convergence not achieved' in text:
        return 'failed'
    if parse_convergence(str(target)):
        return 'converged'
    return 'incomplete'


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
                     timeout_s: int = 86400, restart_incomplete: bool = False) -> dict:
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
        outcome = run_pw(str(input_path), str(output_path), timeout_s=timeout_s)
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
    return {
        'endpoints': endpoints,
        'endpoints_converged': all(x == 'converged' for x in endpoints.values()),
        'neb': neb,
        'frequency': frequency,
        'complete': bool(neb.get('converged') and
                         frequency.get('valid_transition_state')),
    }


def run_methane_neb(calc_dir: str | Path, prefix: str,
                    n_images: int = 7, timeout_s: int = 86400) -> dict:
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
    run_neb(str(input_path), str(output_path), timeout_s=timeout_s)
    return methane_neb_status(root)
