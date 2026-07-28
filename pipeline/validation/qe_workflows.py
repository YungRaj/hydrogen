"""Validated SSSP selection and candidate-specific QE NEB workflows."""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import subprocess
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.mep import NEB
from ase.io import read as ase_read

from pipeline.common.utils import BASE_DIR

SSSP_DIR = BASE_DIR / 'quantum_espresso/sssp/1.3.0-pbe-efficiency'
SSSP_MANIFEST = SSSP_DIR / 'SSSP_1.3.0_PBE_efficiency.json'


@dataclass(frozen=True)
class QEExecutionConfig:
    """Explicit, provenance-recorded local QE resource allocation."""

    mpi_ranks: int = 1
    omp_threads: int = 1
    kpoint_pools: int = 1
    image_groups: int = 1

    def validate(self, *, neb: bool = False) -> None:
        values = (self.mpi_ranks, self.omp_threads, self.kpoint_pools,
                  self.image_groups)
        if any(int(value) < 1 for value in values):
            raise ValueError('QE parallel dimensions must be positive')
        if self.mpi_ranks % self.kpoint_pools:
            raise ValueError('MPI ranks must be divisible by k-point pools')
        if neb and self.mpi_ranks % self.image_groups:
            raise ValueError('MPI ranks must be divisible by NEB image groups')
        if not neb and self.image_groups != 1:
            raise ValueError('image groups are valid only for NEB calculations')

    @classmethod
    def production_default(cls) -> 'QEExecutionConfig':
        return cls(
            mpi_ranks=int(os.environ.get('QE_MPI_RANKS', '4')),
            omp_threads=int(os.environ.get('QE_OMP_THREADS', '1')),
            kpoint_pools=int(os.environ.get('QE_KPOINT_POOLS', '1')),
        )


def build_qe_command(executable: str, input_path: str,
                     config: QEExecutionConfig, *, neb: bool = False) -> list[str]:
    """Build a shell-free MPI/QE command with validated parallel dimensions."""
    config.validate(neb=neb)
    requested = Path(executable)
    resolved = shutil.which(str(executable))
    fallback = Path('/home/ilhanraja/miniconda3/envs/qe-env/bin') / requested.name
    if resolved:
        executable_path = str(Path(resolved).resolve())
    elif requested.is_file():
        executable_path = str(requested.resolve())
    elif fallback.is_file():
        executable_path = str(fallback.resolve())
    else:
        raise RuntimeError(f'QE executable is unavailable: {executable}')
    command = []
    if config.mpi_ranks > 1:
        sibling_mpirun = Path(executable_path).parent / 'mpirun'
        mpirun = str(sibling_mpirun) if sibling_mpirun.is_file() else shutil.which('mpirun')
        if not mpirun:
            raise RuntimeError('MPI ranks requested but mpirun is unavailable')
        command.extend([mpirun, '-np', str(config.mpi_ranks)])
    command.append(executable_path)
    if config.kpoint_pools > 1:
        command.extend(['-nk', str(config.kpoint_pools)])
    if neb and config.image_groups > 1:
        command.extend(['-ni', str(config.image_groups)])
    command.extend(['-inp' if neb else '-in', str(Path(input_path).resolve())])
    return command


def _execution_record(input_path: str, output_path: str, command: list[str],
                      config: QEExecutionConfig, elapsed_s: float,
                      returncode: int, timed_out: bool) -> None:
    source = Path(input_path)
    record = {
        'command': command,
        'configuration': asdict(config),
        'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'elapsed_s': elapsed_s,
        'returncode': returncode,
        'timed_out': timed_out,
    }
    Path(f'{output_path}.execution.json').write_text(
        json.dumps(record, indent=2, sort_keys=True))


def verify_sssp(elements, directory=SSSP_DIR, manifest=SSSP_MANIFEST) -> dict:
    metadata = json.loads(Path(manifest).read_text())
    records, errors = {}, []
    for element in sorted(set(elements)):
        entry = metadata.get(element)
        if not entry:
            errors.append(f'{element}:missing_manifest')
            continue
        path = Path(directory) / entry['filename']
        if not path.is_file():
            errors.append(f'{element}:missing_file:{entry["filename"]}')
            continue
        digest = hashlib.md5(path.read_bytes()).hexdigest()
        if digest != entry['md5']:
            errors.append(f'{element}:checksum_mismatch')
            continue
        records[element] = {'filename': entry['filename'], 'md5': digest,
                            'cutoff_wfc_Ry': float(entry['cutoff_wfc']),
                            'cutoff_rho_Ry': float(entry['cutoff_rho'])}
    return {'valid': not errors, 'records': records, 'errors': errors,
            'ecutwfc_Ry': max((x['cutoff_wfc_Ry'] for x in records.values()), default=None),
            'ecutrho_Ry': max((x['cutoff_rho_Ry'] for x in records.values()), default=None)}


def methane_dissociation_images(slab: Atoms, active_index: int,
                                 n_images: int = 7) -> list[Atoms]:
    """Build candidate-specific CH4(g)+* -> CH3*+H* NEB endpoints/images."""
    if n_images < 5 or active_index < 0 or active_index >= len(slab):
        raise ValueError('invalid NEB image count or active site')
    site = slab.positions[active_index].copy()
    initial = slab.copy()
    ch4 = molecule('CH4')
    ch4.translate(site + np.array([0.0, 0.0, 3.2]) - ch4.get_center_of_mass())
    initial += ch4

    final = slab.copy()
    ch3 = molecule('CH3')
    ch3.translate(site + np.array([0.0, 0.0, 2.0]) - ch3.get_center_of_mass())
    final += ch3
    final += Atoms('H', positions=[site + np.array([1.5, 0.0, 1.2])])
    initial.set_cell(slab.cell); final.set_cell(slab.cell)
    initial.set_pbc(slab.pbc); final.set_pbc(slab.pbc)
    images = [initial] + [initial.copy() for _ in range(n_images - 2)] + [final]
    NEB(images, method='improvedtangent').interpolate(method='idpp')
    return images


def write_qe_neb_input(images: list[Atoms], path: str, prefix: str) -> dict:
    """Write a climbing-image neb.x input using one verified SSSP family."""
    if not images or any(len(x) != len(images[0]) for x in images):
        raise ValueError('NEB images must have identical atom ordering')
    elements = sorted(set(images[0].get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    records = verified['records']
    species = '\n'.join(f" {e} 1.0 {records[e]['filename']}" for e in elements)
    cell = '\n'.join(' '.join(f'{v:.12f}' for v in row) for row in images[0].cell.array)
    positions = []
    for number, image in enumerate(images):
        tag = 'FIRST_IMAGE' if number == 0 else ('LAST_IMAGE' if number == len(images)-1 else 'INTERMEDIATE_IMAGE')
        positions.append(f'{tag}\nATOMIC_POSITIONS angstrom\n' + '\n'.join(
            f'{a.symbol} {a.x:.12f} {a.y:.12f} {a.z:.12f}' for a in image))
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    magnetization = '\n'.join(
        f" starting_magnetization({i})={0.5 if element in magnetic else 0.05}"
        for i, element in enumerate(elements, 1))
    text = f"""BEGIN
BEGIN_PATH_INPUT
&PATH
 num_of_images={len(images)}, opt_scheme='broyden', CI_scheme='auto',
 path_thr=0.05, ds=1.0,
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
&CONTROL
 calculation='scf', prefix='{prefix}', pseudo_dir='{SSSP_DIR}', outdir='./tmp'
/
&SYSTEM
 ibrav=0, nat={len(images[0])}, ntyp={len(elements)},
 ecutwfc={verified['ecutwfc_Ry']}, ecutrho={verified['ecutrho_Ry']},
 occupations='smearing', smearing='mv', degauss=0.02, nspin=2
{magnetization}
/
&ELECTRONS
 conv_thr=1.0d-8, mixing_beta=0.3
/
ATOMIC_SPECIES
{species}
CELL_PARAMETERS angstrom
{cell}
K_POINTS automatic
 2 2 1 0 0 0
BEGIN_POSITIONS
{chr(10).join(positions)}
END_POSITIONS
END_ENGINE_INPUT
END
"""
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return {'path': str(target), 'n_images': len(images), 'sssp': verified}


def run_neb(input_path: str, output_path: str, timeout_s: int = 86400,
            execution: QEExecutionConfig | None = None) -> dict:
    neb = shutil.which('neb.x') or '/home/ilhanraja/miniconda3/envs/qe-env/bin/neb.x'
    workdir = Path(input_path).resolve().parent
    (workdir / 'tmp').mkdir(exist_ok=True)
    config = execution or QEExecutionConfig()
    command = build_qe_command(neb, input_path, config, neb=True)
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = str(config.omp_threads)
    started = time.monotonic()
    timed_out, returncode = False, -1
    with open(output_path, 'w') as sink:
        try:
            proc = subprocess.run(
                command, cwd=str(workdir), env=environment, stdout=sink,
                stderr=subprocess.STDOUT, timeout=timeout_s)
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    elapsed = time.monotonic() - started
    _execution_record(input_path, output_path, command, config, elapsed,
                      returncode, timed_out)
    text = Path(output_path).read_text(errors='replace')
    lower = text.lower()
    converged = returncode == 0 and 'job done' in lower and 'error' not in lower
    return {'converged': converged, 'returncode': returncode,
            'output': str(output_path), 'candidate_specific': True,
            'execution': asdict(config), 'timed_out': timed_out}


def write_qe_relax_input(atoms: Atoms, path: str, prefix: str,
                         kpoints=(2, 2, 1)) -> dict:
    """Write a spin-polarized, force-converged endpoint relaxation."""
    elements = sorted(set(atoms.get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    species = '\n'.join(f"{e} 1.0 {verified['records'][e]['filename']}" for e in elements)
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    mags = '\n'.join(f" starting_magnetization({i})={0.5 if e in magnetic else 0.05}"
                     for i, e in enumerate(elements, 1))
    cell = '\n'.join(' '.join(f'{v:.12f}' for v in row) for row in atoms.cell.array)
    positions = '\n'.join(f'{a.symbol} {a.x:.12f} {a.y:.12f} {a.z:.12f}' for a in atoms)
    text = f"""&CONTROL
 calculation='relax', prefix='{prefix}', pseudo_dir='{SSSP_DIR}', outdir='./tmp',
 forc_conv_thr=1.0d-3, nstep=100, tprnfor=.true.
/
&SYSTEM
 ibrav=0, nat={len(atoms)}, ntyp={len(elements)}, ecutwfc={verified['ecutwfc_Ry']},
 ecutrho={verified['ecutrho_Ry']}, occupations='smearing', smearing='mv', degauss=0.02, nspin=2
{mags}
/
&ELECTRONS
 conv_thr=1.0d-8, mixing_beta=0.3
/
&IONS
 ion_dynamics='bfgs'
/
ATOMIC_SPECIES
{species}
ATOMIC_POSITIONS angstrom
{positions}
CELL_PARAMETERS angstrom
{cell}
K_POINTS automatic
 {kpoints[0]} {kpoints[1]} {kpoints[2]} 0 0 0
"""
    target = Path(path); target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return {'path': str(target), 'sssp': verified}


def run_pw(input_path: str, output_path: str, timeout_s: int = 86400,
           execution: QEExecutionConfig | None = None) -> dict:
    workdir = Path(input_path).resolve().parent
    (workdir / 'tmp').mkdir(exist_ok=True)
    pw = shutil.which('pw.x') or '/home/ilhanraja/miniconda3/envs/qe-env/bin/pw.x'
    config = execution or QEExecutionConfig.production_default()
    command = build_qe_command(pw, input_path, config)
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = str(config.omp_threads)
    started = time.monotonic()
    timed_out, returncode = False, -1
    with open(output_path, 'w') as sink:
        try:
            proc = subprocess.run(
                command, stdout=sink, stderr=subprocess.STDOUT, cwd=str(workdir),
                env=environment, timeout=timeout_s)
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    elapsed = time.monotonic() - started
    _execution_record(input_path, output_path, command, config, elapsed,
                      returncode, timed_out)
    text = Path(output_path).read_text(errors='replace')
    converged = returncode == 0 and 'JOB DONE' in text and 'convergence NOT achieved' not in text
    return {'converged': converged, 'returncode': returncode,
            'output': output_path, 'execution': asdict(config),
            'timed_out': timed_out}


def relaxed_structure(output_path: str) -> Atoms:
    """Read the last geometry only from a cleanly completed QE relaxation."""
    text = Path(output_path).read_text(errors='replace')
    if 'JOB DONE' not in text or 'convergence NOT achieved' in text:
        raise RuntimeError('endpoint relaxation is not converged')
    return ase_read(output_path, format='espresso-out', index=-1)


def parse_neb_result(output_path: str) -> dict:
    text = Path(output_path).read_text(errors='replace')
    energies = [float(x) for x in re.findall(r'activation energy \(->\)\s*=\s*([-+0-9.Ee]+)', text)]
    reverse = [float(x) for x in re.findall(r'activation energy \(<-\)\s*=\s*([-+0-9.Ee]+)', text)]
    errors = [float(x) for x in re.findall(r'path length\s*=\s*([-+0-9.Ee]+)', text)]
    converged = 'JOB DONE' in text and 'neb: convergence achieved' in text.lower()
    return {'converged': converged,
            'forward_barrier_eV': energies[-1] if energies else None,
            'reverse_barrier_eV': reverse[-1] if reverse else None,
            'path_metric': errors[-1] if errors else None,
            'candidate_specific': True}


def partial_hessian(forces_plus: np.ndarray, forces_minus: np.ndarray,
                    displacement_A: float, masses_amu: np.ndarray) -> dict:
    """Construct a mass-weighted partial Hessian from central force differences."""
    plus = np.asarray(forces_plus, float); minus = np.asarray(forces_minus, float)
    if plus.shape != minus.shape or plus.ndim != 3 or displacement_A <= 0:
        raise ValueError('forces must be (3N, N, 3) central-difference arrays')
    n_atoms = plus.shape[1]
    if plus.shape[0] != 3 * n_atoms or len(masses_amu) != n_atoms:
        raise ValueError('partial Hessian dimensions are inconsistent')
    hessian = -(plus - minus).reshape(3*n_atoms, 3*n_atoms).T / (2 * displacement_A)
    hessian = 0.5 * (hessian + hessian.T)
    weights = np.repeat(np.sqrt(np.asarray(masses_amu, float)), 3)
    eigvals, eigvecs = np.linalg.eigh(hessian / np.outer(weights, weights))
    # 1 eV/A^2/amu -> (521.47083 cm^-1)^2
    frequencies = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) * 521.47083
    imaginary = frequencies[frequencies < -50.0]
    return {'frequencies_cm1': frequencies.tolist(),
            'imaginary_count': int(len(imaginary)),
            'valid_transition_state': len(imaginary) == 1,
            'mode_vectors': eigvecs.tolist()}
