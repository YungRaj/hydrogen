"""Validated SSSP selection and candidate-specific QE NEB workflows."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import re
import time
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.geometry import find_mic
from ase.mep import NEB
from ase.io import read as ase_read
from ase.units import Bohr

from pipeline.data_models.quantum import (
    NEBResult,
    SolverRunResult,
    TransitionStateFrequencyResult,
)
from pipeline.utils import BASE_DIR
from pipeline.simulation.executables import resolve_executable, resolve_qe_executable

SSSP_DIR = BASE_DIR / 'quantum_espresso/sssp/1.3.0-pbe-efficiency'
SSSP_MANIFEST = SSSP_DIR / 'SSSP_1.3.0_PBE_efficiency.json'


def _fixed_atom_indices(atoms: Atoms) -> set[int]:
    """Return ASE-constrained atoms that QE must hold fixed in x/y/z."""
    fixed = set()
    for constraint in atoms.constraints:
        if constraint.__class__.__name__ != 'FixAtoms':
            raise ValueError(
                f'unsupported QE relaxation constraint: '
                f'{constraint.__class__.__name__}'
            )
        getter = getattr(constraint, 'get_indices', None)
        if getter is not None:
            fixed.update(int(index) for index in getter())
    return fixed


def _qe_positions(
    atoms: Atoms, include_constraints: bool = True, freeze_all: bool = False
) -> str:
    """Render QE positions with explicit, auditable Cartesian move flags."""
    fixed = (
        set(range(len(atoms)))
        if freeze_all
        else _fixed_atom_indices(atoms) if include_constraints else set()
    )
    rows = []
    for index, atom in enumerate(atoms):
        flags = ' 0 0 0' if index in fixed else ' 1 1 1'
        rows.append(f'{atom.symbol} {atom.x:.12f} {atom.y:.12f} {atom.z:.12f}{flags}')
    return '\n'.join(rows)


@dataclass(frozen=True)
class QEExecutionConfig:
    """Explicit, provenance-recorded local QE resource allocation."""

    mpi_ranks: int = 1
    omp_threads: int = 1
    kpoint_pools: int = 1
    image_groups: int = 1
    # This is provenance, not a caller-selectable policy switch. Production QE
    # is GPU-only and no configuration may weaken that invariant.
    require_gpu: bool = field(default=True, init=False)

    def validate(self, *, neb: bool = False) -> None:
        """Reject invalid configuration before it reaches scientific execution.

        Args:
            neb: Whether to build a NEB rather than pw.x execution command.
        """
        values = (
            self.mpi_ranks,
            self.omp_threads,
            self.kpoint_pools,
            self.image_groups,
        )
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
        """Perform production default.

        Returns:
            A `'QEExecutionConfig'` containing the production default result.
        """
        return cls(
            mpi_ranks=int(os.environ.get('QE_MPI_RANKS', '4')),
            omp_threads=int(os.environ.get('QE_OMP_THREADS', '1')),
            kpoint_pools=int(os.environ.get('QE_KPOINT_POOLS', '1')),
        )


@dataclass(frozen=True)
class NEBPathConfig:
    """Numerical controls for one fidelity tier of an NEB calculation."""

    nstep_path: int = 100
    path_thr_eV_A: float = 0.05
    optimizer: str = 'broyden'
    climbing_scheme: str = 'auto'
    step_size: float = 1.0
    spring_min: float = 0.1
    spring_max: float = 0.1
    electron_maxstep: int = 300
    scf_conv_thr_Ry: float = 1.0e-8
    mixing_beta: float = 0.2
    mixing_ndim: int = 8
    starting_potential: str = 'atomic'
    starting_wavefunctions: str = 'atomic+random'
    diagonalization: str = 'david'
    total_magnetization: float | None = None
    freeze_all_atoms: bool = False
    # QE otherwise interpolates periodic atoms through their literal Cartesian
    # replicas.  Surface paths can then acquire artificial lattice-spanning
    # jumps and atom collisions during mandatory equal-arc initialization.
    minimum_image: bool = True
    # Precondition from-scratch inputs so QE 7.5's one-segment interpolation
    # cursor cannot extrapolate across clusters of short source segments.
    precondition_equal_arc: bool = True
    kpoints: tuple[int, int, int] = (2, 2, 1)
    restart_mode: str = 'from_scratch'

    @classmethod
    def coarse_preconditioner(cls) -> 'NEBPathConfig':
        """Return stable non-climbing controls for a long initial path.

        Returns:
            Coarse path-preconditioning controls.
        """
        return cls(
            nstep_path=150,
            path_thr_eV_A=0.15,
            optimizer='quick-min',
            climbing_scheme='no-CI',
            step_size=0.2,
            # Keep spring physics independent of the optimizer trust step.
            # Scaling k as 1/ds**2 made modest spacing differences dominate
            # the Ni(111) band: the same geometry reported 8.317 eV/A with
            # 9.25--15.42-a.u. springs but 1.212 eV/A with the conventional
            # 0.1-a.u. spring, without changing its energy profile.
            spring_min=0.1,
            spring_max=0.1,
            electron_maxstep=500,
            scf_conv_thr_Ry=1.0e-6,
            # Metallic, spin-polarized Ni images are first preconverged as
            # independent SCFs.  Reuse those charge densities here rather
            # than allowing eleven coupled images to start from atomic
            # superpositions and drift into unrelated magnetic solutions.
            # Independent image SCFs supply stable 3x3x1 densities and
            # wavefunctions.  With those seeds, the conventional Davidson
            # solver and conservative 0.02 mixing preserve the independently
            # validated magnetic branch in the coupled path.  A dense Ni(111)
            # trial with 0.1 mixing drove image 7 from about 8.95 to 4.34 Bohr
            # magnetons during its first coupled SCF, despite a clean seed;
            # the same image is stable with 0.02 mixing.
            mixing_beta=0.02,
            starting_potential='file',
            starting_wavefunctions='file',
            diagonalization='david',
            kpoints=(3, 3, 1),
        )

    @classmethod
    def coarse_fresh_image_seed(cls) -> 'NEBPathConfig':
        """Initialize isolated coarse images without prior QE save data.

        This is an electronic preconditioning tier, not an NEB path result.
        Conservative local-TF mixing is used because stale or geometry-
        incompatible metallic wavefunctions can send Ni images onto a
        catastrophically discontinuous SCF branch.  Each image must converge
        independently before its save directory is admitted to NEB.

        Returns:
            Conservative controls for an independently seeded image SCF.
        """
        return replace(
            cls.coarse_preconditioner(),
            mixing_beta=0.01,
            starting_potential='atomic',
            starting_wavefunctions='atomic+random',
        )

    @classmethod
    def coarse_fixed_spin_image_seed(cls) -> 'NEBPathConfig':
        """Rescue a poor Ni state on the eight-Bohr-magneton branch.

        The Ni(111) image-7 diagnostic established that retaining the
        geometry-compatible charge density while rebuilding the
        wavefunctions converges this metallic state.  Starting from atomic
        charge instead produced large magnetic oscillations, even with CG
        diagonalization and very small mixing.  The spin constraint is only
        an electronic preconditioner; :meth:`coarse_released_spin_image_seed`
        removes it before an image is admitted to the NEB calculation.

        Returns:
            Fixed-spin electronic rescue controls for a Ni path image.
        """
        return replace(
            cls.coarse_fresh_image_seed(),
            mixing_beta=0.01,
            mixing_ndim=24,
            starting_potential='file',
            starting_wavefunctions='atomic+random',
            diagonalization='david',
            total_magnetization=8.0,
        )

    @classmethod
    def coarse_released_spin_image_seed(cls) -> 'NEBPathConfig':
        """Release a fixed-spin seed and reconverge before NEB admission.

        Once the fixed-spin stage has placed the image on the physical
        metallic branch, reuse both density and wavefunctions and restore the
        regular coarse-path mixing.  Retaining the rescue-stage 0.01 damping
        was observed to stall an otherwise stable Ni endpoint a few micro-Ry
        above the admission threshold.

        Returns:
            Unconstrained controls for validating a rescued image state.
        """
        return replace(
            cls.coarse_preconditioner(),
            starting_potential='file',
            starting_wavefunctions='file',
            total_magnetization=None,
        )

    @classmethod
    def coarse_fixed_spin_coupled_seed(cls) -> 'NEBPathConfig':
        """Precondition all coupled images without moving the NEB geometry.

        Standalone image convergence does not by itself guarantee that QE's
        coupled ``from_scratch`` initialization retains the same magnetic
        branch.  The dense Ni(111) validation reproduced a collapse of image
        7 from about 8.95 to below 1 Bohr magneton even after every standalone
        seed had converged.  This non-evidentiary pass constrains the total
        cell magnetization while evaluating every image.  Its deliberately
        loose path threshold makes QE stop after electronic initialization,
        before an optimizer displacement can turn constrained-spin forces
        into a candidate reaction path.

        That constrained state must be followed by an unconstrained coupled
        calculation; neither its forces nor its apparent barrier are valid
        scientific evidence.

        Returns:
            Fixed-spin, no-motion controls for coupled electronic seeding.
        """
        return replace(
            cls.coarse_preconditioner(),
            nstep_path=1,
            path_thr_eV_A=1.0e6,
            total_magnetization=8.0,
            # A fixed total magnetization can increase the required band
            # count.  Reusing unconstrained wavefunctions then makes QE abort
            # in read_collected_wfc; retain the validated density but rebuild
            # wavefunctions for the constrained occupation manifold.
            starting_wavefunctions='atomic+random',
            freeze_all_atoms=True,
        )

    @classmethod
    def coarse_released_spin_coupled_seed(cls) -> 'NEBPathConfig':
        """Release a fixed-spin coupled path without moving its geometry.

        The fixed-spin coupled pass is only an electronic rescue operation.
        Before any NEB force or energy is treated as evidence, every image
        must reconverge without the total-magnetization constraint.  This
        second electronic-only pass retains the rescued charge densities but
        rebuilds wavefunctions because releasing the constraint can reduce
        the number of occupied bands.  All atoms remain frozen and the loose
        path threshold prevents an optimizer displacement.

        Returns:
            Unconstrained, no-motion controls for coupled state validation.
        """
        return replace(
            cls.coarse_fixed_spin_coupled_seed(),
            total_magnetization=None,
            starting_potential='file',
            starting_wavefunctions='atomic+random',
        )

    @classmethod
    def coarse_neighbor_density_image_seed(cls) -> 'NEBPathConfig':
        """Seed one damaged image from a converged neighboring density.

        Adjacent NEB images share cell, composition, and a nearby geometry, so
        their real-space density is a useful recovery seed.  Wavefunctions
        are geometry-specific and are deliberately rebuilt.  Unlike the
        fixed-spin rescue, this tier leaves magnetization unconstrained; that
        avoids forcing an otherwise healthy image away from its natural Ni
        magnetic branch.

        Returns:
            Neighbor-density recovery controls for one path image.
        """
        return replace(
            cls.coarse_preconditioner(),
            starting_potential='file',
            starting_wavefunctions='atomic+random',
            total_magnetization=None,
        )

    @classmethod
    def coarse_refiner(cls) -> 'NEBPathConfig':
        """Resume a quick-min checkpoint with stable quasi-Newton updates.

        Returns:
            Restart controls using the Broyden path optimizer.
        """
        return replace(
            cls.coarse_preconditioner(), optimizer='broyden', restart_mode='restart'
        )

    @classmethod
    def coarse_finisher(cls) -> 'NEBPathConfig':
        """Damp an oscillatory Broyden path before climbing-image refinement.

        Returns:
            Reduced-step Broyden controls with unchanged physical springs.
        """
        return replace(cls.coarse_refiner(), step_size=0.1)

    @classmethod
    def coarse_stabilizer(cls) -> 'NEBPathConfig':
        """Resume a plateaued coarse path with a smaller trust step.

        This tier is deliberately still non-climbing and keeps the inexpensive
        electronic threshold.  It is a path-geometry stabilization pass, not
        evidence for a final transition state or activation barrier.

        Returns:
            Strongly damped non-climbing stabilization controls.
        """
        return replace(cls.coarse_finisher(), step_size=0.05)

    @classmethod
    def coarse_quick_min_history_reset_stabilizer(cls) -> 'NEBPathConfig':
        """Restart an oscillatory coarse quick-min path at half step size.

        The completed ``*.crd`` coordinates must be extracted and written as
        a new ``from_scratch`` path so QE cannot retain the alternating
        quick-min velocity. This profile is valid only when restarting an
        existing quick-min trajectory. It must not be used to hand an SD path
        to quick-min: the Ni(111) reference showed that QE normalizes the first
        SD-to-quick-min move, and reducing ``ds`` from 0.1 to 0.01 left its
        0.210-angstrom maximum displacement unchanged. A diagnostic that
        simultaneously increased the springs by four raised the initial path
        error from 21.8 to 209.1 eV/angstrom, so preserving the spring bounds
        is necessary to isolate optimizer damping. Electronic fidelity and
        the non-climbing coarse force target remain unchanged.

        Returns:
            Fresh-history quick-min controls for an existing quick-min path.
        """
        return replace(
            cls.coarse_preconditioner(),
            nstep_path=250,
            step_size=0.1,
            restart_mode='from_scratch',
            # The interrupted next trial may leave geometry-specific
            # wavefunctions that do not correspond to the last completed
            # ``*.crd`` path. Retain charge-density seeding but rebuild WFCs.
            starting_wavefunctions='atomic+random',
        )

    @classmethod
    def coarse_quick_min_to_sd_stabilizer(cls) -> 'NEBPathConfig':
        """Resume an accepted quick-min checkpoint without its velocity.

        Quick-min can continue improving the lower envelope while its
        retained velocity alternately overshoots the path force.  A native
        restart preserves the exact accepted coordinates (and therefore
        avoids QE's from-scratch spline redistribution), while switching to
        steepest descent makes the next displacement depend only on the
        current NEB force.  Keep the accepted springs and half step unchanged
        so the intervention resets optimizer memory rather than path physics.

        Returns:
            Non-climbing, history-free controls for an exact coarse restart.
        """
        return replace(
            cls.coarse_preconditioner(),
            nstep_path=250,
            optimizer='sd',
            step_size=0.1,
            restart_mode='restart',
        )

    @classmethod
    def coarse_sd_to_broyden_refiner(cls) -> 'NEBPathConfig':
        """Learn curvature after a stable small-step SD smoothing window.

        This handoff preserves the accepted native checkpoint, half step,
        spring model, electronic controls, and sampling.  With no existing
        ``*.broyden`` file, QE's first Broyden displacement has the same
        ``ds**2`` diagonal scaling as steepest descent; later displacements
        can use measured cross-image curvature and are capped internally.
        The tier is appropriate only for a bounded, monitored trial after SD
        improvement has materially slowed.

        Returns:
            Small-step Broyden controls for an exact coarse restart.
        """
        return replace(cls.coarse_quick_min_to_sd_stabilizer(), optimizer='broyden')

    @classmethod
    def coarse_broyden_to_sd_stabilizer(cls) -> 'NEBPathConfig':
        """Retain a productive Broyden geometry but discard its history.

        QE writes ``pathN`` *after* evaluating iteration N and applying that
        iteration's optimizer displacement.  Consequently the geometry whose
        forces were reported at iteration N is held in ``path(N-1)``.  A
        bounded curvature-learning window can discover a lower-force path and
        overshoot in its following trial; callers must restore that preceding
        checkpoint, remove ``*.broyden``, and resume with direct small-step
        forces.  Wavefunctions are rebuilt because a later rejected trial may
        already have overwritten the per-image save directories.  Charge
        density seeding, springs, path sampling, and coarse SCF tolerance are
        retained.

        Returns:
            History-free SD controls for stabilizing a Broyden-discovered path.
        """
        return replace(
            cls.coarse_quick_min_to_sd_stabilizer(),
            starting_wavefunctions='atomic+random',
        )

    @classmethod
    def coarse_alternate_refiner(cls) -> 'NEBPathConfig':
        """Reset a plateaued path with QE's alternate Broyden update.

        Returns:
            Stabilization controls using QE's alternate Broyden update.
        """
        return replace(cls.coarse_stabilizer(), optimizer='broyden2')

    @classmethod
    def coarse_force_escape(cls) -> 'NEBPathConfig':
        """Take conservative steepest-descent steps after a Broyden plateau.

        Returns:
            Restart controls for direct steepest-descent path motion.
        """
        # Direct steepest descent needs a meaningful displacement to escape a
        # quasi-Newton plateau. Keep the independently validated spring model
        # fixed so this changes optimizer motion only.
        return replace(
            cls.coarse_preconditioner(),
            nstep_path=250,
            optimizer='sd',
            step_size=1.0,
            restart_mode='restart',
        )

    @classmethod
    def coarse_force_finisher(cls) -> 'NEBPathConfig':
        """Damp the final direct-force oscillation before strict refinement.

        Returns:
            Reduced-step steepest-descent finishing controls.
        """
        return replace(cls.coarse_force_escape(), step_size=0.5)

    @classmethod
    def coarse_quasi_newton_finisher(cls) -> 'NEBPathConfig':
        """Learn cross-image curvature after direct-force path smoothing.

        Returns:
            Broyden controls initialized from the force-smoothed path.
        """
        return replace(cls.coarse_force_finisher(), optimizer='broyden')

    @classmethod
    def coarse_electronic_stabilizer(cls) -> 'NEBPathConfig':
        """Resume Broyden with damped mixing after a metallic SCF oscillation.

        Path iteration 170 of the Ni(111) reference calculation exposed a
        tail-latency failure in image 7: ``mixing_beta=0.1`` oscillated above
        the coarse electronic threshold after the other ten images had
        converged.  Independent image-7 diagnostics established that 0.02
        local-TF mixing converges this magnetic metallic state.  This tier
        changes only the electronic damping; it preserves the accepted path
        checkpoint and the productive quasi-Newton geometry controls.

        Returns:
            Quasi-Newton controls with stabilized metallic SCF mixing.
        """
        return replace(cls.coarse_quasi_newton_finisher(), mixing_beta=0.02)

    @classmethod
    def coarse_history_reset_refiner(cls) -> 'NEBPathConfig':
        """Start fresh Broyden history from an accepted checkpoint geometry.

        QE restart files retain optimizer state.  A Ni(111) reference retry
        demonstrated that changing ``ds`` on a restart replayed the rejected
        Broyden direction and that stronger springs amplified the path force.
        This configuration must therefore be used with geometries extracted
        by :func:`neb_checkpoint_images` and a newly written ``from_scratch``
        NEB input.  It preserves the spring model under which the checkpoint
        was accepted, damps the first fresh Broyden step, and retains the
        independently validated metallic-SCF stabilization.

        Returns:
            Fresh-history Broyden controls for an extracted checkpoint path.
        """
        return replace(
            cls.coarse_electronic_stabilizer(),
            restart_mode='from_scratch',
            step_size=0.25,
        )

    @classmethod
    def coarse_history_reset_force_stabilizer(cls) -> 'NEBPathConfig':
        """Replace regressing Broyden curvature with small direct-force steps.

        Use this only after consecutive path-force increases with converged
        image SCFs. Geometry is recovered from a chemically valid checkpoint,
        optimizer history is discarded, and the unchanged spring model keeps
        the intervention limited to optimization rather than path physics.

        Returns:
            Small-step direct-force controls with discarded optimizer history.
        """
        return replace(
            cls.coarse_history_reset_refiner(), optimizer='sd', step_size=0.1
        )

    @classmethod
    def coarse_history_reset_force_accelerator(cls) -> 'NEBPathConfig':
        """Scale a verified monotonic direct-force descent to a useful step.

        Returns:
            Fresh-history direct-force controls with a larger trusted step.
        """
        return replace(cls.coarse_history_reset_force_stabilizer(), step_size=0.5)

    @classmethod
    def strict_climbing_refiner(cls) -> 'NEBPathConfig':
        """Refine a converged coarse path without weakening its sampling.

        The native coarse checkpoint is retained because steepest descent has
        no quasi-Newton curvature history to contaminate CI-NEB.  Rewriting
        the same coordinates as a fresh QE path triggers internal spline
        redistribution and can destroy an already converged force profile.
        The validated 3x3x1 k-point grid and saved metallic electronic states
        are retained, while the SCF and path-force thresholds are tightened.
        The validated spring model remains independent of optimizer step size;
        the smaller half-step was observed to leave the climbing force
        effectively stationary.

        Returns:
            Conservative strict climbing-image controls for a coarse path.
        """
        return replace(
            cls.coarse_electronic_stabilizer(),
            nstep_path=500,
            path_thr_eV_A=0.05,
            optimizer='sd',
            climbing_scheme='auto',
            step_size=1.0,
            scf_conv_thr_Ry=1.0e-8,
            restart_mode='restart',
        )

    @classmethod
    def strict_climbing_finisher(cls) -> 'NEBPathConfig':
        """Damp a strict CI path after force transfers to a neighbor image.

        This changes only the optimizer displacement.  The established spring
        constants are deliberately retained: changing them changes the
        elastic-band force and therefore cannot diagnose whether a smaller
        optimization step damps an overshoot.  It is intended for native
        restart from a chemically valid completed path, so no spline
        redistribution or electronic-fidelity change is introduced.  A
        rejected Ni(111) diagnostic confirmed this distinction: following
        QE's rescaled spring suggestion increased the neighboring-image force
        from 0.45 to 1.26 eV/angstrom on the same inherited trial geometry.

        Returns:
            Strict climbing-image controls with a consistently damped step.
        """
        return replace(cls.strict_climbing_refiner(), step_size=0.5)

    @classmethod
    def strict_climbing_quick_min_stabilizer(cls) -> 'NEBPathConfig':
        """Stabilize a strict CI path after SD stalls and Broyden diverges.

        Quick-min retains only velocity components aligned with the current
        NEB force, avoiding the unbounded inverse-Hessian steps observed for
        Broyden on the Ni(111) transition region.  A conservative 0.1 step is
        used while preserving the accepted spring model, strict electronic
        threshold, k-point sampling, and native path checkpoint.

        Returns:
            Strict climbing-image controls using damped quick-min motion.
        """
        return replace(
            cls.strict_climbing_refiner(), optimizer='quick-min', step_size=0.1
        )

    def validate(self) -> None:
        """Reject incomplete or physically unsupported NEB controls.

        Raises:
            ValueError: If a numerical control, execution mode, or physical
                sampling parameter cannot define a supported QE NEB run.
        """
        if len(self.kpoints) != 3:
            raise ValueError('NEB k-point grid must have three dimensions')
        if self.optimizer not in {'quick-min', 'broyden', 'broyden2', 'sd'}:
            raise ValueError(f'unsupported NEB optimizer: {self.optimizer}')
        if self.climbing_scheme not in {'no-CI', 'auto', 'manual'}:
            raise ValueError(
                f'unsupported climbing-image scheme: {self.climbing_scheme}'
            )
        if self.restart_mode not in {'from_scratch', 'restart'}:
            raise ValueError(f'unsupported NEB restart mode: {self.restart_mode}')
        if self.starting_potential not in {'atomic', 'file'}:
            raise ValueError(
                f'unsupported starting potential: {self.starting_potential}'
            )
        if self.starting_wavefunctions not in {
            'atomic',
            'atomic+random',
            'random',
            'file',
        }:
            raise ValueError(
                'unsupported starting wavefunctions: ' f'{self.starting_wavefunctions}'
            )
        if self.diagonalization not in {'david', 'cg', 'paro', 'ppcg'}:
            raise ValueError(f'unsupported diagonalization: {self.diagonalization}')
        if self.total_magnetization is not None and not np.isfinite(
            self.total_magnetization
        ):
            raise ValueError('total magnetization must be finite when set')
        values = (
            self.nstep_path,
            self.path_thr_eV_A,
            self.step_size,
            self.spring_min,
            self.spring_max,
            self.electron_maxstep,
            self.scf_conv_thr_Ry,
            self.mixing_beta,
            self.mixing_ndim,
            *self.kpoints,
        )
        if any(value <= 0 for value in values):
            raise ValueError('NEB path controls must be positive')
        if self.spring_max < self.spring_min:
            raise ValueError('NEB maximum spring must not be below minimum')


def build_qe_command(
    executable: str, input_path: str, config: QEExecutionConfig, *, neb: bool = False
) -> list[str]:
    """Build a shell-free MPI/QE command with validated parallel dimensions.

    Args:
        executable: Executable used by this operation.
        input_path: Filesystem location used for input path.
        config: Configuration controlling this operation.
        neb: Whether to enable neb.

    Returns:
        List of computed or validated records.
    """
    config.validate(neb=neb)
    requested = Path(executable).name
    executable_path = (
        resolve_qe_executable(requested)
        if requested in {'pw.x', 'neb.x'}
        else resolve_executable(executable, required=True)
    )
    command = []
    if requested in {'pw.x', 'neb.x'}:
        mpirun = os.environ.get('MPIEXEC')
        if mpirun:
            mpirun = resolve_executable(mpirun, required=True)
        else:
            raise RuntimeError(
                'MPIEXEC is not set. Source the activate.sh generated by '
                'scripts/install_qe_gpu.sh so QE and its matching NVIDIA HPC-X '
                'MPI runtime cannot be mixed with another MPI installation.'
            )
        command.extend([mpirun, '-np', str(config.mpi_ranks)])
    command.append(executable_path)
    if config.kpoint_pools > 1:
        command.extend(['-nk', str(config.kpoint_pools)])
    if neb and config.image_groups > 1:
        command.extend(['-ni', str(config.image_groups)])
    command.extend(['-inp' if neb else '-in', str(Path(input_path).resolve())])
    return command


def _execution_record(
    input_path: str,
    output_path: str,
    command: list[str],
    config: QEExecutionConfig,
    elapsed_s: float,
    returncode: int,
    timed_out: bool,
    interrupted: bool = False,
    checkpoint_restored: bool = False,
) -> None:
    source = Path(input_path)
    record = {
        'command': command,
        'configuration': asdict(config),
        'input_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
        'elapsed_s': elapsed_s,
        'returncode': returncode,
        'timed_out': timed_out,
        'interrupted': interrupted,
        'checkpoint_restored': checkpoint_restored,
    }
    Path(f'{output_path}.execution.json').write_text(
        json.dumps(record, indent=2, sort_keys=True)
    )


def neb_gpu_accelerated(input_path: str | Path, output_path: str | Path) -> bool:
    """Verify GPU activation in every QE image engine of a NEB calculation.

    ``neb.x`` does not echo the PW GPU banner into its parent output. Each
    image writes the authoritative banner to ``tmp/<prefix>_<image>/PW.out``.
    Requiring all declared image logs prevents one stale or partial GPU marker
    from attesting an otherwise mixed or incomplete calculation.

    Args:
        input_path: NEB input defining the prefix and number of path images.
        output_path: Parent ``neb.x`` output, which may itself attest GPU use.

    Returns:
        Whether the parent output or every declared image log proves that GPU
        acceleration was active.
    """
    parent = Path(output_path)
    if (
        parent.is_file()
        and 'gpu acceleration is active' in parent.read_text(errors='replace').lower()
    ):
        return True
    source = Path(input_path)
    text = source.read_text(errors='replace')
    prefix_match = re.search(r"\bprefix\s*=\s*['\"]([^'\"]+)", text, re.I)
    images_match = re.search(r'\bnum_of_images\s*=\s*(\d+)', text, re.I)
    if not prefix_match or not images_match:
        return False
    prefix = prefix_match.group(1)
    count = int(images_match.group(1))
    logs = [
        source.parent / 'tmp' / f'{prefix}_{index}' / 'PW.out'
        for index in range(1, count + 1)
    ]
    return bool(logs) and all(
        log.is_file()
        and 'gpu acceleration is active' in log.read_text(errors='replace').lower()
        for log in logs
    )


def verify_sssp(elements, directory=SSSP_DIR, manifest=SSSP_MANIFEST) -> dict:
    """Verify required SSSP pseudopotentials against the pinned manifest.

    Args:
        elements: Element symbols whose pseudopotentials must be verified.
        directory: Directory containing the files being verified.
        manifest: Pinned manifest containing expected filenames, checksums, and cutoffs.

    Returns:
        Validated metadata or status; invalid inputs raise an exception.
    """
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
        records[element] = {
            'filename': entry['filename'],
            'md5': digest,
            'cutoff_wfc_Ry': float(entry['cutoff_wfc']),
            'cutoff_rho_Ry': float(entry['cutoff_rho']),
        }
    return {
        'valid': not errors,
        'records': records,
        'errors': errors,
        'ecutwfc_Ry': max((x['cutoff_wfc_Ry'] for x in records.values()), default=None),
        'ecutrho_Ry': max((x['cutoff_rho_Ry'] for x in records.values()), default=None),
    }


def methane_dissociation_images(
    slab: Atoms, active_index: int, n_images: int = 7
) -> list[Atoms]:
    """Build candidate-specific CH4(g)+* -> CH3*+H* NEB endpoints/images.

    Args:
        slab: Slab used by this operation.
        active_index: Active index used by this operation.
        n_images: Number of images to use.

    Returns:
        List of computed or validated records.
    """
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
    initial.set_cell(slab.cell)
    final.set_cell(slab.cell)
    initial.set_pbc(slab.pbc)
    final.set_pbc(slab.pbc)
    images = [initial] + [initial.copy() for _ in range(n_images - 2)] + [final]
    NEB(images, method='improvedtangent').interpolate(method='idpp')
    return images


def neb_checkpoint_images(path: str | Path, templates: list[Atoms]) -> list[Atoms]:
    """Recover geometry only from a QE ``*.path`` restart checkpoint.

    QE stores coordinates in bohr together with energies, gradients, and
    optimizer history. Reusing the restart file also reuses that optimizer
    state. This parser intentionally copies only coordinates onto trusted ASE
    templates, retaining symbols, cells, periodicity, and constraints so a
    newly written ``from_scratch`` input can reset a divergent optimizer.

    Args:
        path: QE NEB ``*.path`` checkpoint.
        templates: Ordered ASE images defining chemical identity and metadata.

    Returns:
        New ASE images carrying checkpoint coordinates in ångströms.

    Raises:
        ValueError: If the checkpoint is incomplete or disagrees with the
            supplied image/atom counts.
    """
    if not templates or any(
        len(image) != len(templates[0])
        or image.get_chemical_symbols() != templates[0].get_chemical_symbols()
        for image in templates
    ):
        raise ValueError('NEB checkpoint templates must share atom ordering')
    lines = Path(path).read_text(errors='replace').splitlines()
    try:
        count_marker = lines.index('NUMBER OF IMAGES')
        image_count = int(lines[count_marker + 1].strip())
        cursor = lines.index('ENERGIES, POSITIONS AND GRADIENTS') + 1
    except (ValueError, IndexError) as exc:
        raise ValueError('invalid QE NEB checkpoint header') from exc
    if image_count != len(templates):
        raise ValueError(
            f'checkpoint has {image_count} images; expected {len(templates)}'
        )

    recovered = []
    atom_count = len(templates[0])
    for expected_index, template in enumerate(templates, 1):
        while cursor < len(lines) and not lines[cursor].lstrip().startswith('Image:'):
            cursor += 1
        if cursor >= len(lines):
            raise ValueError(f'checkpoint is missing image {expected_index}')
        try:
            found_index = int(lines[cursor].split(':', 1)[1])
            cursor += 2  # Skip image marker and energy.
            coordinates = []
            for _ in range(atom_count):
                fields = lines[cursor].split()
                coordinates.append([float(value) for value in fields[:3]])
                cursor += 1
        except (ValueError, IndexError) as exc:
            raise ValueError(
                f'invalid coordinates for checkpoint image {expected_index}'
            ) from exc
        if found_index != expected_index or len(coordinates) != atom_count:
            raise ValueError(f'checkpoint image order mismatch at {expected_index}')
        image = template.copy()
        image.set_positions(np.asarray(coordinates) * Bohr)
        recovered.append(image)
    return recovered


def bound_neb_restart_checkpoint(path: str | Path, final_iteration: int) -> dict:
    """Set an enforceable absolute iteration ceiling in a QE NEB restart.

    On restart, ``neb.x`` takes ``nstep_path`` from the first four lines of
    the ``*.path`` checkpoint and silently overrides the value in the new
    input file.  Updating only the input therefore does not bound GPU use.
    This helper validates the checkpoint header, requires a future ceiling,
    changes only its stored step limit, and records before/after hashes for
    provenance.

    Args:
        path: QE ``*.path`` restart checkpoint to constrain.
        final_iteration: Absolute final path-iteration number permitted.

    Returns:
        Provenance containing current iteration, old/new limits, and hashes.

    Raises:
        ValueError: If the restart header is invalid or the requested ceiling
            does not permit at least one complete future evaluation.
    """
    checkpoint = Path(path)
    original = checkpoint.read_bytes()
    try:
        lines = original.decode(errors='strict').splitlines(keepends=True)
        if lines[0].strip() != 'RESTART INFORMATION':
            raise ValueError
        current_iteration = int(lines[1].strip())
        old_limit = int(lines[2].strip())
        int(lines[3].strip())  # pending-image cursor
    except (IndexError, UnicodeDecodeError, ValueError) as exc:
        raise ValueError('invalid QE NEB restart header') from exc
    limit = int(final_iteration)
    if limit <= current_iteration:
        raise ValueError('NEB restart ceiling must exceed the checkpoint iteration')
    newline = '\r\n' if lines[2].endswith('\r\n') else '\n'
    lines[2] = f'{limit:8d}{newline}'
    updated = ''.join(lines).encode()
    temporary = checkpoint.with_name(f'.{checkpoint.name}.bounded.tmp')
    temporary.write_bytes(updated)
    temporary.replace(checkpoint)
    return {
        'path': str(checkpoint),
        'current_iteration': current_iteration,
        'old_limit': old_limit,
        'new_limit': limit,
        'sha256_before': hashlib.sha256(original).hexdigest(),
        'sha256_after': hashlib.sha256(updated).hexdigest(),
    }


def neb_completed_images(path: str | Path, templates: list[Atoms]) -> list[Atoms]:
    """Recover the last completed NEB geometries from QE's ``*.crd`` file.

    A ``*.path`` restart can contain the optimizer's *next trial* coordinates,
    even when the corresponding force evaluation was interrupted.  QE writes
    ``*.crd`` only after a path iteration has completed, making it the correct
    handoff artifact when changing optimizers or resetting Broyden history.
    Coordinates are in angstrom and are copied onto trusted templates so atom
    identity, cells, periodicity, and constraints remain explicit.

    Args:
        path: QE NEB ``*.crd`` coordinate artifact.
        templates: Ordered ASE images defining chemical identity and metadata.

    Returns:
        New ASE images carrying the last completed path coordinates.

    Raises:
        ValueError: If image count, atom count, symbols, or coordinates are
            incomplete or inconsistent with the trusted templates.
    """
    if not templates or any(
        len(image) != len(templates[0])
        or image.get_chemical_symbols() != templates[0].get_chemical_symbols()
        for image in templates
    ):
        raise ValueError('NEB coordinate templates must share atom ordering')
    lines = Path(path).read_text(errors='replace').splitlines()
    markers = {'FIRST_IMAGE', 'INTERMEDIATE_IMAGE', 'LAST_IMAGE'}
    cursor = 0
    recovered: list[Atoms] = []
    expected_symbols = templates[0].get_chemical_symbols()
    atom_count = len(expected_symbols)
    while cursor < len(lines):
        if lines[cursor].strip() not in markers:
            cursor += 1
            continue
        marker = lines[cursor].strip()
        cursor += 1
        if cursor >= len(lines) or not lines[cursor].strip().upper().startswith(
            'ATOMIC_POSITIONS'
        ):
            raise ValueError(f'missing ATOMIC_POSITIONS after {marker}')
        cursor += 1
        symbols: list[str] = []
        coordinates: list[list[float]] = []
        for _ in range(atom_count):
            if cursor >= len(lines):
                raise ValueError('truncated QE NEB coordinate image')
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) < 4:
                raise ValueError('invalid QE NEB coordinate row')
            try:
                coordinates.append([float(value) for value in fields[1:4]])
            except ValueError as exc:
                raise ValueError('invalid QE NEB coordinate value') from exc
            symbols.append(fields[0])
        if symbols != expected_symbols:
            raise ValueError('QE NEB coordinate atom ordering mismatch')
        if len(recovered) >= len(templates):
            raise ValueError('QE NEB coordinate file has excess images')
        image = templates[len(recovered)].copy()
        image.set_positions(np.asarray(coordinates))
        recovered.append(image)
    if len(recovered) != len(templates):
        raise ValueError(
            f'coordinate file has {len(recovered)} images; '
            f'expected {len(templates)}'
        )
    return recovered


def densify_neb_segment(
    images: list[Atoms], left_index: int, insert_count: int
) -> list[Atoms]:
    """Insert Cartesian interpolation images into one resolved path segment.

    This operation is deliberately local: it preserves every accepted image
    and only adds geometries between a diagnosed pair of neighbors.  It is
    useful when a chemically meaningful coordinate jumps across one NEB
    segment even though the global Cartesian arc length appears uniform.
    The returned geometries are starting points and must be re-evaluated and
    relaxed; interpolation itself is never treated as energy evidence.

    Args:
        images: Ordered, chemically identical NEB images.
        left_index: Zero-based index on the left of the segment to refine.
        insert_count: Number of evenly spaced interior images to insert.

    Returns:
        Copies of the original images with locally interpolated images added.

    Raises:
        ValueError: If the path, segment, or image metadata are inconsistent.
    """
    if len(images) < 2 or insert_count < 1:
        raise ValueError('densification requires a path and inserted images')
    if left_index < 0 or left_index >= len(images) - 1:
        raise ValueError('densification segment index is out of range')
    symbols = images[0].get_chemical_symbols()
    cell = images[0].cell.array
    if any(
        image.get_chemical_symbols() != symbols
        or not np.allclose(image.cell.array, cell)
        for image in images
    ):
        raise ValueError('NEB images must share atom ordering and cell')
    left, right = images[left_index], images[left_index + 1]
    inserted: list[Atoms] = []
    for number in range(1, insert_count + 1):
        fraction = number / (insert_count + 1)
        image = left.copy()
        image.positions = (1.0 - fraction) * left.positions + fraction * right.positions
        inserted.append(image)
    return (
        [image.copy() for image in images[: left_index + 1]]
        + inserted
        + [image.copy() for image in images[left_index + 1 :]]
    )


def redistribute_neb_images_equal_arc(
    images: list[Atoms], image_count: int | None = None
) -> list[Atoms]:
    """Redistribute a periodic path at equal collective Cartesian arc length.

    QE's ``from_scratch`` NEB initialization always redistributes supplied
    images. In QE 7.5, its interpolation cursor advances across at most one
    source segment per target image. A target interval spanning multiple
    short adjacent segments can therefore be extrapolated along the wrong
    segment, creating atom collisions. Supplying an already equi-arc path
    prevents that failure mode while preserving the original polyline.

    Every atom is first unwrapped between source images using its
    minimum-image displacement. Interpolation then follows the full 3N
    Cartesian path and can cross any number of short source segments. These
    geometries remain starting guesses and are never energy evidence.

    Args:
        images: Ordered source images with identical atoms, cell, and PBC.
        image_count: Output image count; defaults to the input count.

    Returns:
        Trusted image copies at equal collective Cartesian arc length.

    Raises:
        ValueError: If metadata is inconsistent, a segment has zero length,
            or fewer than two input/output images are requested.
    """
    if len(images) < 2:
        raise ValueError('equal-arc redistribution requires at least 2 images')
    count = len(images) if image_count is None else int(image_count)
    if count < 2:
        raise ValueError('equal-arc redistribution requires at least 2 outputs')
    symbols = images[0].get_chemical_symbols()
    cell = images[0].cell.array
    pbc = images[0].pbc
    if any(
        image.get_chemical_symbols() != symbols
        or not np.allclose(image.cell.array, cell)
        or not np.array_equal(image.pbc, pbc)
        for image in images
    ):
        raise ValueError('NEB images must share atom ordering, cell, and PBC')

    unwrapped = [np.asarray(images[0].positions, dtype=float).copy()]
    for previous, current in zip(images, images[1:]):
        displacement = current.positions - previous.positions
        minimum_displacement, _ = find_mic(displacement, cell, pbc=pbc)
        unwrapped.append(unwrapped[-1] + minimum_displacement)
    segment_lengths = np.asarray(
        [np.linalg.norm(right - left) for left, right in zip(unwrapped, unwrapped[1:])]
    )
    if not np.all(np.isfinite(segment_lengths)) or np.any(segment_lengths <= 1.0e-12):
        raise ValueError('NEB source path contains a zero-length segment')
    cumulative = np.concatenate(([0.0], np.cumsum(segment_lengths)))
    targets = np.linspace(0.0, cumulative[-1], count)

    redistributed: list[Atoms] = []
    for output_index, target in enumerate(targets):
        if output_index == 0:
            positions = unwrapped[0]
            template = images[0]
        elif output_index == count - 1:
            positions = unwrapped[-1]
            template = images[-1]
        else:
            # Unlike QE 7.5's single increment, searchsorted can cross any
            # number of short segments before locating the target.
            segment = int(np.searchsorted(cumulative, target, side='right') - 1)
            segment = min(segment, len(segment_lengths) - 1)
            fraction = (target - cumulative[segment]) / segment_lengths[segment]
            positions = unwrapped[segment] + fraction * (
                unwrapped[segment + 1] - unwrapped[segment]
            )
            template = images[min(output_index, len(images) - 1)]
        image = template.copy()
        image.set_positions(positions)
        redistributed.append(image)
    return redistributed


def write_qe_neb_input(
    images: list[Atoms],
    path: str,
    prefix: str,
    path_config: NEBPathConfig | None = None,
) -> dict:
    """Write a climbing-image neb.x input using one verified SSSP family.

    Args:
        images: Ordered values supplying images.
        path: Filesystem path to the input or output artifact.
        prefix: Prefix used by this operation.
        path_config: Validated numerical and physical controls for the path;
            defaults to strict climbing-image settings when omitted.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    config = path_config or NEBPathConfig()
    config.validate()
    if not images or any(len(x) != len(images[0]) for x in images):
        raise ValueError('NEB images must have identical atom ordering')
    engine_images = (
        redistribute_neb_images_equal_arc(images)
        if config.restart_mode == 'from_scratch' and config.precondition_equal_arc
        else [image.copy() for image in images]
    )
    elements = sorted(set(engine_images[0].get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    records = verified['records']
    species = '\n'.join(f" {e} 1.0 {records[e]['filename']}" for e in elements)
    cell = '\n'.join(
        ' '.join(f'{v:.12f}' for v in row) for row in engine_images[0].cell.array
    )
    positions = []
    for number, image in enumerate(engine_images):
        tag = (
            'FIRST_IMAGE'
            if number == 0
            else (
                'LAST_IMAGE'
                if number == len(engine_images) - 1
                else 'INTERMEDIATE_IMAGE'
            )
        )
        positions.append(
            f'{tag}\nATOMIC_POSITIONS angstrom\n'
            + _qe_positions(image, freeze_all=config.freeze_all_atoms)
        )
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    magnetization = '\n'.join(
        f" starting_magnetization({i})={0.5 if element in magnetic else 0.05}"
        for i, element in enumerate(elements, 1)
    )
    spin_constraint = (
        ''
        if config.total_magnetization is None
        else f'\n tot_magnetization={config.total_magnetization}'
    )
    text = f"""BEGIN
BEGIN_PATH_INPUT
&PATH
 restart_mode='{config.restart_mode}',
 num_of_images={len(engine_images)}, opt_scheme='{config.optimizer}',
 CI_scheme='{config.climbing_scheme}', nstep_path={config.nstep_path},
 path_thr={config.path_thr_eV_A}, ds={config.step_size},
 minimum_image={'.true.' if config.minimum_image else '.false.'},
 k_min={config.spring_min}, k_max={config.spring_max},
/
END_PATH_INPUT
BEGIN_ENGINE_INPUT
&CONTROL
 calculation='scf', prefix='{prefix}', pseudo_dir='{SSSP_DIR}', outdir='./tmp'
/
&SYSTEM
 ibrav=0, nat={len(engine_images[0])}, ntyp={len(elements)},
 ecutwfc={verified['ecutwfc_Ry']}, ecutrho={verified['ecutrho_Ry']},
 occupations='smearing', smearing='mv', degauss=0.02, nspin=2
{magnetization}{spin_constraint}
/
&ELECTRONS
 conv_thr={config.scf_conv_thr_Ry:.1e}, electron_maxstep={config.electron_maxstep},
 mixing_mode='local-TF', mixing_beta={config.mixing_beta},
 mixing_ndim={config.mixing_ndim}, startingpot='{config.starting_potential}',
 startingwfc='{config.starting_wavefunctions}',
 diagonalization='{config.diagonalization}'
/
ATOMIC_SPECIES
{species}
CELL_PARAMETERS angstrom
{cell}
K_POINTS automatic
 {config.kpoints[0]} {config.kpoints[1]} {config.kpoints[2]} 0 0 0
BEGIN_POSITIONS
{chr(10).join(positions)}
END_POSITIONS
END_ENGINE_INPUT
END
"""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return {
        'path': str(target),
        'n_images': len(engine_images),
        'sssp': verified,
        'equal_arc_preconditioned': (
            config.restart_mode == 'from_scratch' and config.precondition_equal_arc
        ),
    }


def prepare_neb_image_scfs(
    images: list[Atoms],
    directory: str | Path,
    prefix: str,
    path_config: NEBPathConfig | None = None,
) -> dict:
    """Write independent SCF inputs that seed a difficult NEB calculation.

    A spin-polarized metallic path can fail before its first ionic iteration
    when every image initializes from an atomic charge superposition.  These
    jobs converge the *same geometries and electronic controls* independently.
    Their ``<prefix>.save`` directories can then be consumed by an NEB input
    configured with ``startingpot='file'``.  No energy from an unconverged seed
    is accepted as path evidence.

    Args:
        images: Ordered NEB geometries with identical atom ordering.
        directory: Directory receiving ``image_XX.in`` inputs and ``tmp``.
        prefix: Prefix used by the subsequent NEB calculation.
        path_config: Electronic controls shared with that NEB tier.

    Returns:
        Manifest describing every generated image input and output directory.
    """
    config = path_config or NEBPathConfig.coarse_preconditioner()
    config.validate()
    if not images or any(
        len(image) != len(images[0])
        or image.get_chemical_symbols() != images[0].get_chemical_symbols()
        for image in images
    ):
        raise ValueError('NEB images must have identical atom ordering')
    elements = sorted(set(images[0].get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    records = verified['records']
    species = '\n'.join(
        f" {element} 1.0 {records[element]['filename']}" for element in elements
    )
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    magnetization = '\n'.join(
        f" starting_magnetization({index})=" f"{0.5 if element in magnetic else 0.05}"
        for index, element in enumerate(elements, 1)
    )
    spin_constraint = (
        ''
        if config.total_magnetization is None
        else f'\n tot_magnetization={config.total_magnetization}'
    )
    cell = '\n'.join(
        ' '.join(f'{value:.12f}' for value in row) for row in images[0].cell.array
    )
    root = Path(directory)
    root.mkdir(parents=True, exist_ok=True)
    (root / 'tmp').mkdir(exist_ok=True)
    jobs = []
    for index, image in enumerate(images, 1):
        image_outdir = f'./tmp/{prefix}_{index}'
        text = f"""&CONTROL
 calculation='scf', restart_mode='from_scratch', prefix='{prefix}',
 pseudo_dir='{SSSP_DIR}', outdir='{image_outdir}', tprnfor=.true.
/
&SYSTEM
 ibrav=0, nat={len(image)}, ntyp={len(elements)},
 ecutwfc={verified['ecutwfc_Ry']}, ecutrho={verified['ecutrho_Ry']},
 occupations='smearing', smearing='mv', degauss=0.02, nspin=2
{magnetization}{spin_constraint}
/
&ELECTRONS
 conv_thr={config.scf_conv_thr_Ry:.1e},
 electron_maxstep={config.electron_maxstep},
 mixing_mode='local-TF', mixing_beta={config.mixing_beta},
 mixing_ndim={config.mixing_ndim},
 startingpot='{config.starting_potential}',
 startingwfc='{config.starting_wavefunctions}',
 diagonalization='{config.diagonalization}'
/
ATOMIC_SPECIES
{species}
CELL_PARAMETERS angstrom
{cell}
K_POINTS automatic
 {config.kpoints[0]} {config.kpoints[1]} {config.kpoints[2]} 0 0 0
ATOMIC_POSITIONS angstrom
{_qe_positions(image)}
"""
        input_path = root / f'image_{index:02d}.in'
        input_path.write_text(text)
        jobs.append(
            {
                'image': index,
                'input': str(input_path),
                'output': str(root / f'image_{index:02d}.out'),
                'outdir': str(root / 'tmp' / f'{prefix}_{index}'),
            }
        )
    manifest = {
        'prefix': prefix,
        'image_count': len(images),
        'configuration': asdict(config),
        'jobs': jobs,
    }
    (root / 'image_scf_manifest.json').write_text(
        json.dumps(manifest, indent=2, sort_keys=True)
    )
    return manifest


def run_neb(
    input_path: str,
    output_path: str,
    timeout_s: int = 86400,
    execution: QEExecutionConfig | None = None,
    restart_checkpoint: str | Path | None = None,
    final_iteration: int | None = None,
) -> SolverRunResult:
    """Execute a Quantum ESPRESSO NEB calculation and record its provenance.

    Args:
        input_path: Quantum ESPRESSO input file.
        output_path: Destination for solver output and execution provenance.
        timeout_s: Maximum solver wall time in seconds.
        execution: Injected executable-discovery and process-execution services.
        restart_checkpoint: Optional QE ``*.path`` checkpoint whose stored
            iteration limit must be bounded before execution.
        final_iteration: Absolute final iteration allowed for a restart. It is
            required with ``restart_checkpoint`` and must match the input.

    Returns:
        A dictionary containing run neb outputs, status, and supporting metadata.
    """
    if (restart_checkpoint is None) != (final_iteration is None):
        raise ValueError(
            'restart_checkpoint and final_iteration must be supplied together'
        )
    restart_bound = None
    checkpoint_path = (
        Path(restart_checkpoint) if restart_checkpoint is not None else None
    )
    checkpoint_before = (
        checkpoint_path.read_bytes() if checkpoint_path is not None else None
    )
    if restart_checkpoint is not None:
        input_text = Path(input_path).read_text(errors='strict')
        if not re.search(
            r"restart_mode\s*=\s*['\"]restart['\"]", input_text, flags=re.IGNORECASE
        ):
            raise ValueError('bounded NEB execution requires restart_mode=restart')
        match = re.search(r'nstep_path\s*=\s*(\d+)', input_text, flags=re.IGNORECASE)
        if match is None or int(match.group(1)) != int(final_iteration):
            raise ValueError('NEB input nstep_path must match the restart ceiling')
        restart_bound = bound_neb_restart_checkpoint(
            checkpoint_path, int(final_iteration)
        )
        Path(f'{output_path}.restart_bound.json').write_text(
            json.dumps(restart_bound, indent=2, sort_keys=True)
        )

    neb = resolve_qe_executable('neb.x')
    workdir = Path(input_path).resolve().parent
    (workdir / 'tmp').mkdir(exist_ok=True)
    config = execution or QEExecutionConfig()
    command = build_qe_command(neb, input_path, config, neb=True)
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = str(config.omp_threads)
    started = time.monotonic()
    timed_out, returncode = False, -1
    checkpoint_restored = False
    # A sidecar from an earlier attempt must never be mistaken for provenance
    # of a currently running or operator-interrupted calculation.
    Path(f'{output_path}.execution.json').unlink(missing_ok=True)
    try:
        with open(output_path, 'w') as sink:
            proc = subprocess.run(
                command,
                cwd=str(workdir),
                env=environment,
                stdout=sink,
                stderr=subprocess.STDOUT,
                timeout=timeout_s,
            )
            returncode = proc.returncode
    except subprocess.TimeoutExpired:
        timed_out = True
        if checkpoint_path is not None and checkpoint_before is not None:
            temporary = checkpoint_path.with_name(
                f'.{checkpoint_path.name}.pre_timeout.tmp'
            )
            temporary.write_bytes(checkpoint_before)
            temporary.replace(checkpoint_path)
            checkpoint_restored = True
    except BaseException:
        # An operator interrupt can terminate neb.x while it is rewriting the
        # formatted restart file. Preserve the last accepted path just as for
        # a wall-clock timeout, then propagate the interrupt to the caller.
        if checkpoint_path is not None and checkpoint_before is not None:
            temporary = checkpoint_path.with_name(
                f'.{checkpoint_path.name}.pre_interrupt.tmp'
            )
            temporary.write_bytes(checkpoint_before)
            temporary.replace(checkpoint_path)
            checkpoint_restored = True
        elapsed = time.monotonic() - started
        _execution_record(
            input_path,
            output_path,
            command,
            config,
            elapsed,
            -1,
            False,
            interrupted=True,
            checkpoint_restored=checkpoint_restored,
        )
        raise
    elapsed = time.monotonic() - started
    _execution_record(
        input_path,
        output_path,
        command,
        config,
        elapsed,
        returncode,
        timed_out,
        checkpoint_restored=checkpoint_restored,
    )
    text = Path(output_path).read_text(errors='replace')
    lower = text.lower()
    gpu_active = neb_gpu_accelerated(input_path, output_path)
    converged = (
        returncode == 0
        and 'job done' in lower
        and 'neb: convergence achieved' in lower
        and 'error in routine' not in lower
        and gpu_active
    )
    result: SolverRunResult = {
        'converged': converged,
        'returncode': returncode,
        'output': str(output_path),
        'candidate_specific': True,
        'execution': asdict(config),
        'timed_out': timed_out,
        'gpu_accelerated': gpu_active,
    }
    if restart_bound is not None:
        result['restart_bound'] = restart_bound
    if checkpoint_restored:
        result['checkpoint_restored'] = True
    return result


def write_qe_relax_input(
    atoms: Atoms, path: str, prefix: str, kpoints=(2, 2, 1)
) -> dict:
    """Write a spin-polarized, force-converged endpoint relaxation.

    Args:
        atoms: Atomic structure consumed by the calculator.
        path: Filesystem path to the input or output artifact.
        prefix: Prefix used by this operation.
        kpoints: Kpoints used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    elements = sorted(set(atoms.get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    species = '\n'.join(
        f"{e} 1.0 {verified['records'][e]['filename']}" for e in elements
    )
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    mags = '\n'.join(
        f" starting_magnetization({i})={0.5 if e in magnetic else 0.05}"
        for i, e in enumerate(elements, 1)
    )
    cell = '\n'.join(' '.join(f'{v:.12f}' for v in row) for row in atoms.cell.array)
    positions = _qe_positions(atoms)
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
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return {'path': str(target), 'sssp': verified}


def write_qe_force_input(
    atoms: Atoms, path: str, prefix: str, kpoints=(2, 2, 1)
) -> dict:
    """Write a fixed-geometry SCF input that prints Cartesian atomic forces.

    Args:
        atoms: Atomic structure consumed by the calculator.
        path: Filesystem path to the input or output artifact.
        prefix: Prefix used by this operation.
        kpoints: Kpoints used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    elements = sorted(set(atoms.get_chemical_symbols()))
    verified = verify_sssp(elements)
    if not verified['valid']:
        raise RuntimeError(f'SSSP verification failed: {verified["errors"]}')
    species = '\n'.join(
        f"{element} 1.0 {verified['records'][element]['filename']}"
        for element in elements
    )
    magnetic = {'Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Gd', 'Ce', 'Eu'}
    mags = '\n'.join(
        f" starting_magnetization({i})={0.5 if element in magnetic else 0.05}"
        for i, element in enumerate(elements, 1)
    )
    cell = '\n'.join(
        ' '.join(f'{value:.12f}' for value in row) for row in atoms.cell.array
    )
    positions = _qe_positions(atoms, include_constraints=False)
    text = f"""&CONTROL
 calculation='scf', prefix='{prefix}', pseudo_dir='{SSSP_DIR}', outdir='./tmp',
 tprnfor=.true.
/
&SYSTEM
 ibrav=0, nat={len(atoms)}, ntyp={len(elements)}, ecutwfc={verified['ecutwfc_Ry']},
 ecutrho={verified['ecutrho_Ry']}, occupations='smearing', smearing='mv', degauss=0.02, nspin=2
{mags}
/
&ELECTRONS
 conv_thr=1.0d-8, mixing_beta=0.3
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
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(text)
    return {'path': str(target), 'sssp': verified}


def run_pw(
    input_path: str,
    output_path: str,
    timeout_s: int = 86400,
    execution: QEExecutionConfig | None = None,
) -> SolverRunResult:
    """Execute a Quantum ESPRESSO pw.x calculation and record its provenance.

    Args:
        input_path: Quantum ESPRESSO input file.
        output_path: Destination for solver output and execution provenance.
        timeout_s: Maximum solver wall time in seconds.
        execution: Injected executable-discovery and process-execution services.

    Returns:
        A dictionary containing run pw outputs, status, and supporting metadata.
    """
    workdir = Path(input_path).resolve().parent
    (workdir / 'tmp').mkdir(exist_ok=True)
    pw = resolve_qe_executable('pw.x')
    config = execution or QEExecutionConfig.production_default()
    command = build_qe_command(pw, input_path, config)
    environment = os.environ.copy()
    environment['OMP_NUM_THREADS'] = str(config.omp_threads)
    started = time.monotonic()
    timed_out, returncode = False, -1
    with open(output_path, 'w') as sink:
        try:
            proc = subprocess.run(
                command,
                stdout=sink,
                stderr=subprocess.STDOUT,
                cwd=str(workdir),
                env=environment,
                timeout=timeout_s,
            )
            returncode = proc.returncode
        except subprocess.TimeoutExpired:
            timed_out = True
    elapsed = time.monotonic() - started
    _execution_record(
        input_path, output_path, command, config, elapsed, returncode, timed_out
    )
    text = Path(output_path).read_text(errors='replace')
    lower = text.lower()
    gpu_active = 'gpu acceleration is active' in lower
    # ``prefix.EXIT`` requests a clean, restartable QE shutdown. Such a run
    # returns zero and prints ``JOB DONE`` but is not an SCF solution. Require
    # QE's explicit convergence declaration so a graceful checkpoint can never
    # be promoted to scientific evidence.
    converged = (
        returncode == 0
        and 'job done' in lower
        and 'convergence has been achieved' in lower
        and 'convergence not achieved' not in lower
        and gpu_active
    )
    return {
        'converged': converged,
        'returncode': returncode,
        'output': output_path,
        'execution': asdict(config),
        'timed_out': timed_out,
        'gpu_accelerated': gpu_active,
    }


def relaxed_structure(output_path: str, input_path: str | Path | None = None) -> Atoms:
    """Read the last geometry only from a cleanly completed QE relaxation.

    Args:
        output_path: Filesystem destination for generated output.
        input_path: Optional QE input whose constraints are restored on the
            relaxed geometry. QE output does not reliably preserve ``if_pos``
            flags when ASE reads the final structure.

    Returns:
        Computed `Atoms` result.
    """
    text = Path(output_path).read_text(errors='replace')
    if 'JOB DONE' not in text or 'convergence NOT achieved' in text:
        raise RuntimeError('endpoint relaxation is not converged')
    relaxed = ase_read(output_path, format='espresso-out', index=-1)
    if input_path is not None:
        source = ase_read(input_path, format='espresso-in')
        if source.get_chemical_symbols() != relaxed.get_chemical_symbols():
            raise RuntimeError(
                'QE relaxation input/output have different atom ordering'
            )
        relaxed.set_constraint(source.constraints)
    return relaxed


def parse_atomic_forces(
    output_path: str, expected_atoms: int | None = None
) -> np.ndarray:
    """Read the final QE force block and return forces in eV/angstrom.

    Args:
        output_path: Filesystem destination for generated output.
        expected_atoms: Expected atoms used by this operation.

    Returns:
        Computed `np.ndarray` result.
    """
    text = Path(output_path).read_text(errors='replace')
    blocks = re.findall(
        r'Forces acting on atoms[^\n]*\n(.*?)(?=\n\s*Total force|\n\s*!|\Z)',
        text,
        flags=re.IGNORECASE | re.DOTALL,
    )
    if not blocks:
        raise RuntimeError('QE output contains no atomic-force block')
    rows = re.findall(
        r'atom\s+\d+\s+type\s+\d+\s+force\s*=\s*'
        r'([-+0-9.Ee]+)\s+([-+0-9.Ee]+)\s+([-+0-9.Ee]+)',
        blocks[-1],
        flags=re.IGNORECASE,
    )
    if expected_atoms is not None and len(rows) != expected_atoms:
        raise RuntimeError(f'expected {expected_atoms} force rows, found {len(rows)}')
    if not rows:
        raise RuntimeError('QE force block contains no parseable atoms')
    # Quantum ESPRESSO reports forces in Ry/Bohr.
    ry_bohr_to_ev_ang = 13.605693122994 / 0.529177210903
    return np.asarray(rows, dtype=float) * ry_bohr_to_ev_ang


def parse_neb_result(output_path: str) -> NEBResult:
    """Parse a completed NEB output into convergence and barrier evidence.

    Args:
        output_path: Destination for solver output and execution provenance.

    Returns:
        Convergence state, image energies, and activation barrier parsed from output.
    """
    text = Path(output_path).read_text(errors='replace')
    energies = [
        float(x)
        for x in re.findall(r'activation energy \(->\)\s*=\s*([-+0-9.Ee]+)', text)
    ]
    reverse = [
        float(x)
        for x in re.findall(r'activation energy \(<-\)\s*=\s*([-+0-9.Ee]+)', text)
    ]
    errors = [float(x) for x in re.findall(r'path length\s*=\s*([-+0-9.Ee]+)', text)]
    converged = 'JOB DONE' in text and 'neb: convergence achieved' in text.lower()
    return {
        'converged': converged,
        'forward_barrier_eV': energies[-1] if energies else None,
        'reverse_barrier_eV': reverse[-1] if reverse else None,
        'path_metric': errors[-1] if errors else None,
        'candidate_specific': True,
    }


def partial_hessian(
    forces_plus: np.ndarray,
    forces_minus: np.ndarray,
    displacement_A: float,
    masses_amu: np.ndarray,
    reaction_direction: np.ndarray | None = None,
    min_reaction_mode_overlap: float = 0.5,
) -> TransitionStateFrequencyResult:
    """Construct a mass-weighted partial Hessian from central force differences.

    Args:
        forces_plus: Forces plus used by this operation.
        forces_minus: Forces minus used by this operation.
        displacement_A: Displacement in ångströms.
        masses_amu: Masses amu used by this operation.
        reaction_direction: Optional Cartesian NEB tangent for the active atoms.
        min_reaction_mode_overlap: Minimum absolute overlap between the single
            imaginary normal mode and the mass-weighted reaction direction.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    plus = np.asarray(forces_plus, float)
    minus = np.asarray(forces_minus, float)
    if plus.shape != minus.shape or plus.ndim != 3 or displacement_A <= 0:
        raise ValueError('forces must be (3N, N, 3) central-difference arrays')
    n_atoms = plus.shape[1]
    if plus.shape[0] != 3 * n_atoms or len(masses_amu) != n_atoms:
        raise ValueError('partial Hessian dimensions are inconsistent')
    if not 0.0 <= min_reaction_mode_overlap <= 1.0:
        raise ValueError('minimum reaction-mode overlap must be in [0, 1]')
    hessian = -(plus - minus).reshape(3 * n_atoms, 3 * n_atoms).T / (2 * displacement_A)
    hessian = 0.5 * (hessian + hessian.T)
    weights = np.repeat(np.sqrt(np.asarray(masses_amu, float)), 3)
    eigvals, eigvecs = np.linalg.eigh(hessian / np.outer(weights, weights))
    # 1 eV/A^2/amu -> (521.47083 cm^-1)^2
    frequencies = np.sign(eigvals) * np.sqrt(np.abs(eigvals)) * 521.47083
    imaginary_indices = np.flatnonzero(frequencies < -50.0)
    reaction_overlap = None
    reaction_mode_valid = None
    if reaction_direction is not None:
        direction = np.asarray(reaction_direction, float)
        if direction.shape == (n_atoms, 3):
            direction = direction.reshape(-1)
        if direction.shape != (3 * n_atoms,) or not np.all(np.isfinite(direction)):
            raise ValueError('reaction direction must contain 3N finite values')
        mass_weighted = direction * weights
        norm = float(np.linalg.norm(mass_weighted))
        if norm <= np.finfo(float).eps:
            raise ValueError('reaction direction must be nonzero')
        if len(imaginary_indices) == 1:
            reaction_overlap = float(
                abs(np.dot(eigvecs[:, int(imaginary_indices[0])], mass_weighted / norm))
            )
            reaction_mode_valid = reaction_overlap >= min_reaction_mode_overlap
        else:
            reaction_mode_valid = False
    frequency_valid = len(imaginary_indices) == 1
    if reaction_mode_valid is not None:
        frequency_valid = frequency_valid and reaction_mode_valid
    return {
        'frequencies_cm1': frequencies.tolist(),
        'imaginary_count': int(len(imaginary_indices)),
        'valid_transition_state': frequency_valid,
        'mode_vectors': eigvecs.tolist(),
        'reaction_mode_overlap': reaction_overlap,
        'reaction_mode_valid': reaction_mode_valid,
    }
