"""Quantum ESPRESSO NEB path configuration presets, checkpoint recovery, and interpolation."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.build import molecule
from ase.geometry import find_mic
from ase.mep import NEB
from ase.units import Bohr


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
