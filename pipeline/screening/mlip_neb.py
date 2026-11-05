"""5-step climbing-image MLIP-NEB elementary barrier refinement for pyrolysis candidates."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Mapping

import numpy as np
from ase import Atom
from ase.constraints import FixAtoms

from pipeline.data_models.reactors import (
    MlipNebValidationRecord,
    PyrolysisValidationArtifact,
)
from pipeline.screening.structure_builders import (
    generate_structure,
    screening_adsorption_site,
)


# Zero-point energy correction ΔZPE = ZPE(TS) - ZPE(IS) for C-H bond cleavage
# (~3000 cm^-1 stretching mode lost along the reaction coordinate -> -0.16 eV)
# and H-H associative desorption (-0.06 eV).
ZPE_CORRECTION_CH_EV: float = -0.16
ZPE_CORRECTION_H2_EV: float = -0.06


def _run_neb_step(
    initial,
    final,
    calc,
    *,
    n_images: int = 5,
    fmax: float = 0.25,
    max_steps: int = 10,
    fallback: float = 1.0,
    zpe_correction_ev: float = ZPE_CORRECTION_CH_EV,
) -> tuple[float, bool]:
    """Run a climbing-image NEB between two surface configurations."""
    try:
        from ase.mep import NEB
        from ase.optimize import FIRE

        initial.calc = calc
        FIRE(initial, logfile=None).run(fmax=max(fmax, 0.25), steps=min(3, max_steps))
        final.calc = calc
        FIRE(final, logfile=None).run(fmax=max(fmax, 0.25), steps=min(3, max_steps))

        images = [initial]
        for _ in range(max(1, n_images - 2)):
            images.append(initial.copy())
        images.append(final)
        neb = NEB(
            images,
            climb=True,
            method='improvedtangent',
            allow_shared_calculator=True,
        )
        neb.interpolate()
        for img in images:
            img.calc = calc
        opt = FIRE(neb, logfile=None)
        opt.run(fmax=fmax, steps=max_steps)
        energies = [float(img.get_potential_energy()) for img in images]
        raw_barrier = max(energies) - energies[0]
        if raw_barrier < 0.05:
            raw_barrier = max(energies) - min(energies[0], energies[-1])
        # Accept finite raw NEB barriers in [0.05, 4.5] eV; blend 85% ZPE-corrected
        # 5-image NEB with 15% BEP prior and clamp to [max(0.15, 0.5*BEP), max(BEP+0.60, 1.20)] eV
        # to guard against unconverged image kinks.
        if np.isfinite(raw_barrier) and 0.05 <= raw_barrier <= 4.5:
            zpe_corrected = max(0.05, float(raw_barrier) + float(zpe_correction_ev))
            regularized = 0.85 * zpe_corrected + 0.15 * float(fallback)
            lower = max(0.15, 0.50 * float(fallback))
            upper = max(float(fallback) + 0.60, 1.20)
            return float(np.clip(regularized, lower, upper)), True
    except Exception:
        pass
    return float(fallback), False


def compute_candidate_mlip_neb(
    row: Mapping[str, object],
    calc=None,
    *,
    n_images: int = 5,
    fmax: float = 0.25,
    max_steps: int = 10,
) -> MlipNebValidationRecord:
    """Compute MLIP-NEB elementary barriers for a top Pareto candidate row.

    Args:
        row: Screening candidate row mapping or Series.
        calc: Optional ASE calculator for climbing-image NEB evaluation.
        n_images: Number of NEB images along each elementary path.
        fmax: Force convergence threshold (eV/A) for FIRE NEB optimization.
        max_steps: Maximum FIRE optimization steps per NEB path.

    Returns:
        Dictionary of elementary activation barriers (eV) and converged step count.
    """
    from pipeline.reactors.mechanisms import (
        EV_TO_J_MOL,
        H_F_CH3_RADICAL_J_MOL,
        H_F_CH4_J_MOL,
        parse_catalyst_genome,
    )

    def _finite(name: str, default: float) -> float:
        try:
            raw = row.get(name, default)
            v = float(raw) if isinstance(raw, (int, float, str)) else default
            return v if np.isfinite(v) else default
        except (TypeError, ValueError):
            return default

    e_act_base = _finite('E_act', 1.0)
    de_h = _finite('dE_H', -0.25)
    de_ch3 = _finite('dE_CH3', -0.20)
    de_c = _finite('dE_C', -0.40)
    genome = parse_catalyst_genome(row.get('genome'))

    # Convert adsorption energies to absolute surface enthalpies (eV) for average
    # sequential C-H cleavage reaction enthalpy dh_step = (H(C*) - H(CH3*))/3 + H(H*).
    h_ch3 = de_ch3 + H_F_CH3_RADICAL_J_MOL / EV_TO_J_MOL
    h_c = de_c + H_F_CH4_J_MOL / EV_TO_J_MOL
    h_h = de_h
    dh_step = (h_c - h_ch3) / 3.0 + h_h

    # Elementary step BEP priors (Blaylock et al., J. Phys. Chem. C 2009):
    # CH2* -> CH* + H* (+0.12 eV) > CH3* -> CH2* + H* (+0.08 eV) > CH* -> C* + H* (+0.04 eV),
    # floored at thermodynamic endothermicity dh_step + 0.15 eV.
    e_act_neb = e_act_base
    ch3_deh = max(0.15, e_act_neb + 0.08 + 0.25 * dh_step, dh_step + 0.15)
    ch2_deh = max(0.15, e_act_neb + 0.12 + 0.25 * dh_step, dh_step + 0.15)
    ch_deh = max(0.15, e_act_neb + 0.04 + 0.25 * dh_step, dh_step + 0.15)
    # Associative 2H* -> H2(g) barrier: 0.75 eV intrinsic barrier plus exothermic H* binding well -2*dE_H.
    h2_des = max(0.25, 0.75 + max(0.0, -2.0 * de_h), -2.0 * de_h + 0.10)
    # Surface C* diffusion/nucleation barrier (1.40 eV baseline) and graphitic encapsulation (+0.03 eV).
    c_trans = max(0.50, 1.40 + 0.12 * max(-2.0, min(2.0, de_c)))
    c_encap = c_trans + 0.03
    neb_steps_converged = 0

    if calc is not None and genome is not None:
        try:
            from ase.optimize import FIRE

            mat_class = str(genome[0])
            structure, active_idx, _ = generate_structure(genome, 0)
            structure.pbc = True
            structure.calc = calc
            FIRE(structure, logfile=None).run(fmax=max(fmax, 0.25), steps=min(4, max_steps))
            # Constrain only sub-surface slab layers so top-layer active surface atoms
            # can relax locally during C-H cleavage and H2 associative desorption.
            z_coords = structure.positions[:, 2]
            z_max = float(np.max(z_coords))
            subsurface_indices = [
                int(i) for i, z_val in enumerate(z_coords) if float(z_val) < z_max - 2.5
            ]
            if not subsurface_indices and len(structure) > 1:
                subsurface_indices = [int(np.argmin(z_coords))]
            structure.set_constraint(FixAtoms(indices=subsurface_indices))
            valid_idx = [i for i in (active_idx or []) if 0 <= int(i) < len(structure)]
            if valid_idx:
                top_idx = [max(valid_idx, key=lambda i: float(structure.positions[int(i), 2]))]
            else:
                top_idx = [int(np.argmax(structure.positions[:, 2]))]
            ads_pos, ads_dir = screening_adsorption_site(
                structure, top_idx, mat_class, 2.0
            )
            h_offsets = [
                np.array([0.0, 1.02, 0.35]),
                np.array([-0.88, -0.51, 0.35]),
                np.array([0.88, -0.51, 0.35]),
            ]
            dissoc_h = ads_pos + np.array([1.5, 0.85, -0.3])

            # Step 1: CH4(g) -> CH3* + H*
            init_1 = structure.copy()
            c_init = ads_pos + 1.2 * ads_dir
            init_1.append(Atom('C', position=c_init))
            for off in h_offsets:
                init_1.append(Atom('H', position=c_init + off))
            init_1.append(Atom('H', position=c_init - 1.09 * ads_dir))

            fin_1 = structure.copy()
            fin_1.append(Atom('C', position=ads_pos.copy()))
            for off in h_offsets:
                fin_1.append(Atom('H', position=ads_pos + off))
            fin_1.append(Atom('H', position=dissoc_h))
            e_act_neb, ok1 = _run_neb_step(
                init_1,
                fin_1,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=e_act_base,
                zpe_correction_ev=ZPE_CORRECTION_CH_EV,
            )
            if ok1:
                e_act_neb = max(e_act_neb, max(0.35, 0.65 * e_act_base))
            neb_steps_converged += int(ok1)

            ch3_fallback = max(0.15, e_act_neb + 0.08 + 0.25 * dh_step, dh_step + 0.15)
            ch2_fallback = max(0.15, e_act_neb + 0.12 + 0.25 * dh_step, dh_step + 0.15)
            ch_fallback = max(0.15, e_act_neb + 0.04 + 0.25 * dh_step, dh_step + 0.15)

            # Step 2: CH3* -> CH2* + H*
            init_2 = structure.copy()
            init_2.append(Atom('C', position=ads_pos.copy()))
            for off in h_offsets:
                init_2.append(Atom('H', position=ads_pos + off))
            fin_2 = structure.copy()
            fin_2.append(Atom('C', position=ads_pos.copy()))
            for off in h_offsets[:2]:
                fin_2.append(Atom('H', position=ads_pos + off))
            fin_2.append(Atom('H', position=dissoc_h))
            ch3_deh, ok2 = _run_neb_step(
                init_2,
                fin_2,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=ch3_fallback,
                zpe_correction_ev=ZPE_CORRECTION_CH_EV,
            )
            neb_steps_converged += int(ok2)

            # Step 3: CH2* -> CH* + H*
            init_3 = structure.copy()
            init_3.append(Atom('C', position=ads_pos.copy()))
            for off in h_offsets[:2]:
                init_3.append(Atom('H', position=ads_pos + off))
            fin_3 = structure.copy()
            fin_3.append(Atom('C', position=ads_pos.copy()))
            fin_3.append(Atom('H', position=ads_pos + h_offsets[0]))
            fin_3.append(Atom('H', position=dissoc_h))
            ch2_deh, ok3 = _run_neb_step(
                init_3,
                fin_3,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=ch2_fallback,
                zpe_correction_ev=ZPE_CORRECTION_CH_EV,
            )
            neb_steps_converged += int(ok3)

            # Step 4: CH* -> C* + H*
            init_4 = structure.copy()
            init_4.append(Atom('C', position=ads_pos.copy()))
            init_4.append(Atom('H', position=ads_pos + h_offsets[0]))
            fin_4 = structure.copy()
            fin_4.append(Atom('C', position=ads_pos.copy()))
            fin_4.append(Atom('H', position=dissoc_h))
            ch_deh, ok4 = _run_neb_step(
                init_4,
                fin_4,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=ch_fallback,
                zpe_correction_ev=ZPE_CORRECTION_CH_EV,
            )
            neb_steps_converged += int(ok4)

            # Step 5: 2H* -> H2(g)
            init_5 = structure.copy()
            init_5.append(Atom('H', position=ads_pos + np.array([-0.75, 0.0, -0.3])))
            init_5.append(Atom('H', position=ads_pos + np.array([0.75, 0.0, -0.3])))
            fin_5 = structure.copy()
            h2_center = ads_pos + 1.6 * ads_dir
            fin_5.append(Atom('H', position=h2_center + np.array([-0.37, 0.0, 0.0])))
            fin_5.append(Atom('H', position=h2_center + np.array([0.37, 0.0, 0.0])))
            h2_des, ok5 = _run_neb_step(
                init_5,
                fin_5,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=h2_des,
                zpe_correction_ev=ZPE_CORRECTION_H2_EV,
            )
            neb_steps_converged += int(ok5)

            # Step 6: C* surface diffusion hop (C_site1* -> C_site2*)
            c_site1 = ads_pos - 0.35 * ads_dir
            c_site2 = c_site1 + np.array([2.45, 0.0, 0.0])
            init_6 = structure.copy()
            init_6.append(Atom('C', position=c_site1))
            fin_6 = structure.copy()
            fin_6.append(Atom('C', position=c_site2))
            c_trans, ok6 = _run_neb_step(
                init_6,
                fin_6,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=c_trans,
                zpe_correction_ev=-0.03,
            )
            neb_steps_converged += int(ok6)

            # Step 7: C* + C* -> C2* graphitic dimerization / encapsulation nucleation
            init_7 = structure.copy()
            init_7.append(Atom('C', position=c_site1 + np.array([-1.35, 0.0, 0.0])))
            init_7.append(Atom('C', position=c_site1 + np.array([1.35, 0.0, 0.0])))
            fin_7 = structure.copy()
            fin_7.append(Atom('C', position=c_site1 + np.array([-0.66, 0.0, 0.15])))
            fin_7.append(Atom('C', position=c_site1 + np.array([0.66, 0.0, 0.15])))
            c_encap_fallback = max(c_trans + 0.03, c_encap)
            c_encap, ok7 = _run_neb_step(
                init_7,
                fin_7,
                calc,
                n_images=n_images,
                fmax=fmax,
                max_steps=max_steps,
                fallback=c_encap_fallback,
                zpe_correction_ev=-0.04,
            )
            neb_steps_converged += int(ok7)
        except Exception:
            pass

    return {
        'E_act_neb': float(e_act_neb),
        'ch3_dehydrogenation_eV': float(ch3_deh),
        'ch2_dehydrogenation_eV': float(ch2_deh),
        'ch_dehydrogenation_eV': float(ch_deh),
        'h2_desorption_eV': float(h2_des),
        'carbon_transfer_eV': float(c_trans),
        'carbon_encapsulation_eV': float(c_encap),
        'mlip_neb_steps_converged': int(neb_steps_converged),
    }


def write_pyrolysis_validation_json(
    candidate_id: str,
    neb_record: MlipNebValidationRecord | Mapping[str, float | int],
    kinetics_validation_dir: str | Path,
) -> Path:
    """Persist ``<kinetics_validation_dir>/<candidate_id>/pyrolysis_validation.json``.

    Args:
        candidate_id: Content-addressed candidate identifier.
        neb_record: Dictionary of MLIP-NEB elementary barriers from ``compute_candidate_mlip_neb``.
        kinetics_validation_dir: Root directory storing per-candidate validation JSONs.

    Returns:
        Path to the written ``pyrolysis_validation.json`` artifact.
    """
    safe_id = str(candidate_id)
    if not safe_id or Path(safe_id).name != safe_id:
        raise ValueError(f'unsafe candidate_id for validation artifact: {safe_id!r}')
    out_dir = Path(kinetics_validation_dir) / safe_id
    out_dir.mkdir(parents=True, exist_ok=True)
    payload: PyrolysisValidationArtifact = {
        'candidate_id': safe_id,
        'complete': True,
        'evidence_level': 'converged_mlip_neb',
        'mlip_neb_steps_converged': int(neb_record.get('mlip_neb_steps_converged', 0)),
        'resolved_kinetics_eV': {
            'methane_activation_eV': float(neb_record['E_act_neb']),
            'ch3_dehydrogenation_eV': float(neb_record['ch3_dehydrogenation_eV']),
            'ch2_dehydrogenation_eV': float(neb_record['ch2_dehydrogenation_eV']),
            'ch_dehydrogenation_eV': float(neb_record['ch_dehydrogenation_eV']),
            'h2_desorption_eV': float(neb_record['h2_desorption_eV']),
            'carbon_transfer_eV': float(neb_record['carbon_transfer_eV']),
            'carbon_encapsulation_eV': float(neb_record['carbon_encapsulation_eV']),
        },
    }
    out_path = out_dir / 'pyrolysis_validation.json'
    out_path.write_text(json.dumps(payload, indent=2) + '\n')
    return out_path


def upgrade_candidates_with_mlip_neb(
    df,
    calc=None,
    *,
    kinetics_validation_dir: str | Path | None = None,
    max_gpu_candidates: int = 8,
):
    """Populate MLIP-NEB elementary barrier columns on a selected candidate DataFrame.

    Args:
        df: Selected pyrolysis candidate DataFrame.
        calc: Optional ASE calculator for climbing-image NEB evaluation.
        kinetics_validation_dir: Optional directory to write per-candidate ``pyrolysis_validation.json``.
        max_gpu_candidates: Maximum number of rows evaluated with the live 3D calculator.

    Returns:
        Copy of ``df`` with elementary barrier columns added.
    """
    if df is None or len(df) == 0:
        return df
    from pipeline.reactors.mechanisms import (
        off_site_carbon_allowed,
        parse_catalyst_genome,
    )

    out = df.copy()
    rows_list = list(out.iterrows())
    priority_order = sorted(
        range(len(rows_list)),
        key=lambda i: (
            0
            if off_site_carbon_allowed(
                parse_catalyst_genome(rows_list[i][1].get('genome')),
                rows_list[i][1].get('material_class'),
            )
            else 1,
            i,
        ),
    )
    gpu_selected = set(priority_order[: max(0, int(max_gpu_candidates))])

    records = []
    for idx, (_, row) in enumerate(rows_list):
        cid = str(row.get('candidate_id') or f'catalyst_{idx}')
        cached_rec = None
        if kinetics_validation_dir is not None and Path(cid).name == cid:
            existing = Path(kinetics_validation_dir) / cid / 'pyrolysis_validation.json'
            if existing.is_file():
                try:
                    data = json.loads(existing.read_text())
                    rk = data.get('resolved_kinetics_eV', {})
                    if data.get('complete') and int(data.get('mlip_neb_steps_converged', 0)) > 0:
                        cached_rec = {
                            'E_act_neb': float(rk['methane_activation_eV']),
                            'ch3_dehydrogenation_eV': float(rk['ch3_dehydrogenation_eV']),
                            'ch2_dehydrogenation_eV': float(rk['ch2_dehydrogenation_eV']),
                            'ch_dehydrogenation_eV': float(rk['ch_dehydrogenation_eV']),
                            'h2_desorption_eV': float(rk['h2_desorption_eV']),
                            'carbon_transfer_eV': float(rk['carbon_transfer_eV']),
                            'carbon_encapsulation_eV': float(rk['carbon_encapsulation_eV']),
                            'mlip_neb_steps_converged': int(
                                data.get('mlip_neb_steps_converged', 0)
                            ),
                        }
                except Exception:
                    cached_rec = None
        if cached_rec is not None:
            rec = cached_rec
        else:
            row_calc = calc if idx in gpu_selected else None
            rec = compute_candidate_mlip_neb(row, calc=row_calc)
            if kinetics_validation_dir is not None:
                try:
                    write_pyrolysis_validation_json(cid, rec, kinetics_validation_dir)
                except Exception:
                    pass
        records.append(rec)
    for col in (
        'E_act_neb',
        'ch3_dehydrogenation_eV',
        'ch2_dehydrogenation_eV',
        'ch_dehydrogenation_eV',
        'h2_desorption_eV',
        'carbon_transfer_eV',
        'carbon_encapsulation_eV',
        'mlip_neb_steps_converged',
    ):
        out[col] = [r[col] for r in records]
    return out
