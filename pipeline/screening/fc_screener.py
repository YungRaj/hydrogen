#!/usr/bin/env python3
# Fuel-cell candidate screening.
"""
Multi-GPU Meta eSen-SM Screening for ORR Fuel Cell Cathode Catalysts.

For each catalyst candidate from the 21.1B encoded design space, computes:
  1. ΔG_OH* — hydroxyl adsorption free energy
  2. ΔG_O*  — oxygen adsorption free energy
  3. ΔG_OOH* — peroxyl adsorption free energy
  4. ORR overpotential (4-electron CHE method)
  5. Fenton susceptibility index (radical degradation risk)
  6. Dissolution stability estimate

Uses the SAME design space and structure generators as the methane
pyrolysis screener, but evaluates ORR-specific descriptors.

Parallelized across all available GPUs with 2 workers per device.
"""

import numpy as np
import multiprocessing as mp
from typing import List

from ase import Atom
from pipeline.screening.relaxation import relax_with_record, require_relaxation
from pipeline.screening.protocols import ORR_PROTOCOL

from pipeline.utils import (
    FUEL_CELL_DIR,
    setup_logger,
    orr_overpotential,
    abundance_cost_penalty,
    check_element_safety,
    is_valid_for_application,
)
from pipeline.screening.surface_screener import (
    generate_structure,
    screening_adsorption_site,
)
from pipeline.search.design_space import extract_elements as _extract_elements

logger = setup_logger('fc_screener', 'fuel_cell/fc_screening.log')

SCREENING_PROTOCOL_ID = ORR_PROTOCOL.protocol_id


# ═══════════════════════════════════════════════════════════════════════════════
# ORR REFERENCE ENERGIES
# ═══════════════════════════════════════════════════════════════════════════════

# Zero-point energy and entropy corrections for ORR intermediates (eV)
# From standard DFT+CHE literature (Nørskov et al., J. Phys. Chem. B, 2004)
ZPE_CORRECTIONS = {
    'OH': 0.35,  # ZPE(OH*) - 0.5*ZPE(H2O)
    'O': 0.05,  # ZPE(O*)
    'OOH': 0.40,  # ZPE(OOH*)
}
TS_CORRECTIONS = {
    'OH': -0.07,
    'O': 0.00,
    'OOH': -0.10,
}

# Elements with known Fenton reactivity (radical generation in acid)
FENTON_RISK = {
    'Fe': 3,
    'Cu': 2,
    'Co': 1,
    'Mn': 1,
    'Cr': 1,
    'V': 1,
    'Ti': 0,
    'Ni': 0,
    'Zn': 0,
    'Mo': 0,
}

ORR_MAX_BINDING_ABS_EV = 4.0
ORR_MAX_REPORTED_OVERPOTENTIAL_V = 3.0


def _adsorbate_geometry_valid(slab, surface_count: int, kind: str) -> tuple[bool, str]:
    """Reject dissociated/embedded relaxed adsorbates before scoring."""
    positions = slab.positions
    ads = positions[surface_count:]
    if len(ads) == 0 or not np.isfinite(ads).all():
        return False, 'adsorbate_nonfinite'
    from scipy.spatial.distance import pdist
    if len(ads) > 1 and float(np.min(pdist(ads))) < 0.65:
        return False, 'adsorbate_embedded_or_collapsed'
    if kind == 'oh' and len(ads) == 2:
        distances = float(np.linalg.norm(ads[1] - ads[0]))
        if not 0.65 <= distances <= 1.25:
            return False, 'oh_bond_broken'
    if kind == 'ooh' and len(ads) == 3:
        oo = float(np.linalg.norm(ads[1] - ads[0]))
        oh = float(np.linalg.norm(ads[2] - ads[1]))
        if not 0.95 <= oo <= 1.8 or not 0.65 <= oh <= 1.25:
            return False, 'ooh_bond_broken'
    # An adsorbate buried in the slab is not a meaningful surface state.
    if float(np.min(np.linalg.norm(ads[:, None, :] - positions[:surface_count][None, :, :], axis=2))) < 0.65:
        return False, 'adsorbate_embedded'
    return True, ''


def compute_water_ref(calc) -> float:
    """Compute H₂O reference energy.

    Args:
        calc: Atomic calculator used for the evaluation.

    Returns:
        Computed `float` value in the units documented above.
    """
    from ase.build import molecule

    h2o = molecule('H2O')
    h2o.set_cell([10, 10, 10])
    h2o.center()
    h2o.pbc = True
    h2o.calc = calc
    budget = ORR_PROTOCOL.reference
    record = relax_with_record(h2o, 'reference_h2o', budget.fmax_eV_A, budget.steps)
    if not record['relax_reference_h2o_converged']:
        raise RuntimeError('H2O reference relaxation did not converge')
    return h2o.get_potential_energy()


def compute_h2_ref(calc) -> float:
    """Compute H₂ reference energy.

    Args:
        calc: Atomic calculator used for the evaluation.

    Returns:
        Computed `float` value in the units documented above.
    """
    from ase.build import molecule

    h2 = molecule('H2')
    h2.set_cell([10, 10, 10])
    h2.center()
    h2.pbc = True
    h2.calc = calc
    budget = ORR_PROTOCOL.reference
    record = relax_with_record(h2, 'reference_h2', budget.fmax_eV_A, budget.steps)
    if not record['relax_reference_h2_converged']:
        raise RuntimeError('H2 reference relaxation did not converge')
    return h2.get_potential_energy()


# ═══════════════════════════════════════════════════════════════════════════════
# ORR EVALUATION
# ═══════════════════════════════════════════════════════════════════════════════


def evaluate_orr_candidate(genome: tuple, calc, e_h2o: float, e_h2: float) -> dict:
    """
        Evaluate a single catalyst candidate for ORR activity using Meta eSen-SM.

        Computes adsorption free energies for OH*, O*, OOH* intermediates
        and derives the theoretical ORR overpotential via the CHE method.

    Args:
        genome: Encoded catalyst composition and structural configuration.
        calc: Atomic calculator used for the evaluation.
        e_h2o: E h2o used by this operation.
        e_h2: E h2 used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    mat_class = genome[0]
    result = {
        'genome': str(genome),
        'material_class': mat_class,
        'valid': False,
        'screening_protocol': SCREENING_PROTOCOL_ID,
    }

    try:
        structure, active_idx, _ = generate_structure(genome)
        structure.pbc = True  # eSen requires PBC set to True in all dimensions

        # 1. Relax clean surface
        structure.calc = calc
        budget = ORR_PROTOCOL.clean
        if not require_relaxation(
            result, structure, 'clean', budget.fmax_eV_A, budget.steps
        ):
            return result
        e_clean = structure.get_potential_energy()
        result['e_clean'] = e_clean

        # 2. OH* adsorption
        slab_oh = structure.copy()
        oh_pos, ads_direction = screening_adsorption_site(
            structure, active_idx, mat_class, 1.9
        )
        slab_oh.append(Atom('O', position=oh_pos))
        slab_oh.append(Atom('H', position=oh_pos + 0.97 * ads_direction))
        slab_oh.calc = calc
        budget = ORR_PROTOCOL.adsorbate
        if not require_relaxation(
            result, slab_oh, 'oh', budget.fmax_eV_A, budget.steps
        ):
            return result
        valid_geometry, geometry_reason = _adsorbate_geometry_valid(slab_oh, len(structure), 'oh')
        if not valid_geometry:
            result.update({'candidate_disposition': 'validation_required', 'error': geometry_reason,
                           'needs_dft_validation': True})
            return result
        e_oh = slab_oh.get_potential_energy()
        # ΔG_OH* = E(slab+OH) - E(slab) - (E(H2O) - 0.5*E(H2)) + ZPE + TS
        dG_OH = (
            (e_oh - e_clean)
            - (e_h2o - 0.5 * e_h2)
            + ZPE_CORRECTIONS['OH']
            + TS_CORRECTIONS['OH']
        )
        result['dG_OH_eV'] = float(dG_OH)

        # 3. O* adsorption
        slab_o = structure.copy()
        o_pos, _ = screening_adsorption_site(structure, active_idx, mat_class, 1.7)
        slab_o.append(Atom('O', position=o_pos))
        slab_o.calc = calc
        if not require_relaxation(result, slab_o, 'o', budget.fmax_eV_A, budget.steps):
            return result
        e_o = slab_o.get_potential_energy()
        # ΔG_O* = E(slab+O) - E(slab) - (E(H2O) - E(H2)) + ZPE + TS
        dG_O = (
            (e_o - e_clean)
            - (e_h2o - e_h2)
            + ZPE_CORRECTIONS['O']
            + TS_CORRECTIONS['O']
        )
        result['dG_O_eV'] = float(dG_O)

        # 4. OOH* adsorption
        slab_ooh = structure.copy()
        o1_pos, ads_direction = screening_adsorption_site(
            structure, active_idx, mat_class, 1.9
        )
        o2_pos = o1_pos + np.array([1.2, 0.0, 0.6])
        if ads_direction[2] < 0:
            o2_pos[2] = o1_pos[2] - 0.6
        h_pos = o2_pos + 0.97 * ads_direction
        slab_ooh.append(Atom('O', position=o1_pos))
        slab_ooh.append(Atom('O', position=o2_pos))
        slab_ooh.append(Atom('H', position=h_pos))
        slab_ooh.calc = calc
        if not require_relaxation(
            result, slab_ooh, 'ooh', budget.fmax_eV_A, budget.steps
        ):
            return result
        valid_geometry, geometry_reason = _adsorbate_geometry_valid(slab_ooh, len(structure), 'ooh')
        if not valid_geometry:
            result.update({'candidate_disposition': 'validation_required', 'error': geometry_reason,
                           'needs_dft_validation': True})
            return result
        e_ooh = slab_ooh.get_potential_energy()
        # ΔG_OOH* = E(slab+OOH) - E(slab) - (2*E(H2O) - 1.5*E(H2)) + ZPE + TS
        dG_OOH = (
            (e_ooh - e_clean)
            - (2 * e_h2o - 1.5 * e_h2)
            + ZPE_CORRECTIONS['OOH']
            + TS_CORRECTIONS['OOH']
        )
        result['dG_OOH_eV'] = float(dG_OOH)

        # 5. ORR overpotential (4e⁻ CHE method)
        eta, rds = orr_overpotential(dG_OH, dG_O, dG_OOH)
        result['orr_overpotential_V'] = float(eta)
        result['rate_determining_step'] = rds

        # 6. Fenton susceptibility
        elements = _extract_elements(genome)
        fenton_score = sum(FENTON_RISK.get(e, 0) for e in elements)
        result['fenton_stability'] = max(0, 10 - fenton_score)

        # 7. Cost
        result['cost_penalty'] = abundance_cost_penalty(elements)

        # 7b. Safety check
        is_safe, safety_reason = check_element_safety(elements)
        if not is_safe:
            result['valid'] = False
            result['candidate_disposition'] = 'hard_excluded'
            result['error'] = safety_reason
            return result

        # 7c. Application feasibility flag
        result['fc_viable'] = is_valid_for_application(mat_class, 'fuel_cell')

        # 8. Stability estimate (binding strength of active metal)
        result['binding_strength'] = float(abs(dG_OH) + abs(dG_O))
        gc_corr = compute_grand_canonical_orr_corrections(
            dG_OH, dG_O, dG_OOH, genome=genome
        )
        result.update(gc_corr)

        # 9. Physical sanity filters
        SANE_LIMIT = 10.0  # eV
        for key in ('dG_OH_eV', 'dG_O_eV', 'dG_OOH_eV'):
            val = result.get(key, 0)
            if abs(val) > SANE_LIMIT:
                result['valid'] = False
                result['candidate_disposition'] = 'validation_required'
                result['error'] = f'Unphysical {key}={val:.2f} eV'
                result['needs_dft_validation'] = True
                return result

        # Never turn an out-of-domain prediction into a plausible-looking score.
        eta = float(result['orr_overpotential_V'])
        if eta > ORR_MAX_REPORTED_OVERPOTENTIAL_V:
            result.update({'valid': False, 'candidate_disposition': 'validation_required',
                           'error': 'orr_overpotential_capped_at_3V',
                           'overpotential_capped': True, 'needs_dft_validation': True})
            return result
        result['overpotential_capped'] = False
        for key in ('dG_OH_eV', 'dG_O_eV', 'dG_OOH_eV'):
            if abs(float(result[key])) > ORR_MAX_BINDING_ABS_EV:
                result.update({'valid': False, 'candidate_disposition': 'validation_required',
                               'error': f'{key}_outside_4eV_domain', 'needs_dft_validation': True})
                return result

        # 10. OOD confidence — how much we trust this prediction
        from pipeline.screening.ood import compute_model_confidence

        conf = compute_model_confidence(genome, elements)
        result['model_confidence'] = float(conf)
        result['needs_dft_validation'] = conf < 0.5

        result['valid'] = True
        result['candidate_disposition'] = 'quantitative_screening'

    except Exception as e:
        result['error'] = str(e)[:200]
        result['candidate_disposition'] = 'validation_required'
        result['needs_dft_validation'] = True

    return result


STANDARD_REDUCTION_POTENTIALS_V: dict[str, float] = {
    'Au': 1.520, 'Pt': 1.188, 'Ir': 1.156, 'Pd': 0.915, 'Os': 0.840,
    'Ag': 0.800, 'Rh': 0.758, 'Ru': 0.455, 'Cu': 0.340, 'Bi': 0.317,
    'Re': 0.300, 'Sb': 0.212, 'Te': 0.530, 'W': -0.119, 'Pb': -0.126,
    'Sn': -0.138, 'Mo': -0.200, 'Ni': -0.257, 'Co': -0.277, 'In': -0.338,
    'Fe': -0.440, 'Ga': -0.529, 'Ta': -0.600, 'Cr': -0.740, 'Zn': -0.763,
    'Nb': -1.100, 'V': -1.130, 'Mn': -1.185, 'Zr': -1.450, 'Hf': -1.550,
    'Ti': -1.630, 'Al': -1.676, 'Sc': -2.077, 'Ce': -2.336, 'Y': -2.372,
    'Mg': -2.372, 'La': -2.379, 'Ca': -2.868, 'Sr': -2.899,
}

# Surface dipole moments (e*A) and excess polarizabilities (e*A^2/V) for ORR intermediates
_ADSORBATE_DIPOLE_E_A = {'OH*': -0.35, 'O*': 0.18, 'OOH*': -0.42}
_ADSORBATE_ALPHA_E_A2_V = {'OH*': 0.12, 'O*': 0.06, 'OOH*': 0.18}
# Implicit H-bond water stabilization (eV) on hydrophilic metallic/oxide vs carbon hosts
_IMPLICIT_HBOND_EV = {'OH*': -0.30, 'O*': -0.05, 'OOH*': -0.25}


# Coordination/lattice Pourbaix stabilization shift (V) by material class across all 14 classes
POURBAIX_CLASS_COORD_SHIFT_V: dict[str, float] = {
    'MoltenMetal': 0.20,
    'SolidCatalyst': 0.40,
    'SAC': 1.05,
    'DAC': 1.05,
    'MOF': 1.00,
    'COF': 1.05,
    'Perovskite': 1.10,
    'MetalHydride': 0.25,
    'MAXPhase': 0.95,
    'HEA': 0.55,
    'Spinel': 1.10,
    'MXene': 0.95,
    'SAA': 0.55,
    'MetalFreeCarbon': 1.25,
}

# Interfacial H-bond water stabilization scale factor across all 14 classes
HBOND_CLASS_SCALE: dict[str, float] = {
    'MoltenMetal': 0.85,
    'SolidCatalyst': 1.00,
    'SAC': 0.75,
    'DAC': 0.75,
    'MOF': 0.75,
    'COF': 0.75,
    'Perovskite': 1.05,
    'MetalHydride': 0.90,
    'MAXPhase': 0.95,
    'HEA': 1.00,
    'Spinel': 1.05,
    'MXene': 1.00,
    'SAA': 1.00,
    'MetalFreeCarbon': 0.75,
}


from pipeline.data_models.fuel_cells import (
    GrandCanonicalORRRecord,
    PourbaixStabilityRecord,
)


def compute_pourbaix_stability_margin(
    genome,
    *,
    U_V_vs_RHE: float = 0.80,
    pH: float = 0.0,
) -> PourbaixStabilityRecord:
    """Estimate Pourbaix anodic dissolution onset and stability margin at (U, pH).

    Args:
        genome: Encoded catalyst genome tuple or string representation.
        U_V_vs_RHE: Operating electrode potential in V vs RHE.
        pH: Electrolyte pH.

    Returns:
        Dictionary with dissolution onset potential, margin, and acid stability flag.
    """
    from pipeline.reactors.mechanisms import parse_catalyst_genome

    parsed = parse_catalyst_genome(genome) if not isinstance(genome, tuple) else genome
    if not parsed:
        return {
            'pourbaix_dissolution_onset_V': 0.80,
            'pourbaix_dissolution_margin_V': 0.0,
            'pourbaix_acid_stable': True,
        }
    mat_class = str(parsed[0])
    elements = [el for el in _extract_elements(parsed) if el in STANDARD_REDUCTION_POTENTIALS_V]
    if not elements:
        # MetalFreeCarbon / covalent non-metal frameworks resist anodic metal dissolution
        onset = POURBAIX_CLASS_COORD_SHIFT_V.get(mat_class, 1.25)
    else:
        # Active site is governed by the primary metal (first extracted metal)
        e0_primary = STANDARD_REDUCTION_POTENTIALS_V[elements[0]]
        e0_mean = float(np.mean([STANDARD_REDUCTION_POTENTIALS_V[el] for el in elements]))
        base_e0 = 0.65 * e0_primary + 0.35 * e0_mean
        coord_shift = POURBAIX_CLASS_COORD_SHIFT_V.get(mat_class, 0.35)
        if mat_class == 'SolidCatalyst' and len(elements) == 1:
            coord_shift = 0.30
        onset = base_e0 + coord_shift + 0.05916 * float(pH)
    margin = float(onset - U_V_vs_RHE)
    return {
        'pourbaix_dissolution_onset_V': float(onset),
        'pourbaix_dissolution_margin_V': margin,
        'pourbaix_acid_stable': bool(margin >= -0.05),
    }


def compute_grand_canonical_orr_corrections(
    dG_OH: float,
    dG_O: float,
    dG_OOH: float,
    *,
    genome=None,
    U_V_vs_RHE: float = 0.80,
    U_pzc_V: float = 0.30,
    helmholtz_thickness_A: float = 3.0,
    pH: float = 0.0,
) -> GrandCanonicalORRRecord:
    """Compute constant-potential double-layer solvation + electric-field corrections.

    Args:
        dG_OH: Unsolvated OH* adsorption free energy (eV).
        dG_O: Unsolvated O* adsorption free energy (eV).
        dG_OOH: Unsolvated OOH* adsorption free energy (eV).
        genome: Optional catalyst genome tuple or string.
        U_V_vs_RHE: Operating electrode potential in V vs RHE.
        U_pzc_V: Potential of zero charge in V vs RHE.
        helmholtz_thickness_A: Effective compact double-layer thickness in Angstrom.
        pH: Electrolyte pH.

    Returns:
        Dictionary of grand-canonical solvated free energies, overpotential, and Pourbaix stability.
    """
    from pipeline.reactors.mechanisms import parse_catalyst_genome

    parsed = parse_catalyst_genome(genome) if (genome is not None and not isinstance(genome, tuple)) else genome
    mat_class = str(parsed[0]) if parsed else 'SolidCatalyst'
    hbond_scale = HBOND_CLASS_SCALE.get(mat_class, 1.0)
    e_field_V_A = (float(U_V_vs_RHE) - float(U_pzc_V)) / max(1.0, float(helmholtz_thickness_A))

    def _shift(ads: str) -> float:
        mu = _ADSORBATE_DIPOLE_E_A[ads]
        alpha = _ADSORBATE_ALPHA_E_A2_V[ads]
        solv = hbond_scale * _IMPLICIT_HBOND_EV[ads]
        field_term = -mu * e_field_V_A - 0.5 * alpha * (e_field_V_A ** 2)
        return float(solv + field_term)

    dg_oh_gc = float(dG_OH) + _shift('OH*')
    dg_o_gc = float(dG_O) + _shift('O*')
    dg_ooh_gc = float(dG_OOH) + _shift('OOH*')
    eta_gc, rds_gc = orr_overpotential(dg_oh_gc, dg_o_gc, dg_ooh_gc)
    pourbaix = compute_pourbaix_stability_margin(
        parsed, U_V_vs_RHE=U_V_vs_RHE, pH=pH
    )
    return {
        'dG_OH_gc_eV': float(dg_oh_gc),
        'dG_O_gc_eV': float(dg_o_gc),
        'dG_OOH_gc_eV': float(dg_ooh_gc),
        'orr_overpotential_gc_V': float(eta_gc),
        'limiting_potential_gc_V': float(1.23 - eta_gc),
        'rate_determining_step_gc': str(rds_gc),
        **pourbaix,
    }


def upgrade_orr_candidates_with_solvation_pourbaix(
    df,
    *,
    U_V_vs_RHE: float = 0.80,
    pH: float = 0.0,
):
    """Populate grand-canonical solvation and Pourbaix stability columns on an ORR DataFrame.

    Args:
        df: Selected ORR candidate DataFrame.
        U_V_vs_RHE: Operating electrode potential in V vs RHE.
        pH: Electrolyte pH.

    Returns:
        Copy of ``df`` with grand-canonical solvation and Pourbaix columns added.
    """
    if df is None or len(df) == 0:
        return df
    out = df.copy()
    records = []
    for _, row in out.iterrows():
        dg_oh = float(row.get('dG_OH_eV', 0.8))
        dg_o = float(row.get('dG_O_eV', 1.8))
        dg_ooh = float(row.get('dG_OOH_eV', 3.8))
        records.append(
            compute_grand_canonical_orr_corrections(
                dg_oh, dg_o, dg_ooh,
                genome=row.get('genome'),
                U_V_vs_RHE=U_V_vs_RHE,
                pH=pH,
            )
        )
    for col in (
        'dG_OH_gc_eV',
        'dG_O_gc_eV',
        'dG_OOH_gc_eV',
        'orr_overpotential_gc_V',
        'limiting_potential_gc_V',
        'rate_determining_step_gc',
        'pourbaix_dissolution_onset_V',
        'pourbaix_dissolution_margin_V',
        'pourbaix_acid_stable',
    ):
        out[col] = [r[col] for r in records]
    return out


# ═══════════════════════════════════════════════════════════════════════════════
# MULTI-GPU WORKER
# ═══════════════════════════════════════════════════════════════════════════════


def orr_worker(
    worker_id: int,
    gpu_id: int,
    gpu_uuid: str,
    task_queue: mp.Queue,
    result_queue: mp.Queue,
    stop_event,
    candidate_threads: int = 1,
    batched: bool = False,
):
    """Worker process: loads Meta eSen on assigned GPU, evaluates ORR candidates.

    Args:
        worker_id: Worker id used by this operation.
        gpu_id: Gpu id used by this operation.
        gpu_uuid: Gpu uuid used by this operation.
        task_queue: Task queue used by this operation.
        result_queue: Result queue used by this operation.
        stop_event: Stop event used by this operation.
        candidate_threads: Candidate threads used by this operation.
        batched: Whether to enable batched.
    """
    from pipeline.screening.gpu_executor import run_gpu_worker

    def prepare_evaluator(calculator):
        water_reference = compute_water_ref(calculator)
        hydrogen_reference = compute_h2_ref(calculator)
        return lambda genome, thread_calculator: evaluate_orr_candidate(
            genome,
            thread_calculator,
            water_reference,
            hydrogen_reference,
        )

    run_gpu_worker(
        worker_id,
        gpu_id,
        gpu_uuid,
        task_queue,
        result_queue,
        stop_event,
        candidate_threads,
        batched,
        prepare_evaluator,
        SCREENING_PROTOCOL_ID,
        logger,
        'ORR worker',
    )


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN ORR SCREENING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════


def run_orr_screening(
    genomes: List[tuple],
    db_filename: str = "fc_screening.csv",
    workers_per_gpu: int = 2,
    engine: str = 'batched',
) -> 'pd.DataFrame':
    """
        Run parallel Meta eSen-SM ORR screening on a list of catalyst genomes.

        Same interface as surface_screener.run_screening but evaluates
        OH*, O*, OOH* for fuel cell cathode applications.

    Args:
        genomes: Sequence of encoded catalyst candidates.
        db_filename: Db filename used by this operation.
        workers_per_gpu: Workers per gpu used by this operation.
        engine: Engine used by this operation.

    Returns:
        Computed `'pd.DataFrame'` result.
    """
    from pipeline.screening.gpu_executor import ScreeningRunSpec, run_gpu_screening

    df = run_gpu_screening(
        genomes,
        db_filename,
        workers_per_gpu,
        engine,
        orr_worker,
        logger,
        ScreeningRunSpec(
            banner='META ESEN-SM ORR CATHODE CATALYST SCREENING',
            application='fuel_cell_orr',
            manifest_path=FUEL_CELL_DIR / 'orr_worker_health.json',
            output_subdir='fuel_cell',
            start_message='ORR screening {count} candidates...',
            completion_label='ORR screening',
            progress_noun='cand',
            protocol_id=SCREENING_PROTOCOL_ID,
            cache_path=FUEL_CELL_DIR / 'esen_result_cache.sqlite',
        ),
    )

    valid_df = df[df['valid'] == True]
    if len(valid_df) > 0:
        logger.info(f"  Valid: {len(valid_df)}/{len(df)}")
        logger.info(
            f"  Best overpotential: {valid_df['orr_overpotential_V'].min():.4f} V"
        )
        logger.info(
            f"  ΔG_OH range: [{valid_df['dG_OH_eV'].min():.3f}, {valid_df['dG_OH_eV'].max():.3f}] eV"
        )

    return df


if __name__ == '__main__':
    from pipeline.search.indexed_space import deterministic_tree_probes

    pop = deterministic_tree_probes(20)
    df = run_orr_screening(pop, db_filename="test_orr_screening.csv", workers_per_gpu=2)
    print(f"\nORR screening complete: {len(df)} candidates")
