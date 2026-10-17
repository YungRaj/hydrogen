#!/usr/bin/env python3
# Turquoise-hydrogen candidate screening.
"""
Multi-GPU Parallelized Meta eSen-SM Catalyst Screening Engine.

For each catalyst candidate, computes:
  1. Relaxed surface energy
  2. H* adsorption energy (methane activation descriptor)
  3. CH₃* adsorption energy (dissociation intermediate)
  4. C* adsorption energy (coking resistance descriptor)
  5. Reaction energy: CH₄ → CH₃* + H*
  6. Segregation / binding stability energy
  7. NEB transition state barrier (top candidates only)

3D ASE structure generation lives in ``pipeline.screening.structure_builders``
and 5-step MLIP-NEB barrier refinement lives in ``pipeline.screening.mlip_neb``.
Both are re-exported here for backwards compatibility.
"""

import json
import multiprocessing as mp
import os
from typing import Dict, List

import numpy as np
from ase import Atom
from ase.build import molecule

from pipeline.screening.mlip_neb import (
    _run_neb_step,
    compute_candidate_mlip_neb,
    upgrade_candidates_with_mlip_neb,
    write_pyrolysis_validation_json,
)
from pipeline.screening.protocols import PYROLYSIS_PROTOCOL
from pipeline.screening.relaxation import relax_with_record, require_relaxation
from pipeline.screening.structure_builders import (
    LATTICE_CONSTANTS,
    _add_surface_markers,
    _append_marker,
    _generate_hydride_slab,
    _generate_perovskite_slab,
    _generate_spinel_slab,
    _label_fraction,
    _representative_element,
    _stable_seed,
    adsorption_site,
    generate_alloy_slab,
    generate_porphyrin_cluster,
    generate_structure,
    screening_adsorption_site,
)
from pipeline.search.design_space import extract_elements as _extract_elements
from pipeline.utils import (
    SCREENING_DIR,
    abundance_cost_penalty,
    arrhenius_rate,
    bep_activation_energy,
    check_element_safety,
    setup_logger,
)

logger = setup_logger('surface_screener', 'screening/surface_screening.log')

SCREENING_PROTOCOL_ID = PYROLYSIS_PROTOCOL.protocol_id


# ═══════════════════════════════════════════════════════════════════════════════
# MACE / ESEN EVALUATION WORKER
# ═══════════════════════════════════════════════════════════════════════════════


def compute_reference_energies(calc) -> Dict[str, float]:
    """Compute gas-phase reference energies using Meta eSen-SM.

    Args:
        calc: Atomic calculator used for the evaluation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    refs = {}

    # H₂
    h2 = molecule('H2')
    h2.set_cell([10, 10, 10])
    h2.center()
    h2.pbc = True
    h2.calc = calc
    budget = PYROLYSIS_PROTOCOL.reference
    record = relax_with_record(h2, 'reference_h2', budget.fmax_eV_A, budget.steps)
    if not record['relax_reference_h2_converged']:
        raise RuntimeError('H2 reference relaxation did not converge')
    refs['H2'] = h2.get_potential_energy()

    # CH₄
    ch4 = molecule('CH4')
    ch4.set_cell([10, 10, 10])
    ch4.center()
    ch4.pbc = True
    ch4.calc = calc
    record = relax_with_record(ch4, 'reference_ch4', budget.fmax_eV_A, budget.steps)
    if not record['relax_reference_ch4_converged']:
        raise RuntimeError('CH4 reference relaxation did not converge')
    refs['CH4'] = ch4.get_potential_energy()

    # CH₃ radical
    ch3 = molecule('CH3')
    ch3.set_cell([10, 10, 10])
    ch3.center()
    ch3.pbc = True
    ch3.calc = calc
    record = relax_with_record(ch3, 'reference_ch3', budget.fmax_eV_A, budget.steps)
    if not record['relax_reference_ch3_converged']:
        raise RuntimeError('CH3 reference relaxation did not converge')
    refs['CH3'] = ch3.get_potential_energy()

    # Isolated C reference calculated thermodynamically: CH4 - 2*H2
    # This avoids zero-edge errors for single-atom systems in fairchem/eSen.
    refs['C'] = refs['CH4'] - 2.0 * refs['H2']

    return refs


def _evaluate_candidate_replica(
    genome: tuple, calc, refs: dict, structure_seed: int
) -> dict:
    """Evaluate a single catalyst candidate using Meta eSen-SM.

    Args:
        genome: Encoded catalyst composition and structural configuration.
        calc: Atomic calculator used for the evaluation.
        refs: Mapping supplying refs.
        structure_seed: Deterministic structural replica index.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    mat_class = genome[0]
    result = {
        'genome': str(genome),
        'material_class': mat_class,
        'valid': False,
        'screening_protocol': SCREENING_PROTOCOL_ID,
        'structure_seed': structure_seed,
    }

    try:
        structure, active_idx, _ = generate_structure(genome, structure_seed)
        structure.pbc = True

        # 1. Relax clean surface/cluster
        structure.calc = calc
        budget = PYROLYSIS_PROTOCOL.clean
        if not require_relaxation(
            result, structure, 'clean', budget.fmax_eV_A, budget.steps
        ):
            return result
        e_clean = structure.get_potential_energy()
        result['e_clean'] = e_clean

        ads_pos, ads_direction = screening_adsorption_site(
            structure, active_idx, mat_class, 1.8
        )

        # 2. H* adsorption
        slab_h = structure.copy()
        h_pos = ads_pos - 0.3 * ads_direction
        slab_h.append(Atom('H', position=h_pos))
        slab_h.calc = calc
        budget = PYROLYSIS_PROTOCOL.adsorbate
        if not require_relaxation(result, slab_h, 'h', budget.fmax_eV_A, budget.steps):
            return result
        e_h = slab_h.get_potential_energy()
        dE_H = e_h - e_clean - 0.5 * refs['H2']
        result['dE_H'] = dE_H
        h_index = len(slab_h) - 1
        result['h_surface_clearance_A'] = float(
            min(
                slab_h.get_distance(h_index, index, mic=True)
                for index in range(h_index)
            )
        )
        result['h_height_above_surface_A'] = float(
            slab_h[h_index].position[2] - slab_h.positions[:h_index, 2].max()
        )

        # 3. CH₃* adsorption
        slab_ch3 = structure.copy()
        c_pos = ads_pos.copy()
        slab_ch3.append(Atom('C', position=c_pos))
        z_offset = 0.35 * ads_direction
        slab_ch3.append(Atom('H', position=c_pos + np.array([0.0, 1.02, 0.0]) + z_offset))
        slab_ch3.append(
            Atom('H', position=c_pos + np.array([-0.88, -0.51, 0.0]) + z_offset)
        )
        slab_ch3.append(
            Atom('H', position=c_pos + np.array([0.88, -0.51, 0.0]) + z_offset)
        )
        slab_ch3.calc = calc
        if not require_relaxation(
            result, slab_ch3, 'ch3', budget.fmax_eV_A, budget.steps
        ):
            return result
        e_ch3 = slab_ch3.get_potential_energy()
        dE_CH3 = e_ch3 - e_clean - refs['CH3']
        result['dE_CH3'] = dE_CH3
        carbon_index = len(structure)
        ch_bonds = [
            slab_ch3.get_distance(carbon_index, index, mic=True)
            for index in range(carbon_index + 1, carbon_index + 4)
        ]
        result['ch3_min_CH_A'] = float(min(ch_bonds))
        result['ch3_max_CH_A'] = float(max(ch_bonds))
        if result['ch3_min_CH_A'] < 0.7 or result['ch3_max_CH_A'] > 1.5:
            result['error'] = 'CH3 adsorbate lost molecular integrity'
            result['candidate_disposition'] = 'validation_required'
            result['needs_dft_validation'] = True
            return result

        # 4. C* adsorption (coking indicator)
        slab_c = structure.copy()
        c_ads_pos = ads_pos - 0.4 * ads_direction
        slab_c.append(Atom('C', position=c_ads_pos))
        slab_c.calc = calc
        if not require_relaxation(result, slab_c, 'c', budget.fmax_eV_A, budget.steps):
            return result
        e_c = slab_c.get_potential_energy()
        dE_C = e_c - e_clean - refs['C']
        result['dE_C'] = dE_C
        c_index = len(slab_c) - 1
        result['c_surface_clearance_A'] = float(
            min(
                slab_c.get_distance(c_index, index, mic=True)
                for index in range(c_index)
            )
        )
        result['c_height_above_surface_A'] = float(
            slab_c[c_index].position[2] - slab_c.positions[:c_index, 2].max()
        )

        # 5. CH₄ → CH₃* + H* reaction energy
        dE_split = dE_CH3 + dE_H + refs['CH3'] + 0.5 * refs['H2'] - refs['CH4']
        result['dE_split'] = dE_split

        # 6. Activation barrier (BEP correlation)
        E_act_raw = bep_activation_energy(
            dE_split, material_class=mat_class, clamp=False
        )
        E_act = bep_activation_energy(dE_split, material_class=mat_class)
        result['E_act_raw'] = E_act_raw
        result['E_act'] = E_act

        # 7. Coking resistance index (slab descriptor — out of scope for melts)
        from pipeline.search.scope import slab_coking_index_scope

        coking_scope = slab_coking_index_scope(genome)
        result['coking_index_scope'] = coking_scope['status']
        result['coking_index_scope_reason'] = coking_scope.get('reason')

        coking_index = dE_C - 2.0 * dE_H

        py_mode = os.environ.get('PYROLYSIS_MODE', 'thermocatalytic')
        if py_mode == 'ntec' and coking_scope['status'] == 'candidate':
            from pipeline.electrochemistry.ntec import (
                conditions_from_environment,
                ntec_assistance,
            )
            from pipeline.screening.genetic_optimizer import (
                _extract_elements_from_genome,
            )

            assistance = ntec_assistance(conditions_from_environment())
            elements = _extract_elements_from_genome(genome)
            if any(e in {'Ga', 'In', 'Sn', 'Bi'} for e in elements):
                coking_index += assistance['coking_bonus']
                result['ntec_evidence_status'] = assistance['status']
                result['ntec_assistance'] = assistance

        if coking_scope['status'] == 'out_of_scope':
            result['coking_index'] = float('nan')
            result['coking_index_raw_slab'] = coking_index
        else:
            result['coking_index'] = coking_index

        # 8. Stability metric
        if mat_class in ('SolidCatalyst', 'MoltenMetal'):
            e_seg = _compute_segregation_energy(structure, active_idx, calc, e_clean)
            result['segregation_energy'] = e_seg
        else:
            e_bind = _compute_binding_energy(structure, calc, e_clean, genome)
            result['segregation_energy'] = e_bind

        # 9. Rate constant at reference temperature
        T_ref = genome[4] if mat_class == 'MoltenMetal' else 1000.0
        A_prefactor = 1e13
        rate_k = arrhenius_rate(A_prefactor, E_act, T_ref)
        result['rate_constant'] = rate_k
        result['temperature_K'] = T_ref

        # 10. Extract element list for cost scoring
        elements = _extract_elements(genome)
        result['cost_penalty'] = abundance_cost_penalty(elements)

        # 10b. Safety check — reject toxic/radioactive elements
        is_safe, safety_reason = check_element_safety(elements)
        if not is_safe:
            result['valid'] = False
            result['candidate_disposition'] = 'hard_excluded'
            result['error'] = safety_reason
            return result

        # 10c. Encoded-phase admissibility (ADR 0001). Coverage lists stay 14-class.
        from pipeline.search.scope import phase_stable_at_application_T

        result['pyrolysis_viable'] = (
            phase_stable_at_application_T(genome)['status'] == 'candidate'
        )

        # 11. Physical sanity filters
        SANE_LIMIT = 10.0  # eV
        for key in ('dE_H', 'dE_CH3', 'dE_C', 'dE_split'):
            val = result.get(key, 0)
            if abs(val) > SANE_LIMIT:
                result['valid'] = False
                result['candidate_disposition'] = 'validation_required'
                result['error'] = f'Unphysical {key}={val:.2f} eV (|val|>{SANE_LIMIT})'
                result['needs_dft_validation'] = True
                return result

        result['E_act'] = max(0.01, min(result.get('E_act', 5.0), 5.0))
        result['E_act_censored'] = result['E_act'] in (0.01, 5.0)
        if result['E_act_censored']:
            result['needs_dft_validation'] = True

        if 'coking_index' not in result:
            result['valid'] = False
            result['candidate_disposition'] = 'validation_required'
            result['error'] = 'coking_index was never set'
            return result
        coking_val = result['coking_index']
        if coking_val == coking_val:
            result['coking_index'] = max(-20.0, min(float(coking_val), 20.0))

        # 12. OOD confidence
        from pipeline.screening.ood import compute_model_confidence

        conf = compute_model_confidence(genome, elements)
        result['model_confidence'] = float(conf)
        result['needs_dft_validation'] = (
            result.get('needs_dft_validation', False) or conf < 0.5
        )

        result['valid'] = True
        result['candidate_disposition'] = 'quantitative_screening'

    except Exception as e:
        result['error'] = str(e)[:200]
        result['candidate_disposition'] = 'validation_required'

    return result


MOLTEN_ENSEMBLE_SIZE = 5
MOLTEN_ENSEMBLE_MIN_VALID = 4
MOLTEN_ENSEMBLE_FIELDS = (
    'e_clean',
    'dE_H',
    'dE_CH3',
    'dE_C',
    'dE_split',
    'E_act_raw',
    'segregation_energy',
    'h_surface_clearance_A',
    'h_height_above_surface_A',
    'ch3_min_CH_A',
    'ch3_max_CH_A',
    'c_surface_clearance_A',
    'c_height_above_surface_A',
)


def evaluate_candidate(
    genome: tuple, calc, refs: dict, structure_seed: int | None = None
) -> dict:
    """Evaluate one candidate, averaging molten alloys over five replicas.

    Args:
        genome: Encoded catalyst composition and conditions.
        calc: Atomic calculator used for every replica.
        refs: Gas-phase reference-energy mapping.
        structure_seed: Explicit diagnostic replica, or None for production.

    Returns:
        Screening result with ensemble provenance for molten candidates.
    """
    if genome[0] != 'MoltenMetal' or structure_seed is not None:
        return _evaluate_candidate_replica(
            genome, calc, refs, 0 if structure_seed is None else structure_seed
        )

    replicas = [
        _evaluate_candidate_replica(genome, calc, refs, seed)
        for seed in range(MOLTEN_ENSEMBLE_SIZE)
    ]
    valid_replicas = [row for row in replicas if row.get('valid')]
    failures = [row for row in replicas if not row.get('valid')]
    if len(valid_replicas) < MOLTEN_ENSEMBLE_MIN_VALID:
        result = dict(failures[0])
        result['structure_seed'] = None
        result['ensemble_size'] = MOLTEN_ENSEMBLE_SIZE
        result['ensemble_valid_replicas'] = len(valid_replicas)
        result['error'] = 'molten ensemble incomplete: ' + '; '.join(
            str(row.get('error', 'invalid replica')) for row in failures
        )
        return result

    result = dict(valid_replicas[0])
    for field in MOLTEN_ENSEMBLE_FIELDS:
        result[field] = float(
            np.mean([float(row[field]) for row in valid_replicas])
        )
    result['structure_seed'] = None
    result['ensemble_size'] = MOLTEN_ENSEMBLE_SIZE
    result['ensemble_valid_replicas'] = len(valid_replicas)
    result['ensemble_failed_replicas'] = len(failures)
    result['E_act'] = max(0.01, min(float(result['E_act_raw']), 5.0))
    result['E_act_censored'] = result['E_act'] in (0.01, 5.0)
    result['rate_constant'] = arrhenius_rate(
        1e13, result['E_act'], float(result['temperature_K'])
    )
    result['molten_replica_spread_json'] = json.dumps(
        {
            field: {
                'min': min(float(row[field]) for row in valid_replicas),
                'max': max(float(row[field]) for row in valid_replicas),
                'std': float(np.std([float(row[field]) for row in valid_replicas])),
            }
            for field in ('E_act_raw', 'dE_H', 'dE_CH3', 'dE_C')
        },
        sort_keys=True,
    )
    result['molten_replicas_json'] = json.dumps(
        [
            {
                'structure_seed': int(row['structure_seed']),
                'E_act_raw': float(row['E_act_raw']),
                'dE_H': float(row['dE_H']),
                'dE_C': float(row['dE_C']),
            }
            for row in valid_replicas
        ],
        sort_keys=True,
    )
    result['needs_dft_validation'] = True
    return result


def _compute_segregation_energy(slab, active_idx, calc, e_clean):
    """Compute segregation energy: tendency for dopant to stay at surface."""
    try:
        slab_swap = slab.copy()
        z_coords = slab_swap.positions[:, 2]
        bulk_idx = np.argsort(z_coords)[0]
        surf_idx = active_idx[0] if active_idx[0] < len(slab_swap) else 0

        if bulk_idx < len(slab_swap) and surf_idx < len(slab_swap):
            sym_s = slab_swap[surf_idx].symbol
            sym_b = slab_swap[bulk_idx].symbol
            if sym_s != sym_b:
                slab_swap[surf_idx].symbol = sym_b
                slab_swap[bulk_idx].symbol = sym_s
                slab_swap.calc = calc
                e_swap = slab_swap.get_potential_energy()
                return e_clean - e_swap
        return 0.0
    except Exception:
        return 0.0


def _compute_binding_energy(structure, calc, e_clean, genome):
    """Compute metal-cavity binding energy for SAC/MOF."""
    try:
        cavity_struct = structure.copy()
        del cavity_struct[0]
        cavity_struct.calc = calc
        e_empty = cavity_struct.get_potential_energy()

        ref_energies = {
            'Fe': -5.0,
            'Co': -4.8,
            'Ni': -5.5,
            'Mn': -4.2,
            'Cu': -3.8,
            'Zn': -2.5,
            'Mo': -8.0,
            'W': -9.5,
            'V': -7.0,
            'Cr': -6.0,
            'Ru': -6.5,
            'Rh': -5.5,
            'Ti': -6.5,
            'Zr': -7.5,
            'Sn': -3.0,
            'In': -2.5,
            'None': 0.0,
        }
        metal = genome[1]
        e_metal = ref_energies.get(metal, -4.0)
        return (e_clean - e_empty - e_metal) / 10.0
    except Exception:
        return 0.0


# ═══════════════════════════════════════════════════════════════════════════════
# WORKER PROCESS & MAIN SCREENING ENGINE
# ═══════════════════════════════════════════════════════════════════════════════


def eval_worker(
    worker_id: int,
    gpu_id: int,
    gpu_uuid: str,
    task_queue: mp.Queue,
    result_queue: mp.Queue,
    stop_event,
    candidate_threads: int = 1,
    batched: bool = False,
):
    """GPU worker; optionally share one dynamically batched model across threads.

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
        references = compute_reference_energies(calculator)
        return lambda genome, thread_calculator: evaluate_candidate(
            genome, thread_calculator, references
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
        'surface worker',
    )


def run_screening(
    genomes: List[tuple],
    db_filename: str = "surface_screening.csv",
    workers_per_gpu: int = 2,
    engine: str = 'batched',
) -> 'pd.DataFrame':
    """Run parallel Meta eSen-SM screening on a list of catalyst genomes.

    Args:
        genomes: List of catalyst genome tuples
        db_filename: Output CSV filename
        workers_per_gpu: Number of parallel workers per GPU
        engine: Batched production engine or explicitly requested legacy engine.

    Returns:
        DataFrame with all screening results
    """
    from pipeline.screening.gpu_executor import ScreeningRunSpec, run_gpu_screening

    df = run_gpu_screening(
        genomes,
        db_filename,
        workers_per_gpu,
        engine,
        eval_worker,
        logger,
        ScreeningRunSpec(
            banner='META ESEN-SM SURFACE CATALYST SCREENING',
            application='turquoise_hydrogen',
            manifest_path=SCREENING_DIR / 'surface_worker_health.json',
            output_subdir='screening',
            start_message='Screening {count} catalyst candidates...',
            completion_label='Screening',
            protocol_id=SCREENING_PROTOCOL_ID,
            cache_path=SCREENING_DIR / 'esen_result_cache.sqlite',
        ),
    )

    valid_df = df[df['valid'] == True]
    if len(valid_df) > 0:
        logger.info(f"  Valid candidates: {len(valid_df)}/{len(df)}")
        logger.info(f"  Best E_act: {valid_df['E_act'].min():.4f} eV")
        logger.info(f"  Best coking index: {valid_df['coking_index'].max():.4f} eV")
        logger.info(
            f"  dE_H range: [{valid_df['dE_H'].min():.3f}, {valid_df['dE_H'].max():.3f}] eV"
        )

    return df


if __name__ == '__main__':
    from pipeline.search.indexed_space import deterministic_tree_probes

    pop = deterministic_tree_probes(20)
    df = run_screening(pop, db_filename="test_screening.csv", workers_per_gpu=2)
    print(f"\nScreening complete: {len(df)} candidates evaluated")
