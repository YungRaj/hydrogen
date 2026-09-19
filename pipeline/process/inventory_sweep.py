"""B1 solids-inventory sweep on the general staged-sweep harness.

Coarse (gridded) and targeted (ROI) results live under
`results/sweeps/inventory_b1/`. The coarse JSON is write-once.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

import pandas as pd

from pipeline.common.utils import (
    REACTOR_DIR, SCREENING_DIR, setup_logger,
)
from pipeline.process.reactor_mechanisms import (
    CandidateKinetics, write_full_mechanism)
from pipeline.process.reactor_models import (
    INVENTORY_METAL_DISPERSION, INVENTORY_METAL_LOADING, INVENTORY_PARTICLE_MM,
    SINGLE_REACTOR_MODE, ReactorConfig, run_reactor_sweep, simulate_reactor,
)
from pipeline.process.staged_sweep import (
    SweepSpec, cartesian_cells, import_existing_stage, load_stage,
    propose_roi, write_stage,
)

logger = setup_logger('inventory_sweep', 'reactor/inventory_sweep.log')

SWEEP_NAME = 'inventory_b1'
PROBE_INDEX = 9
CONTROL_INDEX = 40
SOLIDS_TYPES = ('PFR', 'Fluidized')
COARSE_TEMPERATURES = (773.15, 1300.0)
TARGETED_TEMPERATURES = (773.15, 1100.0, 1300.0)
LEGACY_COARSE = 'inventory_sweep_b1_2.json'
LEGACY_TARGETED = 'inventory_sweep_b1_2_roi.json'


def inventory_spec() -> SweepSpec:
    return SweepSpec(
        name=SWEEP_NAME,
        levers={
            'catalyst_particle_mm': list(INVENTORY_PARTICLE_MM),
            'metal_loading': list(INVENTORY_METAL_LOADING),
            'metal_dispersion': list(INVENTORY_METAL_DISPERSION),
        },
        score_key='CH4_conversion',
        constraint_key='ergun_ok',
        record_filter=lambda r: (
            r.get('role') == 'probe'
            and r.get('reactor_type') == 'PFR'
            and (r.get('T_K') or 0) >= 1200
        ),
        keep_fraction_of_max=0.30,
        n_targeted_per_lever=6,
        hard_bounds={
            'catalyst_particle_mm': (0.05, 2.0),
            'metal_loading': (0.0, 1.0),
            'metal_dispersion': (0.0, 1.0),
        },
    )


def _row_for_index(idx: int):
    path = SCREENING_DIR / 'ga_full_database.csv'
    df = pd.read_csv(path)
    if idx not in df.index:
        raise KeyError(f'screening row {idx} missing from {path}')
    row = df.loc[idx]
    if not bool(row.get('valid', True)):
        raise ValueError(f'screening row {idx} is not valid')
    return row


def _write_mech(idx: int, row) -> Path:
    """Write the probe/control mechanism from its screening row.

    Uses the same ``CandidateKinetics`` path as the headline YAML sweep and
    the orchestrator, so the inventory probe carries the row's adsorbate
    thermochemistry (dE_H, dE_CH3, dE_C) and B6 class gate. The template
    path (``E_act_CH4=``) has default enthalpies and gave a different
    catalyst under the same name.
    """
    kinetics = CandidateKinetics.from_screening_row(
        row, candidate_id=f'cat_{idx}')
    return write_full_mechanism(f'inv_cat_{idx}', kinetics=kinetics)


def _material_class(row) -> Optional[str]:
    value = row.get('material_class')
    return str(value) if isinstance(value, str) and value else None


def _record(role: str, cat: str, r: dict, cell: dict = None) -> dict:
    src = cell or r
    return {
        'role': role,
        'cat': cat,
        'reactor_type': r.get('reactor_type'),
        'status': r.get('status'),
        'reason': r.get('reason'),
        'material_class': r.get('material_class'),
        'T_K': r.get('T_K'),
        'CH4_conversion': r.get('CH4_conversion'),
        'single_pass_CH4_conversion': r.get(
            'single_pass_CH4_conversion', r.get('CH4_conversion')),
        'conversion_basis': r.get('conversion_basis'),
        'catalyst_particle_mm': src.get('catalyst_particle_mm', r.get('catalyst_particle_mm')),
        'metal_loading': src.get('metal_loading', r.get('metal_loading')),
        'metal_dispersion': src.get('metal_dispersion', r.get('metal_dispersion')),
        'geometric_sv_1_m': r.get('geometric_sv_1_m'),
        'active_sv_1_m': r.get('active_sv_1_m'),
        'WHSV_h-1': r.get('WHSV_h-1'),
        'ergun_delta_p_Pa': r.get('ergun_delta_p_Pa'),
        'ergun_delta_p_bar': r.get('ergun_delta_p_bar'),
        'ergun_ok': r.get('ergun_ok'),
        'catalyst_E_act_eV': r.get('catalyst_E_act_eV'),
        'catalyst_dE_H_eV': r.get('catalyst_dE_H_eV'),
    }


def _evaluate_cells(cells, temperatures, probe, control, probe_mech, control_mech):
    policy = {
        'co2_permitted': False,
        'fluidized_mode': 'circulating',
        'max_regen_cycles': 3,
        'regen_mechanism': 'mechanical',
    }
    records = []
    probe_class = _material_class(probe)
    control_class = _material_class(control)
    for cell in cells:
        # Both solids reactors share the thermocatalytic mode.
        results = run_reactor_sweep(
            f'inv_cat_{PROBE_INDEX}',
            str(probe_mech),
            temperatures=list(temperatures),
            reactor_types=list(SOLIDS_TYPES),
            catalyst_E_act_eV=float(probe.get('E_act', 0.8)),
            pathway_mode='thermocatalytic',
            material_class=probe_class,
            candidate_id=f'cat_{PROBE_INDEX}',
            catalyst_dE_H_eV=float(probe.get('dE_H', 0.0)),
            reactor_config_kwargs={**policy, **cell},
        )
        for r in results:
            records.append(_record('probe', f'cat_{PROBE_INDEX}', r, cell))

    if cells:
        leak_cells = (cells[0], cells[-1])
    else:
        leak_cells = ()
    for cell in leak_cells:
        # MMBCR on a solids probe is not_applicable under the class gate; the
        # record is kept so the leak check is visible, not silently dropped.
        r = simulate_reactor(ReactorConfig(
            T_inlet_K=1300.0,
            reactor_type='MMBCR',
            pathway_mode=SINGLE_REACTOR_MODE['MMBCR'],
            material_class=probe_class,
            mechanism_file=str(probe_mech),
            catalyst_name=f'inv_cat_{PROBE_INDEX}',
            candidate_id=f'cat_{PROBE_INDEX}',
            catalyst_E_act_eV=float(probe.get('E_act', 0.8)),
            catalyst_dE_H_eV=float(probe.get('dE_H', 0.0)),
            **policy,
            **cell,
        ))
        records.append(_record('mmbcr_leak', f'cat_{PROBE_INDEX}', r, cell))
        r = simulate_reactor(ReactorConfig(
            T_inlet_K=1300.0,
            reactor_type='PFR',
            pathway_mode=SINGLE_REACTOR_MODE['PFR'],
            material_class=control_class,
            mechanism_file=str(control_mech),
            catalyst_name=f'inv_cat_{CONTROL_INDEX}',
            candidate_id=f'cat_{CONTROL_INDEX}',
            catalyst_E_act_eV=float(control.get('E_act', 0.01)),
            catalyst_dE_H_eV=float(control.get('dE_H', 0.0)),
            **policy,
            **cell,
        ))
        records.append(_record('h_blocked_control', f'cat_{CONTROL_INDEX}', r, cell))
    return records


def _write_legacy_copy(filename: str, payload: dict) -> None:
    REACTOR_DIR.mkdir(parents=True, exist_ok=True)
    path = REACTOR_DIR / filename
    if path.exists():
        logger.info(f'Leaving existing legacy sweep file in place: {path}')
        return
    path.write_text(json.dumps(payload, indent=2), encoding='utf-8')


def _load_coarse_payload() -> dict:
    archived = load_stage(SWEEP_NAME, 'coarse')
    if archived is not None:
        return archived
    legacy = REACTOR_DIR / LEGACY_COARSE
    if legacy.exists():
        return json.loads(legacy.read_text(encoding='utf-8'))
    raise FileNotFoundError(
        f'no coarse sweep for {SWEEP_NAME}; run stage=coarse first')


def propose_inventory_roi(coarse_payload: dict = None) -> dict:
    payload = coarse_payload if coarse_payload is not None else _load_coarse_payload()
    return propose_roi(payload, inventory_spec())


def run_inventory_sweep(stage: str = 'coarse',
                        overwrite_coarse: bool = False) -> dict:
    spec = inventory_spec()
    if stage == 'propose':
        roi = propose_inventory_roi()
        logger.info(f'Proposed ROI: {roi["reason"]}')
        return {'stage': 'propose', 'roi': roi}

    if stage == 'archive':
        imported = {}
        coarse_src = REACTOR_DIR / LEGACY_COARSE
        targeted_src = REACTOR_DIR / LEGACY_TARGETED
        if coarse_src.exists():
            imported['coarse'] = str(import_existing_stage(
                SWEEP_NAME, 'coarse', coarse_src, overwrite=overwrite_coarse))
        if targeted_src.exists():
            imported['targeted'] = str(import_existing_stage(
                SWEEP_NAME, 'targeted', targeted_src, overwrite=True))
        return {'stage': 'archive', 'imported': imported}

    probe = _row_for_index(PROBE_INDEX)
    control = _row_for_index(CONTROL_INDEX)
    probe_mech = _write_mech(PROBE_INDEX, probe)
    control_mech = _write_mech(CONTROL_INDEX, control)
    roi = None
    if stage == 'coarse':
        cells = cartesian_cells(spec.levers)
        temperatures = COARSE_TEMPERATURES
        levels = {k: list(v) for k, v in spec.levers.items()}
        legacy_name = LEGACY_COARSE
    elif stage == 'targeted':
        coarse = _load_coarse_payload()
        roi = propose_roi(coarse, spec)
        cells = cartesian_cells(roi['levers'])
        temperatures = TARGETED_TEMPERATURES
        levels = roi['levers']
        legacy_name = LEGACY_TARGETED
    else:
        raise ValueError(f'unknown inventory sweep stage {stage!r}')

    records = _evaluate_cells(
        cells, temperatures, probe, control, probe_mech, control_mech)
    payload = {
        'stage': stage,
        'sweep_name': SWEEP_NAME,
        'area_model': 'a_geom * loading * dispersion (both <= 1)',
        'roi': roi,
        'levels': levels,
        'temperatures_K': list(temperatures),
        'n_grid_cells': len(cells),
        'probe': {
            'index': PROBE_INDEX,
            'E_act': float(probe.get('E_act', 0.8)),
            'dE_H': float(probe.get('dE_H', 0.0)),
            'genome': probe.get('genome'),
        },
        'control': {
            'index': CONTROL_INDEX,
            'E_act': float(control.get('E_act', 0.01)),
            'dE_H': float(control.get('dE_H', 0.0)),
            'genome': control.get('genome'),
        },
        'records': records,
    }
    write_stage(SWEEP_NAME, stage, payload, overwrite=overwrite_coarse and stage == 'coarse')
    _write_legacy_copy(legacy_name, payload)
    logger.info(f'{stage} sweep: {len(records)} rows, {len(cells)} cells')
    return payload


if __name__ == '__main__':
    import sys
    stage = sys.argv[1] if len(sys.argv) > 1 else 'coarse'
    run_inventory_sweep(stage=stage)
