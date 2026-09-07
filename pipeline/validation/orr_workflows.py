"""Multi-site ORR validation with explicit CHE environmental corrections."""

from collections import defaultdict
from dataclasses import dataclass, asdict
from typing import Sequence
import numpy as np
from ase import Atoms

from pipeline.common.utils import orr_overpotential


@dataclass(frozen=True)
class ORRCorrections:
    solvation_OH_eV: float = -0.30
    solvation_O_eV: float = 0.00
    solvation_OOH_eV: float = -0.35
    electrode_potential_V: float = 0.0
    pH: float = 0.0
    temperature_K: float = 298.15
    source_id: str = ''


def build_orr_validation_plan(atoms: Atoms,
                              coverages: Sequence[float] = (0.25, 0.5, 1.0),
                              top_tolerance_A: float = 0.75) -> list[dict]:
    """Create deterministic site/coverage/adsorbate tasks for explicit DFT.

    The coverage value is a declared model condition. Structure builders must
    realize it in a suitable supercell; this function never equates a label
    with a physically constructed coverage.
    """
    values = tuple(float(value) for value in coverages)
    if not values or any(not 0 < value <= 1 for value in values):
        raise ValueError('ORR coverages must lie in (0, 1]')
    sites = enumerate_surface_sites(atoms, top_tolerance_A=top_tolerance_A)
    plan = []
    for site_number, site in enumerate(sites):
        site_id = f"{site['kind']}_{site_number:04d}"
        for coverage in sorted(set(values)):
            for adsorbate in ('OH', 'O', 'OOH'):
                plan.append({
                    'task_id': f'{site_id}:cov={coverage:g}:{adsorbate}',
                    'site_id': site_id,
                    'site': site,
                    'coverage_ML': coverage,
                    'adsorbate': adsorbate,
                    'requires_explicit_structure': True,
                })
    return plan


def enumerate_surface_sites(atoms: Atoms, top_tolerance_A: float = 0.75) -> list[dict]:
    """Enumerate symmetry-unreduced atop, bridge and hollow trial sites."""
    z = atoms.positions[:, 2]
    top = np.where(z >= z.max() - top_tolerance_A)[0].tolist()
    sites = [{'kind': 'atop', 'atom_indices': [i], 'position': atoms.positions[i].tolist()}
             for i in top]
    for a, i in enumerate(top):
        for j in top[a+1:]:
            distance = atoms.get_distance(i, j, mic=True)
            if distance < 3.5:
                sites.append({'kind': 'bridge', 'atom_indices': [i, j],
                              'position': ((atoms.positions[i] + atoms.positions[j]) / 2).tolist()})
    if len(top) >= 3:
        for a in range(len(top)-2):
            tri = top[a:a+3]
            sites.append({'kind': 'hollow', 'atom_indices': tri,
                          'position': atoms.positions[tri].mean(axis=0).tolist()})
    return sites


def apply_orr_corrections(dg_oh: float, dg_o: float, dg_ooh: float,
                          corrections: ORRCorrections) -> dict:
    """Apply solvation to states and potential/pH to each one-electron step."""
    if not corrections.source_id:
        raise ValueError('ORR corrections require a traceable source_id')
    kbt_ln10_eV = 8.617333262e-5 * corrections.temperature_K * np.log(10.0)
    pcet_shift = corrections.electrode_potential_V + kbt_ln10_eV * corrections.pH
    values = {
        'dG_OH_eV': dg_oh + corrections.solvation_OH_eV,
        'dG_O_eV': dg_o + corrections.solvation_O_eV,
        'dG_OOH_eV': dg_ooh + corrections.solvation_OOH_eV,
    }
    base_steps = {
        'step_1_OOH': values['dG_OOH_eV'] - 4.92,
        'step_2_O': values['dG_O_eV'] - values['dG_OOH_eV'],
        'step_3_OH': values['dG_OH_eV'] - values['dG_O_eV'],
        'step_4_H2O': -values['dG_OH_eV'],
    }
    operating_steps = {
        name: value + pcet_shift for name, value in base_steps.items()}
    return {
        **values,
        'orr_steps_U0_pH0_eV': base_steps,
        'orr_steps_at_condition_eV': operating_steps,
        'pcet_shift_per_step_eV': float(pcet_shift),
        'corrections': asdict(corrections),
        'evidence_level': 'corrected_DFT',
    }


def select_lowest_site(site_results: list[dict], key: str) -> dict:
    valid = [x for x in site_results if x.get('converged') is True and
             np.isfinite(float(x.get(key, np.nan)))]
    if not valid:
        raise RuntimeError(f'no converged adsorption site for {key}')
    return min(valid, key=lambda x: (float(x[key]), str(x.get('site_id', ''))))


def evaluate_orr_ensemble(site_results: list[dict],
                          correction_models: Sequence[ORRCorrections],
                          expected_cases: int | None = None) -> dict:
    """Evaluate complete site/coverage pathways and expose model uncertainty.

    Each input row is one adsorbate calculation with ``site_id``,
    ``coverage_ML``, ``adsorbate``, ``dG_eV`` and ``converged``. A pathway is
    admitted only when OH, O and OOH all converged for the same site/coverage.
    Every correction model requires provenance. Missing cases keep the result
    fail-closed rather than making the best surviving site look definitive.
    """
    if not correction_models:
        raise ValueError('at least one correction model is required')
    if any(not correction.source_id for correction in correction_models):
        raise ValueError('all ORR correction models require source_id')
    grouped = defaultdict(dict)
    invalid_rows = 0
    for row in site_results:
        try:
            key = (str(row['site_id']), float(row['coverage_ML']))
            adsorbate = str(row['adsorbate'])
            value = float(row['dG_eV'])
        except (KeyError, TypeError, ValueError):
            invalid_rows += 1
            continue
        if adsorbate not in ('OH', 'O', 'OOH') or not row.get('converged') \
                or not np.isfinite(value):
            invalid_rows += 1
            continue
        grouped[key][adsorbate] = value

    cases = []
    for (site_id, coverage), values in sorted(grouped.items()):
        if set(values) != {'OH', 'O', 'OOH'}:
            continue
        for correction in correction_models:
            corrected = apply_orr_corrections(
                values['OH'], values['O'], values['OOH'], correction)
            eta, limiting_step = orr_overpotential(
                corrected['dG_OH_eV'], corrected['dG_O_eV'],
                corrected['dG_OOH_eV'])
            cases.append({
                'site_id': site_id,
                'coverage_ML': coverage,
                'correction_source_id': correction.source_id,
                'orr_overpotential_V': float(eta),
                'limiting_step': limiting_step,
                **corrected,
            })
    anticipated = int(expected_cases) if expected_cases is not None else len(cases)
    if anticipated < 0:
        raise ValueError('expected_cases cannot be negative')
    complete = bool(cases) and len(cases) == anticipated and invalid_rows == 0
    if not cases:
        return {
            'complete': False, 'evidence_level': 'incomplete', 'cases': [],
            'case_count': 0, 'expected_cases': anticipated,
            'invalid_or_incomplete_rows': invalid_rows,
            'orr_overpotential_V': None,
        }
    eta = np.asarray([case['orr_overpotential_V'] for case in cases])
    best = min(cases, key=lambda case: (
        case['orr_overpotential_V'], case['site_id'],
        case['coverage_ML'], case['correction_source_id']))
    return {
        'complete': complete,
        'evidence_level': 'corrected_DFT_ensemble' if complete else 'incomplete',
        'cases': cases,
        'case_count': len(cases),
        'expected_cases': anticipated,
        'invalid_or_incomplete_rows': invalid_rows,
        'best_case': best,
        'orr_overpotential_V': best['orr_overpotential_V'] if complete else None,
        'uncertainty': {
            'definition': 'spread_across_declared_site_coverage_correction_cases',
            'mean_V': float(np.mean(eta)),
            'std_V': float(np.std(eta)),
            'min_V': float(np.min(eta)),
            'max_V': float(np.max(eta)),
            'range_V': float(np.ptp(eta)),
        },
    }
