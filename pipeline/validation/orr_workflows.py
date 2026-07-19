"""Multi-site ORR validation with explicit CHE environmental corrections."""

from dataclasses import dataclass, asdict
import numpy as np
from ase import Atoms


@dataclass(frozen=True)
class ORRCorrections:
    solvation_OH_eV: float = -0.30
    solvation_O_eV: float = 0.00
    solvation_OOH_eV: float = -0.35
    electrode_potential_V: float = 0.0
    pH: float = 0.0
    temperature_K: float = 298.15
    source_id: str = ''


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
    """Apply declared solvation, potential and pH terms; provenance is required."""
    if not corrections.source_id:
        raise ValueError('ORR corrections require a traceable source_id')
    kbt_ln10_eV = 8.617333262e-5 * corrections.temperature_K * np.log(10.0)
    proton_term = corrections.electrode_potential_V + kbt_ln10_eV * corrections.pH
    values = {
        'dG_OH_eV': dg_oh + corrections.solvation_OH_eV - proton_term,
        'dG_O_eV': dg_o + corrections.solvation_O_eV - 2 * proton_term,
        'dG_OOH_eV': dg_ooh + corrections.solvation_OOH_eV - proton_term,
    }
    return {**values, 'corrections': asdict(corrections), 'evidence_level': 'corrected_DFT'}


def select_lowest_site(site_results: list[dict], key: str) -> dict:
    valid = [x for x in site_results if x.get('converged') is True and
             np.isfinite(float(x.get(key, np.nan)))]
    if not valid:
        raise RuntimeError(f'no converged adsorption site for {key}')
    return min(valid, key=lambda x: (float(x[key]), str(x.get('site_id', ''))))
