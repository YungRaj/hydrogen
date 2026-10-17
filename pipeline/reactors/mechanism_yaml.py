"""Cantera YAML mechanism templates and kinetics sidecar formatting."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from pipeline.utils import repo_relative


def _jsonable(obj: Any) -> Any:
    if isinstance(obj, tuple):
        return [_jsonable(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_jsonable(x) for x in obj]
    return obj


# NASA-7 graphite (C(gr)) from Cantera's graphite.yaml / NASA thermo.
_GRAPHITE_SPECIES_YAML = """\
- name: C(gr)
  composition: {C: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 5000.0]
    data:
    - [-0.310872072, 4.40353686e-03, 1.90394118e-06, -6.38546966e-09, 2.98964248e-12,
      -108.650794, 1.11382953]
    - [1.45571829, 1.71702216e-03, -6.97562786e-07, 1.35277032e-10, -9.67590652e-15,
      -695.138814, -8.52583033]
  equation-of-state:
    model: constant-volume
    density: 2.16 g/cm^3
  note: Condensed graphite (fixed-stoichiometry); unit activity, not a gas species
"""

_GAS_SPECIES_YAML = """\
- name: CH4
  composition: {C: 1, H: 4}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11,
      -10246.6, -4.64130376]
    - [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.01815230e-13,
      -9468.34459, 18.437318]
- name: H2
  composition: {H: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12,
      -917.935173, 0.683010238]
    - [2.93286575, 8.26608026e-04, -1.46402364e-07, 1.54100414e-11, -6.888048e-16,
      -813.065581, -1.02432865]
- name: C2H2
  composition: {C: 2, H: 2}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [0.808681094, 0.0233615629, -3.55171815e-05, 2.80152437e-08, -8.50072974e-12,
      26428.9807, 13.9397051]
    - [4.14756964, 5.96166664e-03, -2.37294852e-06, 4.67412171e-10, -3.61235213e-14,
      25935.9992, -1.23028121]
- name: C2H4
  composition: {C: 2, H: 4}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [3.95920148, -7.57052247e-03, 5.70990292e-05, -6.91588753e-08, 2.69884373e-11,
      5089.77593, 4.09733096]
    - [3.99182724, 0.0104833908, -3.71721342e-06, 5.94628366e-10, -3.53630386e-14,
      4268.65851, -0.269081762]
- name: C2H6
  composition: {C: 2, H: 6}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [4.29142492, -5.50154270e-03, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11,
      -11522.2055, 2.66682316]
    - [4.04666411, 0.0153538802, -5.47039485e-06, 8.77826544e-10, -5.23167531e-14,
      -12447.3273, -0.968698313]
- name: Ar
  composition: {Ar: 1}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]
    - [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]
"""


def _gas_graphite_phases_yaml(T_ref: float) -> str:
    return f"""\
- name: gas
  thermo: ideal-gas
  elements: [C, H, Ar]
  species: [CH4, H2, C2H2, C2H4, C2H6, Ar]
  kinetics: gas
  state: {{T: {T_ref:.1f}, P: 1 atm}}

- name: graphite
  thermo: fixed-stoichiometry
  elements: [C]
  species: [C(gr)]
  state: {{T: {T_ref:.1f}, P: 1 atm}}
"""


def _surface_phase_yaml(
    catalyst_name: str, site_density: float, include_off_site: bool
) -> str:
    surf_species = (
        '[site, CH3_s, CH2_s, CH_s, H_s, C_s, C_encap_s]'
        if include_off_site
        else '[site, CH3_s, CH2_s, CH_s, H_s, C_s]'
    )
    adjacent = '[gas, graphite]' if include_off_site else '[gas]'
    return f"""\
- name: {catalyst_name}_surface
  thermo: ideal-surface
  elements: [C, H]
  species: {surf_species}
  kinetics: surface
  reactions: [{catalyst_name}_surface-reactions]
  site-density: {site_density:.3e} mol/cm^2
  adjacent-phases: {adjacent}
"""


def _surface_species_yaml(h0: dict, include_off_site: bool) -> str:
    species = f"""\
- name: site
  composition: {{}}
  thermo:
    model: constant-cp
    h0: 0.0 J/mol
    s0: 0.0 J/mol/K
  sites: 1
- name: CH3_s
  composition: {{C: 1, H: 3}}
  thermo:
    model: constant-cp
    h0: {h0['CH3_s']:.8g} J/mol
    s0: 50.0 J/mol/K
  sites: 1
- name: CH2_s
  composition: {{C: 1, H: 2}}
  thermo:
    model: constant-cp
    h0: {h0['CH2_s']:.8g} J/mol
    s0: 40.0 J/mol/K
  sites: 1
  note: template enthalpy, CH3_s + (C_s - CH3_s)/3
- name: CH_s
  composition: {{C: 1, H: 1}}
  thermo:
    model: constant-cp
    h0: {h0['CH_s']:.8g} J/mol
    s0: 30.0 J/mol/K
  sites: 1
  note: template enthalpy, CH3_s + 2 (C_s - CH3_s)/3
- name: H_s
  composition: {{H: 1}}
  thermo:
    model: constant-cp
    h0: {h0['H_s']:.8g} J/mol
    s0: 20.0 J/mol/K
  sites: 1
- name: C_s
  composition: {{C: 1}}
  thermo:
    model: constant-cp
    h0: {h0['C_s']:.8g} J/mol
    s0: 10.0 J/mol/K
  sites: 1
  note: Surface carbon (Cα); occupies catalytic sites until removed by policy or off-site Cγ
"""
    if include_off_site:
        species += f"""\
- name: C_encap_s
  composition: {{C: 1}}
  thermo:
    model: constant-cp
    h0: {h0['C_s']:.8g} J/mol
    s0: 10.0 J/mol/K
  sites: 1
  note: Encapsulating carbon (Cδ); occupies the site; no off-site return
"""
    return species


def _carbon_removal_reactions_yaml(
    rates: dict, include_off_site: bool, regeneration_rate_1_s: Optional[float]
) -> str:
    """Off-site Cγ/Cδ channels, or the optional diagnostic regeneration step."""
    if include_off_site:
        return f"""\
- equation: C_s => C(gr) + site
  rate-constant: {{A: {rates['A_Cgamma']:.6g}, b: 0.0, Ea: {rates['Ea_Cgamma']:.1f}}}
  note: Cα → Cγ transport-to-edge lump (Baker / Abild-Pedersen); not nucleation
- equation: C_s => C_encap_s
  rate-constant: {{A: {rates['A_Cdelta']:.6g}, b: 0.0, Ea: {rates['Ea_Cdelta']:.1f}}}
  coverage-dependencies:
    C_s: {{a: 0.0, m: 1.0, E: 0.0}}
  note: Cα → Cδ encapsulating, rate ∝ θ_C² (Cantera k·10^(a θ)·θ^m with m=1 times [C_s]); A = A_γ/θ*, θ* = {rates['theta_star']:g}; site stays blocked
"""
    if regeneration_rate_1_s is not None:
        return f"""\
- equation: C_s => C(gr) + site
  rate-constant: {{A: {float(regeneration_rate_1_s):.6g}, b: 0.0, Ea: 0.0}}
  note: Diagnostic mechanically driven carbon regeneration; non-pruning native model
"""
    return ""


def _surface_reactions_yaml(
    catalyst_name: str, rates: dict, carbon_removal: str
) -> str:
    return f"""\
{catalyst_name}_surface-reactions:
- equation: CH4 + 2 site <=> CH3_s + H_s
  sticking-coefficient: {{A: {rates['s0_ch4']:.6g}, b: 0.0, Ea: {rates['Ea_CH4']:.1f}}}
  note: s0 sets the carbon arrival rate; provenance in the .kinetics.json sidecar
- equation: CH3_s + site <=> CH2_s + H_s
  rate-constant: {{A: {rates['A_bimol']:.6g}, b: 0.0, Ea: {rates['Ea_CH3']:.1f}}}
  note: A = 1e13/s / Γ in cm^2/mol/s (bimolecular surface TST)
- equation: CH2_s + site <=> CH_s + H_s
  rate-constant: {{A: {rates['A_bimol']:.6g}, b: 0.0, Ea: {rates['Ea_CH2']:.1f}}}
- equation: CH_s + site <=> C_s + H_s
  rate-constant: {{A: {rates['A_bimol']:.6g}, b: 0.0, Ea: {rates['Ea_CH']:.1f}}}
- equation: 2 H_s <=> H2 + 2 site
  rate-constant: {{A: {rates['A_h2_des']:.6g}, b: 0.0, Ea: {rates['Ea_H2']:.1f}}}
{carbon_removal}"""


def _write_kinetics_sidecar(
    filepath: Path,
    catalyst_name: str,
    values: dict,
    rates: dict,
    h0: dict,
    include_surface_sites: bool,
    include_off_site: bool,
    *,
    surface_thermo_reference: dict,
    ch4_sticking_provenance: str,
    c_gamma_prefactor_provenance: str,
    c_gamma_provenance: str,
    c_delta_form: str,
    theta_star_provenance: str,
    c_delta_provenance: str,
) -> None:
    """Write the provenance sidecar recording every kinetic input and prefactor."""
    provenance = values['provenance']
    record = {
        'schema_version': 1,
        'catalyst_name': catalyst_name,
        'mechanism_file': repo_relative(filepath),
        'inputs': values,
        'carbon_phase_model': (
            'condensed_graphite_plus_surface_C_s_off_site_Cgamma_Cdelta'
            if include_off_site
            else 'condensed_graphite_plus_surface_C_s'
        ),
        'surface_thermo_reference': surface_thermo_reference,
        'surface_prefactors': (
            {
                'bimolecular_cm2_mol_s': rates['A_bimol'],
                'h2_desorption_cm2_mol_s': rates['A_h2_des'],
                'basis': 'k_TST(1/s) / site_density(mol/cm^2)',
                'ch4_sticking_coefficient': rates['s0_ch4'],
                'ch4_sticking_source': provenance.get(
                    'ch4_sticking_coefficient', ch4_sticking_provenance
                ),
            }
            if include_surface_sites
            else {}
        ),
        'surface_enthalpies_J_mol': dict(h0) if include_surface_sites else {},
        'off_site_carbon': include_off_site,
        'coking_index_mapped_to_off_site': False,
        'off_site_channels': (
            {
                'C_gamma': {
                    'equation': 'C_s => C(gr) + site',
                    'barrier_eV': float(values['carbon_transfer_eV']),
                    'preexponential_1_s': rates['A_Cgamma'],
                    'preexponential_source': provenance.get(
                        'carbon_transfer_prefactor_1_s',
                        c_gamma_prefactor_provenance,
                    ),
                    'particle_nm': values.get('carbon_transfer_particle_nm'),
                    'form': 'first_order_theta_C',
                    'kind': 'transport_to_edge',
                    'source': c_gamma_provenance,
                },
                'C_delta': {
                    'equation': 'C_s => C_encap_s',
                    'barrier_eV': float(values['carbon_encapsulation_eV']),
                    'preexponential_1_s': rates['A_Cdelta'],
                    'form': c_delta_form,
                    'crossover_coverage_theta_star': rates['theta_star'],
                    'crossover_coverage_source': provenance.get(
                        'encapsulation_crossover_coverage', theta_star_provenance
                    ),
                    'effective_crossover_note': (
                        'theta_x(T) = theta_star * exp((Ea_delta - Ea_gamma) / kT); '
                        'C_gamma wins below theta_x, C_delta above'
                    ),
                    'kind': 'encapsulating',
                    'source': c_delta_provenance,
                },
            }
            if include_off_site
            else {}
        ),
    }
    filepath.with_suffix('.kinetics.json').write_text(
        json.dumps(_jsonable(record), indent=2, sort_keys=True) + '\n',
        encoding='utf-8',
    )
