#!/usr/bin/env python3
# Reactor mechanism generation.
"""
Cantera-Compatible Mechanism Generator.

Generates YAML mechanism files for methane pyrolysis surface kinetics
using DFT/MACE-derived activation barriers.

Output follows the Cantera 3.x YAML format: top-level sections for
units, phases, species, and reactions.
"""

import numpy as np
import json
from pathlib import Path
from dataclasses import asdict, dataclass, field
from typing import Mapping, Optional

from pipeline.common.utils import (
    eV_to_J, MECHANISMS_DIR, setup_logger,
)

logger = setup_logger('reactor_mechanisms', 'reactor/mechanism_generation.log')

NA = 6.02214076e23  # Avogadro's number
EV_TO_J_MOL = eV_to_J * NA


@dataclass(frozen=True)
class CandidateKinetics:
    """Candidate-specific inputs and honest provenance for one mechanism.

    Screening adsorption energies are used to distinguish the thermochemistry
    of adsorbed H, CH3, and C. They are not silently reinterpreted as activation
    barriers. Missing elementary barriers retain declared template values until
    candidate-specific NEB or measured kinetics replaces them.
    """

    methane_activation_eV: float
    h_adsorption_eV: Optional[float] = None
    ch3_adsorption_eV: Optional[float] = None
    c_adsorption_eV: Optional[float] = None
    ch3_dehydrogenation_eV: Optional[float] = None
    ch2_dehydrogenation_eV: Optional[float] = None
    ch_dehydrogenation_eV: Optional[float] = None
    h2_desorption_eV: Optional[float] = None
    carbon_transfer_eV: Optional[float] = None
    site_density_mol_cm2: float = 2.5e-9
    screening_protocol: str = 'unknown'
    candidate_id: str = 'unknown'
    sources: Mapping[str, str] = field(default_factory=dict)

    @classmethod
    def from_screening_row(cls, row, candidate_id: str = 'unknown'):
        """Build kinetics from a pandas Series or ordinary mapping."""
        def finite(name):
            value = row.get(name)
            try:
                value = float(value)
            except (TypeError, ValueError):
                return None
            return value if np.isfinite(value) else None

        barrier = finite('E_act')
        if barrier is None or barrier <= 0:
            raise ValueError('a finite positive E_act is required')
        protocol = str(row.get('screening_protocol', 'unknown'))
        sources = {'methane_activation_eV': f'screening:{protocol}'}
        mapping = {
            'h_adsorption_eV': 'dE_H',
            'ch3_adsorption_eV': 'dE_CH3',
            'c_adsorption_eV': 'dE_C',
        }
        values = {}
        for target, source in mapping.items():
            values[target] = finite(source)
            if values[target] is not None:
                sources[target] = f'screening:{protocol}:{source}'
        return cls(methane_activation_eV=barrier, candidate_id=candidate_id,
                   screening_protocol=protocol, sources=sources, **values)

    def resolved(self) -> dict:
        """Return numerical values plus whether each was observed or templated."""
        defaults = {
            'ch3_dehydrogenation_eV': self.methane_activation_eV + 0.10,
            'ch2_dehydrogenation_eV': self.methane_activation_eV + 0.15,
            'ch_dehydrogenation_eV': self.methane_activation_eV + 0.05,
            'h2_desorption_eV': 0.8,
            'carbon_transfer_eV': 1.5,
        }
        values = asdict(self)
        provenance = dict(self.sources)
        for name, default in defaults.items():
            if values[name] is None:
                values[name] = default
                provenance[name] = 'template_default'
            else:
                provenance.setdefault(name, 'candidate_specific')
        values['provenance'] = provenance
        values['quantitative_status'] = (
            'candidate_specific' if not any(
                provenance.get(name) == 'template_default' for name in defaults)
            else 'screening_template_incomplete')
        return values


def write_gas_only_mechanism() -> Path:
    """Write a gas-phase-only CH₄ decomposition mechanism."""
    MECHANISMS_DIR.mkdir(parents=True, exist_ok=True)

    yaml_content = """\
units: {length: cm, time: s, quantity: mol, activation-energy: J/mol}

phases:
- name: gas
  thermo: ideal-gas
  elements: [C, H, Ar]
  species: [CH4, H2, C2H2, C2H4, C2H6, Ar]
  kinetics: gas
  state: {T: 1000.0, P: 1 atm}

species:
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

reactions:
- equation: 2 CH4 => C2H6 + H2
  rate-constant: {A: 2.3e+13, b: 0.0, Ea: 356000.0}
- equation: C2H6 => C2H4 + H2
  rate-constant: {A: 4.65e+13, b: 0.0, Ea: 273000.0}
- equation: C2H4 => C2H2 + H2
  rate-constant: {A: 1.0e+14, b: 0.0, Ea: 331000.0}
"""

    filepath = MECHANISMS_DIR / "gri30_ch4_subset.yaml"
    with open(filepath, 'w') as f:
        f.write(yaml_content)
    logger.info(f"Wrote gas-phase mechanism: {filepath}")
    return filepath


def write_full_mechanism(catalyst_name: str, E_act_CH4: float = None,
                          E_act_H_desorb: float = 0.8,
                          E_act_C_diffuse: float = 1.5,
                          site_density: float = 2.5e-9,
                          T_ref: float = 1000.0,
                          kinetics: CandidateKinetics = None) -> Path:
    """
    Write a complete Cantera mechanism file (gas + surface) to disk.
    """
    MECHANISMS_DIR.mkdir(parents=True, exist_ok=True)

    if kinetics is None:
        if E_act_CH4 is None:
            raise ValueError('E_act_CH4 or kinetics is required')
        kinetics = CandidateKinetics(
            methane_activation_eV=float(E_act_CH4),
            h2_desorption_eV=float(E_act_H_desorb),
            carbon_transfer_eV=float(E_act_C_diffuse),
            site_density_mol_cm2=float(site_density),
            sources={'methane_activation_eV': 'legacy_argument',
                     'h2_desorption_eV': 'legacy_argument',
                     'carbon_transfer_eV': 'legacy_argument'})
    values = kinetics.resolved()

    # Convert eV → J/mol. Adsorption energies alter surface enthalpies but are
    # not used as activation barriers.
    Ea_CH4 = values['methane_activation_eV'] * EV_TO_J_MOL
    Ea_CH3 = values['ch3_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_CH2 = values['ch2_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_CH = values['ch_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_H2 = values['h2_desorption_eV'] * EV_TO_J_MOL
    Ea_C = values['carbon_transfer_eV'] * EV_TO_J_MOL
    h0_h = (values['h_adsorption_eV'] * EV_TO_J_MOL
            if values['h_adsorption_eV'] is not None else -25000.0)
    h0_ch3 = (values['ch3_adsorption_eV'] * EV_TO_J_MOL
              if values['ch3_adsorption_eV'] is not None else -20000.0)
    h0_c = (values['c_adsorption_eV'] * EV_TO_J_MOL
            if values['c_adsorption_eV'] is not None else -40000.0)
    site_density = values['site_density_mol_cm2']

    yaml_content = f"""\
units: {{length: cm, time: s, quantity: mol, activation-energy: J/mol}}

phases:
- name: gas
  thermo: ideal-gas
  elements: [C, H, Ar]
  species: [CH4, H2, C2H2, C2H4, C2H6, Ar, C_graphite]
  kinetics: gas
  state: {{T: {T_ref:.1f}, P: 1 atm}}

- name: {catalyst_name}_surface
  thermo: ideal-surface
  elements: [C, H]
  species: [site, CH3_s, CH2_s, CH_s, H_s, C_s]
  kinetics: surface
  reactions: [{catalyst_name}_surface-reactions]
  site-density: {site_density:.3e} mol/cm^2
  adjacent-phases: [gas]

species:
- name: CH4
  composition: {{C: 1, H: 4}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [5.14987613, -0.0136709788, 4.91800599e-05, -4.84743026e-08, 1.66693956e-11,
      -10246.6, -4.64130376]
    - [0.074851495, 0.0133909467, -5.73285809e-06, 1.22292535e-09, -1.01815230e-13,
      -9468.34459, 18.437318]
- name: H2
  composition: {{H: 2}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [2.34433112, 7.98052075e-03, -1.9478151e-05, 2.01572094e-08, -7.37611761e-12,
      -917.935173, 0.683010238]
    - [2.93286575, 8.26608026e-04, -1.46402364e-07, 1.54100414e-11, -6.888048e-16,
      -813.065581, -1.02432865]
- name: C2H2
  composition: {{C: 2, H: 2}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [0.808681094, 0.0233615629, -3.55171815e-05, 2.80152437e-08, -8.50072974e-12,
      26428.9807, 13.9397051]
    - [4.14756964, 5.96166664e-03, -2.37294852e-06, 4.67412171e-10, -3.61235213e-14,
      25935.9992, -1.23028121]
- name: C2H4
  composition: {{C: 2, H: 4}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [3.95920148, -7.57052247e-03, 5.70990292e-05, -6.91588753e-08, 2.69884373e-11,
      5089.77593, 4.09733096]
    - [3.99182724, 0.0104833908, -3.71721342e-06, 5.94628366e-10, -3.53630386e-14,
      4268.65851, -0.269081762]
- name: C2H6
  composition: {{C: 2, H: 6}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [4.29142492, -5.50154270e-03, 5.99438288e-05, -7.08466285e-08, 2.68685771e-11,
      -11522.2055, 2.66682316]
    - [4.04666411, 0.0153538802, -5.47039485e-06, 8.77826544e-10, -5.23167531e-14,
      -12447.3273, -0.968698313]
- name: Ar
  composition: {{Ar: 1}}
  thermo:
    model: NASA7
    temperature-ranges: [200.0, 1000.0, 6000.0]
    data:
    - [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]
    - [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, 4.366]
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
    h0: {h0_ch3:.8g} J/mol
    s0: 50.0 J/mol/K
  sites: 1
- name: CH2_s
  composition: {{C: 1, H: 2}}
  thermo:
    model: constant-cp
    h0: -15000.0 J/mol
    s0: 40.0 J/mol/K
  sites: 1
- name: CH_s
  composition: {{C: 1, H: 1}}
  thermo:
    model: constant-cp
    h0: -10000.0 J/mol
    s0: 30.0 J/mol/K
  sites: 1
- name: H_s
  composition: {{H: 1}}
  thermo:
    model: constant-cp
    h0: {h0_h:.8g} J/mol
    s0: 20.0 J/mol/K
  sites: 1
- name: C_s
  composition: {{C: 1}}
  thermo:
    model: constant-cp
    h0: {h0_c:.8g} J/mol
    s0: 10.0 J/mol/K
  sites: 1
- name: C_graphite
  composition: {{C: 1}}
  thermo:
    model: constant-cp
    h0: 0.0 J/mol
    s0: 5.74 J/mol/K
  note: Solid carbon product (treated as gas-phase trace for balance)

reactions:
- equation: 2 CH4 => C2H6 + H2
  rate-constant: {{A: 2.3e+13, b: 0.0, Ea: 356000.0}}
- equation: C2H6 => C2H4 + H2
  rate-constant: {{A: 4.65e+13, b: 0.0, Ea: 273000.0}}
- equation: C2H4 => C2H2 + H2
  rate-constant: {{A: 1.0e+14, b: 0.0, Ea: 331000.0}}

{catalyst_name}_surface-reactions:
- equation: CH4 + 2 site <=> CH3_s + H_s
  sticking-coefficient: {{A: 0.01, b: 0.0, Ea: {Ea_CH4:.1f}}}
- equation: CH3_s + site <=> CH2_s + H_s
  rate-constant: {{A: 1.0e+13, b: 0.0, Ea: {Ea_CH3:.1f}}}
- equation: CH2_s + site <=> CH_s + H_s
  rate-constant: {{A: 1.0e+13, b: 0.0, Ea: {Ea_CH2:.1f}}}
- equation: CH_s + site <=> C_s + H_s
  rate-constant: {{A: 1.0e+13, b: 0.0, Ea: {Ea_CH:.1f}}}
- equation: 2 H_s <=> H2 + 2 site
  rate-constant: {{A: 5.0e+13, b: 0.0, Ea: {Ea_H2:.1f}}}
- equation: C_s => C_graphite + site
  rate-constant: {{A: 1.0e+10, b: 0.0, Ea: {Ea_C:.1f}}}
"""

    filepath = MECHANISMS_DIR / f"mechanism_{catalyst_name}.yaml"
    with open(filepath, 'w') as f:
        f.write(yaml_content)

    metadata_path = filepath.with_suffix('.kinetics.json')
    metadata_path.write_text(json.dumps({
        'schema_version': 1,
        'catalyst_name': catalyst_name,
        'mechanism_file': str(filepath),
        'inputs': values,
        'carbon_phase_model': 'legacy_gas_tracer',
    }, indent=2, sort_keys=True) + '\n')
    logger.info(
        f"Wrote mechanism: {filepath} "
        f"(E_act={values['methane_activation_eV']:.3f} eV, "
        f"status={values['quantitative_status']})")
    return filepath


# Aliases
write_gri30_subset = write_gas_only_mechanism


if __name__ == '__main__':
    write_gas_only_mechanism()
    write_full_mechanism("NiBi_10pct", E_act_CH4=0.85, T_ref=1000)
    write_full_mechanism("FeNi_graphene", E_act_CH4=0.65, T_ref=900)
    write_full_mechanism("CuSn_20pct", E_act_CH4=1.10, T_ref=1100)
    print("Mechanism files written to:", MECHANISMS_DIR)
