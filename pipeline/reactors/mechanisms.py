#!/usr/bin/env python3
# Reactor mechanism generation.
"""
Cantera-Compatible Mechanism Generator.

Generates YAML mechanism files for methane pyrolysis surface kinetics
using DFT/MACE-derived activation barriers.

Solid carbon is a condensed fixed-stoichiometry graphite phase (C(gr)),
not a gas-phase tracer. Surface carbon remains as C_s (site-blocking)
unless the class gate admits off-site Cγ / Cδ (nanoparticle Ni/Fe/Co).
Output follows the Cantera 3.x YAML format.
"""

import ast
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np

from pipeline.utils import (
    eV_to_J, MECHANISMS_DIR, repo_relative, setup_logger,
)

logger = setup_logger('reactor_mechanisms', 'reactor/mechanism_generation.log')

NA = 6.02214076e23  # Avogadro's number
EV_TO_J_MOL = eV_to_J * NA


@dataclass(frozen=True)
class CandidateKinetics:
    """Candidate-specific inputs and honest provenance for one mechanism.

    Screening adsorption energies distinguish adsorbed H, CH3, and C
    thermochemistry. They are not silently reinterpreted as activation
    barriers. Missing elementary barriers keep declared template values
    until candidate-specific NEB or measured kinetics replaces them.
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
    carbon_encapsulation_eV: Optional[float] = None
    # B6: C_s coverage at which encapsulation (Cδ, ∝ θ_C²) overtakes
    # transport-to-edge (Cγ, ∝ θ_C). Declared, not measured; swept in B6-6.
    encapsulation_crossover_coverage: Optional[float] = None
    # B6: Cγ prefactor (1/s). 1e13 is a single-hop TST label; the real
    # channel is transport + precipitation (D0/L², particle-size dependent)
    # and is swept in B6-6. Cδ inherits A_γ / θ*.
    carbon_transfer_prefactor_1_s: Optional[float] = None
    # B6-6: A_γ may instead be derived as D0 / L² from a metal particle
    # diameter (nm) and a carbon bulk-diffusion prefactor (m²/s). Setting
    # both the particle size and an explicit prefactor is an error.
    carbon_transfer_particle_nm: Optional[float] = None
    carbon_diffusion_prefactor_m2_s: Optional[float] = None
    # B6-6: CH4 dissociative sticking prefactor s0 (dimensionless). Sets the
    # carbon arrival rate at the surface; template 0.01 until swept.
    ch4_sticking_coefficient: Optional[float] = None
    site_density_mol_cm2: float = 2.5e-9
    screening_protocol: str = 'unknown'
    candidate_id: str = 'unknown'
    material_class: Optional[str] = None
    genome: Optional[tuple] = None
    sources: Mapping[str, str] = field(default_factory=dict)

    @classmethod
    def from_screening_row(cls, row, candidate_id: str = 'unknown',
                           validation: Optional[Mapping] = None):
        """Build kinetics from screening and optional converged NEB evidence.

        Validation values replace templates only when the candidate identity
        matches and the complete campaign carries the required evidence level.

        Args:
            row: pandas Series or ordinary mapping from the screening database.
            candidate_id: Stable candidate identifier.
            validation: Candidate-specific validation evidence.

        Returns:
            A ``CandidateKinetics`` with per-field provenance in ``sources``.
        """
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
        genome = parse_catalyst_genome(row.get('genome'))
        material_class = row.get('material_class')
        if material_class is None or (
                isinstance(material_class, float) and not np.isfinite(material_class)):
            material_class = genome[0] if genome else None
        elif material_class is not None:
            material_class = str(material_class)
        if validation is not None:
            if str(validation.get('candidate_id')) != str(candidate_id):
                raise ValueError('kinetics validation candidate_id mismatch')
            if not validation.get('complete') or validation.get(
                    'evidence_level') != 'converged_dft_neb_frequency':
                raise ValueError('kinetics validation campaign is incomplete')
            resolved = validation.get('resolved_kinetics_eV', {})
            allowed = {
                'methane_activation_eV', 'ch3_dehydrogenation_eV',
                'ch2_dehydrogenation_eV', 'ch_dehydrogenation_eV',
                'h2_desorption_eV', 'carbon_transfer_eV',
                'carbon_encapsulation_eV',
            }
            unknown = set(resolved) - allowed
            if unknown:
                raise ValueError(f'unknown validated kinetics fields: {sorted(unknown)}')
            for name, raw_value in resolved.items():
                value = float(raw_value)
                if not np.isfinite(value) or value < 0:
                    raise ValueError(f'invalid validated barrier for {name}')
                if name == 'methane_activation_eV':
                    barrier = value
                else:
                    values[name] = value
                sources[name] = 'candidate_specific:converged_dft_neb_frequency'
        return cls(methane_activation_eV=barrier, candidate_id=candidate_id,
                   screening_protocol=protocol, sources=sources,
                   material_class=material_class, genome=genome, **values)

    def resolved(self) -> dict:
        """Return numerical values plus whether each was observed or templated.

        Returns:
            Dictionary of resolved barriers, per-field ``provenance``, and
            ``quantitative_status`` (``candidate_specific`` or
            ``screening_template_incomplete``).
        """
        defaults = {
            'ch3_dehydrogenation_eV': self.methane_activation_eV + 0.10,
            'ch2_dehydrogenation_eV': self.methane_activation_eV + 0.15,
            'ch_dehydrogenation_eV': self.methane_activation_eV + 0.05,
            'h2_desorption_eV': 0.8,
            'carbon_transfer_eV': 1.5,
            'carbon_encapsulation_eV': 1.53,
        }
        values = asdict(self)
        provenance = dict(self.sources)
        for name, default in defaults.items():
            if values[name] is None:
                values[name] = default
                provenance[name] = 'template_default'
            else:
                provenance.setdefault(name, 'candidate_specific')
        if values['encapsulation_crossover_coverage'] is None:
            values['encapsulation_crossover_coverage'] = (
                DEFAULT_ENCAPSULATION_CROSSOVER_COVERAGE)
            provenance['encapsulation_crossover_coverage'] = THETA_STAR_PROVENANCE
        else:
            provenance.setdefault(
                'encapsulation_crossover_coverage', THETA_STAR_PROVENANCE)
        particle_nm = values['carbon_transfer_particle_nm']
        if particle_nm is not None:
            if values['carbon_transfer_prefactor_1_s'] is not None:
                raise ValueError(
                    'set carbon_transfer_particle_nm or '
                    'carbon_transfer_prefactor_1_s, not both')
            if not particle_nm > 0:
                raise ValueError('carbon_transfer_particle_nm must be positive')
            d0 = values['carbon_diffusion_prefactor_m2_s']
            if d0 is None:
                d0 = CARBON_DIFFUSION_PREFACTOR_M2_S
                values['carbon_diffusion_prefactor_m2_s'] = d0
                provenance['carbon_diffusion_prefactor_m2_s'] = (
                    CARBON_DIFFUSION_PREFACTOR_PROVENANCE)
            elif not d0 > 0:
                raise ValueError('carbon_diffusion_prefactor_m2_s must be positive')
            else:
                provenance.setdefault(
                    'carbon_diffusion_prefactor_m2_s', 'declared_not_measured')
            values['carbon_transfer_prefactor_1_s'] = (
                carbon_transfer_prefactor_from_particle(particle_nm, d0))
            provenance['carbon_transfer_prefactor_1_s'] = (
                f'derived: D0/L^2 with D0={d0:.3g} m^2/s, L={particle_nm:g} nm')
        elif values['carbon_transfer_prefactor_1_s'] is None:
            values['carbon_transfer_prefactor_1_s'] = OFF_SITE_PREEXPONENTIAL_1_S
            provenance['carbon_transfer_prefactor_1_s'] = C_GAMMA_PREFACTOR_PROVENANCE
        else:
            provenance.setdefault(
                'carbon_transfer_prefactor_1_s', 'declared_not_measured')
        if values['ch4_sticking_coefficient'] is None:
            values['ch4_sticking_coefficient'] = DEFAULT_CH4_STICKING_COEFFICIENT
            provenance['ch4_sticking_coefficient'] = CH4_STICKING_PROVENANCE
        else:
            provenance.setdefault('ch4_sticking_coefficient', 'declared_not_measured')
        values['provenance'] = provenance
        # The Cδ (encapsulation) barrier is only written into the mechanism
        # when the B6 class gate admits off-site carbon. For every other
        # candidate it is not part of the kinetics and must not downgrade
        # the status of an otherwise complete converged-NEB set.
        required = set(defaults)
        if not off_site_carbon_allowed(self.genome, self.material_class):
            required.discard('carbon_encapsulation_eV')
        values['quantitative_status'] = (
            'candidate_specific' if not any(
                provenance.get(name) == 'template_default' for name in required)
            else 'screening_template_incomplete')
        return values


# Physical monolayer. ~10^19 atoms/m^2 = 2.5e-9 mol/cm^2.
# Do not raise this to force Damköhler (B1). Extra sites come from
# particle S/V, loading, and dispersion only.
MONOLAYER_SITE_DENSITY_MOL_CM2 = 2.5e-9

# Baker / Helveg cycle needs an extended metal particle (B6-3).
OFF_SITE_CARBON_METALS = frozenset({'Ni', 'Fe', 'Co'})
OFF_SITE_CARBON_DENIED_CLASSES = frozenset({
    'SAC', 'DAC', 'MetalFreeCarbon', 'MoltenMetal',
    'MOF', 'COF', 'Perovskite', 'MetalHydride', 'MXene',
    'MAXPhase', 'Spinel',
})
C_GAMMA_PROVENANCE = 'template_default: Abild-Pedersen 2006 / Baker 1972'
C_DELTA_PROVENANCE = (
    'template_default: Amin IEC Res 2011 50 12460; '
    '147-149 kJ/mol encapsulating-carbon')
# Cδ is second order in θ_C (mean-field island nucleation; Snoeck
# supersaturation picture). θ* is the coverage where k_δ θ_C² = k_γ θ_C at
# equal barriers; A_δ = A_γ / θ*. With Ea_δ − Ea_γ = 0.03 eV the effective
# crossover is θ*·exp(0.03 eV / kT) (~0.73 at 923 K for θ* = 0.5).
DEFAULT_ENCAPSULATION_CROSSOVER_COVERAGE = 0.5
THETA_STAR_PROVENANCE = 'declared_not_measured: swept in B6-6'
C_DELTA_FORM = 'coverage_dependent_theta_C_squared'
OFF_SITE_PREEXPONENTIAL_1_S = 1.0e13
C_GAMMA_PREFACTOR_PROVENANCE = (
    'template_default: single-hop TST 1e13/s; transport+precipitation lump '
    'is D0/L^2 and particle-size dependent; swept in B6-6')
# B6-6: A_γ = D0 / L² from a particle diameter. Carbon bulk diffusion in
# Ni (Lander, Kern & Beach 1952): D = 2.48 cm²/s · exp(−40.2 kcal/mol / RT),
# i.e. D0 = 2.48e-4 m²/s with E ≈ 1.74 eV (close to carbon_transfer_eV
# 1.5). This D0 puts a 10 nm particle at 2.5e12 1/s, within a decade of the
# TST label; A_γ ~ 1e9 corresponds to L ≈ 500 nm on this scale. The lump is
# order-of-magnitude; the sweep, not the derivation, carries the result.
CARBON_DIFFUSION_PREFACTOR_M2_S = 2.48e-4
CARBON_DIFFUSION_PREFACTOR_PROVENANCE = (
    'template_default: Lander J.Appl.Phys. 1952 23 1305; '
    'D0 = 2.48 cm^2/s (ln D = 0.909 - 20200/T)')
# B6-6: CH4 dissociative sticking prefactor. Template order of magnitude
# (Deutschmann-style methane_pox_on_pt uses 0.01 on Pt); Ni(111) molecular-
# beam values span 1e-4 .. 1e-2 at these T. Sets the carbon arrival rate.
DEFAULT_CH4_STICKING_COEFFICIENT = 0.01
CH4_STICKING_PROVENANCE = (
    'template_default: s0 = 0.01 (methane_pox_on_pt order of magnitude); '
    'swept in B6-6')


def carbon_transfer_prefactor_from_particle(particle_nm: float,
                                            d0_m2_s: float) -> float:
    """A_γ = D0 / L² (1/s) for a metal particle of diameter ``particle_nm``.

    Args:
        d0_m2_s: Input controlling d0 m2 s.

    Returns:
        Validated float output for this operation.
    """
    length_m = float(particle_nm) * 1e-9
    return float(d0_m2_s) / (length_m * length_m)

# Surface rate constants. Unimolecular surface steps (C_s => ...) take A in
# 1/s. Bimolecular steps (X_s + site, 2 H_s) are mass-action in surface
# concentrations (mol/cm^2), so A is in cm^2/mol/s and the TST prefactor is
# k_TST / Γ = 1e13 / 2.5e-9 = 4e21 cm^2/mol/s (Deutschmann convention;
# Cantera's methane_pox_on_pt uses 3.7e21). Writing 1e13 cm^2/mol/s gives an
# effective 2.5e4 1/s and freezes the dehydrogenation ladder.
SURFACE_TST_PREFACTOR_1_S = 1.0e13
H2_DESORPTION_PREFACTOR_1_S = 5.0e13


def bimolecular_surface_prefactor_cm2_mol_s(k_1_s: float,
                                            site_density_mol_cm2: float) -> float:
    """Convert a first-order site rate into Cantera surface-reaction units.

    Args:
        k_1_s: Input controlling k 1 s.
        site_density_mol_cm2: Input controlling site density mol cm2.

    Returns:
        Validated float output for this operation.
    """
    return float(k_1_s) / float(site_density_mol_cm2)

# Surface species enthalpies must sit on Cantera's absolute scale, where the
# elements' standard states (H2, graphite) are zero. The screener reports
# adsorption energies against gas references, so each needs the reference's
# formation enthalpy added back:
#   H_s   = dE_H                 (½ H2 reference; h_f = 0)
#   CH3_s = dE_CH3 + h_f(CH3•)   (CH3 radical reference)
#   C_s   = dE_C + h_f(CH4)      (screener C reference is CH4 − 2 H2)
# Values are 298 K standard formation enthalpies (NIST / ATcT).
H_F_CH3_RADICAL_J_MOL = 145700.0
H_F_CH4_J_MOL = -74600.0
SURFACE_THERMO_REFERENCE = {
    'scale': 'cantera_absolute_elements_zero',
    'H_s': 'dE_H (reference 1/2 H2)',
    'CH3_s': 'dE_CH3 + h_f(CH3 radical) = dE_CH3 + 145.7 kJ/mol',
    'C_s': 'dE_C + h_f(CH4) = dE_C - 74.6 kJ/mol (screener C reference is CH4 - 2 H2)',
    'CH2_s': 'template: CH3_s + (C_s - CH3_s)/3',
    'CH_s': 'template: CH3_s + 2 (C_s - CH3_s)/3',
    'C_encap_s': 'same as C_s',
}


def parse_catalyst_genome(raw: Any) -> Optional[tuple]:
    """Parse a serialized or tuple catalyst genome without evaluating arbitrary code.

    Args:
        raw: Input controlling raw.

    Returns:
        Validated Optional[tuple] output for this operation.
    """
    if raw is None:
        return None
    if isinstance(raw, tuple):
        return raw
    if isinstance(raw, list):
        return tuple(raw)
    if isinstance(raw, str):
        try:
            parsed = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            return None
        if isinstance(parsed, tuple):
            return parsed
        if isinstance(parsed, list):
            return tuple(parsed)
    return None


def _nanoparticle_metals(genome: Optional[tuple], material_class: Optional[str]) -> frozenset:
    if not genome:
        return frozenset()
    cls = material_class or genome[0]
    if cls == 'SolidCatalyst' and len(genome) > 1:
        return frozenset({genome[1]})
    if cls == 'SAA' and len(genome) > 2:
        return frozenset({genome[2]})
    if cls == 'HEA' and len(genome) > 1:
        metals = genome[1]
        if isinstance(metals, (tuple, list)):
            return frozenset(metals)
        return frozenset({metals})
    return frozenset()


def off_site_carbon_allowed(genome: Any = None, material_class: Optional[str] = None) -> bool:
    """True only for nanoparticle Ni/Fe/Co (SolidCatalyst / HEA / SAA host).

    Args:
        genome: Input controlling genome.
        material_class: Input controlling material class.

    Returns:
        Validated bool output for this operation.
    """
    parsed = parse_catalyst_genome(genome)
    cls = material_class
    if parsed:
        cls = parsed[0]
    if not cls or cls in OFF_SITE_CARBON_DENIED_CLASSES:
        return False
    return bool(_nanoparticle_metals(parsed, cls) & OFF_SITE_CARBON_METALS)


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


def write_gas_only_mechanism() -> Path:
    """Write gas + condensed graphite (no surface) for equilibrium checks.

    Returns:
        Filesystem path of the written mechanism YAML.
    """
    MECHANISMS_DIR.mkdir(parents=True, exist_ok=True)

    yaml_content = f"""\
units: {{length: cm, time: s, quantity: mol, activation-energy: J/mol}}

phases:
- name: gas
  thermo: ideal-gas
  elements: [C, H, Ar]
  species: [CH4, H2, C2H2, C2H4, C2H6, Ar]
  kinetics: gas
  state: {{T: 1000.0, P: 1 atm}}

- name: graphite
  thermo: fixed-stoichiometry
  elements: [C]
  species: [C(gr)]
  state: {{T: 1000.0, P: 1 atm}}

species:
{_GAS_SPECIES_YAML}
{_GRAPHITE_SPECIES_YAML}
reactions:
- equation: 2 CH4 => C2H6 + H2
  rate-constant: {{A: 2.3e+13, b: 0.0, Ea: 356000.0}}
- equation: C2H6 => C2H4 + H2
  rate-constant: {{A: 4.65e+13, b: 0.0, Ea: 273000.0}}
- equation: C2H4 => C2H2 + H2
  rate-constant: {{A: 1.0e+14, b: 0.0, Ea: 331000.0}}
"""

    filepath = MECHANISMS_DIR / "gri30_ch4_subset.yaml"
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(yaml_content)
    logger.info(f"Wrote gas-phase + graphite mechanism: {repo_relative(filepath)}")
    return filepath


def write_full_mechanism(catalyst_name: str, E_act_CH4: float = None,
                          E_act_H_desorb: float = 0.8,
                          E_act_C_diffuse: float = 1.5,
                          site_density: float = MONOLAYER_SITE_DENSITY_MOL_CM2,
                          T_ref: float = 1000.0,
                          include_surface_sites: bool = True,
                          kinetics: CandidateKinetics = None,
                          off_site_carbon: Optional[bool] = None,
                          output_dir: Optional[Path] = None) -> Path:
    """
    Write a Cantera mechanism (gas + condensed graphite + optional surface).

    Surface carbon remains as C_s unless the class gate admits off-site
    Cγ / Cδ (nanoparticle Ni/Fe/Co only). There is no gas-phase carbon
    product. Condensed C(gr) is for multiphase equilibrium and, when
    gated, the Cγ product.

    Args:
        catalyst_name: Name used for the surface phase and output file.
        E_act_CH4: Legacy CH4 activation barrier (eV); ignored if ``kinetics``.
        E_act_H_desorb: Legacy H2 desorption barrier (eV).
        E_act_C_diffuse: Legacy carbon transfer barrier (eV).
        site_density: Site density (mol/cm^2); must equal the B1 monolayer lock.
        T_ref: Reference temperature written into the phase states (K).
        include_surface_sites: Write the Langmuir surface phase.
        kinetics: Typed candidate kinetics with provenance.
        off_site_carbon: ``None`` = class gate decides; ``True`` on a denied
            class raises; ``False`` suppresses Cγ / Cδ.
        output_dir: Directory for the YAML and kinetics sidecar. Defaults to
            the shared mechanism directory; sweep jobs supply their own.

    Returns:
        Filesystem path of the written mechanism YAML; a ``.kinetics.json``
        sidecar is written alongside it.
    """
    output_dir = MECHANISMS_DIR if output_dir is None else Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
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
    allowed = off_site_carbon_allowed(
        kinetics.genome, kinetics.material_class)
    if off_site_carbon is True and not allowed:
        raise ValueError(
            'off-site carbon is nanoparticle Ni/Fe/Co only '
            '(SolidCatalyst / HEA / SAA host); not SAC/DAC/cat_9')
    include_off_site = (
        include_surface_sites
        and (allowed if off_site_carbon is None else bool(off_site_carbon) and allowed)
    )
    values = kinetics.resolved()
    site_density = float(values['site_density_mol_cm2'])
    if abs(site_density - MONOLAYER_SITE_DENSITY_MOL_CM2) > 1e-15:
        raise ValueError(
            f'site_density={site_density} mol/cm^2; B1 locks Γ at '
            f'{MONOLAYER_SITE_DENSITY_MOL_CM2} mol/cm^2 '
            '(raise particle S/V, loading, or dispersion instead)')

    # Convert eV → J/mol. Adsorption energies alter surface enthalpies
    # but are not used as activation barriers.
    Ea_CH4 = values['methane_activation_eV'] * EV_TO_J_MOL
    Ea_CH3 = values['ch3_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_CH2 = values['ch2_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_CH = values['ch_dehydrogenation_eV'] * EV_TO_J_MOL
    Ea_H2 = values['h2_desorption_eV'] * EV_TO_J_MOL
    Ea_Cgamma = values['carbon_transfer_eV'] * EV_TO_J_MOL
    Ea_Cdelta = values['carbon_encapsulation_eV'] * EV_TO_J_MOL
    theta_star = float(values['encapsulation_crossover_coverage'])
    if not 0.0 < theta_star <= 1.0:
        raise ValueError(
            f'encapsulation_crossover_coverage={theta_star} must be in (0, 1]')
    A_Cgamma = float(values['carbon_transfer_prefactor_1_s'])
    if A_Cgamma <= 0:
        raise ValueError('carbon_transfer_prefactor_1_s must be positive')
    A_Cdelta = A_Cgamma / theta_star
    s0_ch4 = float(values['ch4_sticking_coefficient'])
    if not 0.0 < s0_ch4 <= 1.0:
        raise ValueError(
            f'ch4_sticking_coefficient={s0_ch4} must be in (0, 1]')
    A_bimol = bimolecular_surface_prefactor_cm2_mol_s(
        SURFACE_TST_PREFACTOR_1_S, site_density)
    A_h2_des = bimolecular_surface_prefactor_cm2_mol_s(
        H2_DESORPTION_PREFACTOR_1_S, site_density)
    # Adsorption energies → absolute surface enthalpies (see
    # SURFACE_THERMO_REFERENCE). Template h0 values apply when the screener
    # gave no adsorption energy; they are diagnostic, not candidate data.
    h0_h = (values['h_adsorption_eV'] * EV_TO_J_MOL
            if values['h_adsorption_eV'] is not None else -25000.0)
    h0_ch3 = (values['ch3_adsorption_eV'] * EV_TO_J_MOL + H_F_CH3_RADICAL_J_MOL
              if values['ch3_adsorption_eV'] is not None else -20000.0)
    h0_c = (values['c_adsorption_eV'] * EV_TO_J_MOL + H_F_CH4_J_MOL
            if values['c_adsorption_eV'] is not None else -40000.0)
    # CH2_s / CH_s have no screening descriptor: interpolate the
    # dehydrogenation ladder between CH3_s and C_s on the same scale.
    h0_ch2 = h0_ch3 + (h0_c - h0_ch3) / 3.0
    h0_ch = h0_ch3 + 2.0 * (h0_c - h0_ch3) / 3.0
    E_act_CH4 = float(values['methane_activation_eV'])

    if include_surface_sites:
        surf_species = (
            '[site, CH3_s, CH2_s, CH_s, H_s, C_s, C_encap_s]'
            if include_off_site else
            '[site, CH3_s, CH2_s, CH_s, H_s, C_s]'
        )
        adjacent = '[gas, graphite]' if include_off_site else '[gas]'
        phases_and_surface = f"""\
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

- name: {catalyst_name}_surface
  thermo: ideal-surface
  elements: [C, H]
  species: {surf_species}
  kinetics: surface
  reactions: [{catalyst_name}_surface-reactions]
  site-density: {site_density:.3e} mol/cm^2
  adjacent-phases: {adjacent}
"""
        surface_species = f"""\
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
    h0: {h0_ch2:.8g} J/mol
    s0: 40.0 J/mol/K
  sites: 1
  note: template enthalpy, CH3_s + (C_s - CH3_s)/3
- name: CH_s
  composition: {{C: 1, H: 1}}
  thermo:
    model: constant-cp
    h0: {h0_ch:.8g} J/mol
    s0: 30.0 J/mol/K
  sites: 1
  note: template enthalpy, CH3_s + 2 (C_s - CH3_s)/3
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
  note: Surface carbon (Cα); occupies catalytic sites until removed by policy or off-site Cγ
"""
        if include_off_site:
            surface_species += f"""\
- name: C_encap_s
  composition: {{C: 1}}
  thermo:
    model: constant-cp
    h0: {h0_c:.8g} J/mol
    s0: 10.0 J/mol/K
  sites: 1
  note: Encapsulating carbon (Cδ); occupies the site; no off-site return
"""
        off_site_rxns = ""
        if include_off_site:
            off_site_rxns = f"""\
- equation: C_s => C(gr) + site
  rate-constant: {{A: {A_Cgamma:.6g}, b: 0.0, Ea: {Ea_Cgamma:.1f}}}
  note: Cα → Cγ transport-to-edge lump (Baker / Abild-Pedersen); not nucleation
- equation: C_s => C_encap_s
  rate-constant: {{A: {A_Cdelta:.6g}, b: 0.0, Ea: {Ea_Cdelta:.1f}}}
  coverage-dependencies:
    C_s: {{a: 0.0, m: 1.0, E: 0.0}}
  note: Cα → Cδ encapsulating, rate ∝ θ_C² (Cantera k·10^(a θ)·θ^m with m=1 times [C_s]); A = A_γ/θ*, θ* = {theta_star:g}; site stays blocked
"""
        surface_rxns = f"""\
{catalyst_name}_surface-reactions:
- equation: CH4 + 2 site <=> CH3_s + H_s
  sticking-coefficient: {{A: {s0_ch4:.6g}, b: 0.0, Ea: {Ea_CH4:.1f}}}
  note: s0 sets the carbon arrival rate; provenance in the .kinetics.json sidecar
- equation: CH3_s + site <=> CH2_s + H_s
  rate-constant: {{A: {A_bimol:.6g}, b: 0.0, Ea: {Ea_CH3:.1f}}}
  note: A = 1e13/s / Γ in cm^2/mol/s (bimolecular surface TST)
- equation: CH2_s + site <=> CH_s + H_s
  rate-constant: {{A: {A_bimol:.6g}, b: 0.0, Ea: {Ea_CH2:.1f}}}
- equation: CH_s + site <=> C_s + H_s
  rate-constant: {{A: {A_bimol:.6g}, b: 0.0, Ea: {Ea_CH:.1f}}}
- equation: 2 H_s <=> H2 + 2 site
  rate-constant: {{A: {A_h2_des:.6g}, b: 0.0, Ea: {Ea_H2:.1f}}}
{off_site_rxns}"""
    else:
        phases_and_surface = f"""\
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
        surface_species = ""
        surface_rxns = ""

    yaml_content = f"""\
units: {{length: cm, time: s, quantity: mol, activation-energy: J/mol}}

phases:
{phases_and_surface}
species:
{_GAS_SPECIES_YAML}
{_GRAPHITE_SPECIES_YAML}
{surface_species}
reactions:
- equation: 2 CH4 => C2H6 + H2
  rate-constant: {{A: 2.3e+13, b: 0.0, Ea: 356000.0}}
- equation: C2H6 => C2H4 + H2
  rate-constant: {{A: 4.65e+13, b: 0.0, Ea: 273000.0}}
- equation: C2H4 => C2H2 + H2
  rate-constant: {{A: 1.0e+14, b: 0.0, Ea: 331000.0}}

{surface_rxns}
"""

    filepath = output_dir / f"mechanism_{catalyst_name}.yaml"
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(yaml_content)
    sidecar = filepath.with_suffix('.kinetics.json')
    carbon_model = (
        'condensed_graphite_plus_surface_C_s_off_site_Cgamma_Cdelta'
        if include_off_site else
        'condensed_graphite_plus_surface_C_s'
    )
    sidecar.write_text(json.dumps(_jsonable({
        'schema_version': 1,
        'catalyst_name': catalyst_name,
        'mechanism_file': repo_relative(filepath),
        'inputs': values,
        'carbon_phase_model': carbon_model,
        'surface_thermo_reference': SURFACE_THERMO_REFERENCE,
        'surface_prefactors': ({
            'bimolecular_cm2_mol_s': A_bimol,
            'h2_desorption_cm2_mol_s': A_h2_des,
            'basis': 'k_TST(1/s) / site_density(mol/cm^2)',
            'ch4_sticking_coefficient': s0_ch4,
            'ch4_sticking_source': values['provenance'].get(
                'ch4_sticking_coefficient', CH4_STICKING_PROVENANCE),
        } if include_surface_sites else {}),
        'surface_enthalpies_J_mol': ({
            'H_s': h0_h, 'CH3_s': h0_ch3, 'CH2_s': h0_ch2,
            'CH_s': h0_ch, 'C_s': h0_c,
        } if include_surface_sites else {}),
        'off_site_carbon': include_off_site,
        'coking_index_mapped_to_off_site': False,
        'off_site_channels': ({
            'C_gamma': {
                'equation': 'C_s => C(gr) + site',
                'barrier_eV': float(values['carbon_transfer_eV']),
                'preexponential_1_s': A_Cgamma,
                'preexponential_source': values['provenance'].get(
                    'carbon_transfer_prefactor_1_s', C_GAMMA_PREFACTOR_PROVENANCE),
                'particle_nm': values.get('carbon_transfer_particle_nm'),
                'form': 'first_order_theta_C',
                'kind': 'transport_to_edge',
                'source': C_GAMMA_PROVENANCE,
            },
            'C_delta': {
                'equation': 'C_s => C_encap_s',
                'barrier_eV': float(values['carbon_encapsulation_eV']),
                'preexponential_1_s': A_Cdelta,
                'form': C_DELTA_FORM,
                'crossover_coverage_theta_star': theta_star,
                'crossover_coverage_source': values['provenance'].get(
                    'encapsulation_crossover_coverage', THETA_STAR_PROVENANCE),
                'effective_crossover_note': (
                    'theta_x(T) = theta_star * exp((Ea_delta - Ea_gamma) / kT); '
                    'C_gamma wins below theta_x, C_delta above'),
                'kind': 'encapsulating',
                'source': C_DELTA_PROVENANCE,
            },
        } if include_off_site else {}),
    }), indent=2, sort_keys=True) + '\n', encoding='utf-8')

    logger.info(
        f"Wrote mechanism: {repo_relative(filepath)} "
        f"(E_act={E_act_CH4:.3f} eV, "
        f"off_site_carbon={include_off_site}, "
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
