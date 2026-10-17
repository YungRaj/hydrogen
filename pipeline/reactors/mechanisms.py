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
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np

from pipeline.reactors.mechanism_yaml import (
    _GAS_SPECIES_YAML,
    _GRAPHITE_SPECIES_YAML,
    _carbon_removal_reactions_yaml,
    _gas_graphite_phases_yaml,
    _jsonable,
    _surface_phase_yaml,
    _surface_reactions_yaml,
    _surface_species_yaml,
    _write_kinetics_sidecar as _render_kinetics_sidecar,
)
from pipeline.utils import (
    MECHANISMS_DIR,
    eV_to_J,
    repo_relative,
    setup_logger,
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
    encapsulation_crossover_coverage: Optional[float] = None
    carbon_transfer_prefactor_1_s: Optional[float] = None
    carbon_transfer_particle_nm: Optional[float] = None
    carbon_diffusion_prefactor_m2_s: Optional[float] = None
    ch4_sticking_coefficient: Optional[float] = None
    site_density_mol_cm2: float = 2.5e-9
    screening_protocol: str = 'unknown'
    candidate_id: str = 'unknown'
    material_class: Optional[str] = None
    genome: Optional[tuple] = None
    sources: Mapping[str, str] = field(default_factory=dict)

    @classmethod
    def from_screening_row(
        cls,
        row: Mapping[str, object],
        candidate_id: str = 'unknown',
        validation: Optional[Mapping[str, object]] = None,
    ) -> 'CandidateKinetics':
        """Build kinetics from screening and optional converged NEB evidence.

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
        neb_barrier = finite('E_act_neb')
        if neb_barrier is not None and neb_barrier > 0:
            barrier = neb_barrier
            sources['methane_activation_eV'] = f'mlip_neb:{protocol}:E_act_neb'
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
        step_columns = (
            'ch3_dehydrogenation_eV',
            'ch2_dehydrogenation_eV',
            'ch_dehydrogenation_eV',
            'h2_desorption_eV',
            'carbon_transfer_eV',
            'carbon_encapsulation_eV',
        )
        for step_col in step_columns:
            step_val = finite(step_col)
            if step_val is not None and step_val >= 0:
                values[step_col] = step_val
                sources[step_col] = f'mlip_neb:{protocol}:{step_col}'
        genome = parse_catalyst_genome(row.get('genome'))
        material_class = row.get('material_class')
        if material_class is None or (
            isinstance(material_class, float) and not np.isfinite(material_class)
        ):
            material_class = genome[0] if genome else None
        elif material_class is not None:
            material_class = str(material_class)
        if validation is not None:
            if str(validation.get('candidate_id')) != str(candidate_id):
                raise ValueError('kinetics validation candidate_id mismatch')
            evidence_level = validation.get('evidence_level')
            if not validation.get('complete') or evidence_level not in (
                'converged_dft_neb_frequency',
                'converged_mlip_neb',
            ):
                raise ValueError('kinetics validation campaign is incomplete')
            source_prefix = (
                'candidate_specific'
                if evidence_level == 'converged_dft_neb_frequency'
                else 'mlip_neb'
            )
            resolved = validation.get('resolved_kinetics_eV', {})
            allowed = {
                'methane_activation_eV',
                'ch3_dehydrogenation_eV',
                'ch2_dehydrogenation_eV',
                'ch_dehydrogenation_eV',
                'h2_desorption_eV',
                'carbon_transfer_eV',
                'carbon_encapsulation_eV',
            }
            unknown = set(resolved) - allowed
            if unknown:
                raise ValueError(
                    f'unknown validated kinetics fields: {sorted(unknown)}'
                )
            for name, raw_value in resolved.items():
                value = float(raw_value)
                if not np.isfinite(value) or value < 0:
                    raise ValueError(f'invalid validated barrier for {name}')
                if name == 'methane_activation_eV':
                    barrier = value
                else:
                    values[name] = value
                sources[name] = f'{source_prefix}:{evidence_level}'
        return cls(
            methane_activation_eV=barrier,
            candidate_id=candidate_id,
            screening_protocol=protocol,
            sources=sources,
            material_class=material_class,
            genome=genome,
            **values,
        )

    def resolved(self) -> dict:
        """Return numerical values plus whether each was observed or templated.

        Returns:
            Dictionary of resolved barriers, per-field ``provenance``, and
            ``quantitative_status`` (``candidate_specific`` or
            ``screening_complete``).
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
                DEFAULT_ENCAPSULATION_CROSSOVER_COVERAGE
            )
            provenance['encapsulation_crossover_coverage'] = THETA_STAR_PROVENANCE
        else:
            provenance.setdefault(
                'encapsulation_crossover_coverage', THETA_STAR_PROVENANCE
            )
        particle_nm = values['carbon_transfer_particle_nm']
        if particle_nm is not None:
            if values['carbon_transfer_prefactor_1_s'] is not None:
                raise ValueError(
                    'set carbon_transfer_particle_nm or '
                    'carbon_transfer_prefactor_1_s, not both'
                )
            if not particle_nm > 0:
                raise ValueError('carbon_transfer_particle_nm must be positive')
            d0 = values['carbon_diffusion_prefactor_m2_s']
            if d0 is None:
                d0 = CARBON_DIFFUSION_PREFACTOR_M2_S
                values['carbon_diffusion_prefactor_m2_s'] = d0
                provenance['carbon_diffusion_prefactor_m2_s'] = (
                    CARBON_DIFFUSION_PREFACTOR_PROVENANCE
                )
            elif not d0 > 0:
                raise ValueError('carbon_diffusion_prefactor_m2_s must be positive')
            else:
                provenance.setdefault(
                    'carbon_diffusion_prefactor_m2_s', 'declared_not_measured'
                )
            values['carbon_transfer_prefactor_1_s'] = (
                carbon_transfer_prefactor_from_particle(particle_nm, d0)
            )
            provenance['carbon_transfer_prefactor_1_s'] = (
                f'derived: D0/L^2 with D0={d0:.3g} m^2/s, L={particle_nm:g} nm'
            )
        elif values['carbon_transfer_prefactor_1_s'] is None:
            values['carbon_transfer_prefactor_1_s'] = OFF_SITE_PREEXPONENTIAL_1_S
            provenance['carbon_transfer_prefactor_1_s'] = C_GAMMA_PREFACTOR_PROVENANCE
        else:
            provenance.setdefault(
                'carbon_transfer_prefactor_1_s', 'declared_not_measured'
            )
        if values['ch4_sticking_coefficient'] is None:
            values['ch4_sticking_coefficient'] = DEFAULT_CH4_STICKING_COEFFICIENT
            provenance['ch4_sticking_coefficient'] = CH4_STICKING_PROVENANCE
        else:
            provenance.setdefault('ch4_sticking_coefficient', 'declared_not_measured')
        values['provenance'] = provenance
        required = set(defaults)
        if not off_site_carbon_allowed(self.genome, self.material_class):
            required.discard('carbon_encapsulation_eV')
        if any(provenance.get(name) == 'template_default' for name in required):
            values['quantitative_status'] = 'screening_complete'
        elif any(
            str(provenance.get(name, '')).startswith('mlip_neb:')
            for name in required
        ):
            values['quantitative_status'] = 'mlip_neb_complete'
        else:
            values['quantitative_status'] = 'candidate_specific'
        return values


MONOLAYER_SITE_DENSITY_MOL_CM2 = 2.5e-9
OFF_SITE_CARBON_METALS = frozenset({'Ni', 'Fe', 'Co'})
OFF_SITE_CARBON_DENIED_CLASSES = frozenset(
    {
        'SAC',
        'DAC',
        'MetalFreeCarbon',
        'MoltenMetal',
        'MOF',
        'COF',
        'Perovskite',
        'MetalHydride',
        'MXene',
        'MAXPhase',
        'Spinel',
    }
)
C_GAMMA_PROVENANCE = 'template_default: Abild-Pedersen 2006 / Baker 1972'
C_DELTA_PROVENANCE = (
    'template_default: Amin IEC Res 2011 50 12460; '
    '147-149 kJ/mol encapsulating-carbon'
)
DEFAULT_ENCAPSULATION_CROSSOVER_COVERAGE = 0.5
THETA_STAR_PROVENANCE = 'declared_not_measured: swept in B6-6'
C_DELTA_FORM = 'coverage_dependent_theta_C_squared'
OFF_SITE_PREEXPONENTIAL_1_S = 1.0e13
C_GAMMA_PREFACTOR_PROVENANCE = (
    'template_default: single-hop TST 1e13/s; transport+precipitation lump '
    'is D0/L^2 and particle-size dependent; swept in B6-6'
)
CARBON_DIFFUSION_PREFACTOR_M2_S = 2.48e-4
CARBON_DIFFUSION_PREFACTOR_PROVENANCE = (
    'template_default: Lander J.Appl.Phys. 1952 23 1305; '
    'D0 = 2.48 cm^2/s (ln D = 0.909 - 20200/T)'
)
DEFAULT_CH4_STICKING_COEFFICIENT = 0.01
CH4_STICKING_PROVENANCE = (
    'template_default: s0 = 0.01 (methane_pox_on_pt order of magnitude); '
    'swept in B6-6'
)
SURFACE_TST_PREFACTOR_1_S = 1.0e13
H2_DESORPTION_PREFACTOR_1_S = 5.0e13
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


def carbon_transfer_prefactor_from_particle(
    particle_nm: float, d0_m2_s: float
) -> float:
    """A_γ = D0 / L² (1/s) for a metal particle of diameter ``particle_nm``.

    Args:
        particle_nm: Particle diameter in nanometers.
        d0_m2_s: Input controlling d0 m2 s.

    Returns:
        Validated float output for this operation.
    """
    length_m = float(particle_nm) * 1e-9
    return float(d0_m2_s) / (length_m * length_m)


def bimolecular_surface_prefactor_cm2_mol_s(
    k_1_s: float, site_density_mol_cm2: float
) -> float:
    """Convert a first-order site rate into Cantera surface-reaction units.

    Args:
        k_1_s: Input controlling k 1 s.
        site_density_mol_cm2: Input controlling site density mol cm2.

    Returns:
        Validated float output for this operation.
    """
    return float(k_1_s) / float(site_density_mol_cm2)


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


def _nanoparticle_metals(
    genome: Optional[tuple], material_class: Optional[str]
) -> frozenset:
    if not genome:
        return frozenset()
    cls = material_class or genome[0]
    if cls in ('SolidCatalyst', 'SupportedNanoparticle') and len(genome) > 1:
        return frozenset({genome[1]})
    if cls == 'SAA' and len(genome) > 2:
        return frozenset({genome[2]})
    if cls == 'HEA' and len(genome) > 1:
        metals = genome[1]
        if isinstance(metals, (tuple, list)):
            return frozenset(metals)
        return frozenset({metals})
    return frozenset()


def off_site_carbon_allowed(
    genome: Any = None, material_class: Optional[str] = None
) -> bool:
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


def _surface_enthalpies(values: dict) -> dict:
    """Absolute surface enthalpies (J/mol) on the SURFACE_THERMO_REFERENCE scale."""
    def adsorbed(key: str, offset: float | None, template: float) -> float:
        if values[key] is None:
            return template
        energy = values[key] * EV_TO_J_MOL
        return energy if offset is None else energy + offset

    h0_h = adsorbed('h_adsorption_eV', None, -25000.0)
    h0_ch3 = adsorbed('ch3_adsorption_eV', H_F_CH3_RADICAL_J_MOL, -20000.0)
    h0_c = adsorbed('c_adsorption_eV', H_F_CH4_J_MOL, -40000.0)
    return {
        'H_s': h0_h,
        'CH3_s': h0_ch3,
        'CH2_s': h0_ch3 + (h0_c - h0_ch3) / 3.0,
        'CH_s': h0_ch3 + 2.0 * (h0_c - h0_ch3) / 3.0,
        'C_s': h0_c,
    }


def _surface_rate_parameters(values: dict, site_density: float) -> dict:
    """Validate kinetic inputs and convert barriers to J/mol and prefactors."""
    theta_star = float(values['encapsulation_crossover_coverage'])
    if not 0.0 < theta_star <= 1.0:
        raise ValueError(
            f'encapsulation_crossover_coverage={theta_star} must be in (0, 1]'
        )
    A_Cgamma = float(values['carbon_transfer_prefactor_1_s'])
    if A_Cgamma <= 0:
        raise ValueError('carbon_transfer_prefactor_1_s must be positive')
    s0_ch4 = float(values['ch4_sticking_coefficient'])
    if not 0.0 < s0_ch4 <= 1.0:
        raise ValueError(f'ch4_sticking_coefficient={s0_ch4} must be in (0, 1]')
    h0 = _surface_enthalpies(values)
    dh_ch4 = max(0.0, h0['CH3_s'] + h0['H_s'] - H_F_CH4_J_MOL)
    dh_ch3 = max(0.0, h0['CH2_s'] + h0['H_s'] - h0['CH3_s'])
    dh_ch2 = max(0.0, h0['CH_s'] + h0['H_s'] - h0['CH2_s'])
    dh_ch = max(0.0, h0['C_s'] + h0['H_s'] - h0['CH_s'])
    dh_h2 = max(0.0, -2.0 * h0['H_s'])
    return {
        'Ea_CH4': max(values['methane_activation_eV'] * EV_TO_J_MOL, dh_ch4),
        'Ea_CH3': max(values['ch3_dehydrogenation_eV'] * EV_TO_J_MOL, dh_ch3),
        'Ea_CH2': max(values['ch2_dehydrogenation_eV'] * EV_TO_J_MOL, dh_ch2),
        'Ea_CH': max(values['ch_dehydrogenation_eV'] * EV_TO_J_MOL, dh_ch),
        'Ea_H2': max(values['h2_desorption_eV'] * EV_TO_J_MOL, dh_h2),
        'Ea_Cgamma': values['carbon_transfer_eV'] * EV_TO_J_MOL,
        'Ea_Cdelta': values['carbon_encapsulation_eV'] * EV_TO_J_MOL,
        'theta_star': theta_star,
        'A_Cgamma': A_Cgamma,
        'A_Cdelta': A_Cgamma / theta_star,
        's0_ch4': s0_ch4,
        'A_bimol': bimolecular_surface_prefactor_cm2_mol_s(
            SURFACE_TST_PREFACTOR_1_S, site_density
        ),
        'A_h2_des': bimolecular_surface_prefactor_cm2_mol_s(
            H2_DESORPTION_PREFACTOR_1_S, site_density
        ),
    }


def _write_kinetics_sidecar(
    filepath: Path,
    catalyst_name: str,
    values: dict,
    rates: dict,
    h0: dict,
    include_surface_sites: bool,
    include_off_site: bool,
) -> None:
    """Write the provenance sidecar recording every kinetic input and prefactor."""
    _render_kinetics_sidecar(
        filepath,
        catalyst_name,
        values,
        rates,
        h0,
        include_surface_sites,
        include_off_site,
        surface_thermo_reference=SURFACE_THERMO_REFERENCE,
        ch4_sticking_provenance=CH4_STICKING_PROVENANCE,
        c_gamma_prefactor_provenance=C_GAMMA_PREFACTOR_PROVENANCE,
        c_gamma_provenance=C_GAMMA_PROVENANCE,
        c_delta_form=C_DELTA_FORM,
        theta_star_provenance=THETA_STAR_PROVENANCE,
        c_delta_provenance=C_DELTA_PROVENANCE,
    )


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


def write_full_mechanism(
    catalyst_name: str,
    E_act_CH4: float = None,
    E_act_H_desorb: float = 0.8,
    E_act_C_diffuse: float = 1.5,
    site_density: float = MONOLAYER_SITE_DENSITY_MOL_CM2,
    T_ref: float = 1000.0,
    include_surface_sites: bool = True,
    kinetics: CandidateKinetics = None,
    off_site_carbon: Optional[bool] = None,
    carbon_regeneration_rate_1_s: Optional[float] = None,
    output_dir: Optional[Path] = None,
) -> Path:
    """Write a Cantera mechanism (gas + condensed graphite + optional surface).

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
        carbon_regeneration_rate_1_s: Optional mechanically estimated C_s
            regeneration rate.
        output_dir: Directory for the YAML and kinetics sidecar.

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
            sources={
                'methane_activation_eV': 'legacy_argument',
                'h2_desorption_eV': 'legacy_argument',
                'carbon_transfer_eV': 'legacy_argument',
            },
        )
    allowed = off_site_carbon_allowed(kinetics.genome, kinetics.material_class)
    if off_site_carbon is True and not allowed:
        raise ValueError(
            'off-site carbon is nanoparticle Ni/Fe/Co only '
            '(SolidCatalyst / HEA / SAA host); not SAC/DAC/cat_9'
        )
    include_off_site = include_surface_sites and (
        allowed if off_site_carbon is None else bool(off_site_carbon) and allowed
    )
    values = kinetics.resolved()
    site_density = float(values['site_density_mol_cm2'])
    if abs(site_density - MONOLAYER_SITE_DENSITY_MOL_CM2) > 1e-15:
        raise ValueError(
            f'site_density={site_density} mol/cm^2; B1 locks Γ at '
            f'{MONOLAYER_SITE_DENSITY_MOL_CM2} mol/cm^2 '
            '(raise particle S/V, loading, or dispersion instead)'
        )

    rates = _surface_rate_parameters(values, site_density)
    h0 = _surface_enthalpies(values)
    E_act_CH4 = float(values['methane_activation_eV'])

    if include_surface_sites:
        if (
            carbon_regeneration_rate_1_s is not None
            and not include_off_site
            and carbon_regeneration_rate_1_s <= 0.0
        ):
            raise ValueError('carbon regeneration rate must be positive')
        phases = _gas_graphite_phases_yaml(T_ref) + '\n' + _surface_phase_yaml(
            catalyst_name, site_density, include_off_site
        )
        surface_species = _surface_species_yaml(h0, include_off_site)
        surface_rxns = _surface_reactions_yaml(
            catalyst_name,
            rates,
            _carbon_removal_reactions_yaml(
                rates, include_off_site, carbon_regeneration_rate_1_s
            ),
        )
    else:
        phases = _gas_graphite_phases_yaml(T_ref)
        surface_species = ""
        surface_rxns = ""

    yaml_content = f"""\
units: {{length: cm, time: s, quantity: mol, activation-energy: J/mol}}

phases:
{phases}
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
    _write_kinetics_sidecar(
        filepath,
        catalyst_name,
        values,
        rates,
        h0,
        include_surface_sites,
        include_off_site,
    )

    logger.info(
        f"Wrote mechanism: {repo_relative(filepath)} "
        f"(E_act={E_act_CH4:.3f} eV, "
        f"off_site_carbon={include_off_site}, "
        f"status={values['quantitative_status']})"
    )
    return filepath


# Aliases
write_gri30_subset = write_gas_only_mechanism


if __name__ == '__main__':
    write_gas_only_mechanism()
    write_full_mechanism("NiBi_10pct", E_act_CH4=0.85, T_ref=1000)
    write_full_mechanism("FeNi_graphene", E_act_CH4=0.65, T_ref=900)
    write_full_mechanism("CuSn_20pct", E_act_CH4=1.10, T_ref=1100)
    print("Mechanism files written to:", MECHANISMS_DIR)
