"""Canonical methane-conversion pathway modes and reactor routing.

Mode selection is intentionally separate from catalyst scoring.  A mode chooses
the physical reactor/pathway model that is allowed to consume a candidate; it
must never silently grant a candidate a performance bonus.
"""

from __future__ import annotations

from dataclasses import dataclass


DEFAULT_MODE = 'thermocatalytic'
MODE_CHOICES = (
    'thermocatalytic',
    'thermocatalytic_pfr',
    'thermocatalytic_fluidized',
    'mmbcr',
    'ntec',
    'electrochemical',
)


@dataclass(frozen=True)
class PathwayMode:
    name: str
    reactor_types: tuple[str, ...]
    requires_specialized_validation: bool = False


@dataclass(frozen=True)
class ReactorModel:
    """Physical interpretation and current fidelity of one dispatch target."""
    reactor_type: str
    bed_or_interface: str
    cantera_model: str | None
    reaction_domain: str
    compatible_material_classes: frozenset[str] | None
    fidelity: str


_MODES = {
    'thermocatalytic': PathwayMode(
        'thermocatalytic', ('PFR', 'Fluidized')),
    'thermocatalytic_pfr': PathwayMode(
        'thermocatalytic_pfr', ('PFR',)),
    'thermocatalytic_fluidized': PathwayMode(
        'thermocatalytic_fluidized', ('Fluidized',)),
    'mmbcr': PathwayMode('mmbcr', ('MMBCR',)),
    'ntec': PathwayMode('ntec', ('NTEC',), True),
    'electrochemical': PathwayMode(
        'electrochemical', ('Electrochemical',), True),
}

_SOLID_CLASSES = frozenset({
    'SolidCatalyst', 'SAC', 'DAC', 'MOF', 'COF', 'Perovskite',
    'MetalHydride', 'MAXPhase', 'HEA', 'Spinel', 'MXene', 'SAA',
    'MetalFreeCarbon',
})

REACTOR_MODELS = {
    'PFR': ReactorModel(
        'PFR', 'fixed_packed_bed', 'staged_lagrangian_ideal_gas_reactors',
        'heterogeneous_gas_solid_surface', _SOLID_CLASSES,
        'screening_approximation'),
    'Fluidized': ReactorModel(
        'Fluidized', 'gas_solid_fluidized_bed',
        'emulsion_reactor_plus_bubble_bypass',
        'heterogeneous_gas_solid_surface', _SOLID_CLASSES,
        'screening_approximation'),
    'MMBCR': ReactorModel(
        'MMBCR', 'gas_bubbles_in_molten_metal',
        'steady_cstrs_in_series_with_interfacial_surface',
        'gas_liquid_interface_represented_as_ideal_surface_proxy',
        frozenset({'MoltenMetal'}), 'screening_approximation'),
    'NTEC': ReactorModel(
        'NTEC', 'mechanically_agitated_liquid_solid_interface', None,
        'coupled_mechanical_electrical_liquid_solid_pathway', None,
        'validated_multiphysics_artifact_required'),
    'Electrochemical': ReactorModel(
        'Electrochemical', 'electrode_electrolyte_interface', None,
        'electrode_electrolyte_charge_transfer_pathway', None,
        'validated_multiphysics_artifact_required'),
}


def resolve_pathway_mode(mode: str | None) -> PathwayMode:
    """Return a validated mode definition; absent mode means thermocatalytic."""
    normalized = (mode or DEFAULT_MODE).strip().lower()
    try:
        return _MODES[normalized]
    except KeyError as exc:
        raise ValueError(
            f'unknown pathway mode {mode!r}; expected one of {MODE_CHOICES}') from exc


def reactor_types_for_mode(mode: str | None) -> tuple[str, ...]:
    return resolve_pathway_mode(mode).reactor_types


def validate_mode_reactors(mode: str | None,
                           reactor_types: tuple[str, ...] | list[str]) -> None:
    """Reject cross-wiring such as an MMBCR mode dispatched into a PFR."""
    selected = resolve_pathway_mode(mode)
    requested = tuple(reactor_types)
    unknown = set(requested) - set(REACTOR_MODELS)
    if unknown:
        raise ValueError(f'unknown reactor type(s): {sorted(unknown)}')
    if requested != selected.reactor_types:
        raise ValueError(
            f'mode {selected.name!r} requires reactors '
            f'{selected.reactor_types}, received {requested}')


def reactor_applicability(reactor_type: str,
                          material_class: str | None) -> tuple[bool, str | None]:
    """Return physical bed/material compatibility without excluding candidates."""
    spec = REACTOR_MODELS[reactor_type]
    if spec.compatible_material_classes is None:
        return True, None
    if not material_class:
        return False, 'material_class_required_for_reactor_applicability'
    if material_class not in spec.compatible_material_classes:
        return False, (
            f'{material_class}_is_not_compatible_with_{spec.bed_or_interface}')
    return True, None


def is_ntec_mode(mode: str | None) -> bool:
    return resolve_pathway_mode(mode).name == 'ntec'
