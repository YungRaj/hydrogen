"""Shared application-specific physical scope checks for encoded candidates."""


DIRECT_PEMFC_CATHODE_CLASSES = frozenset({
    'SolidCatalyst', 'HEA', 'SAC', 'DAC', 'SAA', 'Perovskite', 'Spinel',
    'MOF', 'COF', 'MetalFreeCarbon', 'MAXPhase', 'MXene',
})
OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES = frozenset({'MetalHydride', 'MoltenMetal'})


def pemfc_cathode_scope(genome: tuple) -> dict:
    """Reject classes with no physically defined solid PEMFC cathode realization."""
    material_class = genome[0]
    if material_class in OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES:
        return {'status': 'out_of_scope', 'reason':
                f'{material_class} has no encoded solid catalyst-layer realization'}
    if material_class not in DIRECT_PEMFC_CATHODE_CLASSES:
        return {'status': 'unknown', 'reason': 'unrecognized material class'}
    return {'status': 'candidate', 'reason': None}
