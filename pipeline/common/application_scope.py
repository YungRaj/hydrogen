"""Shared application-specific physical scope checks for encoded candidates.

Pyrolysis reactor admissibility is phase_stable_at_application_T (ADR 0001).
The check is applied at candidate selection only. It does not change the
indexed-space coverage certificate or the 21.1B denominator.
"""

import ast


DIRECT_PEMFC_CATHODE_CLASSES = frozenset({
    'SolidCatalyst', 'HEA', 'SAC', 'DAC', 'SAA', 'Perovskite', 'Spinel',
    'MOF', 'COF', 'MetalFreeCarbon', 'MAXPhase', 'MXene',
})
OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES = frozenset({'MetalHydride', 'MoltenMetal'})

# Slab coking index ΔE_C − 2ΔE_H is defined on a solid surface slab.
# MoltenMetal has no slab; the descriptor must not rank that class.
SLAB_COKING_INDEX_CLASSES = frozenset({
    'SolidCatalyst', 'HEA', 'SAC', 'DAC', 'SAA', 'Perovskite', 'Spinel',
    'MOF', 'COF', 'MetalFreeCarbon', 'MAXPhase', 'MXene', 'MetalHydride',
})
OUT_OF_SCOPE_SLAB_COKING_CLASSES = frozenset({'MoltenMetal'})

APPLICATION_PYROLYSIS = 'turquoise_pyrolysis'
APPLICATION_PEMFC = 'pemfc_cathode'
TURQUOISE_PYROLYSIS_MIN_T_K = 700.0

# Dead under both applications: keep enumerated for coverage; no reserved
# eSen/DFT validation slots (ADR 0001).
VALIDATION_QUOTA_EXEMPT_CLASSES = frozenset({'MetalHydride'})

# Encoded phase is not present at pyrolysis T. MetalHydride is first, not only.
# Basis: docs/adr/0001-pyrolysis-phase-admissibility.md
_PYROLYSIS_PHASE_STABILITY = {
    'MetalHydride': {
        'status': 'out_of_scope',
        'reason': (
            'encoded hydride phase is not present under flowing CH4 at '
            f'>=~{int(TURQUOISE_PYROLYSIS_MIN_T_K)} K (P_H2~0; residual '
            'metal carburizes; leftover metal is not this genome)'),
    },
    'MOF': {
        'status': 'out_of_scope',
        'reason': (
            'encoded MOF framework does not survive 773-1300 K; residue is a '
            'carbonized MOF-derived composite (precursor, not the genome)'),
    },
    'COF': {
        'status': 'out_of_scope',
        'reason': (
            'encoded COF does not survive 773-1300 K; residue is carbonaceous, '
            'not the covalent framework genome'),
    },
    'MXene': {
        'status': 'out_of_scope',
        'reason': (
            'MXene terminations are lost and Ti3C2Tx converts toward TiC; '
            'the encoded MXene is not the high-T solid'),
    },
    'Perovskite': {
        'status': 'out_of_scope',
        'reason': (
            'under reducing CH4 the perovskite exsolves metal; working '
            'catalyst is not the encoded ABO3 genome'),
    },
}

# Backward-compatible alias used by older call sites / logs.
OUT_OF_SCOPE_TURQUOISE_PYROLYSIS_CLASSES = frozenset(
    cls for cls, spec in _PYROLYSIS_PHASE_STABILITY.items()
    if spec['status'] == 'out_of_scope'
)


def pemfc_cathode_scope(genome: tuple) -> dict:
    """Reject classes with no physically defined solid PEMFC cathode realization.

    Args:
        genome: Encoded catalyst composition and structural configuration.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    material_class = genome[0]
    if material_class in OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES:
        return {'status': 'out_of_scope', 'reason':
                f'{material_class} has no encoded solid catalyst-layer realization'}
    if material_class not in DIRECT_PEMFC_CATHODE_CLASSES:
        return {'status': 'unknown', 'reason': 'unrecognized material class'}
    return {'status': 'candidate', 'reason': None}


def slab_coking_index_scope(genome: tuple) -> dict:
    """Reject slab coking_index for classes with no solid-slab realization."""
    material_class = genome[0] if genome else None
    if material_class in OUT_OF_SCOPE_SLAB_COKING_CLASSES:
        return {
            'status': 'out_of_scope',
            'reason': (
                f'{material_class} has no slab; coking_index=ΔE_C-2ΔE_H '
                'is not a melt carbon-separation descriptor'),
        }
    if material_class not in SLAB_COKING_INDEX_CLASSES:
        return {'status': 'unknown', 'reason': 'unrecognized material class'}
    return {'status': 'candidate', 'reason': None}


def phase_stable_at_application_T(genome: tuple, T_K: float = None,
                                  application: str = APPLICATION_PYROLYSIS) -> dict:
    """
    Encoded-phase admissibility at the application's operating temperature.

    See ADR 0001. Applied at candidate selection only — not to the coverage
    certificate denominator.
    """
    material_class = genome[0] if genome else None
    if application == APPLICATION_PEMFC:
        return pemfc_cathode_scope(genome if genome else (None,))
    if application != APPLICATION_PYROLYSIS:
        return {'status': 'unknown', 'reason': f'unrecognized application {application!r}'}

    spec = _PYROLYSIS_PHASE_STABILITY.get(material_class)
    if spec is not None:
        return {'status': spec['status'], 'reason': spec['reason']}
    if T_K is not None and T_K < TURQUOISE_PYROLYSIS_MIN_T_K:
        return {
            'status': 'unknown',
            'reason': f'T={T_K} K below turquoise pyrolysis band',
        }
    if not material_class:
        return {'status': 'unknown', 'reason': 'missing material class'}
    return {'status': 'candidate', 'reason': None}


def turquoise_pyrolysis_scope(genome: tuple, T_K: float = None) -> dict:
    """Admissibility for turquoise methane-pyrolysis ranking / Phase 2."""
    return phase_stable_at_application_T(genome, T_K, APPLICATION_PYROLYSIS)


def parse_encoded_genome(raw):
    """Parse a screening-row genome (tuple or literal string) or return None."""
    if isinstance(raw, str):
        try:
            raw = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            return None
    if isinstance(raw, list):
        raw = tuple(raw)
    if not isinstance(raw, tuple) or not raw:
        return None
    return raw


def is_turquoise_pyrolysis_candidate(raw, T_K: float = None) -> bool:
    """True only if the genome is admissible for turquoise pyrolysis reactors."""
    genome = parse_encoded_genome(raw)
    if genome is None:
        return False
    return phase_stable_at_application_T(genome, T_K)['status'] == 'candidate'


def select_turquoise_pyrolysis_candidates(df, top_k=None, genome_col: str = 'genome'):
    """
    Drop phase-unstable rows then optionally take top_k by E_act.

    Candidate-selection filter only. Does not alter coverage certificates.
    """
    if df is None:
        raise ValueError('screening frame is required for turquoise pyrolysis scope')
    if genome_col not in df.columns:
        raise ValueError(
            f'{genome_col!r} is required to apply phase_stable_at_application_T')
    scoped = df[df[genome_col].map(is_turquoise_pyrolysis_candidate)].copy()
    if top_k is None or len(scoped) == 0:
        return scoped
    if 'E_act' in scoped.columns:
        return scoped.nsmallest(int(top_k), 'E_act')
    return scoped.head(int(top_k))


def scope_pyrolysis_pool(df, genome_col: str = 'genome'):
    """Admissibility pool for slate drawing, tolerant of genome-less frames.

    Returns ``(pool, note)``. Rows whose encoded phase is unstable at the
    pyrolysis T band are dropped (ADR 0001). Validity is *not* applied here:
    the validation route may still rescue invalid rows, and the caller
    applies ``valid`` for the reactor route. A frame without a genome column
    (synthetic tables, fixtures) is returned unchanged with the reason
    recorded so the omission is visible rather than silent.
    """
    if df is None:
        raise ValueError('screening frame is required for turquoise pyrolysis scope')
    if genome_col not in df.columns:
        return df, {'filter': None,
                    'reason': f'{genome_col} column absent; admissibility not applied'}
    pool = select_turquoise_pyrolysis_candidates(df, top_k=None, genome_col=genome_col)
    return pool, {'filter': 'phase_stable_at_application_T',
                  'admissible_count': int(len(pool)),
                  'dropped_count': int(len(df) - len(pool))}


def is_validation_quota_class(material_class: str, application: str = None) -> bool:
    """Reserved eSen/DFT slots: dead classes stay enumerated, they do not get a quota.

    MetalHydride is exempt on every application. Pyrolysis additionally
    exempts every class whose encoded phase is out of scope (ADR 0001).
    """
    if material_class in VALIDATION_QUOTA_EXEMPT_CLASSES:
        return False
    if application in (APPLICATION_PYROLYSIS, 'pyrolysis', 'turquoise_hydrogen'):
        return material_class not in OUT_OF_SCOPE_TURQUOISE_PYROLYSIS_CLASSES
    return True


def validation_quota_class_count(material_classes, application: str = None) -> int:
    """How many classes consume --min-validation-per-class reserved slots."""
    return sum(
        1 for cls in material_classes
        if is_validation_quota_class(cls, application)
    )
