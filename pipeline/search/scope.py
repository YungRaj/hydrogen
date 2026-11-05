"""Shared application-specific physical scope checks for encoded candidates.

Pyrolysis reactor admissibility is phase_stable_at_application_T (ADR 0001).
The check is applied at candidate selection only. It does not change the
indexed-space coverage certificate or the 21.1B denominator.
"""

import ast
from typing import TypeVar

from pipeline.data_models.discovery import AdmissibilitySummary


FrameT = TypeVar('FrameT')


DIRECT_PEMFC_CATHODE_CLASSES = frozenset(
    {
        'SolidCatalyst',
        'HEA',
        'SAC',
        'DAC',
        'SAA',
        'Perovskite',
        'Spinel',
        'MOF',
        'COF',
        'MetalFreeCarbon',
        'MAXPhase',
        'MXene',
    }
)
OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES = frozenset({'MetalHydride', 'MoltenMetal'})

# Slab coking index ΔE_C − 2ΔE_H is defined on a solid surface slab.
# MoltenMetal has no slab; the descriptor must not rank that class.
SLAB_COKING_INDEX_CLASSES = DIRECT_PEMFC_CATHODE_CLASSES | {
    'MetalHydride',
    'SupportedNanoparticle',
}
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
            'metal carburizes; leftover metal is not this genome)'
        ),
    },
    'MOF': {
        'status': 'out_of_scope',
        'reason': (
            'encoded MOF framework does not survive 773-1300 K; residue is a '
            'carbonized MOF-derived composite (precursor, not the genome)'
        ),
    },
    'COF': {
        'status': 'out_of_scope',
        'reason': (
            'encoded COF does not survive 773-1300 K; residue is carbonaceous, '
            'not the covalent framework genome'
        ),
    },
    'MXene': {
        'status': 'out_of_scope',
        'reason': (
            'MXene terminations are lost and Ti3C2Tx converts toward TiC; '
            'the encoded MXene is not the high-T solid'
        ),
    },
    'Perovskite': {
        'status': 'out_of_scope',
        'reason': (
            'under reducing CH4 the perovskite exsolves metal; working '
            'catalyst is not the encoded ABO3 genome'
        ),
    },
}

# Backward-compatible alias used by older call sites / logs.
OUT_OF_SCOPE_TURQUOISE_PYROLYSIS_CLASSES = frozenset(
    cls
    for cls, spec in _PYROLYSIS_PHASE_STABILITY.items()
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
        return {
            'status': 'out_of_scope',
            'reason': f'{material_class} has no encoded solid catalyst-layer realization',
        }
    if material_class not in DIRECT_PEMFC_CATHODE_CLASSES:
        return {'status': 'unknown', 'reason': 'unrecognized material class'}
    host_idx = 2 if material_class == 'SAA' else 1
    if (
        material_class in ('SolidCatalyst', 'SupportedNanoparticle', 'SAA')
        and len(genome) > host_idx
        and str(genome[host_idx]) in LOW_MELTING_HOST_METALS_K
    ):
        host = str(genome[host_idx])
        return {
            'status': 'out_of_scope',
            'reason': (
                f'{material_class} host {host} is liquid or rapidly acid-soluble '
                'under PEMFC cathode operating conditions (353.15 K, pH ~ 0)'
            ),
        }
    return {'status': 'candidate', 'reason': None}


def slab_coking_index_scope(genome: tuple) -> dict:
    """Reject slab coking_index for classes with no solid-slab realization.

    Args:
        genome: Input controlling genome.

    Returns:
        Validated dict output for this operation.
    """
    material_class = genome[0] if genome else None
    if material_class in OUT_OF_SCOPE_SLAB_COKING_CLASSES:
        return {
            'status': 'out_of_scope',
            'reason': (
                f'{material_class} has no slab; coking_index=ΔE_C-2ΔE_H '
                'is not a melt carbon-separation descriptor'
            ),
        }
    if material_class not in SLAB_COKING_INDEX_CLASSES:
        return {'status': 'unknown', 'reason': 'unrecognized material class'}
    return {'status': 'candidate', 'reason': None}


LOW_MELTING_HOST_METALS_K = {
    'Ga': 302.9,
    'In': 429.8,
    'Sn': 505.1,
    'Bi': 544.7,
    'Pb': 600.6,
    'Zn': 692.7,
    'Mg': 923.0,
    'Al': 933.5,
    'Ba': 1000.0,
    'Sr': 1050.0,
    'Ce': 1068.0,
    'Yb': 1097.0,
    'Eu': 1099.0,
    'Ca': 1115.0,
    'La': 1193.0,
    'Pr': 1208.0,
}

PYROLYSIS_VOLATILE_ELEMENTS = frozenset(
    {
        'Zn',
        'Te',
        'Cd',
        'Hg',
        'Se',
        'As',
        'Br',
        'Cl',
        'F',
        'I',
        'Li',
        'Na',
        'K',
        'Rb',
        'Cs',
    }
)

PYROLYSIS_DEMETALLATING_SINGLE_ATOMS = frozenset(
    {'Zn', 'Cd', 'Hg', 'Ga', 'In', 'Tl', 'Pb', 'Bi'}
)


def phase_stable_at_application_T(
    genome: tuple, T_K: float = None, application: str = APPLICATION_PYROLYSIS
) -> dict:
    """
    Encoded-phase admissibility at the application's operating temperature.

    See ADR 0001. Applied at candidate selection only — not to the coverage
    certificate denominator.

    Args:
        genome: Input controlling genome.
        T_K: Input controlling T K.

    Returns:
        Validated dict output for this operation.
    """
    material_class = genome[0] if genome else None
    if application == APPLICATION_PEMFC:
        return pemfc_cathode_scope(genome if genome else (None,))
    if application != APPLICATION_PYROLYSIS:
        return {
            'status': 'unknown',
            'reason': f'unrecognized application {application!r}',
        }

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

    if T_K is not None:
        effective_T_K = float(T_K)
    elif (
        material_class in ('MoltenMetal', 'HEA')
        and len(genome) > 4
        and isinstance(genome[4], (int, float))
    ):
        effective_T_K = max(float(genome[4]), 1073.15)
    else:
        effective_T_K = 1073.15

    host_idx = 2 if material_class == 'SAA' else 1
    if (
        material_class in ('SolidCatalyst', 'SupportedNanoparticle', 'SAA')
        and len(genome) > host_idx
        and str(genome[host_idx]) in LOW_MELTING_HOST_METALS_K
        and max(effective_T_K, 1073.15) > LOW_MELTING_HOST_METALS_K[str(genome[host_idx])]
    ):
        host = str(genome[host_idx])
        tm_k = LOW_MELTING_HOST_METALS_K[host]
        return {
            'status': 'out_of_scope',
            'reason': (
                f'{material_class} host {host} melts at {tm_k} K < pyrolysis T '
                f'({max(effective_T_K, 1073.15):.1f} K; use MoltenMetal pathway)'
            ),
        }

    if material_class in ('SolidCatalyst', 'SupportedNanoparticle'):
        dop_idx = 5 if material_class == 'SolidCatalyst' else 6
        nsub_idx = 6 if material_class == 'SolidCatalyst' else 7
        n_sub = (
            int(genome[nsub_idx])
            if len(genome) > nsub_idx and isinstance(genome[nsub_idx], (int, float))
            else 0
        )
        dopants = (
            tuple(str(d) for d in genome[dop_idx])
            if len(genome) > dop_idx and isinstance(genome[dop_idx], (tuple, list))
            else ()
        )
        if n_sub > 0:
            volatile_dops = [d for d in dopants if d in PYROLYSIS_VOLATILE_ELEMENTS]
            if volatile_dops:
                return {
                    'status': 'out_of_scope',
                    'reason': (
                        f'{material_class} dopant(s) {volatile_dops} volatilize or '
                        'desorb as halide/metal vapor at pyrolysis temperatures'
                    ),
                }

    if (
        material_class == 'SAA'
        and len(genome) > 1
        and str(genome[1]) in PYROLYSIS_VOLATILE_ELEMENTS
    ):
        return {
            'status': 'out_of_scope',
            'reason': f'SAA trace element {genome[1]} is volatile at pyrolysis temperatures',
        }

    if (
        material_class == 'HEA'
        and len(genome) > 1
        and isinstance(genome[1], (tuple, list))
    ):
        comps = [str(c) for c in genome[1]]
        volatile_comps = [
            c
            for c in comps
            if c in PYROLYSIS_VOLATILE_ELEMENTS or c in LOW_MELTING_HOST_METALS_K
        ]
        if volatile_comps:
            return {
                'status': 'out_of_scope',
                'reason': (
                    f'HEA component(s) {volatile_comps} melt or boil below solid-solution '
                    'pyrolysis operating temperatures (e.g. Zn boils at 1180 K)'
                ),
            }

    if material_class == 'MoltenMetal' and len(genome) >= 3:
        host_m = str(genome[1])
        prom_m = str(genome[2])
        if host_m in PYROLYSIS_VOLATILE_ELEMENTS:
            return {
                'status': 'out_of_scope',
                'reason': (
                    f'MoltenMetal host {host_m} boils/volatilizes or forms toxic '
                    'hydride gas (H2Te/Zn vapor) at pyrolysis temperatures'
                ),
            }
        if prom_m not in ('None', 'none') and prom_m in PYROLYSIS_VOLATILE_ELEMENTS:
            return {
                'status': 'out_of_scope',
                'reason': f'MoltenMetal promoter {prom_m} is volatile at pyrolysis temperatures',
            }

    if material_class == 'SAC':
        if len(genome) > 1 and str(genome[1]) in PYROLYSIS_DEMETALLATING_SINGLE_ATOMS:
            return {
                'status': 'out_of_scope',
                'reason': (
                    f'SAC metal {genome[1]} demetallates and volatilizes from M-Nx '
                    'coordination sites under high-T reducing CH4/H2'
                ),
            }
        if len(genome) > 2 and 'Se' in str(genome[2]):
            return {
                'status': 'out_of_scope',
                'reason': f'SAC seleno-coordination {genome[2]} volatilizes as H2Se at pyrolysis temperatures',
            }
        if len(genome) > 3 and str(genome[3]) == 'Phosphorene':
            return {
                'status': 'out_of_scope',
                'reason': 'SAC Phosphorene support sublimes/decomposes above 650 K',
            }
        if len(genome) > 4 and str(genome[4]) in (
            PYROLYSIS_VOLATILE_ELEMENTS
            | {'CO', 'NO', 'CN', 'OH', 'OOH', 'SH', 'NH2', 'H2O', 'py'}
        ):
            return {
                'status': 'out_of_scope',
                'reason': f'SAC axial ligand {genome[4]} desorbs at pyrolysis temperatures',
            }

    if material_class == 'DAC':
        if len(genome) > 2:
            m1, m2 = str(genome[1]), str(genome[2])
            if (
                m1 in PYROLYSIS_DEMETALLATING_SINGLE_ATOMS
                or m2 in PYROLYSIS_DEMETALLATING_SINGLE_ATOMS
            ):
                return {
                    'status': 'out_of_scope',
                    'reason': (
                        f'DAC pair {m1}-{m2} contains low-melting/volatile metal that '
                        'demetallates at pyrolysis temperatures'
                    ),
                }
        if len(genome) > 4 and str(genome[4]) == 'Phosphorene':
            return {
                'status': 'out_of_scope',
                'reason': 'DAC Phosphorene support sublimes/decomposes above 650 K',
            }

    if (
        material_class == 'MetalFreeCarbon'
        and len(genome) > 5
        and str(genome[5]) in PYROLYSIS_VOLATILE_ELEMENTS
    ):
        return {
            'status': 'out_of_scope',
            'reason': f'MetalFreeCarbon halogen/volatile dopant {genome[5]} desorbs at pyrolysis temperatures',
        }

    return {'status': 'candidate', 'reason': None}


def parse_encoded_genome(raw):
    """Parse a screening-row genome (tuple or literal string) or return None.

    Args:
        raw: Input controlling raw.
    """
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
    """True only if the genome is admissible for turquoise pyrolysis reactors.

    Args:
        raw: Input controlling raw.
        T_K: Input controlling T K.

    Returns:
        Validated bool output for this operation.
    """
    genome = parse_encoded_genome(raw)
    if genome is None:
        return False
    return phase_stable_at_application_T(genome, T_K)['status'] == 'candidate'


def select_turquoise_pyrolysis_candidates(df, top_k=None, genome_col: str = 'genome'):
    """
    Drop phase-unstable rows then optionally take top_k by E_act.

    Candidate-selection filter only. Does not alter coverage certificates.

    Args:
        df: Input controlling df.
        genome_col: Input controlling genome col.

    Returns:
        A dataframe containing only phase-admissible candidate rows.
    """
    if df is None:
        raise ValueError('screening frame is required for turquoise pyrolysis scope')
    if genome_col not in df.columns:
        raise ValueError(
            f'{genome_col!r} is required to apply phase_stable_at_application_T'
        )
    scoped = df[df[genome_col].map(is_turquoise_pyrolysis_candidate)].copy()
    if top_k is None or len(scoped) == 0:
        return scoped
    if 'E_act' in scoped.columns:
        return scoped.nsmallest(int(top_k), 'E_act')
    return scoped.head(int(top_k))


def scope_pyrolysis_pool(
    df: FrameT, genome_col: str = 'genome'
) -> tuple[FrameT, AdmissibilitySummary]:
    """Admissibility pool for slate drawing, tolerant of genome-less frames.

    Returns ``(pool, note)``. Rows whose encoded phase is unstable at the
    pyrolysis T band are dropped (ADR 0001). Validity is *not* applied here:
    the validation route may still rescue invalid rows, and the caller
    applies ``valid`` for the reactor route. A frame without a genome column
    (synthetic tables, fixtures) is returned unchanged with the reason
    recorded so the omission is visible rather than silent.

    Args:
        df: Input controlling df.
        genome_col: Input controlling genome col.

    Returns:
        The type-preserving admissible pool and a named filter summary.
    """
    if df is None:
        raise ValueError('screening frame is required for turquoise pyrolysis scope')
    if genome_col not in df.columns:
        return df, {
            'filter': None,
            'reason': f'{genome_col} column absent; admissibility not applied',
        }
    pool = select_turquoise_pyrolysis_candidates(df, top_k=None, genome_col=genome_col)
    return pool, {
        'filter': 'phase_stable_at_application_T',
        'admissible_count': int(len(pool)),
        'dropped_count': int(len(df) - len(pool)),
    }


def scope_pemfc_pool(
    df: FrameT, genome_col: str = 'genome'
) -> tuple[FrameT, AdmissibilitySummary]:
    """Return candidates with a physically defined direct PEMFC cathode.

    Search coverage remains unchanged: every encoded material class is still
    traversed and certified.  This boundary only prevents classes such as
    molten metals and metal hydrides from being interpreted as solid PEMFC
    catalyst layers.  Genome-less fixture and legacy tables are retained with
    an explicit note so callers do not fail silently.

    Args:
        df: Candidate table entering PEMFC validation or stack modeling.
        genome_col: Column containing the encoded catalyst genome.

    Returns:
        The type-preserving admissible pool and a named filter summary.
    """
    if df is None:
        raise ValueError('screening frame is required for PEMFC cathode scope')
    if genome_col not in df.columns:
        return df, {
            'filter': None,
            'reason': f'{genome_col} column absent; admissibility not applied',
        }
    mask = df[genome_col].map(
        lambda raw: (
            pemfc_cathode_scope(genome)['status'] == 'candidate'
            if (genome := parse_encoded_genome(raw)) is not None
            else False
        )
    )
    pool = df[mask].copy()
    return pool, {
        'filter': 'pemfc_cathode_scope',
        'admissible_count': int(len(pool)),
        'dropped_count': int(len(df) - len(pool)),
    }


def is_validation_quota_class(material_class: str, application: str = None) -> bool:
    """Reserved eSen/DFT slots: dead classes stay enumerated, they do not get a quota.

    MetalHydride is exempt on every application. Pyrolysis and the PEMFC
    cathode (ORR) workflows additionally exempt every class that is out of
    their scope (ADR 0001), so no reserved slot goes to a class the
    application's objectives can never accept.

    Args:
        material_class: Input controlling material class.
        application: Workflow application name, or None for no scope filter.

    Returns:
        Validated bool output for this operation.
    """
    if material_class in VALIDATION_QUOTA_EXEMPT_CLASSES:
        return False
    if application in (APPLICATION_PYROLYSIS, 'pyrolysis', 'turquoise_hydrogen'):
        return material_class not in OUT_OF_SCOPE_TURQUOISE_PYROLYSIS_CLASSES
    if application in (APPLICATION_PEMFC, 'fuel_cell_orr', 'fuel_cell'):
        return material_class not in OUT_OF_SCOPE_PEMFC_CATHODE_CLASSES
    return True


def validation_quota_class_count(material_classes, application: str = None) -> int:
    """How many classes consume --min-validation-per-class reserved slots.

    Args:
        material_classes: Input controlling material classes.
        application: Input controlling application.

    Returns:
        Validated int output for this operation.
    """
    return sum(
        1 for cls in material_classes if is_validation_quota_class(cls, application)
    )
