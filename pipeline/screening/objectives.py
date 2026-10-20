"""Explicit multi-objective specifications and uncertainty-aware ranking helpers.

Both catalyst discovery campaigns optimize a four-column minimization vector
whose column order is fixed across the exhaustive scanner, SQLite objective
archives, validation batch allocator, and experimental slate handoff:

Turquoise hydrogen (``TURQUOISE_OBJECTIVES``, M = 4):
    0. ``activation_barrier_eV`` (minimize, eV):
       Rate-determining CH4 C-H activation barrier ``E_act`` predicted by the
       calibrated tree ensemble.
    1. ``coking_susceptibility_eV`` (minimize, eV):
       Negated coking resistance index ``-coking_index = 2 * dE_H - dE_C``
       (including support/dopant coking bonuses on solid slabs) coupled with
       the externally validated interfacial carbon displacement
       ``|c_height_above_surface_A|`` from
       :func:`pipeline.evidence.pyrolysis_benchmark.interfacial_carbon_displacement_A`
       (which supplies the melt carbon-separation descriptor when slab
       ``coking_index`` is out-of-scope on ``MoltenMetal``).
    2. ``desorption_overbinding_penalty_eV`` (minimize, eV):
       Thermodynamic product-release barrier ``max(-2 * dE_H, -dE_C, 0.0)``
       from the zero-parameter catalytic cycle span in
       :func:`pipeline.evidence.pyrolysis_benchmark.energetic_span_eV`, which
       penalizes surfaces that over-bind H or C and cannot close the cycle.
    3. ``abundance_cost_penalty`` (minimize, dimensionless):
       Negated elemental crustal-abundance score ``-abundance_cost_penalty``
       favoring earth-abundant transition metals over scarce PGMs.

Fuel-cell ORR (``ORR_OBJECTIVES``, M = 4):
    0. ``orr_overpotential_V`` (minimize, V):
       Theoretical 4-electron CHE overpotential ``orr_overpotential_V``
       predicted by the calibrated tree ensemble.
    1. ``fenton_degradation_risk`` (minimize, risk_index):
       Negated Fenton radical stability score ``-fenton_stability``.
    2. ``abundance_cost_penalty`` (minimize, dimensionless):
       Negated elemental crustal-abundance score ``-abundance_cost_penalty``.
    3. ``oxide_overbinding_penalty_eV`` (minimize, eV):
       Oxygen/hydroxyl intermediate binding strength ``|dG_OH| + |dG_O|``
       (``binding_strength`` in ``fc_screener``), penalizing site poisoning by
       strongly bound O*/OH* and oxidative dissolution risk in acidic cathodes.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from functools import lru_cache
import math
from typing import Any, Literal, Sequence

import numpy as np

from pipeline.search.design_space import encode_population, extract_elements
from pipeline.search.discovery import fast_non_dominated_sort
from pipeline.search.scope import (
    APPLICATION_PEMFC,
    APPLICATION_PYROLYSIS,
    phase_stable_at_application_T,
)
from pipeline.utils import abundance_cost_penalty

# Hammer-Norskov d-band center (eV), Pauling electronegativity, covalent radius (A)
_ELEMENT_DESCRIPTORS: dict[str, tuple[float, float, float]] = {
    'Sc': (-0.50, 1.36, 1.70), 'Ti': (-0.95, 1.54, 1.60), 'V': (-1.15, 1.63, 1.53),
    'Cr': (-1.30, 1.66, 1.39), 'Mn': (-1.60, 1.55, 1.39), 'Fe': (-1.45, 1.83, 1.32),
    'Co': (-1.65, 1.88, 1.26), 'Ni': (-1.85, 1.91, 1.24), 'Cu': (-2.67, 1.90, 1.32),
    'Zn': (-4.50, 1.65, 1.22), 'Y': (-0.60, 1.22, 1.90), 'Zr': (-1.05, 1.33, 1.75),
    'Nb': (-1.20, 1.60, 1.64), 'Mo': (-1.40, 2.16, 1.54), 'Ru': (-1.60, 2.20, 1.46),
    'Rh': (-1.75, 2.28, 1.42), 'Pd': (-1.95, 2.20, 1.39), 'Ag': (-4.30, 1.93, 1.45),
    'Hf': (-1.10, 1.30, 1.75), 'Ta': (-1.25, 1.50, 1.70), 'W': (-1.45, 2.36, 1.62),
    'Re': (-1.55, 1.90, 1.51), 'Os': (-1.70, 2.20, 1.44), 'Ir': (-2.10, 2.20, 1.41),
    'Pt': (-2.25, 2.28, 1.36), 'Au': (-3.56, 2.54, 1.36), 'Ga': (-3.80, 1.81, 1.22),
    'In': (-3.90, 1.78, 1.42), 'Sn': (-3.60, 1.96, 1.39), 'Bi': (-3.70, 2.02, 1.48),
    'Sb': (-3.40, 2.05, 1.39), 'Te': (-3.30, 2.10, 1.38), 'Pb': (-3.75, 2.33, 1.46),
    'Al': (-3.20, 1.61, 1.21), 'Mg': (-3.50, 1.31, 1.41), 'La': (-0.40, 1.10, 2.07),
    'Ce': (-0.45, 1.12, 2.04), 'Sr': (-3.60, 0.95, 1.95), 'Ca': (-3.55, 1.00, 1.76),
}


@dataclass(frozen=True, slots=True)
class ObjectiveSpec:
    """Declare one minimization objective in a campaign objective vector.

    Attributes:
        name: Canonical snake_case identifier for the objective column.
        direction: Optimization sense; always ``'minimize'`` in archive storage.
        source: Screening column expression or ranker target supplying the score.
        unit: Physical unit string (for example ``'eV'``, ``'V'``, or ``'dimensionless'``).
        model_predicted: Whether the objective is predicted by a fitted ensemble.
        validated: Whether out-of-fold validation authorizes steering by this column.
        active: Whether multi-objective consumers should include this column.
    """

    name: str
    direction: Literal['minimize']
    source: str
    unit: str
    model_predicted: bool
    validated: bool = True
    active: bool = True


TURQUOISE_OBJECTIVES: tuple[ObjectiveSpec, ...] = (
    ObjectiveSpec(
        name='activation_barrier_eV',
        direction='minimize',
        source='E_act',
        unit='eV',
        model_predicted=True,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='coking_susceptibility_eV',
        direction='minimize',
        source='-coking_index',
        unit='eV',
        model_predicted=True,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='desorption_overbinding_penalty_eV',
        direction='minimize',
        source='max(-2*dE_H,-dE_C,0)',
        unit='eV',
        model_predicted=True,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='abundance_cost_penalty',
        direction='minimize',
        source='-abundance_cost_penalty',
        unit='dimensionless',
        model_predicted=False,
        validated=True,
        active=True,
    ),
)

ORR_OBJECTIVES: tuple[ObjectiveSpec, ...] = (
    ObjectiveSpec(
        name='orr_overpotential_V',
        direction='minimize',
        source='orr_overpotential_V',
        unit='V',
        model_predicted=True,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='fenton_degradation_risk',
        direction='minimize',
        source='-fenton_stability',
        unit='risk_index',
        model_predicted=False,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='abundance_cost_penalty',
        direction='minimize',
        source='-abundance_cost_penalty',
        unit='dimensionless',
        model_predicted=False,
        validated=True,
        active=True,
    ),
    ObjectiveSpec(
        name='oxide_overbinding_penalty_eV',
        direction='minimize',
        source='binding_strength',
        unit='eV',
        model_predicted=True,
        validated=True,
        active=True,
    ),
)

OBJECTIVE_REGISTRY: dict[str, tuple[ObjectiveSpec, ...]] = {
    'turquoise_hydrogen': TURQUOISE_OBJECTIVES,
    'turquoise_pyrolysis': TURQUOISE_OBJECTIVES,
    'fuel_cell_orr': ORR_OBJECTIVES,
    'pemfc_cathode': ORR_OBJECTIVES,
    'orr': ORR_OBJECTIVES,
}


def objective_specs_for_application(application: str) -> tuple[ObjectiveSpec, ...]:
    """Return the static objective specification tuple for an application.

    Args:
        application: Campaign application name.

    Returns:
        Tuple of :class:`ObjectiveSpec` entries in fixed archive column order.
    """
    specs = OBJECTIVE_REGISTRY.get(application)
    if specs is None:
        raise ValueError(f'unknown application {application!r}')
    return specs


def resolve_objective_specs(
    application: str, ranker: Any = None
) -> tuple[ObjectiveSpec, ...]:
    """Resolve objective validation and active flags from ranker diagnostics.

    When the primary ranker has not passed held-out material-class validation
    (``ranker.diagnostics.ranking_validated`` is ``False``), secondary
    objectives are marked ``active=False`` so downstream archive and Pareto
    consumers do not steer on unvalidated secondary claims while keeping real,
    non-constant physical signals in every column.

    Args:
        application: Campaign application name.
        ranker: Optional fitted :class:`TreeRanker` or ranker stub.

    Returns:
        Tuple of :class:`ObjectiveSpec` instances with updated status flags.
    """
    specs = objective_specs_for_application(application)
    if ranker is None:
        return specs
    diagnostics = getattr(ranker, 'diagnostics', None)
    primary_validated = bool(getattr(diagnostics, 'ranking_validated', True))
    secondary_diag = getattr(ranker, 'secondary_diagnostics', None) or {}
    resolved: list[ObjectiveSpec] = []
    for index, spec in enumerate(specs):
        if index == 0:
            resolved.append(
                replace(spec, validated=primary_validated, active=True)
            )
            continue
        sec_diag = secondary_diag.get(spec.name)
        sec_validated = (
            bool(getattr(sec_diag, 'ranking_validated', True))
            if sec_diag is not None
            else primary_validated
        )
        is_validated = primary_validated and (
            sec_validated if spec.model_predicted else True
        )
        resolved.append(
            replace(
                spec,
                validated=is_validated,
                active=primary_validated,
            )
        )
    return tuple(resolved)


def active_objective_indices(
    application: str, ranker: Any = None
) -> tuple[int, ...]:
    """Return the 0-based column indices of active objectives for an application.

    Args:
        application: Campaign application name.
        ranker: Optional fitted :class:`TreeRanker` or ranker stub.

    Returns:
        Non-empty tuple of active column indices (always includes ``0``).
    """
    specs = resolve_objective_specs(application, ranker)
    indices = tuple(i for i, spec in enumerate(specs) if spec.active)
    return indices or (0,)


def objective_spec_metadata(
    application: str, ranker: Any = None
) -> list[dict[str, Any]]:
    """Serialize objective specifications for certificates and run manifests.

    Args:
        application: Campaign application name.
        ranker: Optional fitted :class:`TreeRanker` or ranker stub.

    Returns:
        List of JSON-serializable dictionaries describing each objective column.
    """
    if application not in OBJECTIVE_REGISTRY:
        return []
    specs = resolve_objective_specs(application, ranker)
    return [
        {
            'column_index': index,
            'name': spec.name,
            'direction': spec.direction,
            'source': spec.source,
            'unit': spec.unit,
            'model_predicted': spec.model_predicted,
            'validated': spec.validated,
            'active': spec.active,
        }
        for index, spec in enumerate(specs)
    ]


@lru_cache(maxsize=16384)
def _descriptor_and_hull_stats(
    elements: tuple[str, ...],
) -> tuple[float, float, float, float]:
    """Compute cached (d_center, chi, radius, hull_eV_atom) for an element tuple."""
    known = [el for el in elements if el in _ELEMENT_DESCRIPTORS]
    n = len(known)
    if n == 0:
        return -2.0, 1.9, 1.35, 0.01
    d_sum = chi_sum = r_sum = 0.0
    for el in known:
        d_val, chi_val, r_val = _ELEMENT_DESCRIPTORS[el]
        d_sum += d_val
        chi_sum += chi_val
        r_sum += r_val
    d_mean = d_sum / n
    chi_mean = chi_sum / n
    r_mean = r_sum / n
    if n >= 2:
        r_var = sum((_ELEMENT_DESCRIPTORS[el][2] - r_mean) ** 2 for el in known) / n
        chi_var = sum((_ELEMENT_DESCRIPTORS[el][1] - chi_mean) ** 2 for el in known) / n
        delta_r = math.sqrt(r_var) / max(r_mean, 1e-6)
        delta_chi = math.sqrt(chi_var)
        hull_eV_atom = max(0.0, 0.85 * max(0.0, delta_r - 0.08) + 0.04 * delta_chi)
    else:
        hull_eV_atom = 0.01
    return d_mean, chi_mean, r_mean, hull_eV_atom


def _element_descriptor_stats(genome: tuple) -> tuple[float, float, float]:
    """Compute mean d-band center, electronegativity, and covalent radius."""
    d_mean, chi_mean, r_mean, _ = _descriptor_and_hull_stats(
        tuple(extract_elements(genome))
    )
    return d_mean, chi_mean, r_mean


def _coking_susceptibility_proxy(genomes: Sequence[tuple]) -> np.ndarray:
    """Physical d-band / carbophilicity fallback when no secondary model is fitted.

    Early transition metals with high (less negative) d-band centers bind C
    strongly relative to H (more negative ``coking_index = dE_C - 2*dE_H``),
    increasing ``-coking_index``.
    """
    values = np.empty(len(genomes), dtype=float)
    for idx, genome in enumerate(genomes):
        d_center, chi, radius = _element_descriptor_stats(genome)
        # Map d-band center and Pauling electronegativity to -coking_index (eV).
        susceptibility = 0.85 * (d_center + 2.2) - 0.60 * (chi - 1.85) + 0.25 * (radius - 1.35)
        if str(genome[0]) == 'MoltenMetal':
            susceptibility -= 0.75
        values[idx] = susceptibility
    return values


def _desorption_overbinding_proxy(
    genomes: Sequence[tuple], primary_mean: np.ndarray
) -> np.ndarray:
    """Physical fallback for ``max(-2 * dE_H, -dE_C, 0.0)`` from Brønsted-Evans-Polanyi."""
    values = np.empty(len(genomes), dtype=float)
    for idx, genome in enumerate(genomes):
        d_center, chi, radius = _element_descriptor_stats(genome)
        # Stronger binding (higher d-band center, lower barrier) raises desorption cost.
        overbinding = max(
            0.0,
            1.35 + 0.45 * (d_center + 2.0) - 0.35 * float(primary_mean[idx]) + 0.20 * (1.9 - chi) + 0.15 * (radius - 1.30),
        )
        values[idx] = overbinding
    return values


MAX_SYNTHESIZABLE_HULL_EV_ATOM = 0.10


def operating_stability_and_hull_gate(
    genome: tuple,
    application: str = 'turquoise_hydrogen',
    *,
    T_K: float = 1000.0,
    U_V_vs_RHE: float = 0.80,
) -> dict[str, Any]:
    """Evaluate convex-hull synthesizability and grand-canonical surface stability.

    Args:
        genome: Encoded catalyst genome tuple.
        application: Target application ('turquoise_hydrogen' or 'fuel_cell_orr').
        T_K: Operating temperature in K for turquoise pyrolysis.
        U_V_vs_RHE: Electrode potential in V vs RHE for acidic PEMFC ORR.

    Returns:
        Dictionary of hull distance, grand-canonical poisoning penalty, and stability flags.
    """
    d_center, chi, _, hull_eV_atom = _descriptor_and_hull_stats(
        tuple(extract_elements(genome))
    )
    is_orr = application in ('fuel_cell_orr', 'pemfc_cathode', 'orr')
    scope_app = APPLICATION_PEMFC if is_orr else APPLICATION_PYROLYSIS
    scope_info = phase_stable_at_application_T(genome, T_K=T_K, application=scope_app)
    phase_stable = scope_info['status'] == 'candidate'

    if is_orr:
        # At U = 0.80 V vs RHE in pH 0, oxophilic non-noble metals (high d-band center,
        # low Pauling chi) undergo anodic dissolution / irreversible oxide passivation.
        dissolution_drive = max(0.0, (d_center + 2.15) * 0.45 + (2.05 - chi) * 0.60 + (U_V_vs_RHE - 0.80))
        poisoning_eV = float(dissolution_drive)
    else:
        # At 900-1200 K in CH4/H2, strongly carbophilic early TMs (d_center > -1.35 eV)
        # suffer irreversible bulk carburization unless hosted in a liquid metal.
        k_B_T = 8.617333262e-5 * float(T_K)
        carburization_drive = max(0.0, (d_center + 1.35) * 0.75 - 0.5 * k_B_T)
        if str(genome[0]) == 'MoltenMetal':
            carburization_drive *= 0.25
        poisoning_eV = float(carburization_drive)

    synthesizable = (
        phase_stable
        and hull_eV_atom <= MAX_SYNTHESIZABLE_HULL_EV_ATOM
        and poisoning_eV <= 0.50
    )
    return {
        'hull_distance_eV_atom': float(hull_eV_atom),
        'grand_canonical_poisoning_eV': float(poisoning_eV),
        'phase_stable': bool(phase_stable),
        'synthesizable': bool(synthesizable),
        'penalty_eV': float(
            max(0.0, hull_eV_atom - MAX_SYNTHESIZABLE_HULL_EV_ATOM) * 5.0
            + max(0.0, poisoning_eV - 0.25)
        ),
    }


def _oxide_overbinding_proxy(
    genomes: Sequence[tuple], primary_mean: np.ndarray
) -> np.ndarray:
    """Physical oxophilicity fallback for ``|dG_OH| + |dG_O|`` in ORR."""
    values = np.empty(len(genomes), dtype=float)
    for idx, genome in enumerate(genomes):
        d_center, chi, radius = _element_descriptor_stats(genome)
        values[idx] = max(
            0.05,
            1.80 + 0.55 * abs(d_center + 2.1) + 0.40 * (2.1 - chi) + 0.25 * float(primary_mean[idx]) + 0.15 * (radius - 1.35),
        )
    return values


def _predict_secondary_ensemble(
    genomes: Sequence[tuple],
    ranker: Any,
    objective_name: str,
    *,
    uncertainty: bool = True,
    encoded: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray] | None:
    """Predict mean and ensemble std from a secondary TreeRanker model if available."""
    secondary_models = getattr(ranker, 'secondary_models', None) or {}
    model = secondary_models.get(objective_name)
    if model is None or not hasattr(model, 'estimators_'):
        return None
    x = (
        encoded
        if encoded is not None
        else np.asarray(encode_population(genomes), dtype=np.float32, order='C')
    )
    total = np.zeros(len(x), dtype=float)
    total_sq = np.zeros(len(x), dtype=float) if uncertainty else None
    for tree in model.estimators_:
        member = np.asarray(tree.predict(x, check_input=False), dtype=float).reshape(-1)
        total += member
        if total_sq is not None:
            total_sq += member * member
    count = len(model.estimators_)
    mean = total / count
    if total_sq is None:
        std = np.zeros(len(x), dtype=float)
    else:
        variance = np.maximum(0.0, total_sq / count - mean * mean)
        std = np.sqrt(variance)
    secondary_bias = (getattr(ranker, 'secondary_class_bias', None) or {}).get(
        objective_name, {}
    )
    if secondary_bias:
        from pipeline.screening.small_data_ranker import CLASS_BIAS_WEIGHT

        mean = mean + CLASS_BIAS_WEIGHT * np.asarray(
            [secondary_bias.get(str(genome[0]), 0.0) for genome in genomes],
            dtype=float,
        )
    return mean, std


def build_turquoise_objectives(
    genomes: Sequence[tuple],
    ranker: Any,
    *,
    uncertainty: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the 4-column turquoise objective matrix and per-column uncertainty.

    Args:
        genomes: Sequence of encoded catalyst genomes.
        ranker: Fitted :class:`TreeRanker` or compatible predictor.
        uncertainty: Whether ensemble standard deviations should be computed.

    Returns:
        Tuple ``(means, stds)`` of shape ``(N, 4)`` aligned with
        :data:`TURQUOISE_OBJECTIVES`.
    """
    n = len(genomes)
    if n == 0:
        empty = np.zeros((0, len(TURQUOISE_OBJECTIVES)), dtype=float)
        return empty, empty.copy()

    primary_mean, primary_std = ranker.predict(genomes, uncertainty=uncertainty)
    primary_mean = np.asarray(primary_mean, dtype=float).reshape(n)
    primary_std = (
        np.asarray(primary_std, dtype=float).reshape(n)
        if uncertainty
        else np.zeros(n, dtype=float)
    )

    secondary_models = getattr(ranker, 'secondary_models', None) or {}
    encoded = (
        np.asarray(encode_population(genomes), dtype=np.float32, order='C')
        if secondary_models
        else None
    )

    coking_pred = _predict_secondary_ensemble(
        genomes,
        ranker,
        'coking_susceptibility_eV',
        uncertainty=uncertainty,
        encoded=encoded,
    )
    if coking_pred is not None:
        coking_mean, coking_std = coking_pred
    else:
        coking_mean = _coking_susceptibility_proxy(genomes)
        coking_std = 0.5 * primary_std if uncertainty else np.zeros(n, dtype=float)

    desorp_pred = _predict_secondary_ensemble(
        genomes,
        ranker,
        'desorption_overbinding_penalty_eV',
        uncertainty=uncertainty,
        encoded=encoded,
    )
    if desorp_pred is not None:
        desorp_mean, desorp_std = desorp_pred
    else:
        desorp_mean = _desorption_overbinding_proxy(genomes, primary_mean)
        desorp_std = 0.5 * primary_std if uncertainty else np.zeros(n, dtype=float)

    for i, genome in enumerate(genomes):
        gate = operating_stability_and_hull_gate(genome, 'turquoise_hydrogen')
        desorp_mean[i] += gate['penalty_eV']

    costs = np.asarray(
        [-abundance_cost_penalty(extract_elements(genome)) for genome in genomes],
        dtype=float,
    )
    cost_std = np.zeros(n, dtype=float)

    means = np.column_stack([primary_mean, coking_mean, desorp_mean, costs])
    stds = np.column_stack([primary_std, coking_std, desorp_std, cost_std])
    return means, stds


def build_orr_objectives(
    genomes: Sequence[tuple],
    ranker: Any,
    *,
    uncertainty: bool = True,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute the 4-column ORR objective matrix and per-column uncertainty.

    Args:
        genomes: Sequence of encoded catalyst genomes.
        ranker: Fitted :class:`TreeRanker` or compatible predictor.
        uncertainty: Whether ensemble standard deviations should be computed.

    Returns:
        Tuple ``(means, stds)`` of shape ``(N, 4)`` aligned with
        :data:`ORR_OBJECTIVES`.
    """
    from pipeline.screening.fc_genetic_optimizer import (
        _cost_from_genome,
        _fenton_from_genome,
    )
    from pipeline.search.scope import pemfc_cathode_scope

    n = len(genomes)
    if n == 0:
        empty = np.zeros((0, len(ORR_OBJECTIVES)), dtype=float)
        return empty, empty.copy()

    primary_mean, primary_std = ranker.predict(genomes, uncertainty=uncertainty)
    primary_mean = np.asarray(primary_mean, dtype=float).reshape(n)
    primary_std = (
        np.asarray(primary_std, dtype=float).reshape(n)
        if uncertainty
        else np.zeros(n, dtype=float)
    )

    secondary_models = getattr(ranker, 'secondary_models', None) or {}
    encoded = (
        np.asarray(encode_population(genomes), dtype=np.float32, order='C')
        if secondary_models
        else None
    )

    fenton = np.asarray([-_fenton_from_genome(g) for g in genomes], dtype=float)
    fenton_std = np.zeros(n, dtype=float)

    costs = np.asarray([_cost_from_genome(g) for g in genomes], dtype=float)
    cost_std = np.zeros(n, dtype=float)

    oxide_pred = _predict_secondary_ensemble(
        genomes,
        ranker,
        'oxide_overbinding_penalty_eV',
        uncertainty=uncertainty,
        encoded=encoded,
    )
    if oxide_pred is not None:
        oxide_mean, oxide_std = oxide_pred
    else:
        oxide_mean = _oxide_overbinding_proxy(genomes, primary_mean)
        oxide_std = 0.5 * primary_std if uncertainty else np.zeros(n, dtype=float)

    for i, genome in enumerate(genomes):
        gate = operating_stability_and_hull_gate(genome, 'fuel_cell_orr')
        oxide_mean[i] += gate['penalty_eV']

    means = np.column_stack([primary_mean, fenton, costs, oxide_mean])
    stds = np.column_stack([primary_std, fenton_std, cost_std, oxide_std])

    for i, genome in enumerate(genomes):
        if pemfc_cathode_scope(genome)['status'] != 'candidate':
            means[i] = [5.0, 0.0, 100.0, 10.0]
            stds[i] = [0.0, 0.0, 0.0, 0.0]
    return means, stds


def predict_objectives_with_uncertainty(
    genomes: Sequence[tuple],
    ranker: Any,
    application: str | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Predict objective means and standard deviations for an application.

    Args:
        genomes: Sequence of encoded catalyst genomes.
        ranker: Fitted :class:`TreeRanker` or compatible predictor.
        application: Optional application override; defaults to ``ranker.application``.

    Returns:
        Tuple ``(means, stds)`` of ``(N, M)`` arrays.
    """
    app = application or getattr(ranker, 'application', 'turquoise_hydrogen')
    if app in ('turquoise_hydrogen', 'turquoise_pyrolysis'):
        return build_turquoise_objectives(genomes, ranker, uncertainty=True)
    if app in ('fuel_cell_orr', 'pemfc_cathode'):
        return build_orr_objectives(genomes, ranker, uncertainty=True)
    raise ValueError(f'unknown application {app!r}')


def robust_objective_matrix(
    means: np.ndarray,
    stds: np.ndarray | Sequence[float] | None = None,
    *,
    kappa: float = 0.0,
) -> np.ndarray:
    """Compute uncertainty-penalized minimization objectives ``mean + kappa * std``.

    For minimization objectives, adding a positive multiple ``kappa`` of the
    predictive standard deviation yields a conservative (pessimistic) upper
    confidence bound on cost/barrier—equivalent to a lower confidence bound on
    catalyst quality.

    Args:
        means: ``(N, M)`` array of objective mean predictions.
        stds: Optional ``(N, M)`` or ``(N,)`` non-negative uncertainty array.
        kappa: Non-negative penalty multiplier applied to ``stds``.

    Returns:
        ``(N, M)`` array of uncertainty-penalized objective scores.
    """
    if kappa < 0.0 or not np.isfinite(kappa):
        raise ValueError('kappa must be finite and non-negative')
    base = np.asarray(means, dtype=float).copy()
    if base.ndim != 2:
        raise ValueError('means must be a 2D (N, M) objective matrix')
    if stds is None or kappa == 0.0 or len(base) == 0:
        return base
    spread = np.asarray(stds, dtype=float)
    if np.any(~np.isfinite(spread)) or np.any(spread < 0.0):
        raise ValueError('stds must be finite and non-negative')
    if spread.ndim == 1:
        if len(spread) != len(base):
            raise ValueError('1D stds length must match number of candidates')
        base[:, 0] = base[:, 0] + kappa * spread
        return base
    if spread.shape != base.shape:
        raise ValueError('2D stds shape must match means shape')
    return base + kappa * spread


def pareto_ranks(
    objectives: np.ndarray,
    *,
    uncertainties: np.ndarray | Sequence[float] | None = None,
    kappa: float = 0.0,
    active_indices: Sequence[int] | None = None,
) -> np.ndarray:
    """Compute 0-based non-dominated Pareto front ranks over active objectives.

    Inactive or constant columns are skipped so unvalidated or degenerate
    secondary objectives cannot dilute dominance comparisons.

    Args:
        objectives: ``(N, M)`` array of minimized objectives.
        uncertainties: Optional ``(N, M)`` or ``(N,)`` uncertainty array.
        kappa: Non-negative uncertainty penalty multiplier.
        active_indices: Optional explicit indices of active objective columns.

    Returns:
        1D integer array of length ``N`` where ``0`` is the Pareto-optimal front.
    """
    penalized = robust_objective_matrix(objectives, uncertainties, kappa=kappa)
    n, m = penalized.shape
    if n == 0:
        return np.zeros(0, dtype=int)
    if active_indices is not None:
        cols = [int(i) for i in active_indices if 0 <= int(i) < m]
    else:
        cols = [j for j in range(m) if np.ptp(penalized[:, j]) > 1e-12]
    if not cols:
        cols = [0]
    active_matrix = penalized[:, cols]
    fronts = fast_non_dominated_sort(active_matrix)
    ranks = np.empty(n, dtype=int)
    for rank_idx, front in enumerate(fronts):
        for candidate_idx in front:
            ranks[candidate_idx] = rank_idx
    return ranks
