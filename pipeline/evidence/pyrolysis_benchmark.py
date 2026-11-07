"""Fail-closed contracts for the external pyrolysis descriptor benchmark."""

from __future__ import annotations

import ast
import hashlib
import json
import math
from pathlib import Path
import random
from typing import Any

from pipeline.search.discovery import candidate_id
from pipeline.utils import K_B_EV_K

MIN_PER_LABEL = 12
LABELS = ('active', 'poor_or_poisoning')
REQUIRED_TEXT_FIELDS = (
    'condition_id',
    'catalyst_identity',
    'reactor_context',
    'source_id',
    'source_family',
    'label_basis',
    'source_url',
)
TRUE_VALUES = (True, 'True', 'true', 1, '1')
EV_PER_KJ_MOL = 1.0 / 96.4853321233  # 1 eV/molecule = F J/mol = 96.4853321233 kJ/mol
# NIST-JANAF integrated sensible enthalpy H°(298.15 K) - H°(0 K) = 8.468 kJ/mol for H2(g),
# converting Shomate H°(T) - H°(298.15 K) to a 0 K DFT reference state.
H2_H_298_MINUS_H_0_KJ_MOL = 8.468
# NIST WebBook Shomate polynomial parameters (A, B, C, D, E, F, G) for H2(g), t = T / 1000 K:
# H°(T) - H°(298.15) = A*t + B*t²/2 + C*t³/3 + D*t⁴/4 - E/t + F [kJ/mol]
# S°(T) = A*ln(t) + B*t + C*t²/2 + D*t³/3 - E/(2*t²) + G [J/(mol·K)]
H2_SHOMATE = (
    (298.0, 1000.0, (33.066178, -11.363417, 11.432816, -2.772874, -0.158558, -9.980797, 172.707974)),
    (1000.0, 2500.0, (18.563083, 12.257357, -2.859786, 0.268238, 1.977990, -1.147438, 156.288133)),
)


def effective_barrier_eV(barriers: list[float], temperature_K: float) -> float:
    """Return the barrier whose Arrhenius factor equals the mean site rate.

    Args:
        barriers: Non-negative replica barriers in eV.
        temperature_K: Absolute condition temperature.

    Returns:
        Numerically stable rate-averaged effective barrier in eV.
    """
    if not barriers or temperature_K <= 0.0:
        raise ValueError('rate averaging requires barriers and positive temperature')
    values = [float(value) for value in barriers]
    if not all(math.isfinite(value) and value >= 0.0 for value in values):
        raise ValueError('rate averaging requires finite non-negative barriers')
    # Log-sum-exp ensemble rate average: exp(-E_eff / k_B T) = (1/N) * sum_i exp(-E_i / k_B T),
    # shifted by E_min for numerical stability: E_eff = E_min - k_B*T*ln((1/N)*sum_i exp(-(E_i - E_min)/k_B T)).
    minimum = min(values)
    mean_relative_rate = sum(
        math.exp(-(value - minimum) / (K_B_EV_K * temperature_K))
        for value in values
    ) / len(values)
    return minimum - K_B_EV_K * temperature_K * math.log(mean_relative_rate)


def h2_chemical_potential_eV(
    temperature_K: float, pressure_bar: float = 1.0
) -> float:
    """Return H2 ideal-gas chemical potential relative to 0 K at 1 bar.

    Args:
        temperature_K: Absolute gas temperature within 298--2500 K.
        pressure_bar: Ideal-gas H2 partial pressure in bar.

    Returns:
        Thermal enthalpy, entropy, and pressure correction in eV per H2.
    """
    if pressure_bar <= 0.0:
        raise ValueError('H2 pressure must be positive')
    selected = next(
        (coefficients for low, high, coefficients in H2_SHOMATE
         if low <= temperature_K <= high),
        None,
    )
    if selected is None:
        raise ValueError('H2 Shomate correction requires 298 <= T <= 2500 K')
    a, b, c, d, e, f, g = selected
    t = temperature_K / 1000.0
    enthalpy_kj_mol = (
        a * t + b * t**2 / 2.0 + c * t**3 / 3.0 + d * t**4 / 4.0
        - e / t + f + H2_H_298_MINUS_H_0_KJ_MOL
    )
    entropy_j_mol_k = (
        a * math.log(t) + b * t + c * t**2 / 2.0 + d * t**3 / 3.0
        - e / (2.0 * t**2) + g
    )
    standard_eV = (
        enthalpy_kj_mol - temperature_K * entropy_j_mol_k / 1000.0
    ) * EV_PER_KJ_MOL
    return standard_eV + K_B_EV_K * temperature_K * math.log(pressure_bar)


def _parse_molten_replicas(
    row: dict[str, Any],
) -> tuple[float, list[tuple[float, float, float]]]:
    """Validate a molten alloy row and return (temperature_K, [(E_act_raw, dG_H2_rel, dG_C_rel), ...])."""
    try:
        genome = ast.literal_eval(str(row['genome']))
        temperature_K = float(row['temperature_K'])
        replicas = json.loads(str(row['molten_replicas_json']))
    except (KeyError, TypeError, ValueError, SyntaxError, json.JSONDecodeError) as exc:
        raise ValueError('missing molten rate/free-energy input') from exc
    if genome[0] != 'MoltenMetal' or genome[2] == 'None':
        raise ValueError('pure-host molten BEP barrier is outside validated scope')
    h2_mu = h2_chemical_potential_eV(temperature_K)
    parsed: list[tuple[float, float, float]] = []
    for replica in replicas:
        try:
            raw_act = float(replica['E_act_raw'])
            hydrogen_release = -2.0 * float(replica['dE_H']) + h2_mu
            carbon_release = -float(replica['dE_C'])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError('missing molten rate/free-energy input') from exc
        parsed.append((raw_act, hydrogen_release, carbon_release))
    return temperature_K, parsed


def rate_free_energy_span_eV(row: dict[str, Any]) -> float:
    """Return the rate-averaged, H2-free-energy-corrected molten score.

    Args:
        row: Completed ensemble screening row for a promoted molten alloy.

    Returns:
        Effective turnover-limiting free-energy span in eV.

    Raises:
        ValueError: If ensemble data are absent or the BEP scope is unsupported.
    """
    temperature_K, parsed_replicas = _parse_molten_replicas(row)
    spans: list[float] = []
    for raw_act, hydrogen_release, carbon_release in parsed_replicas:
        span = max(max(0.0, raw_act), hydrogen_release, carbon_release, 0.0)
        if not math.isfinite(span):
            raise ValueError('missing molten rate/free-energy input')
        spans.append(span)
    return effective_barrier_eV(spans, temperature_K)


def dft_anchor_offset_eV(
    reference_barriers_eV: list[float],
    predicted_barriers_eV: list[float],
    temperature_K: float,
) -> float:
    """Compute the single-anchor DFT minus MLIP/BEP rate-averaged barrier offset.

    Args:
        reference_barriers_eV: Independent DFT NEB transition-state barriers in eV.
        predicted_barriers_eV: Uncalibrated MLIP/BEP replica barriers in eV.
        temperature_K: Absolute condition temperature in kelvin.

    Returns:
        Additive barrier correction in eV restoring rate-consistent anchor agreement.
    """
    ref_eff = effective_barrier_eV(reference_barriers_eV, temperature_K)
    pred_eff = effective_barrier_eV(predicted_barriers_eV, temperature_K)
    return float(ref_eff - pred_eff)


def calibrated_rate_free_energy_span_eV(
    row: dict[str, Any],
    *,
    anchor_offset_eV: float = 0.0,
    anchor_source: str = 'uncalibrated',
) -> dict[str, Any]:
    """Return the DFT-anchor-calibrated molten rate/free-energy span and diagnostics.

    Args:
        row: Completed ensemble screening row for a promoted molten alloy.
        anchor_offset_eV: Additive DFT minus MLIP/BEP activation-barrier offset in eV.
        anchor_source: Provenance identifier for the DFT anchor.

    Returns:
        Dictionary with calibrated and uncalibrated spans, calibrated activation
        barrier, dominant limiting step, and anchor provenance.
    """
    if not math.isfinite(float(anchor_offset_eV)):
        raise ValueError('anchor_offset_eV must be finite')
    temperature_K, parsed_replicas = _parse_molten_replicas(row)
    if not str(anchor_source).strip():
        raise ValueError('anchor_source must be non-empty')
    uncalibrated_spans = rate_free_energy_span_eV(row)
    cal_spans: list[float] = []
    cal_activations: list[float] = []
    replica_steps: list[tuple[float, str]] = []
    for raw_act, hydrogen_release, carbon_release in parsed_replicas:
        cal_act = max(0.0, raw_act + float(anchor_offset_eV))
        terms = (
            ('methane_activation', cal_act),
            ('hydrogen_release', hydrogen_release),
            ('carbon_release', carbon_release),
        )
        step_name, step_val = max(terms, key=lambda item: item[1])
        span = max(step_val, 0.0)
        if not math.isfinite(span):
            raise ValueError('missing molten rate/free-energy input')
        cal_spans.append(span)
        cal_activations.append(cal_act)
        replica_steps.append((span, step_name))
    fastest_step = min(replica_steps, key=lambda item: item[0])[1]
    return {
        'rate_free_energy_span_eV': effective_barrier_eV(cal_spans, temperature_K),
        'uncalibrated_span_eV': float(uncalibrated_spans),
        'calibrated_activation_eV': effective_barrier_eV(
            cal_activations, temperature_K
        ),
        'anchor_offset_eV': float(anchor_offset_eV),
        'anchor_source': str(anchor_source),
        'dominant_limiting_step': fastest_step,
    }


def _row_h2_mu_eV(row: dict[str, Any]) -> float:
    """Return H2 chemical potential when condition temperature is present."""
    temp_raw = row.get('temperature_K')
    if temp_raw is None and 'genome' in row:
        try:
            genome = ast.literal_eval(str(row['genome']))
            if isinstance(genome, tuple) and len(genome) >= 5:
                temp_raw = genome[4]
        except (SyntaxError, ValueError, TypeError):
            temp_raw = None
    if temp_raw is None:
        return 0.0
    try:
        temperature_K = float(temp_raw)
    except (TypeError, ValueError):
        return 0.0
    if 298.0 <= temperature_K <= 2500.0:
        return h2_chemical_potential_eV(temperature_K)
    return 0.0


def energetic_span_eV(row: dict[str, Any]) -> float:
    """Return the largest uphill step in the minimal CH4 -> C + 2 H2 cycle.

    The score has no fitted coefficients.  It compares methane activation with
    the thermodynamic costs of releasing two adsorbed H atoms as H2 and
    forming or releasing one adsorbed C atom relative to the CH4 - 2 H2
    reference state at the condition temperature.  Lower values are better.

    Args:
        row: Completed pyrolysis screening result.

    Returns:
        Two-sided energetic-span score in eV.

    Raises:
        ValueError: If a required value is missing, non-finite, or censored.
    """
    names = ('E_act_raw', 'dE_H', 'dE_C')
    try:
        e_act, d_e_h, d_e_c = (float(row[name]) for name in names)
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('missing finite energetic-span input') from exc
    if not all(math.isfinite(value) for value in (e_act, d_e_h, d_e_c)):
        raise ValueError('missing finite energetic-span input')
    h2_mu = _row_h2_mu_eV(row)
    dg_c = d_e_c + 2.0 * h2_mu
    return max(e_act, -2.0 * d_e_h + h2_mu, dg_c, -dg_c, 0.0)


# Literature DFT transition-state anchors vs. eSEN-sm-conserving-all MLIP/BEP barriers:
# - Solid Ni(111) terrace CH4 -> CH3* + H* TS: 101 kJ/mol = 1.05 eV (Bengaard et al., J. Catal. 209, 2002)
#   vs. 0.386 eV eSEN climbing-image NEB (docs/evidence/molten_neb_calibration_lock.json; Delta E = +0.664 eV).
# - Molten Cu45Bi55 alloy CH4 TS at 1273 K: [2.53, 2.60] eV (Palmer et al., ACS Catal. 9, 2019; J. Phys. Chem. C 2020)
#   vs. 5-replica eSEN barriers on Cu45Bi55 at 1273.15 K (docs/evidence/pyrolysis_descriptor_benchmark.json).
SOLID_NI111_DFT_BARRIER_EV = 1.05
SOLID_NI111_ESEN_NEB_BARRIER_EV = 0.386
MOLTEN_CU45BI55_DFT_BARRIERS_EV = (2.53, 2.60)
MOLTEN_CU45BI55_ESEN_REPLICA_BARRIERS_EV = (
    1.6259,
    1.2933,
    1.4835,
    1.3954,
    1.3462,
)


def class_dft_barrier_offset_eV(
    material_class: str, temperature_K: float = 1273.15
) -> dict[str, Any]:
    """Return the literature DFT minus MLIP/BEP activation-barrier offset for a catalyst class.

    Uses the Palmer et al. (2019/2020) $\\text{Cu}_{45}\\text{Bi}_{55}$ ab initio MD NEB
    anchor ($[2.53, 2.60]\\text{ eV}$) for ``'MoltenMetal'`` and the Bengaard et al. (2002)
    $\\text{Ni}(111)$ terrace DFT NEB anchor ($1.05\\text{ eV}$ vs. $0.386\\text{ eV}$ eSEN
    CI-NEB) for solid catalyst classes.

    Args:
        material_class: Catalyst material class identifier (e.g., ``'MoltenMetal'`` or ``'SolidCatalyst'``).
        temperature_K: Absolute temperature in K used for Boltzmann rate averaging of replica anchors.

    Returns:
        Dictionary containing ``material_class``, ``anchor_offset_eV``, ``reference_barrier_eV``,
        ``predicted_barrier_eV``, and ``anchor_source``.
    """
    if not str(material_class).strip():
        raise ValueError('material_class must be non-empty')
    if temperature_K <= 0.0 or not math.isfinite(float(temperature_K)):
        raise ValueError('temperature_K must be positive and finite')
    cls_name = str(material_class).strip()
    if cls_name == 'MoltenMetal':
        ref_barriers = list(MOLTEN_CU45BI55_DFT_BARRIERS_EV)
        pred_barriers = list(MOLTEN_CU45BI55_ESEN_REPLICA_BARRIERS_EV)
        source = 'palmer_2019_2020_cu45bi55_aimd_neb'
    else:
        ref_barriers = [SOLID_NI111_DFT_BARRIER_EV]
        pred_barriers = [SOLID_NI111_ESEN_NEB_BARRIER_EV]
        source = 'bengaard_2002_ni111_terrace_neb'
    ref_eff = effective_barrier_eV(ref_barriers, float(temperature_K))
    pred_eff = effective_barrier_eV(pred_barriers, float(temperature_K))
    return {
        'material_class': cls_name,
        'anchor_offset_eV': float(ref_eff - pred_eff),
        'reference_barrier_eV': float(ref_eff),
        'predicted_barrier_eV': float(pred_eff),
        'anchor_source': source,
    }


def interfacial_carbon_displacement_A(row: dict[str, Any]) -> float:
    """Return the absolute vertical displacement of C* from the top melt/surface plane.

    In molten-metal and liquid-interface methane pyrolysis, active alloys pin
    nascent C* near the gas-liquid interface (``|z_C - z_surf| -> 0``) where
    carbon nucleates and floats off, whereas inactive p-block hosts either
    expel C* into the gas phase (``z_C - z_surf >> 0``, e.g. Sn, Pb) or sink
    C* into subsurface traps (``z_C - z_surf << 0``, e.g. Bi, single-phase
    liquid Ni33Sn67). Lower values are better and no parameters are fitted.

    Args:
        row: Completed pyrolysis screening result containing ``c_height_above_surface_A``.

    Returns:
        Non-negative interfacial carbon displacement ``|c_height_above_surface_A|`` in Å.

    Raises:
        ValueError: If ``c_height_above_surface_A`` is missing or non-finite.
    """
    try:
        c_height = float(row['c_height_above_surface_A'])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('missing finite interfacial carbon displacement input') from exc
    if not math.isfinite(c_height):
        raise ValueError('missing finite interfacial carbon displacement input')
    return abs(c_height)


def interfacial_pyrolysis_score_eV(row: dict[str, Any]) -> float:
    """Return the zero-parameter coupled activation and interfacial-carbon score.

    Combines the physical forward activation barrier ``max(E_act_raw, 0)`` (eV)
    with the unit-stiffness (1 eV/Å) interfacial carbon displacement penalty
    ``|c_height_above_surface_A|`` (Å). Lower values are better.

    Args:
        row: Completed pyrolysis screening result.

    Returns:
        Coupled activation + interfacial carbon retention score in eV-equivalent units.

    Raises:
        ValueError: If required inputs are missing or non-finite.
    """
    try:
        e_act = float(row['E_act_raw'])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError('missing finite interfacial pyrolysis score input') from exc
    if not math.isfinite(e_act):
        raise ValueError('missing finite interfacial pyrolysis score input')
    return max(e_act, 0.0) + interfacial_carbon_displacement_A(row)


def calibrated_interfacial_pyrolysis_score_eV(row: dict[str, Any]) -> dict[str, Any]:
    """Return the class-anchored DFT-calibrated activation + interfacial carbon score.

    Applies the literature DFT anchor offset $\\Delta E_{\\text{anchor}}$ for the
    row's ``material_class`` (Palmer $\\text{Cu}_{45}\\text{Bi}_{55}$ for ``'MoltenMetal'``,
    Bengaard $\\text{Ni}(111)$ for solid classes) to ``E_act_raw`` before adding the
    unit-stiffness ($1\\text{ eV/\\AA}$) interfacial carbon displacement penalty
    $|z_{\\text{C}} - z_{\\text{surf}}|$.

    Args:
        row: Completed pyrolysis screening result containing ``E_act_raw`` and
            ``c_height_above_surface_A``.

    Returns:
        Dictionary with ``calibrated_interfacial_score_eV``,
        ``uncalibrated_interfacial_score_eV``, ``calibrated_activation_eV``,
        ``interfacial_carbon_displacement_A``, ``anchor_offset_eV``, and ``anchor_source``.
    """
    uncalibrated = interfacial_pyrolysis_score_eV(row)
    displacement = interfacial_carbon_displacement_A(row)
    e_act_raw = float(row['E_act_raw'])
    mat_class = str(row.get('material_class') or '').strip()
    if not mat_class and 'genome' in row:
        try:
            parsed = ast.literal_eval(str(row['genome']))
            if isinstance(parsed, tuple) and parsed:
                mat_class = str(parsed[0])
        except (SyntaxError, ValueError, TypeError):
            mat_class = ''
    if not mat_class:
        mat_class = 'SolidCatalyst'
    temp_raw = row.get('temperature_K')
    temperature_K = 1273.15
    if temp_raw is not None:
        try:
            cand_t = float(temp_raw)
            if math.isfinite(cand_t) and cand_t > 0.0:
                temperature_K = cand_t
        except (TypeError, ValueError):
            pass
    anchor = class_dft_barrier_offset_eV(mat_class, temperature_K=temperature_K)
    cal_act = max(0.0, e_act_raw + float(anchor['anchor_offset_eV']))
    return {
        'calibrated_interfacial_score_eV': float(cal_act + displacement),
        'uncalibrated_interfacial_score_eV': float(uncalibrated),
        'calibrated_activation_eV': float(cal_act),
        'interfacial_carbon_displacement_A': float(displacement),
        'anchor_offset_eV': float(anchor['anchor_offset_eV']),
        'anchor_source': str(anchor['anchor_source']),
    }


def roc_auc(labels: list[bool], scores: list[float]) -> float:
    """Compute tie-aware ROC AUC where larger scores predict ``True``.

    Args:
        labels: Binary ground-truth values.
        scores: Aligned prediction scores.

    Returns:
        Area under the receiver-operating-characteristic curve.
    """
    if len(labels) != len(scores):
        raise ValueError('ROC AUC requires aligned positive and negative examples')
    positives = [score for label, score in zip(labels, scores) if label]
    negatives = [score for label, score in zip(labels, scores) if not label]
    if not positives or not negatives:
        raise ValueError('ROC AUC requires aligned positive and negative examples')
    wins = sum(
        float(positive > negative) + 0.5 * float(positive == negative)
        for positive in positives
        for negative in negatives
    )
    return wins / (len(positives) * len(negatives))


def grouped_bootstrap_auc(
    rows: list[dict[str, Any]],
    score_key: str,
    resamples: int = 10_000,
    seed: int = 20261010,
) -> dict[str, float]:
    """Estimate AUC and a percentile interval by resampling source families.

    Args:
        rows: Labelled score rows containing source family and label.
        score_key: Lower-is-better score to evaluate.
        resamples: Number of deterministic bootstrap resamples.
        seed: Reproducible pseudo-random seed.

    Returns:
        Point AUC and 95% percentile bounds.
    """
    groups: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault(str(row['source_family']), []).append(row)
    group_rows = list(groups.values())

    def auc(sample: list[dict[str, Any]]) -> float:
        labels = [row['label'] == 'active' for row in sample]
        scores = [-float(row[score_key]) for row in sample]
        return roc_auc(labels, scores)

    point = auc(rows)
    generator = random.Random(seed)
    estimates: list[float] = []
    while len(estimates) < resamples:
        sample = [
            row
            for _ in group_rows
            for row in generator.choice(group_rows)
        ]
        try:
            estimates.append(auc(sample))
        except ValueError:
            continue
    estimates.sort()
    return {
        'auc': point,
        'ci95_low': estimates[int(0.025 * resamples)],
        'ci95_high': estimates[int(0.975 * resamples) - 1],
    }


def dataset_digest(document: dict[str, Any]) -> str:
    """Return the content identity of a benchmark document without its digest.

    Args:
        document: Parsed benchmark dataset.

    Returns:
        Lowercase SHA-256 digest of canonical JSON content.
    """
    from pipeline.io.hashing import canonical_json_sha256

    body = {key: value for key, value in document.items() if key != 'dataset_sha256'}
    return canonical_json_sha256(body)


def validate_benchmark(document: dict[str, Any], require_locked: bool = True) -> list[str]:
    """Validate labels, provenance, independence, and exact genome mappings.

    Args:
        document: Parsed external benchmark dataset.
        require_locked: Whether readiness requires a matching immutable digest.

    Returns:
        Human-readable contract violations; an empty list means scoreable.
    """
    errors: list[str] = []
    if document.get('schema_version') != 1:
        errors.append('schema_version must be 1')
    records = document.get('records')
    if not isinstance(records, list):
        return errors + ['records must be a list']
    condition_ids: set[str] = set()
    sources: dict[str, set[str]] = {label: set() for label in LABELS}
    for index, row in enumerate(records):
        prefix = f'record {index}'
        if not isinstance(row, dict):
            errors.append(f'{prefix} must be an object')
            continue
        label = str(row.get('label', ''))
        if label not in sources:
            errors.append(f'{prefix} has invalid label')
        if row.get('exact_encoding') is not True:
            errors.append(f'{prefix} is not an exact encoding')
        for field in REQUIRED_TEXT_FIELDS:
            if not isinstance(row.get(field), str) or not row[field].strip():
                errors.append(f'{prefix} missing {field}')
        condition_id = row.get('condition_id')
        if isinstance(condition_id, str) and condition_id.strip():
            if condition_id in condition_ids:
                errors.append(f'{prefix} duplicates condition_id')
            condition_ids.add(condition_id)
        try:
            temperature_K = float(row['temperature_K'])
            if not math.isfinite(temperature_K) or temperature_K < 700.0:
                raise ValueError
        except (KeyError, TypeError, ValueError):
            errors.append(f'{prefix} temperature_K must be finite and >= 700 K')
        if not str(row.get('source_url', '')).startswith(('https://doi.org/', 'https://')):
            errors.append(f'{prefix} source_url is not stable HTTPS provenance')
        source_id = str(row.get('source_id', ''))
        if source_id.startswith('doi:'):
            expected_url = 'https://doi.org/' + source_id.removeprefix('doi:')
            if row.get('source_url') != expected_url:
                errors.append(f'{prefix} source_id and source_url disagree')
        try:
            genome = ast.literal_eval(str(row['genome']))
            candidate_id(genome)
        except (KeyError, SyntaxError, ValueError, TypeError):
            errors.append(f'{prefix} has invalid genome')
            continue
        if genome[0] == 'MoltenMetal':
            from pipeline.search.design_space import validate_molten_metal

            if not validate_molten_metal(genome):
                errors.append(f'{prefix} molten genome is outside the design space')
        if label in sources and isinstance(row.get('source_id'), str):
            sources[label].add(row['source_id'])
    for label, independent_sources in sources.items():
        if len(independent_sources) < MIN_PER_LABEL:
            errors.append(
                f'{label} has {len(independent_sources)} independent primary sources; '
                f'{MIN_PER_LABEL} required'
            )
    if require_locked:
        if document.get('locked') is not True:
            errors.append('dataset is not locked')
        if document.get('dataset_sha256') != dataset_digest(document):
            errors.append('dataset digest mismatch')
    return errors


def load_scoreable_benchmark(path: str | Path) -> dict[str, Any]:
    """Load a benchmark only after every preregistered curation gate passes.

    Args:
        path: JSON benchmark dataset to load.

    Returns:
        Parsed, checksum-verified benchmark document.

    Raises:
        ValueError: If any provenance, independence, encoding, or lock gate fails.
    """
    document = json.loads(Path(path).read_text())
    errors = validate_benchmark(document)
    if errors:
        raise ValueError('pyrolysis benchmark is not scoreable: ' + '; '.join(errors))
    return document


def verify_analysis_lock(path: str | Path, root: str | Path = '.') -> list[str]:
    """Return mismatches between an analysis lock and its frozen files.

    Args:
        path: JSON lock document.
        root: Repository root against which locked paths are resolved.

    Returns:
        Missing or content-mismatched locked-file errors.
    """
    lock = json.loads(Path(path).read_text())
    errors: list[str] = []
    for relative, expected in lock.get('sha256', {}).items():
        source = Path(root) / relative
        if not source.is_file():
            errors.append(f'missing locked file: {relative}')
            continue
        actual = hashlib.sha256(source.read_bytes()).hexdigest()
        if actual != expected:
            errors.append(f'locked file changed: {relative}')
    if not lock.get('sha256'):
        errors.append('lock contains no file hashes')
    return errors


def analyze_screening_results(
    benchmark: dict[str, Any], screening_rows: list[dict[str, Any]],
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    """Join locked conditions to screening results and run the frozen analysis.

    Args:
        benchmark: Checksum-verified benchmark document.
        screening_rows: Screening results keyed by their serialized genome.
        allow_incomplete: Exploratory only: score the resolvable conditions
            and list the rest instead of raising. Such a report never passes.

    Returns:
        Complete benchmark report including the preregistered pass verdict.

    Raises:
        ValueError: If any condition lacks an eligible screening result.
    """
    by_genome = {str(row.get('genome')): row for row in screening_rows}
    scored: list[dict[str, Any]] = []
    failures: list[str] = []
    for condition in benchmark['records']:
        result = by_genome.get(condition['genome'])
        if result is None:
            failures.append(f"{condition['condition_id']}: missing screening result")
            continue
        if result.get('valid') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: invalid screening result")
            continue
        if result.get('pyrolysis_viable') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: non-viable screening result")
            continue
        row_with_temp = {
            **result,
            'temperature_K': result.get('temperature_K', condition.get('temperature_K')),
        }
        try:
            span = energetic_span_eV(row_with_temp)
            current = float(result['E_act_raw'])
        except (KeyError, TypeError, ValueError) as exc:
            failures.append(f"{condition['condition_id']}: {exc}")
            continue
        h2_mu = _row_h2_mu_eV(row_with_temp)
        dg_c = float(result['dE_C']) + 2.0 * h2_mu
        scored.append({
            'condition_id': condition['condition_id'],
            'source_family': condition['source_family'],
            'label': condition['label'],
            'genome': condition['genome'],
            'E_act_eV': current,
            'energetic_span_eV': span,
            'limiting_step': max(
                (
                    ('methane_activation', current),
                    ('hydrogen_release', -2.0 * float(result['dE_H']) + h2_mu),
                    ('carbon_formation', dg_c),
                    ('carbon_release', -dg_c),
                ),
                key=lambda item: item[1],
            )[0],
        })
    if failures and not allow_incomplete:
        raise ValueError('benchmark screening is incomplete: ' + '; '.join(failures))
    current_stats = grouped_bootstrap_auc(scored, 'E_act_eV')
    span_stats = grouped_bootstrap_auc(scored, 'energetic_span_eV')
    improvement = span_stats['auc'] - current_stats['auc']
    passed = (
        span_stats['auc'] >= 0.70
        and span_stats['ci95_low'] > 0.50
        and improvement >= 0.10
    )
    return {
        'schema_version': 1,
        'dataset_sha256': benchmark['dataset_sha256'],
        'conditions_scored': len(scored),
        'current_E_act': current_stats,
        'energetic_span': span_stats,
        'auc_improvement': improvement,
        'acceptance_passed': passed and not failures,
        'unscored_conditions': failures,
        'records': scored,
    }


def analyze_interfacial_pyrolysis_benchmark(
    benchmark: dict[str, Any],
    screening_rows: list[dict[str, Any]],
    allow_incomplete: bool = False,
) -> dict[str, Any]:
    """Evaluate zero-parameter interfacial carbon retention and coupled pyrolysis descriptors.

    Args:
        benchmark: Checksum-verified benchmark document.
        screening_rows: Screening results keyed by their serialized genome.
        allow_incomplete: Exploratory only: score resolvable conditions instead of raising.

    Returns:
        Complete benchmark report comparing ``interfacial_carbon_displacement_A``
        and ``interfacial_pyrolysis_score_eV`` against ``E_act_raw``.
    """
    by_genome = {str(row.get('genome')): row for row in screening_rows}
    scored: list[dict[str, Any]] = []
    failures: list[str] = []
    for condition in benchmark['records']:
        result = by_genome.get(condition['genome'])
        if result is None:
            failures.append(f"{condition['condition_id']}: missing screening result")
            continue
        if result.get('valid') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: invalid screening result")
            continue
        if result.get('pyrolysis_viable') not in TRUE_VALUES:
            failures.append(f"{condition['condition_id']}: non-viable screening result")
            continue
        try:
            current = float(result['E_act_raw'])
            c_disp = interfacial_carbon_displacement_A(result)
            coupled = interfacial_pyrolysis_score_eV(result)
        except (KeyError, TypeError, ValueError) as exc:
            failures.append(f"{condition['condition_id']}: {exc}")
            continue
        scored.append({
            'condition_id': condition['condition_id'],
            'source_family': condition['source_family'],
            'label': condition['label'],
            'genome': condition['genome'],
            'E_act_eV': current,
            'interfacial_carbon_displacement_A': c_disp,
            'interfacial_pyrolysis_score_eV': coupled,
        })
    if failures and not allow_incomplete:
        raise ValueError('benchmark screening is incomplete: ' + '; '.join(failures))
    current_stats = grouped_bootstrap_auc(scored, 'E_act_eV')
    disp_stats = grouped_bootstrap_auc(scored, 'interfacial_carbon_displacement_A')
    coupled_stats = grouped_bootstrap_auc(scored, 'interfacial_pyrolysis_score_eV')
    improvement = coupled_stats['auc'] - current_stats['auc']
    passed = (
        coupled_stats['auc'] >= 0.70
        and coupled_stats['ci95_low'] > 0.50
        and improvement >= 0.10
    )
    return {
        'schema_version': 1,
        'dataset_sha256': benchmark['dataset_sha256'],
        'conditions_scored': len(scored),
        'current_E_act': current_stats,
        'interfacial_carbon_displacement': disp_stats,
        'interfacial_pyrolysis_score': coupled_stats,
        'auc_improvement': improvement,
        'acceptance_passed': passed and not failures,
        'unscored_conditions': failures,
        'records': scored,
    }


def main() -> None:
    """Screen or analyze the locked external benchmark."""
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('screen', 'analyze'))
    parser.add_argument(
        '--dataset', default='docs/evidence/pyrolysis_benchmark_dataset.json'
    )
    parser.add_argument(
        '--screening-csv', default='results/screening/pyrolysis_benchmark.csv'
    )
    parser.add_argument(
        '--report', default='results/evidence/pyrolysis_descriptor_benchmark.json'
    )
    parser.add_argument(
        '--lock', default='docs/evidence/pyrolysis_descriptor_benchmark_lock.json'
    )
    parser.add_argument(
        '--exploratory', action='store_true',
        help='allow code that differs from the lock (e.g. another screening '
             'model); the report is labelled exploratory and cannot pass')
    args = parser.parse_args()
    lock_errors = verify_analysis_lock(args.lock)
    if lock_errors and not args.exploratory:
        parser.error(
            '; '.join(lock_errors) + ' -- this analysis is frozen; re-run it from '
            'the commit recorded in docs/evidence/lock_provenance.json'
        )
    benchmark = load_scoreable_benchmark(args.dataset)
    if args.action == 'screen':
        from pipeline.screening.surface_screener import run_screening

        genomes = list(dict.fromkeys(
            ast.literal_eval(row['genome']) for row in benchmark['records']
        ))
        run_screening(genomes, Path(args.screening_csv).name)
        return

    import pandas as pd

    rows = pd.read_csv(args.screening_csv).to_dict(orient='records')
    report = analyze_screening_results(
        benchmark, rows, allow_incomplete=args.exploratory)
    report['screening_models'] = sorted(
        {str(row['screening_model']) for row in rows if row.get('screening_model')})
    if lock_errors:
        # Not the preregistered analysis: record why and never claim a pass.
        report['exploratory_modified_implementation'] = True
        report['lock_errors'] = lock_errors
        report['acceptance_passed'] = False
    target = Path(args.report)
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n')
    temporary.replace(target)
    print(json.dumps({key: value for key, value in report.items() if key != 'records'}, indent=2))


if __name__ == '__main__':
    main()
