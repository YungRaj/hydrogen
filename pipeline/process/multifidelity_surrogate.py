"""Portable, fail-closed surrogates for validated multiphysics artifacts.

These models approximate reactor-scale transport closures.  They do not replace
Cantera kinetics and never authorize candidate exclusion.  Predictions outside
the calibrated domain are explicitly referred to a full-physics calculation.
"""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
from typing import Iterable, Mapping

import numpy as np

from pipeline.process.multiphysics_contract import REQUIRED_OUTPUTS


MODEL_SCHEMA_VERSION = 1
TRANSPORT_CLOSURE_TARGETS = {
    'Fluidized': frozenset({
        'gas_velocity_m_s', 'u_mf_m_s', 'bubble_fraction'}),
    'MMBCR': frozenset({
        'gas_velocity_m_s', 'gas_holdup_fraction', 'bubble_diameter_mm'}),
    # Existing NTEC/electrochemical artifacts report overall performance, not
    # separable transport closures.  Enabling those modes requires explicit
    # closure fields and candidate-specific kinetic features first.
}


def _physical_closure(reactor: str, predictions: Mapping[str, float],
                      features: Mapping[str, float]) -> bool:
    """Check identities and bounds that hold independently of calibration."""
    if reactor == 'Fluidized':
        velocity = predictions.get('gas_velocity_m_s')
        minimum = predictions.get('u_mf_m_s')
        fraction = predictions.get('bubble_fraction')
        if velocity is not None and velocity <= 0:
            return False
        if minimum is not None and minimum <= 0:
            return False
        if velocity is not None and minimum is not None and velocity <= minimum:
            return False
        if fraction is not None and not 0 < fraction < 1:
            return False
    elif reactor == 'MMBCR':
        velocity = predictions.get('gas_velocity_m_s')
        holdup = predictions.get('gas_holdup_fraction')
        diameter_mm = predictions.get('bubble_diameter_mm')
        if velocity is not None and velocity <= 0:
            return False
        if holdup is not None and not 0 < holdup < 1:
            return False
        if diameter_mm is not None and diameter_mm <= 0:
            return False
        column_m = features.get('geometry.column_diameter_m')
        if (diameter_mm is not None and column_m is not None and
                diameter_mm / 1000 >= column_m):
            return False
    return True


@dataclass(frozen=True)
class PhysicsRecord:
    """One accepted full-physics result, expressed without solver dependencies."""

    case_id: str
    pathway_mode: str
    reactor_type: str
    features: Mapping[str, float]
    outputs: Mapping[str, float]


def record_from_artifact(artifact: Mapping, *, case_id: str | None = None) -> PhysicsRecord:
    """Extract a training row from an already validated artifact.

    Callers loading files should use ``load_validated_artifact`` first.  This
    function independently checks the fields relevant to model training so a
    malformed or legacy artifact cannot silently become a label.
    """
    if artifact.get('complete') is not True or artifact.get('convergence', {}).get(
            'converged') is not True:
        raise ValueError('surrogate labels require a complete, converged artifact')
    reactor = str(artifact.get('reactor_type', ''))
    snapshot = artifact.get('surrogate_inputs', {})
    if (snapshot.get('schema_version') != 1 or
            snapshot.get('reactor_type') != reactor or
            snapshot.get('units_in_field_names') is not True):
        raise ValueError('artifact lacks a compatible surrogate-input snapshot')
    features = _finite_mapping(snapshot.get('values'), 'features')
    outputs = _finite_mapping(artifact.get('outputs'), 'outputs')
    missing = REQUIRED_OUTPUTS.get(reactor, set()).difference(outputs)
    if missing:
        raise ValueError(f'artifact lacks required outputs: {sorted(missing)}')
    identity = case_id or str(artifact.get('provenance', {}).get('input_sha256', ''))
    if not identity:
        raise ValueError('a stable case identity is required')
    return PhysicsRecord(identity, str(artifact.get('pathway_mode', '')),
                         reactor, features, outputs)


def _finite_mapping(value, label: str) -> dict[str, float]:
    if not isinstance(value, Mapping) or not value:
        raise ValueError(f'{label} must be a non-empty mapping')
    result = {}
    for name, number in value.items():
        if (not isinstance(name, str) or not isinstance(number, (int, float)) or
                isinstance(number, bool) or not math.isfinite(float(number))):
            raise ValueError(f'{label} must contain finite numeric values')
        result[name] = float(number)
    return result


@dataclass(frozen=True)
class TransportSurrogate:
    """JSON-serializable ridge ensemble with calibration and domain limits."""

    pathway_mode: str
    reactor_type: str
    feature_names: tuple[str, ...]
    target_names: tuple[str, ...]
    center: np.ndarray
    scale: np.ndarray
    feature_min: np.ndarray
    feature_max: np.ndarray
    coefficients: np.ndarray
    validation_rmse: np.ndarray
    uncertainty_limits: np.ndarray
    training_case_ids: tuple[str, ...]
    validation_case_ids: tuple[str, ...]

    def __post_init__(self) -> None:
        feature_count, target_count = len(self.feature_names), len(self.target_names)
        if not self.pathway_mode or not self.reactor_type:
            raise ValueError('transport surrogate identity is required')
        if (feature_count < 1 or target_count < 1 or
                len(set(self.feature_names)) != feature_count or
                len(set(self.target_names)) != target_count):
            raise ValueError('model feature and target schemas must be nonempty and unique')
        unsupported = set(self.target_names).difference(
            TRANSPORT_CLOSURE_TARGETS.get(self.reactor_type, frozenset()))
        if unsupported:
            raise ValueError('serialized model contains unsupported transport targets')
        arrays = {
            'center': (self.center, (feature_count,)),
            'scale': (self.scale, (feature_count,)),
            'feature_min': (self.feature_min, (feature_count,)),
            'feature_max': (self.feature_max, (feature_count,)),
            'validation_rmse': (self.validation_rmse, (target_count,)),
            'uncertainty_limits': (self.uncertainty_limits, (target_count,)),
        }
        for name, (value, shape) in arrays.items():
            if value.shape != shape or not np.all(np.isfinite(value)):
                raise ValueError(f'model {name} has an invalid shape or value')
        design_size = 1 + 2 * feature_count
        if (self.coefficients.ndim != 3 or self.coefficients.shape[0] < 2 or
                self.coefficients.shape[1:] != (design_size, target_count) or
                not np.all(np.isfinite(self.coefficients))):
            raise ValueError('model coefficients have an invalid shape or value')
        if (np.any(self.scale <= 0) or np.any(self.feature_min > self.feature_max) or
                np.any(self.validation_rmse < 0) or
                np.any(self.uncertainty_limits < 0)):
            raise ValueError('model scales, domains, and errors must be physical')
        training, validation = set(self.training_case_ids), set(self.validation_case_ids)
        if (not training or not validation or
                len(training) != len(self.training_case_ids) or
                len(validation) != len(self.validation_case_ids) or
                training.intersection(validation)):
            raise ValueError('model case identities must be unique and split-disjoint')

    def sha256(self) -> str:
        """Return a deterministic identity for the complete serialized model."""
        payload = json.dumps(
            self.to_dict(), sort_keys=True, separators=(',', ':')).encode()
        return hashlib.sha256(payload).hexdigest()

    def predict(self, features: Mapping[str, float]) -> dict:
        try:
            actual = _finite_mapping(features, 'features')
            if set(actual) != set(self.feature_names):
                raise ValueError('feature schema does not match the trained model')
            raw = np.asarray([actual[name] for name in self.feature_names])
            z = (raw - self.center) / self.scale
            design = np.concatenate(([1.0], z, z * z))
            ensemble = np.einsum('edf,d->ef', self.coefficients, design)
            mean = ensemble.mean(axis=0)
            uncertainty = ensemble.std(axis=0)
            prediction_map = dict(zip(self.target_names, mean.tolist()))
            outside = bool(np.any(raw < self.feature_min) or
                           np.any(raw > self.feature_max))
            uncertain = bool(np.any(uncertainty > self.uncertainty_limits))
            finite = bool(np.all(np.isfinite(mean)) and np.all(np.isfinite(uncertainty)))
            physical = finite and _physical_closure(
                self.reactor_type, prediction_map, actual)
        except (TypeError, ValueError) as exc:
            return {'usable': False, 'decision': 'full_physics_required',
                    'reason': str(exc)}
        usable = finite and physical and not outside and not uncertain
        return {
            'usable': usable,
            'decision': 'surrogate_closure' if usable else 'full_physics_required',
            'reason': None if usable else (
                'outside_calibrated_domain' if outside else
                'ensemble_uncertainty_exceeded' if uncertain else
                'nonphysical_closure' if not physical else
                'nonfinite_prediction'),
            'predictions': prediction_map,
            'uncertainty_1sigma': dict(zip(self.target_names, uncertainty.tolist())),
            'validation_rmse': dict(zip(self.target_names,
                                        self.validation_rmse.tolist())),
            'candidate_exclusion_authorized': False,
        }

    def to_dict(self) -> dict:
        return {
            'schema_version': MODEL_SCHEMA_VERSION,
            'pathway_mode': self.pathway_mode,
            'reactor_type': self.reactor_type,
            'feature_names': list(self.feature_names),
            'target_names': list(self.target_names),
            'center': self.center.tolist(), 'scale': self.scale.tolist(),
            'feature_min': self.feature_min.tolist(),
            'feature_max': self.feature_max.tolist(),
            'coefficients': self.coefficients.tolist(),
            'validation_rmse': self.validation_rmse.tolist(),
            'uncertainty_limits': self.uncertainty_limits.tolist(),
            'training_case_ids': list(self.training_case_ids),
            'validation_case_ids': list(self.validation_case_ids),
            'scope': 'transport_closure_only',
            'candidate_exclusion_authorized': False,
        }

    def save(self, path: str | Path) -> Path:
        target = Path(path)
        target.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + '\n')
        return target

    @classmethod
    def load(cls, path: str | Path) -> 'TransportSurrogate':
        value = json.loads(Path(path).read_text())
        if (value.get('schema_version') != MODEL_SCHEMA_VERSION or
                value.get('scope') != 'transport_closure_only' or
                value.get('candidate_exclusion_authorized') is not False):
            raise ValueError('unsupported or unsafe transport-surrogate model')
        return cls(
            pathway_mode=value['pathway_mode'], reactor_type=value['reactor_type'],
            feature_names=tuple(value['feature_names']),
            target_names=tuple(value['target_names']),
            center=np.asarray(value['center'], dtype=float),
            scale=np.asarray(value['scale'], dtype=float),
            feature_min=np.asarray(value['feature_min'], dtype=float),
            feature_max=np.asarray(value['feature_max'], dtype=float),
            coefficients=np.asarray(value['coefficients'], dtype=float),
            validation_rmse=np.asarray(value['validation_rmse'], dtype=float),
            uncertainty_limits=np.asarray(value['uncertainty_limits'], dtype=float),
            training_case_ids=tuple(value['training_case_ids']),
            validation_case_ids=tuple(value['validation_case_ids']))


def fit_transport_surrogate(
        training: Iterable[PhysicsRecord], validation: Iterable[PhysicsRecord], *,
        targets: Iterable[str], validation_rmse_limits: Mapping[str, float],
        ensemble_size: int = 8, ridge: float = 1e-8,
        random_seed: int = 0) -> TransportSurrogate:
    """Fit one mode-local model and require a disjoint passing holdout."""
    train, holdout = list(training), list(validation)
    if not isinstance(ensemble_size, int) or ensemble_size < 2:
        raise ValueError('ensemble_size must be an integer of at least two')
    if not math.isfinite(float(ridge)) or ridge < 0:
        raise ValueError('ridge must be finite and nonnegative')
    if len(train) < 4 or len(holdout) < 2:
        raise ValueError('at least four training and two validation cases are required')
    if {row.case_id for row in train}.intersection(row.case_id for row in holdout):
        raise ValueError('training and validation case identities must be disjoint')
    identity = {(row.pathway_mode, row.reactor_type) for row in train + holdout}
    if len(identity) != 1:
        raise ValueError('a surrogate may cover exactly one pathway mode and reactor type')
    mode, reactor = next(iter(identity))
    features = tuple(sorted(train[0].features))
    target_names = tuple(targets)
    if not target_names:
        raise ValueError('at least one target is required')
    allowed_targets = TRANSPORT_CLOSURE_TARGETS.get(reactor, frozenset())
    unsupported = set(target_names).difference(allowed_targets)
    if unsupported:
        raise ValueError(
            f'outputs are not established transport closures for {reactor}: '
            f'{sorted(unsupported)}')
    if set(validation_rmse_limits) != set(target_names):
        raise ValueError('every target needs exactly one validation RMSE limit')
    for row in train + holdout:
        if set(row.features) != set(features):
            raise ValueError('all records must share an exact feature schema')
        if not set(target_names).issubset(row.outputs):
            raise ValueError('all records must contain every requested target')
        _finite_mapping(row.features, 'features')
        _finite_mapping({name: row.outputs[name] for name in target_names},
                        'outputs')
    x = np.asarray([[row.features[name] for name in features] for row in train])
    y = np.asarray([[row.outputs[name] for name in target_names] for row in train])
    center, scale = x.mean(axis=0), x.std(axis=0)
    scale[scale <= 1e-15] = 1.0
    z = (x - center) / scale
    design = np.column_stack((np.ones(len(x)), z, z * z))
    rng = np.random.default_rng(random_seed)
    coefficients = []
    penalty = np.eye(design.shape[1]) * float(ridge)
    penalty[0, 0] = 0.0
    for _ in range(ensemble_size):
        sample = rng.integers(0, len(train), len(train))
        xd, yd = design[sample], y[sample]
        coefficients.append(np.linalg.pinv(xd.T @ xd + penalty) @ xd.T @ yd)
    coefficients = np.asarray(coefficients)
    xv = np.asarray([[row.features[name] for name in features] for row in holdout])
    yv = np.asarray([[row.outputs[name] for name in target_names] for row in holdout])
    zv = (xv - center) / scale
    dv = np.column_stack((np.ones(len(xv)), zv, zv * zv))
    predicted = np.einsum('edf,vd->evf', coefficients, dv).mean(axis=0)
    rmse = np.sqrt(np.mean((predicted - yv) ** 2, axis=0))
    limits = np.asarray([float(validation_rmse_limits[name]) for name in target_names])
    if np.any(~np.isfinite(limits)) or np.any(limits < 0):
        raise ValueError('validation RMSE limits must be finite and nonnegative')
    if np.any(rmse > limits):
        failed = {name: float(error) for name, error, limit in
                  zip(target_names, rmse, limits) if error > limit}
        raise ValueError(f'holdout validation failed: {failed}')
    # A prediction must have no more ensemble spread than the independently
    # accepted holdout error.  Zero limits use a small numerical floor.
    uncertainty_limits = np.maximum(limits, 1e-12)
    return TransportSurrogate(
        mode, reactor, features, target_names, center, scale,
        x.min(axis=0), x.max(axis=0), coefficients, rmse,
        uncertainty_limits, tuple(row.case_id for row in train),
        tuple(row.case_id for row in holdout))
