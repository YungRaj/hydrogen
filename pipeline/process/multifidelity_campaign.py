"""Composable controller for one complete multi-fidelity campaign iteration."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from typing import Callable, Iterable, Mapping

from pipeline.process.campaign_ledger import CampaignLedger
from pipeline.process.full_physics_scheduler import (
    FullPhysicsRequest, schedule_full_physics_cases)
from pipeline.process.multiphysics_runner import run_backend
from pipeline.process.transport_model_registry import TransportModelRegistry
from pipeline.process.transport_training import CaseArtifactReference
from pipeline.process.transport_training import train_and_publish_transport_model


@dataclass(frozen=True)
class MultiFidelityCampaignServices:
    """All stateful or scientific operations required by the controller."""

    execute_case: Callable[[Mapping], CaseArtifactReference]
    train_and_publish: Callable[[list[CaseArtifactReference]], Mapping]
    load_published_model: Callable[[str, str], object]
    append_lineage: Callable[[str, Mapping], object]
    schedule_referrals: Callable[..., list[dict]] = schedule_full_physics_cases


def default_campaign_services(
        *, results_dir, registry_dir, ledger_path, targets,
        validation_rmse_limits, ensemble_size: int = 8, ridge: float = 1e-8,
        random_seed: int = 0, solver_runner: Callable = run_backend,
        trainer: Callable = train_and_publish_transport_model
        ) -> MultiFidelityCampaignServices:
    """Bind the controller to the production solver, validator, and registry.

    Args:
        results_dir: Directory containing or receiving calculation results.
        registry_dir: Directory used for registry dir.
        ledger_path: Filesystem location used for ledger path.
        targets: Output properties learned or evaluated by the model.
        validation_rmse_limits: Validation rmse limits used by this operation.
        ensemble_size: Number of ensemble size to use.
        ridge: Ridge used by this operation.
        random_seed: Seed controlling deterministic sampling or fitting.
        solver_runner: Injected callable used to perform solver runner.
        trainer: Injected callable used to perform trainer.

    Returns:
        Computed `MultiFidelityCampaignServices` result.
    """
    registry = TransportModelRegistry(registry_dir)
    ledger = CampaignLedger(ledger_path)

    def execute(plan: Mapping) -> CaseArtifactReference:
        required = {
            'candidate_id', 'pathway_mode', 'reactor_type', 'temperature_K',
            'case_dir', 'model_source', 'partition'}
        missing = required.difference(plan)
        if missing:
            raise ValueError(f'case plan lacks production inputs: {sorted(missing)}')
        solver_runner(
            mode=plan['pathway_mode'], reactor_type=plan['reactor_type'],
            candidate_id=plan['candidate_id'],
            temperature_K=float(plan['temperature_K']),
            case_dir=plan['case_dir'], results_dir=results_dir,
            model_source=plan['model_source'],
            fenics_model=plan.get('fenics_model'),
            timeout_s=int(plan.get('timeout_s', 86400)),
            max_coupling_iterations=int(
                plan.get('max_coupling_iterations', 20)))
        return CaseArtifactReference(
            candidate_id=str(plan['candidate_id']),
            pathway_mode=str(plan['pathway_mode']),
            reactor_type=str(plan['reactor_type']),
            temperature_K=float(plan['temperature_K']),
            partition=str(plan['partition']))

    def train(references: list[CaseArtifactReference]) -> Mapping:
        return trainer(
            references, results_dir=results_dir, model_registry=registry,
            targets=tuple(targets),
            validation_rmse_limits=dict(validation_rmse_limits),
            ensemble_size=ensemble_size, ridge=ridge,
            random_seed=random_seed)

    return MultiFidelityCampaignServices(
        execute_case=execute, train_and_publish=train,
        load_published_model=registry.load, append_lineage=ledger.append)


def _metric(query: Mapping, name: str, default: float = 0.0) -> float:
    value = query.get(name, default)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError(f'query {name} must be numeric')
    return float(value)


def _mapping_sha256(value: Mapping) -> str:
    try:
        encoded = json.dumps(value, sort_keys=True, separators=(',', ':')).encode()
    except (TypeError, ValueError) as exc:
        raise ValueError('campaign inputs must be JSON serializable') from exc
    return hashlib.sha256(encoded).hexdigest()


def run_multifidelity_iteration(
        *, designed_cases: Iterable[Mapping], screening_queries: Iterable[Mapping],
        services: MultiFidelityCampaignServices, referral_budget: int,
        minimum_per_region: int = 1,
        prior_references: Iterable[CaseArtifactReference] = ()) -> dict:
    """Execute, train, screen, refer, and record one fail-closed iteration.

    Args:
        designed_cases: Mapping supplying designed cases.
        screening_queries: Mapping supplying screening queries.
        services: Services used by this operation.
        referral_budget: Referral budget used by this operation.
        minimum_per_region: Minimum per region used by this operation.
        prior_references: Prior references used by this operation.

    Returns:
        Dictionary containing the computed values, status, and supporting metadata.
    """
    plans = [dict(case) for case in designed_cases]
    if not plans:
        raise ValueError('a campaign iteration requires designed cases')
    if not isinstance(referral_budget, int) or referral_budget < 1:
        raise ValueError('referral_budget must be a positive integer')
    queries = [dict(query) for query in screening_queries]
    services.append_lineage('multifidelity_iteration_started', {
        'designed_cases': [{
            'case_id': plan.get('case_id'), 'partition': plan.get('partition'),
            'input_sha256': _mapping_sha256(plan)} for plan in plans],
        'screening_queries': [{
            'query_id': query.get('query_id'),
            'input_sha256': _mapping_sha256(query)} for query in queries],
    })
    references = []
    for plan in plans:
        if plan.get('partition') not in {'training', 'validation'}:
            raise ValueError('case partitions must be preassigned before execution')
        reference = services.execute_case(plan)
        if not isinstance(reference, CaseArtifactReference):
            raise TypeError('case executor must return a CaseArtifactReference')
        if reference.partition != plan['partition']:
            raise ValueError('executed case changed its preassigned partition')
        references.append(reference)
    historical = list(prior_references)
    if any(not isinstance(value, CaseArtifactReference) for value in historical):
        raise TypeError('prior references must be CaseArtifactReference values')
    training = dict(services.train_and_publish([*historical, *references]))
    if training.get('status') != 'published':
        raise RuntimeError('campaign model was not published')
    mode = str(training.get('pathway_mode', ''))
    reactor = str(training.get('reactor_type', ''))
    model = services.load_published_model(mode, reactor)
    if model is None or model.sha256() != training.get('model_sha256'):
        raise RuntimeError('published model cannot be reloaded with matching identity')
    accepted, requests, decisions = [], [], []
    query_ids = set()
    for query in queries:
        query_id, region = query.get('query_id'), query.get('region')
        if (not isinstance(query_id, str) or not query_id or
                not isinstance(region, str) or not region):
            raise ValueError('every screening query needs query_id and region')
        if query_id in query_ids:
            raise ValueError('screening query identities must be unique')
        query_ids.add(query_id)
        prediction = model.predict(query.get('features', {}))
        decision = {
            'query_id': query_id, 'region': region,
            'input_sha256': _mapping_sha256(query), **prediction}
        decisions.append(decision)
        if prediction.get('usable') is True:
            accepted.append(decision)
        else:
            requests.append(FullPhysicsRequest(
                case_id=query_id, region=region,
                expected_improvement=_metric(query, 'expected_improvement'),
                uncertainty=_metric(query, 'uncertainty'),
                calibration_error=_metric(query, 'calibration_error'),
                disagreement=_metric(query, 'disagreement'),
                unproductive_history=_metric(query, 'unproductive_history')))
    if requests:
        budget = min(referral_budget, len(requests))
        referrals = services.schedule_referrals(
            requests, total_budget=budget,
            minimum_per_region=minimum_per_region)
    else:
        referrals = []
    summary = {
        'pathway_mode': mode,
        'reactor_type': reactor,
        'model_sha256': training['model_sha256'],
        'executed_case_ids': [reference.candidate_id for reference in references],
        'prior_case_count': len(historical),
        'training_case_ids': training.get('training_case_ids', []),
        'validation_case_ids': training.get('validation_case_ids', []),
        'screening_decisions': decisions,
        'accepted_surrogate_count': len(accepted),
        'full_physics_required_count': len(requests),
        'scheduled_referrals': referrals,
        'candidate_exclusion_authorized': False,
    }
    event = services.append_lineage('multifidelity_iteration_complete', summary)
    event_hash = event.get('event_sha256') if isinstance(event, Mapping) else None
    return {**summary, 'lineage_event_sha256': event_hash}
