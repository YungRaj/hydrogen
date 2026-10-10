#!/usr/bin/env python3
"""One append-only prospective campaign for catalyst-guided search.

Historical v2-v8 files are evidence inputs only; batch IDs never select code.
"""
from __future__ import annotations

import argparse, ast, hashlib, json, os, subprocess, tempfile
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

from pipeline.search.branch_search import _probe_indices
from pipeline.search.design_space import encode_population
from pipeline.search.discovery import candidate_id
from pipeline.search.indexed_space import CLASS_ORDER, CLASS_SIZES, candidate_at_class
from pipeline.search.scope import phase_stable_at_application_T, pemfc_cathode_scope
from pipeline.screening.protocols import ORR_PROTOCOL, PYROLYSIS_PROTOCOL

REPO_ROOT = Path(__file__).resolve().parent
ROOT = Path('results/prospective_search')
HISTORY = Path('docs/evidence/legacy_pilot_rounds.jsonl')
POLICY_LOCK = Path('docs/evidence/prospective_policy_lock.json')
APPS = (('turquoise_hydrogen', 'E_act', 'pyrolysis'),
        ('fuel_cell_orr', 'orr_overpotential_V', 'orr'))
PROTOCOLS = {
    'turquoise_hydrogen': PYROLYSIS_PROTOCOL.protocol_id,
    'fuel_cell_orr': ORR_PROTOCOL.protocol_id,
}
LEGACY_PROTOCOLS = {
    'turquoise_hydrogen': 'esen-sm-conserving-all-oc25:relax-v3:pyrolysis-v2',
    'fuel_cell_orr': 'esen-sm-conserving-all-oc25:relax-v3:orr-che-v2',
}
BASELINE_OUTCOMES = {
    'turquoise_hydrogen': Path('results/screening/pilot/divide_conquer_v8_pyrolysis.csv'),
    'fuel_cell_orr': Path('results/fuel_cell/pilot/divide_conquer_v8_orr.csv'),
}
HIT_FRACTION = 0.20
POLICY_SOURCE_PATHS = (
    Path('run_divide_conquer_pilot.py'),
    Path('pipeline/screening/small_data_ranker.py'),
    Path('pipeline/search/scope.py'),
    Path('pipeline/screening/protocols.py'),
    Path('pipeline/screening/surface_screener.py'),
    Path('pipeline/screening/fc_screener.py'),
)


def _coverage_select(
    scores: np.ndarray, genomes: list[tuple], budget: int
) -> np.ndarray:
    """Apply one class-floor slot before filling the budget by score."""
    by_class: dict[str, list[int]] = {}
    for index, genome in enumerate(genomes):
        by_class.setdefault(str(genome[0]), []).append(index)
    selected = [
        max(
            indices,
            key=lambda index: (
                float(scores[index]),
                candidate_id(genomes[index]),
            ),
        )
        for _, indices in sorted(by_class.items())
    ]
    selected.sort(
        key=lambda index: (-float(scores[index]), candidate_id(genomes[index]))
    )
    chosen = set(selected)
    remainder = sorted(
        (index for index in range(len(genomes)) if index not in chosen),
        key=lambda index: (-float(scores[index]), candidate_id(genomes[index])),
    )
    return np.asarray((selected + remainder)[:budget], dtype=int)


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _policy_source_hash() -> str:
    """Return a checkout-independent digest of every frozen source file."""
    digest = hashlib.sha256()
    for relative_path in POLICY_SOURCE_PATHS:
        normalized_path = relative_path.as_posix()
        content = (REPO_ROOT / relative_path).read_bytes().replace(b'\r\n', b'\n')
        digest.update(normalized_path.encode())
        digest.update(b'\0')
        digest.update(content)
        digest.update(b'\0')
    return digest.hexdigest()


def _require_locked_policy() -> str:
    """Fail unless current policy sources match the preregistered digest."""
    if not POLICY_LOCK.is_file():
        raise RuntimeError(f'missing prospective policy lock: {POLICY_LOCK}')
    lock = json.loads(POLICY_LOCK.read_text())
    expected = str(lock.get('policy_source_sha256', ''))
    actual = _policy_source_hash()
    if not expected or actual != expected:
        raise RuntimeError(
            f'prospective policy hash mismatch: expected {expected}, got {actual}'
        )
    return actual


def _atomic_json(path: Path, payload: dict) -> None:
    """Replace one JSON artifact atomically after flushing it to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f'.{path.name}.', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w') as stream:
            json.dump(payload, stream, indent=2)
            stream.write('\n')
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _atomic_text(path: Path, content: str) -> None:
    """Replace one text artifact atomically after flushing it to disk."""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(prefix=f'.{path.name}.', dir=path.parent)
    try:
        with os.fdopen(descriptor, 'w') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    except BaseException:
        Path(temporary).unlink(missing_ok=True)
        raise


def _root(batch: str) -> Path:
    if not batch or not all(c.isalnum() or c in '-_' for c in batch):
        raise ValueError('invalid batch ID')
    return ROOT / 'batches' / batch


def _outcome(batch: str, app: str) -> Path:
    locations = {
        'turquoise_hydrogen': ('screening', 'pyrolysis'),
        'fuel_cell_orr': ('fuel_cell', 'orr'),
    }
    if app not in locations:
        raise ValueError(f'unknown application {app}')
    folder, suffix = locations[app]
    return Path('results') / folder / 'prospective' / f'{batch}_{suffix}.csv'


def _batch_ledger_rows() -> list[dict]:
    """Return batch records in append-only campaign chronology."""
    ledger = ROOT / 'batches.jsonl'
    if not ledger.exists():
        return []
    rows = [json.loads(line) for line in ledger.read_text().splitlines() if line.strip()]
    batch_ids = [str(row['batch_id']) for row in rows]
    if len(batch_ids) != len(set(batch_ids)):
        raise RuntimeError('batch ledger contains duplicate batch IDs')
    return rows


def _ordered_finalized_batches() -> list[tuple[str, Path]]:
    """Return finalized analyses in append-only batch-ledger order."""
    finalized = []
    for row in _batch_ledger_rows():
        batch = str(row['batch_id'])
        manifest = _root(batch) / 'manifest.json'
        if not manifest.is_file() or _hash(manifest) != row['manifest_sha256']:
            raise RuntimeError(f'batch ledger manifest mismatch: {manifest}')
        analysis = _root(batch) / 'analysis.json'
        if analysis.is_file():
            finalized.append((batch, analysis))
    return finalized


def _completed_outcomes(app: str, before_batch: str | None = None) -> list[Path]:
    """Return checksum-verified outcomes from finalized prospective batches."""
    paths = []
    for batch, analysis in _ordered_finalized_batches():
        report = json.loads(analysis.read_text())
        if str(report['batch_id']) != batch:
            raise RuntimeError(f'analysis batch ID mismatch: {analysis}')
        if batch == before_batch:
            break
        path = _outcome(batch, app)
        expected = report.get('outcome_sha256', {}).get(app)
        if expected is None:
            continue
        if not path.is_file() or _hash(path) != expected:
            raise RuntimeError(f'finalized outcome checksum mismatch: {path}')
        paths.append(path)
    return paths


def _training(
    app: str, valid_only: bool = True, before_batch: str | None = None
) -> pd.DataFrame:
    """Load only tracked, checksum-verifiable outcomes from one protocol."""
    if app not in PROTOCOLS:
        raise ValueError(f'unknown application {app}')
    baseline = BASELINE_OUTCOMES[app]
    if not baseline.is_file():
        raise RuntimeError(
            f'missing tracked {app} training baseline: {baseline}; restore the file from git'
        )
    paths = [baseline, *_completed_outcomes(app, before_batch)]
    frames = []
    for path in paths:
        frame = pd.read_csv(path)
        if 'screening_protocol' not in frame.columns:
            raise RuntimeError(f'training data has no screening_protocol: {path}')
        allowed = {PROTOCOLS[app], LEGACY_PROTOCOLS[app]}
        mismatched = ~frame.screening_protocol.isin(allowed)
        if bool(mismatched.any()):
            protocols = sorted(
                frame.loc[mismatched, 'screening_protocol'].astype(str).unique()
            )
            raise RuntimeError(f'incompatible screening protocol in {path}: {protocols}')
        legacy = frame.screening_protocol.eq(LEGACY_PROTOCOLS[app])
        if bool(legacy.any()):
            # Relax-v4 changes SAC/DAC geometry and adsorption placement. Other
            # material classes are structurally identical and remain compatible.
            changed_geometry = frame.material_class.isin(('SAC', 'DAC'))
            frame = frame.loc[~(legacy & changed_geometry)].copy()
            frame.loc[:, 'screening_protocol'] = PROTOCOLS[app]
        frames.append(frame)
    frame = pd.concat(frames, ignore_index=True, sort=False)
    if valid_only:
        from pipeline.screening.small_data_ranker import training_metric_eligibility
        frame = frame.loc[training_metric_eligibility(frame, app)].copy()
    frame['_id'] = [candidate_id(ast.literal_eval(x)) for x in frame.genome]
    return frame.drop_duplicates('_id', keep='last').drop(columns='_id')


def _select_policies(
    pool: list[tuple],
    app: str,
    budget: int,
    train: pd.DataFrame,
    all_train: pd.DataFrame,
    random_seed: int,
) -> dict:
    """Fit the frozen policy and select equal-budget candidates from one pool."""
    from dataclasses import asdict
    from sklearn.ensemble import ExtraTreesClassifier
    from pipeline.screening.small_data_ranker import (
        CLASS_SUCCESS_WEIGHT,
        class_success_probability,
        fit_tree_ranker,
        orr_catalyst_acquisition,
        orr_tree_objectives,
        turquoise_catalyst_acquisition,
        turquoise_tree_objectives,
    )

    eligible = [
        index for index, genome in enumerate(pool)
        if (
            phase_stable_at_application_T(genome)['status'] == 'candidate'
            if app == 'turquoise_hydrogen'
            else pemfc_cathode_scope(genome)['status'] == 'candidate'
        )
    ]
    ranker = fit_tree_ranker(train, app)
    objectives = (
        turquoise_tree_objectives(pool, ranker)
        if app == 'turquoise_hydrogen'
        else orr_tree_objectives(pool, ranker)
    )
    _, uncertainty = ranker.predict(pool)
    genomes = [ast.literal_eval(raw) for raw in all_train.genome]
    validity = ExtraTreesClassifier(
        n_estimators=256,
        min_samples_leaf=2,
        max_features=1.0,
        class_weight='balanced',
        random_state=20260722,
        n_jobs=-1,
    )
    validity.fit(encode_population(genomes), all_train.valid.eq(True).to_numpy(int))
    validity_score = validity.predict_proba(encode_population(pool))[
        :, list(validity.classes_).index(1)
    ]
    candidates = [pool[index] for index in eligible]
    selected_budget = min(len(candidates), budget)
    catalyst_score = (
        turquoise_catalyst_acquisition(
            objectives[eligible, 0], validity_score[eligible]
        )
        if app == 'turquoise_hydrogen'
        else orr_catalyst_acquisition(
            objectives[eligible, 0], uncertainty[eligible], validity_score[eligible]
        )
    )
    catalyst_score += CLASS_SUCCESS_WEIGHT * class_success_probability(
        all_train, candidates, app
    )
    selections = {
        'catalyst': _coverage_select(catalyst_score, candidates, selected_budget),
        'uncertainty': _coverage_select(
            uncertainty[eligible], candidates, selected_budget
        ),
        'validity': _coverage_select(
            validity_score[eligible], candidates, selected_budget
        ),
    }
    ids = [candidate_id(genome) for genome in pool]
    return {
        'application': app,
        'training_rows': len(train),
        'training_digest': hashlib.sha256(
            '\n'.join(
                sorted(candidate_id(ast.literal_eval(raw)) for raw in train.genome)
            ).encode()
        ).hexdigest(),
        'eligible_ids': [ids[index] for index in eligible],
        'budget': selected_budget,
        'policy_selected_ids': {
            name: [ids[eligible[index]] for index in selected]
            for name, selected in selections.items()
        },
        'ranker_diagnostics': asdict(ranker.diagnostics),
        'random_seed': random_seed,
        'random_trials': 50000,
        'random_one_sided_alpha': 0.05,
    }


def _prior() -> set[str]:
    prior = set()
    paths = list(Path('results/pilot').glob('divide_conquer_v*/manifest.json'))
    paths += list(ROOT.glob('batches/*/manifest.json'))
    for path in paths:
        records = json.loads(path.read_text()).get('records', [])
        if records:
            prior.update(candidate_id(ast.literal_eval(x)) for x in records[0]['pool'])
    return prior


def _pool(per_class: int) -> list[tuple]:
    prior, pool = _prior(), []
    for cls in CLASS_ORDER:
        count, choices = max(32, per_class * 8), []
        while len(choices) < per_class:
            candidates = [candidate_at_class(cls, i) for i in _probe_indices(0, CLASS_SIZES[cls], min(count, CLASS_SIZES[cls]))]
            fresh: dict[str, tuple] = {}
            for genome in candidates:
                identity = candidate_id(genome)
                if identity not in prior and identity not in fresh:
                    fresh[identity] = genome
            choices = list(fresh.values())[:per_class]
            if count >= CLASS_SIZES[cls]: break
            count = min(count * 2, CLASS_SIZES[cls])
        if len(choices) != per_class: raise RuntimeError(f'fresh pool exhausted for {cls}')
        pool += choices
    return pool


def prepare(
    batch: str,
    per_class: int = 10,
    budget: int = 20,
    evidence_role: str = 'confirmatory',
) -> Path:
    """Lock one standard-policy batch before outcomes exist.

    Args:
        batch: Unique append-only batch identifier.
        per_class: Fresh pool members per material class.
        budget: Fixed application evaluation budget.
        evidence_role: ``confirmatory`` or non-confirmatory ``smoke`` evidence.

    Returns:
        Locked manifest path.
    """
    policy_source_sha256 = _require_locked_policy()
    lock = json.loads(POLICY_LOCK.read_text())
    if evidence_role not in ('confirmatory', 'smoke'):
        raise ValueError(f'unknown evidence role {evidence_role}')
    if evidence_role == 'confirmatory' and not batch.startswith('prereg-'):
        raise ValueError('confirmatory batch IDs must start with prereg-')
    if evidence_role == 'confirmatory':
        expected_per_class = int(lock['candidates_per_class'])
        expected_budget = int(lock['selection_budget'])
        if (per_class, budget) != (expected_per_class, expected_budget):
            raise ValueError(
                'confirmatory design mismatch: expected '
                f'per_class={expected_per_class}, budget={expected_budget}'
            )
    root, manifest = _root(batch), _root(batch) / 'manifest.json'
    if manifest.exists(): raise FileExistsError(manifest)
    if any(_outcome(batch, app).exists() for app, _, _ in APPS): raise RuntimeError('outcomes predate lock')
    pool, records = _pool(per_class), []
    ids = [candidate_id(g) for g in pool]
    applications = (
        APPS if evidence_role == 'smoke' else
        tuple(item for item in APPS if item[0] in lock['confirmatory_applications'])
    )
    for app, _, _ in applications:
        train, all_train = _training(app), _training(app, False)
        random_seed = int.from_bytes(
            hashlib.sha256(f'{batch}:{app}'.encode()).digest()[:8], 'big'
        )
        selection = _select_policies(
            pool, app, budget, train, all_train, random_seed
        )
        records.append({
            **selection,
            'pool': [repr(genome) for genome in pool],
            'candidate_ids': ids,
        })
    ROOT.mkdir(parents=True, exist_ok=True)
    campaign = ROOT / 'manifest.json'
    if not campaign.exists():
        _atomic_json(campaign, {'schema_version':1,'campaign':'prospective_catalyst_search','append_only':True,'historical_evidence':str(HISTORY)})
    payload = {'schema_version':3,'batch_id':batch,'evidence_role':evidence_role,'locked_before_outcomes':True,
      'created_utc':datetime.now(timezone.utc).isoformat(),
      'git_commit':subprocess.run(['git','rev-parse','HEAD'],check=True,capture_output=True,text=True).stdout.strip(),
      'policy_source_sha256':policy_source_sha256,
      'design': {'candidates_per_class': per_class, 'selection_budget': budget},
      'acceptance':{
          'confirmatory_applications':list(lock['confirmatory_applications']),
          'exploratory_applications':list(lock.get('exploratory_applications', [])),
          'catalyst_must_beat':['uncertainty','validity'],
          'random_one_sided_alpha':0.05,
      },'records':records}
    root.mkdir(parents=True, exist_ok=True); _atomic_json(manifest, payload)
    _commit_batch_manifest(batch, _hash(manifest))
    return manifest


def evaluate(batch: str, app: str) -> pd.DataFrame:
    """Evaluate every eligible candidate in a locked batch.

    Args:
        batch: Locked batch identifier.
        app: Application to screen.

    Returns:
        Complete screening table.
    """
    outcome = _outcome(batch, app)
    if outcome.exists():
        raise FileExistsError(f'prospective outcome is immutable: {outcome}')
    data=json.loads((_root(batch)/'manifest.json').read_text()); record=next(x for x in data['records'] if x['application']==app)
    allowed=set(record['eligible_ids']); pool=[ast.literal_eval(x) for x,i in zip(record['pool'],record['candidate_ids']) if i in allowed]
    if app=='turquoise_hydrogen':
        from pipeline.screening.surface_screener import run_screening
        return run_screening(pool,db_filename=f'prospective/{batch}_pyrolysis.csv',workers_per_gpu=1)
    from pipeline.screening.fc_screener import run_orr_screening
    return run_orr_screening(pool,db_filename=f'prospective/{batch}_orr.csv',workers_per_gpu=1)


def _score_record(
    record: dict, frame: pd.DataFrame, batch: str
) -> tuple[dict, list[dict], np.ndarray]:
    """Score one locked application record with deterministic top-k hits."""
    from pipeline.screening.small_data_ranker import screening_metric_eligibility

    app = str(record['application'])
    column = 'E_act' if app == 'turquoise_hydrogen' else 'orr_overpotential_V'
    required = {'genome', 'valid', column}
    if app == 'turquoise_hydrogen':
        required |= {'E_act_censored', 'pyrolysis_viable'}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f'{app} outcome missing columns: {sorted(missing)}')
    identities = [candidate_id(ast.literal_eval(raw)) for raw in frame.genome]
    if len(identities) != len(set(identities)):
        raise ValueError(f'{app} outcome contains duplicate candidates')
    expected = set(record['eligible_ids'])
    observed = set(identities)
    if observed != expected:
        raise ValueError(
            f'{app} outcome candidate mismatch: '
            f'missing={len(expected - observed)}, extra={len(observed - expected)}'
        )
    eligible = screening_metric_eligibility(frame, app)
    ranked = sorted(
        (
            (float(frame.iloc[position][column]), identity)
            for position, identity in enumerate(identities)
            if eligible[position]
        ),
        key=lambda item: (item[0], item[1]),
    )
    if not ranked:
        raise ValueError(f'{app} outcome contains no eligible metrics')
    hit_count = max(1, int(np.ceil(HIT_FRACTION * len(ranked))))
    hit_ids = {identity for _, identity in ranked[:hit_count]}
    cutoff = ranked[hit_count - 1][0]
    count_hits = lambda ids: sum(identity in hit_ids for identity in ids)
    policy = {
        name: count_hits(selected)
        for name, selected in record['policy_selected_ids'].items()
    }
    genomes = {
        identity: ast.literal_eval(raw)
        for identity, raw in zip(identities, frame.genome)
    }
    by_class: dict[str, list[str]] = {}
    for identity in record['eligible_ids']:
        by_class.setdefault(str(genomes[identity][0]), []).append(identity)
    rng = np.random.default_rng(record['random_seed'])
    random_hits = np.empty(record['random_trials'], dtype=np.int16)
    for trial in range(record['random_trials']):
        chosen = [str(rng.choice(by_class[name])) for name in sorted(by_class)]
        remaining = [
            identity for identity in record['eligible_ids'] if identity not in chosen
        ]
        chosen.extend(
            str(identity) for identity in rng.choice(
                remaining, record['budget'] - len(chosen), replace=False
            )
        )
        random_hits[trial] = count_hits(chosen)
    alpha = float(record.get('random_one_sided_alpha', .025))
    upper = float(np.quantile(random_hits, 1.0 - alpha))
    passed = (
        policy['catalyst'] > policy['uncertainty']
        and policy['catalyst'] > policy['validity']
        and policy['catalyst'] > upper
    )
    observations = []
    for position, identity in enumerate(identities):
        value = float(frame.iloc[position][column]) if eligible[position] else None
        observations.append({
            'batch_id': batch,
            'application': app,
            'candidate_id': identity,
            'valid': bool(eligible[position]),
            'outcome': value,
        })
    result = {
        'application': app,
        'eligible': len(record['eligible_ids']),
        'valid': len(ranked),
        'budget': record['budget'],
        'hit_definition': f'deterministic_top_{HIT_FRACTION:g}_eligible',
        'hit_count': hit_count,
        'hit_cutoff': cutoff,
        'policy_hits': policy,
        'policy_matched_random_mean_hits': float(np.mean(random_hits)),
        'policy_matched_random_quantiles': [
            float(np.quantile(random_hits, alpha)), upper
        ],
        'random_one_sided_alpha': alpha,
        'acceptance_passed': bool(passed),
    }
    return result, observations, random_hits


def _commit_observations(observations: list[dict]) -> None:
    """Idempotently commit candidate-keyed observations to the ledger."""
    ledger = ROOT / 'outcomes.jsonl'
    existing_rows = (
        [json.loads(line) for line in ledger.read_text().splitlines()]
        if ledger.exists()
        else []
    )
    keyed = {
        (row['batch_id'], row['application'], row['candidate_id']): row
        for row in existing_rows
    }
    for row in observations:
        key = (row['batch_id'], row['application'], row['candidate_id'])
        if key in keyed and keyed[key] != row:
            raise RuntimeError(f'conflicting observation already exists: {key}')
        keyed[key] = row
    ordered = sorted(
        keyed.values(),
        key=lambda row: (row['batch_id'], row['application'], row['candidate_id']),
    )
    _atomic_text(
        ledger, ''.join(json.dumps(row, sort_keys=True) + '\n' for row in ordered)
    )


def _commit_batch_manifest(batch: str, manifest_sha256: str) -> None:
    """Idempotently record one immutable manifest in the batch ledger."""
    ledger = ROOT / 'batches.jsonl'
    rows = (
        [json.loads(line) for line in ledger.read_text().splitlines()]
        if ledger.exists()
        else []
    )
    row = {'batch_id': batch, 'manifest_sha256': manifest_sha256}
    existing = next(
        (item for item in rows if str(item['batch_id']) == batch), None
    )
    if existing is not None and existing != row:
        raise RuntimeError(f'conflicting batch manifest already exists: {batch}')
    if existing is None:
        rows.append(row)
    _atomic_text(
        ledger,
        ''.join(json.dumps(item, sort_keys=True) + '\n' for item in rows),
    )


def analyze(batch: str) -> Path:
    """Analyze a completed batch once and atomically commit its evidence.

    Args:
        batch: Completed locked batch identifier.

    Returns:
        Path to the immutable batch analysis.
    """
    root = _root(batch)
    target = root / 'analysis.json'
    if target.exists():
        raise FileExistsError(target)
    manifest = root / 'manifest.json'
    data = json.loads(manifest.read_text())
    results, observations = [], []
    for record in data['records']:
        result, rows, _ = _score_record(
            record, pd.read_csv(_outcome(batch, record['application'])), batch
        )
        results.append(result)
        observations.extend(rows)
    report = {
        'schema_version': 2,
        'batch_id': batch,
        'manifest_sha256': _hash(manifest),
        'outcome_sha256': {
            str(record['application']): _hash(
                _outcome(batch, str(record['application']))
            )
            for record in data['records']
        },
        'results': results,
        'acceptance_passed': all(
            result['acceptance_passed']
            for result in results
            if result['application'] in data.get('acceptance', {}).get(
                'confirmatory_applications', [app for app, _, _ in APPS]
            )
        ),
    }
    # The ledger is idempotent, so a crash before the analysis rename is safe
    # to retry.  Publishing analysis first could permanently strand the batch.
    _commit_observations(observations)
    _atomic_json(target, report)
    reports = [
        json.loads(path.read_text())
        for _, path in _ordered_finalized_batches()
    ]
    _atomic_json(ROOT / 'analysis.json', {
        'schema_version': 1,
        'batches': [item['batch_id'] for item in reports],
        'batches_passed': sum(item['acceptance_passed'] for item in reports),
        'latest': reports[-1],
    })
    return target


def rescore_campaign() -> Path:
    """Write a retrospective corrected analysis without changing batch evidence.

    Returns:
        Path to the new campaign-level retrospective artifact.
    """
    batch_reports, random_by_app = [], {app: [] for app, _, _ in APPS}
    for batch, analysis in _ordered_finalized_batches():
        immutable = json.loads(analysis.read_text())
        if str(immutable['batch_id']) != batch:
            raise RuntimeError(f'analysis batch ID mismatch: {analysis}')
        manifest = _root(batch) / 'manifest.json'
        if _hash(manifest) != immutable['manifest_sha256']:
            raise RuntimeError(f'finalized manifest checksum mismatch: {manifest}')
        data = json.loads(manifest.read_text())
        results = []
        for record in data['records']:
            app = str(record['application'])
            outcome = _outcome(batch, app)
            if _hash(outcome) != immutable['outcome_sha256'][app]:
                raise RuntimeError(f'finalized outcome checksum mismatch: {outcome}')
            result, _, random_hits = _score_record(record, pd.read_csv(outcome), batch)
            results.append(result)
            random_by_app[app].append(random_hits)
        batch_reports.append({'batch_id': batch, 'results': results})
    pooled = []
    for app, _, _ in APPS:
        app_results = [
            next(result for result in batch['results'] if result['application'] == app)
            for batch in batch_reports
        ]
        random_sum = np.sum(np.vstack(random_by_app[app]), axis=0)
        policy_hits = {
            policy: int(sum(result['policy_hits'][policy] for result in app_results))
            for policy in ('catalyst', 'uncertainty', 'validity')
        }
        observed = policy_hits['catalyst']
        pooled.append({
            'application': app,
            'policy_hits': policy_hits,
            'policy_matched_random_mean_hits': float(np.mean(random_sum)),
            'policy_matched_random_95pct': [
                float(np.quantile(random_sum, .025)),
                float(np.quantile(random_sum, .975)),
            ],
            'random_one_sided_p_value': float(
                (1 + np.count_nonzero(random_sum >= observed)) / (len(random_sum) + 1)
            ),
        })
    output = ROOT / 'rescored_analysis.json'
    _atomic_json(output, {
        'schema_version': 1,
        'interpretation': 'retrospective_exploratory_not_confirmatory',
        'hit_definition': f'deterministic_top_{HIT_FRACTION:g}_eligible',
        'batches': batch_reports,
        'pooled': pooled,
    })
    return output


def replay_and_power() -> Path:
    """Replay the fixed policy chronologically and bootstrap ten-batch power.

    Returns:
        Path to the CPU-only replay and power artifact.
    """
    replayed, random_by_app = [], {app: [] for app, _, _ in APPS}
    for batch, analysis in _ordered_finalized_batches():
        immutable = json.loads(analysis.read_text())
        if str(immutable['batch_id']) != batch:
            raise RuntimeError(f'analysis batch ID mismatch: {analysis}')
        manifest = _root(batch) / 'manifest.json'
        data = json.loads(manifest.read_text())
        pool = [ast.literal_eval(raw) for raw in data['records'][0]['pool']]
        results = []
        for app, _, _ in APPS:
            train = _training(app, before_batch=batch)
            all_train = _training(app, False, before_batch=batch)
            seed = int.from_bytes(
                hashlib.sha256(f'replay:{batch}:{app}'.encode()).digest()[:8],
                'big',
            )
            record = _select_policies(pool, app, 20, train, all_train, seed)
            outcome = _outcome(batch, app)
            if _hash(outcome) != immutable['outcome_sha256'][app]:
                raise RuntimeError(f'finalized outcome checksum mismatch: {outcome}')
            frame = pd.read_csv(outcome)
            identities = [candidate_id(ast.literal_eval(raw)) for raw in frame.genome]
            frame = frame.loc[
                [identity in set(record['eligible_ids']) for identity in identities]
            ].reset_index(drop=True)
            result, _, random_hits = _score_record(record, frame, batch)
            results.append(result)
            random_by_app[app].append(random_hits)
        replayed.append({'batch_id': batch, 'results': results})

    rng = np.random.default_rng(20261008)
    bootstrap_replicates = 2000
    null_trials = 1000
    planned_batches = int(json.loads(POLICY_LOCK.read_text())['planned_batches'])
    power = []
    for app, _, _ in APPS:
        app_results = [
            next(result for result in batch['results'] if result['application'] == app)
            for batch in replayed
        ]
        def estimate(indices: list[int], planned_batches: int) -> float:
            passes = 0
            for _ in range(bootstrap_replicates):
                sampled = rng.choice(indices, size=planned_batches, replace=True)
                catalyst = sum(
                    app_results[index]['policy_hits']['catalyst'] for index in sampled
                )
                uncertainty = sum(
                    app_results[index]['policy_hits']['uncertainty'] for index in sampled
                )
                validity = sum(
                    app_results[index]['policy_hits']['validity'] for index in sampled
                )
                null = np.zeros(null_trials, dtype=np.int16)
                for index in sampled:
                    draws = random_by_app[app][index]
                    null += draws[rng.integers(0, len(draws), size=null_trials)]
                p_value = (
                    1 + np.count_nonzero(null >= catalyst)
                ) / (null_trials + 1)
                passes += (
                    catalyst > uncertainty
                    and catalyst > validity
                    and p_value <= .05
                )
            return passes / bootstrap_replicates

        all_indices = list(range(len(app_results)))
        ten_batch_power = estimate(all_indices, 10)
        without_powered = [
            index for index, batch in enumerate(replayed)
            if batch['batch_id'] != '20260930-catalyst-04-powered'
        ]
        estimated_power = estimate(all_indices, planned_batches)
        leave_batch_04_out_power = estimate(without_powered, planned_batches)
        power.append({
            'application': app,
            'bootstrap_replicates': bootstrap_replicates,
            'planned_batches': planned_batches,
            'estimated_pass_probability': estimated_power,
            'ten_batch_pass_probability': ten_batch_power,
            'leave_batch_04_out_pass_probability': leave_batch_04_out_power,
            'leave_batch_04_out_power_curve': {
                str(batch_count): estimate(without_powered, batch_count)
                for batch_count in (10, 15, 20, 25, 30, 35, 40)
            },
            'confirmatory_recommendation': leave_batch_04_out_power >= .80,
        })
    output = ROOT / 'fixed_policy_replay.json'
    _atomic_json(output, {
        'schema_version': 1,
        'interpretation': (
            'retrospective chronological replay; power bootstrap is indicative '
            'because historical pools are smaller than the planned pool'
        ),
        'batches': replayed,
        'power': power,
    })
    return output


def confirmatory_analysis() -> Path:
    """Evaluate only fully collected batches registered under the policy lock.

    Before all planned batches exist, this writes progress metadata without a
    statistic, p-value, or verdict. The final analysis is immutable.

    Returns:
        Path to either the interim progress or final confirmatory artifact.
    """
    lock = json.loads(POLICY_LOCK.read_text())
    locked_hash = _require_locked_policy()
    planned = int(lock['planned_batches'])
    confirmatory_apps = tuple(lock['confirmatory_applications'])
    exploratory_apps = tuple(lock.get('exploratory_applications', ()))
    included: list[tuple[str, dict, dict]] = []
    excluded = []
    for batch, analysis_path in _ordered_finalized_batches():
        analysis = json.loads(analysis_path.read_text())
        manifest_path = _root(batch) / 'manifest.json'
        manifest = json.loads(manifest_path.read_text())
        reason = None
        if manifest.get('evidence_role') != 'confirmatory':
            reason = 'not_confirmatory_role'
        elif manifest.get('policy_source_sha256') != locked_hash:
            reason = 'policy_hash_mismatch'
        elif not batch.startswith('prereg-'):
            reason = 'batch_id_not_preregistered'
        elif analysis.get('manifest_sha256') != _hash(manifest_path):
            raise RuntimeError(f'confirmatory manifest checksum mismatch: {batch}')
        if reason:
            excluded.append({'batch_id': batch, 'reason': reason})
        else:
            included.append((batch, manifest, analysis))
    if len(included) > planned:
        raise RuntimeError(
            f'confirmatory batch count exceeds preregistration: {len(included)} > {planned}'
        )
    if len(included) < planned:
        progress = ROOT / 'confirmatory_progress.json'
        _atomic_json(progress, {
            'schema_version': 1,
            'status': 'collecting',
            'policy_source_sha256': locked_hash,
            'planned_batches': planned,
            'completed_batches': len(included),
            'remaining_batches': planned - len(included),
            'batch_ids': [batch for batch, _, _ in included],
            'excluded': excluded,
            'inferential_statistics_withheld': True,
        })
        return progress

    output = ROOT / 'confirmatory_analysis.json'
    if output.exists():
        raise FileExistsError(f'confirmatory analysis is immutable: {output}')
    result_apps = tuple(
        str(record['application']) for record in included[0][1]['records']
    )
    if set(result_apps) != set(confirmatory_apps):
        raise RuntimeError('confirmatory manifests do not match locked applications')
    results_by_app: dict[str, list[tuple[dict, np.ndarray]]] = {
        app: [] for app in result_apps
    }
    for batch, manifest, analysis in included:
        records = {
            str(record['application']): record for record in manifest['records']
        }
        if set(records) != set(result_apps):
            raise RuntimeError('confirmatory application set changed between batches')
        for app in result_apps:
            outcome = _outcome(batch, app)
            if _hash(outcome) != analysis['outcome_sha256'][app]:
                raise RuntimeError(f'confirmatory outcome checksum mismatch: {outcome}')
            result, _, random_hits = _score_record(
                records[app], pd.read_csv(outcome), batch
            )
            results_by_app[app].append((result, random_hits))
    pooled = []
    for app in result_apps:
        scored = results_by_app[app]
        random_sum = np.sum(
            np.vstack([random_hits for _, random_hits in scored]), axis=0
        )
        policy_hits = {
            policy: int(sum(result['policy_hits'][policy] for result, _ in scored))
            for policy in ('catalyst', 'uncertainty', 'validity')
        }
        catalyst = policy_hits['catalyst']
        p_value = float(
            (1 + np.count_nonzero(random_sum >= catalyst)) / (len(random_sum) + 1)
        )
        alpha = float(lock['random_one_sided_alpha'])
        is_confirmatory = app in confirmatory_apps
        passed = (
            catalyst > policy_hits['uncertainty']
            and catalyst > policy_hits['validity']
            and p_value <= alpha
        ) if is_confirmatory else None
        pooled.append({
            'application': app,
            'evidence_role': (
                'confirmatory' if is_confirmatory else 'exploratory'
            ),
            'policy_hits': policy_hits,
            'policy_matched_random_mean_hits': float(np.mean(random_sum)),
            'random_one_sided_p_value': p_value,
            'random_one_sided_alpha': alpha if is_confirmatory else None,
            'acceptance_passed': passed,
        })
    report = {
        'schema_version': 1,
        'status': 'complete',
        'policy_source_sha256': locked_hash,
        'planned_batches': planned,
        'batch_ids': [batch for batch, _, _ in included],
        'confirmatory_applications': list(confirmatory_apps),
        'exploratory_applications': list(exploratory_apps),
        'pooled': pooled,
        'acceptance_passed': all(
            result['acceptance_passed']
            for result in pooled
            if result['evidence_role'] == 'confirmatory'
        ),
        'excluded': excluded,
    }
    _atomic_json(output, report)
    return output


def archive_history() -> Path:
    """Write one checksum ledger for immutable v2-v8 evidence.

    Returns:
        Historical ledger path.
    """
    rows=[]
    for manifest in sorted(Path('results/pilot').glob('divide_conquer_v*/manifest.json')):
        version=manifest.parent.name.removeprefix('divide_conquer_'); artifacts={'manifest':manifest}
        for name,path in [('analysis',manifest.parent/'analysis.json'),('pyrolysis',Path(f'results/screening/pilot/divide_conquer_{version}_pyrolysis.csv')),('orr',Path(f'results/fuel_cell/pilot/divide_conquer_{version}_orr.csv'))]:
            if path.exists(): artifacts[name]=path
        rows.append({'historical_round':version,'read_only':True,'artifacts':{k:{'path':str(v),'sha256':_hash(v)} for k,v in artifacts.items()}})
    HISTORY.parent.mkdir(parents=True,exist_ok=True); HISTORY.write_text(''.join(json.dumps(x,sort_keys=True)+'\n' for x in rows)); return HISTORY


def main() -> None:
    """Dispatch the prospective campaign command-line interface."""
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=('prepare','evaluate-pyrolysis','evaluate-orr','analyze','rescore','replay-power','confirmatory-analysis','archive-history')); parser.add_argument('--batch'); parser.add_argument('--per-class',type=int,default=10); parser.add_argument('--budget',type=int,default=20); parser.add_argument('--evidence-role',choices=('confirmatory','smoke'),default='confirmatory'); args=parser.parse_args()
    if args.action=='archive-history': print(archive_history()); return
    if args.action=='rescore': print(rescore_campaign()); return
    if args.action=='replay-power': print(replay_and_power()); return
    if args.action=='confirmatory-analysis': print(confirmatory_analysis()); return
    if not args.batch: parser.error('--batch is required')
    if args.action=='prepare': print(prepare(args.batch,args.per_class,args.budget,args.evidence_role))
    elif args.action=='evaluate-pyrolysis': print(evaluate(args.batch,'turquoise_hydrogen').shape)
    elif args.action=='evaluate-orr': print(evaluate(args.batch,'fuel_cell_orr').shape)
    else: print(analyze(args.batch))


if __name__=='__main__': main()
