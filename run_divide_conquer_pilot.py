#!/usr/bin/env python3
"""One append-only prospective campaign for catalyst-guided search.

Historical v2-v8 files are evidence inputs only; batch IDs never select code.
"""
from __future__ import annotations

import argparse, ast, hashlib, json, subprocess
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd

from pipeline.evidence.pilot_benchmark import default_specs, load_legacy_outcomes
from pipeline.evidence.search_policy_benchmark import _coverage_select
from pipeline.search.branch_search import _probe_indices
from pipeline.search.design_space import encode_population
from pipeline.search.discovery import candidate_id
from pipeline.search.indexed_space import CLASS_ORDER, CLASS_SIZES, candidate_at_class
from pipeline.search.scope import pemfc_cathode_scope

ROOT = Path('results/prospective_search')
HISTORY = Path('docs/evidence/legacy_pilot_rounds.jsonl')
APPS = (('turquoise_hydrogen', 'E_act', 'pyrolysis'),
        ('fuel_cell_orr', 'orr_overpotential_V', 'orr'))


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


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


def _completed_outcomes(app: str) -> list[Path]:
    """Return checksum-verified outcomes from finalized prospective batches."""
    paths = []
    for analysis in sorted(ROOT.glob('batches/*/analysis.json')):
        report = json.loads(analysis.read_text())
        batch = str(report['batch_id'])
        path = _outcome(batch, app)
        expected = report.get('outcome_sha256', {}).get(app)
        if expected is None:
            continue
        if not path.is_file() or _hash(path) != expected:
            raise RuntimeError(f'finalized outcome checksum mismatch: {path}')
        paths.append(path)
    return paths


def _training(app: str, valid_only: bool = True) -> pd.DataFrame:
    spec = next(x for x in default_specs() if x.application == app)
    frames = [load_legacy_outcomes(spec).drop(columns=['parsed_genome', 'candidate_id', 'replicates'])]
    seed = Path('results/screening/pilot/prospective_pyrolysis.csv' if app == 'turquoise_hydrogen' else 'results/fuel_cell/pilot/prospective_orr.csv')
    frames.append(pd.read_csv(seed))
    folder, suffix = ('screening', 'pyrolysis') if app == 'turquoise_hydrogen' else ('fuel_cell', 'orr')
    paths = sorted(Path(f'results/{folder}/pilot').glob(f'divide_conquer_v*_{suffix}.csv'))
    paths += _completed_outcomes(app)
    frames += [pd.read_csv(path) for path in paths]
    frame = pd.concat(frames, ignore_index=True, sort=False)
    if valid_only:
        frame = frame[frame.valid.eq(True)].copy()
    frame['_id'] = [candidate_id(ast.literal_eval(x)) for x in frame.genome]
    return frame.drop_duplicates('_id', keep='last').drop(columns='_id')


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
            choices = [g for g in candidates if candidate_id(g) not in prior][:per_class]
            if count >= CLASS_SIZES[cls]: break
            count = min(count * 2, CLASS_SIZES[cls])
        if len(choices) != per_class: raise RuntimeError(f'fresh pool exhausted for {cls}')
        pool += choices
    return pool


def prepare(batch: str, per_class: int = 4, extra_slots: int = 6) -> Path:
    """Lock one standard-policy batch before outcomes exist.

    Args:
        batch: Unique append-only batch identifier.
        per_class: Fresh pool members per material class.
        extra_slots: Budget remaining after class floors.

    Returns:
        Locked manifest path.
    """
    from dataclasses import asdict
    from sklearn.ensemble import ExtraTreesClassifier
    from pipeline.screening.small_data_ranker import (
        fit_tree_ranker,
        orr_catalyst_acquisition,
        orr_tree_objectives,
        turquoise_tree_objectives,
    )
    root, manifest = _root(batch), _root(batch) / 'manifest.json'
    if manifest.exists(): raise FileExistsError(manifest)
    if any(_outcome(batch, app).exists() for app, _, _ in APPS): raise RuntimeError('outcomes predate lock')
    pool, records = _pool(per_class), []
    ids = [candidate_id(g) for g in pool]
    for app, _, _ in APPS:
        train, all_train = _training(app), _training(app, False)
        eligible = list(range(len(pool))) if app == 'turquoise_hydrogen' else [i for i,g in enumerate(pool) if pemfc_cathode_scope(g)['status'] == 'candidate']
        ranker = fit_tree_ranker(train, app)
        obj = turquoise_tree_objectives(pool, ranker) if app == 'turquoise_hydrogen' else orr_tree_objectives(pool, ranker)
        _, uncertainty = ranker.predict(pool)
        genomes = [ast.literal_eval(x) for x in all_train.genome]
        validity = ExtraTreesClassifier(n_estimators=256, min_samples_leaf=2, max_features=1.0, class_weight='balanced', random_state=20260722, n_jobs=-1)
        validity.fit(encode_population(genomes), all_train.valid.eq(True).to_numpy(int))
        vscore = validity.predict_proba(encode_population(pool))[:, list(validity.classes_).index(1)]
        candidates = [pool[i] for i in eligible]
        budget = min(len(candidates), len({g[0] for g in candidates}) + extra_slots)
        catalyst_score = (
            -obj[eligible, 0]
            if app == 'turquoise_hydrogen'
            else orr_catalyst_acquisition(
                obj[eligible, 0], uncertainty[eligible], vscore[eligible]
            )
        )
        local = {'catalyst': _coverage_select(catalyst_score, candidates, budget),
                 'uncertainty': _coverage_select(uncertainty[eligible], candidates, budget),
                 'validity': _coverage_select(vscore[eligible], candidates, budget)}
        records.append({'application': app, 'training_rows': len(train),
          'training_digest': hashlib.sha256('\n'.join(sorted(candidate_id(ast.literal_eval(x)) for x in train.genome)).encode()).hexdigest(),
          'pool': [repr(g) for g in pool], 'candidate_ids': ids,
          'eligible_ids': [ids[i] for i in eligible], 'budget': budget,
          'policy_selected_ids': {name: [ids[eligible[i]] for i in ix] for name,ix in local.items()},
          'ranker_diagnostics': asdict(ranker.diagnostics), 'random_seed': 20260720, 'random_trials': 50000})
    ROOT.mkdir(parents=True, exist_ok=True)
    campaign = ROOT / 'manifest.json'
    if not campaign.exists():
        campaign.write_text(json.dumps({'schema_version':1,'campaign':'prospective_catalyst_search','append_only':True,'historical_evidence':str(HISTORY)}, indent=2)+'\n')
    payload = {'schema_version':2,'batch_id':batch,'locked_before_outcomes':True,
      'created_utc':datetime.now(timezone.utc).isoformat(),
      'git_commit':subprocess.run(['git','rev-parse','HEAD'],check=True,capture_output=True,text=True).stdout.strip(),
      'policy_source_sha256':hashlib.sha256(Path(__file__).read_bytes()+Path('pipeline/screening/small_data_ranker.py').read_bytes()).hexdigest(),
      'acceptance':{'catalyst_must_beat':['uncertainty','validity','policy_matched_random_97_5pct'],'required_independently_for_each_application':True},'records':records}
    root.mkdir(parents=True, exist_ok=False); manifest.write_text(json.dumps(payload,indent=2)+'\n')
    with (ROOT/'batches.jsonl').open('a') as f: f.write(json.dumps({'batch_id':batch,'manifest_sha256':_hash(manifest)})+'\n')
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


def analyze(batch: str) -> Path:
    """Analyze a completed batch once and append its observations.

    Args:
        batch: Completed batch identifier.

    Returns:
        Immutable batch analysis path.
    """
    root=_root(batch); target=root/'analysis.json'
    if target.exists(): raise FileExistsError(target)
    manifest=root/'manifest.json'; data=json.loads(manifest.read_text()); results=[]; observations=[]
    for record in data['records']:
        app=record['application']; column='E_act' if app=='turquoise_hydrogen' else 'orr_overpotential_V'
        frame=pd.read_csv(_outcome(batch,app)); values={}; genomes={}
        required = {'genome', 'valid', column}
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
        for _,row in frame.iterrows():
            genome=ast.literal_eval(row.genome); identity=candidate_id(genome); genomes[identity]=genome
            valid: bool = bool(row.valid) and bool(np.isfinite(float(row[column])))
            value = float(row[column]) if valid else None
            if valid: values[identity]=value
            observations.append({'batch_id':batch,'application':app,'candidate_id':identity,'valid':valid,'outcome':value})
        cutoff=float(np.quantile(list(values.values()),.2))
        hits=lambda ids:sum(i in values and values[i]<=cutoff for i in ids)
        policy={name:hits(ids) for name,ids in record['policy_selected_ids'].items()}
        by_class={}
        for identity in record['eligible_ids']: by_class.setdefault(genomes[identity][0],[]).append(identity)
        rng=np.random.default_rng(record['random_seed']); random=[]
        for _ in range(record['random_trials']):
            chosen=[rng.choice(by_class[name]) for name in sorted(by_class)]; remaining=[i for i in record['eligible_ids'] if i not in chosen]
            chosen.extend(rng.choice(remaining,record['budget']-len(chosen),replace=False)); random.append(hits(chosen))
        upper=float(np.quantile(random,.975)); passed=policy['catalyst']>policy['uncertainty'] and policy['catalyst']>policy['validity'] and policy['catalyst']>upper
        results.append({'application':app,'eligible':len(record['eligible_ids']),'valid':len(values),'budget':record['budget'],'hit_cutoff':cutoff,'policy_hits':policy,'policy_matched_random_mean_hits':float(np.mean(random)),'policy_matched_random_95pct':[float(np.quantile(random,.025)),upper],'acceptance_passed':bool(passed)})
    report={'schema_version':1,'batch_id':batch,'manifest_sha256':_hash(manifest),'outcome_sha256':{app:_hash(_outcome(batch,app)) for app,_,_ in APPS},'results':results,'acceptance_passed':all(x['acceptance_passed'] for x in results)}
    serialized_observations = ''.join(
        json.dumps(row, sort_keys=True) + '\n' for row in observations
    )
    target.write_text(json.dumps(report,indent=2)+'\n')
    with (ROOT/'outcomes.jsonl').open('a') as f:
        f.write(serialized_observations)
    reports=[json.loads(path.read_text()) for path in sorted(ROOT.glob('batches/*/analysis.json'))]
    (ROOT/'analysis.json').write_text(json.dumps({'schema_version':1,'batches':[x['batch_id'] for x in reports],'batches_passed':sum(x['acceptance_passed'] for x in reports),'latest':reports[-1]},indent=2)+'\n')
    return target


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
    parser=argparse.ArgumentParser(); parser.add_argument('action',choices=('prepare','evaluate-pyrolysis','evaluate-orr','analyze','archive-history')); parser.add_argument('--batch'); parser.add_argument('--per-class',type=int,default=4); parser.add_argument('--extra-slots',type=int,default=6); args=parser.parse_args()
    if args.action=='archive-history': print(archive_history()); return
    if not args.batch: parser.error('--batch is required')
    if args.action=='prepare': print(prepare(args.batch,args.per_class,args.extra_slots))
    elif args.action=='evaluate-pyrolysis': print(evaluate(args.batch,'turquoise_hydrogen').shape)
    elif args.action=='evaluate-orr': print(evaluate(args.batch,'fuel_cell_orr').shape)
    else: print(analyze(args.batch))


if __name__=='__main__': main()
