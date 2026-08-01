"""Prospective time-split evidence benchmark for discovery claims."""

import hashlib
import json

from pipeline.search.discovery import candidate_id, discovery_region


def _manifest_hash(records) -> str:
    normalized = sorted(
        (candidate_id(tuple(item['genome'])), int(item['publication_year']),
         str(item['source_id']))
        for item in records)
    return hashlib.sha256(
        json.dumps(normalized, separators=(',', ':')).encode()).hexdigest()


def time_split_recovery(ranked_candidates, held_out_discoveries,
                        cutoff_year: int, training_records, k=100) -> dict:
    """Measure whether a ranking recovers discoveries hidden by publication year.

    `held_out_discoveries` must contain mappings with genome, publication_year,
    source_id, and a non-empty citation. Missing provenance fails closed.
    """
    if k <= 0:
        raise ValueError('k must be positive')
    training_records = list(training_records)
    held_out_discoveries = list(held_out_discoveries)
    all_records = training_records + held_out_discoveries
    malformed = [x for x in all_records if not x.get('source_id') or
                 not x.get('citation') or not isinstance(x.get('publication_year'), int)]
    if malformed:
        return {'valid': False, 'reason': 'held-out records lack provenance'}
    if any(x['publication_year'] > cutoff_year for x in training_records):
        return {'valid': False, 'reason': 'training record occurs after cutoff'}
    if any(x['publication_year'] <= cutoff_year for x in held_out_discoveries):
        return {'valid': False, 'reason': 'held-out record does not occur after cutoff'}
    training_ids = {candidate_id(tuple(x['genome'])) for x in training_records}
    held_out_ids = {candidate_id(tuple(x['genome'])) for x in held_out_discoveries}
    if training_ids & held_out_ids:
        return {'valid': False, 'reason': 'candidate leakage across time split'}
    top = list(ranked_candidates[:k])
    top_ids = {candidate_id(g) for g in top}
    top_regions = {discovery_region(g) for g in top}
    exact = sum(candidate_id(tuple(x['genome'])) in top_ids for x in held_out_discoveries)
    regional = sum(discovery_region(tuple(x['genome'])) in top_regions
                   for x in held_out_discoveries)
    n = len(held_out_discoveries)
    return {'valid': n > 0, 'n_held_out': n, 'k': k,
            'cutoff_year': int(cutoff_year),
            'training_manifest_sha256': _manifest_hash(training_records),
            'held_out_manifest_sha256': _manifest_hash(held_out_discoveries),
            'exact_recall_at_k': exact / n if n else 0.0,
            'region_recall_at_k': regional / n if n else 0.0,
            'exact_recovered': exact, 'region_recovered': regional}
