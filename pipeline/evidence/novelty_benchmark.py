"""Prospective time-split evidence benchmark for discovery claims."""

from pipeline.search.discovery import candidate_id, discovery_region


def time_split_recovery(ranked_candidates, held_out_discoveries, k=100) -> dict:
    """Measure whether a ranking recovers discoveries hidden by publication year.

    `held_out_discoveries` must contain mappings with genome, publication_year,
    source_id, and a non-empty citation. Missing provenance fails closed.
    """
    if k <= 0:
        raise ValueError('k must be positive')
    malformed = [x for x in held_out_discoveries if not x.get('source_id') or
                 not x.get('citation') or not isinstance(x.get('publication_year'), int)]
    if malformed:
        return {'valid': False, 'reason': 'held-out records lack provenance'}
    top = list(ranked_candidates[:k])
    top_ids = {candidate_id(g) for g in top}
    top_regions = {discovery_region(g) for g in top}
    exact = sum(candidate_id(tuple(x['genome'])) in top_ids for x in held_out_discoveries)
    regional = sum(discovery_region(tuple(x['genome'])) in top_regions
                   for x in held_out_discoveries)
    n = len(held_out_discoveries)
    return {'valid': n > 0, 'n_held_out': n, 'k': k,
            'exact_recall_at_k': exact / n if n else 0.0,
            'region_recall_at_k': regional / n if n else 0.0,
            'exact_recovered': exact, 'region_recovered': regional}
