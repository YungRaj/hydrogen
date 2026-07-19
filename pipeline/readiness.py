"""Fail-closed preflight for a final discovery campaign."""

import json
from pathlib import Path

from pipeline.indexed_space import TOTAL_SIZE
from pipeline.prior_art import PriorArtRegistry


def campaign_readiness(coverage_certificate: str, prior_art_db: str,
                       require_complete_coverage: bool = True) -> dict:
    failures, warnings = [], []
    cert_path = Path(coverage_certificate)
    if not cert_path.exists():
        failures.append('coverage_certificate_missing')
        certificate = {}
    else:
        certificate = json.loads(cert_path.read_text())
        if certificate.get('declared_encoded_population') != TOTAL_SIZE:
            failures.append('coverage_denominator_mismatch')
        if require_complete_coverage and not certificate.get('complete', False):
            failures.append('coverage_incomplete')
    registry = PriorArtRegistry(prior_art_db)
    prior_count = registry.count()
    if prior_count == 0:
        failures.append('prior_art_registry_empty')
    warnings.append('industrial_viability_requires_reactor_or_stack_measurements')
    warnings.append('surrogate_coverage_is_not_experimental_validation')
    return {'ready': not failures, 'failures': failures, 'warnings': warnings,
            'prior_art_records': prior_count, 'coverage': certificate}
