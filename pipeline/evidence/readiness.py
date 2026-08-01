"""Fail-closed evidence preflight for a final discovery campaign."""

import json
from pathlib import Path

from pipeline.search.indexed_space import TOTAL_SIZE
from pipeline.evidence.prior_art import PriorArtRegistry
from pipeline.evidence.manifest import verify_evidence_manifest


def campaign_readiness(coverage_certificate: str, prior_art_db: str,
                       require_complete_coverage: bool = True,
                       evidence_manifest: str | None = None,
                       application: str | None = None,
                       pyrolysis_mode: str = 'ntec') -> dict:
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
    evidence = {}
    if evidence_manifest is not None:
        verification = verify_evidence_manifest(evidence_manifest)
        evidence = verification['counts']
        failures.extend(verification['errors'])
        if verification['valid']:
            required = {
                'turquoise_hydrogen': ('converged_dft', 'measured_reactor',
                                       'measured_deactivation'),
                'fuel_cell': ('converged_orr_dft', 'measured_mea',
                              'measured_durability', 'hydrogen_impurity_test',
                              'time_split_benchmark', 'curated_prior_art_source'),
            }.get(application, ())
            if application == 'turquoise_hydrogen' and pyrolysis_mode == 'ntec':
                required += ('ntec_control_pair',)
            for key in required:
                if int(evidence.get(key, 0) or 0) < 1:
                    failures.append(f'evidence_missing:{key}')
    warnings.append('industrial_viability_requires_reactor_or_stack_measurements')
    warnings.append('surrogate_coverage_is_not_experimental_validation')
    return {'ready': not failures, 'failures': failures, 'warnings': warnings,
            'prior_art_records': prior_count, 'coverage': certificate,
            'evidence': evidence}
