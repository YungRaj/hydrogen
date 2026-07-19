"""Fail-closed preflight for a final discovery campaign."""

import json
from pathlib import Path

from pipeline.indexed_space import TOTAL_SIZE
from pipeline.prior_art import PriorArtRegistry


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
        path = Path(evidence_manifest)
        if not path.exists():
            failures.append('evidence_manifest_missing')
        else:
            evidence = json.loads(path.read_text())
            required = {
                'turquoise_hydrogen': ('converged_dft_count', 'measured_reactor_count',
                                       'measured_deactivation_count'),
                'fuel_cell': ('converged_orr_dft_count', 'measured_mea_count',
                              'measured_durability_count', 'hydrogen_impurity_test_count',
                              'time_split_benchmark_count', 'curated_prior_art_sources'),
            }.get(application, ())
            if application == 'turquoise_hydrogen' and pyrolysis_mode == 'ntec':
                required += ('ntec_control_pair_count',)
            for key in required:
                if int(evidence.get(key, 0) or 0) < 1:
                    failures.append(f'evidence_missing:{key}')
    warnings.append('industrial_viability_requires_reactor_or_stack_measurements')
    warnings.append('surrogate_coverage_is_not_experimental_validation')
    return {'ready': not failures, 'failures': failures, 'warnings': warnings,
            'prior_art_records': prior_count, 'coverage': certificate,
            'evidence': evidence}
