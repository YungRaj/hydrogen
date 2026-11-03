"""Fail-closed evidence preflight for a final discovery campaign."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from pipeline.evidence.manifest import verify_evidence_manifest
from pipeline.evidence.prior_art import PriorArtRegistry
from pipeline.search.indexed_space import TOTAL_SIZE


def campaign_readiness(
    coverage_certificate: str,
    prior_art_db: str,
    require_complete_coverage: bool = True,
    evidence_manifest: str | None = None,
    application: str | None = None,
    pyrolysis_mode: str = 'thermocatalytic',
) -> dict[str, Any]:
    """Evaluate whether all required discovery evidence gates have passed.

    Args:
        coverage_certificate: Search certificate describing visited intervals and gaps.
        prior_art_db: Path to the curated prior-art registry.
        require_complete_coverage: Whether incomplete search coverage blocks readiness.
        evidence_manifest: Manifest containing required experimental and computational evidence.
        application: Scientific application or objective family.
        pyrolysis_mode: Methane-conversion pathway used for evaluation.

    Returns:
        A readiness record listing passed and missing scientific gates.
    """
    failures: list[str] = []
    warnings: list[str] = []
    cert_path = Path(coverage_certificate)
    certificate: dict[str, Any]
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
    prior_count = int(registry.count())
    if prior_count == 0:
        failures.append('prior_art_registry_empty')
    evidence: dict[str, int] = {}
    if evidence_manifest is not None:
        verification = verify_evidence_manifest(evidence_manifest)
        evidence = dict(verification['counts'])
        failures.extend(verification['errors'])
        if verification['valid']:
            required: tuple[str, ...] = {
                'turquoise_hydrogen': (
                    'converged_dft',
                    'measured_reactor',
                    'measured_deactivation',
                ),
                'fuel_cell': (
                    'converged_orr_dft',
                    'measured_mea',
                    'measured_durability',
                    'hydrogen_impurity_test',
                    'time_split_benchmark',
                    'curated_prior_art_source',
                ),
            }.get(application or '', ())
            if application == 'turquoise_hydrogen' and pyrolysis_mode == 'ntec':
                required += ('ntec_control_pair',)
            for key in required:
                # Require at least 1 hash-verified artifact record per mandatory evidence category
                if int(evidence.get(key, 0) or 0) < 1:
                    failures.append(f'evidence_missing:{key}')
    warnings.append('industrial_viability_requires_reactor_or_stack_measurements')
    warnings.append('surrogate_coverage_is_not_experimental_validation')
    return {
        'ready': not failures,
        'failures': failures,
        'warnings': warnings,
        'prior_art_records': prior_count,
        'coverage': certificate,
        'evidence': evidence,
    }
