"""Machine-readable evidence status for scientific acceptance criteria."""

import json
from pathlib import Path
from pipeline.search.indexed_space import TOTAL_SIZE
from pipeline.evidence.manifest import verify_evidence_manifest


def assess_campaign(results_dir='results', pyrolysis_mode='ntec') -> dict:
    root = Path(results_dir)
    coverage = []
    for path in (root / 'screening/turquoise_hydrogen_coverage_certificate.json',
                 root / 'fuel_cell/coverage_certificate.json'):
        value = json.loads(path.read_text()) if path.exists() else {}
        coverage.append(bool(value.get('complete')) and
                        value.get('declared_encoded_population') == TOTAL_SIZE)
    path = root / 'evidence_manifest.json'
    verification = verify_evidence_manifest(path)
    evidence = verification['counts']
    criteria = {
        'complete_search': all(coverage),
        'validated_champions': evidence.get('converged_dft', 0) > 0 and
                               evidence.get('converged_orr_dft', 0) > 0,
        'validated_reactor': evidence.get('measured_reactor', 0) > 0 and
                             evidence.get('measured_deactivation', 0) > 0,
        'validated_pemfc': evidence.get('measured_mea', 0) > 0 and
                           evidence.get('measured_durability', 0) > 0 and
                           evidence.get('hydrogen_impurity_test', 0) > 0,
        'defensible_novelty': evidence.get('time_split_benchmark', 0) > 0 and
                              evidence.get('curated_prior_art_source', 0) > 0,
    }
    if pyrolysis_mode == 'ntec':
        criteria['calibrated_ntec'] = evidence.get('ntec_control_pair', 0) > 0
    return {'ready': all(criteria.values()), 'criteria': criteria,
            'missing': [k for k, passed in criteria.items() if not passed],
            'evidence_manifest_valid': verification['valid'],
            'evidence_errors': verification['errors']}
