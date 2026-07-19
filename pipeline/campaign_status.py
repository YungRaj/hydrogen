"""Machine-readable status for the six scientific acceptance criteria."""

import json
from pathlib import Path
from pipeline.indexed_space import TOTAL_SIZE


def assess_campaign(results_dir='results', pyrolysis_mode='ntec') -> dict:
    root = Path(results_dir)
    coverage = []
    for path in (root / 'screening/turquoise_hydrogen_coverage_certificate.json',
                 root / 'fuel_cell/coverage_certificate.json'):
        value = json.loads(path.read_text()) if path.exists() else {}
        coverage.append(bool(value.get('complete')) and
                        value.get('declared_encoded_population') == TOTAL_SIZE)
    path = root / 'evidence_manifest.json'
    evidence = json.loads(path.read_text()) if path.exists() else {}
    criteria = {
        'complete_search': all(coverage),
        'validated_champions': evidence.get('converged_dft_count', 0) > 0 and
                               evidence.get('converged_orr_dft_count', 0) > 0,
        'validated_reactor': evidence.get('measured_reactor_count', 0) > 0 and
                             evidence.get('measured_deactivation_count', 0) > 0,
        'validated_pemfc': evidence.get('measured_mea_count', 0) > 0 and
                           evidence.get('measured_durability_count', 0) > 0 and
                           evidence.get('hydrogen_impurity_test_count', 0) > 0,
        'defensible_novelty': evidence.get('time_split_benchmark_count', 0) > 0 and
                              evidence.get('curated_prior_art_sources', 0) > 0,
    }
    if pyrolysis_mode == 'ntec':
        criteria['calibrated_ntec'] = evidence.get('ntec_control_pair_count', 0) > 0
    return {'ready': all(criteria.values()), 'criteria': criteria,
            'missing': [k for k, passed in criteria.items() if not passed]}
