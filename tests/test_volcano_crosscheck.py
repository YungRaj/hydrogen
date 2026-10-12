from pipeline.evidence.volcano_crosscheck import rank_descriptors, volcano_rate_log10


def test_volcano_is_maximal_near_optimum():
    assert volcano_rate_log10(0.86) > volcano_rate_log10(0.2)
    assert volcano_rate_log10(0.86) > volcano_rate_log10(1.5)


def test_rank_does_not_mutate_source_rows():
    rows = [{'candidate_id': 'a', 'dG_OH_eV': 0.86}, {'candidate_id': 'b', 'dG_OH_eV': 0.2}]
    ranked = rank_descriptors(rows)
    assert [row['candidate_id'] for row in ranked] == ['a', 'b']
    assert 'volcano_log10_rate' not in rows[0]
