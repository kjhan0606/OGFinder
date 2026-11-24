import copy
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prepare_s23b_s16a_shape_qa_pilot as pilot


def test_pilot_is_deterministic_and_keeps_each_support_stratum():
    rows = []
    for band in ('i', 'i2'):
        for positive, size in ((5, 8), (3, 4), (0, 1)):
            for k in range(size):
                rows.append(dict(matched_lens_id=f'{band}-{positive}-{k}', icl_filter=band,
                                 target_raw_count=1000+k, n_positive=positive))
    original = copy.deepcopy(rows)
    result = pilot.pilot_lenses(rows)
    assert result == pilot.pilot_lenses(list(reversed(rows)))
    assert rows == original
    assert len(result) == 10  # Two singleton strata contribute one each.
    five_i = [r for r in result if r['pilot_stratum'] == 'i:five_positive']
    assert [r['stratum_rank_zero_based'] for r in five_i] == [2, 6]
    assert {r['pilot_stratum'] for r in result} == {
        f'{b}:{c}' for b in ('i', 'i2') for c in ('five_positive', 'partial_positive', 'all_zero')}
    with pytest.raises(ValueError, match='duplicate'):
        pilot.pilot_lenses(rows+[rows[0]])


def center_fixture():
    lenses = [dict(matched_lens_id='S23B_1', icl_filter='i', z_lens=.36,
                   target_raw_count=2000, lens_pdr3_field='SPRING')]
    targets = {'S23B_1': dict(ra_deg=20, dec_deg=0, z_cl=.36)}
    randoms = [dict(matched_lens_id='S23B_1', random_center_id=f'R{k}', control_slot=k,
                    ra_deg=20+k*.1, dec_deg=0, target_raw_count=2000,
                    random_raw_count=0 if k < 3 else 2000) for k in range(1, 6)]
    return lenses, randoms, targets


def test_zero_centers_preserved_in_ledger_but_not_queried():
    result = pilot.make_centers(*center_fixture())
    assert len(result) == 6
    assert sum(r['query_planned'] for r in result) == 4
    assert sum(r['expected_raw_pairs'] for r in result) == 8000
    assert result[0]['role'] == 'lens'
    assert result[0]['center_id'] == 'LC-S23B_1'
    assert result[1]['expected_raw_pairs'] == result[2]['expected_raw_pairs'] == 0
    assert result[1]['query_planned'] is result[2]['query_planned'] is False


@pytest.mark.parametrize('case', ['missing_slot', 'duplicate_position', 'redshift', 'target_count'])
def test_invalid_pilot_centers_are_refused(case):
    lenses, randoms, targets = center_fixture()
    if case == 'missing_slot': randoms[0]['control_slot'] = 2
    elif case == 'duplicate_position': randoms[0]['ra_deg'] = targets['S23B_1']['ra_deg']
    elif case == 'redshift': lenses[0]['z_lens'] = .4
    elif case == 'target_count': randoms[0]['target_raw_count'] = 2001
    with pytest.raises(ValueError):
        pilot.make_centers(lenses, randoms, targets)


def test_sql_is_limited_to_documented_twelve_columns_and_same_cone_join():
    center = pilot.make_centers(*center_fixture())[0]
    sql = pilot.query(center)
    assert len(pilot.OUTPUT_COLUMNS) == 12
    assert sql.startswith("SELECT DISTINCT 'LC-S23B_1' AS center_id,")
    assert 's16a_wide.meas2 AS a' in sql
    assert 's16a_wide.weaklensing_hsm_regauss AS b' in sql
    assert 's16a_wide.meas AS c' in sql
    assert sql.count(' AS ') == 15  # Twelve output aliases plus three tables.
    assert 'coneSearch(c.icoord,' in sql
    for unwanted in ('s19a', 'photoz', 'zbin', '*', 'psf_shape'):
        assert unwanted not in sql
    for expr in pilot.COLUMN_EXPRESSIONS.values():
        assert expr in sql
