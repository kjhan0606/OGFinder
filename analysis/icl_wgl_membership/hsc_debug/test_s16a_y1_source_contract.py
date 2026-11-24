from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import s16a_y1_source_contract as c
import prepare_s23b_s16a_quality_supplement as p


def good():
    return dict(is_primary='t', **{k: 'f' for k in c.BAD_FLAGS},
                raw_sigma=.4, resolution=.3, extendedness=1, blendedness=.1,
                i_flux=100, i_flux_err=10, i_mag=24.7, a_i=.2,
                e1=.2, e2=-.3,
                **{f'{b}_{k}': v for b in c.BANDS for k, v in [('flux', 5), ('flux_err', 1)]})


def test_all_published_inclusive_boundaries_pass_without_science_promotion():
    row = good()
    result = c.select_row(row)
    assert result['table4_diagnostic'] == 'pass'
    assert result['missing_or_invalid_inputs'] == []
    assert result['science_measurement'] is False
    row.update(raw_sigma=0, e1=1.99, e2=0)
    assert c.select_row(row)['table4_diagnostic'] == 'pass'


@pytest.mark.parametrize('field,value,criterion', [
    ('raw_sigma', -.001, 'raw_sigma_0_to_0p4'),
    ('raw_sigma', .400001, 'raw_sigma_0_to_0p4'),
    ('resolution', .299999, 'resolution_ge_0p3'),
    ('e1', 2, 'distortion_lt_2'),
    ('blendedness', c.BLEND_MAX, 'blend_lt_limit'),
    ('i_flux', 99.999, 'i_snr_ge_10'),
    ('i_mag', 24.70001, 'dereddened_i_le_24p5'),
    ('extendedness', 0, 'extended'),
])
def test_rejects_beyond_limits_and_strict_upper_boundaries(field, value, criterion):
    row = good(); row.update(e2=0); row[field] = value
    assert c.select_row(row)['criteria'][criterion] == 'fail'


def test_raw_sigma_is_never_substituted_by_calibrated_noise_or_m_tail_clipped():
    row = good(); row.update(sigma_e=.63, m=-1.1)
    assert c.select_row(row)['table4_diagnostic'] == 'pass'
    del row['raw_sigma']
    assert c.select_row(row)['table4_diagnostic'] == 'unknown'
    assert 'raw_sigma' in c.select_row(row)['missing_or_invalid_inputs']
    row.update(raw_sigma=.5, sigma_e=.1)
    assert c.select_row(row)['table4_diagnostic'] == 'fail'


def test_boolean_strings_do_not_use_python_truthiness_and_null_never_passes():
    row = good()
    row['pixel_edge'] = 'false'
    assert c.select_row(row)['table4_diagnostic'] == 'pass'
    row['pixel_edge'] = 'true'
    assert c.select_row(row)['table4_diagnostic'] == 'fail'
    for invalid in ('null', '', 'yes', None, 0):
        row['pixel_edge'] = invalid
        assert c.select_row(row)['table4_diagnostic'] == 'unknown'


def test_two_of_four_detection_not_all_four_with_explicit_ambiguity():
    row = good(); row.update(z_flux=None, y_flux=0)
    assert c.select_row(row)['criteria']['two_other_bands_ge_5'] == 'pass'
    row.update(r_flux=0)
    assert c.select_row(row)['criteria']['two_other_bands_ge_5'] == 'unknown'
    row.update(z_flux=0)
    assert c.select_row(row)['criteria']['two_other_bands_ge_5'] == 'fail'


@pytest.mark.parametrize('field,value', [('i_flux_err', 0), ('i_flux_err', -1),
                                      ('resolution', 'nan'), ('e1', 'inf'),
                                      ('a_i', True), ('raw_sigma', None)])
def test_invalid_numeric_inputs_cannot_pass(field, value):
    row = good(); row[field] = value
    assert c.select_row(row)['table4_diagnostic'] == 'unknown'


def test_existing_pilot_columns_do_not_make_up_retained_counts():
    result = c.select_row(dict(e1=.1, e2=.2, sigma_e=.1, m=-.02))
    assert result['table4_diagnostic'] == 'unknown'
    assert len(result['missing_or_invalid_inputs']) == len(c.REQUIRED)-2


def test_exact_large_ids_and_pair_reuse_survive_join():
    pilot = [dict(center_id='L', object_id='9007199254740993', e1=.1),
             dict(center_id='R1', object_id='9007199254740993', e1=.1),
             dict(center_id='R1', object_id='9007199254740992', e1=.2)]
    supplement = [dict(object_id='9007199254740992', photoz_best=None),
                  dict(object_id='9007199254740993', photoz_best=.8)]
    merged = c.merge_supplement(pilot, supplement)
    assert len(merged) == 3 and [r['e1'] for r in merged] == [.1, .1, .2]
    assert [r['photoz_best'] for r in merged] == [.8, .8, None]
    assert 'photoz_best' not in pilot[0]
    for bad in (supplement[:1], supplement + supplement[:1],
                supplement + [dict(object_id='123')],
                [dict(object_id='9007199254740993', e1=0), supplement[0]]):
        with pytest.raises(ValueError): c.merge_supplement(pilot, bad)
    with pytest.raises(ValueError): c.merge_supplement(pilot+pilot[:1], supplement)


@pytest.mark.parametrize('value', [9007199254740992., '1.0', '01', '-1', True,
                                  str(2**63), '1); SELECT 1;--'])
def test_no_float_rounding_or_sql_injection_through_ids(value):
    with pytest.raises(ValueError): c.object_id(value)
    with pytest.raises(ValueError): p.query([value])


def test_quality_sql_has_only_fixed_ids_and_no_shape_reexport_or_lossy_join():
    sql = p.query(['9007199254740992', '9007199254740993'])
    assert sql.count('::bigint)') == 2
    assert 'WHERE' not in sql and 'LIMIT' not in sql and 'DISTINCT' not in sql
    assert sql.count('LEFT JOIN') == 4
    assert 'weaklensing_hsm_regauss' not in sql and 'coneSearch' not in sql
    assert 'derived_' not in sql and 'ira' not in sql and 'idec' not in sql
    assert 'm.idetect_is_primary AS is_primary' in sql
    assert 'a.ishape_hsm_regauss_sigma AS raw_sigma' in sql
    assert 'm.gcmodel_flux AS g_flux' in sql
    assert 'f.a_i AS a_i' in sql and 'z.photoz_best AS photoz_best' in sql
    assert len(c.EXPRESSIONS) == 31
    with pytest.raises(ValueError): p.query([])
    with pytest.raises(ValueError): p.query(['1', '1'])


def test_contract_separates_selection_photoz_and_measurement_gates():
    contract = c.contract()
    assert contract['photoz']['background_cut'] is None
    assert not contract['photoz']['calibrated_beta_available']
    assert not contract['science_measurement']
    assert len(contract['unresolved_measurement_gates']) == 6
