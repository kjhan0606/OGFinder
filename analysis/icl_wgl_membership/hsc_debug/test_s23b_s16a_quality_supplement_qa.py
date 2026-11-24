from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_s23b_s16a_quality_supplement as qa
from s16a_y1_source_contract import BAD_FLAGS


def centers():
    return [dict(center_id='L' if i == 0 else f'R{i}', matched_lens_id='A',
                 role='lens' if i == 0 else 'random', icl_filter='i', ra_deg=0., dec_deg=0.,
                 z_lens=.3, expected_raw_pairs=1, query_planned=True) for i in range(6)]


def row(ident, obj='9007199254740993'):
    return dict(center_id=ident, object_id=obj, source_ra='0', source_dec='0',
                is_primary='t', **{k:'f' for k in BAD_FLAGS}, raw_sigma='.2', resolution='.5',
                extendedness='1', blendedness='.1', i_mag='24', a_i='.1',
                i_flux='100', i_flux_err='10', e1='.2', e2='.1', sigma_e='.5', rms_e='.3',
                weight='2', m='-1.1', c1='0', c2='0', photoz_best='', photoz_risk_best='',
                **{f'{b}_{suffix}':value for b in 'grzy' for suffix,value in [('flux','5'),('flux_err','1')]})


def test_known_zero_control_is_retained_and_not_assigned_zero_shear():
    cs = centers();cs[-1].update(expected_raw_pairs=0, query_planned=False)
    result = qa.analyze(cs, [row(c['center_id']) for c in cs[:-1]])
    assert len(result['centers']) == 6 and len(result['radial']) == 48
    zero = result['centers'][-1]
    assert zero['raw_count'] == zero['selected_count'] == zero['effective_number_shape_weight'] == 0
    assert zero['source_origin'] == 'prior_count_only_zero_not_requeried'
    assert all(not c['raw_common_support'] and not c['selected_common_support'] for c in result['common'])
    assert all('shear' not in k for k in zero)


def test_reused_objects_and_nonexclusive_failures_are_counted_at_both_levels():
    cs = centers();cs[0]['expected_raw_pairs'] = 2
    records = [row(c['center_id']) for c in cs]
    failed = row('L','9007199254740992');failed.update(resolution='.2',pixel_edge='t')
    result = qa.analyze(cs, [*records,failed])
    assert result['unique_status'] == {'pass':1,'fail':1}
    assert result['pair_status'] == {'pass':6,'fail':1}
    criteria = {r['criterion']:r for r in result['criteria']}
    assert criteria['resolution_ge_0p3']['unique_fail'] == 1
    assert criteria['pixel_edge']['unique_fail'] == 1
    assert sum(r['selected_count'] for r in result['radial']) == 6
    assert result['common'][0]['selected_common_support']
    assert len(result['selection_rows']) == len(result['tails']) == 2
    # Diagnostic m/derived-sigma tails do not override the raw-sigma object cut.
    assert sum(r['table4_diagnostic']=='pass' for r in result['tails']) == 1


def test_missing_point_photoz_does_not_apply_an_unapproved_background_cut():
    cs = centers(); records = [row(c['center_id']) for c in cs]
    result = qa.analyze(cs, records)
    assert result['pair_status'] == {'pass':6}
    assert result['photoz']['photoz_best_all_unique']['missing_or_nonfinite'] == 1
    assert result['photoz']['photoz_best_table4_pass_unique']['total_objects'] == 1


def test_unknown_quality_is_preserved_without_fabricating_support():
    cs = centers(); records = [row(c['center_id']) for c in cs]
    for r in records:r['raw_sigma']=''
    result = qa.analyze(cs, records)
    assert result['pair_status'] == {'unknown':6}
    assert result['common'][0]['raw_common_support']
    assert not any(r['selected_common_support'] for r in result['common'])


def test_source_inconsistency_duplicate_pair_and_raw_count_mismatch_abort():
    cs = centers(); records = [row(c['center_id']) for c in cs]
    with pytest.raises(ValueError, match='duplicate'):
        qa.analyze(cs, records+records[:1])
    records[-1]['e1']='.3'
    with pytest.raises(ValueError, match='properties differ'):
        qa.analyze(cs, records)
    with pytest.raises(ValueError, match='raw count'):
        qa.analyze(cs, records[:-1])
