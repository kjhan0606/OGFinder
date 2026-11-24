import copy
import importlib.util
from pathlib import Path
import sys

import pytest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
spec = importlib.util.spec_from_file_location('s16a_random_support', HERE/'audit_s23b_s16a_wl_random_support.py')
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


def fixture_counts(values=(0, 1000, 2000, 3000, 4000)):
    selected, counts = [], {}
    targets = {'L': dict(ra_deg=10., dec_deg=0., z_cl=.36)}
    for i, n in enumerate(values, 1):
        ra = 10+i*.1
        selected.append(dict(random_center_id=f'R{i}', matched_lens_id='L',
            target_s16a_raw_count='2000', z_lens='.36', ra_deg=str(ra), dec_deg='0',
            icl_filter='i', random_object_id=f'O{i}', control_slot=str(i),
            offset_from_target_mpc='7',
            offset_from_target_arcmin=str(audit.angular_separation_arcmin(10, 0, ra, 0)),
            match_score='.2', delta_depth_mag='.1', seeing_ratio='1.1', delta_ebv_sfd='.001',
            cluster_overlap_2mpc='0', tract='1', patch='2'))
        counts[f'R{i}'] = n
    return selected, counts, {'L': 2000}, targets


def test_zero_counts_stay_in_fixed_five_center_mean():
    inputs = fixture_counts()
    before = copy.deepcopy(inputs)
    rows = audit.join_centers(*inputs)
    lens = audit.lens_rows(rows)[0]
    assert lens['n_positive'] == lens['n_ge1000'] == 4
    assert lens['mean_random_count'] == 2000  # Not 2500 from dropping the zero.
    assert lens['mean_random_to_target'] == 1
    assert inputs == before
    summary = audit.summarize(rows, [lens])
    assert summary[0]['zero_centers'] == 1
    assert summary[0]['lenses_with_five_positive'] == 0
    assert summary[2]['n_lenses'] == 0
    assert summary[2]['center_ratio']['median'] is None


def test_all_zero_group_is_preserved_and_not_science_support():
    rows = audit.join_centers(*fixture_counts((0, 0, 0, 0, 0)))
    lens = audit.lens_rows(rows)[0]
    assert lens['n_positive'] == 0
    assert lens['mean_random_to_target'] == 0
    assert all(r['n_matching'] == 0 for r in audit.sensitivities([lens]))


@pytest.mark.parametrize('change', ['missing_count', 'extra_count', 'wrong_target', 'duplicate_slot', 'reused_object', 'wrong_separation', 'negative', 'noninteger'])
def test_invalid_join_is_refused(change):
    selected, counts, target, targets = fixture_counts()
    if change == 'missing_count': counts.pop('R1')
    elif change == 'extra_count': counts['other'] = 1
    elif change == 'wrong_target': target['L'] = 2001
    elif change == 'duplicate_slot': selected[1]['control_slot'] = '1'
    elif change == 'reused_object': selected[1]['random_object_id'] = 'O1'
    elif change == 'wrong_separation': selected[0]['offset_from_target_arcmin'] = '99'
    elif change == 'negative': counts['R1'] = -1
    elif change == 'noninteger': counts['R1'] = 1.5
    with pytest.raises(ValueError):
        audit.join_centers(selected, counts, target, targets)


def test_positive_does_not_imply_count_floor_or_ratio_balance():
    rows = audit.join_centers(*fixture_counts((500, 1000, 2000, 4000, 5000)))
    lens = audit.lens_rows(rows)[0]
    assert lens['n_positive'] == 5
    assert lens['n_ge1000'] == 4
    assert lens['n_ratio_0p5_to_2'] == 3  # Includes both diagnostic boundaries.
    result = {r['diagnostic']: r for r in audit.sensitivities([lens]) if r['icl_filter'] == 'all'}
    assert result['five_positive']['n_matching'] == 1
    assert result['five_ge1000']['n_matching'] == 0
    assert result['five_ratios_between_0p5_and_2']['n_matching'] == 0
    assert all(r['diagnostic_only'] for r in result.values())


def test_spherical_separation_wraps_right_ascension():
    assert audit.angular_separation_arcmin(359.9, 0, .1, 0) == pytest.approx(12)
    assert audit.angular_separation_arcmin(12, 80, 12, 80) == 0


def test_nonfinite_metadata_and_partial_group_are_refused():
    selected, counts, target, targets = fixture_counts()
    selected[0]['seeing_ratio'] = 'nan'
    with pytest.raises(ValueError, match='metadata'):
        audit.join_centers(selected, counts, target, targets)
    rows = audit.join_centers(*fixture_counts())
    with pytest.raises(ValueError, match='five-center'):
        audit.lens_rows(rows[:-1])
