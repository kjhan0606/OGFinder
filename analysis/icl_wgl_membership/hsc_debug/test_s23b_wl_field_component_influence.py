import sys
from pathlib import Path
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from analyze_s23b_wl_field_component_influence import decompose


def row(field, lens, control):
    return {'field': field, 'lens_delta_sigma_msun_mpc2': str(lens),
            'equal_center_random_delta_sigma_msun_mpc2': str(control),
            'lens_cross_msun_mpc2': str(2*lens), 'equal_center_random_cross_msun_mpc2': str(2*control)}


def test_omission_separates_lens_and_subtracted_control():
    rows = [row('a', 10, 3), row('a', 6, 1), row('b', 4, 2)]
    value = decompose(rows, 'a', 'tangential')
    assert value['n_retained'] == 1
    assert value['full_profile'] == pytest.approx(14 / 3)
    assert value['lens_component_shift'] == pytest.approx(-8 / 3)
    assert value['minus_control_component_shift'] == pytest.approx(0)
    assert value['lens_only_omission_profile'] == pytest.approx(2)
    assert value['lens_only_omission_shift'] == pytest.approx(-8 / 3)
    cross = decompose(rows, 'a', 'cross')
    assert cross['lens_only_omission_shift'] == pytest.approx(2 * value['lens_only_omission_shift'])


def test_rejects_empty_or_nonfinite_group():
    with pytest.raises(ValueError, match='invalid'):
        decompose([row('a', 1, 2)], 'a', 'tangential')
    bad = row('b', 1, 2)
    bad['lens_delta_sigma_msun_mpc2'] = 'nan'
    with pytest.raises(ValueError, match='nonfinite'):
        decompose([row('a', 1, 2), bad], 'a', 'tangential')
