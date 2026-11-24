import sys
from pathlib import Path
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent))
from analyze_s23b_fixed52_tract_region_influence import group_exact_regions_by_field


def test_groups_parent_tracts_and_preserves_unassigned():
    exact = {'lens_center_lookup': {'l1': 0, 'l2': 2}, 'control_center_lookup': {'c': 1},
             'lens_center': np.array([0, 2, -1]), 'lens_source': np.array([1, 0, -1]),
             'control_center': np.array([1, 2]), 'control_source': np.array([2, 0])}
    audit = {'region_labels': {'0': 'SPRING:tract-1', '1': 'SPRING:tract-2', '2': 'AUTUMN:tract-3'}}
    result, details = group_exact_regions_by_field(exact, audit)
    assert details['region_labels'] == {'0': 'AUTUMN', '1': 'SPRING'}
    assert result['lens_center'].tolist() == [1, 0, -1]
    assert result['lens_source'].tolist() == [1, 1, -1]
    assert result['control_center_lookup'] == {'c': 1}
    assert details['tract_region_labels'] == audit['region_labels']
    assert exact['lens_center'].tolist() == [0, 2, -1]


def test_rejects_unknown_or_malformed_tract_codes():
    exact = {'lens_center_lookup': {'l': 9}, 'control_center_lookup': {},
             'lens_center': np.array([9]), 'lens_source': np.array([]),
             'control_center': np.array([]), 'control_source': np.array([])}
    with pytest.raises(ValueError, match='unknown'):
        group_exact_regions_by_field(exact, {'region_labels': {'0': 'SPRING:tract-1'}})
    with pytest.raises(ValueError, match='invalid'):
        group_exact_regions_by_field(exact, {'region_labels': {'0': 'SPRING'}})
