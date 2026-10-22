import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).with_name("compare_icl_s23b_control_exclusion.py")
SPEC = importlib.util.spec_from_file_location("compare_icl_s23b_control_exclusion", SCRIPT)
compare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compare)


def test_area_weighted_mean_uses_annular_area_and_exact_bin_coverage():
    rows = [
        {"r_min_kpc": "50", "r_max_kpc": "100", "value": "2"},
        {"r_min_kpc": "100", "r_max_kpc": "200", "value": "4"},
    ]
    expected = (2 * (100**2 - 50**2) + 4 * (200**2 - 100**2)) / (200**2 - 50**2)
    assert compare.area_weighted_mean(rows, "value", 50, 200) == pytest.approx(expected)


def test_area_weighted_mean_rejects_incomplete_radial_bins():
    rows = [{"r_min_kpc": "50", "r_max_kpc": "100", "value": "2"}]
    with pytest.raises(ValueError, match="do not exactly tile"):
        compare.area_weighted_mean(rows, "value", 50, 200)
