from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_public_psf_target_control_overlap import fit_window_scales, pair_by_actual_filter


def test_fit_windows_recover_constant_scale():
    edges = np.array([1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3])
    public = np.array([1.0, 0.8, 0.6, 0.4, 0.3, 0.2])
    local = 0.73 * public
    counts = np.array([20, 25, 30, 35, 40, 45])
    result = fit_window_scales(local, public, counts, edges)
    assert set(result) == {
        "inner_1p5_2p4_arcsec", "middle_1p8_2p7_arcsec",
        "outer_2p4_3p3_arcsec", "full_1p5_3p3_arcsec",
    }
    assert all(value == pytest.approx(0.73) for value in result.values())


def test_fit_windows_reject_inconsistent_radial_dimensions():
    with pytest.raises(ValueError, match="inconsistent"):
        fit_window_scales(np.ones(5), np.ones(6), np.ones(6), np.arange(7))


def test_pair_by_actual_filter_excludes_mismatched_controls():
    targets = {
        "a": {"filter": "i"},
        "b": {"filter": "i2"},
        "c": {"filter": "i"},
    }
    controls = {
        "a": {"filter": "i"},
        "b": {"filter": "i"},
        "c": {"filter": "i2"},
    }
    matched, mismatched = pair_by_actual_filter(targets, controls)
    assert [item[0] for item in matched] == ["a"]
    assert mismatched == [
        {"target_cluster_id": "b", "target_filter": "i2", "control_filter": "i"},
        {"target_cluster_id": "c", "target_filter": "i", "control_filter": "i2"},
    ]
