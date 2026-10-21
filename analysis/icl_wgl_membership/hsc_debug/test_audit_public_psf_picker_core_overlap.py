from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_public_psf_picker_core_overlap import fit_profile_scale, radial_mean_profile


def test_overlap_scale_recovers_exact_amplitude_without_shape_residual():
    y, x = np.indices((43, 43), dtype=float)
    public = np.exp(-((x - 21) ** 2 + (y - 21) ** 2) / (2 * 8.0 ** 2))
    local = 2.4 * public[1:42, 1:42]
    edges = np.array([1.5, 1.8, 2.1, 2.4, 2.7, 3.0, 3.3])
    p, counts = radial_mean_profile(public, edges)
    l, local_counts = radial_mean_profile(local, edges)
    assert np.array_equal(counts, local_counts)
    fit = fit_profile_scale(l, p, counts)
    assert fit["amplitude_scale_local_over_public"] == pytest.approx(2.4)
    assert fit["weighted_fractional_shape_residual_rms"] < 1e-12


def test_overlap_scale_rejects_too_few_bins():
    with pytest.raises(ValueError):
        fit_profile_scale(np.array([1.0, np.nan]), np.array([1.0, 2.0]), np.array([4, 4]))


def test_radial_profile_rejects_even_stamp():
    with pytest.raises(ValueError):
        radial_mean_profile(np.ones((4, 5)))
