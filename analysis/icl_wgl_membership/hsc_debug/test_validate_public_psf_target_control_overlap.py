from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from validate_public_psf_target_control_overlap import _independent_fit, _independent_profile


def test_independent_profile_and_scale_recover_known_factor():
    y, x = np.indices((43, 43), dtype=float)
    image = np.exp(-((x - 21) ** 2 + (y - 21) ** 2) / 100)
    edges = np.array([0.0, 1.0, 2.0, 3.0])
    profile, counts = _independent_profile(image, edges, 0.168)
    scaled, scaled_counts = _independent_profile(0.61 * image, edges, 0.168)
    assert np.array_equal(counts, scaled_counts)
    scale, residual, _ = _independent_fit(scaled, profile, counts, edges, (0.0, 3.0))
    assert scale == pytest.approx(0.61)
    assert residual < 1e-12


import pytest
