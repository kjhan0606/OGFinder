import numpy as np
import pytest

from audit_wl_published_nz_shift_sensitivity import (
    PUBLISHED_HIGH_Z_SHIFT_INTERVALS,
    PUBLISHED_SHIFT_SCENARIOS,
    translate_nz_curve,
)


def test_published_shift_scenarios_cover_central_and_interval_corners():
    assert len(PUBLISHED_SHIFT_SCENARIOS) == 5
    assert PUBLISHED_SHIFT_SCENARIOS["central"][3] == pytest.approx(0.043)
    assert PUBLISHED_SHIFT_SCENARIOS["central"][4] == pytest.approx(0.050)
    assert PUBLISHED_HIGH_Z_SHIFT_INTERVALS == {
        3: (0.021, 0.063),
        4: (0.038, 0.062),
    }
    corners = [tuple(value[b] for b in (3, 4))
               for name, value in PUBLISHED_SHIFT_SCENARIOS.items()
               if name != "central"]
    assert len(corners) == 4
    for expected in ((0.021, 0.038), (0.021, 0.062),
                     (0.063, 0.038), (0.063, 0.062)):
        assert any(actual == pytest.approx(expected) for actual in corners)


def test_translation_moves_distribution_and_renormalizes():
    redshift = np.linspace(0.0, 5.0, 501)
    nz = np.exp(-0.5 * ((redshift - 1.1) / 0.25) ** 2)
    shifted, details = translate_nz_curve(redshift, nz, 0.043)
    assert np.trapezoid(shifted, redshift) == pytest.approx(1.0, abs=1e-12)
    assert details["realized_mean_shift_after_grid_truncation"] == pytest.approx(0.043, abs=1e-6)
    assert details["mean_redshift_after"] > details["mean_redshift_before"]


def test_translation_rejects_nonpositive_result():
    redshift = np.linspace(0.0, 1.0, 101)
    nz = np.exp(-0.5 * ((redshift - 0.1) / 0.03) ** 2)
    with pytest.raises(ValueError, match="no positive finite integral"):
        translate_nz_curve(redshift, nz, 2.0)
