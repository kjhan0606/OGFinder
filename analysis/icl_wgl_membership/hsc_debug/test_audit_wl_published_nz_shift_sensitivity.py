import numpy as np
import pytest

from audit_wl_published_nz_shift_sensitivity import translate_nz_curve


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
