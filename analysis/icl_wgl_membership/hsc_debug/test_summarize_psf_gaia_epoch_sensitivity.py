import math

from summarize_psf_gaia_epoch_sensitivity import (
    propagate_tangent_motion,
    separation_arcsec,
)


def test_zero_motion_keeps_sky_position():
    moved = propagate_tangent_motion(123.4, -22.5, 0.0, 0.0, 10.0)
    assert moved == (123.4, -22.5)


def test_ra_motion_uses_cos_dec_tangent_component_and_wraps_ra():
    ra, dec = propagate_tangent_motion(359.9999, 0.0, 1000.0, 0.0, 1.0)
    assert 0.0 <= ra < 0.001
    assert math.isclose(dec, 0.0, abs_tol=1e-12)
    assert math.isclose(separation_arcsec(359.9999, 0.0, ra, dec), 1.0, abs_tol=1e-7)


def test_dec_motion_and_negative_epoch_reverse_on_sphere():
    ra_plus, dec_plus = propagate_tangent_motion(17.0, 63.0, 0.0, 2500.0, 2.0)
    ra_minus, dec_minus = propagate_tangent_motion(17.0, 63.0, 0.0, 2500.0, -2.0)
    assert math.isclose(separation_arcsec(17.0, 63.0, ra_plus, dec_plus), 5.0, abs_tol=1e-7)
    assert math.isclose(separation_arcsec(17.0, 63.0, ra_minus, dec_minus), 5.0, abs_tol=1e-7)
    assert math.isclose(dec_plus, 63.0 + 5.0 / 3600.0, abs_tol=1e-8)
    assert math.isclose(dec_minus, 63.0 - 5.0 / 3600.0, abs_tol=1e-8)
