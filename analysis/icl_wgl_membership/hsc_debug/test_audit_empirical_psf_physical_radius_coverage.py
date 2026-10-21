import pytest
from astropy.cosmology import FlatLambdaCDM

from audit_empirical_psf_physical_radius_coverage import (
    inscribed_stamp_radius_pixels,
    kpc_per_pixel,
)


def test_inscribed_radius_uses_nearest_stamp_edge():
    assert inscribed_stamp_radius_pixels(200.25, 199.75) == pytest.approx(199.75)
    assert inscribed_stamp_radius_pixels(5.0, 200.0) == pytest.approx(5.0)
    with pytest.raises(ValueError, match="outside"):
        inscribed_stamp_radius_pixels(-0.1, 200.0)


def test_proper_kpc_pixel_scale_matches_astropy_distance():
    cosmology = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0)
    scale = kpc_per_pixel(0.4, 0.168, cosmology)
    expected = cosmology.angular_diameter_distance(0.4).to_value("kpc") * 0.168 / 206264.80624709636
    assert scale == pytest.approx(expected, rel=1e-12)
    with pytest.raises(ValueError, match="positive"):
        kpc_per_pixel(0.0, 0.168, cosmology)
