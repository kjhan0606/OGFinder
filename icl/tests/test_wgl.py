"""Scaled ICL stack: apparent size and estimated mass set the common radius."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit.wlcluster import Cosmology, rho_crit  # noqa: E402
from icl.wgl import (  # noqa: E402
    mean_stack, r_delta_mpc, scaled_radius, stack_clusters, theta_delta_arcsec,
    warp_surface_brightness)


def test_mass_radius_closes_at_two_redshifts():
    cosmo = Cosmology()
    for z, mass in ((0.3, 1.0e15), (0.8, 2.5e14)):
        radius = r_delta_mpc(mass, z, cosmo, delta=200.0)
        rho = rho_crit(z, cosmo.H0, cosmo.Om)
        enclosed = 4.0 * np.pi / 3.0 * 200.0 * rho * radius ** 3
        assert enclosed == pytest.approx(mass, rel=1e-12)
    assert r_delta_mpc(8.0e15, 0.3, cosmo) == pytest.approx(2.0 * r_delta_mpc(1.0e15, 0.3, cosmo), rel=1e-12)


def test_apparent_angle_is_the_mass_radius_over_distance():
    cosmo = Cosmology()
    z, mass = 0.4, 8.0e14
    distance = float(np.asarray(cosmo.angular_diameter(z)))
    theta = theta_delta_arcsec(mass, z, cosmo)
    proper = theta * (np.pi / 648000.0) * distance
    assert proper == pytest.approx(r_delta_mpc(mass, z, cosmo), rel=1e-12)


def test_scaled_radius_follows_apparent_size_and_mass():
    cosmo = Cosmology()
    z = 0.35
    theta = 120.0
    base = float(scaled_radius(theta, 1.0e15, z, cosmo))
    assert float(scaled_radius(2.0 * theta, 1.0e15, z, cosmo)) == pytest.approx(2.0 * base, rel=1e-12)
    assert float(scaled_radius(theta, 8.0e15, z, cosmo)) == pytest.approx(0.5 * base, rel=1e-12)
    assert float(scaled_radius(theta_delta_arcsec(1.0e15, z, cosmo), 1.0e15, z, cosmo)) == pytest.approx(1.0)


def _paint(z, mass, pixels_per_radius, n, s0, u_scale, cosmo):
    theta = theta_delta_arcsec(mass, z, cosmo)
    pixel_scale = theta / pixels_per_radius
    center = ((n - 1) / 2.0, (n - 1) / 2.0)
    yy, xx = np.mgrid[0:n, 0:n]
    apparent = np.hypot(xx - center[0], yy - center[1]) * pixel_scale
    u = apparent / theta
    sb = s0 * np.exp(-u / u_scale)
    return sb * pixel_scale ** 2, center, pixel_scale


def test_constant_surface_brightness_does_not_depend_on_the_mass_scale():
    cosmo = Cosmology()
    mass, z = 1.0e15, 0.3
    theta = theta_delta_arcsec(mass, z, cosmo)
    pixel_scale = theta / 20.0
    image = np.full((81, 81), 2.0)
    frame = warp_surface_brightness(image, (40, 40), pixel_scale, mass, z, cosmo, u_max=0.8, n_side=17)
    assert np.all(frame['weight'] == 1)
    assert frame['sb'] == pytest.approx(2.0 / pixel_scale ** 2)


def test_two_clusters_recover_one_scaled_profile_and_a_mask_is_excluded():
    cosmo = Cosmology()
    s0, u_scale = 2.5, 0.45
    image_a, center_a, scale_a = _paint(0.25, 1.2e15, 30, 91, s0, u_scale, cosmo)
    image_b, center_b, scale_b = _paint(0.62, 3.0e14, 24, 91, s0, u_scale, cosmo)
    common = dict(cosmo=cosmo, u_max=1.1, n_side=41, delta=200.0)
    clean = warp_surface_brightness(image_a, center_a, scale_a, 1.2e15, 0.25, **common)
    other = warp_surface_brightness(image_b, center_b, scale_b, 3.0e14, 0.62, **common)
    uy, ux = np.meshgrid(clean['u'], clean['u'], indexing='ij')
    model = s0 * np.exp(-np.hypot(ux, uy) / u_scale)
    assert clean['weight'].all() and other['weight'].all()
    assert clean['sb'] == pytest.approx(model, rel=0.01, abs=0.01)
    assert other['sb'] == pytest.approx(model, rel=0.01, abs=0.01)

    spiked = image_a.copy()
    mask = np.zeros(image_a.shape, bool)
    cx, cy = int(round(center_a[0])), int(round(center_a[1]))
    spiked[cy - 2:cy + 3, cx + 8:cx + 13] = 1.0e9
    mask[cy - 2:cy + 3, cx + 8:cx + 13] = True
    hidden = warp_surface_brightness(spiked, center_a, scale_a, 1.2e15, 0.25, mask=mask, **common)
    assert np.any(hidden['weight'] == 0)
    stacked = mean_stack([hidden, other])
    assert stacked['sb'][hidden['weight'] == 0] == pytest.approx(other['sb'][hidden['weight'] == 0], rel=1e-12, abs=1e-12)
    assert np.all(stacked['weight'][hidden['weight'] == 0] == 1)
    both = mean_stack([clean, other])
    assert both['sb'] == pytest.approx(model, rel=0.01, abs=0.01)
    assert np.all(both['weight'] == 2)

    leaked = warp_surface_brightness(spiked, center_a, scale_a, 1.2e15, 0.25, **common)
    assert np.nanmax(leaked['sb']) > 10.0 * s0


def test_one_cluster_stack_matches_its_own_warp():
    cosmo = Cosmology()
    image, center, pixel_scale = _paint(0.3, 1.0e15, 28, 81, 1.5, 0.5, cosmo)
    record = dict(image=image, center=center, pixel_scale=pixel_scale, mass=1.0e15, z=0.3)
    stacked = stack_clusters([record], cosmo=cosmo, u_max=1.0, n_side=31)
    alone = warp_surface_brightness(image, center, pixel_scale, 1.0e15, 0.3, cosmo, u_max=1.0, n_side=31)
    assert stacked['n'] == 1
    assert stacked['sb'] == pytest.approx(alone['sb'], nan_ok=True)
    assert np.array_equal(stacked['weight'], alone['weight'])


def test_rejects_a_non_positive_distance_mass_or_scale():
    cosmo = Cosmology()
    image = np.ones((9, 9))
    with pytest.raises(ValueError):
        theta_delta_arcsec(1.0e15, 0.0, cosmo)
    with pytest.raises(ValueError):
        r_delta_mpc(-1.0, 0.3, cosmo)
    with pytest.raises(ValueError):
        warp_surface_brightness(image, (4, 4), 0.0, 1.0e15, 0.3, cosmo, n_side=5, u_max=0.2)
    with pytest.raises(ValueError):
        warp_surface_brightness(image, (4, 4), 0.1, 1.0e15, 0.3, cosmo, mask=np.zeros((4, 4), bool), n_side=5, u_max=0.2)
    with pytest.raises(ValueError):
        mean_stack([])
