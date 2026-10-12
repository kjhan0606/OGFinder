"""Scaled ICL stack: apparent size and estimated mass set the common radius."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit.wlcluster import ARCSEC, Cosmology, rho_crit  # noqa: E402
from icl.wgl import (  # noqa: E402
    area_corrected_residual, forward_scaled_profile, gaussian_ln_mass_nodes,
    mass_interpretation, mean_stack, named_light_fraction, r_delta_mpc,
    radial_sample, radius_log_shift, scaled_radius, stack_clusters,
    survey_checklist, theta_delta_arcsec, warp_metric, warp_surface_brightness)


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


def _paint_kpc(z, cosmo, kpc_max, s0, scale_kpc, half, cover=1.5):
    distance = float(np.asarray(cosmo.angular_diameter(z)))
    kpc_per_arcsec = distance * ARCSEC * 1000.0
    kpc_per_pixel = (kpc_max * cover) / half
    pixel_scale = kpc_per_pixel / kpc_per_arcsec
    n = 2 * half + 1
    center = ((n - 1) / 2.0, (n - 1) / 2.0)
    yy, xx = np.mgrid[0:n, 0:n]
    radius = np.hypot(xx - center[0], yy - center[1]) * kpc_per_pixel
    sb = s0 * np.exp(-radius / scale_kpc)
    return sb * pixel_scale ** 2, center, pixel_scale


def test_common_kpc_grid_recovers_one_profile_at_two_redshifts():
    cosmo = Cosmology()
    kpc_max, s0, scale_kpc = 80.0, 2.5, 30.0
    n_side = 41
    records = []
    frames = []
    for z, mass in ((0.25, 1.2e15), (0.62, 3.0e14)):
        image, center, pixel_scale = _paint_kpc(z, cosmo, kpc_max, s0, scale_kpc, half=48)
        frame = warp_metric(image, center, pixel_scale, z, cosmo, kpc_max=kpc_max, n_side=n_side)
        assert np.all(frame['weight'] == 1)
        frames.append(frame)
        records.append(dict(image=image, center=center, pixel_scale=pixel_scale, mass=mass, z=z, id='z%s' % z))
    averaged = mean_stack(frames)
    ky, kx = np.meshgrid(averaged['kpc'], averaged['kpc'], indexing='ij')
    model = s0 * np.exp(-np.hypot(kx, ky) / scale_kpc)
    assert averaged['sb'] == pytest.approx(model, rel=0.02, abs=0.02)
    assert 'u' not in averaged
    stacked = stack_clusters(records, cosmo, kpc_max=kpc_max, u_max=0.4, n_side=n_side)
    assert stacked['product'] == 'bcg+diffuse'
    assert 'BCG+diffuse' in stacked['description']
    assert 'counts per arcsec^2' in stacked['surface_brightness_note']
    assert stacked['kpc']['product'] == 'bcg+diffuse'
    assert stacked['kpc']['kpc_max_source'] == 'caller'
    assert stacked['kpc']['sb'] == pytest.approx(model, rel=0.02, abs=0.02)
    assert stacked['survey']['science_measurement'] is False
    assert stacked['interpretation']['independent_evidence'] is False
    assert 'M_200c and M_200m' in stacked['interpretation']['definition_note']


def test_mass_error_shifts_the_scaled_radius_by_one_third():
    assert radius_log_shift(3.0 * np.log(2.0)) == pytest.approx(np.log(2.0))
    u = np.linspace(0.0, 4.0, 2001)
    sb = np.exp(-u)
    one = forward_scaled_profile(u, sb, 3.0 * np.log(2.0))
    assert one == pytest.approx(np.exp(-0.5 * u), rel=1e-4, abs=1e-4)
    assert forward_scaled_profile(u, sb, 0.0) == pytest.approx(sb)
    shifts = np.array([3.0 * np.log(2.0), -3.0 * np.log(2.0)])
    weights = np.array([1.0, 3.0])
    averaged = forward_scaled_profile(u, sb, shifts, weights)
    manual = (
        np.interp(0.5 * u, u, sb, left=np.nan, right=np.nan)
        + 3.0 * np.interp(2.0 * u, u, sb, left=np.nan, right=np.nan)) / 4.0
    assert averaged == pytest.approx(manual, nan_ok=True)
    nodes, node_weights = gaussian_ln_mass_nodes(0.3, n=7)
    assert node_weights.sum() == pytest.approx(1.0)
    assert nodes == pytest.approx(-nodes[::-1])


def test_outer_radial_bin_drops_a_cluster_that_ends_inside():
    axis = np.linspace(-1.0, 1.0, 21)
    yy, xx = np.meshgrid(axis, axis, indexing='ij')
    radius = np.hypot(xx, yy)
    wide = {'id': 'A', 'u': axis, 'sb': np.ones((21, 21)), 'weight': np.ones((21, 21))}
    inner = {'id': 'B', 'u': axis, 'sb': np.full((21, 21), 3.0), 'weight': (radius <= 0.3).astype(float)}
    table = radial_sample([wide, inner], n_bins=4)
    assert table['bins'][0]['ids'] == ['A', 'B']
    assert table['bins'][0]['n_clusters'] == 2
    assert table['bins'][-1]['ids'] == ['A']
    assert table['sample_changes_with_radius'] is True
    assert 'outermost' in table['note']
    stacked = mean_stack([wide, inner])
    assert stacked['radial']['sample_changes_with_radius'] is True
    assert stacked['ids'] == ['A', 'B']
    assert 'inverse_variance_sb' not in stacked


def test_inverse_variance_mean_stays_off_the_default_mean():
    axis = np.linspace(-1.0, 1.0, 5)
    shape = (5, 5)
    faint = {
        'id': 'faint', 'u': axis, 'sb': np.full(shape, 2.0),
        'weight': np.ones(shape), 'variance': np.ones(shape)}
    bright = {
        'id': 'bright', 'u': axis, 'sb': np.full(shape, 8.0),
        'weight': np.ones(shape), 'variance': np.full(shape, 3.0)}
    stacked = mean_stack([faint, bright])
    assert stacked['sb'] == pytest.approx(5.0)
    assert stacked['inverse_variance_sb'] == pytest.approx(3.5)
    assert stacked['inverse_variance_sb'] != pytest.approx(stacked['sb'])
    missing = dict(bright)
    del missing['variance']
    with pytest.raises(ValueError):
        mean_stack([faint, missing])


def test_open_area_correction_keeps_the_unmasked_extra_light():
    diffuse = np.full((8, 8), 1.0)
    observed = diffuse.copy()
    observed[0, :] += 0.3
    observed[1, 1] += 0.7
    observed[1, 2] += 0.7
    observed[4, 4] += 5.0
    mask = np.zeros((8, 8), bool)
    mask[6:8, 6:8] = True
    got = area_corrected_residual(diffuse, observed, mask)
    assert got['residual'] == pytest.approx(got['extra_light_mean'])
    assert got['residual'] == pytest.approx(8.8 / 60.0)
    assert got['masked_diffuse_sum'] == pytest.approx(4.0)
    assert got['measured'] == pytest.approx(got['diffuse_on_open_pixels'] + got['residual'])
    assert 'open area' in got['note']


def test_survey_gate_stays_down_until_the_checklist_is_filled():
    required = (
        'luminosity_unit', 'band', 'extinction', 'sky_residual', 'enlarged_mask',
        'psf', 'radial_sample', 'covariance')
    empty = survey_checklist()
    assert empty['science_measurement'] is False
    assert 'rest_frame' in empty['empty']
    for key in required:
        assert key in empty['empty']
    blank = {key: '  ' for key in required}
    assert survey_checklist(blank, rest_frame='observed')['science_measurement'] is False
    filled = {key: 'recorded' for key in required}
    assert survey_checklist(filled, rest_frame='observed')['science_measurement'] is True
    rest = survey_checklist(filled, rest_frame='rest')
    assert rest['science_measurement'] is False
    assert 'k_correction' in rest['empty']
    assert 'surface_brightness_dimming' in rest['empty']
    filled['k_correction'] = 'applied'
    filled['surface_brightness_dimming'] = 'applied'
    done = survey_checklist(filled, rest_frame='rest')
    assert done['science_measurement'] is True
    assert done['items']['rest_frame'] == 'rest'


def test_light_fraction_names_the_callers_definition():
    included = named_light_fraction(
        10.0, 40.0, definition='centre_exclude', centre_exclude_radius=50.0,
        numerator_includes_bcg=True)
    assert included['product'] == 'bcg+icl'
    assert 'BCG' in included['description']
    assert included['ratio'] == pytest.approx(0.25)
    assert included['definition_parameters']['centre_exclude_radius'] == pytest.approx(50.0)
    diffuse = named_light_fraction(definition='surface_brightness', sb_limit=26.5)
    assert diffuse['product'] == 'diffuse'
    assert diffuse['ratio'] is None
    assert diffuse['definition_parameters']['sb_limit'] == pytest.approx(26.5)
    decomposed = named_light_fraction(definition='model_decomposition', decomposition='double Sersic')
    assert decomposed['definition_parameters']['decomposition'] == 'double Sersic'
    with pytest.raises(ValueError):
        named_light_fraction(definition='centre_exclude')
    with pytest.raises(ValueError):
        named_light_fraction(definition='surface_brightness')
    with pytest.raises(ValueError):
        named_light_fraction(definition='model_decomposition')
    with pytest.raises(ValueError):
        named_light_fraction(definition='icl')


def test_same_lens_mass_is_not_recorded_as_independent_evidence():
    cov = np.eye(2)
    u = np.linspace(0.0, 2.0, 21)
    model = np.exp(-u)
    same = mass_interpretation(
        mass_origin='same_lens', delta_ln_mass=0.0, lensing_covariance=cov,
        u=u, model_sb=model, stack_sb=model)
    assert same['independent_evidence'] is False
    assert same['covariance_propagated'] is True
    assert same['smeared_model'] == pytest.approx(model)
    assert same['residual'] == pytest.approx(np.zeros_like(model))
    assert any('same result' in note for note in same['notes'])
    richness = mass_interpretation(
        mass_origin='richness', splits_sample=True, lensing_covariance=cov)
    assert richness['independent_evidence'] is False
    assert richness['covariance_propagated'] is False
    assert any('denominator' in note for note in richness['notes'])
    assert any('Eddington' in note for note in richness['notes'])
    external = mass_interpretation(mass_origin='external', mass_definition='M500c')
    assert any('Record the conversion' in note for note in external['notes'])
    converted = mass_interpretation(
        mass_origin='external', mass_definition='M200m', conversion='NFW c=4, M200m to M200c')
    assert converted['conversion'] == 'NFW c=4, M200m to M200c'
    assert converted['independent_evidence'] is False
    assert not any('Record the conversion' in note for note in converted['notes'])


def test_warped_variance_is_counts_variance_per_area_squared():
    cosmo = Cosmology()
    mass, z = 1.0e15, 0.3
    theta = theta_delta_arcsec(mass, z, cosmo)
    pixel_scale = theta / 20.0
    image = np.full((41, 41), 4.0)
    variance = np.full((41, 41), 9.0)
    frame = warp_surface_brightness(
        image, (20, 20), pixel_scale, mass, z, cosmo, u_max=0.5, n_side=11, variance=variance)
    assert np.all(frame['weight'] == 1)
    assert frame['sb'] == pytest.approx(4.0 / pixel_scale ** 2)
    assert frame['variance'] == pytest.approx(9.0 / pixel_scale ** 4)


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
