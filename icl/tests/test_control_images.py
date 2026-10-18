from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from validate_control_images import control_masks, nJy_factor, paired_row, ExactCosmology


def test_random_center_source_is_masked_and_science_image_is_unchanged():
    y, x = np.indices((256, 256))
    image = np.random.default_rng(22).normal(0, 1, (256, 256))
    for cx, cy in [(128, 128), (60, 60)]:
        image += 200*np.exp(-((x-cx)**2+(y-cy)**2)/8)
    image = image.astype('f4'); original = image.copy()
    var = np.ones_like(image); var[20, 30] = -1; var[40, 50] = np.nan
    bits = np.zeros(image.shape, 'i4'); bits[100, 100] = 512; bits[90, 90] = 1
    base, masks, seg, _, result = control_masks(image, var, bits, (128, 128), .168)
    assert result['center_segment'] > 0
    for mask in masks:
        assert mask[128, 128] and mask[60, 60]
        assert mask[20, 30] and mask[40, 50] and mask[100, 100] and mask[90, 90]
        assert np.all(mask[seg > 0]) and np.all(mask[base])
    assert np.all(~masks[0] | masks[1]) and np.all(~masks[1] | masks[2])
    np.testing.assert_array_equal(image, original)


def test_all_invalid_control_has_no_fabricated_detection():
    image = np.ones((64, 64), 'f4')
    with pytest.raises(ValueError, match='Too few valid pixels'):
        control_masks(image, -image, np.zeros(image.shape, 'i4'), (32, 32), .168)


def test_flux_units_and_missing_pair_statistic():
    info = dict(input_unit='instrumental_flux/pixel', fluxmag0=10**(27/2.5), pixel_scale_arcsec=.168)
    np.testing.assert_allclose(nJy_factor(info), 3.631e12*10**(-27/2.5)/.168**2)
    info['input_unit'] = 'nJy/pixel'
    np.testing.assert_allclose(nJy_factor(info), 1/.168**2)
    info['input_unit'] = 'counts/second'
    with pytest.raises(ValueError, match='Unknown input unit'):
        nJy_factor(info)
    control = dict(control_id='random1', median_tile_mean=None, status='insufficient_tiles',
                   n_tiles=0, outer_valid_fraction=0, plane_residual_mad_over_nominal=None)
    cluster = dict(cluster_id='1', cohort='area_qualified', kron_factor='2.5', median_tile_mean='-2',
                   n_tiles='12', status='diagnostic_only', outer_valid_fraction='.4', plane_residual_mad_over_nominal='2')
    row = paired_row(control, cluster, 3, 5, False)
    assert row['cluster_median_tile_mean_njy_arcsec2'] == -10
    assert row['cluster_minus_control_median_tile_mean_njy_arcsec2'] is None
    control.update(median_tile_mean=-3, status='diagnostic_only')
    assert paired_row(control, cluster, 3, 5, False)['cluster_minus_control_median_tile_mean_njy_arcsec2'] == -1


def test_radial_scale_matches_current_astropy_outer_diagnostics():
    from astropy.cosmology import FlatLambdaCDM
    from ogfkit.wlcluster import ARCSEC
    old_interface = ExactCosmology(70, .3)
    reference = FlatLambdaCDM(H0=70, Om0=.3, Tcmb0=0)
    for z in (.3, .42, .55):
        np.testing.assert_allclose(old_interface.angular_diameter(z)*ARCSEC*1000,
                                   reference.kpc_proper_per_arcmin(z).value/60, rtol=1e-12)
