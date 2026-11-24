import numpy as np
import pytest

from ogfkit.wl_y1 import rotate_sky_components, source_frame_angle, y1_bin_delta_sigma, y1_shape_response


def test_source_frame_wrap_and_tangential_sign():
    phi = source_frame_angle(.001, 0, 359.999, 0)
    assert rotate_sky_components(-.2, 0., phi) == pytest.approx([.2, 0.])
    phi = source_frame_angle(22, 6.1, 22, 6.)
    assert rotate_sky_components(.2, 0., phi) == pytest.approx([.2, 0.])
    with pytest.raises(ValueError, match='undefined'):
        source_frame_angle(0., 0., 360., 0.)


def test_known_delta_sigma_with_variable_weights_and_additive_calibration():
    beta = np.array([1e-15, 2e-15, 0.])
    weights = np.array([1., 3., 400.])
    lens = np.array([2., 1., 1.])
    rms = np.array([.3, .4, np.nan])
    m = np.array([-.97, -.1, np.nan])
    R = 1-(2*.3**2+12*.4**2)/14
    mean_m = (2*-.97+12*-.1)/14
    expected = [2e13, -3e12]
    additive = np.array([[.002, .003], [-.004, .001], [np.nan, np.nan]])
    e = 2*R*(beta[:, None]*np.asarray(expected)*(1+mean_m)*1.02*.92 + additive)
    result = y1_bin_delta_sigma(e, additive, weights, beta, rms, m,
                               lens_weight=lens, selection_m=.02, photoz_response=.92)
    assert result['n_pairs'] == 2
    assert result['responsivity'] == pytest.approx(R)
    assert result['mean_m'] == pytest.approx(mean_m)
    assert [result['tangential'], result['cross']] == pytest.approx(expected)


def test_ensemble_handles_individual_low_m_without_dividing_by_it():
    assert y1_shape_response([1, 4], [.3, .4], [-1.1, -.1])['mean_m'] == pytest.approx(-.3)
    with pytest.raises(ValueError, match='ensemble calibration'):
        y1_shape_response([1, 0], [.3, .4], [-1.1, -.1])
