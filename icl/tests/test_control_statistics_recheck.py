from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from recheck_control_statistics import fit_plane, lag_statistic, unique


def test_independent_weighted_plane_recovers_known_gradient():
    rows = []
    for y in (-400, -100, 200, 500):
        for x in (-500, -250, 0, 250):
            rows.append(dict(x_kpc=x, y_kpc=y, raw_mean=-.3+.7*x/1000-.2*y/1000,
                             nominal_mean_sigma=.1+(x+500)/1000))
    fit, residual = fit_plane(rows)
    np.testing.assert_allclose([fit[k] for k in ('plane_constant', 'plane_x_per_mpc', 'plane_y_per_mpc')],
                               [-.3, .7, -.2], atol=1e-12)
    assert max(abs(residual)) < 1e-12
    for row in rows:
        row['y_kpc'] = 0
    with pytest.raises(ValueError, match='independent tile positions'):
        fit_plane(rows)


def test_pairs_match_explicit_coordinates_without_crossing_tile_edges():
    image = np.random.default_rng(17).normal(size=(9, 11))
    valid = np.random.default_rng(22).uniform(size=image.shape) > .2
    for axis in (0, 1):
        for lag in (1, 2, 3):
            pairs = []
            for y in range(image.shape[0]):
                for x in range(image.shape[1]):
                    other = [y, x]; other[axis] += lag
                    if (other[axis] >= image.shape[axis] or (y, x)[axis]//4 != other[axis]//4
                            or not valid[y, x] or not valid[tuple(other)]):
                        continue
                    pairs.append([image[tuple(other)], image[y, x]])
            pairs = np.array(pairs)
            result = lag_statistic(image, valid, lag, axis, tile_size=4)
            assert result['n_pairs'] == len(pairs)
            np.testing.assert_allclose(result['correlation'], np.corrcoef(pairs.T)[0, 1], atol=1e-14)
            np.testing.assert_allclose(result['variance_normalized_covariance'], np.cov(pairs.T, bias=True)[0, 1], atol=1e-14)


def test_empty_or_constant_pairs_have_no_fabricated_correlation():
    image = np.ones((8, 8))
    result = lag_statistic(image, np.zeros(image.shape, bool), 1, 0, 4)
    assert result == dict(n_pairs=0, variance_normalized_covariance=None, correlation=None)
    result = lag_statistic(image, np.ones(image.shape, bool), 1, 0, 4)
    assert result['n_pairs'] == 48 and result['correlation'] is None
    assert result['variance_normalized_covariance'] == 0


def test_duplicate_method_rows_cannot_silently_replace_each_other():
    with pytest.raises(ValueError, match='Duplicate row'):
        unique([dict(id='1', factor='2.5'), dict(id='1', factor='2.5')], ['id', 'factor'])
