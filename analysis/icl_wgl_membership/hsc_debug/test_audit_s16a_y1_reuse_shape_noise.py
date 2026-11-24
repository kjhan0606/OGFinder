import numpy as np
import pytest

from audit_s16a_y1_reuse_shape_noise import covariance


def test_shared_source_transfers_noise_between_tangential_and_cross():
    cov = covariance(np.array([0, 2]), ['object', 'object'], np.array([0., np.pi/4]), np.ones(2), 4)
    assert np.diag(cov) == pytest.approx(np.ones(4))
    assert cov[0, 3] == pytest.approx(-1.)
    assert cov[1, 2] == pytest.approx(1.)
    separate = covariance(np.array([0, 2]), ['pair1', 'pair2'], np.array([0., np.pi/4]), np.ones(2), 4)
    assert separate == pytest.approx(np.eye(4))


def test_reused_noise_cancels_when_same_angle_appears_with_opposite_weights():
    cov = covariance(np.array([0, 0]), ['object', 'object'], np.zeros(2), np.array([1., -1.]), 2)
    assert cov == pytest.approx(np.zeros((2, 2)))
