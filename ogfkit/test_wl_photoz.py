import numpy as np
import pytest

from ogfkit.wl_photoz import point_redshift_response, stack_photoz_response
from ogfkit.wlcluster import Cosmology


def test_weighted_foreground_dilution():
    # Identical estimated efficiencies. The reference foreground has 6/8
    # of the usable weight, and the estimated foreground never enters the fit.
    r = point_redshift_response(.3, [.8, .8, .2], [.8, .2, 1.2],
                                [2, 3, 100], [1, 2, 100], Cosmology())
    assert r['used_pairs'] == 2
    assert r['response'] == pytest.approx(.25)
    assert r['reference_foreground_weight_fraction'] == pytest.approx(.75)
    assert r['delta_sigma_multiplier'] == pytest.approx(4.)


def test_unequal_estimated_efficiencies_weight_the_estimator_denominator():
    cosmology = Cosmology()
    estimate = np.array([.5, 1.])
    truth = np.array([1., .6])
    inverse_est = 1/cosmology.sigma_crit(.3, estimate)
    inverse_true = 1/cosmology.sigma_crit(.3, truth)
    response = point_redshift_response(.3, estimate, truth, [3, 1], [1, 2], cosmology)
    numerator = np.sum(np.array([3, 2])*inverse_est*inverse_true)
    denominator = np.sum(np.array([3, 2])*inverse_est**2)
    assert response['response'] == pytest.approx(numerator/denominator)


def test_stack_weights_and_unsupported_cases():
    r = stack_photoz_response([.5, 1.2, np.nan], [2, 3, 0])
    assert r['response'] == pytest.approx(.92)
    assert r['delta_sigma_multiplier'] == pytest.approx(1/.92)
    with pytest.raises(ValueError, match='positive calibrated response'):
        stack_photoz_response([np.nan], [1])
    with pytest.raises(ValueError, match='no positive calibration weight'):
        point_redshift_response(.3, [.3], [1.], [1], [1], Cosmology())
    with pytest.raises(ValueError, match='reference response is zero'):
        point_redshift_response(.3, [.8], [.2], [1], [1], Cosmology())
