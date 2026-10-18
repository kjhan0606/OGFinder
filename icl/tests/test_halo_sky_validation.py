from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'analysis/icl_wgl_membership/hsc_debug'))
from validate_halo_sky import operator, solve, profile_result, make_columns


def fixture():
    data = dict(base=np.zeros((91, 91), bool), center=np.array([45., 45.]),
                members=np.array([[66., 49.]]))
    data['base'][5:9, :] = True
    cols = make_columns(data, 1., 10., True)
    variance = np.full(data['base'].shape, 2.)
    return data, cols, variance


def test_augmented_fit_matches_independent_scaled_svd_and_residual():
    data, cols, variance = fixture()
    truth = np.array([4000, 1500, 2000, 1000, .03, .02, -.01, .05, 0, .05])
    rng = np.random.default_rng(42)
    image = sum(c*t for c, t in zip(cols, truth)) + rng.normal(size=variance.shape)*np.sqrt(variance)
    anchor = truth[-6:] + rng.normal(size=6)*.01
    op, status = operator(cols, variance, data['base'], 6, .01)
    assert status == 'ok'
    result = solve(op, image, variance, anchor)
    valid = ~data['base']
    matrix = np.column_stack([c[valid]/np.sqrt(variance[valid]) for c in cols])
    constraint = np.zeros((6, len(cols))); constraint[:, -6:] = np.eye(6)/.01
    matrix = np.vstack([matrix, constraint])
    y = np.r_[image[valid]/np.sqrt(variance[valid]), anchor/.01]
    norms = np.linalg.norm(matrix, axis=0)
    x, _, rank, _ = np.linalg.lstsq(matrix/norms, y, rcond=None)
    assert rank == len(cols)
    np.testing.assert_allclose(result['coef'], x/norms, rtol=1e-9, atol=1e-7)
    np.testing.assert_allclose(result['chi2'], np.sum((y-matrix@(x/norms))**2), rtol=1e-11)
    _, s, vt = np.linalg.svd(matrix/norms, full_matrices=False)
    covariance = ((vt.T/s**2)@vt)/norms[:, None]/norms[None, :]
    np.testing.assert_allclose(result['covariance'], covariance, rtol=1e-8, atol=1e-8)


def test_noiseless_recovery_and_biased_sky_constraint_changes_diffuse():
    data, cols, variance = fixture()
    truth = np.array([4000, 1500, 2000, 1000, .03, .02, -.01, .05, 0, .05])
    image = sum(c*t for c, t in zip(cols, truth))
    op, _ = operator(cols, variance, data['base'], 6, .01)
    fitted = solve(op, image, variance, truth[-6:])
    np.testing.assert_allclose(fitted['coef'], truth, atol=1e-7)
    biased = truth[-6:].copy(); biased[[3, 5]] += .02
    wrong = solve(op, image, variance, biased)
    assert abs(wrong['coef'][1]-truth[1]) > 1
    assert wrong['chi2'] > fitted['chi2']


def test_profile_envelope_includes_ambiguity_but_not_unsupported_scales():
    def record(scale, dflux, chi2):
        return dict(scale=scale, coef=np.array([1., dflux, 1., 2.]),
                    covariance=np.eye(4)*4, chi2=chi2, condition=1., min_visible_fraction=.8)
    out = profile_result([record(8, 10, 20), record(10, 14, 20.75), record(12, 1000, 22)], 1)
    assert out['best_scale_kpc'] == 8
    assert out['supported_scale_count'] == 2
    assert out['profile_low'] == 8 and out['profile_high'] == 15
    assert out['valid_scale_count'] == 3
    assert profile_result([], 1)['status'] == 'no_valid_scale'


def test_invalid_geometry_and_anchor_are_not_silently_accepted():
    image = np.ones((20, 20)); mask = np.zeros_like(image, bool)
    assert operator([image, image], image, mask, 1)[1] == 'ill_conditioned'
    assert operator([image], image, ~mask, 1)[1] == 'too_few_pixels'
    with pytest.raises(ValueError):
        operator([image], image, mask, 1, 0.)
    op, _ = operator([image], image, mask, 1, .01)
    with pytest.raises(ValueError):
        solve(op, image, image)


def test_full_rank_hidden_component_still_gets_information_veto():
    fit = dict(scale=8., coef=np.ones(4), covariance=np.eye(4), chi2=0.,
               condition=2., min_visible_fraction=1e-10)
    result = profile_result([fit], 1)
    assert result['status'] == 'ok'
    assert result['information_status'] == 'blocked_hidden_component'


def test_unconstrained_fit_retains_signed_flux_and_variance_scaling():
    data, cols, variance = fixture()
    truth = np.array([4000, -500, 2000, 1000, .03, .02, -.01, .05, 0, .05])
    image = sum(c*t for c, t in zip(cols, truth))
    op, _ = operator(cols, variance, data['base'], 6)
    louder, _ = operator(cols, 4*variance, data['base'], 6)
    np.testing.assert_allclose(solve(op, image, variance)['coef'], truth, atol=1e-7)
    np.testing.assert_allclose(louder['covariance'], 4*op['covariance'], rtol=1e-9)
