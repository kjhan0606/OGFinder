import math

import numpy as np

import audit_s16a_afterburner_mask_effects as a


def test_tangent_plane_matches_small_angle_offsets():
    x, y = a.tangent_arcsec([150.01], [2.0], 150.0, 2.0)
    assert abs(x[0]-36*math.cos(math.radians(2.0))) < 1e-3 and abs(y[0]) < 1e-2
    x, y = a.tangent_arcsec([150.0], [2.01], 150.0, 2.0)
    assert abs(x[0]) < 1e-9 and abs(y[0]-36) < 1e-3


def test_grid_area_approximates_disc():
    x, y = a.grid(400.0)
    assert abs(len(x)*a.GRID_ARCSEC**2/(math.pi*400.0**2)-1) < 0.01


def test_ratio_and_jackknife():
    assert a.ratio(10, 2, 20, 8) == 2.0 and a.ratio(0, 1, 1, 1) is None and a.ratio(1, 0, 1, 1) is None
    assert a.jackknife([1.0, None]) is None
    values = [1.0, 2.0, 3.0, 4.0]
    assert abs(a.jackknife(values)-math.sqrt(3/4*sum((v-2.5)**2 for v in values))) < 1e-12


def test_stratified_test_detects_shift_and_ignores_stratum_offsets():
    rng = np.random.default_rng(0)
    group = np.repeat(np.arange(20), 100)
    label = np.tile(np.r_[np.ones(30, bool), np.zeros(70, bool)], 20)
    base = group*1.0+rng.normal(0, 0.1, len(group))
    null = a.stratified_test(base+0.0, group, label, rng)
    assert null['permutation_p'] > 0.01 and abs(null['stratified_mean_difference']) < 0.03
    shifted = a.stratified_test(base+0.1*label, group, label, rng)
    assert shifted['permutation_p'] < 0.01 and abs(shifted['stratified_mean_difference']-0.1) < 0.03


def test_stratified_test_drops_strata_without_both_groups():
    rng = np.random.default_rng(1)
    group = np.r_[np.zeros(50, int), np.ones(50, int)]
    label = np.r_[np.ones(50, bool), np.r_[np.ones(10, bool), np.zeros(40, bool)]]
    result = a.stratified_test(np.arange(100.0), group, label, rng)
    assert result['n_strata'] == 1 and result['n_true'] == 10 and result['n_false'] == 40
