from pathlib import Path
import sys

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_s16a_seeing_transfer_field_tests as a


def terms(n=400, nbins=4, seed=0, ratio=None):
    rng = np.random.default_rng(seed)
    fwhm = 0.7+0.2*rng.random(n)
    return dict(weight=np.ones(n), ratio=np.ones(n) if ratio is None else ratio(fwhm), fwhm=fwhm,
                index=np.arange(n) % nbins, nbins=nbins)


def lens(nbins=4):
    return dict(bin_eff=np.ones(nbins), fwhm_eff=1.05)


def test_weighted_quantile_matches_unweighted_median():
    x = np.arange(1, 102, dtype=float)
    assert np.isclose(a.weighted_quantile(x, np.ones_like(x), 0.5), 51, atol=0.6)


def test_bin_stats_returns_zero_for_empty_bins():
    t = terms(nbins=5)
    t['index'] = np.minimum(t['index'], 3)
    mu, f, mass = a.bin_stats(t, t['weight'])
    assert mass[4] == 0 and mu[4] == 0 and f[4] == 0 and np.all(np.isfinite(mu))


def test_split_response_is_zero_without_seeing_dependence():
    t = terms()
    high = t['fwhm'] > np.median(t['fwhm'])
    out = a.split_response(t, lens(), np.ones(len(high)), high)
    assert abs(out['change']) < 1e-12 and out['F_hi'] > out['F_lo'] and out['fallback_eff_fraction'] == 0


def test_split_response_recovers_a_linear_seeing_trend():
    # mu = 1 + 0.1 (FWHM - 0.8): the multiplier 1/mu changes by about -0.1 per arcsec.
    t = terms(n=4000, ratio=lambda f: 1+0.1*(f-0.8))
    high = t['fwhm'] > np.median(t['fwhm'])
    out = a.split_response(t, lens(), np.ones(len(high)), high)
    assert out['slope'] == pytest.approx(-0.1, rel=0.05)
    assert out['extrapolated'] == pytest.approx(out['slope']*(1.05-out['F_all']))


def test_split_response_refuses_unsupported_bin_unless_pooled_fallback_given():
    t = terms(nbins=4)
    high = t['fwhm'] > np.median(t['fwhm'])
    mult = (t['index'] != 2).astype(float)
    with pytest.raises(ValueError):
        a.split_response(t, lens(), mult, high)
    pooled = a.bin_stats(t, t['weight'])[:2]
    assert a.split_response(t, lens(), mult, high, pooled)['empty_pooled_bins']


def test_weighted_label_test_centres_the_composition_term():
    # Labels concentrate in a high-value stratum, with no within-stratum difference: the pooled difference is
    # large but the within-stratum effect is about zero and the test does not reject.
    rng = np.random.default_rng(1)
    strata = np.repeat([0, 1], 300)
    values = np.where(strata == 1, 2.0, 1.0)+0.01*rng.standard_normal(600)
    label = np.zeros(600, bool)
    label[300:450] = True
    label[:15] = True
    out = a.weighted_label_test(values, np.ones(600), strata, label, rng)
    assert out['pooled_difference'] > 0.5 and abs(out['within_strata_difference']) < 0.01 and out['permutation_p'] > 0.05


def test_weighted_label_test_recovers_a_uniform_within_stratum_shift():
    # Uneven label fractions (0.1 and 0.9) with a true within-stratum shift of 0.3: the direct estimator
    # returns 0.3, while the pooled-minus-null-mean column is attenuated.
    rng = np.random.default_rng(2)
    strata = np.repeat([0, 1], 500)
    label = np.zeros(1000, bool)
    label[:50] = True
    label[500:950] = True
    values = np.where(strata == 1, 1.0, 0.0)+0.3*label+0.01*rng.standard_normal(1000)
    out = a.weighted_label_test(values, np.ones(1000), strata, label, rng)
    assert out['within_strata_difference'] == pytest.approx(0.3, abs=0.005) and out['permutation_p'] < 0.01
    assert out['pooled_minus_null_mean'] < 0.2
