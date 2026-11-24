from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import audit_s16a_cmodel_photoz_missingness as a


def test_weighted_sums_and_summary_match_direct_formulas():
    rng = np.random.default_rng(1)
    w, z, eff = rng.random(50), 1+rng.random(50), rng.random(50)
    keep = np.arange(50) % 3 != 0
    sums = {'all': a.weighted_sums(w, z, eff), 'available': a.weighted_sums(w*keep, z, eff)}
    out = a.summarise(sums, None)
    direct = lambda m: ((m*w*eff**2*z).sum()/(m*w*eff**2).sum(), (m*w*eff**3).sum()/(m*w*eff**2).sum())
    (za, ea), (zb, eb) = direct(np.ones(50)), direct(keep)
    assert np.isclose(out['mean_z_available_minus_all'], zb-za) and np.isclose(out['mean_eff_available_over_all_minus_1'], eb/ea-1)
    assert np.isclose(out['retained_eff2_weight_fraction'], (keep*w*eff**2).sum()/(w*eff**2).sum())


def test_bootstrap_sums_equal_count_weighted_sums():
    rng = np.random.default_rng(2)
    w, z, eff = rng.random(20), rng.random(20), rng.random(20)
    counts = [rng.poisson(1.0, size=(3, 20)).astype(float) for _ in range(2)]
    boot = a.weighted_sums(w, z, eff, counts).reshape(-1, 3)
    assert boot.shape == (6, 3)
    assert np.allclose(boot[4], a.weighted_sums(w*counts[1][1], z, eff))


def test_identical_thinning_gives_zero_difference_and_unit_retention():
    w, z, eff = np.ones(4), np.array([1., 2, 3, 4]), np.full(4, .5)
    sums = {'all': a.weighted_sums(w, z, eff), 'available': a.weighted_sums(w, z, eff)}
    out = a.summarise(sums, None)
    assert out['mean_z_available_minus_all'] == 0 and out['retained_eff2_weight_fraction'] == 1


def test_compare_reports_offset_scatter_and_outliers():
    rng = np.random.default_rng(3)
    ref = np.full(400, 1.0)
    code = ref+0.02
    code[:40] = 2.0
    out = a.compare(code, ref, rng)
    assert np.isclose(out['median_dz'], 0.02) and np.isclose(out['outlier_fraction'], 0.1)
    assert out['nmad_scaled'] == 0 and a.compare(code, None, rng).keys() >= {'n', 'median_z_code'}


def test_tract_jackknife_sees_field_level_scatter_that_source_bootstrap_misses():
    rng = np.random.default_rng(5)
    tract = np.repeat(np.arange(6), 400)
    fwhm = 0.8+0.05*tract+0.01*rng.standard_normal(len(tract))
    dz = 0.05*rng.standard_normal(6)[tract]+0.01*rng.standard_normal(len(tract))  # no true trend, field offsets
    out = a.slope(fwhm, dz, dz, tract, rng)
    assert out['n_tracts'] == 6 and out['slope_tract_jackknife_se'] > 3*out['slope_source_iid_bootstrap_se']
    assert out['slope_leave_one_tract_out_min'] <= out['slope'] <= out['slope_leave_one_tract_out_max']


def test_fixed_weight_stack_removes_the_lens_mix_term():
    # Lens A loses half its weight with no within-lens change; lens B is untouched. The pooled stack moves,
    # the fixed-weight stack does not.
    rng = np.random.default_rng(6)
    parts = {}
    for name, z, keep in (('A', 1.0, np.arange(10) < 5), ('B', 2.0, np.ones(10, bool))):
        w, zs, eff = np.ones(10), np.full(10, z), np.full(10, 0.5)
        counts = [rng.poisson(1.0, size=(4, 10)).astype(float)]
        sums = {'all': a.weighted_sums(w, zs, eff), 'available': a.weighted_sums(w*keep, zs, eff)}
        boot = {'all': a.weighted_sums(w, zs, eff, counts).reshape(-1, 3),
                'available': a.weighted_sums(w*keep, zs, eff, counts).reshape(-1, 3)}
        parts[name] = (sums, boot)
    pooled = a.summarise({s: parts['A'][0][s]+parts['B'][0][s] for s in ('all', 'available')}, None)
    fixed = a.fixed_weight_stack(parts)
    assert pooled['mean_z_available_minus_all'] > 0.1
    assert abs(fixed['mean_z_available_minus_all']) < 1e-12 and abs(fixed['mean_eff_available_over_all_minus_1']) < 1e-12


def test_untestable_when_a_group_has_no_finite_value_and_holm_is_monotone():
    rng = np.random.default_rng(4)
    values = np.array([np.nan, np.nan, 1., 2.])
    label = np.array([True, True, False, False])
    assert a.tested(values, np.zeros(4, int), label, rng)['testable'] is False
    rows = [dict(family='f', testable=True, permutation_p=p) for p in (0.01, 0.04, 0.03)] + [dict(family='f', testable=False)]
    a.holm(rows)
    assert [r['holm_p_within_family'] for r in rows[:3]] == [0.03, 0.06, 0.06] and 'holm_p_within_family' not in rows[3]
