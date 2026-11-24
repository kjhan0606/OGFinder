import math

import numpy as np
from scipy.spatial import cKDTree

import audit_s16a_afterburner_mask_effects as old
import audit_s16a_boost_mask_field_tests as a


def part(n, area, post=None, unmasked=None):
    out = {key: np.zeros(a.NBINS) for key in a.KEYS}
    out['n_pre'][:] = out['w_pre'][:] = out['w2_pre'][:] = n
    out['n_post'][:] = out['w_post'][:] = out['w2_post'][:] = n if post is None else post
    out['a_fp'][:] = area
    out['a_unm'][:] = area if unmasked is None else unmasked
    return out


def test_grid_matches_recorded_runner_at_4_arcsec():
    x, y = old.grid(400.0)
    assert np.array_equal(a.grid(400.0, 4.0), np.c_[x, y])
    assert abs(len(a.grid(400.0, 1.0))/(math.pi*400.0**2)-1) < 1e-3


def test_rotation_identity_is_exact_and_quarter_turn_swaps_axes():
    xy = np.array([[1.5, -2.25], [0.0, 3.0]])
    assert np.array_equal(a.rotate(xy, 0.0), xy)
    assert np.allclose(a.rotate(xy, math.pi/2), np.array([[2.25, 1.5], [-3.0, 0.0]]))


def test_donor_mask_uses_nearest_donor_source_within_footprint():
    xy = np.array([[0.0, 0.0], [100.0, 0.0]])
    donor = dict(xy=xy, tree=cKDTree(xy), missing=np.array([False, True]), n=2)
    query = np.array([[95.0, 0.0], [5.0, 0.0], [300.0, 0.0], [0.0, 50.0]])
    mask, covered = a.donor_query(donor, query, 1.0, 0.0, 30.0)
    assert mask.tolist() == [True, False, False, False] and covered.tolist() == [True, True, False, False]
    assert a.donor_mask(donor, np.array([[50.0, 0.0]]), 2.0, 0.0, 30.0).tolist() == [True]
    assert a.donor_mask(donor, np.array([[0.0, -50.0]]), 2.0, math.pi/2, 30.0).tolist() == [True]


def test_parts_area_labels_weights_and_geometric_extra_mask():
    xy = np.array([[0.0, 0.0], [40.0, 0.0]])
    entry = dict(n=2, xy=xy, tree=cKDTree(xy), bins=np.zeros(2, int), weight=np.array([1.0, 2.0]), missing=np.array([False, True]))
    lab = dict(xy=np.array([[0.0, 0.0], [39.0, 0.0], [500.0, 0.0]]), bins=np.zeros(3, int), step=60.0)
    lab['dist'], lab['index'] = entry['tree'].query(lab['xy'])
    p = a.parts(entry, lab, 30.0, entry['missing'])
    assert p['a_fp'][0] == 2 and p['a_unm'][0] == 1 and p['n_post'][0] == 1 and p['w_post'][0] == 1.0 and p['w2_pre'][0] == 5.0
    g = a.parts(entry, lab, 30.0, np.array([True, True]), label_mask=entry['missing'], grid_extra=np.array([True, False, False]))
    assert g['a_unm'][0] == 0 and g['n_post'][0] == 0


def test_zero_count_is_boost_zero_with_one_inclusion_rule():
    lp, rs = part(0, 10), [part(30, 10)]
    assert a.density([lp], 'n', 'pre', (0,)) == 0.0 and a.single('pooled', lp, rs, 0, 'n', 'pre', ()) == 0.0
    value, used = a.stack('expected', [(lp, rs), (part(30, 10), [part(30, 10)])], 0, 'n', 'pre', ())
    assert used == [0, 1] and abs(value-0.5) < 1e-12
    assert a.single('pooled', part(5, 0), rs, 0, 'n', 'pre', ()) is None
    assert a.ratio_defined(0, 1, 2, 1) == 0.0 and a.ratio(0, 1, 2, 1) is None


def test_expected_count_stack_pairs_each_lens_with_its_own_field():
    items = [(part(10, 10), [part(40, 40)]), (part(90, 30), [part(60, 20)])]
    assert abs(a.stack('expected', items, 0, 'n', 'pre', ())[0]-1) < 1e-12
    assert abs(a.stack('pooled', items, 0, 'n', 'pre', ())[0]-1) > 0.1
    assert a.single('pooled', *items[0], 0, 'n', 'pre', ()) == a.single('expected', *items[0], 0, 'n', 'pre', ())


def test_normalised_definitions_and_correlated_poisson_se():
    lp, rs = part(20, 10), [part(30, 10)]
    lp['n_pre'][3] = 40
    assert a.single('outer_normalised', lp, rs, 7, 'n', 'pre', (7,)) == 1.0
    assert abs(a.single('outer_normalised', lp, rs, 3, 'n', 'pre', (7,))-2.0) < 1e-12
    assert abs(a.single('self_normalised', lp, [], 3, 'n', 'pre', (7,))-2.0) < 1e-12
    assert abs(a.stack('self_normalised', [(lp, [])], 3, 'n', 'pre', (7,))[0]-2.0) < 1e-12
    assert a.poisson_se('self_normalised', lp, [], 3, (7,), 'n', 'pre') == math.sqrt(1/40+1/20)
    assert a.poisson_se('self_normalised', lp, [], 7, (7,), 'n', 'pre') == 0.0
    assert abs(a.poisson_se('outer_normalised', lp, rs, 6, (6, 7), 'n', 'pre')-math.sqrt((1/20-1/40)+(1/30-1/60))) < 1e-12


def test_weighted_poisson_se_uses_effective_counts():
    lp, rs = part(4, 1), [part(4, 1)]
    lp['w_pre'][:], lp['w2_pre'][:] = 4.0, 8.0   # effective count 2
    assert abs(a.poisson_se('pooled', lp, rs, 0, (), 'w', 'pre')-math.sqrt(1/2+1/4)) < 1e-12


def test_decomposition_is_exact_and_separates_composition():
    # Two centres at densities 1 and 2; the dense one is half masked with no density change inside it.
    parts = [part(10, 10), part(20, 10, post=10, unmasked=5)]
    d = a.decompose(parts, 0)
    assert abs((1+d['composition'])*(1+d['within'])-(1+d['total'])) < 1e-12
    assert abs(d['within']) < 1e-12 and abs(d['composition']-((10+10)/15/1.5-1)) < 1e-12


def test_unmasked_density_change_handles_full_masking():
    assert a.unmasked_density_change(part(10, 2, post=0, unmasked=1), 0) == -1.0
    assert a.unmasked_density_change(part(10, 2, post=0, unmasked=0), 0) is None


def test_jackknife_strict_and_stack_jackknife_use_entering_items():
    assert a.jackknife_strict([1.0, None]) is None and a.jackknife_strict([1.0]) is None
    items = [(part(10, 10), [part(10, 10)]), (part(20, 10), [part(10, 10)]), (part(5, 0), [part(10, 10)])]
    value, used, se, replicas = a.stack_jackknife('expected', items, 0, 'n', 'pre', ())
    assert used == [0, 1] and len(replicas) == 2 and abs(value-1.5) < 1e-12 and se is not None


def test_null_stats_combined_test_and_holm():
    out = a.null_stats(3.0, [-1.0, 0.0, 1.0, 0.0])
    assert out['null_mean'] == 0.0 and out['p_two_sided'] == 1/5 and out['z'] > 3
    rng = np.random.default_rng(0)
    null = rng.normal(size=(2000, 3))
    assert a.combined_test([0.0, 0.0, 0.0], null)['sum_z2_p'] > 0.9
    assert a.combined_test([4.0, 4.0, 4.0], null)['mahalanobis_p'] < 0.01
    rows = [dict(p_two_sided=0.01), dict(p_two_sided=0.04), dict(p_two_sided=None)]
    a.holm(rows)
    assert rows[0]['holm_p'] == 0.02 and rows[1]['holm_p'] == 0.04 and rows[2]['holm_p'] is None and rows[0]['family_size'] == 2


def test_own_source_rule_drops_dilated_area_without_sources():
    lp, rs = part(10, 10), [part(30, 10), part(0, 10)]
    assert a.density(rs, 'n', 'pre', (0,)) == 1.5 and a.density(rs, 'n', 'pre', (0,), own=True) == 3.0
    empty = part(0, 4)
    assert a.single('pooled', empty, rs, 0, 'n', 'pre', ()) == 0.0
    assert a.single('pooled', empty, rs, 0, 'n', 'pre', (), own=True) is None
    masked = part(10, 10, post=0, unmasked=5)
    assert a.single('pooled', masked, [part(30, 10)], 0, 'n', 'area', (), own=True) == 0.0


def test_weighted_se_inside_normalisation_is_exact_delta_form():
    lp = part(10, 1)
    lp['w_pre'][:], lp['w2_pre'][:] = 20.0, 50.0
    w_k, w2_k, w_n, w2_n = 20.0, 50.0, 40.0, 100.0
    expected = math.sqrt(w2_k/w_k**2+w2_n/w_n**2-2*w2_k/(w_k*w_n))
    assert abs(a.poisson_se('self_normalised', lp, [], 6, (6, 7), 'w', 'pre')-expected) < 1e-12


def test_stack_only_restricts_items():
    items = [(part(10, 10), [part(10, 10)]), (part(30, 10), [part(10, 10)])]
    assert a.stack('expected', items, 0, 'n', 'pre', (), only={0}) == (1.0, [0])
