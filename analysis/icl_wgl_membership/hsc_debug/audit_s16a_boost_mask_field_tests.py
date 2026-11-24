#!/usr/bin/env python3
"""Offline follow-up of the S16A afterburner-mask boost diagnostics (2026-11-24 candidates b and c).

b1. Full-field transplant: the S23B_3620 lens and random masks are placed on the lens and on the
    randoms of every other lens group with randoms. The unmasked area comes either from the donor
    mask evaluated at each grid point (geometric, on top of the recipient's own labels) or from the
    recipient's nearest-source labels, as in the real analysis.
b2. Rotation null: each centre's own mask is rotated about the centre. This keeps the radial mask
    profile (where the rotated mask stays inside the footprint) and removes any link to the local
    source density. The area-normalised residual is split exactly into a composition term (pooled
    centres reweighted by unmasked area) and a within-centre term; real and rotated values are compared.
b3. Grid resolution: per-annulus boosts on 1, 2 and 4 arcsec grids. The 4-arcsec values must
    reproduce the recorded mask-effects table, and the unrotated lens-side transplant its table.
c.  Boost normalisation: pooled, per-field expected-count, outer-annulus normalised and lens-only
    self-normalised definitions, on the lenses and on random centres used as pseudo-lenses.

Reads verified local inputs only. Descriptive sensitivities: no shear, mass, correction or selection
change. Constants, realisation counts and seeds are fixed before the run. Revision 2 follows the
independent review of the first run (see summary.json 'revision').
"""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_s16a_afterburner_mask_effects import DONOR_LENS, FOOTPRINTS_ARCSEC, PRIMARY_FOOTPRINT, annulus, ratio, tangent_arcsec
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs, sha, write_table
from audit_s16a_y1_reuse_shape_noise import GEOMETRY, GEOMETRY_SHA
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES

PRIOR = Path('/home/kjhan/BACKUP/Astrafek/s16a_afterburner_mask_effects_20261124_02')
PRIOR_SUMMARY_SHA = 'f658b60508794d72f6bc421e1bb8f80da6eda284b053cbb57096f38bdf4c9676'
PRIOR_RUNNER_SHA = 'daf79c9c508c00c7092b1d518d89be7e590b721f8929ca01be3c7f6a50b334e8'
GRIDS_ARCSEC = (1.0, 2.0, 4.0)
NULL_GRID_ARCSEC = 2.0
ANCHOR_GRID_ARCSEC = 4.0
TRANSPLANT_REALISATIONS = 256
ROTATION_REALISATIONS = 2000
ROTATION_FOOTPRINTS = (20.0, 30.0, 45.0)
COVERAGE_FLAG = 0.8
OWN = True  # normalisation section: a centre enters annulus k only with an own pre-mask source there
NORMALISATION = {'2-3': (7,), '1.5-3': (6, 7)}
SEED = 2026112503
WORKERS = 16
NBINS = len(RADIAL_EDGES)-1
KEYS = ('a_fp', 'a_unm', 'n_pre', 'n_post', 'w_pre', 'w_post', 'w2_pre', 'w2_post')
DEFINITIONS = ('pooled', 'expected', 'outer_normalised', 'self_normalised')
REVISION = ['Zero counts over a positive area are boost 0 (well defined), with one inclusion rule for per-centre and stacked '
            'values; the first run dropped them from per-centre summaries only.',
            'Poisson SE inside the normalisation annuli uses the correlated form 1/n_k - 1/N; weighted rows use effective counts.',
            'Per-centre rotation deltas are -1 when every source is masked over a positive unmasked area (first run: dropped).',
            'Pseudo-lens stacks are reported only for the self-normalised definition; for pooled, expected and outer-normalised '
            'definitions leave-one-out closure within a group pins them near 1 by construction.',
            'Jackknives run over the items that enter the annulus and are undefined if a replicate is undefined.',
            'Rotation null: 2000 realisations, footprints 20/30/45 arcsec, exact composition/within-centre split, combined '
            'tests across annuli, Holm only over rows with masked sources, per-centre rotated coverage.',
            'Transplant: 256 realisations with Monte Carlo SE of the null means; grid table adds expected-count stacks.',
            'Revision 3 (second review): normalisation uses the own-source rule (a centre enters annulus k only with an own '
            'pre-mask source there; zero-count cells in revision 2 were footprint-dilation artefacts); exact weighted SE inside '
            'the normalisation annuli; SD/Poisson over one set; STACK_expected pre and area over the same groups; '
            'jackknife guards; group-level rotated coverage with edge-affected groups (a masked centre below the coverage '
            'flag) excluded from Holm, combined tests and STACK_core; composition combined tests reported as informational; '
            'nulls on a 2-arcsec grid (anchors stay at 4 arcsec).']


def grid(radius_arcsec, step):
    n = int(math.ceil(radius_arcsec/step))
    axis = np.arange(-n, n+1)*step
    x, y = np.meshgrid(axis, axis)
    keep = x**2+y**2 <= radius_arcsec**2
    return np.c_[x[keep], y[keep]]


def rotate(xy, theta):
    c, s = math.cos(theta), math.sin(theta)
    return np.c_[c*xy[:, 0]-s*xy[:, 1], s*xy[:, 0]+c*xy[:, 1]]


def add(parts):
    total = {k: np.zeros(NBINS) for k in KEYS}
    for part in parts:
        for k in KEYS:
            total[k] = total[k]+part[k]
    return total


def ratio_defined(nl, al, nr, ar):
    """As the recorded ratio, but a zero lens count over positive areas is 0 rather than undefined."""
    return (nl/al)/(nr/ar) if al > 0 and ar > 0 and nr > 0 else None


def boost(lp, rp, k, kind='n', state='area', strict=False):
    f = ratio if strict else ratio_defined
    if state == 'pre':
        return f(lp[kind+'_pre'][k], lp['a_fp'][k], rp[kind+'_pre'][k], rp['a_fp'][k])
    area = 'a_fp' if state == 'naive' else 'a_unm'
    return f(lp[kind+'_post'][k], lp[area][k], rp[kind+'_post'][k], rp[area][k])


def over(a, b):
    return a/b-1 if a is not None and b is not None and b != 0 else None


def jackknife_strict(replicas):
    if len(replicas) < 2 or any(r is None for r in replicas):
        return None
    v = np.asarray(replicas, float)
    return float(math.sqrt((len(v)-1)/len(v)*((v-v.mean())**2).sum()))


def counts(p, kind, state):
    return p[kind+('_pre' if state == 'pre' else '_post')], p['a_fp' if state == 'pre' else 'a_unm']


def density(parts, kind, state, ks, own=False):
    """Pooled density over annuli ks. With own=True a centre contributes annulus k only if it has at least one
    of its own pre-mask sources there: label area dilated in from neighbouring annuli alone is not counted."""
    num = area = 0.0
    for p in parts:
        n, a = counts(p, kind, state)
        for k in ks:
            if own and p['n_pre'][k] <= 0:
                continue
            num, area = num+n[k], area+a[k]
    return num/area if area > 0 else None


def positive(x):
    return x is not None and x > 0


def single(definition, lp, randoms, k, kind, state, norm, own=False):
    """Boost of one lens-like centre in annulus k; randoms is a list of parts of the same field. A zero lens
    count over a positive area gives 0; undefined (None) when an area or a reference density vanishes, or
    (own=True) when the lens-like centre has no own pre-mask source in annulus k."""
    if own and lp['n_pre'][k] <= 0:
        return None
    rho_l = density([lp], kind, state, (k,), own)
    if rho_l is None:
        return None
    if definition == 'self_normalised':
        ref = density([lp], kind, state, norm, own)
        return rho_l/ref if positive(ref) else None
    rho_r = density(randoms, kind, state, (k,), own)
    if not positive(rho_r):
        return None
    if definition in ('pooled', 'expected'):
        return rho_l/rho_r
    sl, sr = density([lp], kind, state, norm, own), density(randoms, kind, state, norm, own)
    return (rho_l/rho_r)/(sl/sr) if positive(sl) and positive(sr) else None


def stack(definition, items, k, kind, state, norm, own=False, only=None):
    """items: (lens parts, list of random parts). An item enters annulus k exactly when single() is defined
    (and, if given, its index is in only). Returns (value, indices used)."""
    used, num, expected, pooled_l, pooled_r = [], 0.0, 0.0, [], []
    for i, (lp, randoms) in enumerate(items):
        if (only is not None and i not in only) or single(definition, lp, randoms, k, kind, state, norm, own) is None:
            continue
        used.append(i)
        n_l, a_l = (v[k] for v in counts(lp, kind, state))
        if definition == 'pooled':
            pooled_l.append(lp)
            pooled_r.extend(randoms)
        elif definition == 'self_normalised':
            num, expected = num+n_l, expected+density([lp], kind, state, norm, own)*a_l
        else:
            scale = 1.0
            if definition == 'outer_normalised':
                scale = density([lp], kind, state, norm, own)/density(randoms, kind, state, norm, own)
            num, expected = num+n_l, expected+density(randoms, kind, state, (k,), own)*scale*a_l
    if not used:
        return None, used
    if definition == 'pooled':
        return density(pooled_l, kind, state, (k,), own)/density(pooled_r, kind, state, (k,), own), used
    return (num/expected if expected > 0 else None), used


def stack_jackknife(definition, items, k, kind, state, norm, labels=None, own=False):
    """Delete-one jackknife SE of (stack - 1) over the items (or label groups) that enter annulus k."""
    value, used = stack(definition, items, k, kind, state, norm, own)
    labels = list(range(len(items))) if labels is None else labels
    drops = sorted({labels[i] for i in used})
    replicas = [over(stack(definition, [it for j, it in enumerate(items) if labels[j] != d], k, kind, state, norm, own)[0], 1.0) for d in drops]
    return value, used, jackknife_strict(replicas), replicas


def sums(parts, kind, state, ks, own=False):
    """Sum of weights W and of squared weights W2 (counts for kind 'n') over annuli ks, with the own rule."""
    s, w, w2 = ('_pre' if state == 'pre' else '_post'), 0.0, 0.0
    for p in parts:
        for k in ks:
            if own and p['n_pre'][k] <= 0:
                continue
            w += p[('n' if kind == 'n' else 'w')+s][k]
            w2 += p[('n' if kind == 'n' else 'w2')+s][k]
    return w, w2


def effective(parts, kind, state, ks, own=False):
    w, w2 = sums(parts, kind, state, ks, own)
    return w*w/w2 if w2 > 0 else 0.0


def poisson_se(definition, lp, randoms, k, norm, kind, state, own=False):
    """Delta-method relative SE of the boost. A side contributes W2_k/W_k^2; a normalised side adds
    W2_N/W_N^2, minus 2 W2_k/(W_k W_N) when annulus k lies inside the normalisation annuli (shared sources).
    For counts this is 1/n_k + 1/N, or 1/n_k - 1/N inside."""
    sides = [[lp]] if definition == 'self_normalised' else [[lp], randoms]
    var = 0.0
    for parts in sides:
        w_k, w2_k = sums(parts, kind, state, (k,), own)
        if w_k <= 0:
            return None
        var += w2_k/w_k**2
        if definition in ('outer_normalised', 'self_normalised'):
            w_n, w2_n = sums(parts, kind, state, norm, own)
            if w_n <= 0:
                return None
            var += w2_n/w_n**2-(2*w2_k/(w_k*w_n) if k in norm else 0.0)
    return math.sqrt(max(var, 0.0))


def decompose(parts, k, kind='n'):
    """Exact split of a pooled set's unmasked/pre-mask density ratio: (1+total) = (1+composition)(1+within).
    composition: pre-mask densities reweighted by unmasked area; within: unmasked density over that baseline."""
    a = np.array([p['a_fp'][k] for p in parts])
    u = np.array([p['a_unm'][k] for p in parts])
    pre = np.array([p[kind+'_pre'][k] for p in parts])
    post = np.array([p[kind+'_post'][k] for p in parts])
    if a.sum() <= 0 or u.sum() <= 0 or pre.sum() <= 0:
        return None
    rho = np.divide(pre, a, out=np.zeros_like(pre), where=a > 0)
    baseline = (rho*u).sum()/u.sum()
    if baseline <= 0:
        return None
    full = pre.sum()/a.sum()
    return dict(total=(post.sum()/u.sum())/full-1, composition=baseline/full-1, within=(post.sum()/u.sum())/baseline-1)


def unmasked_density_change(p, k):
    """Per-centre (unmasked density)/(pre-mask density) - 1; -1 when every source is masked over positive area."""
    if p['n_pre'][k] <= 0 or p['a_fp'][k] <= 0 or p['a_unm'][k] <= 0:
        return None
    return (p['n_post'][k]/p['a_unm'][k])/(p['n_pre'][k]/p['a_fp'][k])-1


def load_centres():
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    if sha(GEOMETRY/'summary.json') != GEOMETRY_SHA:
        raise ValueError('geometry summary changed')
    geometry = json.loads((GEOMETRY/'summary.json').read_text())
    if sha(GEOMETRY/'pair_geometry.tsv') != geometry['output_sha256']['pair_geometry.tsv']:
        raise ValueError('geometry table changed')
    by_center = defaultdict(list)
    with (GEOMETRY/'pair_geometry.tsv').open() as stream:
        for pair, geo in zip(pairs, csv.DictReader(stream, delimiter='\t'), strict=True):
            if (pair['center_id'], pair['object_id']) != (geo['center_id'], geo['object_id']):
                raise ValueError('pair/geometry mismatch')
            if pair['object_id'] in selected:
                by_center[pair['center_id']].append((pair, float(geo['radius_proper_mpc']), int(geo['radial_bin'])))
    data, radius_error = {}, 0.0
    for c in centers:
        rows = by_center[c['center_id']]
        entry = dict(center=c, id=c['center_id'], role=c['role'], group=c['matched_lens_id'], radius=float(c['radius_arcsec']),
                     ids=[p['object_id'] for p, _, _ in rows], n=len(rows))
        if rows:
            x, y = tangent_arcsec([float(p['source_ra']) for p, _, _ in rows], [float(p['source_dec']) for p, _, _ in rows],
                                  c['ra_deg'], c['dec_deg'])
            radius_error = max(radius_error, float(np.max(np.abs(np.hypot(x, y)*3.0/entry['radius']-np.array([r for _, r, _ in rows])))))
            entry.update(xy=np.c_[x, y], bins=np.array([b for _, _, b in rows]), weight=np.array([float(p['weight']) for p, _, _ in rows]),
                         missing=np.array([o not in available for o in entry['ids']]))
            entry['tree'] = cKDTree(entry['xy'])
        data[c['center_id']] = entry
    if radius_error > 5e-3:
        raise ValueError(f'tangent-plane radius differs from pair geometry by {radius_error} Mpc')
    return centers, data, provenance, radius_error


def labels(entry, step):
    xy = grid(entry['radius'], step)
    bins = annulus(np.hypot(xy[:, 0], xy[:, 1])*3.0/entry['radius'])
    if not entry['n']:
        return dict(xy=xy, bins=bins, dist=np.full(len(xy), np.inf), index=np.zeros(len(xy), int), step=step)
    dist, index = entry['tree'].query(xy, workers=WORKERS)
    return dict(xy=xy, bins=bins, dist=dist, index=index, step=step)


def parts(entry, lab, footprint, source_mask, label_mask=None, grid_extra=None):
    """Per-annulus footprint and unmasked areas (arcmin^2), counts, weights and squared weights before/after
    source_mask. Unmasked grid points are labelled by the nearest source's label_mask (default source_mask);
    a geometric grid_extra mask is ORed on top when given."""
    out = {k: np.zeros(NBINS) for k in KEYS}
    if not entry['n']:
        return out
    inside = lab['dist'] <= footprint
    masked = (source_mask if label_mask is None else label_mask)[lab['index']]
    if grid_extra is not None:
        masked = masked | grid_extra
    cell = lab['step']**2/3600
    w, keep = entry['weight'], ~source_mask
    out['a_fp'] = np.bincount(lab['bins'][inside], minlength=NBINS)*cell
    out['a_unm'] = np.bincount(lab['bins'][inside & ~masked], minlength=NBINS)*cell
    out['n_pre'] = np.bincount(entry['bins'], minlength=NBINS).astype(float)
    out['n_post'] = np.bincount(entry['bins'][keep], minlength=NBINS).astype(float)
    out['w_pre'] = np.bincount(entry['bins'], weights=w, minlength=NBINS)
    out['w_post'] = np.bincount(entry['bins'][keep], weights=w[keep], minlength=NBINS)
    out['w2_pre'] = np.bincount(entry['bins'], weights=w*w, minlength=NBINS)
    out['w2_post'] = np.bincount(entry['bins'][keep], weights=w[keep]**2, minlength=NBINS)
    return out


def donor_query(donor, xy, scale, theta, footprint):
    """Donor's nearest-source mask at recipient positions rotated by theta and scaled to the donor frame,
    and whether each position lies within the donor footprint at all."""
    dist, index = donor['tree'].query(rotate(xy, theta)*scale, workers=WORKERS)
    covered = dist <= footprint
    return covered & donor['missing'][np.minimum(index, donor['n']-1)], covered


def donor_mask(donor, xy, scale, theta, footprint):
    return donor_query(donor, xy, scale, theta, footprint)[0]


def null_stats(real, null):
    null = np.asarray([v for v in null if v is not None], float)
    if real is None or len(null) < 2:
        return dict(null_n=len(null), null_mean=None, null_sd=None, z=None, p_two_sided=None)
    mean, sd = float(null.mean()), float(null.std(ddof=1))
    return dict(null_n=len(null), null_mean=mean, null_sd=sd, z=(real-mean)/sd if sd > 0 else None,
                p_two_sided=float((1+(np.abs(null-mean) >= abs(real-mean)).sum())/(1+len(null))))


def combined_test(real, null):
    """real: vector over annuli; null: realisations x annuli (complete rows only). Sum of z^2 and Mahalanobis
    distance, each with an empirical p from the same null."""
    real, null = np.asarray(real, float), np.asarray(null, float)
    sd = null.std(axis=0, ddof=1)
    use = sd > 0
    if use.sum() == 0 or len(null) < 10:
        return dict(annuli_used=int(use.sum()), sum_z2=None, sum_z2_p=None, mahalanobis=None, mahalanobis_p=None)
    mean = null[:, use].mean(axis=0)
    z_null, z_real = (null[:, use]-mean)/sd[use], (real[use]-mean)/sd[use]
    s_null, s_real = (z_null**2).sum(1), float((z_real**2).sum())
    inv = np.linalg.pinv(np.cov(z_null, rowvar=False).reshape(use.sum(), use.sum()))
    m_null = np.einsum('ij,jk,ik->i', z_null, inv, z_null)
    m_real = float(z_real@inv@z_real)
    return dict(annuli_used=int(use.sum()), sum_z2=s_real, sum_z2_p=float((1+(s_null >= s_real).sum())/(1+len(null))),
                mahalanobis=m_real, mahalanobis_p=float((1+(m_null >= m_real).sum())/(1+len(null))))


def holm(rows, key='p_two_sided'):
    usable = sorted((r for r in rows if r[key] is not None), key=lambda r: r[key])
    running = 0.0
    for rank, r in enumerate(usable):
        running = max(running, min(1.0, (len(usable)-rank)*r[key]))
        r['holm_p'] = running
    for r in rows:
        r.setdefault('holm_p', None)
        r['family_size'] = len(usable)


def num(text):
    return None if text == '' else float(text)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    if sha(PRIOR/'summary.json') != PRIOR_SUMMARY_SHA:
        raise ValueError('prior mask-effects summary changed')
    prior = json.loads((PRIOR/'summary.json').read_text())
    if prior['runner_sha256'] != PRIOR_RUNNER_SHA:
        raise ValueError('prior runner differs')
    for name, digest in prior['output_sha256'].items():
        if sha(PRIOR/name) != digest:
            raise ValueError('prior table changed '+name)
    rng = np.random.default_rng(SEED)
    centers, data, provenance, radius_error = load_centres()
    lenses = {c['matched_lens_id']: data[c['center_id']] for c in centers if c['role'] == 'lens'}
    randoms = defaultdict(list)
    for c in centers:
        if c['role'] == 'random':
            randoms[c['matched_lens_id']].append(data[c['center_id']])
    nonempty = {g: sorted((r for r in rs if r['n']), key=lambda r: r['id']) for g, rs in randoms.items()}
    fp = PRIMARY_FOOTPRINT

    # b3. Grid resolution, with exact reproduction of the recorded 4-arcsec table.
    real, grid_area = {}, {}
    for step in GRIDS_ARCSEC:
        footprints = FOOTPRINTS_ARCSEC if step == 4.0 else (fp,)
        cache = {cid: labels(e, step) for cid, e in data.items()}
        for cid, e in data.items():
            disc = np.bincount(cache[cid]['bins'], minlength=NBINS)*step**2/3600
            exact = np.pi*(np.square(RADIAL_EDGES[1:])-np.square(RADIAL_EDGES[:-1]))*(e['radius']/3.0/60)**2
            grid_area[(step, cid)] = disc/exact-1
        for footprint in footprints:
            real[(step, footprint)] = {cid: parts(e, cache[cid], footprint, e.get('missing')) for cid, e in data.items()}
        if step == NULL_GRID_ARCSEC:
            lab_null = cache
        if step == ANCHOR_GRID_ARCSEC:
            lab_anchor = cache
        del cache
    supported = [g for g in lenses if add([real[(4.0, fp)][r['id']] for r in randoms[g]])['a_fp'].sum() > 0]

    def group_parts(step, footprint, g):
        table = real[(step, footprint)]
        return table[lenses[g]['id']], add([table[r['id']] for r in randoms[g]])

    def items_for(step, footprint, subset):
        table = real[(step, footprint)]
        return [(table[lenses[g]['id']], [table[r['id']] for r in nonempty[g]]) for g in subset]

    def matched(step, footprint, subset, k):
        use = [g for g in subset if all(p['a_fp'][k] > 0 for p in group_parts(step, footprint, g))]
        return add([group_parts(step, footprint, g)[0] for g in use]), add([group_parts(step, footprint, g)[1] for g in use]), len(use)

    reproduction = dict(rows=0, rows_with_defined_boost=0, max_relative_difference=0.0)
    with (PRIOR/'annulus_boost.tsv').open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            footprint, g, k = float(row['footprint_arcsec']), row['lens_id'], int(row['annulus'])
            if g.startswith('STACK:'):
                subset = supported if g == 'STACK:all_supported' else [x for x in supported if x != DONOR_LENS]
                lp, rp, _ = matched(4.0, footprint, subset, k)
            else:
                lp, rp = group_parts(4.0, footprint, g)
            mine = dict(lens_n_pre=lp['n_pre'][k], lens_n_post=lp['n_post'][k], lens_area_fp=lp['a_fp'][k], lens_area_unmasked=lp['a_unm'][k],
                        random_n_pre=rp['n_pre'][k], random_n_post=rp['n_post'][k], random_area_fp=rp['a_fp'][k],
                        random_area_unmasked=rp['a_unm'][k])
            for kind, suffix in (('n', ''), ('w', '_weighted')):
                for state in ('pre', 'naive', 'area'):
                    mine[f'boost_{state}{suffix}'] = boost(lp, rp, k, kind, state, strict=True)
            for key, value in mine.items():
                recorded = num(row[key])
                if (recorded is None) != (value is None):
                    raise ValueError(f'4-arcsec reproduction differs: {g} {k} {key}')
                if recorded is not None:
                    reproduction['max_relative_difference'] = max(reproduction['max_relative_difference'],
                                                                  abs(value-recorded)/max(abs(recorded), 1e-300))
            reproduction['rows'] += 1
            reproduction['rows_with_defined_boost'] += mine['boost_area'] is not None
    if reproduction['max_relative_difference'] > 1e-12:
        raise ValueError('4-arcsec boost table not reproduced')

    grid_rows = []
    for step in GRIDS_ARCSEC:
        for g in [*supported, 'STACK:all_supported', 'STACK:without_'+DONOR_LENS]:
            for k in range(NBINS):
                expected = {}
                if g.startswith('STACK:'):
                    subset = supported if g == 'STACK:all_supported' else [x for x in supported if x != DONOR_LENS]
                    lp, rp, n_groups = matched(step, fp, subset, k)
                    disc_error = None
                    for kind, suffix in (('n', ''), ('w', '_weighted')):
                        pre, area = (stack('expected', items_for(step, fp, subset), k, kind, s, ())[0] for s in ('pre', 'area'))
                        expected.update({'boost_pre_expected'+suffix: pre, 'boost_area_expected'+suffix: area,
                                         'area_minus_pre_expected'+suffix: over(area, pre)})
                else:
                    lp, rp = group_parts(step, fp, g)
                    n_groups, disc_error = None, float(grid_area[(step, lenses[g]['id'])][k])
                row = dict(grid_arcsec=step, footprint_arcsec=fp, lens_id=g, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                           stacked_groups=n_groups, lens_grid_disc_area_relative_error=disc_error,
                           lens_area_fp=lp['a_fp'][k], lens_area_unmasked=lp['a_unm'][k], random_area_fp=rp['a_fp'][k],
                           random_area_unmasked=rp['a_unm'][k])
                for kind, suffix in (('n', ''), ('w', '_weighted')):
                    values = {s: boost(lp, rp, k, kind, s) for s in ('pre', 'naive', 'area')}
                    row.update({f'boost_{s}{suffix}': v for s, v in values.items()})
                    row['area_minus_pre'+suffix] = over(values['area'], values['pre'])
                row.update(expected)
                grid_rows.append(row)
    reference = {(r['lens_id'], r['annulus']): r for r in grid_rows if r['grid_arcsec'] == 1.0}
    for r in grid_rows:
        base = reference[(r['lens_id'], r['annulus'])]
        for key in ('area_minus_pre', 'boost_pre', 'area_minus_pre_expected', 'boost_pre_expected'):
            if key in r:
                r[key+'_minus_1arcsec'] = r[key]-base[key] if r[key] is not None and base[key] is not None else None

    # b1. Full-field transplant of the donor lens and random masks.
    real_null = real[(NULL_GRID_ARCSEC, fp)]
    donor_lens, donor_randoms = lenses[DONOR_LENS], nonempty[DONOR_LENS]
    recipients = [g for g in supported if g != DONOR_LENS]

    def transplanted(entry, donor, theta, lab_table=None):
        scale = donor['radius']/entry['radius']
        lab = (lab_null if lab_table is None else lab_table)[entry['id']]
        source_t = donor_mask(donor, entry['xy'], scale, theta, fp)
        grid_t = donor_mask(donor, lab['xy'], scale, theta, fp)
        combined = entry['missing'] | source_t
        return (parts(entry, lab, fp, combined, label_mask=entry['missing'], grid_extra=grid_t),
                parts(entry, lab, fp, combined))

    # Anchor: unrotated lens-side-only geometric transplant reproduces the recorded transplant table.
    anchor = dict(rows=0, rows_with_defined_value=0, max_relative_difference=0.0,
                  scope='geometric, lens side only, unweighted area_after_over_area_before at theta=0, donor lens:'+DONOR_LENS)
    with (PRIOR/'transplant_boost.tsv').open() as stream:
        recorded = [r for r in csv.DictReader(stream, delimiter='\t') if r['donor'] == 'lens:'+DONOR_LENS and not r['recipient'].startswith('STACK:')]
    lens_only = {}
    for r in recorded:
        g, k = r['recipient'], int(r['annulus'])
        if g not in lens_only:
            lens_only[g] = transplanted(lenses[g], donor_lens, 0.0, lab_anchor)[0]
        lp0, rp = group_parts(ANCHOR_GRID_ARCSEC, fp, g)
        after, before = boost(lens_only[g], rp, k, strict=True), boost(lp0, rp, k, strict=True)
        value = after/before if after and before else None
        if (num(r['area_after_over_area_before']) is None) != (value is None):
            raise ValueError('lens-side transplant anchor differs '+g)
        if value is not None:
            anchor['max_relative_difference'] = max(anchor['max_relative_difference'], abs(value/num(r['area_after_over_area_before'])-1))
            anchor['rows_with_defined_value'] += 1
        anchor['rows'] += 1
    if anchor['max_relative_difference'] > 1e-12:
        raise ValueError('lens-side transplant table not reproduced')

    transplant_rows = []
    for realisation in range(TRANSPLANT_REALISATIONS+1):
        after = {}
        for g in recipients:
            rs = nonempty[g]
            if realisation == 0:
                thetas, assign = np.zeros(1+len(rs)), [i % len(donor_randoms) for i in range(len(rs))]
            else:
                thetas, assign = rng.uniform(0, 2*np.pi, 1+len(rs)), rng.permutation(len(donor_randoms))[:len(rs)]
            lens_g, lens_n = transplanted(lenses[g], donor_lens, thetas[0])
            rand = [transplanted(r, donor_randoms[a], t) for r, a, t in zip(rs, assign, thetas[1:])]
            after[g] = dict(geometric=(lens_g, add([x[0] for x in rand])), nearest=(lens_n, add([x[1] for x in rand])))
        for g in [*recipients, 'STACK']:
            for k in range(NBINS):
                if g == 'STACK':
                    use = [x for x in recipients if all(p['a_fp'][k] > 0 for p in group_parts(NULL_GRID_ARCSEC, fp, x))]
                    before = (add([group_parts(NULL_GRID_ARCSEC, fp, x)[0] for x in use]), add([group_parts(NULL_GRID_ARCSEC, fp, x)[1] for x in use]))
                    modes = {m: (add([after[x][m][0] for x in use]), add([after[x][m][1] for x in use])) for m in ('geometric', 'nearest')}
                else:
                    use, before, modes = [g], group_parts(NULL_GRID_ARCSEC, fp, g), after[g]
                nl, nr = modes['nearest'][0]['n_post'][k], modes['nearest'][1]['n_post'][k]
                row = dict(realisation=realisation, recipient=g, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], stacked_groups=len(use),
                           lens_transplanted_count_fraction=1-nl/before[0]['n_post'][k] if before[0]['n_post'][k] else None,
                           random_transplanted_count_fraction=1-nr/before[1]['n_post'][k] if before[1]['n_post'][k] else None,
                           poisson_relative_se_after=math.sqrt(1/nl+1/nr) if nl > 0 and nr > 0 else None)
                for kind, suffix in (('n', ''), ('w', '_weighted')):
                    b0 = boost(*before, k, kind)
                    for m in ('geometric', 'nearest'):
                        row[f'{m}_after_over_before_minus_1{suffix}'] = over(boost(*modes[m], k, kind), b0)
                    gv, nv = row[f'geometric_after_over_before_minus_1{suffix}'], row[f'nearest_after_over_before_minus_1{suffix}']
                    row['nearest_over_geometric_minus_1'+suffix] = (1+nv)/(1+gv)-1 if gv is not None and nv is not None and gv != -1 else None
                transplant_rows.append(row)
    transplant_summary = []
    for g in [*recipients, 'STACK']:
        for k in range(NBINS):
            rows = [r for r in transplant_rows if r['recipient'] == g and r['annulus'] == k]
            primary, null = rows[0], rows[1:]
            out = dict(recipient=g, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], stacked_groups=primary['stacked_groups'],
                       lens_transplanted_count_fraction=primary['lens_transplanted_count_fraction'],
                       random_transplanted_count_fraction=primary['random_transplanted_count_fraction'],
                       poisson_relative_se_after=primary['poisson_relative_se_after'])
            for suffix in ('', '_weighted'):
                for base_key in ('geometric_after_over_before_minus_1', 'nearest_after_over_before_minus_1', 'nearest_over_geometric_minus_1'):
                    key = base_key+suffix
                    values = np.asarray([r[key] for r in null if r[key] is not None], float)
                    out[key+'_primary'] = primary[key]
                    out[key+'_null_n'] = len(values)
                    out[key+'_null_mean'] = float(values.mean()) if len(values) else None
                    out[key+'_null_sd'] = float(values.std(ddof=1)) if len(values) > 1 else None
                    out[key+'_null_mean_mc_se'] = float(values.std(ddof=1)/math.sqrt(len(values))) if len(values) > 1 else None
            transplant_summary.append(out)

    # b2. Rotation null of each centre's own mask (nearest-source areas, as in the real analysis).
    rotation_centres = {g: [lenses[g], *nonempty[g]] for g in supported}
    flat = [e for es in rotation_centres.values() for e in es]
    thetas = rng.uniform(0, 2*np.pi, (ROTATION_REALISATIONS, len(flat)))
    self_query_mismatch = 0
    for e in flat:
        _, index = e['tree'].query(e['xy'], workers=WORKERS)
        self_query_mismatch += int((e['missing'][index] != e['missing']).sum())
    rot = {}
    coverage = {}
    for footprint in ROTATION_FOOTPRINTS:
        base = {e['id']: parts(e, lab_null[e['id']], footprint, e['missing']) for e in flat}
        rot[(footprint, None)] = base
        for i, e in enumerate(flat):
            series, cover = [], []
            for r in range(ROTATION_REALISATIONS):
                mask, covered = donor_query(e, e['xy'], 1.0, thetas[r, i], footprint)
                series.append(parts(e, lab_null[e['id']], footprint, mask))
                cover.append(covered.mean())
            rot[(footprint, e['id'])] = series
            coverage[(footprint, e['id'])] = float(np.mean(cover))

    def sets(footprint, g, which):
        es = rotation_centres[g]
        table = (lambda e: rot[(footprint, None)][e['id']]) if which is None else (lambda e: rot[(footprint, e['id'])][which])
        return table(es[0]), [table(e) for e in es[1:]]

    def statistics(footprint, groups, which, kind, k):
        """Residual, composition and within-centre terms for a pooled set of groups (one group or the stack)."""
        lens_parts, random_parts = [], []
        for g in groups:
            lp, rs = sets(footprint, g, which)
            lens_parts.append(lp)
            random_parts.extend(rs)
        dl, dr = decompose(lens_parts, k, kind), decompose(random_parts, k, kind)
        if dl is None or dr is None:
            return None
        return dict(residual=(1+dl['total'])/(1+dr['total'])-1, composition=(1+dl['composition'])/(1+dr['composition'])-1,
                    within=(1+dl['within'])/(1+dr['within'])-1, lens_composition=dl['composition'], random_composition=dr['composition'])

    def expected_residual(footprint, groups, which, kind, k):
        """Expected-count area over pre-mask stack - 1, both over the groups that enter both stacks."""
        items = [sets(footprint, g, which) for g in groups]
        common = set(stack('expected', items, k, kind, 'pre', ())[1]) & set(stack('expected', items, k, kind, 'area', ())[1])
        pre = stack('expected', items, k, kind, 'pre', (), only=common)[0]
        area = stack('expected', items, k, kind, 'area', (), only=common)[0]
        return over(area, pre)

    def edge_groups(footprint):
        """Groups with a centre that has masked sources and rotated coverage below the flag."""
        return {g for g, es in rotation_centres.items() if any(e['missing'].any() and coverage[(footprint, e['id'])] < COVERAGE_FLAG for e in es)}

    rotation_rows, combined_rows, edge_record = [], [], {}
    for footprint in ROTATION_FOOTPRINTS:
        edge = edge_groups(footprint)
        edge_record[str(footprint)] = sorted(edge)
        core = [x for x in supported if x not in edge]
        for kind in ('n', 'w'):
            family = []
            for g in [*supported, 'STACK', 'STACK_core', 'STACK_expected', 'STACK_expected_core']:
                vectors = defaultdict(lambda: dict(real=[], null=[]))
                pool = core if g.endswith('_core') else supported
                members = [g] if not g.startswith('STACK') else pool
                edge_affected = any(x in edge for x in members)
                min_cov = min(coverage[(footprint, e['id'])] for x in members for e in rotation_centres[x] if e['missing'].any()) \
                    if any(e['missing'].any() for x in members for e in rotation_centres[x]) else None
                for k in range(NBINS):
                    entering = [x for x in members if all(p['a_fp'][k] > 0 for p in (sets(footprint, x, None)[0], add(sets(footprint, x, None)[1])))]
                    use = entering
                    if g.startswith('STACK_expected'):
                        real_stats = {'residual': expected_residual(footprint, use, None, kind, k)}
                        null_stats_list = [{'residual': expected_residual(footprint, use, i, kind, k)} for i in range(ROTATION_REALISATIONS)]
                    else:
                        real_stats = statistics(footprint, use, None, kind, k) if use else None
                        null_stats_list = [statistics(footprint, use, i, kind, k) if use else None for i in range(ROTATION_REALISATIONS)]
                    lp, rs = (add([sets(footprint, x, None)[0] for x in use]), add([p for x in use for p in sets(footprint, x, None)[1]])) if use else (None, None)
                    mask_l = 1-lp[kind+'_post'][k]/lp[kind+'_pre'][k] if use and lp[kind+'_pre'][k] else None
                    mask_r = 1-rs[kind+'_post'][k]/rs[kind+'_pre'][k] if use and rs[kind+'_pre'][k] else None
                    for statistic in ('residual', 'composition', 'within'):
                        if g.startswith('STACK_expected') and statistic != 'residual':
                            continue
                        real_value = real_stats.get(statistic) if real_stats else None
                        null_values = [s.get(statistic) if s else None for s in null_stats_list]
                        row = dict(footprint_arcsec=footprint, weighting=kind, lens_id=g, statistic=statistic, annulus=k,
                                   r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], stacked_groups=len(use),
                                   lens_masked_fraction=mask_l, random_masked_fraction=mask_r, real=real_value,
                                   real_lens_composition=real_stats.get('lens_composition') if real_stats else None,
                                   real_random_composition=real_stats.get('random_composition') if real_stats else None,
                                   **null_stats(real_value, null_values))
                        row.update(masked_sources_present=bool((mask_l or 0) > 0 or (mask_r or 0) > 0), edge_affected=edge_affected,
                                   min_rotated_coverage_masked_centres=min_cov)
                        rotation_rows.append(row)
                        if statistic in ('residual', 'within') and row['masked_sources_present'] and row['null_sd'] and not edge_affected:
                            family.append(row)
                        if real_value is not None and len([v for v in null_values if v is not None]) == ROTATION_REALISATIONS:
                            vectors[statistic]['real'].append(real_value)
                            vectors[statistic]['null'].append(null_values)
                for statistic, v in vectors.items():
                    if not v['real']:
                        continue
                    combined_rows.append(dict(footprint_arcsec=footprint, weighting=kind, lens_id=g, statistic=statistic,
                                              edge_affected=edge_affected, informational_only=statistic == 'composition',
                                              **combined_test(v['real'], np.asarray(v['null']).T)))
            holm(family)
    for r in rotation_rows:
        r.setdefault('holm_p', None)
        r.setdefault('family_size', None)

    centre_rows = []
    for footprint in ROTATION_FOOTPRINTS:
        family = []
        for g, es in rotation_centres.items():
            for e in es:
                p0 = rot[(footprint, None)][e['id']]
                for k in range(NBINS):
                    null = [unmasked_density_change(p, k) for p in rot[(footprint, e['id'])]]
                    fractions = [1-p['n_post'][k]/p['n_pre'][k] for p in rot[(footprint, e['id'])] if p['n_pre'][k] > 0]
                    row = dict(footprint_arcsec=footprint, center_id=e['id'], role=e['role'], lens_id=g, annulus=k,
                               r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], n_pre=int(p0['n_pre'][k]),
                               masked_count_fraction=1-p0['n_post'][k]/p0['n_pre'][k] if p0['n_pre'][k] else None,
                               rotated_masked_count_fraction_mean=float(np.mean(fractions)) if fractions else None,
                               rotated_source_coverage=coverage[(footprint, e['id'])],
                               edge_affected=coverage[(footprint, e['id'])] < COVERAGE_FLAG,
                               real_unmasked_density_minus_pre=unmasked_density_change(p0, k), **null_stats(unmasked_density_change(p0, k), null))
                    centre_rows.append(row)
                    if (row['masked_count_fraction'] or 0) >= 0.01 and not row['edge_affected']:
                        family.append(row)
        holm(family)
    for r in centre_rows:
        r.setdefault('holm_p', None)
        r.setdefault('family_size', None)

    # c. Boost normalisation on the lenses and on random centres used as pseudo-lenses (1-arcsec areas).
    fine = real[(1.0, fp)]
    lens_items = {g: (fine[lenses[g]['id']], [fine[r['id']] for r in nonempty[g]]) for g in lenses}
    pseudo = []
    for g in supported:
        rs = nonempty[g]
        if len(rs) >= 2:
            for r in rs:
                pseudo.append((r['id'], g, fine[r['id']], [fine[x['id']] for x in rs if x['id'] != r['id']]))
    overlap_rows = []
    for g in supported:
        es = [lenses[g], *nonempty[g]]
        id_sets = {e['id']: set(e['ids']) for e in es}
        for i, a in enumerate(es):
            for b in es[i+1:]:
                shared = len(id_sets[a['id']] & id_sets[b['id']])
                overlap_rows.append(dict(lens_id=g, center_a=a['id'], role_a=a['role'], n_a=a['n'], center_b=b['id'], role_b=b['role'], n_b=b['n'],
                                         shared_sources=shared, shared_over_smaller=shared/min(a['n'], b['n'])))
    norm_rows, norm_summary = [], []
    variants = [(d, None) for d in ('pooled', 'expected')]+[(d, n) for d in ('outer_normalised', 'self_normalised') for n in NORMALISATION]
    pseudo_items = [(lp, rs) for _, _, lp, rs in pseudo]
    pseudo_groups = [g for _, g, _, _ in pseudo]
    for definition, norm_name in variants:
        norm = NORMALISATION.get(norm_name, ())
        for state in ('pre', 'area'):
            for kind in ('n', 'w'):
                for k in range(NBINS):
                    base = dict(definition=definition, normalisation_annuli=norm_name, state=state, weighting=kind, annulus=k,
                                r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], inside_normalisation=k in norm)
                    values, paired, no_own, zero_post = [], [], 0, 0
                    for cid, g, lp, rs in pseudo:
                        v = single(definition, lp, rs, k, kind, state, norm, OWN)
                        se = poisson_se(definition, lp, rs, k, norm, kind, state, OWN)
                        no_own += lp['n_pre'][k] <= 0 and lp['a_fp'][k] > 0
                        zero_post += v == 0
                        norm_rows.append(dict(**base, sample='pseudo_lens', center_id=cid, lens_id=g, boost=v, poisson_relative_se=se))
                        if v is not None:
                            values.append(v-1)
                            if se is not None:
                                paired.append((v-1, se))
                    for g, (lp, rs) in lens_items.items():
                        if definition != 'self_normalised' and g not in supported:
                            continue
                        norm_rows.append(dict(**base, sample='lens', center_id=lenses[g]['id'], lens_id=g,
                                              boost=single(definition, lp, rs, k, kind, state, norm, OWN),
                                              poisson_relative_se=poisson_se(definition, lp, rs, k, norm, kind, state, OWN)))
                    v = np.asarray(values, float)
                    pv, pse = (np.asarray([x[i] for x in paired], float) for i in (0, 1))
                    out = dict(**base, pseudo_n=len(v), pseudo_excluded_no_own_source=int(no_own), pseudo_zero_post_count=int(zero_post),
                               pseudo_mean_minus_1=float(v.mean()) if len(v) else None,
                               pseudo_sd=float(v.std(ddof=1)) if len(v) > 1 else None,
                               pseudo_robust_sd=float(1.4826*np.median(np.abs(v-np.median(v)))) if len(v) else None,
                               pseudo_n_with_poisson=len(pv), pseudo_rms_poisson_se=float(np.sqrt(np.mean(pse**2))) if len(pse) else None)
                    rms = out['pseudo_rms_poisson_se']
                    # Ratios over one set: centres with a defined Poisson SE.
                    out['pseudo_sd_over_poisson'] = float(pv.std(ddof=1))/rms if len(pv) > 1 and rms else None
                    out['pseudo_robust_sd_over_poisson'] = float(1.4826*np.median(np.abs(pv-np.median(pv))))/rms if len(pv) and rms else None
                    if definition == 'self_normalised':
                        value, used, se, _ = stack_jackknife(definition, pseudo_items, k, kind, state, norm, pseudo_groups, OWN)
                        out.update(pseudo_stack_minus_1=over(value, 1.0), pseudo_stack_group_jackknife_se=se, pseudo_stack_groups=len({pseudo_groups[i] for i in used}))
                    for name, subset in (('supported', supported), ('all', list(lenses))):
                        if name == 'all' and definition != 'self_normalised':
                            continue
                        items = [lens_items[g] for g in subset]
                        value, used, se, replicas = stack_jackknife(definition, items, k, kind, state, norm, own=OWN)
                        out[f'lens_stack_{name}_minus_1'] = over(value, 1.0)
                        out[f'lens_stack_{name}_jackknife_se'] = se
                        out[f'lens_stack_{name}_n'] = len(used)
                        if name == 'supported' and definition != 'expected':
                            ref, ref_used = stack('expected', items, k, kind, state, (), OWN)
                            if value is not None and ref is not None and set(used) == set(ref_used):
                                diffs = []
                                for d in used:
                                    rest = [it for j, it in enumerate(items) if j != d]
                                    a_ = over(stack(definition, rest, k, kind, state, norm, OWN)[0], 1.0)
                                    b_ = over(stack('expected', rest, k, kind, state, (), OWN)[0], 1.0)
                                    diffs.append(a_-b_ if a_ is not None and b_ is not None else None)
                                out['lens_stack_minus_expected'] = value-ref
                                out['lens_stack_minus_expected_jackknife_se'] = jackknife_strict(diffs)
                    norm_summary.append(out)
    outer_excess = {}
    for norm_name, norm in NORMALISATION.items():
        for state in ('pre', 'area'):
            def stacked_norm(items):
                """Expected-count stack over the normalisation annuli; returns (value - 1, entering indices)."""
                n = e = 0.0
                entering = []
                for i, (lp, rs) in enumerate(items):
                    rho = density(rs, 'n', state, norm, OWN)
                    count, area = (sum(counts(lp, 'n', state)[j][k] for k in norm if not OWN or lp['n_pre'][k] > 0) for j in (0, 1))
                    if positive(rho) and area > 0:
                        entering.append(i)
                        n, e = n+count, e+rho*area
                return (n/e-1 if e > 0 else None), entering
            lens_list = [lens_items[g] for g in supported]
            value, entering = stacked_norm(lens_list)
            outer_excess[f'{norm_name}:{state}'] = dict(
                lens_expected_stack_minus_1=value, lenses_entering=len(entering),
                lens_jackknife_se=jackknife_strict([stacked_norm([it for j, it in enumerate(lens_list) if j != d])[0] for d in entering]))

    output.mkdir(parents=True)
    tables = {'grid_resolution.tsv': grid_rows, 'transplant_field_realisations.tsv': transplant_rows,
              'transplant_field_summary.tsv': transplant_summary, 'rotation_null_groups.tsv': rotation_rows,
              'rotation_null_combined.tsv': combined_rows, 'rotation_null_centres.tsv': centre_rows,
              'normalisation_centres.tsv': norm_rows, 'normalisation_summary.tsv': norm_summary, 'centre_overlap.tsv': overlap_rows}
    for name, rows in tables.items():
        write_table(output/name, rows, list(dict.fromkeys(k for r in rows for k in r)))
    pairs = [r for r in overlap_rows if r['role_a'] == 'random' and r['role_b'] == 'random']
    report = dict(
        science_measurement=False, calibration_applied=False, selection_changed=False, revision=REVISION,
        constants=dict(grids_arcsec=GRIDS_ARCSEC, null_grid_arcsec=NULL_GRID_ARCSEC, anchor_grid_arcsec=ANCHOR_GRID_ARCSEC, own_source_rule_in_normalisation=OWN, primary_footprint_arcsec=fp,
                       transplant_realisations=TRANSPLANT_REALISATIONS, rotation_realisations=ROTATION_REALISATIONS,
                       rotation_footprints_arcsec=ROTATION_FOOTPRINTS, coverage_flag=COVERAGE_FLAG,
                       normalisation_annuli={k: list(v) for k, v in NORMALISATION.items()}, seed=SEED, donor_lens=DONOR_LENS,
                       radial_edges_mpc=RADIAL_EDGES.tolist()),
        supported_lens_groups=supported, nonempty_randoms_per_group={g: len(nonempty[g]) for g in supported},
        transplant_recipients=recipients, donor_randoms=[r['id'] for r in donor_randoms],
        pseudo_lens_centres=len(pseudo), pseudo_lens_groups=sorted(set(pseudo_groups)),
        reproduction_4arcsec_table=reproduction, reproduction_lens_side_transplant=anchor,
        rotation_edge_affected_groups=edge_record, rotation_self_query_label_mismatches=self_query_mismatch, max_tangent_radius_error_mpc=radius_error,
        centre_overlap=dict(lens_random_pairs_sharing=sum(1 for r in overlap_rows if 'lens' in (r['role_a'], r['role_b']) and r['shared_sources']),
                            random_pairs=len(pairs), random_pairs_sharing=sum(1 for r in pairs if r['shared_sources']),
                            random_pairs_sharing_at_least_half=sum(1 for r in pairs if r['shared_over_smaller'] >= 0.5),
                            max_shared_over_smaller=max((r['shared_over_smaller'] for r in pairs), default=None)),
        outer_normalisation_excess=outer_excess,
        definitions=dict(
            mask='Selected source lacks Ephor_AB best/risk (identical to the afterburner NaN mask, verified by the prior runner).',
            area='Grid points within the footprint of the nearest selected source of the same centre; unmasked if that source is unmasked. '
                 'In the geometric transplant mode the donor mask at each grid point is ORed with the recipient own labels.',
            transplant='Recipient positions rotated by theta and scaled to the donor frame in proper Mpc; a position is masked when the '
                       'donor nearest source within the footprint is masked. Realisation 0: no rotation, recipient random i takes donor '
                       'random i (sorted ids); null realisations: independent uniform rotation per centre and a random donor-random '
                       'assignment. The truth is the recipient boost before the transplant (same estimator).',
            rotation_null='Each centre of a lens group with randoms keeps its sources; its own mask is rotated by an independent uniform '
                          'angle (the same angles for every footprint) and applied with nearest-source areas. Rotated positions with no '
                          'own source within the footprint are unmasked; rotated_source_coverage records that fraction. '
                          'residual = area-normalised over pre-mask boost - 1 = (1+composition)(1+within) - 1 exactly, where for each '
                          'side composition reweights the pooled pre-mask densities by unmasked area and within is the unmasked '
                          'density over that baseline. STACK pools the groups (lens side too); STACK_expected uses the expected-count '
                          'stack. Combined tests: sum of z^2 and Mahalanobis over annuli with an empirical p from the same null.',
            normalisation=dict(
                pooled='(sum lens n / sum lens area) / (sum random n / sum random area), randoms pooled across groups (recorded definition).',
                expected='sum lens n / sum over groups of (own random density x lens area); identical to pooled for one lens.',
                outer_normalised='expected-count boost with each group random density rescaled by its lens/random density ratio in the '
                                 'normalisation annuli.',
                self_normalised='sum lens n / sum over lenses of (own lens density in the normalisation annuli x lens area); no randoms.',
                pseudo_lens='Each non-empty random of a group with at least two, against the other non-empty randoms of the same group; '
                            'the truth is 1. Stacks of pooled, expected and outer-normalised pseudo-lenses are pinned near 1 by '
                            'leave-one-out closure and are not reported.',
                inclusion='One rule for per-centre and stacked values: a centre enters annulus k only with at least one own pre-mask '
                          'source there (own-source rule; dilated label area alone does not count), and wherever its areas and '
                          'reference densities are positive; a zero post-mask count over positive area gives boost 0.')),
        limitations=['Transplant and rotation nulls make the mask independent of the recipient density by construction; they test '
                     'the estimator, not whether the real mask is redshift-independent.',
                     'One heavily masked group (S23B_3620) carries the rotation test and supplies the only donor; its two ~50% '
                     'masked randoms share 59% of their sources, so there are about three independent mask patterns.',
                     'Rotated masks can fall outside the footprint (see rotated_source_coverage; flagged below the coverage flag).',
                     'Random centres of a group often share sources (centre_overlap.tsv); pseudo-lens scatter is then a lower bound '
                     'on field variance and rotations of overlapping centres are independent although their sky masks coincide.',
                     'Normalised definitions remove any real excess in the normalisation annuli (members, two-halo term); the size '
                     'of that excess is not measured here.',
                     'Poisson errors are delta-method count errors (effective counts for weights) and ignore clustering and reuse; '
                     '8 groups give noisy jackknife errors.',
                     'Transplant and rotation nulls use a 2-arcsec grid; the grid and normalisation tables use 1 arcsec; the anchors 4 arcsec.',
                     'The 30-arcsec nearest-source footprint dilates label area across annulus edges, holes and survey edges; the '
                     'grid, transplant and rotation sections keep this recorded convention, only the normalisation section applies '
                     'the own-source rule.',
                     'Only one donor and one masked field exist; no second masked field is available locally.'],
        provenance=dict(**provenance, geometry_sha256=GEOMETRY_SHA, prior_summary_sha256=PRIOR_SUMMARY_SHA, prior_runner_sha256=PRIOR_RUNNER_SHA),
        runner_sha256=sha(Path(__file__)), output_sha256={name: sha(output/name) for name in tables})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(reproduction=reproduction, anchor=anchor, supported=supported, pseudo=len(pseudo),
                          self_query_mismatch=self_query_mismatch), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
