#!/usr/bin/env python3
"""Candidate (c): is the S16A central source loss shape dependent? Boost of primary sources per |e| class. Local only.

The Y1 cuts remove only 6 of 37,029 lens shape-table rows, so the lost sources have no shapes; the test is indirect.
Fixed before the run: primary sources (audit_s16a_primary_boost) split by raw |e| = hypot(e1, e2) tertiles of the
pooled random primary sources, for all colours and for the blue g-i tertile (as audit_s16a_round_member_test).
Boost per class with the primary-boost areas, randoms, stack and area noise. A shape-independent loss gives equal
boosts in all |e| classes; a loss of elongated sources gives a lower boost in the high-|e| class.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_background_cut_boosts import masked
from audit_s16a_fdfc_areas import VALIDATION, in_fdfc
from audit_s16a_new_inner_point_areas import areas, new_inner_rows, with_area
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, label_fdfc, load_new, load_new_centres, receipt, rows
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from audit_s16a_random_point_areas import INNER, inner_ok, i_null, load_rows, unit
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

VARIANT = 'null_duplicates + fdfc'
REALISATIONS = 1000
SEED = 2026112603
OUTLIER = 'S23B_9307'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fdfc = load_fdfc()
    acc = json.loads((VALIDATION/'user_acceptance.json').read_text())
    if not (acc['accepted'] and acc['grid_summary_sha256'] == fdfc['summary_sha256']):
        raise ValueError('FDFC grid not accepted')
    wide_check = merge_wide_fdfc(fdfc['grids'], load_wide())
    centres, selection_sha = load_new_centres()
    new, counts, provenance, _ = load_new(centres)
    pilot_centres, pilot, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in pilot_centres if c['role'] == 'lens'}
    groups = sorted(lens_row)
    new_rand = defaultdict(list)
    for cid, c in centres.items():
        new_rand[c['matched_lens_id']].append(cid)
    _, prow = load_rows(INNER)
    _, nrow = new_inner_rows()
    by_c = defaultdict(list)
    for r in prow+nrow:
        by_c[r['center_id']].append(r)
    pool, seen = [], set()
    for r in prow+nrow:
        if not i_null(r) and inner_ok(r) and r['object_id'] not in seen:
            seen.add(r['object_id']); pool.append((float(r['ra']), float(r['dec'])))
    tree = cKDTree(unit(*map(np.array, zip(*pool))))
    base, entries, rp = {}, {}, {}
    for g in groups:
        lens, grid = lens_row[g], fdfc['grids'][g]
        fh = (lambda ra, dec, lens=lens, grid=grid: ~in_fdfc(lens, ra, dec, grid))
        lid = lens['center_id']
        for cid, e, c in [(lid, pilot[lid], lens), *((r, new[r], centres[r]) for r in new_rand[g])]:
            base[cid], entries[cid] = label_fdfc(e, lens, grid), e
            rp[cid] = areas(c, float(e['radius']), by_c[cid], tree, fh)[VARIANT]
    pq, nq = pilot_quality(), new_quality()
    pz, pz_sha = load_pz()
    _, ppairs, _, _, _, _ = load_inputs(PLAN, OBJECT_QA)
    res, summ = receipt(SOURCE_PLAN, SOURCE_SHA)
    npairs = [{k: r[k] for k in PAIR_FIELDS} for r in rows(res, sorted(n for n in summ['output_sha256'] if 'pairs' in n))]
    pair = {(p['center_id'], p['object_id']): p for p in ppairs+npairs}
    props = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        best, lo95, imag = g('photoz_best', q)-zl, g('ephor_ab_photoz_err95_min', pz)-zl, g('i_mag', q)-g('a_i', q)
        ae = np.array([np.hypot(num(pair[(cid, o)]['e1']), num(pair[(cid, o)]['e2'])) for o in e['ids']])
        props[cid] = dict(primary=(best > .2) & (lo95 > 0) & (imag > 23), abs_e=ae, gi=np.array([colour(q[o]) for o in e['ids']]),
                          ok=~e['missing'])
    rand_cids = [r for g in groups for r in new_rand[g]]
    t1 = float(np.nanpercentile(np.concatenate([props[c]['gi'][props[c]['ok']] for c in rand_cids]), 100/3))
    e1t, e2t = np.nanpercentile(np.concatenate([props[c]['abs_e'][props[c]['ok'] & props[c]['primary']] for c in rand_cids]), [100/3, 200/3])
    masks = {}
    for cid, d in props.items():
        lo, mid, hi = d['abs_e'] <= e1t, (d['abs_e'] > e1t) & (d['abs_e'] <= e2t), d['abs_e'] > e2t
        blue = d['gi'] <= t1
        masks[cid] = {'primary all': d['primary'], 'primary |e| low': d['primary'] & lo, 'primary |e| mid': d['primary'] & mid,
                      'primary |e| high': d['primary'] & hi, 'primary blue': d['primary'] & blue,
                      'primary blue |e| low': d['primary'] & blue & lo, 'primary blue |e| mid': d['primary'] & blue & mid,
                      'primary blue |e| high': d['primary'] & blue & hi}
    rng = np.random.default_rng(SEED)
    draws = {cid: [(v['annulus'], v['a_geo'], v['n_present'], v['n_pass']) for v in vals] for cid, vals in rp.items()}
    realised = [{cid: [dict(annulus=k, area=(a*rng.binomial(n, m/n)/n if n else 0.0)) for k, a, n, m in d] for cid, d in draws.items()}
                for _ in range(REALISATIONS)]
    out = []
    for dname in masks[next(iter(masks))]:
        src = {cid: masked(base[cid], entries[cid], masks[cid][dname]) for cid in base}
        std = {cid: with_area(src[cid], rp[cid]) for cid in src}
        for gname, gs in (('all groups', groups),):
            items = lambda s: [(s[lens_row[g]['center_id']], [s[r] for r in new_rand[g]]) for g in gs]
            for kind in ('n',):
                for k in range(bt.NBINS):
                    value, used, jk, _ = bt.stack_jackknife('expected', items(std), k, kind, 'area', (), None, False)
                    mc = None
                    if k <= 2:
                        vals = [v for real in realised
                                if (v := bt.stack('expected', items({c: with_area(src[c], real[c]) for c in src}), k, kind, 'area', ())[0]) is not None]
                        mc = float(np.std(vals, ddof=1))
                    out.append(dict(definition=dname, groups=gname, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                    lens_groups=len(used), boost=value, group_jackknife_se=jk, area_binomial_se=mc,
                                    total_se=None if jk is None else float(np.hypot(jk, mc or 0.0)),
                                    n_lens=float(sum(lp['n_post'][k] for lp, _ in items(std)))))
    output.mkdir(parents=True)
    with (output/'abs_e_class_boosts.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)
    summary = dict(schema='astrafek-s16a-shape-dependent-loss-v1', science_measurement=False, new_query=False,
                   abs_e_tertiles=[float(e1t), float(e2t)], blue_gi_max=t1, realisations=REALISATIONS, seed=SEED,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, groups=groups, counts=counts,
                   limitations=['Indirect: the lost sources have no shapes; a lower high-|e| boost could also come from a lens-side |e| shift.',
                                'Raw |e| classes; tertile edges from the random primary sources.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'abs_e_class_boosts.tsv': sha(output/'abs_e_class_boosts.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_shape_dependent_loss_20261126_01'))
    f = lambda x: '—' if x is None else f'{x:.3f}'
    for r in run(parser.parse_args().output):
        if r['weighting'] == 'n':
            print(f"{r['definition']:22s} {r['r_lo']}-{r['r_hi']:<4} B={f(r['boost'])} ± {f(r['total_se'])} "
                  f"(area {f(r['area_binomial_se'])}) N={r['n_lens']:.0f} g={r['lens_groups']}")
