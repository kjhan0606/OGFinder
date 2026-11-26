#!/usr/bin/env python3
"""Cause of the 0-0.4 Mpc deficit of clean background sources around the S16A lenses; local only.

Clean background: standard sources with photoz_best - z_l > 0.2 and Ephor_AB err95_min > z_l (receipt be6bd1a9).
Setup as audit_s16a_pz_dilution (standard areas, 76 new randoms, expected count, group jackknife). Fixed before run:
A. Magnitude: boost per dereddened i bin (<=23, 23-24, 24-24.5) per annulus and pooled 0-0.4 Mpc. Magnification
   prediction from the pooled random clean counts: n'(bin) = n(m1+d, m2+d)/mu with d = 2.5 log10 mu, for mu = 1.05,
   1.1, 1.2 (flux magnification only; size-dependent selection is ignored), and the count slope s at 24.5.
B. Crowding: mean blendedness and fraction with blendedness > 0.05 of clean sources, and fraction of clean sources
   with another shape-table source within 3 arcsec, lens and random per annulus (pooled 0-0.4 Mpc too).
C. Per lens (groups with 0-0.4 Mpc area): Spearman of the clean-background 0-0.4 Mpc boost against the member-like
   (red tertile, photoz_best - z_l <= 0.1) 0-0.4 Mpc boost.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_afterburner_mask_effects import tangent_arcsec
from audit_s16a_background_cut_boosts import masked
from audit_s16a_fdfc_areas import VALIDATION, in_fdfc
from audit_s16a_new_inner_point_areas import areas, new_inner_rows, with_area
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import label_fdfc, load_new, load_new_centres
from audit_s16a_random_point_areas import INNER, inner_ok, i_null, load_rows, unit
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

VARIANT = 'null_duplicates + fdfc'
MAG_EDGES = (-np.inf, 23., 24., 24.5)
MUS = (1.05, 1.1, 1.2)
NEIGHBOUR_ARCSEC = 3.0
BLEND = 0.05
POOL = (0, 1, 2)


def pooled(p):
    """Parts with annulus 2 replaced by the sum over 0-0.4 Mpc (used only as 'annulus 2')."""
    q = {k: v.copy() for k, v in p.items()}
    for key in bt.KEYS:
        q[key][2] = p[key][list(POOL)].sum()
    return q


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
    parts, entries = {}, {}
    for g in groups:
        lens, grid = lens_row[g], fdfc['grids'][g]
        fh = (lambda ra, dec, lens=lens, grid=grid: ~in_fdfc(lens, ra, dec, grid))
        lid = lens['center_id']
        for cid, e, c in [(lid, pilot[lid], lens), *((r, new[r], centres[r]) for r in new_rand[g])]:
            parts[cid] = with_area(label_fdfc(e, lens, grid), areas(c, float(e['radius']), by_c[cid], tree, fh)[VARIANT])
            entries[cid] = e
    pq, nq = pilot_quality(), new_quality()
    pz, pz_sha = load_pz()
    arr = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        arr[cid] = dict(best=g('photoz_best', q)-zl, lo95=g('ephor_ab_photoz_err95_min', pz)-zl,
                        imag=g('i_mag', q)-g('a_i', q), blend=g('blendedness', q),
                        gi=np.array([colour(q[o]) for o in e['ids']]),
                        nn=(e['tree'].query(e['xy'], k=2)[0][:, 1] if e['n'] > 1 else np.full(e['n'], np.inf)))
    clean = {cid: (a['best'] > .2) & (a['lo95'] > 0) & ~entries[cid]['missing'] for cid, a in arr.items()}
    rand_ids = [r for g in groups for r in new_rand[g]]
    lens_ids = [lens_row[g]['center_id'] for g in groups]

    def stack_rows(label, mask_fn):
        src = {cid: masked(parts[cid], entries[cid], mask_fn(cid)) for cid in parts}
        out = []
        for name, s in (('annulus', src), ('pooled 0-0.4', {c: pooled(p) for c, p in src.items()})):
            items = [(s[lens_row[g]['center_id']], [s[r] for r in new_rand[g]]) for g in groups]
            for k in (range(bt.NBINS) if name == 'annulus' else (2,)):
                value, used, se, _ = bt.stack_jackknife('expected', items, k, 'n', 'area', (), None, False)
                out.append(dict(subset=label, bin=name, annulus=k, r_lo=RADIAL_EDGES[k] if name == 'annulus' else 0.0, r_hi=RADIAL_EDGES[k+1],
                                boost=value, group_jackknife_se=se, lens_groups=len(used),
                                n_lens=float(sum(lp['n_post'][k] for lp, _ in items))))
        return out
    mag = stack_rows('clean all', lambda c: clean[c])
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        mag += stack_rows(f'clean i_dered ({lo},{hi}]', lambda c, lo=lo, hi=hi: clean[c] & (arr[c]['imag'] > lo) & (arr[c]['imag'] <= hi))

    m = np.sort(np.concatenate([arr[c]['imag'][clean[c]] for c in rand_ids]))
    cum = lambda x: np.searchsorted(m, x, side='right')
    pred = []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        base = cum(hi)-cum(lo)
        row = dict(bin=f'({lo},{hi}]', n_random=int(base))
        for mu in MUS:
            d = 2.5*np.log10(mu)
            row[f'pred_mu_{mu}'] = float((cum(min(hi+d, 24.5+d) if hi == 24.5 else hi+d)-cum(lo+d))/mu/base) if base else None
        pred.append(row)
    total = dict(bin='all (<=24.5)', n_random=int(cum(24.5)))
    for mu in MUS:
        total[f'pred_mu_{mu}'] = float(cum(24.5+2.5*np.log10(mu))/mu/cum(24.5))
    pred.append(total)
    slope = float((np.log10(cum(24.5)) - np.log10(cum(24.0)))/0.5)

    crowd = []
    for role, cids in (('lens', lens_ids), ('random', rand_ids)):
        for k, ks in [*((k, (k,)) for k in range(bt.NBINS)), ('0-0.4', POOL)]:
            sel_b, sel_n = [], []
            for c in cids:
                s = clean[c] & np.isin(entries[c]['bins'], ks)
                sel_b.append(arr[c]['blend'][s]); sel_n.append(arr[c]['nn'][s])
            b, nn = np.concatenate(sel_b), np.concatenate(sel_n)
            crowd.append(dict(role=role, annulus=str(k), n=len(b), mean_blendedness=float(np.nanmean(b)) if len(b) else None,
                              frac_blend_gt=float(np.mean(b > BLEND)) if len(b) else None,
                              frac_neighbour_within=float(np.mean(nn <= NEIGHBOUR_ARCSEC)) if len(nn) else None))
    t2 = float(np.nanpercentile(np.concatenate([arr[c]['gi'][~entries[c]['missing']] for c in rand_ids]), 200/3))
    member = {c: (arr[c]['gi'] > t2) & (arr[c]['best'] <= .1) & ~entries[c]['missing'] for c in arr}
    per = []
    for g in groups:
        it = lambda mk: [(pooled(masked(parts[lens_row[g]['center_id']], entries[lens_row[g]['center_id']], mk[lens_row[g]['center_id']])),
                          [pooled(masked(parts[r], entries[r], mk[r])) for r in new_rand[g]])]
        b_clean = bt.stack('expected', it(clean), 2, 'n', 'area', ())[0]
        b_mem = bt.stack('expected', it(member), 2, 'n', 'area', ())[0]
        if b_clean is not None and b_mem is not None:
            per.append(dict(lens_id=g, clean_boost_0_0p4=b_clean, member_boost_0_0p4=b_mem))
    rho = spearmanr([p['clean_boost_0_0p4'] for p in per], [p['member_boost_0_0p4'] for p in per])
    output.mkdir(parents=True)
    tables = {'magnitude_boosts.tsv': mag, 'magnification_prediction.tsv': pred, 'crowding.tsv': crowd, 'per_lens.tsv': per}
    for name, tab in tables.items():
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in tab for k in r)), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-clean-deficit-v1', science_measurement=False, selection_change=False, new_query=False,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check,
                   count_slope_s_at_24p5=slope, red_tertile=t2, per_lens_spearman=dict(rho=float(rho.statistic), p=float(rho.pvalue), n=len(per)),
                   limitations=['Magnification prediction uses flux only; resolution/size selection response is ignored.',
                                'Neighbours are shape-table sources only, not all detections.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in tables})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return tables, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_clean_deficit_20261126_01'))
    tables, s = run(parser.parse_args().output)
    f = lambda x: '—' if x is None else f'{x:.3f}'
    for r in tables['magnitude_boosts.tsv']:
        if r['bin'] == 'pooled 0-0.4' or r['annulus'] in (0, 1, 2, 3, 7):
            print(f"{r['subset']:28s} {r['bin'][:6]} {r['r_lo']}-{r['r_hi']:<4} B={f(r['boost'])} ± {f(r['group_jackknife_se'])} N={r['n_lens']:.0f} g={r['lens_groups']}")
    for t in ('magnification_prediction.tsv', 'crowding.tsv', 'per_lens.tsv'):
        for r in tables[t]:
            print(t[:5], {k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print(json.dumps({k: s[k] for k in ('count_slope_s_at_24p5', 'per_lens_spearman')}))
