#!/usr/bin/env python3
"""Axis convention of the S16A blendedness moments and the tangential shear bias they allow. Local only (receipt 6851a9e9).

Fixed before the run. Four candidate maps of the i blendedness moments (Q11 - Q22, 2 Q12)/(Q11 + Q22) onto the shape
e1/e2 frame: (s1, s2) in {+1, -1}^2 multiplying the two components.
A1. Child moments versus the catalogue shape: weighted correlation of s1*m1 with e1 and of s2*m2 with e2 over all selected
    sources with moments. The right map gives two positive correlations.
A2. Sibling geometry: children (parent_id != 0) with exactly one other member of the same parent in the same disc. The
    extra light (parent - child moments) must be elongated toward the sibling: under the right map the tangential
    component about the sibling position (ogfkit tangential, xy in arcsec, x east, y north) is negative.
B. With the map fixed by A, lens - random (group jackknife) in angular and proper-radius bins, for all and blue primary
   sources, of the weighted tangential extra-light ellipticity et_extra and the catalogue tangential ellipticity et.
   et_extra is the shift of the parent-footprint distortion over the child one; the deblended child shape need not carry
   it, so et_extra / et is an upper bound on the fraction of the lens signal that this light could induce.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_local_background_abs_e import load_lb
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wlcluster import tangential  # noqa: E402

THETA = (0, 30, 60, 90, 150, 300, 600)
RPROP = (0, 0.2, 0.4, 0.7, 1.0, 1.5, 3.0)
MAPS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centres, selection_sha = load_new_centres()
    new, _, provenance, _ = load_new(centres)
    pilot_centres, pilot, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in pilot_centres if c['role'] == 'lens'}
    groups = sorted(lens_row)
    new_rand = defaultdict(list)
    for cid, c in centres.items():
        new_rand[c['matched_lens_id']].append(cid)
    _, ppairs, _, _, _, _ = load_inputs(PLAN, OBJECT_QA)
    res, summ = receipt(SOURCE_PLAN, SOURCE_SHA)
    npairs = [{k: r[k] for k in PAIR_FIELDS} for r in rows(res, sorted(n for n in summ['output_sha256'] if 'pairs' in n))]
    pair = {(p['center_id'], p['object_id']): p for p in ppairs+npairs}
    pq, nq = pilot_quality(), new_quality()
    pz, pz_sha = load_pz()
    lb, lb_sha = load_lb()
    lens_ids = {g: lens_row[g]['center_id'] for g in groups}
    entries = {**{lens_ids[g]: pilot[lens_ids[g]] for g in groups}, **{r: new[r] for g in groups for r in new_rand[g]}}
    data = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        keep = ~e['missing']
        prim = ((g('photoz_best', q)-zl > .2) & (g('ephor_ab_photoz_err95_min', pz)-zl > 0) & (g('i_mag', q)-g('a_i', q) > 23))[keep]
        ids = np.array(e['ids'])[keep]
        P = [pair[(cid, o)] for o in ids]
        e1, e2, w = (np.array([num(p[k]) for p in P]) for k in ('e1', 'e2', 'weight'))
        L = lambda key: np.array([num(lb[o][key]) for o in ids])
        c11, c22, c12 = (L('i_blendedness_abs_shape_child_'+k) for k in ('11', '22', '12'))
        d11, d22, d12 = (L('i_blendedness_abs_shape_parent_'+k)-v for k, v in (('11', c11), ('22', c22), ('12', c12)))
        tr = c11+c22
        ok = tr > 0
        xy = e['xy'][keep]
        data[cid] = dict(ids=ids, w=w, e1=e1, e2=e2, x=xy[:, 0], y=xy[:, 1], theta=np.hypot(*xy.T),
                         rprop=np.hypot(*xy.T)*3.0/float(e['radius']), prim=prim, gi=np.array([colour(q[o]) for o in ids]),
                         m1=np.where(ok, (c11-c22)/np.where(ok, tr, 1), np.nan), m2=np.where(ok, 2*c12/np.where(ok, tr, 1), np.nan),
                         d1=np.where(ok, (d11-d22)/np.where(ok, tr, 1), np.nan), d2=np.where(ok, 2*d12/np.where(ok, tr, 1), np.nan),
                         pid=np.array([lb[o]['parent_id'].strip() for o in ids]))
    allc = list(data)

    def wcorr(a, b, w):
        ok = np.isfinite(a) & np.isfinite(b) & (w > 0)
        a, b, w = a[ok], b[ok], w[ok]
        am, bm = np.average(a, weights=w), np.average(b, weights=w)
        return float(np.average((a-am)*(b-bm), weights=w)/np.sqrt(np.average((a-am)**2, weights=w)*np.average((b-bm)**2, weights=w)))
    cat = lambda cids, k: np.concatenate([data[c][k] for c in cids])
    W, E1, E2, M1, M2 = (cat(allc, k) for k in ('w', 'e1', 'e2', 'm1', 'm2'))
    a1 = [dict(s1=s1, s2=s2, corr_e1=wcorr(s1*M1, E1, W), corr_e2=wcorr(s2*M2, E2, W), n=int(np.isfinite(M1).sum())) for s1, s2 in MAPS]
    sib = defaultdict(list)
    for cid, d in data.items():
        fam = defaultdict(list)
        for i, p in enumerate(d['pid']):
            if p not in ('', '0'):
                fam[p].append(i)
        for idx in fam.values():
            if len(idx) == 2:
                for i, j in (idx, idx[::-1]):
                    if np.isfinite(d['d1'][i]):
                        sib['dx'].append(d['x'][j]-d['x'][i]); sib['dy'].append(d['y'][j]-d['y'][i])
                        sib['d1'].append(d['d1'][i]); sib['d2'].append(d['d2'][i])
    S = {k: np.array(v) for k, v in sib.items()}
    sep = np.hypot(S['dx'], S['dy'])
    a2 = []
    for s1, s2 in MAPS:
        et, ex = tangential(S['dx'], S['dy'], s1*S['d1'], s2*S['d2'])
        rad = -et/np.maximum(np.hypot(S['d1'], S['d2']), 1e-12)
        a2.append(dict(s1=s1, s2=s2, n=len(et), median_sep_arcsec=float(np.median(sep)), mean_et=float(et.mean()),
                       se_et=float(et.std(ddof=1)/np.sqrt(len(et))), mean_ex=float(ex.mean()), se_ex=float(ex.std(ddof=1)/np.sqrt(len(ex))),
                       frac_radial=float((rad > 0).mean()), mean_cos2_alignment=float(rad.mean())))
    win1 = [r for r in a1 if r['corr_e1'] > 0 and r['corr_e2'] > 0]
    win2 = min(a2, key=lambda r: r['mean_et'])
    agree = len(win1) == 1 and (win1[0]['s1'], win1[0]['s2']) == (win2['s1'], win2['s2'])
    s1, s2 = (win1[0]['s1'], win1[0]['s2']) if agree else (1, 1)
    for d in data.values():
        d['et'], d['ex'] = tangential(d['x'], d['y'], d['e1'], d['e2'])
        ee1, ee2 = s1*d['d1'], s2*d['d2']
        et, ex = tangential(d['x'], d['y'], np.nan_to_num(ee1), np.nan_to_num(ee2))
        d['et_extra'], d['ex_extra'] = np.where(np.isfinite(ee1), et, np.nan), np.where(np.isfinite(ee1), ex, np.nan)
    rand_cids = [r for g in groups for r in new_rand[g]]
    t1 = float(np.nanpercentile(cat(rand_cids, 'gi')[cat(rand_cids, 'prim')], 100/3))
    sel = {'all': lambda d: d['prim'], 'blue': lambda d: d['prim'] & (d['gi'] <= t1)}

    def wm(cids, key, var, lo, hi, cname):
        num_ = den = 0.0; n = 0
        for c in cids:
            d = data[c]
            s = sel[cname](d) & (d[var] >= lo) & (d[var] < hi) & np.isfinite(d[key]) & np.isfinite(d['d1'])
            num_ += float((d['w'][s]*d[key][s]).sum()); den += float(d['w'][s].sum()); n += int(s.sum())
        return (num_/den if den > 0 else np.nan), n

    def diff(gs, key, var, lo, hi, cname):
        l, n = wm([lens_ids[g] for g in gs], key, var, lo, hi, cname)
        r, _ = wm([r for g in gs for r in new_rand[g]], key, var, lo, hi, cname)
        return l, r, n

    def jk(fn):
        reps = np.array([fn([g for g in groups if g != d]) for d in groups])
        reps = reps[np.isfinite(reps)]
        return float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum()))
    out = []
    for var, edges in (('theta', THETA), ('rprop', RPROP)):
        for cname in sel:
            for lo, hi in zip(edges[:-1], edges[1:]):
                row = dict(binning=var, lo=lo, hi=hi, colour=cname)
                for key in ('et', 'et_extra', 'ex_extra'):
                    l, r, n = diff(groups, key, var, lo, hi, cname)
                    row.update({f'{key}_lens': l, f'{key}_random': r, f'{key}_diff': l-r,
                                f'{key}_se': jk(lambda gs: np.subtract(*diff(gs, key, var, lo, hi, cname)[:2])), 'n_lens': n})
                row['ratio_extra_over_et'] = row['et_extra_diff']/row['et_diff'] if row['et_diff'] else np.nan
                row['ratio_se'] = jk(lambda gs: (lambda a, b: (a[0]-a[1])/(b[0]-b[1]))(diff(gs, 'et_extra', var, lo, hi, cname),
                                                                                       diff(gs, 'et', var, lo, hi, cname)))
                out.append(row)
    output.mkdir(parents=True)
    with (output/'bias_profiles.tsv').open('w', newline='') as stream:
        wr = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        wr.writeheader(); wr.writerows(out)
    summary = dict(schema='astrafek-s16a-extra-light-axes-v1', science_measurement=False, selection_change=False, new_query=False,
                   a1_child_moment_vs_shape=a1, a2_sibling_alignment=a2, maps_agree=agree, chosen_map=dict(s1=s1, s2=s2),
                   blue_gi_max=t1, theta_edges_arcsec=THETA, rprop_edges_mpc=RPROP,
                   local_background_receipt_summary_sha256=lb_sha, pz_receipt_summary_sha256=pz_sha,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, groups=groups,
                   limitations=['et_extra is the parent-footprint distortion shift; the deblended child shape may carry little of it, '
                                'so the ratio is an upper bound.',
                                'Smooth light removed by sky subtraction is not in the parent moments.',
                                'Sibling test uses only families with exactly two members in the selected set.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'bias_profiles.tsv': sha(output/'bias_profiles.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return a1, a2, agree, (s1, s2), out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_extra_light_axes_20261127_01'))
    a1, a2, agree, m, out = run(parser.parse_args().output)
    print(json.dumps(dict(a1=a1, a2=a2, agree=agree, map=m), indent=1, default=float))
    f = lambda v: f'{v:+.4f}'
    for r in out:
        print(f"{r['binning']:5s} {r['lo']:>5}-{r['hi']:<5} {r['colour']:4s} n={r['n_lens']:5d} "
              f"et {f(r['et_diff'])}±{r['et_se']:.4f}  et_x {f(r['et_extra_diff'])}±{r['et_extra_se']:.4f}  "
              f"ex_x {f(r['ex_extra_diff'])}±{r['ex_extra_se']:.4f}  ratio {r['ratio_extra_over_et']:+.3f}±{r['ratio_se']:.3f}")
