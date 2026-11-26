#!/usr/bin/env python3
"""Candidate (a): does centre (BCG/ICL) light explain the lower |e| of S16A primary sources near the lenses? Local only.

Primary sources as audit_s16a_primary_boost; centre = lens catalogue centre (assumed at the BCG). Fixed before the run:
A. Weighted |e| lens - random (group jackknife) in angular bins (arcsec) and in proper-radius bins, for all primary
   sources and the blue g-i tertile, over all lens groups and separately for low-z (z_l < 0.45) and high-z groups.
   Light predicts a deficit ordered in angle (stronger at low z, less surface-brightness dimming); cluster structure
   predicts one ordered in proper radius.
B. Local-noise proxy from the existing columns (no sky column was received): q_b = flux_err_b / median random
   flux_err_b in the same cell of dereddened i (0.25 mag) x raw_sigma decile, for b = g, i. Lens - random of the
   weighted mean log10 q_b by angular bin; fraction with any of pixel_interpolated_center, pixel_suspect_center,
   pixel_clipped_any; blendedness. Random |e| in log10 q_i quintiles and the q_i-matched lens - random |e| at
   theta < 60 arcsec.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from select_s16a_random_centres import sha

THETA = (0, 30, 60, 90, 150, 300, 600)
RPROP = (0, 0.2, 0.4, 0.7, 1.0, 1.5, 3.0)
ZSPLIT = 0.45
FLAGS = ('pixel_interpolated_center', 'pixel_suspect_center', 'pixel_clipped_any')


def truth(v):
    return str(v).strip().lower() in ('true', 't', '1')


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centres, selection_sha = load_new_centres()
    new, counts, provenance, _ = load_new(centres)
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
    lens_ids = {g: lens_row[g]['center_id'] for g in groups}
    entries = {**{lens_ids[g]: pilot[lens_ids[g]] for g in groups}, **{r: new[r] for g in groups for r in new_rand[g]}}
    data = {}
    for cid, e in entries.items():
        c = e['center']
        q, zl = (pq if cid.startswith('LC-') else nq), float(c['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        imag = g('i_mag', q)-g('a_i', q)
        keep = ~e['missing'] & (g('photoz_best', q)-zl > .2) & (g('ephor_ab_photoz_err95_min', pz)-zl > 0) & (imag > 23)
        ids = np.array(e['ids'])[keep]
        P = [pair[(cid, o)] for o in ids]
        e1, e2, w = (np.array([num(p[k]) for p in P]) for k in ('e1', 'e2', 'weight'))
        theta = np.hypot(*e['xy'][keep].T)
        data[cid] = dict(w=w, abs_e=np.hypot(e1, e2), theta=theta, rprop=theta*3.0/float(e['radius']),
                         gi=np.array([colour(q[o]) for o in ids]), imag=imag[keep], raw_sigma=g('raw_sigma', q)[keep],
                         blend=g('blendedness', q)[keep], gerr=g('g_flux_err', q)[keep], ierr=g('i_flux_err', q)[keep],
                         flag=np.array([any(truth(q[o][f]) for f in FLAGS) for o in ids], float), zl=np.full(keep.sum(), zl))
    rand_cids = [r for g in groups for r in new_rand[g]]

    def cat(cids):
        return {k: np.concatenate([data[c][k] for c in cids]) for k in data[cids[0]]}
    R = cat(rand_cids)
    t1 = float(np.nanpercentile(R['gi'], 100/3))
    rs_edges = np.nanpercentile(R['raw_sigma'], np.linspace(0, 100, 11))[1:-1]
    cell = lambda A: np.floor(A['imag']/0.25).astype(int)*100 + np.digitize(A['raw_sigma'], rs_edges)
    med = {}
    for band in ('g', 'i'):
        cr = cell(R)
        for k in np.unique(cr):
            x = R[band+'err'][(cr == k) & (R[band+'err'] > 0)]
            if len(x) >= 20:
                med[(band, k)] = float(np.median(x))
    for d in data.values():
        cd = cell(d)
        for band in ('g', 'i'):
            m = np.array([med.get((band, k), np.nan) for k in cd])
            d['lq'+band] = np.where(d[band+'err'] > 0, np.log10(d[band+'err']/m), np.nan)
    R = cat(rand_cids)
    sel = {'all': lambda A: np.ones(len(A['w']), bool), 'blue': lambda A: A['gi'] <= t1}
    zsets = {'all z': groups, f'z < {ZSPLIT}': [g for g in groups if float(lens_row[g]['z_lens']) < ZSPLIT],
             f'z >= {ZSPLIT}': [g for g in groups if float(lens_row[g]['z_lens']) >= ZSPLIT]}

    def wm(A, key, s):
        ok = s & np.isfinite(A[key])
        return float(np.average(A[key][ok], weights=A['w'][ok])) if A['w'][ok].sum() > 0 else np.nan

    def diff(gs, key, var, lo, hi, cname):
        L, Rr = cat([lens_ids[g] for g in gs]), cat([r for g in gs for r in new_rand[g]])
        sl = sel[cname](L) & (L[var] >= lo) & (L[var] < hi)
        sr = sel[cname](Rr) & (Rr[var] >= lo) & (Rr[var] < hi)
        return wm(L, key, sl), wm(Rr, key, sr), int(sl.sum())

    def jk(gs, *a):
        reps = []
        for d in gs:
            l, r, _ = diff([g for g in gs if g != d], *a)
            reps.append(l-r)
        reps = np.array(reps); reps = reps[np.isfinite(reps)]
        return float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum())) if len(reps) > 1 else np.nan
    out = []
    for var, edges in (('theta', THETA), ('rprop', RPROP)):
        for zname, gs in zsets.items():
            for cname in sel:
                for lo, hi in zip(edges[:-1], edges[1:]):
                    for key in (('abs_e', 'lqi', 'lqg', 'flag', 'blend') if zname == 'all z' else ('abs_e',)):
                        l, r, n = diff(gs, key, var, lo, hi, cname)
                        out.append(dict(binning=var, lo=lo, hi=hi, z_groups=zname, colour=cname, quantity=key, n_lens=n,
                                        lens=l, random=r, lens_minus_random=l-r, group_jackknife_se=jk(gs, key, var, lo, hi, cname)))
    qe = np.nanpercentile(R['lqi'], [20, 40, 60, 80])
    curve = []
    for k in range(5):
        s = np.isfinite(R['lqi']) & (np.digitize(R['lqi'], qe) == k)
        curve.append(dict(quintile=k, mean_log10_q_i=wm(R, 'lqi', s), abs_e=wm(R, 'abs_e', s), n=int(s.sum())))

    def matched(gs, cname):
        L, Rr = cat([lens_ids[g] for g in gs]), cat([r for g in gs for r in new_rand[g]])
        sl = sel[cname](L) & (L['theta'] < 60) & np.isfinite(L['lqi'])
        sr = sel[cname](Rr) & np.isfinite(Rr['lqi'])
        cl, cr = np.digitize(L['lqi'][sl], qe), np.digitize(Rr['lqi'][sr], qe)
        pred = wt = 0.0
        for k in range(5):
            if (cr == k).any() and (cl == k).any():
                ww = L['w'][sl][cl == k].sum()
                pred, wt = pred+ww*np.average(Rr['abs_e'][sr][cr == k], weights=Rr['w'][sr][cr == k]), wt+ww
        return wm(L, 'abs_e', sl), pred/wt
    match = {}
    for cname in sel:
        o, p_ = matched(groups, cname)
        reps = np.array([np.subtract(*matched([g for g in groups if g != d], cname)) for d in groups])
        match[cname] = dict(lens_theta_lt_60=o, random_q_i_matched=p_, difference=o-p_,
                            group_jackknife_se=float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum())))
    output.mkdir(parents=True)
    for name, tab in (('profiles.tsv', out), ('noise_curve.tsv', curve)):
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-bcg-light-abs-e-v1', science_measurement=False, selection_change=False, new_query=False,
                   theta_edges_arcsec=THETA, rprop_edges_mpc=RPROP, z_split=ZSPLIT, blue_gi_max=t1, noise_quintile_edges=list(qe),
                   noise_matched=match, z_groups={k: v for k, v in zsets.items()}, pz_receipt_summary_sha256=pz_sha,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, groups=groups,
                   limitations=['No sky or BCG-light column: flux_err in (i, raw_sigma) cells is a local-noise proxy only.',
                                'Lens z spans 0.345-0.542, so arcsec per Mpc varies by only ~30%; angle and radius are weakly separable.',
                                'Centre assumed to coincide with the BCG.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('profiles.tsv', 'noise_curve.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out, curve, match


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_bcg_light_abs_e_20261126_01'))
    out, curve, match = run(parser.parse_args().output)
    for r in out:
        print(f"{r['binning']:5s} {r['lo']:>5}-{r['hi']:<5} {r['z_groups']:9s} {r['colour']:4s} {r['quantity']:6s} n={r['n_lens']:5d} "
              f"L={r['lens']:.4f} R={r['random']:.4f} d={r['lens_minus_random']:+.4f}±{r['group_jackknife_se']:.4f}")
    for r in curve:
        print(r)
    print(json.dumps(match, indent=1, default=float))
