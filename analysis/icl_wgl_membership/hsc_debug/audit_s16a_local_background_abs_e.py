#!/usr/bin/env python3
"""Candidate (a) with the received local-background columns (receipt 6851a9e9): does BCG/ICL or neighbour light at the
source position explain the lower |e| of S16A primary sources near the lens centres? Local only.

Primary sources as audit_s16a_primary_boost; blue = g-i <= random tertile (as audit_s16a_bcg_light_abs_e). Fixed before the run:
A. Lens - random (group jackknife) in angular bins, all and blue, of:
   lv_b  = log10 variance_b minus the random median in the same dereddened i (0.25 mag) x raw_sigma decile cell (b = g, i);
   icountinputs; extra = flux_parent/flux_child - 1 (i, child-matched Gaussian weight);
   e_extra = (dQ11 - dQ22, 2 dQ12)/(child Q11 + Q22) with dQ = parent - child moments: |e_extra| and its tangential and
   cross components about the centre (moment axes assumed as the shape e1/e2 axes; convention not independently checked);
   child fraction (parent_id != 0) and any i deblend flag.
B. Random |e| in quintiles of extra and of |e_extra|, and the lens - random |e| at theta < 60 arcsec matched to the lens
   distribution of each proxy (quintiles of random primary sources).
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
from audit_s16a_outer_excess_split import colour, new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from prepare_s16a_local_background_rows import OUTPUT as LB_PLAN
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from s16a_y1_source_contract import boolean
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wlcluster import tangential  # noqa: E402

LB_SHA = '6851a9e9c795895a5b37c8c29523403831972a122d4aadfa9c7955b0660150b2'
THETA = (0, 30, 60, 90, 150, 300, 600)
FLAGS = ('i_deblend_has_stray_flux', 'i_deblend_too_many_peaks', 'i_deblend_parent_too_big', 'i_deblend_masked')
QUANTITIES = ('abs_e', 'lvi', 'lvg', 'icount', 'extra', 'abs_e_extra', 'et_extra', 'ex_extra', 'child', 'flag')


def load_lb():
    res, summary = receipt(LB_PLAN, LB_SHA)
    return {r['object_id']: r for r in rows(res, sorted(summary['output_sha256']))}, sha(res/'download_summary.json')


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
    lb, lb_sha = load_lb()
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
        L = lambda key: np.array([num(lb[o][key]) for o in ids])
        xy = e['xy'][keep]
        c11, c22, c12 = (L('i_blendedness_abs_shape_child_'+k) for k in ('11', '22', '12'))
        d11, d22, d12 = (L('i_blendedness_abs_shape_parent_'+k)-v for k, v in (('11', c11), ('22', c22), ('12', c12)))
        tr = c11+c22
        ee1, ee2 = np.where(tr > 0, (d11-d22)/tr, np.nan), np.where(tr > 0, 2*d12/tr, np.nan)
        et, ex = tangential(xy[:, 0], xy[:, 1], np.nan_to_num(ee1), np.nan_to_num(ee2))
        fc, fp = L('i_blendedness_abs_flux_child'), L('i_blendedness_abs_flux_parent')
        pid = [lb[o]['parent_id'].strip() for o in ids]
        data[cid] = dict(w=w, abs_e=np.hypot(e1, e2), theta=np.hypot(*xy.T), gi=np.array([colour(q[o]) for o in ids]), imag=imag[keep],
                         raw_sigma=g('raw_sigma', q)[keep], vi=L('i_variance'), vg=L('g_variance'), icount=L('i_countinputs'),
                         extra=np.where(fc > 0, fp/fc-1, np.nan), abs_e_extra=np.hypot(ee1, ee2),
                         et_extra=np.where(np.isfinite(ee1), et, np.nan), ex_extra=np.where(np.isfinite(ee1), ex, np.nan),
                         child=np.array([np.nan if p == '' else float(p != '0') for p in pid]),
                         flag=np.array([float(any(boolean(lb[o][f].strip()) is True for f in FLAGS)) for o in ids]))
    rand_cids = [r for g in groups for r in new_rand[g]]

    def cat(cids):
        return {k: np.concatenate([data[c][k] for c in cids]) for k in data[cids[0]]}
    R = cat(rand_cids)
    t1 = float(np.nanpercentile(R['gi'], 100/3))
    rs_edges = np.nanpercentile(R['raw_sigma'], np.linspace(0, 100, 11))[1:-1]
    cell = lambda A: np.floor(A['imag']/0.25).astype(int)*100 + np.digitize(A['raw_sigma'], rs_edges)
    med = {}
    cr = cell(R)
    for band in ('g', 'i'):
        for k in np.unique(cr):
            x = R['v'+band][(cr == k) & (R['v'+band] > 0)]
            if len(x) >= 20:
                med[(band, k)] = float(np.median(x))
    for d in data.values():
        cd = cell(d)
        for band in ('g', 'i'):
            m = np.array([med.get((band, k), np.nan) for k in cd])
            d['lv'+band] = np.where(d['v'+band] > 0, np.log10(d['v'+band]/m), np.nan)
    R = cat(rand_cids)
    sel = {'all': lambda A: np.ones(len(A['w']), bool), 'blue': lambda A: A['gi'] <= t1}

    def wm(A, key, s):
        ok = s & np.isfinite(A[key])
        return float(np.average(A[key][ok], weights=A['w'][ok])) if A['w'][ok].sum() > 0 else np.nan

    def diff(gs, key, lo, hi, cname):
        L, Rr = cat([lens_ids[g] for g in gs]), cat([r for g in gs for r in new_rand[g]])
        sl = sel[cname](L) & (L['theta'] >= lo) & (L['theta'] < hi)
        sr = sel[cname](Rr) & (Rr['theta'] >= lo) & (Rr['theta'] < hi)
        return wm(L, key, sl), wm(Rr, key, sr), int(sl.sum())

    def jk(fn):
        reps = np.array([fn([g for g in groups if g != d]) for d in groups])
        reps = reps[np.isfinite(reps)]
        return float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum())) if len(reps) > 1 else np.nan
    out = []
    for cname in sel:
        for lo, hi in zip(THETA[:-1], THETA[1:]):
            for key in QUANTITIES:
                l, r, n = diff(groups, key, lo, hi, cname)
                out.append(dict(colour=cname, lo=lo, hi=hi, quantity=key, n_lens=n, lens=l, random=r, lens_minus_random=l-r,
                                group_jackknife_se=jk(lambda gs: np.subtract(*diff(gs, key, lo, hi, cname)[:2]))))
    curves, match = [], {}
    for proxy in ('extra', 'abs_e_extra', 'lvi'):
        qe = np.nanpercentile(R[proxy], [20, 40, 60, 80])
        for k in range(5):
            s = np.isfinite(R[proxy]) & (np.digitize(R[proxy], qe) == k)
            curves.append(dict(proxy=proxy, quintile=k, mean_proxy=wm(R, proxy, s), abs_e=wm(R, 'abs_e', s), n=int(s.sum())))

        def matched(gs, cname, proxy=proxy, qe=qe):
            L, Rr = cat([lens_ids[g] for g in gs]), cat([r for g in gs for r in new_rand[g]])
            sl = sel[cname](L) & (L['theta'] < 60) & np.isfinite(L[proxy])
            sr = sel[cname](Rr) & np.isfinite(Rr[proxy])
            cl, crr = np.digitize(L[proxy][sl], qe), np.digitize(Rr[proxy][sr], qe)
            pred = wt = 0.0
            for k in range(5):
                if (crr == k).any() and (cl == k).any():
                    ww = L['w'][sl][cl == k].sum()
                    pred, wt = pred+ww*np.average(Rr['abs_e'][sr][crr == k], weights=Rr['w'][sr][crr == k]), wt+ww
            return wm(L, 'abs_e', sl), pred/wt
        for cname in sel:
            o, p_ = matched(groups, cname)
            match[f'{proxy} {cname}'] = dict(lens_theta_lt_60=o, random_matched=p_, difference=o-p_,
                                             group_jackknife_se=jk(lambda gs: np.subtract(*matched(gs, cname))))
    output.mkdir(parents=True)
    for name, tab in (('profiles.tsv', out), ('proxy_curves.tsv', curves)):
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-local-background-abs-e-v1', science_measurement=False, selection_change=False, new_query=False,
                   theta_edges_arcsec=THETA, blue_gi_max=t1, proxy_matched=match, local_background_receipt_summary_sha256=lb_sha,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   groups=groups,
                   limitations=['Variance includes source photon noise; matched in i x raw_sigma cells only.',
                                'Moment axes assumed equal to the shape e1/e2 axes for the tangential split; |e_extra| is axis-free.',
                                'Centre assumed to coincide with the BCG.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('profiles.tsv', 'proxy_curves.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out, curves, match


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_local_background_abs_e_20261127_01'))
    out, curves, match = run(parser.parse_args().output)
    for r in out:
        print(f"{r['colour']:4s} {r['lo']:>4}-{r['hi']:<4} {r['quantity']:11s} n={r['n_lens']:5d} L={r['lens']:.4f} R={r['random']:.4f} "
              f"d={r['lens_minus_random']:+.4f}±{r['group_jackknife_se']:.4f}")
    for r in curves:
        print(r)
    print(json.dumps(match, indent=1, default=float))
