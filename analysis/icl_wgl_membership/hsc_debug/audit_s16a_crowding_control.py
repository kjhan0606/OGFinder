#!/usr/bin/env python3
"""Random-field crowding control for the S16A primary-source shape changes near the lens centres; local only.

Neighbours: every shape-table pair row (before the Y1 cuts) of the same centre. For each primary source
(audit_s16a_primary_boost definition): n10 = other shape-table sources within 10 arcsec, the nearest neighbour,
and the weighted centroid direction of the neighbours within 10 arcsec. Fixed before the run:
A. Random sources by n10 class (0, 1, 2, >=3): weighted |e|, e_t/e_x about the nearest neighbour (nn <= 5 arcsec)
   and about the local neighbour centroid (n10 >= 2). Lens sources in 0-0.4, 0.4-1, 1-3 Mpc, same metrics.
B. Crowding match: random |e| reweighted to the lens 0-0.4 Mpc n10-class distribution, compared with the lens value;
   the same for e_t about the nearest neighbour.
Shape-noise SE sqrt(sum w^2 x^2)/sum w; for |e| differences a group jackknife over lens groups is added.
Convention: ogfkit.wlcluster.tangential, x east, y north, e2 as given (confirmed in audit_s16a_loss_shear_bias).
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_afterburner_mask_effects import tangent_arcsec
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_outer_excess_split import new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wlcluster import tangential  # noqa: E402

R10, NN_MAX = 10.0, 5.0
CLASSES = ((0, 0), (1, 1), (2, 2), (3, 10**9))
RANGES = {'0-0.4': (0, 1, 2), '0.4-1': (3, 4), '1-3': (5, 6, 7)}


def mean_se(x, w):
    ok = np.isfinite(x)
    x, w = x[ok], w[ok]
    sw = w.sum()
    return (float((w*x).sum()/sw), float(np.sqrt(((w*x)**2).sum())/sw), int(ok.sum())) if sw > 0 else (None, None, 0)


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
    allrows = defaultdict(list)
    for p in ppairs+npairs:
        allrows[p['center_id']].append(p)
    pair = {(p['center_id'], p['object_id']): p for p in ppairs+npairs}
    pq, nq = pilot_quality(), new_quality()
    pz, pz_sha = load_pz()
    lens_ids = {g: lens_row[g]['center_id'] for g in groups}
    entries = {**{lens_ids[g]: pilot[lens_ids[g]] for g in groups}, **{r: new[r] for g in groups for r in new_rand[g]}}
    data = {}
    for cid, e in entries.items():
        c = e['center']
        nb = allrows[cid]
        nx, ny = tangent_arcsec([float(p['source_ra']) for p in nb], [float(p['source_dec']) for p in nb], c['ra_deg'], c['dec_deg'])
        nxy = np.c_[nx, ny]
        tree = cKDTree(nxy)
        q, zl = (pq if cid.startswith('LC-') else nq), float(c['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        best, lo95, imag = g('photoz_best', q)-zl, g('ephor_ab_photoz_err95_min', pz)-zl, g('i_mag', q)-g('a_i', q)
        keep = (best > .2) & (lo95 > 0) & (imag > 23) & ~e['missing']
        xy = e['xy'][keep]
        P = [pair[(cid, o)] for o in np.array(e['ids'])[keep]]
        e1, e2, w = (np.array([num(p[k]) for p in P]) for k in ('e1', 'e2', 'weight'))
        lists = tree.query_ball_point(xy, R10)
        n10, cx, cy = np.zeros(len(xy), int), np.full(len(xy), np.nan), np.full(len(xy), np.nan)
        for i, L in enumerate(lists):
            others = [j for j in L if np.hypot(*(nxy[j]-xy[i])) > 0.01]
            n10[i] = len(others)
            if others:
                cx[i], cy[i] = nxy[others].mean(0)
        dist, idx = tree.query(xy, k=2)
        self0 = dist[:, 0] <= 0.01
        nnd = np.where(self0, dist[:, 1], dist[:, 0])
        nni = np.where(self0, idx[:, 1], idx[:, 0])
        et_nn, ex_nn = tangential(xy[:, 0], xy[:, 1], e1, e2, nxy[nni, 0], nxy[nni, 1])
        et_c, ex_c = tangential(xy[:, 0], xy[:, 1], e1, e2, np.nan_to_num(cx), np.nan_to_num(cy))
        data[cid] = dict(bins=e['bins'][keep], w=w, abs_e=np.hypot(e1, e2), n10=n10, nnd=nnd,
                         et_nn=np.where(nnd <= NN_MAX, et_nn, np.nan), ex_nn=np.where(nnd <= NN_MAX, ex_nn, np.nan),
                         et_c=np.where(n10 >= 2, et_c, np.nan), ex_c=np.where(n10 >= 2, ex_c, np.nan))

    def gather(cids, ks=None):
        out = defaultdict(list)
        for c in cids:
            s = np.ones(len(data[c]['w']), bool) if ks is None else np.isin(data[c]['bins'], ks)
            for k, v in data[c].items():
                out[k].append(v[s])
        return {k: np.concatenate(v) for k, v in out.items()}
    rand_cids = [r for g in groups for r in new_rand[g]]
    rows_out = []

    def metrics(label, A):
        row = dict(sample=label, n=len(A['w']), mean_n10=float(np.average(A['n10'], weights=A['w'])),
                   frac_nn_le_3=float(np.average(A['nnd'] <= 3, weights=A['w'])))
        for k in ('abs_e', 'et_nn', 'ex_nn', 'et_c', 'ex_c'):
            m, s, n = mean_se(A[k], A['w'])
            row.update({k: m, k+'_se': s, k+'_n': n})
        rows_out.append(row)
    R = gather(rand_cids)
    for lo, hi in CLASSES:
        s = (R['n10'] >= lo) & (R['n10'] <= hi)
        metrics(f'random n10 {lo}-{hi if hi < 10**9 else "inf"}', {k: v[s] for k, v in R.items()})
    metrics('random all', R)
    for rname, ks in RANGES.items():
        metrics(f'lens {rname} Mpc', gather(list(lens_ids.values()), ks))
        metrics(f'random {rname} Mpc', gather(rand_cids, ks))

    def matched(gs, key):
        L, Rr = gather([lens_ids[g] for g in gs], RANGES['0-0.4']), gather([r for g in gs for r in new_rand[g]])
        wl, pred = 0.0, 0.0
        for lo, hi in CLASSES:
            sl, sr = (L['n10'] >= lo) & (L['n10'] <= hi), (Rr['n10'] >= lo) & (Rr['n10'] <= hi)
            m = mean_se(Rr[key][sr], Rr['w'][sr])[0]
            if m is not None and L['w'][sl].sum() > 0:
                okl = np.isfinite(L[key][sl])
                ww = L['w'][sl][okl].sum()
                pred, wl = pred+m*ww, wl+ww
        return mean_se(L[key], L['w'])[0], pred/wl
    match = {}
    for key in ('abs_e', 'et_nn', 'et_c'):
        obs, pred = matched(groups, key)
        reps = [np.subtract(*matched([g for g in groups if g != d], key)) for d in groups]
        reps = np.array(reps)
        jk = float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum()))
        match[key] = dict(lens_0_0p4=obs, random_crowding_matched=pred, difference=obs-pred, group_jackknife_se=jk)
    output.mkdir(parents=True)
    with (output/'crowding_metrics.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows_out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows_out)
    summary = dict(schema='astrafek-s16a-crowding-control-v1', science_measurement=False, selection_change=False, new_query=False,
                   radius_arcsec=R10, nn_max_arcsec=NN_MAX, classes=[list(c) for c in CLASSES], crowding_matched=match,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   groups=groups, counts=counts,
                   limitations=['Neighbours are shape-table sources only; bright members and stars outside the shape table are not counted.',
                                'Raw weighted ellipticities; shape-noise SE for single means.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'crowding_metrics.tsv': sha(output/'crowding_metrics.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return rows_out, match


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_crowding_control_20261126_01'))
    rows_out, match = run(parser.parse_args().output)
    for r in rows_out:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items() if not k.endswith('_n')})
    print(json.dumps(match, indent=1, default=float))
