#!/usr/bin/env python3
"""Does the higher resolution of the S16A primary sources near the lenses explain their lower |e|? Local only.

Primary sources as audit_s16a_primary_boost. Fixed before the run:
A. Random primary sources: weighted |e| in resolution deciles (decile edges from the random primary set) and the
   weighted linear slope d|e|/d resolution.
B. Matched prediction: random |e| reweighted to the lens resolution-decile distribution per radial range, for all
   primary sources and per g-i tertile (deciles within the tertile); lens - matched with group jackknife.
   Also matched jointly in resolution decile x dereddened i bin (<=24, >24).
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

RANGES = {'0-0.4': (0, 1, 2), '0.4-1': (3, 4), '1-3': (5, 6, 7)}


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
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        P = [pair[(cid, o)] for o in e['ids']]
        e1, e2, w = (np.array([num(p[k]) for p in P]) for k in ('e1', 'e2', 'weight'))
        imag = g('i_mag', q)-g('a_i', q)
        keep = ~e['missing'] & (g('photoz_best', q)-zl > .2) & (g('ephor_ab_photoz_err95_min', pz)-zl > 0) & (imag > 23)
        gi = np.array([colour(q[o]) for o in e['ids']])
        data[cid] = dict(bins=e['bins'][keep], w=w[keep], abs_e=np.hypot(e1, e2)[keep], res=g('resolution', q)[keep],
                         imag=imag[keep], gi=gi[keep], all_gi=gi[~e['missing']])
    rand_cids = [r for g in groups for r in new_rand[g]]
    t1, t2 = np.nanpercentile(np.concatenate([data[c]['all_gi'] for c in rand_cids]), [100/3, 200/3])

    def cat(cids, ks=None):
        out = defaultdict(list)
        for c in cids:
            s = np.ones(len(data[c]['w']), bool) if ks is None else np.isin(data[c]['bins'], ks)
            for k in ('w', 'abs_e', 'res', 'imag', 'gi'):
                out[k].append(data[c][k][s])
        return {k: np.concatenate(v) for k, v in out.items()}
    R = cat(rand_cids)
    edges = np.nanpercentile(R['res'], np.linspace(0, 100, 11))
    edges[0], edges[-1] = -np.inf, np.inf
    dec = lambda x: np.clip(np.digitize(x, edges[1:-1]), 0, 9)
    curve = []
    for k in range(10):
        s = dec(R['res']) == k
        curve.append(dict(decile=k, res_lo=float(edges[k]), res_hi=float(edges[k+1]), n=int(s.sum()),
                          mean_res=float(np.average(R['res'][s], weights=R['w'][s])), abs_e=float(np.average(R['abs_e'][s], weights=R['w'][s]))))
    ok = np.isfinite(R['res']) & np.isfinite(R['abs_e'])
    slope = float(np.polyfit(R['res'][ok], R['abs_e'][ok], 1, w=np.sqrt(R['w'][ok]))[0])
    tert = {'all': lambda A: np.ones(len(A['w']), bool), 'blue': lambda A: A['gi'] <= t1,
            'mid': lambda A: (A['gi'] > t1) & (A['gi'] <= t2), 'red': lambda A: A['gi'] > t2}

    def matched(gs, ks, tname, joint):
        L = cat([lens_ids[g] for g in gs], ks)
        Rr = cat([r for g in gs for r in new_rand[g]])
        sl, sr = tert[tname](L), tert[tname](Rr)
        cell = lambda A: dec(A['res'])*2 + ((A['imag'] > 24).astype(int) if joint else 0)
        cl, cr = cell(L)[sl], cell(Rr)[sr]
        wl, el, wr, er = L['w'][sl], L['abs_e'][sl], Rr['w'][sr], Rr['abs_e'][sr]
        pred = num_ = 0.0
        for c in np.unique(cl):
            m = cr == c
            if m.any() and wr[m].sum() > 0:
                ww = wl[cl == c].sum()
                pred, num_ = pred+ww*np.average(er[m], weights=wr[m]), num_+ww
        obs = float(np.average(el, weights=wl)) if wl.sum() > 0 else np.nan
        return obs, (pred/num_ if num_ else np.nan), float(np.average(er, weights=wr)), int(sl.sum()), \
            float(np.average(L['res'][sl], weights=wl)) - float(np.average(Rr['res'][sr], weights=wr))
    out = []
    for rname, ks in RANGES.items():
        for tname in tert:
            for joint in (False, True):
                obs, pred, unm, n, dres = matched(groups, ks, tname, joint)
                reps = []
                for d in groups:
                    o, p_, _, _, _ = matched([g for g in groups if g != d], ks, tname, joint)
                    reps.append(o-p_)
                reps = np.array(reps); reps = reps[np.isfinite(reps)]
                jk = float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum()))
                out.append(dict(range=rname, colour=tname, match='resolution x i' if joint else 'resolution', n_lens=n,
                                lens_abs_e=obs, random_unmatched=unm, random_matched=pred, lens_minus_matched=obs-pred,
                                lens_minus_unmatched=obs-unm, group_jackknife_se=jk, lens_minus_random_resolution=dres,
                                slope_prediction=slope*dres))
    output.mkdir(parents=True)
    for name, tab in (('resolution_curve.tsv', curve), ('matched.tsv', out)):
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-resolution-abs-e-v1', science_measurement=False, selection_change=False, new_query=False,
                   slope_abs_e_per_resolution=slope, resolution_decile_edges=[float(x) for x in edges[1:-1]], gi_tertiles=[float(t1), float(t2)],
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, groups=groups,
                   limitations=['Raw |e|; the resolution-|e| relation is measured on random fields and assumed to hold near the lenses.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('resolution_curve.tsv', 'matched.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return curve, out, slope


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_resolution_abs_e_20261126_01'))
    curve, out, slope = run(parser.parse_args().output)
    for r in curve:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})
    print('slope', round(slope, 4))
    for r in out:
        print(f"{r['range']:6s} {r['colour']:5s} {r['match']:15s} n={r['n_lens']:5d} L={r['lens_abs_e']:.3f} Rm={r['random_matched']:.3f} "
              f"L-Rm={r['lens_minus_matched']:+.3f}±{r['group_jackknife_se']:.3f} L-R={r['lens_minus_unmatched']:+.3f} dres={r['lens_minus_random_resolution']:+.3f} slope_pred={r['slope_prediction']:+.4f}")
