#!/usr/bin/env python3
"""Is the lower |e| of S16A primary sources at 0-0.4 Mpc residual contamination by round members? Local only.

Primary sources as audit_s16a_primary_boost. Fixed before the run:
A. Weighted |e| lens and random (0-0.4, 0.4-1, 1-3 Mpc), lens - random with group jackknife, for primary sources split by
   observed g-i tertile (pooled new-random tertiles of all selected sources, as audit_s16a_outer_excess_split) and by
   photoz_best - z_l in (0.2,0.5], (0.5,inf).
B. Weighted |e| of member-like sources (red tertile, photoz_best - z_l <= 0.1, any magnitude and i > 23), lens and random.
C. Contamination needed: f = (|e|_random - |e|_lens)/(|e|_random - |e|_member-like lens 0-0.4 Mpc) for the primary set.
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
        data[cid] = dict(bins=e['bins'], w=w, abs_e=np.hypot(e1, e2), ok=~e['missing'], best=g('photoz_best', q)-zl,
                         lo95=g('ephor_ab_photoz_err95_min', pz)-zl, imag=g('i_mag', q)-g('a_i', q),
                         gi=np.array([colour(q[o]) for o in e['ids']]))
    rand_cids = [r for g in groups for r in new_rand[g]]
    t1, t2 = np.nanpercentile(np.concatenate([data[c]['gi'][data[c]['ok']] for c in rand_cids]), [100/3, 200/3])
    primary = lambda d: d['ok'] & (d['best'] > .2) & (d['lo95'] > 0) & (d['imag'] > 23)
    subsets = {'primary all': primary,
               f'primary blue (g-i <= {t1:.3f})': lambda d: primary(d) & (d['gi'] <= t1),
               'primary mid': lambda d: primary(d) & (d['gi'] > t1) & (d['gi'] <= t2),
               f'primary red (g-i > {t2:.3f})': lambda d: primary(d) & (d['gi'] > t2),
               'primary dz (0.2,0.5]': lambda d: primary(d) & (d['best'] <= .5),
               'primary dz > 0.5': lambda d: primary(d) & (d['best'] > .5),
               'member-like (red, dz <= 0.1), all i': lambda d: d['ok'] & (d['gi'] > t2) & (d['best'] <= .1),
               'member-like (red, dz <= 0.1), i > 23': lambda d: d['ok'] & (d['gi'] > t2) & (d['best'] <= .1) & (d['imag'] > 23)}

    def wm(cids, f, ks):
        x, w = [], []
        for c in cids:
            d = data[c]
            s = f(d) & np.isin(d['bins'], ks)
            x.append(d['abs_e'][s]); w.append(d['w'][s])
        x, w = np.concatenate(x), np.concatenate(w)
        return (float(np.sum(w*x)/np.sum(w)) if w.sum() > 0 else np.nan), len(x), (float(np.sqrt(np.sum((w*(x-np.sum(w*x)/np.sum(w)))**2))/np.sum(w)) if w.sum() > 0 else np.nan)
    out = []
    for name, f in subsets.items():
        for rname, ks in RANGES.items():
            L, nl, sl = wm([lens_ids[g] for g in groups], f, ks)
            R, nr, sr = wm(rand_cids, f, ks)
            reps = np.array([wm([lens_ids[g] for g in groups if g != d], f, ks)[0]-wm([r for g in groups if g != d for r in new_rand[g]], f, ks)[0]
                             for d in groups])
            reps = reps[np.isfinite(reps)]
            jk = float(np.sqrt((len(reps)-1)/len(reps)*((reps-reps.mean())**2).sum())) if len(reps) > 1 else None
            out.append(dict(subset=name, range=rname, lens=L, lens_se=sl, n_lens=nl, random=R, random_se=sr, n_random=nr,
                            lens_minus_random=L-R, group_jackknife_se=jk))
    get = lambda s, r, k: next(x[k] for x in out if x['subset'] == s and x['range'] == r)
    e_r, e_l = get('primary all', '0-0.4', 'random'), get('primary all', '0-0.4', 'lens')
    e_m = get('member-like (red, dz <= 0.1), i > 23', '0-0.4', 'lens')
    need = dict(random_primary=e_r, lens_primary=e_l, member_like_lens_i_gt_23=e_m,
                member_like_random_i_gt_23=get('member-like (red, dz <= 0.1), i > 23', '0-0.4', 'random'),
                contamination_needed=(e_r-e_l)/(e_r-e_m) if e_r > e_m else None)
    output.mkdir(parents=True)
    with (output/'abs_e_splits.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)
    summary = dict(schema='astrafek-s16a-round-member-test-v1', science_measurement=False, selection_change=False, new_query=False,
                   gi_tertiles=[float(t1), float(t2)], contamination=need, pz_receipt_summary_sha256=pz_sha,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, groups=groups,
                   limitations=['Raw |e| (no responsivity); observed g-i without extinction term.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'abs_e_splits.tsv': sha(output/'abs_e_splits.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out, need


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_round_member_test_20261126_01'))
    out, need = run(parser.parse_args().output)
    f = lambda x: '—' if x is None or (isinstance(x, float) and np.isnan(x)) else f'{x:.3f}'
    for r in out:
        print(f"{r['subset'][:38]:38s} {r['range']:6s} L={f(r['lens'])}({r['n_lens']}) R={f(r['random'])}({r['n_random']}) "
              f"d={f(r['lens_minus_random'])} ± {f(r['group_jackknife_se'])}")
    print(json.dumps(need, default=float))
