#!/usr/bin/env python3
"""Do the crowding-lost sources bias the shear of the remaining S16A primary sources? Local only.

Primary sources (audit_s16a_primary_boost): standard selection + photoz_best - z_l > 0.2 + err95_min > z_l + i > 23,
11 lenses and the 76 new randoms. No areas are needed. Fixed before the run, in radial ranges 0-0.4, 0.4-1, 1-3 Mpc:
A. Weighted means (lens, random, lens - random with group-jackknife SE over lens groups) of: resolution, raw_sigma,
   blendedness, dereddened i, weight, sigma_e, rms_e, m, the responsivity 1 - rms_e^2 and |e|.
   Two-sample KS of lens vs random (0-0.4 Mpc) for resolution, blendedness, i.
B. Weighted tangential and cross ellipticity (ogfkit.wlcluster.tangential on the gnomonic x = east, y = north
   plane) and the same with e2 sign flipped; the convention with positive lens tangential signal is reported as
   primary. Random-centre cross and tangential means are nulls. Shape-noise SE: sqrt(sum w^2 e^2)/sum w.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import ks_2samp

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_outer_excess_split import new_quality, pilot_quality
from audit_s16a_pz_dilution import load_pz, num
from audit_s16a_random_extension_boosts import PAIR_FIELDS, SOURCE_SHA, load_new, load_new_centres, receipt, rows
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wlcluster import tangential  # noqa: E402

RANGES = {'0-0.4': (0, 1, 2), '0.4-1': (3, 4), '1-3': (5, 6, 7)}
PROPS = ('resolution', 'raw_sigma', 'blendedness', 'imag', 'weight', 'sigma_e', 'rms_e', 'm', 'responsivity', 'abs_e')


def wmean(x, w):
    ok = np.isfinite(x) & (w > 0)
    return float(np.sum(w[ok]*x[ok])/np.sum(w[ok])) if ok.any() else np.nan


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
        P = [pair[(cid, o)] for o in e['ids']]
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        f = lambda key: np.array([num(p[key]) for p in P])
        best, lo95, imag = g('photoz_best', q)-zl, g('ephor_ab_photoz_err95_min', pz)-zl, g('i_mag', q)-g('a_i', q)
        keep = (best > .2) & (lo95 > 0) & (imag > 23) & ~e['missing']
        e1, e2, rms = f('e1'), f('e2'), f('rms_e')
        d = dict(bins=e['bins'][keep], w=f('weight')[keep], e1=e1[keep], e2=e2[keep], x=e['xy'][keep, 0], y=e['xy'][keep, 1],
                 resolution=g('resolution', q)[keep], raw_sigma=g('raw_sigma', q)[keep], blendedness=g('blendedness', q)[keep],
                 imag=imag[keep], weight=f('weight')[keep], sigma_e=f('sigma_e')[keep], rms_e=rms[keep], m=f('m')[keep],
                 responsivity=(1-rms**2)[keep], abs_e=np.hypot(e1, e2)[keep])
        data[cid] = d

    def gather(cids, ks):
        out = defaultdict(list)
        for c in cids:
            s = np.isin(data[c]['bins'], ks)
            for k, v in data[c].items():
                out[k].append(v[s])
        return {k: np.concatenate(v) for k, v in out.items()}

    def cids(role, gs):
        return [lens_ids[g] for g in gs] if role == 'lens' else [r for g in gs for r in new_rand[g]]
    prop_rows = []
    for rname, ks in RANGES.items():
        for p in PROPS:
            vals = {role: wmean(gather(cids(role, groups), ks)[p], gather(cids(role, groups), ks)['w']) for role in ('lens', 'random')}
            reps = []
            for dg in groups:
                gs = [g for g in groups if g != dg]
                L, R = gather(cids('lens', gs), ks), gather(cids('random', gs), ks)
                reps.append(wmean(L[p], L['w'])-wmean(R[p], R['w']))
            reps = np.array(reps)
            se = float(np.sqrt((len(reps)-1)/len(reps)*np.nansum((reps-np.nanmean(reps))**2)))
            prop_rows.append(dict(range=rname, property=p, lens=vals['lens'], random=vals['random'], lens_minus_random=vals['lens']-vals['random'],
                                  group_jackknife_se=se, n_lens=len(gather(cids('lens', groups), ks)['w']),
                                  n_random=len(gather(cids('random', groups), ks)['w'])))
    L, R = gather(cids('lens', groups), RANGES['0-0.4']), gather(cids('random', groups), RANGES['0-0.4'])
    ks_rows = [dict(property=p, ks_statistic=float(ks_2samp(L[p][np.isfinite(L[p])], R[p][np.isfinite(R[p])]).statistic),
                    p_value=float(ks_2samp(L[p][np.isfinite(L[p])], R[p][np.isfinite(R[p])]).pvalue)) for p in ('resolution', 'blendedness', 'imag')]
    shear_rows = []
    for conv, sgn in (('e2 as given', 1.0), ('e2 sign flipped', -1.0)):
        for role in ('lens', 'random'):
            for rname, ks in RANGES.items():
                A = gather(cids(role, groups), ks)
                gt, gx = tangential(A['x'], A['y'], A['e1'], sgn*A['e2'])
                sw = np.sum(A['w'])
                shear_rows.append(dict(convention=conv, role=role, range=rname, n=len(A['w']),
                                       e_t=float(np.sum(A['w']*gt)/sw), e_t_se=float(np.sqrt(np.sum((A['w']*gt)**2))/sw),
                                       e_x=float(np.sum(A['w']*gx)/sw), e_x_se=float(np.sqrt(np.sum((A['w']*gx)**2))/sw)))
    lens_signal = {r['convention']: r['e_t']/r['e_t_se'] for r in shear_rows if r['role'] == 'lens' and r['range'] == '0.4-1'}
    primary = max(lens_signal, key=lens_signal.get)
    output.mkdir(parents=True)
    tables = {'properties.tsv': prop_rows, 'ks.tsv': ks_rows, 'shear_nulls.tsv': shear_rows}
    for name, tab in tables.items():
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-loss-shear-bias-v1', science_measurement=False, selection_change=False, new_query=False,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   convention_with_positive_lens_signal=primary, lens_signal_to_noise_0p4_1=lens_signal, groups=groups, counts=counts,
                   limitations=['Raw weighted ellipticities: no responsivity or m correction applied to e_t/e_x.',
                                'Shape-noise SE only for the shear means (no cosmic variance or source-position covariance).'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in tables})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return tables, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_loss_shear_bias_20261126_01'))
    tables, s = run(parser.parse_args().output)
    for t, tab in tables.items():
        for r in tab:
            print(t[:5], {k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()})
    print(json.dumps({k: s[k] for k in ('convention_with_positive_lens_signal', 'lens_signal_to_noise_0p4_1')}))
