#!/usr/bin/env python3
"""Where does the S16A 2-3 Mpc boost of ~1.10 (new randoms, 2026-11-26) come from? Local only.

A. Property split. The lens and new-random entries, areas (label_fdfc) and stack (expected count, group
   jackknife) are those of audit_s16a_random_extension_boosts; only the counted sources are restricted to a
   subset of the selected, photo-z-available sources. Subsets, fixed before the run:
   dereddened i (i_mag - a_i) in (-inf,22], (22,23], (23,24], (24,24.5];
   photoz_best - z_lens in (-inf,0.1], (0.1,0.4], (0.4,inf);
   observed g-i (cmodel fluxes, no colour extinction term) tertiles of the pooled new-random sources, plus
   'g-i undefined' (non-positive flux);
   'member-like': red tertile and photoz_best - z_lens <= 0.1.
   Excess share = (N_lens,S - E_S)/(N_lens - E) in the annulus, E the expected count from the randoms.
B. Selection efficiency. For every pair row in the shape table (before the Y1 cuts), the Y1 decision and each
   criterion. Pass fraction lens/random per annulus (groups pooled, group jackknife), the failure fraction of
   each criterion, and the boost of all shape-table rows over the FDFC footprint area.
No new query, no selection change, no science measurement.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, load_inputs
from audit_s16a_fdfc_areas import VALIDATION
from audit_s16a_random_extension_boosts import (ID_SHA, PAIR_FIELDS, SOURCE_SHA, label_fdfc, load_new, load_new_centres,
                                                receipt, rows)
from audit_s23b_s16a_quality_supplement import numeric_calibration
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES, radial_bins, separations_radians
from prepare_s16a_random_id_supplement import PARTS, output as id_plan
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from s16a_y1_source_contract import finite, select_row
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

MAG_EDGES = (-np.inf, 22., 23., 24., 24.5)
DZ_EDGES = (-np.inf, .1, .4, np.inf)
FOCUS = (6, 7)


def colour(q):
    g, i = finite(q['g_flux']), finite(q['i_flux'])
    return -2.5*np.log10(g/i) if g is not None and i is not None and g > 0 and i > 0 else np.nan


def props(q, z_lens):
    return dict(imag=float(finite(q['i_mag'])-finite(q['a_i'])), dz=float(finite(q['photoz_best']))-z_lens if finite(q['photoz_best']) is not None else np.nan,
                gi=colour(q))


def pilot_quality():
    out = {}
    results = PLAN/'results'
    for r in rows(results, sorted(p.name for p in results.glob('s16a_quality_supplement_*.csv'))):
        out[r['object_id']] = r
    return out


def new_quality():
    out = {}
    for part in range(1, PARTS+1):
        res, s = receipt(id_plan(part), ID_SHA[part])
        for r in rows(res, sorted(n for n in s['output_sha256'] if n.startswith('quality_'))):
            out[r['object_id']] = r
    return out


def subset_parts(p, e, mask):
    q = dict(p)
    keep = ~e['missing'] & mask
    q['n_post'] = np.bincount(e['bins'][keep], minlength=bt.NBINS).astype(float)
    return q


def stack_rows(items, label):
    out = []
    for k in range(bt.NBINS):
        value, used, se, _ = bt.stack_jackknife('expected', items, k, 'n', 'area', (), None, False)
        n_l = sum(lp['n_post'][k] for lp, _ in items)
        e = sum(bt.density(rs, 'n', 'area', (k,))*lp['a_unm'][k] for lp, rs in items if bt.density(rs, 'n', 'area', (k,)))
        out.append(dict(subset=label, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], boost=value, group_jackknife_se=se,
                        n_lens=n_l, expected=e, excess=n_l-e, lens_groups=len(used)))
    return out


def all_pair_rows(pairs, quality, centres, cosmo):
    """Per pair: centre id, annulus, Y1 decision and criteria."""
    out = defaultdict(list)
    for p in pairs:
        c = centres.get(p['center_id'])
        if c is None:
            continue
        merged = {**{k: v for k, v in quality[p['object_id']].items() if k != 'object_id'}, **p}
        sel = select_row(merged)
        ok = sel['table4_diagnostic'] == 'pass' and numeric_calibration(merged)
        out[p['center_id']].append((float(p['source_ra']), float(p['source_dec']), ok, sel['criteria'], numeric_calibration(merged)))
    res = {}
    for cid, rs in out.items():
        c = centres[cid]
        r = separations_radians(np.array([x[0] for x in rs]), np.array([x[1] for x in rs]), c['ra_deg'], c['dec_deg']) \
            * cosmo.angular_diameter_distance(c['z_lens']).value
        res[cid] = [(b, ok, crit, cal) for b, (_, _, ok, crit, cal) in zip(radial_bins(r), rs)]
    return res


def efficiency(per_centre, lens_of, randoms_of, groups):
    def tally(cids):
        t = defaultdict(Counter)
        for cid in cids:
            for b, ok, crit, cal in per_centre.get(cid, []):
                t[b]['all'] += 1
                t[b]['pass'] += ok
                t[b]['calibration_fail'] += not cal
                for name, v in crit.items():
                    t[b]['fail_'+name] += v == 'fail'
                    t[b]['unknown_'+name] += v not in ('pass', 'fail')
        return t
    lens = {g: tally([lens_of[g]]) for g in groups}
    rand = {g: tally(randoms_of[g]) for g in groups}

    def frac(gs, side, k, key):
        a = sum(side[g][k]['all'] for g in gs)
        return sum(side[g][k][key] for g in gs)/a if a else None
    keys = sorted({key for side in (lens, rand) for g in groups for k in side[g] for key in side[g][k] if key != 'all'})
    out = []
    for k in range(bt.NBINS):
        for key in keys:
            fl, fr = frac(groups, lens, k, key), frac(groups, rand, k, key)
            ratio = fl/fr if fl is not None and fr else None
            reps = []
            for d in groups:
                gs = [g for g in groups if g != d]
                a, b = frac(gs, lens, k, key), frac(gs, rand, k, key)
                reps.append(a/b if a is not None and b else None)
            se = bt.jackknife_strict([None if r is None else r-1 for r in reps])
            out.append(dict(annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], quantity=key, lens_fraction=fl, random_fraction=fr,
                            lens_over_random=ratio, group_jackknife_se=se,
                            n_lens_rows=sum(lens[g][k]['all'] for g in groups), n_random_rows=sum(rand[g][k]['all'] for g in groups)))
    return out


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fdfc = load_fdfc()
    acc = json.loads((VALIDATION/'user_acceptance.json').read_text())
    if not (acc['accepted'] and acc['grid_summary_sha256'] == fdfc['summary_sha256']):
        raise ValueError('FDFC grid not accepted')
    wide_check = merge_wide_fdfc(fdfc['grids'], load_wide())
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    centres, selection_sha = load_new_centres()
    new, counts, provenance, _ = load_new(centres)
    pilot_centres, pilot, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in pilot_centres if c['role'] == 'lens'}
    groups = sorted(lens_row)
    new_rand = defaultdict(list)
    for cid, c in centres.items():
        new_rand[c['matched_lens_id']].append(cid)
    pq, nq = pilot_quality(), new_quality()

    parts, entries = {}, {}
    for g in groups:
        lid = lens_row[g]['center_id']
        parts[lid], entries[lid] = label_fdfc(pilot[lid], lens_row[g], fdfc['grids'][g]), pilot[lid]
        for cid in new_rand[g]:
            parts[cid], entries[cid] = label_fdfc(new[cid], lens_row[g], fdfc['grids'][g]), new[cid]
    zl = {cid: float(e['center']['z_lens']) for cid, e in entries.items()}
    P = {cid: [props((pq if cid.startswith('LC-') else nq)[o], zl[cid]) for o in e['ids']] for cid, e in entries.items()}
    arr = {cid: {k: np.array([x[k] for x in v], float) for k in ('imag', 'dz', 'gi')} for cid, v in P.items()}
    pool = np.concatenate([arr[c]['gi'][~entries[c]['missing']] for g in groups for c in new_rand[g]])
    t1, t2 = np.nanpercentile(pool, [100/3, 200/3])
    subsets = [('all', lambda a: np.ones(len(a['imag']), bool))]
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        subsets.append((f'i_dered ({lo},{hi}]', lambda a, lo=lo, hi=hi: (a['imag'] > lo) & (a['imag'] <= hi)))
    for lo, hi in zip(DZ_EDGES[:-1], DZ_EDGES[1:]):
        subsets.append((f'zphot-zlens ({lo},{hi}]', lambda a, lo=lo, hi=hi: (a['dz'] > lo) & (a['dz'] <= hi)))
    subsets += [(f'g-i <= {t1:.3f}', lambda a: a['gi'] <= t1), (f'g-i ({t1:.3f},{t2:.3f}]', lambda a: (a['gi'] > t1) & (a['gi'] <= t2)),
                (f'g-i > {t2:.3f}', lambda a: a['gi'] > t2), ('g-i undefined', lambda a: np.isnan(a['gi'])),
                ('member-like: red tertile and zphot-zlens <= 0.1', lambda a: (a['gi'] > t2) & (a['dz'] <= .1))]
    split = []
    for label, f in subsets:
        sub = {cid: subset_parts(parts[cid], entries[cid], f(arr[cid])) for cid in entries}
        items = [(sub[lens_row[g]['center_id']], [sub[r] for r in new_rand[g]]) for g in groups]
        split += stack_rows(items, label)
    total = {r['annulus']: r['excess'] for r in split if r['subset'] == 'all'}
    for r in split:
        r['excess_share'] = r['excess']/total[r['annulus']] if total[r['annulus']] else None

    # B: selection efficiency over all shape-table pair rows
    _, ppairs, _, _, _, _ = load_inputs(PLAN, OBJECT_QA)
    lens_by_id = {lens_row[g]['center_id']: lens_row[g] for g in groups}
    per = all_pair_rows(ppairs, pq, lens_by_id, cosmo)
    res, summ = receipt(SOURCE_PLAN, SOURCE_SHA)
    npairs = [{k: r[k] for k in PAIR_FIELDS} for r in rows(res, sorted(n for n in summ['output_sha256'] if 'pairs' in n))]
    per.update(all_pair_rows(npairs, nq, centres, cosmo))
    eff = efficiency(per, {g: lens_row[g]['center_id'] for g in groups}, new_rand, groups)
    allrows = {}
    for cid, e in entries.items():
        q = dict(parts[cid])
        q['n_post'] = np.bincount([b for b, *_ in per[cid]], minlength=bt.NBINS).astype(float)
        q['a_unm'] = q['a_fp']
        allrows[cid] = q
    all_boost = stack_rows([(allrows[lens_row[g]['center_id']], [allrows[r] for r in new_rand[g]]) for g in groups],
                           'all shape-table rows over FDFC footprint area')

    output.mkdir(parents=True)
    for fname, tab in (('property_split.tsv', split), ('selection_efficiency.tsv', eff), ('all_rows_boost.tsv', all_boost)):
        with (output/fname).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(dict.fromkeys(k for r in tab for k in r)), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-outer-excess-split-v1', science_measurement=False, selection_change=False, new_query=False,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check,
                   gi_tertiles_from_new_randoms=[float(t1), float(t2)], mag_edges=[float(x) for x in MAG_EDGES], dz_edges=[float(x) for x in DZ_EDGES],
                   groups=groups, counts=counts,
                   limitations=['g-i is observed (no colour extinction correction); tertiles are pooled over all lens redshifts.',
                                'Pass fractions use shape-table rows (weaklensing_hsm_regauss join), not all detections.',
                                'Areas are label_fdfc of the full selected sample for every subset.'],
                   runner_sha256=sha(Path(__file__)),
                   output_sha256={n: sha(output/n) for n in ('property_split.tsv', 'selection_efficiency.tsv', 'all_rows_boost.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return split, eff, all_boost


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_outer_excess_split_20261126_01'))
    split, eff, all_boost = run(parser.parse_args().output)
    for r in split:
        if r['annulus'] in FOCUS or r['annulus'] in (0, 1, 2):
            print(f"{r['subset']:48s} {r['r_lo']}-{r['r_hi']:<4} B={r['boost'] if r['boost'] is None else round(r['boost'],3)} "
                  f"± {r['group_jackknife_se'] if r['group_jackknife_se'] is None else round(r['group_jackknife_se'],3)} "
                  f"N={r['n_lens']:.0f} share={r['excess_share'] if r['excess_share'] is None else round(r['excess_share'],2)}")
    for r in eff:
        if r['annulus'] in FOCUS and r['quantity'] in ('pass',) or (r['annulus'] == 7 and r['lens_over_random'] and abs(r['lens_over_random']-1) > .05 and (r['lens_fraction'] or 0) > .01):
            print(f"eff {r['annulus']} {r['quantity']:32s} L={r['lens_fraction']:.4f} R={r['random_fraction']:.4f} ratio={r['lens_over_random']:.3f} ± {r['group_jackknife_se']}")
    for r in all_boost:
        print(f"allrows {r['r_lo']}-{r['r_hi']} B={r['boost']:.3f} ± {r['group_jackknife_se']}")
