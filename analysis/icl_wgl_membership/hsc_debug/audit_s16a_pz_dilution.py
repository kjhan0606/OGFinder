#!/usr/bin/env python3
"""P(z)-summary dilution check of the S16A standard boost; local only.

Inputs: the standard-boost setup (audit_s16a_standard_boost: 76 new randoms, random-point + FDFC areas in
0-0.4 Mpc, label_fdfc beyond, expected count, group jackknife) and the verified receipt of the Ephor_AB
PDF-summary plan be6bd1a9 (photoz_mc, interval bounds). Fixed before the run:
A. Boost of source subsets within the standard background set (photoz_best - z_l > 0.2):
   PDF draw near the lens (|photoz_mc - z_l| <= 0.1), draw in front (photoz_mc <= z_l + 0.1), draw behind (> z_l + 0.2),
   and err95_min <= z_l (95 per cent interval reaches the lens).
B. Alternative source cuts on all selected sources: photoz_mc > z_l + 0.2; err95_min > z_l + 0.1;
   photoz_best > z_l + 0.2 and err95_min > z_l.
C. Fractions of the background set with photoz_mc <= z_l + 0.1 and err95_min <= z_l, lens and randoms, per annulus.
Area noise is not redrawn here (it was <= 0.025 for the standard set).
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
from audit_s16a_outer_excess_split import new_quality, pilot_quality
from audit_s16a_random_extension_boosts import label_fdfc, load_new, load_new_centres, rows
from audit_s16a_random_point_areas import INNER, inner_ok, i_null, load_rows, unit
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from prepare_s16a_photoz_interval_rows import OUTPUT as PZ_PLAN
from s16a_y1_source_contract import finite
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

PZ_SHA = 'be6bd1a95bee49638a5d5337fa4c63f3cc93f8ac77d40a37e3fafe4a558c60cb'
VARIANT = 'null_duplicates + fdfc'


def load_pz():
    res = PZ_PLAN/'results'
    status = json.loads((res/'status.json').read_text())
    summary = json.loads((res/'download_summary.json').read_text())
    if not (status['status'] == 'complete' and summary['approved_scope_plan_sha256'] == PZ_SHA == status['plan_sha256']
            and sha(res/'download_summary.json') == status['summary_sha256'] and summary['credentials_saved'] is False):
        raise ValueError('photo-z interval receipt incomplete or changed')
    for name, digest in summary['output_sha256'].items():
        if sha(res/name) != digest:
            raise ValueError('CSV changed '+name)
    out = {}
    for r in rows(res, sorted(summary['output_sha256'])):
        out[r['object_id']] = r
    return out, sha(res/'download_summary.json')


def num(x):
    v = finite(x)
    return np.nan if v is None else v


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
        q, zl, ok = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens']), ~e['missing']
        for o in np.array(e['ids'])[ok]:
            if o not in pz or pz[o]['ephor_ab_object_id'].strip() == '':
                raise ValueError('no PDF-summary row for selected source '+o)
        g = lambda key, src: np.array([num(src[o][key]) if o in src else np.nan for o in e['ids']])
        arr[cid] = dict(best=g('photoz_best', q)-zl, mc=g('ephor_ab_photoz_mc', pz)-zl, lo95=g('ephor_ab_photoz_err95_min', pz)-zl)
    bgr = lambda a: a['best'] > .2
    subsets = [('A standard background (best > zl+0.2)', bgr),
               ('A  ... and |mc - zl| <= 0.1', lambda a: bgr(a) & (np.abs(a['mc']) <= .1)),
               ('A  ... and mc <= zl + 0.1', lambda a: bgr(a) & (a['mc'] <= .1)),
               ('A  ... and mc > zl + 0.2', lambda a: bgr(a) & (a['mc'] > .2)),
               ('A  ... and err95_min <= zl', lambda a: bgr(a) & (a['lo95'] <= 0)),
               ('B mc > zl + 0.2', lambda a: a['mc'] > .2),
               ('B err95_min > zl + 0.1', lambda a: a['lo95'] > .1),
               ('B best > zl + 0.2 and err95_min > zl', lambda a: bgr(a) & (a['lo95'] > 0))]
    out = []
    for label, f in subsets:
        src = {cid: masked(parts[cid], entries[cid], f(arr[cid])) for cid in parts}
        items = [(src[lens_row[g]['center_id']], [src[r] for r in new_rand[g]]) for g in groups]
        for kind in ('n', 'w'):
            for k in range(bt.NBINS):
                value, used, se, _ = bt.stack_jackknife('expected', items, k, kind, 'area', (), None, False)
                out.append(dict(subset=label, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], boost=value,
                                group_jackknife_se=se, lens_groups=len(used), n_lens=float(sum(lp['n_post'][k] for lp, _ in items))))
    frac = []
    for k in range(bt.NBINS):
        for role in ('lens', 'random'):
            cids = [lens_row[g]['center_id'] for g in groups] if role == 'lens' else [r for g in groups for r in new_rand[g]]
            n = {key: 0 for key in ('bg', 'mc_front', 'lo95_reach')}
            for cid in cids:
                a, e = arr[cid], entries[cid]
                sel = ~e['missing'] & bgr(a) & (e['bins'] == k)
                n['bg'] += int(sel.sum()); n['mc_front'] += int((sel & (a['mc'] <= .1)).sum()); n['lo95_reach'] += int((sel & (a['lo95'] <= 0)).sum())
            frac.append(dict(annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], role=role, n_background=n['bg'],
                             frac_mc_le_zl_plus_0p1=n['mc_front']/n['bg'] if n['bg'] else None,
                             frac_err95min_le_zl=n['lo95_reach']/n['bg'] if n['bg'] else None))
    output.mkdir(parents=True)
    for name, tab in (('boosts.tsv', out), ('fractions.tsv', frac)):
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-pz-dilution-v1', science_measurement=False, selection_change=False, new_query=False,
                   pz_receipt_summary_sha256=pz_sha, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, groups=groups, counts=counts,
                   limitations=['photoz_mc is a single PDF draw per source: subsets by mc are noisy but unbiased in the stack.',
                                'Area noise not redrawn; areas as the standard boost.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('boosts.tsv', 'fractions.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out, frac


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_pz_dilution_20261126_01'))
    out, frac = run(parser.parse_args().output)
    f = lambda x: '—' if x is None else f'{x:.3f}'
    for r in out:
        if r['weighting'] == 'n':
            print(f"{r['subset']:40s} {r['r_lo']}-{r['r_hi']:<4} B={f(r['boost'])} ± {f(r['group_jackknife_se'])} N={r['n_lens']:.0f} g={r['lens_groups']}")
    for r in frac:
        print(f"frac {r['r_lo']}-{r['r_hi']} {r['role']:6s} n={r['n_background']} mc<=zl+0.1 {f(r['frac_mc_le_zl_plus_0p1'])} lo95<=zl {f(r['frac_err95min_le_zl'])}")
