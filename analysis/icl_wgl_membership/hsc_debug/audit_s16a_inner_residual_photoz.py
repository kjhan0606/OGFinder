#!/usr/bin/env python3
"""Is the S16A central (0-0.4 Mpc) boost left after the dz > 0.2 background cut photo-z leakage? Local only.

Setup as audit_s16a_background_cut_boosts (label_fdfc areas, expected count, group jackknife, 76 new randoms).
Base set: selected, photo-z-available sources with Ephor_AB photoz_best - z_lens > 0.2. Subsets, fixed before
the run (the risk split point is the median over the base set of the new randoms):
  risk low / high (photoz_risk_best <= / > median);
  Mizuki also background (mizuki_photoz_best - z_lens > 0.2) / Mizuki not background (finite, <= 0.2) / no Mizuki;
  strict: risk low and Mizuki also background.
Mizuki values come from the cmodel photo-z receipts (pilot plan 35265ee4, new ID supplements).
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_background_cut_boosts import masked
from audit_s16a_cmodel_photoz_missingness import load_cmodel
from audit_s16a_fdfc_areas import VALIDATION
from audit_s16a_outer_excess_split import new_quality, pilot_quality
from audit_s16a_random_extension_boosts import ID_SHA, label_fdfc, load_new, load_new_centres, receipt, rows
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from prepare_s16a_random_id_supplement import PARTS, output as id_plan
from s16a_y1_source_contract import finite
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

BASE_DZ = 0.2


def new_mizuki():
    out = {}
    for part in range(1, PARTS+1):
        res, s = receipt(id_plan(part), ID_SHA[part])
        for r in rows(res, sorted(n for n in s['output_sha256'] if n.startswith('cmodel_'))):
            out[r['object_id']] = finite(r['mizuki_photoz_best']) if r['mizuki_object_id'].strip() else None
    return out


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
    pq, nq = pilot_quality(), new_quality()
    pcm, cm_prov = load_cmodel()
    pmiz = {o: (r['mizuki_photoz_best'] if r['has_mizuki'] else None) for o, r in pcm.items()}
    nmiz = new_mizuki()
    parts, entries = {}, {}
    for g in groups:
        lid = lens_row[g]['center_id']
        parts[lid], entries[lid] = label_fdfc(pilot[lid], lens_row[g], fdfc['grids'][g]), pilot[lid]
        for cid in new_rand[g]:
            parts[cid], entries[cid] = label_fdfc(new[cid], lens_row[g], fdfc['grids'][g]), new[cid]
    arr = {}
    for cid, e in entries.items():
        lens = cid.startswith('LC-')
        q, mz, zl = (pq if lens else nq), (pmiz if lens else nmiz), float(e['center']['z_lens'])
        def f(x):
            v = finite(x)
            return np.nan if v is None else v
        arr[cid] = dict(dz=np.array([f(q[o]['photoz_best']) for o in e['ids']])-zl,
                        risk=np.array([f(q[o]['photoz_risk_best']) for o in e['ids']]),
                        dzm=np.array([np.nan if mz.get(o) is None else mz[o] for o in e['ids']], float)-zl)
    base = {cid: a['dz'] > BASE_DZ for cid, a in arr.items()}
    pool = np.concatenate([arr[c]['risk'][base[c] & ~entries[c]['missing']] for g in groups for c in new_rand[g]])
    rmed = float(np.nanmedian(pool))
    subsets = [(f'dz > {BASE_DZ} (base)', lambda a: a['dz'] > BASE_DZ),
               (f'base, risk <= {rmed:.4f}', lambda a: (a['dz'] > BASE_DZ) & (a['risk'] <= rmed)),
               (f'base, risk > {rmed:.4f}', lambda a: (a['dz'] > BASE_DZ) & (a['risk'] > rmed)),
               ('base, Mizuki also background', lambda a: (a['dz'] > BASE_DZ) & (a['dzm'] > BASE_DZ)),
               ('base, Mizuki not background', lambda a: (a['dz'] > BASE_DZ) & (a['dzm'] <= BASE_DZ)),
               ('base, no Mizuki', lambda a: (a['dz'] > BASE_DZ) & np.isnan(a['dzm'])),
               ('strict: risk low and Mizuki background', lambda a: (a['dz'] > BASE_DZ) & (a['risk'] <= rmed) & (a['dzm'] > BASE_DZ))]
    out = []
    for label, fn in subsets:
        sub = {cid: masked(parts[cid], entries[cid], fn(arr[cid])) for cid in entries}
        items = [(sub[lens_row[g]['center_id']], [sub[r] for r in new_rand[g]]) for g in groups]
        for kind in ('n', 'w'):
            for k in range(bt.NBINS):
                value, used, se, _ = bt.stack_jackknife('expected', items, k, kind, 'area', (), None, False)
                out.append(dict(subset=label, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], boost=value,
                                group_jackknife_se=se, lens_groups=len(used), n_lens=float(sum(lp['n_post'][k] for lp, _ in items)),
                                n_random=float(sum(r['n_post'][k] for _, rs in items for r in rs))))
    output.mkdir(parents=True)
    with (output/'boosts.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)
    summary = dict(schema='astrafek-s16a-inner-residual-photoz-v1', science_measurement=False, selection_change=False, new_query=False,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, pilot_cmodel=cm_prov,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, base_dz=BASE_DZ,
                   risk_median_new_randoms=rmed, groups=groups, counts=counts,
                   limitations=['Point photo-z values only; no P(z).', 'Areas are label_fdfc of the full selected sample.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'boosts.tsv': sha(output/'boosts.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_inner_residual_photoz_20261126_01'))
    for r in run(parser.parse_args().output):
        if r['weighting'] == 'n':
            b, s = r['boost'], r['group_jackknife_se']
            print(f"{r['subset']:42s} {r['r_lo']}-{r['r_hi']:<4} B={'—' if b is None else f'{b:.3f}'} ± {'—' if s is None else f'{s:.3f}'} "
                  f"Nl={r['n_lens']:.0f} Nr={r['n_random']:.0f}")
