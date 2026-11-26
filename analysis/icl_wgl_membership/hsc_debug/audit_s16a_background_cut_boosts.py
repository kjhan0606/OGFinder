#!/usr/bin/env python3
"""S16A boost profiles after photo-z background cuts, lens vs the 76 new randoms; local only.

Same entries, label_fdfc areas and expected-count stack (group jackknife) as audit_s16a_outer_excess_split;
the counted sources are the selected, photo-z-available sources with photoz_best - z_lens > dz_min, for
dz_min in DZ_CUTS fixed before the run, plus dz > 0.4 without the red g-i tertile. Counts (n) and weights (w).
Descriptive: the areas are those of the full selected sample; no lensing efficiency or new query.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_fdfc_areas import VALIDATION
from audit_s16a_outer_excess_split import new_quality, pilot_quality, props
from audit_s16a_random_extension_boosts import label_fdfc, load_new, load_new_centres
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

DZ_CUTS = (None, .1, .2, .3, .4, .6)


def masked(p, e, mask):
    q, keep = dict(p), ~e['missing'] & mask
    for kind, w in (('n', np.ones(e['n'])), ('w', e['weight'])):
        q[kind+'_post'] = np.bincount(e['bins'][keep], weights=w[keep], minlength=bt.NBINS)
    return q


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
    parts, entries = {}, {}
    for g in groups:
        lid = lens_row[g]['center_id']
        parts[lid], entries[lid] = label_fdfc(pilot[lid], lens_row[g], fdfc['grids'][g]), pilot[lid]
        for cid in new_rand[g]:
            parts[cid], entries[cid] = label_fdfc(new[cid], lens_row[g], fdfc['grids'][g]), new[cid]
    arr = {}
    for cid, e in entries.items():
        P = [props((pq if cid.startswith('LC-') else nq)[o], float(e['center']['z_lens'])) for o in e['ids']]
        arr[cid] = {k: np.array([x[k] for x in P], float) for k in ('dz', 'gi')}
    pool = np.concatenate([arr[c]['gi'][~entries[c]['missing']] for g in groups for c in new_rand[g]])
    t2 = float(np.nanpercentile(pool, 200/3))
    cuts = [('all' if d is None else f'dz > {d}', (lambda a, d=d: np.ones(len(a['dz']), bool) if d is None else a['dz'] > d))
            for d in DZ_CUTS]
    cuts.append((f'dz > 0.4 and g-i <= {t2:.3f}', lambda a: (a['dz'] > .4) & ~(a['gi'] > t2)))
    out = []
    for label, f in cuts:
        sub = {cid: masked(parts[cid], entries[cid], f(arr[cid])) for cid in entries}
        items = [(sub[lens_row[g]['center_id']], [sub[r] for r in new_rand[g]]) for g in groups]
        for kind in ('n', 'w'):
            for k in range(bt.NBINS):
                value, used, se, _ = bt.stack_jackknife('expected', items, k, kind, 'area', (), None, False)
                out.append(dict(cut=label, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], boost=value,
                                group_jackknife_se=se, lens_groups=len(used), n_lens=float(sum(lp['n_post'][k] for lp, _ in items)),
                                n_random=float(sum(r['n_post'][k] for _, rs in items for r in rs))))
    output.mkdir(parents=True)
    with (output/'boosts.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)
    summary = dict(schema='astrafek-s16a-background-cut-boosts-v1', science_measurement=False, selection_change=False, new_query=False,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance, fdfc_grid_summary_sha256=fdfc['summary_sha256'],
                   wide_fdfc_check=wide_check, dz_cuts=list(DZ_CUTS), red_tertile_from_new_randoms=t2, counts=counts, groups=groups,
                   limitations=['Point photo-z (Ephor_AB photoz_best) cut; no P(z), no lensing efficiency.',
                                'Areas are label_fdfc of the full selected sample.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'boosts.tsv': sha(output/'boosts.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_background_cut_boosts_20261126_01'))
    for r in run(parser.parse_args().output):
        if r['weighting'] == 'n':
            b, s = r['boost'], r['group_jackknife_se']
            print(f"{r['cut']:30s} {r['r_lo']}-{r['r_hi']:<4} B={'—' if b is None else f'{b:.3f}'} ± {'—' if s is None else f'{s:.3f}'} "
                  f"Nl={r['n_lens']:.0f} g={r['lens_groups']}")
