#!/usr/bin/env python3
"""S16A standard boost definition (user approval 2026-11-26) with area binomial noise propagated; local only.

Definition:
- sources: Y1 diagnostic selection with numeric calibration, finite Ephor_AB photoz_best and risk, and the
  background cut photoz_best - z_lens > 0.2;
- randoms: the 76 designed random centres (selection 2026-11-26_01), each matched to its lens;
- areas: 0-0.4 Mpc random-point areas with FDFC (null-duplicate rule, audit_s16a_new_inner_point_areas);
  0.4-3 Mpc label_fdfc (2-arcsec lattice, 30-arcsec footprint, FDFC cells of the accepted grid);
- stack: expected count over lens groups, group jackknife.
Area noise: each random-point area is redrawn as a_geo x Binomial(n_present, n_pass/n_present)/n_present, for
all centres at once, R realisations with a fixed seed; the stack spread over realisations is added in quadrature
to the group jackknife. Also reported: the same with all sources (no background cut).
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
from audit_s16a_random_extension_boosts import label_fdfc, load_new, load_new_centres
from audit_s16a_random_point_areas import INNER, inner_ok, i_null, load_rows, unit
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from s16a_y1_source_contract import finite
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

BASE_DZ = 0.2
VARIANT = 'null_duplicates + fdfc'
REALISATIONS = 2000
SEED = 2026112601


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
    parts, entries, rp = {}, {}, {}
    for g in groups:
        lens, grid = lens_row[g], fdfc['grids'][g]
        fh = (lambda ra, dec, lens=lens, grid=grid: ~in_fdfc(lens, ra, dec, grid))
        lid = lens['center_id']
        for cid, e, c in [(lid, pilot[lid], lens), *((r, new[r], centres[r]) for r in new_rand[g])]:
            parts[cid], entries[cid] = label_fdfc(e, lens, grid), e
            rp[cid] = areas(c, float(e['radius']), by_c[cid], tree, fh)[VARIANT]
    pq, nq = pilot_quality(), new_quality()
    bg = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        bg[cid] = np.array([(v := finite(q[o]['photoz_best'])) is not None and v-zl > BASE_DZ for o in e['ids']], bool)
    srcsets = {'background dz > 0.2 (standard)': {cid: masked(parts[cid], entries[cid], bg[cid]) for cid in parts},
               'all sources': {cid: masked(parts[cid], entries[cid], np.ones(entries[cid]['n'], bool)) for cid in parts}}

    def items(src):
        return [(src[lens_row[g]['center_id']], [src[r] for r in new_rand[g]]) for g in groups]
    rng = np.random.default_rng(SEED)
    draws = {cid: [(v['annulus'], v['a_geo'], v['n_present'], v['n_pass']) for v in vals] for cid, vals in rp.items()}
    realised = []
    for _ in range(REALISATIONS):
        realised.append({cid: [dict(annulus=k, area=(a*rng.binomial(n, m/n)/n if n else 0.0)) for k, a, n, m in d] for cid, d in draws.items()})
    out = []
    for sname, src in srcsets.items():
        std = {cid: with_area(src[cid], rp[cid]) for cid in src}
        for kind in ('n', 'w'):
            for k in range(bt.NBINS):
                value, used, jk, _ = bt.stack_jackknife('expected', items(std), k, kind, 'area', (), None, False)
                mc = None
                if k <= 2:
                    vals = []
                    for real in realised:
                        s = {cid: with_area(src[cid], real[cid]) for cid in src}
                        v = bt.stack('expected', items(s), k, kind, 'area', ())[0]
                        if v is not None:
                            vals.append(v)
                    mc = float(np.std(vals, ddof=1))
                tot = None if jk is None else float(np.hypot(jk, mc or 0.0))
                out.append(dict(sources=sname, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                area='random_point + fdfc' if k <= 2 else 'label_fdfc', lens_groups=len(used), boost=value,
                                group_jackknife_se=jk, area_binomial_se=mc, total_se=tot,
                                n_lens=float(sum(lp['n_post'][k] for lp, _ in items(std)))))
    output.mkdir(parents=True)
    with (output/'standard_boosts.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(out[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(out)
    summary = dict(schema='astrafek-s16a-standard-boost-v1', science_measurement=False, new_query=False,
                   definition=dict(sources='Y1 diagnostic + numeric calibration + finite Ephor_AB photoz_best/risk + photoz_best - z_lens > 0.2',
                                   randoms='76 designed random centres, selection 2026-11-26_01',
                                   areas='0-0.4 Mpc random-point + FDFC (null-duplicate rule); 0.4-3 Mpc label_fdfc',
                                   stack='expected count, group jackknife', approved='user, 2026-11-26'),
                   realisations=REALISATIONS, seed=SEED, selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, groups=groups, counts=counts,
                   limitations=['Point photo-z background cut; no P(z) or lensing efficiency.',
                                'Area noise and group jackknife are combined in quadrature (treated as independent).',
                                'FDFC enters the random-point areas through 1-arcmin cells; 0-0.2 Mpc is supported by 7 groups.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={'standard_boosts.tsv': sha(output/'standard_boosts.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_standard_boost_20261126_01'))
    for r in run(parser.parse_args().output):
        if r['weighting'] == 'n' or r['annulus'] <= 2:
            f = lambda x: '—' if x is None else f'{x:.3f}'
            print(f"{r['sources'][:10]:10s} {r['weighting']} {r['r_lo']}-{r['r_hi']:<4} B={f(r['boost'])} jk={f(r['group_jackknife_se'])} "
                  f"area={f(r['area_binomial_se'])} tot={f(r['total_se'])} g={r['lens_groups']}")
