#!/usr/bin/env python3
"""0-0.4 Mpc random-point areas for the 11 lenses and the 76 new random centres, and the inner boosts; local only.

Random points: lens centres from the pilot inner plan (e8b05d6e), new centres from the inner batches of the
source plan (c1dfe5a4). The hole rule, the NULL-row overlap-duplicate probe and the area estimate
(geometric annulus x passing/present points) are those of audit_s16a_random_point_areas. New here: an FDFC
variant that also counts as holes the points outside the FDFC cells of the accepted lens grid (the four wide
lenses use their 180-arcmin grid), as label_fdfc does for its lattice. Sources and counts are the label_fdfc parts
of audit_s16a_random_extension_boosts; the stack is expected count, group jackknife, no own-source rule.
Variants with the dz > 0.2 background cut are included. Binomial area noise is reported per centre only.
"""
from __future__ import annotations

import argparse
import copy
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_afterburner_mask_effects import tangent_arcsec
from audit_s16a_background_cut_boosts import masked
from audit_s16a_fdfc_areas import VALIDATION, in_fdfc
from audit_s16a_outer_excess_split import new_quality, pilot_quality
from audit_s16a_random_extension_boosts import SOURCE_SHA, label_fdfc, load_new, load_new_centres
from audit_s16a_random_point_areas import (DENSITY, DEFICIT_SIGMA, DUP_MIN_FRACTION, DUP_PROBE_ARCSEC, INNER, INNER_BINS,
                                           hole, i_null, inner_ok, load_rows, unit)
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from s16a_y1_source_contract import finite
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

BASE_DZ = 0.2


def new_inner_rows():
    summary = json.loads((SOURCE_PLAN/'results'/'download_summary.json').read_text())
    if summary['approved_scope_plan_sha256'] != SOURCE_SHA or summary.get('credentials_saved') is not False:
        raise ValueError('source receipt differs')
    out = []
    for name, digest in sorted(summary['output_sha256'].items()):
        if 'inner' not in name:
            continue
        path = SOURCE_PLAN/'results'/name
        if sha(path) != digest:
            raise ValueError('csv hash differs '+name)
        with path.open() as stream:
            reader = csv.reader(stream)
            head = [h.lstrip('# ').strip() for h in next(reader)]
            out.extend(dict(zip(head, r)) for r in reader)
    return summary, out


def areas(centre, radius, rows, tree, fdfc_hole):
    chord = 2*math.sin(math.radians(DUP_PROBE_ARCSEC/3600)/2)
    nominal = DENSITY*math.pi*(DUP_PROBE_ARCSEC/60)**2
    ra, dec = np.array([float(r['ra']) for r in rows]), np.array([float(r['dec']) for r in rows])
    x, y = tangent_arcsec(ra, dec, float(centre['ra_deg']), float(centre['dec_deg']))
    mpc = np.hypot(x, y)*3.0/radius
    nulls = np.array([i_null(r) for r in rows], bool)
    dup = np.zeros(len(rows), bool)
    if nulls.any():
        local = np.array([len(v) for v in tree.query_ball_point(unit(ra[nulls], dec[nulls]), chord)])
        dup[np.flatnonzero(nulls)] = local >= DUP_MIN_FRACTION*nominal
    keep = np.array([inner_ok(r) for r in rows], bool) | nulls
    holes = np.array([hole(r) for r in rows], bool)
    out = {}
    for variant, present, hl in (('preregistered', keep, holes), ('null_duplicates', keep & ~dup, holes),
                                 ('null_duplicates + fdfc', keep & ~dup, holes | fdfc_hole(ra, dec))):
        vals = []
        for k in INNER_BINS:
            lo, hi = RADIAL_EDGES[k], RADIAL_EDGES[k+1]
            sel = present & (mpc >= lo) & (mpc < hi)
            a_geo = math.pi*(hi**2-lo**2)*(radius/3.0/60)**2
            n_pres, n_pass = int(sel.sum()), int((sel & ~hl).sum())
            vals.append(dict(annulus=k, a_geo=a_geo, n_present=n_pres, n_pass=n_pass, deficit_sigma=(a_geo*DENSITY-n_pres)/math.sqrt(a_geo*DENSITY),
                             area=a_geo*n_pass/n_pres if n_pres else 0.0,
                             area_binomial_se=a_geo*math.sqrt(n_pass*(n_pres-n_pass)/n_pres**3) if n_pres else None))
        out[variant] = vals
    return out


def with_area(p, vals):
    q = copy.deepcopy(p)
    for v in vals:
        q['a_unm'][v['annulus']] = v['area']
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

    pilot_summary, prow = load_rows(INNER)
    new_summary, nrow = new_inner_rows()
    by_c = defaultdict(list)
    for r in prow+nrow:
        by_c[r['center_id']].append(r)
    pool, seen = [], set()
    for r in prow+nrow:
        if not i_null(r) and inner_ok(r) and r['object_id'] not in seen:
            seen.add(r['object_id']); pool.append((float(r['ra']), float(r['dec'])))
    tree = cKDTree(unit(*map(np.array, zip(*pool))))

    parts, entries, rp, group_of = {}, {}, {}, {}
    for g in groups:
        lens, grid = lens_row[g], fdfc['grids'][g]
        fh = (lambda ra, dec, lens=lens, grid=grid: ~in_fdfc(lens, ra, dec, grid))
        lid = lens['center_id']
        for cid, e, c in [(lid, pilot[lid], lens), *((r, new[r], centres[r]) for r in new_rand[g])]:
            parts[cid], entries[cid], group_of[cid] = label_fdfc(e, lens, grid), e, g
            if not by_c.get(cid):
                raise ValueError('no random points for '+cid)
            rp[cid] = areas(c, float(e['radius']), by_c[cid], tree, fh)

    pq, nq = pilot_quality(), new_quality()
    bg = {}
    for cid, e in entries.items():
        q, zl = (pq if cid.startswith('LC-') else nq), float(e['center']['z_lens'])
        bg[cid] = np.array([(v := finite(q[o]['photoz_best'])) is not None and v-zl > BASE_DZ for o in e['ids']], bool)

    def items(source):
        return [(source[lens_row[g]['center_id']], [source[r] for r in new_rand[g]]) for g in groups]
    sources = {'all': parts, f'dz > {BASE_DZ}': {cid: masked(parts[cid], entries[cid], bg[cid]) for cid in parts}}
    boost_rows = []
    for sname, src in sources.items():
        configs = [('label_fdfc', src)] + [(f'random_point ({v})', {cid: with_area(src[cid], rp[cid][v]) for cid in src})
                                           for v in ('preregistered', 'null_duplicates', 'null_duplicates + fdfc')]
        for aname, source in configs:
            for kind in ('n', 'w'):
                for k in INNER_BINS:
                    value, used, se, _ = bt.stack_jackknife('expected', items(source), k, kind, 'area', (), None, False)
                    boost_rows.append(dict(sources=sname, area=aname, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                           lens_groups=len(used), boost=value, group_jackknife_se=se))
    centre_rows = []
    for cid in parts:
        for v in rp[cid]['null_duplicates + fdfc']:
            k = v['annulus']
            lab = parts[cid]['a_unm'][k]
            centre_rows.append(dict(center_id=cid, role='lens' if cid.startswith('LC-') else 'random', lens_id=group_of[cid], **v,
                                    area_no_fdfc=rp[cid]['null_duplicates'][k]['area'], label_fdfc_area=lab, n_post=parts[cid]['n_post'][k],
                                    random_point_over_label=v['area']/lab if lab > 0 else None,
                                    missing_patch_flag=v['deficit_sigma'] > DEFICIT_SIGMA))
    ratio = defaultdict(list)
    for r in centre_rows:
        if r['random_point_over_label'] is not None:
            ratio[(r['role'], r['annulus'])].append(r['random_point_over_label'])
    ratio_summary = {f'{role}_{k}': dict(n=len(v), median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))
                     for (role, k), v in sorted(ratio.items())}

    output.mkdir(parents=True)
    tables = {'inner_boosts.tsv': boost_rows, 'centre_annulus_areas.tsv': centre_rows}
    for name, tab in tables.items():
        with (output/name).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-new-inner-point-areas-v1', science_measurement=False, selection_change=False, new_query=False,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   pilot_inner_receipt_sha256=sha(INNER/'results'/'download_summary.json'),
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, groups=groups, counts=counts,
                   n_points=dict(pilot=len(prow), new=len(nrow)), base_dz=BASE_DZ,
                   random_point_over_label_area=ratio_summary,
                   missing_patch_flags=sorted({r['center_id'] for r in centre_rows if r['missing_patch_flag']}),
                   limitations=['FDFC enters the random-point areas only geometrically (1-arcmin cells of the accepted grid).',
                                'Binomial area noise is not propagated into the group jackknife.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in tables})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return boost_rows, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_new_inner_point_areas_20261126_01'))
    rows, s = run(parser.parse_args().output)
    for r in rows:
        if r['weighting'] == 'n':
            b, e = r['boost'], r['group_jackknife_se']
            print(f"{r['sources']:9s} {r['area']:38s} {r['r_lo']}-{r['r_hi']:<4} B={'—' if b is None else f'{b:.3f}'} ± {'—' if e is None else f'{e:.3f}'} g={r['lens_groups']}")
    print(json.dumps({k: s[k] for k in ('n_points', 'random_point_over_label_area', 'missing_patch_flags')}, indent=1))
