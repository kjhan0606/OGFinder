#!/usr/bin/env python3
"""S16A boosts with FDFC-restricted areas (grid 019c8228 accepted by the user 2026-11-25); local only.

Areas per annulus, 0-3 Mpc, 66 pilot centres:
- label: 2-arcsec lattice within 30 arcsec of a source, unmasked by the recorded S16A afterburner rule;
- label_fdfc: the same lattice points that also lie in an FDFC cell of the lens grid;
- random_point_fdfc (0-0.4 Mpc): inner random-point area with points outside FDFC cells counted as holes.
Stacks are the expected-count stack over the supported lens groups (own-source rule only for 'label'), with a
group jackknife; randoms are the non-empty ones, and for the FDFC areas also all randoms (empty ones get zero
FDFC area and so drop out by themselves).
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np

import audit_s16a_boost_mask_field_tests as bt
import audit_s16a_random_point_areas as ra
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from select_s16a_random_centres import load_fdfc, sha

VALIDATION = Path('/home/kjhan/BACKUP/Astrafek/s16a_fdfc_grid_validation_20261125_01')
STEP, FOOTPRINT = 2.0, 30.0


def accepted_grid():
    fdfc = load_fdfc()
    acc = json.loads((VALIDATION/'user_acceptance.json').read_text())
    if acc['grid_summary_sha256'] != fdfc['summary_sha256'] or not acc['accepted']:
        raise ValueError('FDFC grid not accepted')
    return fdfc


def in_fdfc(lens, ra_deg, dec_deg, grid):
    ra0, dec0 = float(lens['ra_deg']), float(lens['dec_deg'])
    cx = np.floor((np.asarray(ra_deg)-ra0)*math.cos(math.radians(dec0))*60).astype(int)
    cy = np.floor((np.asarray(dec_deg)-dec0)*60).astype(int)
    return np.array([grid.get((i, j), False) for i, j in zip(cx, cy)], bool)


def untangent(x, y, ra0, dec0):
    x, y = np.radians(np.asarray(x)/3600), np.radians(np.asarray(y)/3600)
    a0, d0 = math.radians(ra0), math.radians(dec0)
    rho = np.hypot(x, y); c = np.arctan(rho)
    with np.errstate(invalid='ignore', divide='ignore'):
        dec = np.arcsin(np.cos(c)*math.sin(d0)+np.where(rho > 0, y*np.sin(c)*math.cos(d0)/rho, 0))
    ra = a0+np.arctan2(x*np.sin(c), rho*math.cos(d0)*np.cos(c)-y*math.sin(d0)*np.sin(c))
    return np.degrees(ra), np.degrees(dec)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fdfc = accepted_grid()
    centers, data, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in centers if c['role'] == 'lens'}
    grids = {g: fdfc['grids'][g] for g in lens_row}

    label, label_fd, frac = {}, {}, []
    for c in centers:
        e, g = data[c['center_id']], c['matched_lens_id']
        lab = bt.labels(e, STEP)
        r, d = untangent(lab['xy'][:, 0], lab['xy'][:, 1], float(c['ra_deg']), float(c['dec_deg']))
        fd = in_fdfc(lens_row[g], r, d, grids[g])
        missing = e.get('missing', np.zeros(0, bool))
        p = bt.parts(e, lab, FOOTPRINT, missing)
        q = bt.parts(e, lab, FOOTPRINT, missing, grid_extra=~fd)
        if e['n']:
            q['a_fp'] = np.bincount(lab['bins'][(lab['dist'] <= FOOTPRINT) & fd], minlength=bt.NBINS)*STEP**2/3600
        label[c['center_id']], label_fd[c['center_id']] = p, q
        geo = np.bincount(lab['bins'], minlength=bt.NBINS).astype(float)
        frac.append(dict(center_id=c['center_id'], role=c['role'], lens_id=g, n_sources=e['n'],
                         **{f'fdfc_{k}': float(np.bincount(lab['bins'][fd], minlength=bt.NBINS)[k]/geo[k]) for k in range(bt.NBINS)},
                         **{f'label_unm_{k}': float(p['a_unm'][k]) for k in range(bt.NBINS)},
                         **{f'label_fdfc_unm_{k}': float(q['a_unm'][k]) for k in range(bt.NBINS)}))

    def extra(c, r):
        return not in_fdfc(lens_row[c['matched_lens_id']], [float(r['ra'])], [float(r['dec'])], grids[c['matched_lens_id']])[0]
    _, areas, _, _, _ = ra.part_areas(centers, data, extra)
    rp = {cid: ra.with_area(ra.label_parts(data[cid]), areas[cid]['null_duplicates']) for cid in data}

    lenses = {g: c['center_id'] for g, c in lens_row.items()}
    rand = {}
    for c in centers:
        if c['role'] == 'random':
            rand.setdefault(c['matched_lens_id'], []).append(c['center_id'])
    nonempty = {g: sorted(r for r in rs if data[r]['n']) for g, rs in rand.items()}
    allr = {g: sorted(rs) for g, rs in rand.items()}
    groups = sorted(g for g in lenses if sum(label[r]['a_fp'].sum() for r in rand[g]) > 0)

    def items(src, randoms):
        return [(src[lenses[g]], [src[r] for r in randoms[g]]) for g in groups]
    configs = [('label (2 arcsec, own rule)', label, True, nonempty, range(bt.NBINS)),
               ('label_fdfc', label_fd, False, nonempty, range(bt.NBINS)),
               ('label_fdfc + empty randoms', label_fd, False, allr, range(bt.NBINS)),
               ('random_point_fdfc', rp, False, nonempty, ra.INNER_BINS),
               ('random_point_fdfc + empty randoms', rp, False, allr, ra.INNER_BINS)]
    rows = []
    for name, src, own, randoms, ks in configs:
        for kind in ('n', 'w'):
            for k in ks:
                value, used, se, _ = bt.stack_jackknife('expected', items(src, randoms), k, kind, 'area', (), None, own)
                rows.append(dict(area=name, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                 lens_groups=len(used), boost=value, group_jackknife_se=se))
    output.mkdir(parents=True)
    for fname, tab in (('boosts.tsv', rows), ('centre_fdfc_areas.tsv', frac)):
        with (output/fname).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-fdfc-areas-v1', science_measurement=False, selection_change=False,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], acceptance_sha256=sha(VALIDATION/'user_acceptance.json'),
                   step_arcsec=STEP, footprint_arcsec=FOOTPRINT, groups=groups,
                   limitations=['FDFC cells are 1 arcmin; inner annuli (<0.1 Mpc ~ 0.3 arcmin) inherit cell-boundary coarseness.',
                                'Label lattice at 2 arcsec, not the recorded 1 arcsec.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('boosts.tsv', 'centre_fdfc_areas.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return rows


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_fdfc_areas_20261125_01'))
    for r in run(parser.parse_args().output):
        if r['weighting'] == 'n':
            print(f"{r['area']:36s} {r['r_lo']:4}-{r['r_hi']:<4} B={r['boost']:.3f} ± {r['group_jackknife_se']:.3f} (g={r['lens_groups']})")
