#!/usr/bin/env python3
"""S16A boosts with the 76 designed random centres (selection 2026-11-26_01); local only.

Inputs are the verified receipts of plans c1dfe5a4 (source pairs) and aeda229e/bce80af8/7d7ed612 (ID supplements).
Each new pair goes through the pilot rules unchanged: Y1 diagnostic selection (select_row table4 pass and finite,
positive numeric calibration), photo-z availability (finite Ephor_AB photoz_best and risk), proper radius at the
lens redshift. Lens centres are the pilot entries. Areas are label_fdfc as in audit_s16a_fdfc_areas (2-arcsec
lattice within 30 arcsec of a source, sources without photo-z masked, lattice points outside FDFC cells of the
accepted lens grid removed; the four wide lenses use their checked 180-arcmin grid). Stack: expected-count,
group jackknife over lens groups. Stacks: new randoms (primary), pilot randoms (reproduction check), both.
Depth/PSF check: per-centre selected-source density per FDFC footprint area and mean weight relative to the lens,
against offset; outer-annulus (1.5-3 Mpc) boost per lens group.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

import audit_s16a_boost_mask_field_tests as bt
from audit_s16a_afterburner_mask_effects import tangent_arcsec
from audit_s16a_fdfc_areas import FOOTPRINT, STEP, VALIDATION, in_fdfc, untangent
from audit_s23b_s16a_quality_supplement import numeric_calibration
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES, radial_bins, separations_radians
from fetch_ancillary_catalogues import catalogue_reader
from prepare_s16a_random_id_supplement import PARTS, output as id_plan
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN, SELECTION
from s16a_y1_source_contract import finite, select_row
from select_s16a_random_centres import load_fdfc, load_wide, merge_wide_fdfc, sha

SOURCE_SHA = 'c1dfe5a4860ceb830281e1444b90b438a99e93702726490518fa32a90073b211'
ID_SHA = {1: 'aeda229e886582c75074c3c58e19bbb727a364fc5dff3fbc44463893cd801fe7',
          2: 'bce80af866047312eeb1d8327a3d6378cc5366ff45a10a13575296f22f73f022',
          3: '7d7ed612bea5678ec115e9d8687a6ce997fd28826c177dbd97b9ef3843ee2cca'}
PAIR_FIELDS = ('center_id', 'object_id', 'source_ra', 'source_dec', 'e1', 'e2', 'sigma_e', 'rms_e', 'weight', 'm', 'c1', 'c2')
OUTER = (6, 7)


def receipt(directory, approved):
    results = directory/'results'
    status = json.loads((results/'status.json').read_text())
    summary = json.loads((results/'download_summary.json').read_text())
    if not (status['status'] == 'complete' and summary['approved_scope_plan_sha256'] == approved
            and sha(results/'download_summary.json') == status['summary_sha256'] and summary['credentials_saved'] is False):
        raise ValueError('receipt incomplete or changed: '+str(directory))
    for name, digest in summary['output_sha256'].items():
        if sha(results/name) != digest:
            raise ValueError('CSV changed '+name)
    return results, summary


def rows(results, names):
    for name in names:
        with (results/name).open(newline='', encoding='utf-8-sig') as stream:
            yield from catalogue_reader(stream)


def load_new_centres():
    summary = json.loads((SELECTION/'summary.json').read_text())
    if sha(SELECTION/'selected_centres.tsv') != summary['selected_centres_sha256']:
        raise ValueError('selection changed')
    with (SELECTION/'selected_centres.tsv').open(newline='') as stream:
        sel = list(csv.DictReader(stream, delimiter='\t'))
    return {r['center_id']: dict(center_id=r['center_id'], role='random', matched_lens_id=r['lens_id'], ra_deg=float(r['ra_deg']),
                                 dec_deg=float(r['dec_deg']), z_lens=float(r['z_lens']), radius_arcsec=float(r['radius_arcsec']),
                                 offset_arcmin=float(r['offset_arcmin']), psf_fwhm=float(r['psf_fwhm']), sky=float(r['sky']),
                                 icountinputs=float(r['icountinputs']), fdfc_fraction=float(r['fdfc_fraction']),
                                 camira_clearance_mpc=float(r['camira_clearance_mpc']))
            for r in sel}, sha(SELECTION/'summary.json')


def load_new(centres):
    results, summary = receipt(SOURCE_PLAN, SOURCE_SHA)
    pairs = [{k: r[k] for k in PAIR_FIELDS} for r in rows(results, sorted(n for n in summary['output_sha256'] if 'pairs' in n))]
    quality, provenance = {}, {'source': sha(results/'download_summary.json')}
    for part in range(1, PARTS+1):
        res, s = receipt(id_plan(part), ID_SHA[part])
        provenance[f'id_part_{part}'] = sha(res/'download_summary.json')
        for r in rows(res, sorted(n for n in s['output_sha256'] if n.startswith('quality_'))):
            if r['object_id'] in quality:
                raise ValueError('duplicate quality row')
            quality[r['object_id']] = r
    status = {}
    for p in pairs:
        q = quality[p['object_id']]
        merged = {**{k: v for k, v in q.items() if k != 'object_id'}, **p}
        status[p['object_id']] = (select_row(merged)['table4_diagnostic'] == 'pass' and numeric_calibration(merged),
                                  finite(q['photoz_best']) is not None and finite(q['photoz_risk_best']) is not None)
    if len(status) != len(quality) or len(status) != len(pairs):
        raise ValueError('pair/quality ID sets differ')
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    by = defaultdict(list)
    for p in pairs:
        if status[p['object_id']][0]:
            by[p['center_id']].append(p)
    data, radius_error = {}, 0.0
    for cid, c in centres.items():
        rs = by[cid]
        ra, dec = np.array([float(p['source_ra']) for p in rs]), np.array([float(p['source_dec']) for p in rs])
        radius = separations_radians(ra, dec, c['ra_deg'], c['dec_deg'])*cosmo.angular_diameter_distance(c['z_lens']).value
        x, y = tangent_arcsec(ra, dec, c['ra_deg'], c['dec_deg'])
        radius_error = max(radius_error, float(np.max(np.abs(np.hypot(x, y)*3.0/c['radius_arcsec']-radius))))
        e = dict(center=c, id=cid, role='random', group=c['matched_lens_id'], radius=c['radius_arcsec'],
                 ids=[p['object_id'] for p in rs], n=len(rs), xy=np.c_[x, y], bins=np.asarray(radial_bins(radius)),
                 weight=np.array([float(p['weight']) for p in rs]), missing=np.array([not status[p['object_id']][1] for p in rs]))
        e['tree'] = cKDTree(e['xy'])
        data[cid] = e
    if radius_error > 5e-3:
        raise ValueError(f'tangent-plane radius differs by {radius_error} Mpc')
    counts = dict(pairs=len(pairs), selected=sum(v[0] for v in status.values()),
                  selected_with_photoz=sum(v[0] and v[1] for v in status.values()))
    return data, counts, provenance, radius_error


def label_fdfc(e, lens, grid):
    lab = bt.labels(e, STEP)
    r, d = untangent(lab['xy'][:, 0], lab['xy'][:, 1], float(e['center']['ra_deg']), float(e['center']['dec_deg']))
    fd = in_fdfc(lens, r, d, grid)
    q = bt.parts(e, lab, FOOTPRINT, e.get('missing', np.zeros(0, bool)), grid_extra=~fd)
    if e['n']:
        q['a_fp'] = np.bincount(lab['bins'][(lab['dist'] <= FOOTPRINT) & fd], minlength=bt.NBINS)*STEP**2/3600
    return q


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fdfc = load_fdfc()
    acc = json.loads((VALIDATION/'user_acceptance.json').read_text())
    if not (acc['accepted'] and acc['grid_summary_sha256'] == fdfc['summary_sha256']):
        raise ValueError('FDFC grid not accepted')
    wide = load_wide()
    wide_check = merge_wide_fdfc(fdfc['grids'], wide)
    centres, selection_sha = load_new_centres()
    new, counts, provenance, radius_error = load_new(centres)
    pilot_centres, pilot, _, _ = bt.load_centres()
    lens_row = {c['matched_lens_id']: c for c in pilot_centres if c['role'] == 'lens'}
    old_rand = defaultdict(list)
    for c in pilot_centres:
        if c['role'] == 'random' and pilot[c['center_id']]['n']:
            old_rand[c['matched_lens_id']].append(c['center_id'])
    new_rand = defaultdict(list)
    for cid, c in centres.items():
        new_rand[c['matched_lens_id']].append(cid)
    groups = sorted(lens_row)
    if sorted(new_rand) != groups:
        raise ValueError('lens groups differ')

    parts = {}
    for g in groups:
        grid, lens = fdfc['grids'][g], lens_row[g]
        for cid in [lens['center_id'], *old_rand[g]]:
            parts[cid] = label_fdfc(pilot[cid], lens, grid)
        for cid in new_rand[g]:
            parts[cid] = label_fdfc(new[cid], lens, grid)

    def items(randoms, subset=groups):
        return [(parts[lens_row[g]['center_id']], [parts[r] for r in randoms[g]]) for g in subset if randoms[g]]
    both = {g: old_rand[g]+new_rand[g] for g in groups}
    stacks = [('new randoms (76)', new_rand), ('pilot non-empty randoms (reproduction)', old_rand), ('pilot + new randoms', both)]
    rows_out = []
    for name, randoms in stacks:
        for kind in ('n', 'w'):
            for k in range(bt.NBINS):
                value, used, se, _ = bt.stack_jackknife('expected', items(randoms), k, kind, 'area', (), None, False)
                rows_out.append(dict(randoms=name, weighting=kind, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                     lens_groups=len(used), boost=value, group_jackknife_se=se))
    per_group = []
    for g in groups:
        for k in range(bt.NBINS):
            value, _ = bt.stack('expected', items(new_rand, [g]), k, 'n', 'area', (), False)
            per_group.append(dict(lens_id=g, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], n_randoms=len(new_rand[g]), boost=value))

    depth = []
    for g in groups:
        lid = lens_row[g]['center_id']
        lp = parts[lid]
        l_area = lp['a_unm'][list(OUTER)].sum()
        l_den = lp['n_post'][list(OUTER)].sum()/l_area if l_area > 0 else np.nan
        l_w = lp['w_post'].sum()/lp['n_post'].sum()
        for cid in [*new_rand[g], *old_rand[g]]:
            p = parts[cid]
            c = centres.get(cid) or dict(offset_arcmin=None, psf_fwhm=None, sky=None, icountinputs=None, camira_clearance_mpc=None)
            area = p['a_unm'].sum()
            depth.append(dict(center_id=cid, lens_id=g, set='new' if cid in centres else 'pilot', offset_arcmin=c['offset_arcmin'],
                              psf_fwhm=c['psf_fwhm'], sky=c['sky'], camira_clearance_mpc=c['camira_clearance_mpc'],
                              icountinputs=c['icountinputs'], n_selected=int(p['n_pre'].sum()), n_with_photoz=int(p['n_post'].sum()),
                              fdfc_unmasked_area_arcmin2=float(area),
                              density_per_arcmin2=float(p['n_post'].sum()/area) if area > 0 else None,
                              density_ratio_to_lens_outer=float(p['n_post'].sum()/area/l_den) if area > 0 and l_den > 0 else None,
                              mean_weight_ratio_to_lens=float(p['w_post'].sum()/p['n_post'].sum()/l_w) if p['n_post'].sum() else None))
    nd = [d for d in depth if d['set'] == 'new']
    dr = np.array([d['density_ratio_to_lens_outer'] for d in nd], float)
    wr = np.array([d['mean_weight_ratio_to_lens'] for d in nd], float)
    off = np.array([d['offset_arcmin'] for d in nd], float)
    clear = np.array([d['camira_clearance_mpc'] for d in nd], float)
    paired = []
    for g in groups:
        a = [d['density_per_arcmin2'] for d in depth if d['lens_id'] == g and d['set'] == 'new' and d['density_per_arcmin2']]
        b = [d['density_per_arcmin2'] for d in depth if d['lens_id'] == g and d['set'] == 'pilot' and d['density_per_arcmin2']]
        if a and b:
            paired.append(dict(lens_id=g, new_over_pilot_mean_density=float(np.mean(a)/np.mean(b)), n_new=len(a), n_pilot=len(b)))
    corr = dict(density_ratio_vs_offset=spearmanr(off, dr, nan_policy='omit')._asdict(),
                weight_ratio_vs_offset=spearmanr(off, wr, nan_policy='omit')._asdict(),
                density_ratio_median=float(np.nanmedian(dr)), density_ratio_range=[float(np.nanmin(dr)), float(np.nanmax(dr))],
                weight_ratio_median=float(np.nanmedian(wr)), weight_ratio_range=[float(np.nanmin(wr)), float(np.nanmax(wr))],
                density_ratio_vs_camira_clearance=spearmanr(clear, dr, nan_policy='omit')._asdict(),
                new_over_pilot_random_density=paired)

    output.mkdir(parents=True)
    for fname, tab in (('boosts.tsv', rows_out), ('per_group_boosts.tsv', per_group), ('random_depth_check.tsv', depth)):
        with (output/fname).open('w', newline='') as stream:
            fields = list(dict.fromkeys(k for r in tab for k in r))
            w = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary = dict(schema='astrafek-s16a-random-extension-boosts-v1', science_measurement=False, selection_change=False,
                   selection_summary_sha256=selection_sha, receipt_summary_sha256=provenance,
                   fdfc_grid_summary_sha256=fdfc['summary_sha256'], wide_fdfc_check=wide_check, counts=counts,
                   tangent_radius_max_error_mpc=radius_error, step_arcsec=STEP, footprint_arcsec=FOOTPRINT, groups=groups,
                   depth_correlation=corr, outer_annuli=list(OUTER),
                   limitations=['label_fdfc areas only; the inner random-point FDFC area (0-0.4 Mpc) for the new centres is not computed here.',
                                'FDFC cells are 1 arcmin; the 0-0.1 Mpc annulus inherits cell-boundary coarseness.',
                                'Density ratio uses the lens 1.5-3 Mpc annuli as reference, which still contain some cluster members.'],
                   runner_sha256=sha(Path(__file__)),
                   output_sha256={n: sha(output/n) for n in ('boosts.tsv', 'per_group_boosts.tsv', 'random_depth_check.tsv')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return rows_out, corr, counts


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_random_extension_boosts_20261126_02'))
    rows_out, corr, counts = run(parser.parse_args().output)
    print(json.dumps(counts))
    for r in rows_out:
        if r['weighting'] == 'n':
            print(f"{r['randoms']:42s} {r['r_lo']:4}-{r['r_hi']:<4} B={r['boost']:.3f} ± {r['group_jackknife_se']:.3f} (g={r['lens_groups']})")
    print(json.dumps(corr, indent=1, default=float))
