#!/usr/bin/env python3
"""Offline diagnostics of the S16A afterburner (Ephor_AB availability) source mask.

1. Per-annulus boosts before the mask, after it without and with mask-area normalisation.
   Mask area comes from a 4-arcsec grid labelled by the nearest selected source.
2. Same-cell colour proxy test of redshift-independent missingness (cmodel colours).
3. Near-threshold survivors: photo-z multiplier change and same-cell z_best test after
   removing sources whose afterburner PSF FWHM exceeds 1.0 arcsec.
4. Lens jackknife and mask-transplant null tests for the boost and the multiplier.

Reads verified local receipts only. Descriptive sensitivities: no shear, mass, correction or
selection change. All thresholds and seeds are fixed before looking at the outputs.
"""
import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_s16a_cosmos_lens_response import INTAKE, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from audit_s16a_missingness_geometry import CELL_DEG, MIN_CELL
from audit_s16a_y1_reuse_shape_noise import GEOMETRY, GEOMETRY_SHA
from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from fetch_ancillary_catalogues import catalogue_reader
from prepare_s16a_afterburner_photoz_rows import BANDS, OUTPUT as AB_PLAN
from submit_s16a_afterburner_photoz_rows import load_plan as load_ab_plan, verify_complete as verify_ab
from ogfkit.wlcluster import Cosmology

CALIB = PLAN.parent/'s16a_cosmos_bin_support_2026-11-21_01'
PIXEL_ARCSEC = 0.168
FWHM_PER_SIGMA = 2*math.sqrt(2*math.log(2))
NEAR_FWHM_ARCSEC = 1.0
GAP = 0.2
GRID_ARCSEC = 4.0
FOOTPRINTS_ARCSEC = (20.0, 30.0, 45.0)
PRIMARY_FOOTPRINT = 30.0
BOOTSTRAP = 1000
THINNING = 500
PERMUTATIONS = 2000
SEED = 2026112401
DONOR_LENS = 'S23B_3620'
WINSOR = (1, 99)
# Added after the trial run: the near-threshold cut removed 94-99.9% of the weight in whole lens
# fields, where a per-lens change is not defined. Below this retained fraction it is not reported.
MIN_RETAINED_WEIGHT = 0.5


def num(value):
    value = value.strip()
    return math.nan if value == '' else float(value)


def tangent_arcsec(ra, dec, ra0, dec0):
    ra, dec = np.radians(np.asarray(ra, float)), np.radians(np.asarray(dec, float))
    ra0, dec0 = math.radians(ra0), math.radians(dec0)
    cosc = math.sin(dec0)*np.sin(dec)+math.cos(dec0)*np.cos(dec)*np.cos(ra-ra0)
    x = np.cos(dec)*np.sin(ra-ra0)/cosc
    y = (math.cos(dec0)*np.sin(dec)-math.sin(dec0)*np.cos(dec)*np.cos(ra-ra0))/cosc
    return np.degrees(x)*3600, np.degrees(y)*3600


def grid(radius_arcsec):
    n = int(math.ceil(radius_arcsec/GRID_ARCSEC))
    axis = np.arange(-n, n+1)*GRID_ARCSEC
    x, y = np.meshgrid(axis, axis)
    keep = x**2+y**2 <= radius_arcsec**2
    return x[keep], y[keep]


def annulus(radius_mpc):
    return np.minimum(np.searchsorted(RADIAL_EDGES, radius_mpc, side='right')-1, len(RADIAL_EDGES)-2)


def ratio(nl, al, nr, ar):
    return (nl/al)/(nr/ar) if al > 0 and ar > 0 and nr > 0 and nl > 0 else None


def jackknife(values):
    values = np.asarray([v for v in values if v is not None], float)
    n = len(values)
    return float(math.sqrt((n-1)/n*((values-values.mean())**2).sum())) if n > 1 else None


def load_afterburner():
    plan = load_ab_plan(AB_PLAN)
    results = AB_PLAN/'results'
    state = json.loads((results/'status.json').read_text())
    receipt = verify_ab(plan, AB_PLAN, state)
    if state['status'] != 'complete' or sha(results/'download_summary.json') != state['summary_sha256']:
        raise ValueError('afterburner receipt is incomplete or changed')
    seeing, unavailable = {}, set()
    for name in receipt['output_sha256']:
        with (results/name).open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                seeing[row['object_id']] = [num(row[b+'parent_flux_convolved_seeing']) for b in BANDS]
                if not all(math.isfinite(num(row[b+'parent_flux_convolved_2_1'])) for b in BANDS):
                    unavailable.add(row['object_id'])
    return seeing, unavailable, dict(afterburner_plan_sha256=plan['query_plan_sha256'],
                                     afterburner_receipt_sha256=sha(results/'download_summary.json'))


def load_supplement():
    receipt = json.loads((PLAN/'results'/'download_summary.json').read_text())
    props = {}
    for name, digest in receipt['output_sha256'].items():
        path = PLAN/'results'/name
        if sha(path) != digest:
            raise ValueError('supplement digest mismatch')
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                props[row['object_id']] = row
    return props


def calibration(cosmo):
    intake = json.loads((INTAKE/'summary.json').read_text())
    cal = json.loads((CALIB/'summary.json').read_text())
    for folder, report in ((INTAKE, intake), (CALIB, cal)):
        for name, expected in report['output_sha256'].items():
            if sha(folder/name) != expected:
                raise ValueError('changed calibration input '+name)
    if cal['intake_sha256'] != sha(INTAKE/'summary.json'):
        raise ValueError('COSMOS lineage mismatch')
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    if len(rows) != intake['eligible_rows']:
        raise ValueError('COSMOS eligibility mismatch')
    photo, ref, shape, som = np.asarray([[float(r[k]) for k in ('photo_z', 'reference_z', 'shape_weight', 'som_weight')]
                                         for r in rows]).T
    with (CALIB/'lens_photoz_reweight.tsv').open() as stream:
        baseline = {r['lens_id']: float(r['reweighted']) for r in csv.DictReader(stream, delimiter='\t')}
    edges = np.asarray(cal['edges'])

    def means(zl):
        chi = float(cosmo.comoving(zl))
        keep = photo > zl+GAP
        p, r = photo[keep], ref[keep]
        est = 1-chi/cosmo.comoving(p)
        true = np.zeros_like(r)
        true[r > zl] = 1-chi/cosmo.comoving(r[r > zl])
        w = shape[keep]*som[keep]*est**2
        index = np.searchsorted(edges, p, side='right')-1
        mass = np.bincount(index, weights=w, minlength=len(edges)-1)
        total = np.bincount(index, weights=w*true/est, minlength=len(edges)-1)
        return np.divide(total, mass, out=np.zeros_like(mass), where=mass > 0), mass
    return means, edges, baseline, dict(intake_sha256=sha(INTAKE/'summary.json'), calibration_sha256=sha(CALIB/'summary.json'))


def stratified_test(values, group, label, rng):
    """Stratum-size-weighted mean difference (label True minus False) of winsorised values,
    with a within-stratum label permutation null."""
    ok = np.isfinite(values)
    values, group, label = values[ok], group[ok], label[ok]
    lo, hi = np.percentile(values, WINSOR)
    values = np.clip(values, lo, hi)
    _, group = np.unique(group, return_inverse=True)
    both = np.bincount(group, weights=label, minlength=group.max()+1)
    sizes = np.bincount(group, minlength=group.max()+1)
    usable = (both > 0) & (both < sizes)
    keep = usable[group]
    values, group, label = values[keep], group[keep], label[keep]
    _, group = np.unique(group, return_inverse=True)
    sizes = np.bincount(group).astype(float)

    def statistic(lab):
        n1 = np.bincount(group, weights=lab, minlength=len(sizes))
        s1 = np.bincount(group, weights=values*lab, minlength=len(sizes))
        s0 = np.bincount(group, weights=values*(~lab), minlength=len(sizes))
        return float((sizes*(s1/n1-s0/(sizes-n1))).sum()/sizes.sum())
    observed = statistic(label)
    order = np.argsort(group, kind='stable')
    null = np.empty(PERMUTATIONS)
    for k in range(PERMUTATIONS):
        shuffled = np.empty_like(label)
        shuffled[order] = label[np.lexsort((rng.random(len(label)), group))]
        null[k] = statistic(shuffled)
    a, b = values[label], values[~label]
    return dict(n_true=int(label.sum()), n_false=int((~label).sum()), n_strata=len(sizes),
                median_true=float(np.median(a)), median_false=float(np.median(b)),
                iqr_true=float(np.subtract(*np.percentile(a, [75, 25]))), iqr_false=float(np.subtract(*np.percentile(b, [75, 25]))),
                stratified_mean_difference=observed, null_sd=float(null.std(ddof=1)),
                permutation_p=float((1+(np.abs(null) >= abs(observed)).sum())/(1+PERMUTATIONS)))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    rng = np.random.default_rng(SEED)
    cosmo = Cosmology()
    seeing, unavailable, ab_provenance = load_afterburner()
    props = load_supplement()
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    point_z = finite_point_redshifts(available)
    if set(seeing) != set(props) or {o for o in selected if o not in available} != unavailable & selected:
        raise ValueError('afterburner mask differs from the photo-z availability ledger')
    if sha(GEOMETRY/'summary.json') != GEOMETRY_SHA:
        raise ValueError('geometry summary changed')
    geometry = json.loads((GEOMETRY/'summary.json').read_text())
    if sha(GEOMETRY/'pair_geometry.tsv') != geometry['output_sha256']['pair_geometry.tsv']:
        raise ValueError('geometry table changed')
    means_for, edges, baseline, cal_provenance = calibration(cosmo)

    by_center = defaultdict(list)
    with (GEOMETRY/'pair_geometry.tsv').open() as stream:
        for pair, geo in zip(pairs, csv.DictReader(stream, delimiter='\t'), strict=True):
            if (pair['center_id'], pair['object_id']) != (geo['center_id'], geo['object_id']):
                raise ValueError('pair/geometry mismatch')
            if pair['object_id'] in selected:
                by_center[pair['center_id']].append((pair, float(geo['radius_proper_mpc']), int(geo['radial_bin'])))
    fwhm = {o: np.asarray(s)*FWHM_PER_SIGMA*PIXEL_ARCSEC for o, s in seeing.items()}

    data, radius_error = {}, 0.0
    for c in centers:
        rows = by_center[c['center_id']]
        ids = [p['object_id'] for p, _, _ in rows]
        entry = dict(center=c, ids=ids, n=len(ids))
        if ids:
            x, y = tangent_arcsec([float(p['source_ra']) for p, _, _ in rows], [float(p['source_dec']) for p, _, _ in rows],
                                  c['ra_deg'], c['dec_deg'])
            scale = 3.0/float(c['radius_arcsec'])
            radius_error = max(radius_error, float(np.max(np.abs(np.hypot(x, y)*scale-np.array([r for _, r, _ in rows])))))
            entry.update(xy=np.c_[x, y], bins=np.array([b for _, _, b in rows]), weight=np.array([float(p['weight']) for p, _, _ in rows]),
                         missing=np.array([o not in available for o in ids]),
                         z=np.array([point_z.get(o, math.nan) for o in ids]),
                         fwhm_r=np.array([fwhm[o][BANDS.index('r')] for o in ids]),
                         fwhm_max=np.array([np.nanmax(fwhm[o]) if np.isfinite(fwhm[o]).any() else math.nan for o in ids]))
            entry['tree'] = cKDTree(entry['xy'])
        gx, gy = grid(float(c['radius_arcsec']))
        entry.update(grid=np.c_[gx, gy], grid_bins=annulus(np.hypot(gx, gy)*3.0/float(c['radius_arcsec'])))
        data[c['center_id']] = entry
    if radius_error > 5e-3:
        raise ValueError(f'tangent-plane radius differs from pair geometry by {radius_error} Mpc')
    nbins = len(RADIAL_EDGES)-1
    lenses = [c for c in centers if c['role'] == 'lens']
    randoms = defaultdict(list)
    for c in centers:
        if c['role'] == 'random':
            randoms[c['matched_lens_id']].append(c['center_id'])

    def area_labels(entry, footprint, extra_mask=None):
        """Per-annulus footprint and unmasked areas (arcmin^2), plus counts and weights before/after the mask."""
        out = dict(a_fp=np.zeros(nbins), a_unm=np.zeros(nbins), n_pre=np.zeros(nbins), n_post=np.zeros(nbins),
                   w_pre=np.zeros(nbins), w_post=np.zeros(nbins))
        if not entry['n']:
            return out
        dist, index = entry['tree'].query(entry['grid'])
        inside = dist <= footprint
        masked = entry['missing'][index]
        source_masked = entry['missing'].copy()
        if extra_mask is not None:
            masked = masked | extra_mask['grid']
            source_masked = source_masked | extra_mask['sources']
        cell = GRID_ARCSEC**2/3600
        out['a_fp'] = np.bincount(entry['grid_bins'][inside], minlength=nbins)*cell
        out['a_unm'] = np.bincount(entry['grid_bins'][inside & ~masked], minlength=nbins)*cell
        out['n_pre'] = np.bincount(entry['bins'], minlength=nbins).astype(float)
        out['n_post'] = np.bincount(entry['bins'][~source_masked], minlength=nbins).astype(float)
        out['w_pre'] = np.bincount(entry['bins'], weights=entry['weight'], minlength=nbins)
        out['w_post'] = np.bincount(entry['bins'][~source_masked], weights=entry['weight'][~source_masked], minlength=nbins)
        return out

    def boosts(lens_part, rand_part, k):
        row = {}
        for kind in ('n', 'w'):
            pre = ratio(lens_part[kind+'_pre'][k], lens_part['a_fp'][k], rand_part[kind+'_pre'][k], rand_part['a_fp'][k])
            naive = ratio(lens_part[kind+'_post'][k], lens_part['a_fp'][k], rand_part[kind+'_post'][k], rand_part['a_fp'][k])
            area = ratio(lens_part[kind+'_post'][k], lens_part['a_unm'][k], rand_part[kind+'_post'][k], rand_part['a_unm'][k])
            suffix = '' if kind == 'n' else '_weighted'
            row.update({'boost_pre'+suffix: pre, 'boost_naive'+suffix: naive, 'boost_area'+suffix: area,
                        'naive_minus_pre'+suffix: naive/pre-1 if pre and naive else None,
                        'area_minus_pre'+suffix: area/pre-1 if pre and area else None})
        nl, nr = lens_part['n_post'][k], rand_part['n_post'][k]
        row['poisson_relative_se'] = math.sqrt(1/nl+1/nr) if nl > 0 and nr > 0 else None
        return row

    def add(parts):
        total = {key: np.zeros(nbins) for key in ('a_fp', 'a_unm', 'n_pre', 'n_post', 'w_pre', 'w_post')}
        for part in parts:
            for key in total:
                total[key] = total[key]+part[key]
        return total

    # 1. Area-normalised boosts per lens group and annulus, and the stacked boost with lens jackknife.
    boost_rows, stacked = [], {}
    parts_by_fp = {}
    for footprint in FOOTPRINTS_ARCSEC:
        parts = {}
        for lens in lenses:
            lid = lens['matched_lens_id']
            lp = area_labels(data[lens['center_id']], footprint)
            rp = add([area_labels(data[r], footprint) for r in randoms[lid]])
            parts[lid] = (lp, rp)
            for k in range(nbins):
                boost_rows.append(dict(footprint_arcsec=footprint, lens_id=lid, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                                       lens_n_pre=int(lp['n_pre'][k]), lens_n_post=int(lp['n_post'][k]),
                                       lens_area_fp=lp['a_fp'][k], lens_area_unmasked=lp['a_unm'][k],
                                       random_n_pre=int(rp['n_pre'][k]), random_n_post=int(rp['n_post'][k]),
                                       random_area_fp=rp['a_fp'][k], random_area_unmasked=rp['a_unm'][k],
                                       lens_masked_count_fraction=1-lp['n_post'][k]/lp['n_pre'][k] if lp['n_pre'][k] else None,
                                       lens_masked_area_fraction=1-lp['a_unm'][k]/lp['a_fp'][k] if lp['a_fp'][k] else None,
                                       random_masked_count_fraction=1-rp['n_post'][k]/rp['n_pre'][k] if rp['n_pre'][k] else None,
                                       random_masked_area_fraction=1-rp['a_unm'][k]/rp['a_fp'][k] if rp['a_fp'][k] else None,
                                       **boosts(lp, rp, k)))
        parts_by_fp[footprint] = parts
        supported = [lid for lid, (lp, rp) in parts.items() if rp['a_fp'].sum() > 0]

        def matched(subset, k):
            # Pool a lens group in annulus k only where both its lens and its randoms have footprint there.
            use = [l for l in subset if parts[l][0]['a_fp'][k] > 0 and parts[l][1]['a_fp'][k] > 0]
            return add([parts[l][0] for l in use]), add([parts[l][1] for l in use]), len(use)
        for subset_name, subset in (('all_supported', supported), ('without_'+DONOR_LENS, [l for l in supported if l != DONOR_LENS])):
            for k in range(nbins):
                lp, rp, n_groups = matched(subset, k)
                row = dict(footprint_arcsec=footprint, lens_id='STACK:'+subset_name, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                           stacked_groups=n_groups,
                           lens_n_pre=int(lp['n_pre'][k]), lens_n_post=int(lp['n_post'][k]), lens_area_fp=lp['a_fp'][k],
                           lens_area_unmasked=lp['a_unm'][k], random_n_pre=int(rp['n_pre'][k]), random_n_post=int(rp['n_post'][k]),
                           random_area_fp=rp['a_fp'][k], random_area_unmasked=rp['a_unm'][k],
                           lens_masked_count_fraction=1-lp['n_post'][k]/lp['n_pre'][k] if lp['n_pre'][k] else None,
                           lens_masked_area_fraction=1-lp['a_unm'][k]/lp['a_fp'][k] if lp['a_fp'][k] else None,
                           random_masked_count_fraction=1-rp['n_post'][k]/rp['n_pre'][k] if rp['n_pre'][k] else None,
                           random_masked_area_fraction=1-rp['a_unm'][k]/rp['a_fp'][k] if rp['a_fp'][k] else None,
                           **boosts(lp, rp, k))
                if subset_name == 'all_supported':
                    for key in ('naive_minus_pre', 'area_minus_pre'):
                        replicas = []
                        for drop in subset:
                            rest_lens, rest_random, _ = matched([l for l in subset if l != drop], k)
                            replicas.append(boosts(rest_lens, rest_random, k)[key])
                        row[key+'_jackknife_se'] = jackknife(replicas)
                boost_rows.append(row)
        stacked[footprint] = supported

    # 4b. Mask transplant onto other lenses: boost null (true answer is the recipient's own boost).
    donors = {'lens:'+DONOR_LENS: next(c['center_id'] for c in lenses if c['matched_lens_id'] == DONOR_LENS)}
    random_fraction = {cid: data[cid]['missing'].mean() for cid in data if data[cid]['n'] and data[cid]['center']['role'] == 'random'}
    donors['random:most_masked'] = max(random_fraction, key=random_fraction.get)

    def transplant(donor_id, entry, footprint):
        donor = data[donor_id]
        scale = float(donor['center']['radius_arcsec'])/float(entry['center']['radius_arcsec'])
        result = {}
        for key, xy in (('grid', entry['grid']), ('sources', entry.get('xy', np.zeros((0, 2))))):
            dist, index = donor['tree'].query(xy*scale)
            result[key] = (dist <= footprint) & donor['missing'][np.minimum(index, donor['n']-1)]
            result[key+'_covered'] = float((dist <= footprint).mean()) if len(xy) else None
        return result

    transplant_boost, transplant_m = [], []
    footprint = PRIMARY_FOOTPRINT

    def transplant_row(dname, recipient, k, lp0, lp1, rp, coverage):
        before, after = boosts(lp0, rp, k), boosts(lp1, rp, k)
        return dict(donor=dname, recipient=recipient, annulus=k, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1],
                    transplanted_count_fraction=1-lp1['n_post'][k]/lp0['n_post'][k] if lp0['n_post'][k] else None,
                    donor_grid_coverage=coverage,
                    area_after_over_area_before=after['boost_area']/before['boost_area'] if after['boost_area'] and before['boost_area'] else None,
                    naive_after_over_naive_before=after['boost_naive']/before['boost_naive'] if after['boost_naive'] and before['boost_naive'] else None,
                    poisson_relative_se=after['poisson_relative_se'])
    for dname, did in donors.items():
        pooled = []
        for lens in lenses:
            lid = lens['matched_lens_id']
            # Recipients are the lightly masked lenses; the donor lens carries its own real mask.
            if lens['center_id'] == did or lid == DONOR_LENS or not data[lens['center_id']]['n']:
                continue
            entry = data[lens['center_id']]
            extra = transplant(did, entry, footprint)
            lp0, rp = parts_by_fp[footprint][lid]
            lp1 = area_labels(entry, footprint, extra)
            for k in range(nbins):
                transplant_boost.append(transplant_row(dname, lid, k, lp0, lp1, rp, extra['grid_covered']))
            if rp['a_fp'].sum() > 0:
                pooled.append((lp0, lp1, rp))
            transplant_m.append(dict(donor=dname, recipient=lid, extra=extra))
        for k in range(nbins):
            use = [p for p in pooled if p[0]['a_fp'][k] > 0 and p[2]['a_fp'][k] > 0]
            lp0, lp1, rp = (add([p[i] for p in use]) for i in range(3))
            row = transplant_row(dname, f'STACK:{len(pooled)}_recipients_with_randoms', k, lp0, lp1, rp, None)
            row['stacked_groups'] = len(use)
            transplant_boost.append(row)

    # 3/4. Photo-z multiplier: all sources, near-threshold cuts, lens jackknife and transplant nulls.
    lens_m = {}
    for lens in lenses:
        lid, zl = lens['matched_lens_id'], float(lens['z_lens'])
        entry = data[lens['center_id']]
        means, mass = means_for(zl)
        z = entry['z']
        use = np.isfinite(z) & (z > zl+GAP)
        chi = float(cosmo.comoving(zl))
        eff = np.zeros(len(z))
        eff[use] = entry['weight'][use]*(1-chi/cosmo.comoving(z[use]))**2
        index = np.zeros(len(z), int)
        index[use] = np.searchsorted(edges, z[use], side='right')-1
        if np.any(mass[index[use]] <= 0):
            raise ValueError('observed redshift without COSMOS support')
        mu = np.where(use, means[index], 0.0)
        lens_m[lid] = dict(use=use, eff=eff, mu=mu, fwhm_r=entry['fwhm_r'], fwhm_max=entry['fwhm_max'])
        m_all = eff.sum()/(eff*mu).sum()
        if abs(m_all/baseline[lid]-1) > 1e-9:
            raise ValueError(f'all-radius multiplier differs from the recorded calibration for {lid}')

    def m_of(lid, keep, counts=None):
        d = lens_m[lid]
        w = d['eff']*keep if counts is None else d['eff']*keep*counts
        return w.sum(), (w*d['mu']).sum()

    cuts = {'r_fwhm_gt_1p0': lambda d: ~(d['fwhm_r'] > NEAR_FWHM_ARCSEC), 'any_band_fwhm_gt_1p0': lambda d: ~(d['fwhm_max'] > NEAR_FWHM_ARCSEC)}
    m_rows, boot = [], {}
    for lid, d in lens_m.items():
        counts = rng.poisson(1.0, size=(BOOTSTRAP, len(d['eff'])))
        s, t = m_of(lid, np.ones(len(d['eff']), bool))
        r_fwhm = d['fwhm_r'][d['use']]
        row = dict(lens_id=lid, sources_with_pointz=int(d['use'].sum()), multiplier_all=s/t,
                   r_fwhm_median=float(np.nanmedian(r_fwhm)), r_fwhm_q05=float(np.nanpercentile(r_fwhm, 5)),
                   r_fwhm_q95=float(np.nanpercentile(r_fwhm, 95)))
        boot[lid] = dict(all=((counts*d['eff']).sum(1), (counts*d['eff']*d['mu']).sum(1)))
        for name, rule in cuts.items():
            keep = rule(d)
            s1, t1 = m_of(lid, keep)
            retained = float((d['eff']*keep).sum()/d['eff'].sum())
            bs, bt = (counts*d['eff']*keep).sum(1), (counts*d['eff']*keep*d['mu']).sum(1)
            boot[lid][name] = (bs, bt)
            supported = retained >= MIN_RETAINED_WEIGHT
            row[name+'_removed_weight_fraction'] = 1-retained
            row[name+'_per_lens_supported'] = supported
            row[name+'_multiplier'] = s1/t1 if supported else None
            row[name+'_change'] = (s1/t1)/(s/t)-1 if supported else None
            row[name+'_change_bootstrap_se'] = float(np.std((bs/bt)/(boot[lid]['all'][0]/boot[lid]['all'][1])-1, ddof=1)) if supported else None
        m_rows.append(row)

    def stack(lids, name, replica=None):
        if replica is None:
            parts = [m_of(l, np.ones(len(lens_m[l]['eff']), bool) if name == 'all' else cuts[name](lens_m[l])) for l in lids]
            return sum(p[0] for p in parts)/sum(p[1] for p in parts)
        return sum(boot[l][name][0] for l in lids)/sum(boot[l][name][1] for l in lids)
    all_lids = list(lens_m)
    for subset_name, lids in (('all', all_lids), ('without_'+DONOR_LENS, [l for l in all_lids if l != DONOR_LENS])):
        row = dict(lens_id='STACK:'+subset_name, sources_with_pointz=int(sum(lens_m[l]['use'].sum() for l in lids)), multiplier_all=stack(lids, 'all'))
        row['multiplier_all_jackknife_se'] = jackknife([stack([l for l in lids if l != drop], 'all') for drop in lids])
        for name in cuts:
            change = stack(lids, name)/stack(lids, 'all')-1
            row[name+'_multiplier'] = stack(lids, name)
            row[name+'_change'] = change
            row[name+'_change_bootstrap_se'] = float(np.std(stack(lids, name, True)/stack(lids, 'all', True)-1, ddof=1))
            row[name+'_change_jackknife_se'] = jackknife([stack([l for l in lids if l != drop], name)/stack([l for l in lids if l != drop], 'all')-1 for drop in lids])
        m_rows.append(row)
    jackknife_rows = [dict(dropped_lens=drop, stacked_multiplier_all=stack([l for l in all_lids if l != drop], 'all'),
                           **{name+'_change': stack([l for l in all_lids if l != drop], name)/stack([l for l in all_lids if l != drop], 'all')-1 for name in cuts})
                      for drop in all_lids]

    transplant_rows = []
    for item in transplant_m:
        lid = item['recipient']
        d = lens_m[lid]
        keep = ~item['extra']['sources']
        s, t = m_of(lid, np.ones(len(d['eff']), bool))
        s1, t1 = m_of(lid, keep)
        removed = int((d['use'] & ~keep).sum())
        thin = []
        candidates = np.flatnonzero(d['use'])
        for _ in range(THINNING):
            drop = np.zeros(len(d['eff']), bool)
            drop[rng.choice(candidates, size=removed, replace=False)] = True
            s2, t2 = m_of(lid, ~drop)
            thin.append((s2/t2)/(s/t)-1)
        change = (s1/t1)/(s/t)-1
        sd = float(np.std(thin, ddof=1)) if removed else None
        transplant_rows.append(dict(donor=item['donor'], recipient=lid, donor_source_coverage=item['extra']['sources_covered'],
                                    removed_pointz_sources=removed,
                                    removed_weight_fraction=float(1-(d['eff']*keep).sum()/d['eff'].sum()),
                                    multiplier_change=change, random_thinning_sd=sd,
                                    change_over_thinning_sd=change/sd if sd else None))

    # 2/3. Same-cell proxy tests on unique selected objects.
    pos = {}
    for p in pairs:
        if p['object_id'] in selected:
            pos.setdefault(p['object_id'], (float(p['source_ra']), float(p['source_dec'])))
    ids = sorted(pos, key=int)
    cell = np.array([int(math.floor(pos[o][0]/CELL_DEG))*100000+int(math.floor(pos[o][1]/CELL_DEG))+50000 for o in ids], np.int64)
    missing = np.array([o not in available for o in ids])
    mag = {b: np.array([-2.5*math.log10(f) if (f := num(props[o][b+'_flux'])) > 0 else math.nan for o in ids]) for b in BANDS}
    imag = np.array([num(props[o]['i_mag'])-num(props[o]['a_i']) for o in ids])
    imag_bin = np.where(np.isfinite(imag), np.floor(imag/0.5), -1).astype(np.int64)
    sizes = defaultdict(int)
    frac = defaultdict(int)
    for c, m in zip(cell, missing):
        sizes[c] += 1
        frac[c] += m
    mixed = np.array([sizes[c] >= MIN_CELL and 0.1 <= frac[c]/sizes[c] <= 0.9 for c in cell])
    proxy_rows = []
    colours = {'g-i': mag['g']-mag['i'], 'i-z': mag['i']-mag['z'], 'z-y': mag['z']-mag['y'], 'r-i': mag['r']-mag['i']}
    for name, colour in colours.items():
        for strata, group in (('cell', cell), ('cell_x_imag0.5', cell*1000+imag_bin)):
            result = stratified_test(colour[mixed], group[mixed], missing[mixed], rng)
            proxy_rows.append(dict(test='colour_'+name, groups='photoz_missing_vs_finite', scope='mixed_cells', strata=strata, **result))
    zb = np.array([point_z.get(o, math.nan) for o in ids])
    near = np.array([fwhm[o][BANDS.index('r')] > NEAR_FWHM_ARCSEC for o in ids])
    finite = ~missing
    for strata, group in (('cell', cell), ('cell_x_imag0.5', cell*1000+imag_bin)):
        result = stratified_test(zb[finite], group[finite], near[finite], rng)
        proxy_rows.append(dict(test='z_best', groups='r_fwhm_gt_1p0_vs_le_1p0', scope='finite_photoz_cells_with_both', strata=strata, **result))
        for name, colour in colours.items():
            result = stratified_test(colour[finite], group[finite], near[finite], rng)
            proxy_rows.append(dict(test='colour_'+name, groups='r_fwhm_gt_1p0_vs_le_1p0', scope='finite_photoz_cells_with_both', strata=strata, **result))

    for family in {r['groups'] for r in proxy_rows}:
        rows = sorted((r for r in proxy_rows if r['groups'] == family), key=lambda r: r['permutation_p'])
        running = 0.0
        for rank, r in enumerate(rows):
            running = max(running, min(1.0, (len(rows)-rank)*r['permutation_p']))
            r['holm_p_within_family'] = running
            r['family_size'] = len(rows)

    output.mkdir(parents=True)
    tables = {'annulus_boost.tsv': boost_rows, 'transplant_boost.tsv': transplant_boost, 'multiplier_sensitivity.tsv': m_rows,
              'multiplier_lens_jackknife.tsv': jackknife_rows, 'transplant_multiplier.tsv': transplant_rows, 'proxy_tests.tsv': proxy_rows}
    for name, rows in tables.items():
        fields = list(dict.fromkeys(k for r in rows for k in r))
        write_table(output/name, rows, fields)
    report = dict(
        science_measurement=False, calibration_applied=False, selection_changed=False,
        constants=dict(pixel_arcsec=PIXEL_ARCSEC, deconv_fwhm_arcsec_inferred=6.5*PIXEL_ARCSEC, near_threshold_fwhm_arcsec=NEAR_FWHM_ARCSEC,
                       gap=GAP, grid_arcsec=GRID_ARCSEC, footprints_arcsec=FOOTPRINTS_ARCSEC, primary_footprint_arcsec=PRIMARY_FOOTPRINT,
                       bootstrap=BOOTSTRAP, thinning=THINNING, permutations=PERMUTATIONS, seed=SEED, winsor_percentiles=WINSOR,
                       radial_edges_mpc=RADIAL_EDGES.tolist(), cell_deg=CELL_DEG, min_cell=MIN_CELL),
        donors=donors, lenses_with_random_support={str(k): v for k, v in stacked.items()},
        max_tangent_radius_error_mpc=radius_error, mixed_cell_objects=int(mixed.sum()),
        definitions=dict(
            mask='Selected source lacks Ephor_AB best/risk; identical to any grizy afterburner parent_flux_convolved_2_1 non-finite.',
            area='Grid points within the footprint radius of the nearest selected source of the same centre; a point is masked if that source is masked.',
            boost='(lens count / lens area) / (pooled random count / pooled random area) per annulus; pre uses all selected sources and footprint area, '
                  'naive uses unmasked sources over footprint area, area uses unmasked sources over unmasked area. Stacks pool a lens '
                  'group in an annulus only where both its lens and randoms have footprint (matched pooling, after independent review).',
            multiplier='1 / sum(w eff^2 mu_bin) / sum(w eff^2) with COSMOS bin means at gap 0.2, reproducing the recorded per-lens reweighted value.',
            fwhm='afterburner parent_flux_convolved_seeing x 2.3548 x 0.168 arcsec; unit inferred from the deconv flag boundary.',
            transplant='Recipient positions scaled to the donor frame in proper Mpc (same orientation); masked where the donor nearest source is masked.',
            proxy_test='Stratum-size-weighted mean difference of winsorised values with a within-stratum permutation null; '
                       'Holm adjustment within each comparison family.',
            near_threshold=f'Per-lens change reported only if at least {MIN_RETAINED_WEIGHT} of the weight is retained (rule added after the '
                           'trial run). The stacked change also includes lenses that lose most sources, so it mixes within-field '
                           'removal with a change of lens composition.'),
        limitations=['Grid areas use nearest-source labels; mask structure finer than the source spacing is not resolved.',
                     'The 4-arcsec grid under-resolves the 0-0.2 Mpc discs (area errors up to about 7% per lens); an independent '
                     '1-arcsec check moved the stacked area-minus-pre residual there by about 0.005-0.007.',
                     'Poisson boost SE is the error of the boost, not of the post-minus-pre difference, and field-to-field '
                     'lens-random density mismatch (per-lens pre-mask boosts 0.87-1.21 at 2-3 Mpc) is not included.',
                     'The transplant masks only the lens side geometrically; random-side masking and nearest-source area '
                     'inference in the real field are not tested by it.',
                     'Poisson boost errors ignore clustering and source reuse across centres; jackknife over 8-11 lenses is noisy.',
                     'The multiplier uses fixed COSMOS bin means; COSMOS sampling and reference-z errors are not propagated.',
                     'Colour and z_best tests are proxies for redshift independence, not a redshift measurement of masked sources.'],
        provenance=dict(**provenance, **ab_provenance, **cal_provenance, geometry_sha256=GEOMETRY_SHA),
        runner_sha256=sha(Path(__file__)), output_sha256={name: sha(output/name) for name in tables})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(radius_error=radius_error, donors=donors, mixed=int(mixed.sum())), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
