#!/usr/bin/env python3
"""Follow-ups to the S16A cmodel photo-z audit, offline.

(b) Seeing transfer of the COSMOS multiplier calibration. The S17A Wide-COSMOS calibration objects
    have r PSF FWHM (sdss PSF-moment proxy) of about 0.7-0.94 arcsec, while many lens sources sit at
    1.0-1.092 arcsec (afterburner seeing). Public CAS has no COSMOS wide-depth photo-z at other seeing,
    so the test stays inside COSMOS: split the calibration objects at the weighted median proxy FWHM,
    rebuild the bin means from each half, apply both to every lens's own sources and express the
    multiplier difference per arcsec (COSMOS spatial-block bootstrap). The linear extrapolation to the
    lens's own FWHM is illustrative only; the proxy and afterburner seeing scales are not cross-checked.
(c) Field-level checks of the Ephor_AB mask: same-cell missing-vs-finite tests per tract, and per lens
    the eff^2-weighted mean code z of missing minus available sources with a within-cell permutation null.

All constants are fixed before the run. Descriptive only: no selection change or correction.
Added after independent review (no earlier number changes except the per-lens effect columns):
two-sided COSMOS support fractions and support-grouped stacks; position, other-band and random-block
placebo splits (separate random stream); a direct within-stratum weighted difference with its own
permutation null as the lens test; and lens comparisons without the z cut, because the used-subset
statistic conditions on the code z it compares.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

from astropy.io import fits
import numpy as np

from audit_s16a_afterburner_mask_effects import (BOOTSTRAP, FWHM_PER_SIGMA, GAP, NEAR_FWHM_ARCSEC, PERMUTATIONS, PIXEL_ARCSEC,
                                                 calibration, load_afterburner, load_supplement, num, stratified_test)
from audit_s16a_cmodel_photoz_missingness import CODES, holm, load_cmodel
from audit_s16a_cosmos_lens_response import INTAKE, OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from audit_s16a_missingness_geometry import CELL_DEG, MIN_CELL
from prepare_s16a_afterburner_photoz_rows import BANDS
from ogfkit.wlcluster import Cosmology

COMPANION = PLAN.parent/'s16a_cosmos_companion_plan_2026-11-24_01'
COMPANION_PLAN_SHA = '29e3cb624a420b44244ce793674b1561fe1c8041b2d3cce01ed64f5fb4fa8788'
SEED = 2026112403
PLACEBO_SEED = 2026112404
PLACEBO_DRAWS = 500
BLOCK_DEG = 0.2
MIN_MISSING_USED = 20


def support_group(fwhm, q):
    return 'above_cosmos_q99' if fwhm > q['q99'] else 'below_cosmos_q01' if fwhm < q['q01'] else 'inside_cosmos_q01_q99'


def load_cosmos():
    """Eligible COSMOS calibration rows (same filter and order as calibration()) with r PSF-moment FWHM."""
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    objid = np.array([int(r['s17a_objid']) for r in rows], np.int64)
    photo, ref, shape, som = np.asarray([[float(r[k]) for k in ('photo_z', 'reference_z', 'shape_weight', 'som_weight')]
                                         for r in rows]).T
    receipt = json.loads((COMPANION/'results'/'receipt.json').read_text())
    if receipt['plan_sha256'] != COMPANION_PLAN_SHA:
        raise ValueError('companion receipt refers to another plan')
    found = {}
    for name, entry in receipt['files'].items():
        path = COMPANION/'results'/name
        if entry['status'] != 'complete' or sha(path) != entry['sha256']:
            raise ValueError('companion file changed: '+name)
        with fits.open(path, memmap=True) as hdul:
            d = hdul[1].data
            oid = np.asarray(d['object_id'], np.int64)
            keep = np.isin(oid, objid)
            bands = [FWHM_PER_SIGMA*np.sqrt((np.asarray(d[b+'shape_sdss_psf_11'][keep], float)
                                             + np.asarray(d[b+'shape_sdss_psf_22'][keep], float))/2) for b in 'rgizy']
            for k, (o, ra, dec, t) in enumerate(zip(oid[keep], np.asarray(d['ra'][keep], float), np.asarray(d['dec'][keep], float),
                                                     np.asarray(d['tract'][keep]))):
                if o in found:
                    raise ValueError('duplicate companion object')
                found[int(o)] = (*(float(f[k]) for f in bands), ra, dec, int(t))
    if set(found) != set(objid.tolist()):
        raise ValueError('calibration objects missing from the companion files')
    fwhm, fg, fi, fz, fy, ra, dec, tract = np.array([found[o] for o in objid.tolist()]).T
    if not np.all(np.isfinite(fwhm)):
        raise ValueError('non-finite COSMOS PSF moments')
    block = np.unique(np.c_[tract, np.floor(ra/BLOCK_DEG), np.floor(dec/BLOCK_DEG)], axis=0, return_inverse=True)[1].ravel()
    return dict(photo=photo, ref=ref, base=shape*som, fwhm=fwhm, block=block, ra=ra, dec=dec,
                other_bands=dict(g=fg, i=fi, z=fz, y=fy)), dict(companion_receipt_sha256=sha(COMPANION/'results'/'receipt.json'))


def weighted_quantile(x, w, q):
    order = np.argsort(x)
    c = np.cumsum(w[order])
    return float(np.interp(q*c[-1], c, x[order]))


def cosmos_terms(cosmo, zl, cos, edges):
    chi = float(cosmo.comoving(zl))
    keep = cos['photo'] > zl+GAP
    p, r = cos['photo'][keep], cos['ref'][keep]
    est = 1-chi/cosmo.comoving(p)
    true = np.zeros_like(r)
    true[r > zl] = 1-chi/cosmo.comoving(r[r > zl])
    return dict(keep=keep, weight=cos['base'][keep]*est**2, ratio=true/est, fwhm=cos['fwhm'][keep],
                index=np.searchsorted(edges, p, side='right')-1, nbins=len(edges)-1)


def bin_stats(terms, w):
    mass = np.bincount(terms['index'], weights=w, minlength=terms['nbins'])
    safe = np.where(mass > 0, mass, 1.0)
    mu = np.where(mass > 0, np.bincount(terms['index'], weights=w*terms['ratio'], minlength=terms['nbins'])/safe, 0.0)
    f = np.where(mass > 0, np.bincount(terms['index'], weights=w*terms['fwhm'], minlength=terms['nbins'])/safe, 0.0)
    return mu, f, mass


def split_response(terms, lens, mult, high, pooled=None):
    """Multiplier with each half's bin means; bins empty in either half use the pooled value (no contribution).
    In a bootstrap replicate a lens bin can lose all COSMOS objects; it then takes the full-sample pooled values."""
    w = terms['weight']*mult
    mu, f, mass = bin_stats(terms, w)
    e = lens['bin_eff']
    empty = (mass <= 0) & (e > 0)
    if empty.any():
        if pooled is None:
            raise ValueError('lens source bin without COSMOS support')
        mu, f = np.where(empty, pooled[0], mu), np.where(empty, pooled[1], f)
    mu_lo, f_lo, m_lo = bin_stats(terms, w*~high)
    mu_hi, f_hi, m_hi = bin_stats(terms, w*high)
    both = (m_lo > 0) & (m_hi > 0)
    mu_lo, mu_hi = np.where(both, mu_lo, mu), np.where(both, mu_hi, mu)
    f_lo, f_hi = np.where(both, f_lo, f), np.where(both, f_hi, f)
    s = e.sum()
    out = dict(S=s, T_all=(e*mu).sum(), T_lo=(e*mu_lo).sum(), T_hi=(e*mu_hi).sum(), F_all=(e*f).sum()/s,
               F_lo=(e*f_lo).sum()/s, F_hi=(e*f_hi).sum()/s, fallback_eff_fraction=float(e[~both].sum()/s),
               empty_pooled_bins=bool(empty.any()))
    out['m_all'], out['m_lo'], out['m_hi'] = s/out['T_all'], s/out['T_lo'], s/out['T_hi']
    out['change'] = out['m_hi']/out['m_lo']-1
    out['slope'] = out['change']/(out['F_hi']-out['F_lo'])
    out['extrapolated'] = out['slope']*(lens['fwhm_eff']-out['F_all'])
    return out


def weighted_label_test(values, weights, strata, label, rng):
    """Weighted difference of label True minus False.

    Primary: the direct within-stratum difference, i.e. per-stratum weighted mean difference averaged over
    strata that hold both labels with the label weight as stratum weight; its within-stratum permutation
    null is centred on zero by construction. Secondary: the pooled difference, whose permutation null is
    centred on a between-stratum composition term; observed minus that centre is not a within-stratum
    effect (it responds to a uniform within-stratum shift by less than one), so it is kept only for reference.
    """
    _, s = np.unique(strata, return_inverse=True)
    n = s.max()+1

    def pooled(lab, w):
        return float((w*values*lab).sum()/(w*lab).sum()-(w*values*~lab).sum()/(w*~lab).sum())

    def direct(lab, w):
        wl, wo = np.bincount(s, w*lab, n), np.bincount(s, w*~lab, n)
        ok = (wl > 0) & (wo > 0)
        if not ok.any():
            return math.nan
        diff = np.bincount(s, w*values*lab, n)[ok]/wl[ok]-np.bincount(s, w*values*~lab, n)[ok]/wo[ok]
        return float((wl[ok]*diff).sum()/wl[ok].sum())
    observed, within = pooled(label, weights), direct(label, weights)
    order = np.argsort(strata, kind='stable')
    null, null_direct = np.empty(PERMUTATIONS), np.empty(PERMUTATIONS)
    for k in range(PERMUTATIONS):
        shuffled = np.empty_like(label)
        shuffled[order] = label[np.lexsort((rng.random(len(label)), strata))]
        null[k], null_direct[k] = pooled(shuffled, weights), direct(shuffled, weights)
    boot, boot_direct = [], []
    for _ in range(BOOTSTRAP):
        c = rng.poisson(1.0, len(values))
        lab_w, other_w = weights*c*label, weights*c*~label
        if lab_w.sum() > 0 and other_w.sum() > 0:
            boot.append((lab_w*values).sum()/lab_w.sum()-(other_w*values).sum()/other_w.sum())
            boot_direct.append(direct(label, weights*c))
    centre = float(null.mean())
    sizes, with_label = Counter(strata.tolist()), Counter(strata[label].tolist())
    return dict(within_strata_difference=within, within_strata_null_sd=float(null_direct.std(ddof=1)),
                within_strata_bootstrap_se=float(np.nanstd(boot_direct, ddof=1)),
                permutation_p=float((1+(np.abs(null_direct) >= abs(within)).sum())/(1+PERMUTATIONS)),
                pooled_difference=observed, pooled_null_mean_composition=centre, pooled_minus_null_mean=observed-centre,
                pooled_null_sd=float(null.std(ddof=1)), pooled_bootstrap_se=float(np.std(boot, ddof=1)),
                pooled_centred_permutation_p=float((1+(np.abs(null-centre) >= abs(observed-centre)).sum())/(1+PERMUTATIONS)),
                strata_with_both=sum(1 for s, n in with_label.items() if n < sizes[s]),
                label_in_strata_with_both=sum(n for s, n in with_label.items() if n < sizes[s]))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    rng = np.random.default_rng(SEED)
    cosmo = Cosmology()
    means_for, edges, baseline, cal_provenance = calibration(cosmo)
    cos, cos_provenance = load_cosmos()
    cm, cm_provenance = load_cmodel()
    seeing, unavailable, ab_provenance = load_afterburner()
    props = load_supplement()
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    point_z = finite_point_redshifts(available)
    if {o for o in selected if o not in available} != unavailable & selected:
        raise ValueError('afterburner mask differs from the photo-z availability ledger')
    fwhm_r = {o: s[BANDS.index('r')]*FWHM_PER_SIGMA*PIXEL_ARCSEC for o, s in seeing.items()}

    # Unique selected objects: cells, tracts, i-mag and code redshifts (as in the cmodel audit).
    pos = {}
    for p in pairs:
        if p['object_id'] in selected:
            pos.setdefault(p['object_id'], (float(p['source_ra']), float(p['source_dec'])))
    ids = sorted(pos, key=int)
    index_of = {o: i for i, o in enumerate(ids)}
    missing = np.array([o not in available for o in ids])
    ra, dec = np.array([pos[o] for o in ids]).T
    cell = np.floor(ra/CELL_DEG).astype(np.int64)*100000+np.floor(dec/CELL_DEG).astype(np.int64)+50000
    tract = np.array([int(o) >> 42 for o in ids])
    imag = np.array([num(props[o]['i_mag'])-num(props[o]['a_i']) for o in ids])
    imag_bin = np.where(np.isfinite(imag), np.floor(imag/0.5), -1).astype(np.int64)
    z = {code: np.array([cm[o][col] for o in ids]) for code, col in CODES.items()}
    z_ab = np.array([point_z.get(o, math.nan) for o in ids])
    sizes, frac = Counter(cell.tolist()), Counter(cell[missing].tolist())
    mixed = np.array([sizes[c] >= MIN_CELL and 0.1 <= frac[c]/sizes[c] <= 0.9 for c in cell.tolist()])
    seeing_ids = np.array([seeing[o] for o in ids], float)*FWHM_PER_SIGMA*PIXEL_ARCSEC
    r_deconv = seeing_ids[:, BANDS.index('r')] > 6.5*PIXEL_ARCSEC
    y_deconv = seeing_ids[:, BANDS.index('y')] > 6.5*PIXEL_ARCSEC

    by_center = defaultdict(list)
    for p in pairs:
        if p['object_id'] in selected:
            by_center[p['center_id']].append(p)
    lenses = [c for c in centers if c['role'] == 'lens']

    # (b) Lens sources with Ephor_AB, as in the multiplier: eff = w (1 - chi_l/chi_s)^2 for z_AB > z_l + gap.
    q = {k: weighted_quantile(cos['fwhm'], cos['base'], v) for k, v in (('q01', .01), ('q50', .5), ('q95', .95), ('q99', .99))}
    high = cos['fwhm'] > q['q50']
    lens_data = {}
    for lens in lenses:
        lid, zl = lens['matched_lens_id'], float(lens['z_lens'])
        rows = by_center[lens['center_id']]
        o = [p['object_id'] for p in rows]
        zs = np.array([point_z.get(i, math.nan) for i in o])
        use = np.isfinite(zs) & (zs > zl+GAP)
        w = np.array([float(p['weight']) for p in rows])[use]
        eff = w*(1-float(cosmo.comoving(zl))/cosmo.comoving(zs[use]))**2
        idx = np.searchsorted(edges, zs[use], side='right')-1
        f = np.array([fwhm_r[i] for i in np.asarray(o)[use]])
        terms = cosmos_terms(cosmo, zl, cos, edges)
        mu_ref, _ = means_for(zl)
        mu_check, f_all, _ = bin_stats(terms, terms['weight'])
        if not np.allclose(mu_check, mu_ref, rtol=0, atol=1e-12):
            raise ValueError('COSMOS bin means differ from the recorded calibration')
        lens_data[lid] = dict(zl=zl, terms=terms, bin_eff=np.bincount(idx, weights=eff, minlength=len(edges)-1),
                              fwhm_eff=float((eff*f).sum()/eff.sum()), eff=eff, fwhm=f, n=int(use.sum()), pooled=(mu_check, f_all))
        m = eff.sum()/(eff*mu_ref[idx]).sum()
        if abs(m/baseline[lid]-1) > 1e-9:
            raise ValueError('multiplier differs from the recorded calibration for '+lid)

    ones = np.ones(len(cos['photo']))
    blocks = cos['block'].max()+1
    replicas = []
    for _ in range(BOOTSTRAP):
        counts = np.bincount(rng.integers(0, blocks, blocks), minlength=blocks)[cos['block']].astype(float)
        replicas.append({lid: split_response(d['terms'], d, counts[d['terms']['keep']], high[d['terms']['keep']], d['pooled'])
                         for lid, d in lens_data.items()})
    seeing_rows = []
    for lid, d in lens_data.items():
        point = split_response(d['terms'], d, ones[d['terms']['keep']], high[d['terms']['keep']])
        rep = {k: np.array([r[lid][k] for r in replicas]) for k in ('change', 'slope', 'extrapolated')}
        e = d['eff']
        seeing_rows.append(dict(
            lens_id=lid, z_lens=d['zl'], sources_used=d['n'], replicas_with_empty_lens_bin=sum(r[lid]['empty_pooled_bins'] for r in replicas),
            multiplier_all=point['m_all'], multiplier_low_fwhm_half=point['m_lo'],
            multiplier_high_fwhm_half=point['m_hi'], high_over_low_minus_1=point['change'],
            high_over_low_minus_1_block_bootstrap_se=float(rep['change'].std(ddof=1)),
            cosmos_fwhm_low_half=point['F_lo'], cosmos_fwhm_high_half=point['F_hi'], cosmos_fwhm_all=point['F_all'],
            change_per_arcsec=point['slope'], change_per_arcsec_block_bootstrap_se=float(rep['slope'].std(ddof=1)),
            lens_source_fwhm_eff_weighted=d['fwhm_eff'], extrapolated_change_to_lens_fwhm=point['extrapolated'],
            extrapolated_change_block_bootstrap_se=float(rep['extrapolated'].std(ddof=1)),
            extrapolated_change_q16=float(np.percentile(rep['extrapolated'], 16)),
            extrapolated_change_q84=float(np.percentile(rep['extrapolated'], 84)),
            fallback_eff_fraction=point['fallback_eff_fraction'],
            eff_fraction_r_fwhm_gt_1p0=float(e[d['fwhm'] > NEAR_FWHM_ARCSEC].sum()/e.sum()),
            eff_fraction_above_cosmos_q99=float(e[d['fwhm'] > q['q99']].sum()/e.sum()),
            eff_fraction_below_cosmos_q01=float(e[d['fwhm'] < q['q01']].sum()/e.sum()),
            support_group=support_group(d['fwhm_eff'], q)))

    def stack(responses):
        s = sum(r['S'] for r in responses.values())
        lo, hi = s/sum(r['T_lo'] for r in responses.values()), s/sum(r['T_hi'] for r in responses.values())
        df = sum(r['S']*(r['F_hi']-r['F_lo']) for r in responses.values())/s
        extrap = sum(r['S']*r['extrapolated'] for r in responses.values())/s
        return hi/lo-1, (hi/lo-1)/df, extrap
    points = {lid: split_response(d['terms'], d, ones[d['terms']['keep']], high[d['terms']['keep']]) for lid, d in lens_data.items()}
    groups = {'STACK:all': list(lens_data)}
    for name in ('above_cosmos_q99', 'inside_cosmos_q01_q99', 'below_cosmos_q01'):
        members = [lid for lid, d in lens_data.items() if support_group(d['fwhm_eff'], q) == name]
        if members:
            groups['STACK:'+name] = members
    for name, members in groups.items():
        point = stack({l: points[l] for l in members})
        rep = np.array([stack({l: r[l] for l in members}) for r in replicas])
        total_eff = sum(lens_data[l]['eff'].sum() for l in members)
        seeing_rows.append(dict(
            lens_id=name, stacked_lenses=','.join(members), high_over_low_minus_1=point[0],
            high_over_low_minus_1_block_bootstrap_se=float(rep[:, 0].std(ddof=1)),
            change_per_arcsec=point[1], change_per_arcsec_block_bootstrap_se=float(rep[:, 1].std(ddof=1)),
            extrapolated_change_to_lens_fwhm=point[2], extrapolated_change_block_bootstrap_se=float(rep[:, 2].std(ddof=1)),
            extrapolated_change_q16=float(np.percentile(rep[:, 2], 16)), extrapolated_change_q84=float(np.percentile(rep[:, 2], 84)),
            eff_fraction_r_fwhm_gt_1p0=float(sum(lens_data[l]['eff'][lens_data[l]['fwhm'] > NEAR_FWHM_ARCSEC].sum() for l in members)/total_eff),
            eff_fraction_above_cosmos_q99=float(sum(lens_data[l]['eff'][lens_data[l]['fwhm'] > q['q99']].sum() for l in members)/total_eff),
            eff_fraction_below_cosmos_q01=float(sum(lens_data[l]['eff'][lens_data[l]['fwhm'] < q['q01']].sum() for l in members)/total_eff)))

    # Placebo splits (added after review, separate random stream): position, position-controlled r, other bands,
    # random halves of whole sky blocks. If these move the stacked multiplier as much as the r split, the r split
    # does not isolate seeing.
    prng = np.random.default_rng(PLACEBO_SEED)
    block_mean = np.bincount(cos['block'], cos['base']*cos['fwhm'])/np.bincount(cos['block'], cos['base'])
    variables = {'r_fwhm': cos['fwhm'], 'ra': cos['ra'], 'dec': cos['dec'], 'r_fwhm_minus_block_mean': cos['fwhm']-block_mean[cos['block']],
                 **{b+'_fwhm': v for b, v in cos['other_bands'].items()}}

    def stacked_change(split):
        return stack({lid: split_response(d['terms'], d, ones[d['terms']['keep']], split[d['terms']['keep']])
                      for lid, d in lens_data.items()})[0]
    placebo_rows = []
    for name, v in variables.items():
        split = v > weighted_quantile(v, cos['base'], 0.5)
        placebo_rows.append(dict(split=name, stacked_high_over_low_minus_1=stacked_change(split),
                                 label_agreement_with_r_split=float(np.average(split == high, weights=cos['base'])),
                                 r_fwhm_high_minus_low=float(np.average(cos['fwhm'][split], weights=cos['base'][split])
                                                             - np.average(cos['fwhm'][~split], weights=cos['base'][~split]))))
    random_changes = np.array([stacked_change(prng.random(blocks)[cos['block']] < 0.5) for _ in range(PLACEBO_DRAWS)])
    observed_change = placebo_rows[0]['stacked_high_over_low_minus_1']
    placebo_rows.append(dict(split=f'random_block_halves_{PLACEBO_DRAWS}_draws', stacked_high_over_low_minus_1_sd=float(random_changes.std(ddof=1)),
                             fraction_abs_ge_r_split=float(np.mean(np.abs(random_changes) >= abs(observed_change)))))

    # (c1) Same-cell tests per tract (mixed cells, as the primary cmodel test).
    tract_rows = []
    for t, n in Counter(tract[mixed].tolist()).most_common():
        scope = mixed & (tract == t)
        for code in CODES:
            ok = np.isfinite(z[code][scope])
            row = dict(family='per_tract_mixed_cells', tract=t, test=code+'_z_best', objects=int(scope.sum()),
                       missing=int((scope & missing).sum()))
            if (missing[scope] & ok).sum() and (~missing[scope] & ok).sum():
                row.update(stratified_test(z[code][scope], cell[scope], missing[scope], rng), testable=True)
            else:
                row['testable'] = False
            tract_rows.append(row)
    holm(tract_rows)

    # (c2) Per lens with missing sources used: eff^2-weighted mean code z, missing minus available.
    lens_rows = []
    for code in CODES:
        for lens in lenses:
            lid, zl = lens['matched_lens_id'], float(lens['z_lens'])
            rows = by_center[lens['center_id']]
            idx = np.array([index_of[p['object_id']] for p in rows])
            zs = z[code][idx]
            use = np.isfinite(zs) & (zs > zl+GAP)
            lab = missing[idx][use]
            if lab.sum() < MIN_MISSING_USED:
                continue
            w = np.array([float(p['weight']) for p in rows])[use]
            weff = w*(1-float(cosmo.comoving(zl))/cosmo.comoving(zs[use]))**2
            zu, iu = zs[use], idx[use]
            row = dict(family='per_lens_'+code, lens_id=lid, z_lens=zl, missing_used=int(lab.sum()), available_used=int((~lab).sum()),
                       missing_r_deconv=int(r_deconv[iu][lab].sum()), missing_y_deconv=int(y_deconv[iu][lab].sum()),
                       median_z_missing=float(np.median(zu[lab])), median_z_available=float(np.median(zu[~lab])),
                       fraction_z_gt_2_missing=float(np.mean(zu[lab] > 2)), fraction_z_gt_2_available=float(np.mean(zu[~lab] > 2)),
                       median_imag_missing=float(np.nanmedian(imag[iu][lab])), median_imag_available=float(np.nanmedian(imag[iu][~lab])),
                       median_r_fwhm_missing=float(np.median([fwhm_r[ids[i]] for i in iu[lab]])),
                       median_r_fwhm_available=float(np.median([fwhm_r[ids[i]] for i in iu[~lab]])))
            for strata_name, strata in (('cell', cell[iu]), ('cell_x_imag0.5', cell[iu]*1000+imag_bin[iu])):
                result = weighted_label_test(zu, weff, strata, lab, rng)
                lens_rows.append(dict(row, test='eff2_weighted_mean_z_missing_minus_available', strata=strata_name, testable=True, **result))
    # Added after review: the used-subset statistic conditions on the code z it compares, so a wider photo-z
    # scatter among missing sources shows up as a mean difference. Compare all lens sources with finite code z,
    # weighted by w only, on z itself, the fractions below the cut and above z = 2, and the per-source lensing
    # efficiency (zero for sources at or below the cut).
    for code in CODES:
        for lens in lenses:
            lid, zl = lens['matched_lens_id'], float(lens['z_lens'])
            rows = by_center[lens['center_id']]
            idx = np.array([index_of[p['object_id']] for p in rows])
            zs = z[code][idx]
            if (missing[idx] & np.isfinite(zs) & (zs > zl+GAP)).sum() < MIN_MISSING_USED:
                continue
            ok = np.isfinite(zs)
            zs, iu, lab = zs[ok], idx[ok], missing[idx][ok]
            w = np.array([float(p['weight']) for p in rows])[ok]
            eff = np.where(zs > zl+GAP, 1-float(cosmo.comoving(zl))/cosmo.comoving(np.maximum(zs, zl+GAP)), 0.0)
            quantities = {'mean_z_all_sources': zs, 'fraction_at_or_below_cut': (zs <= zl+GAP).astype(float),
                          'fraction_z_gt_2': (zs > 2).astype(float), 'mean_efficiency_per_source': eff}
            for qname, values in quantities.items():
                for strata_name, strata in (('cell', cell[iu]), ('cell_x_imag0.5', cell[iu]*1000+imag_bin[iu])):
                    result = weighted_label_test(values, w, strata, lab, rng)
                    lens_rows.append(dict(family='no_cut_'+code, lens_id=lid, z_lens=zl, missing_used=int(lab.sum()),
                                          available_used=int((~lab).sum()), test=qname+'_missing_minus_available',
                                          strata=strata_name, testable=True, **result))
    holm(lens_rows)

    output.mkdir(parents=True)
    tables = {'cosmos_seeing_split.tsv': seeing_rows, 'cosmos_placebo_splits.tsv': placebo_rows, 'tract_tests.tsv': tract_rows,
              'lens_weighted_tests.tsv': lens_rows}
    for name, rows in tables.items():
        write_table(output/name, rows, list(dict.fromkeys(k for r in rows for k in r)))
    report = dict(
        science_measurement=False, calibration_applied=False, selection_changed=False,
        constants=dict(seed=SEED, block_deg=BLOCK_DEG, bootstrap=BOOTSTRAP, permutations=PERMUTATIONS, gap=GAP,
                       near_threshold_fwhm_arcsec=NEAR_FWHM_ARCSEC, min_missing_used=MIN_MISSING_USED, cell_deg=CELL_DEG, min_cell=MIN_CELL),
        cosmos_calibration_objects=len(cos['photo']), cosmos_blocks=int(blocks),
        cosmos_r_fwhm_proxy_weighted_quantiles=q, cosmos_split_threshold_arcsec=q['q50'],
        cosmos_r_fwhm_proxy_max=float(cos['fwhm'].max()),
        missing_in_mixed_cells=int((mixed & missing).sum()),
        public_cas_check=('No public COSMOS wide-depth seeing-variant schema (pdr1_cosmos_widedepth_best/median/worst, '
                          'pdr2_cosmos_wide_depth_*, pdr3/s19a *_dud_wide_depth_*) has a photo-z or afterburner table; no '
                          'cosmos/wide_depth schema has a photo-z table, and the only afterburner table among them is the deep '
                          'pdr1_ssp_uh_cosmos stack (schema browser schemata.json, retrieved 2026-11-24, table names only, '
                          'columns not checked). A CAS Ephor_AB recalibration at worse seeing is not possible.'),
        definitions=dict(
            proxy_fwhm='2.3548 x sqrt of the mean of rshape_sdss_psf_11 and _22 in the S17A companion files; the arcsec^2 unit '
                       'is inferred (no TUNIT; pixel^2 would give an unphysical 0.14 arcsec).',
            support='Lens groups by eff-weighted afterburner r FWHM relative to the COSMOS weighted q01/q99 of the proxy. '
                    'STACK:all nets upward and downward extrapolations; read the support-grouped stacks instead.',
            placebo='Stacked change for median splits on RA, Dec, r FWHM minus its block mean, and g/i/z/y proxy FWHM, plus '
                    'random halves of whole sky blocks (separate random stream).',
            split='COSMOS objects above vs at-or-below the shape x SOM weighted median proxy FWHM. Bin means per half with '
                  'the calibration weights; bins empty in either half use the pooled value for both halves.',
            change_per_arcsec='(m_high / m_low - 1) / (effective FWHM high - low), effective FWHM weighted by the lens source bins.',
            extrapolated='change_per_arcsec x (lens eff-weighted afterburner r FWHM - COSMOS effective FWHM); illustrative.',
            bootstrap='COSMOS objects resampled by 0.2-degree sky blocks with the split threshold fixed; lens sources fixed.',
            lens_test='Used sources (code z > z_l + 0.2) of the lens centre; eff^2 = w (1 - chi_l/chi_s)^2 from the same code. '
                      'within_strata_difference: per-cell weighted mean difference averaged with the missing weight over cells '
                      'holding both groups; its within-cell permutation null gives permutation_p (Holm within each family; '
                      'the two strata variants count as separate members, which is conservative). The pooled_* columns are '
                      'the earlier pooled statistic, kept for reference.',
            no_cut_test='All lens sources with finite code z, weight w: z, P(z <= z_l + 0.2), P(z > 2) and per-source '
                        'efficiency (0 at or below the cut), missing minus available, same estimator and null.'),
        limitations=['The proxy FWHM (sdss PSF moments, S17A) and the afterburner seeing (S16A) are not cross-calibrated; the '
                     'extrapolation assumes the same scale and a linear trend well beyond the COSMOS range.',
                     'COSMOS seeing varies across the field, so the split also splits sky position; block bootstrap covers part of that.',
                     'Lens-side sampling noise is not included in the seeing-split errors.',
                     'Afterburner photometry is PSF-matched to a common target, so seeing acts mainly through noise; this test '
                     'cannot separate that from other position-dependent COSMOS properties.',
                     'In COSMOS the r proxy FWHM tracks RA; position and other-band placebo splits show whether the r split '
                     'isolates seeing (cosmos_placebo_splits.tsv). Only one calibration object is above 1.0 arcsec, so an '
                     'onset near the 1.092-arcsec threshold is untestable.',
                     'Blocks are keyed on tract as well as sky cell, so co-located tract 9812/9813 blocks resample independently.',
                     'The used-subset lens statistic conditions on the code z it compares; read it with the no-cut rows.',
                     'Tract tests cover only missing sources in mixed cells; tracts other than 10048 have little power.'],
        provenance=dict(**provenance, **ab_provenance, **cm_provenance, **cal_provenance, **cos_provenance),
        runner_sha256=sha(Path(__file__)), output_sha256={name: sha(output/name) for name in tables})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(cosmos=len(cos['photo']), blocks=int(blocks), split=q['q50'], lenses=len(lens_data)), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
