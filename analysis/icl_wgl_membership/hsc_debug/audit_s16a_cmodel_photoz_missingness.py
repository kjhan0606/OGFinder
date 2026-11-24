#!/usr/bin/env python3
"""Redshifts of S16A sources without Ephor_AB, from photo-z codes that do not need afterburner fluxes.

Analyses fixed in the approved CAS plan (35265ee4..., analysis_plan) before the download:
1. Ephor(cmodel) and Mizuki minus Ephor_AB z_best on sources with both, by afterburner r-band FWHM.
   6575 of 6671 missing sources lie above the deconvolution threshold (r FWHM > 1.092 arcsec), outside
   the range of sources with Ephor_AB, so the near-threshold bins are reported and no calibration is applied.
2. Same-cell stratified permutation tests of Ephor(cmodel) and Mizuki z_best, missing vs finite.
3. Per lens, lensing-efficiency-weighted mean z of all selected sources vs Ephor_AB-available
   sources only, with the same photo-z code and redshift cut for both sets.
Not in the hashed plan: the binning constants, the secondary comparison (missing vs finite with r FWHM
> 1.0), the Mizuki per-lens rows and the stacks were set while the download ran, before any analysis
output was seen. After independent review: tract-jackknife slope errors, a core (non-outlier) slope and
a fixed-lens-weight stack were added without changing any other number.
Descriptive only: no selection change, correction or science promotion.
"""
import argparse
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

import numpy as np

from audit_s16a_afterburner_mask_effects import (BOOTSTRAP, FWHM_PER_SIGMA, GAP, NEAR_FWHM_ARCSEC, PERMUTATIONS, PIXEL_ARCSEC,
                                                 WINSOR, jackknife, load_afterburner, load_supplement, num, stratified_test)
from audit_s16a_cosmos_lens_response import OBJECT_QA, PLAN, finite_point_redshifts, load_inputs, sha, write_table
from audit_s16a_missingness_geometry import CELL_DEG, MIN_CELL
from fetch_ancillary_catalogues import catalogue_reader
from prepare_s16a_afterburner_photoz_rows import BANDS
from prepare_s16a_cmodel_photoz_rows import OUTPUT as CM_PLAN
from submit_s16a_cmodel_photoz_rows import load_plan as load_cm_plan, verify_complete as verify_cm
from ogfkit.wlcluster import Cosmology

SEED = 2026112402
DECONV_FWHM_ARCSEC = 6.5*PIXEL_ARCSEC
FWHM_EDGES = (0.0, 0.6, 0.7, 0.8, 0.9, 1.0, 1.05, DECONV_FWHM_ARCSEC, 1.12, 1.15, math.inf)
SLOPE_RANGE = (0.8, DECONV_FWHM_ARCSEC)
OUTLIER = 0.15
MIN_BIN = 20
CODES = {'ephor_cmodel': 'ephor_photoz_best', 'mizuki': 'mizuki_photoz_best'}
DIAGNOSTICS = ('ephor_photoz_risk_best', 'ephor_photoz_std_best', 'mizuki_photoz_risk_best', 'mizuki_photoz_std_best',
               'mizuki_reduced_chisq', 'mizuki_prob_star', 'mizuki_prob_gal', 'mizuki_prob_qso')
CHUNK = 100


def load_cmodel():
    plan = load_cm_plan(CM_PLAN)
    results = CM_PLAN/'results'
    state = json.loads((results/'status.json').read_text())
    if state['status'] != 'complete' or sha(results/'download_summary.json') != state['summary_sha256']:
        raise ValueError('cmodel photo-z receipt is incomplete or changed')
    receipt = verify_cm(plan, CM_PLAN, state)
    rows = {}
    for name in receipt['output_sha256']:
        with (results/name).open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                rows[row['object_id']] = dict(
                    has_ephor=row['ephor_object_id'].strip() != '', has_mizuki=row['mizuki_object_id'].strip() != '',
                    **{k: num(row[k]) for k in (*CODES.values(), *DIAGNOSTICS)})
    return rows, dict(cmodel_plan_sha256=plan['query_plan_sha256'], cmodel_receipt_sha256=sha(results/'download_summary.json'))


def nmad(x):
    return float(1.4826*np.median(np.abs(x-np.median(x)))) if len(x) else None


def compare(z_code, z_ref, rng):
    """Bin statistics of z_code - z_ref (None-safe when z_ref is absent)."""
    out = dict(n=len(z_code), median_z_code=float(np.median(z_code)),
               z_code_q16=float(np.percentile(z_code, 16)), z_code_q84=float(np.percentile(z_code, 84)))
    if z_ref is None:
        return out
    dz = z_code-z_ref
    scaled = dz/(1+z_ref)
    index = rng.integers(0, len(dz), size=(BOOTSTRAP, len(dz)))
    out.update(median_z_ref=float(np.median(z_ref)), median_dz=float(np.median(dz)),
               median_dz_bootstrap_se=float(np.std(np.median(dz[index], axis=1), ddof=1)),
               nmad_scaled=nmad(scaled), outlier_fraction=float(np.mean(np.abs(scaled) > OUTLIER)))
    return out


def fit_slope(fwhm, dz, winsor=True):
    if winsor:
        dz = np.clip(dz, *np.percentile(dz, WINSOR))
    return float(np.polyfit(fwhm, dz, 1)[0])


def tract_jackknife(fwhm, dz, tract, winsor=True):
    """Leave-one-tract-out slopes: FWHM is mostly a field property, so tracts are the resampling unit."""
    replicas = [fit_slope(fwhm[tract != t], dz[tract != t], winsor) for t in np.unique(tract)]
    return dict(se=jackknife(replicas), n_tracts=len(replicas), lo=min(replicas), hi=max(replicas))


def slope(fwhm, dz, scaled, tract, rng):
    clipped = np.clip(dz, *np.percentile(dz, WINSOR))
    fit = lambda f, d: float(np.polyfit(f, d, 1)[0])
    index = rng.integers(0, len(dz), size=(BOOTSTRAP, len(dz)))
    out = dict(slope=fit(fwhm, clipped), slope_source_iid_bootstrap_se=float(np.std([fit(fwhm[i], clipped[i]) for i in index], ddof=1)))
    jk = tract_jackknife(fwhm, dz, tract)
    out.update(slope_tract_jackknife_se=jk['se'], n_tracts=jk['n_tracts'], slope_leave_one_tract_out_min=jk['lo'],
               slope_leave_one_tract_out_max=jk['hi'])
    core = np.abs(scaled) <= OUTLIER
    jk = tract_jackknife(fwhm[core], dz[core], tract[core], winsor=False)
    out.update(core_n=int(core.sum()), core_slope=fit_slope(fwhm[core], dz[core], winsor=False), core_slope_tract_jackknife_se=jk['se'])
    return out


def fixed_weight_stack(parts):
    """Within-lens changes averaged with each lens's all-source sum w eff^2 held fixed (no lens-mix term)."""
    lids = list(parts)
    weight = np.array([parts[l][0]['all'][2] for l in lids])

    def changes(sums):
        a, b = sums['all'], sums['available']
        return np.array([b[0]/b[2]-a[0]/a[2], (b[1]/b[2])/(a[1]/a[2])-1])
    point = np.array([changes(parts[l][0]) for l in lids])
    boot = np.array([changes({s: parts[l][1][s].T for s in ('all', 'available')}) for l in lids])
    mean = lambda use: (weight[use, None]*point[use]).sum(0)/weight[use].sum()
    value = mean(np.ones(len(lids), bool))
    replicas = np.array([mean(np.arange(len(lids)) != k) for k in range(len(lids))])
    spread = (weight[:, None, None]*boot).sum(0)/weight.sum()
    return dict(mean_z_available_minus_all=float(value[0]), mean_eff_available_over_all_minus_1=float(value[1]),
                mean_z_available_minus_all_bootstrap_se=float(np.std(spread[0], ddof=1)),
                mean_eff_available_over_all_minus_1_bootstrap_se=float(np.std(spread[1], ddof=1)),
                mean_z_available_minus_all_jackknife_se=jackknife(replicas[:, 0]),
                mean_eff_available_over_all_minus_1_jackknife_se=jackknife(replicas[:, 1]))


def tested(values, group, label, rng):
    ok = np.isfinite(values)
    a, b = int((ok & label).sum()), int((ok & ~label).sum())
    if a == 0 or b == 0:
        return dict(n_true=a, n_false=b, testable=False)
    result = stratified_test(values, group, label, rng)
    return dict(result, testable=True)


def holm(rows, key='family'):
    for family in {r[key] for r in rows}:
        members = sorted((r for r in rows if r[key] == family and r.get('testable')), key=lambda r: r['permutation_p'])
        running = 0.0
        for rank, r in enumerate(members):
            running = max(running, min(1.0, (len(members)-rank)*r['permutation_p']))
            r['holm_p_within_family'] = running
            r['family_size'] = len(members)


def weighted_sums(w, z, eff, counts=None):
    """Sums for the efficiency-weighted mean z and mean efficiency: (sum w eff^2 z, sum w eff^3, sum w eff^2)."""
    terms = np.c_[w*eff**2*z, w*eff**3, w*eff**2]
    if counts is None:
        return terms.sum(0)
    return np.concatenate([c @ terms for c in counts])


def run(output):
    if output.exists():
        raise FileExistsError(output)
    rng = np.random.default_rng(SEED)
    cosmo = Cosmology()
    cm, cm_provenance = load_cmodel()
    seeing, unavailable, ab_provenance = load_afterburner()
    props = load_supplement()
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    point_z = finite_point_redshifts(available)
    if set(cm) != set(seeing) or set(cm) != set(props):
        raise ValueError('cmodel, afterburner and supplement object sets differ')
    if {o for o in selected if o not in available} != unavailable & selected:
        raise ValueError('afterburner mask differs from the photo-z availability ledger')

    pos = {}
    for p in pairs:
        if p['object_id'] in selected:
            pos.setdefault(p['object_id'], (float(p['source_ra']), float(p['source_dec'])))
    ids = sorted(pos, key=int)
    if set(ids) != selected:
        raise ValueError('selected objects without a selected pair')
    missing = np.array([o not in available for o in ids])
    fwhm_r = np.array([seeing[o][BANDS.index('r')]*FWHM_PER_SIGMA*PIXEL_ARCSEC for o in ids])
    z_ab = np.array([point_z.get(o, math.nan) for o in ids])
    z = {code: np.array([cm[o][col] for o in ids]) for code, col in CODES.items()}
    has = {'ephor_cmodel': np.array([cm[o]['has_ephor'] for o in ids]), 'mizuki': np.array([cm[o]['has_mizuki'] for o in ids])}
    r_deconv = fwhm_r > DECONV_FWHM_ARCSEC
    seeing_all = np.array([seeing[o] for o in ids], float)
    any_deconv = np.any(seeing_all*FWHM_PER_SIGMA*PIXEL_ARCSEC > DECONV_FWHM_ARCSEC, axis=1)
    any_seeing_nan = ~np.all(np.isfinite(seeing_all), axis=1)
    tract = np.array([int(o) >> 42 for o in ids])

    # Row existence by group and missingness cause.
    groups = {'ephor_ab_missing': missing, 'ephor_ab_missing_r_fwhm_gt_threshold': missing & r_deconv,
              'ephor_ab_missing_r_fwhm_le_threshold': missing & ~r_deconv, 'ephor_ab_finite': ~missing,
              'ephor_ab_finite_r_fwhm_gt_1p0': ~missing & (fwhm_r > NEAR_FWHM_ARCSEC)}
    existence = []
    for gname, g in groups.items():
        row = dict(group=gname, n=int(g.sum()), r_fwhm_min=float(np.nanmin(fwhm_r[g])), r_fwhm_median=float(np.nanmedian(fwhm_r[g])),
                   r_fwhm_max=float(np.nanmax(fwhm_r[g])), any_band_fwhm_gt_threshold=int((g & any_deconv).sum()),
                   any_band_seeing_nan=int((g & any_seeing_nan).sum()))
        for code in CODES:
            row[code+'_row'] = int((g & has[code]).sum())
            row[code+'_finite_z'] = int((g & np.isfinite(z[code])).sum())
            row[code+'_negative_z'] = int((g & (z[code] < 0)).sum())
        existence.append(row)

    # 1. Code minus Ephor_AB by r FWHM, and the code redshift of missing sources in the same bins.
    fwhm_rows = []
    for code in CODES:
        both = ~missing & np.isfinite(z[code]) & np.isfinite(z_ab)
        for k in range(len(FWHM_EDGES)-1):
            lo, hi = FWHM_EDGES[k], FWHM_EDGES[k+1]
            in_bin = (fwhm_r > lo) & (fwhm_r <= hi)
            for gname, g in (('both_available', both), ('ephor_ab_missing', missing & np.isfinite(z[code]))):
                sel = g & in_bin
                row = dict(code=code, group=gname, r_fwhm_lo=lo, r_fwhm_hi=hi, n=int(sel.sum()))
                if sel.sum() >= MIN_BIN:
                    row.update(compare(z[code][sel], z_ab[sel] if gname == 'both_available' else None, rng))
                fwhm_rows.append(row)
        near = both & (fwhm_r > NEAR_FWHM_ARCSEC) & (fwhm_r <= DECONV_FWHM_ARCSEC)
        fwhm_rows.append(dict(code=code, group='both_available_near_threshold', r_fwhm_lo=NEAR_FWHM_ARCSEC,
                              r_fwhm_hi=DECONV_FWHM_ARCSEC, **compare(z[code][near], z_ab[near], rng)))
        fit = both & (fwhm_r > SLOPE_RANGE[0]) & (fwhm_r <= SLOPE_RANGE[1])
        dz = (z[code]-z_ab)[fit]
        fwhm_rows.append(dict(code=code, group='both_available_slope_dz_per_arcsec', r_fwhm_lo=SLOPE_RANGE[0],
                              r_fwhm_hi=SLOPE_RANGE[1], n=int(fit.sum()),
                              **slope(fwhm_r[fit], dz, dz/(1+z_ab[fit]), tract[fit], rng)))

    # Diagnostics of the fits by group.
    diag_rows = []
    for gname in ('ephor_ab_missing', 'ephor_ab_finite_r_fwhm_gt_1p0', 'ephor_ab_finite'):
        g = groups[gname]
        for col in DIAGNOSTICS:
            v = np.array([cm[o][col] for o in np.asarray(ids)[g]])
            v = v[np.isfinite(v)]
            diag_rows.append(dict(group=gname, column=col, n_finite=len(v),
                                  **(dict(median=float(np.median(v)), q16=float(np.percentile(v, 16)), q84=float(np.percentile(v, 84)))
                                     if len(v) else {})))

    # 2. Same-cell tests, missing vs finite (primary) and missing vs near-threshold finite (secondary).
    ra, dec = np.array([pos[o] for o in ids]).T
    cell = (np.floor(ra/CELL_DEG).astype(np.int64)*100000+np.floor(dec/CELL_DEG).astype(np.int64)+50000)
    imag = np.array([num(props[o]['i_mag'])-num(props[o]['a_i']) for o in ids])
    imag_bin = np.where(np.isfinite(imag), np.floor(imag/0.5), -1).astype(np.int64)
    sizes, frac = Counter(cell.tolist()), Counter(cell[missing].tolist())
    mixed = np.array([sizes[c] >= MIN_CELL and 0.1 <= frac[c]/sizes[c] <= 0.9 for c in cell.tolist()])
    test_rows = []
    near_finite = ~missing & (fwhm_r > NEAR_FWHM_ARCSEC)
    for code in CODES:
        for strata, group in (('cell', cell), ('cell_x_imag0.5', cell*1000+imag_bin)):
            test_rows.append(dict(family='primary_missing_vs_finite', test=code+'_z_best', scope='mixed_cells', strata=strata,
                                  **tested(z[code][mixed], group[mixed], missing[mixed], rng)))
            scope = missing | near_finite
            test_rows.append(dict(family='secondary_missing_vs_finite_r_fwhm_gt_1p0', test=code+'_z_best',
                                  scope='cells_with_both_groups', strata=strata,
                                  **tested(z[code][scope], group[scope], missing[scope], rng)))
    holm(test_rows)

    # 3. Per-lens efficiency-weighted mean z, all selected vs Ephor_AB-available only.
    by_center = defaultdict(list)
    for p in pairs:
        if p['object_id'] in selected:
            by_center[p['center_id']].append(p)
    index_of = {o: i for i, o in enumerate(ids)}
    lenses = [c for c in centers if c['role'] == 'lens']
    lens_rows, stack_parts = [], defaultdict(dict)
    for code in CODES:
        for lens in lenses:
            lid, zl = lens['matched_lens_id'], float(lens['z_lens'])
            rows = by_center[lens['center_id']]
            if not rows:
                continue
            idx = np.array([index_of[p['object_id']] for p in rows])
            w = np.array([float(p['weight']) for p in rows])
            zs = z[code][idx]
            use = np.isfinite(zs) & (zs > zl+GAP)
            eff = np.zeros(len(zs))
            eff[use] = 1-float(cosmo.comoving(zl))/cosmo.comoving(zs[use])
            zs = np.where(use, zs, 0.0)
            avail = ~missing[idx]
            counts = [rng.poisson(1.0, size=(CHUNK, len(w))).astype(float) for _ in range(BOOTSTRAP//CHUNK)]
            sums = {s: weighted_sums(w*use*keep, zs, eff) for s, keep in (('all', np.ones(len(w))), ('available', avail))}
            boot = {s: weighted_sums(w*use*keep, zs, eff, counts).reshape(-1, 3) for s, keep in (('all', np.ones(len(w))), ('available', avail))}
            stack_parts[code][lid] = (sums, boot)
            lens_rows.append(dict(code=code, lens_id=lid, z_lens=zl, sources=len(w), sources_used=int(use.sum()),
                                  missing_sources_used=int((use & ~avail).sum()), **summarise(sums, boot)))
        lids = list(stack_parts[code])
        total = lambda subset, part: {s: sum(stack_parts[code][l][part][s] for l in subset) for s in ('all', 'available')}
        row = dict(code=code, lens_id='STACK:all', sources_used=None, **summarise(total(lids, 0), total(lids, 1)))
        for key in ('mean_z_available_minus_all', 'mean_eff_available_over_all_minus_1'):
            row[key+'_jackknife_se'] = jackknife([summarise(total([l for l in lids if l != d], 0), None)[key] for d in lids])
        lens_rows.append(row)
        lens_rows.append(dict(code=code, lens_id='STACK:fixed_all_source_lens_weight', **fixed_weight_stack(stack_parts[code])))

    output.mkdir(parents=True)
    tables = {'row_existence.tsv': existence, 'fwhm_comparison.tsv': fwhm_rows, 'fit_diagnostics.tsv': diag_rows,
              'stratified_tests.tsv': test_rows, 'lens_mean_redshift.tsv': lens_rows}
    for name, rows in tables.items():
        write_table(output/name, rows, list(dict.fromkeys(k for r in rows for k in r)))
    report = dict(
        science_measurement=False, calibration_applied=False, selection_changed=False,
        constants=dict(seed=SEED, fwhm_edges_arcsec=[e if math.isfinite(e) else 'inf' for e in FWHM_EDGES],
                       deconv_fwhm_arcsec=DECONV_FWHM_ARCSEC, near_threshold_fwhm_arcsec=NEAR_FWHM_ARCSEC, slope_range_arcsec=SLOPE_RANGE,
                       outlier_scaled_dz=OUTLIER, min_bin=MIN_BIN, gap=GAP, bootstrap=BOOTSTRAP, permutations=PERMUTATIONS,
                       winsor_percentiles=WINSOR, cell_deg=CELL_DEG, min_cell=MIN_CELL),
        unique_selected_objects=len(ids), ephor_ab_missing=int(missing.sum()), mixed_cell_objects=int(mixed.sum()),
        mixed_cell_objects_by_tract={str(t): n for t, n in Counter(tract[mixed].tolist()).most_common()},
        missing_by_tract={str(t): n for t, n in Counter(tract[missing].tolist()).most_common()},
        definitions=dict(
            fwhm='afterburner r parent_flux_convolved_seeing x 2.3548 x 0.168 arcsec; deconvolution threshold 6.5 px = 1.092 arcsec.',
            dz='code z_best minus Ephor_AB z_best; scaled dz divides by 1 + Ephor_AB z; outlier |scaled dz| > 0.15.',
            tests='Stratum-size-weighted mean difference (missing minus finite) of winsorised z_best, within-stratum permutation '
                  'null, Holm within each family. Primary: mixed cells (>= MIN_CELL objects, missing fraction 0.1-0.9). '
                  'Secondary: missing vs finite with r FWHM > 1.0 arcsec in any stratum holding both.',
            lens_mean_z='Sources of each lens with code z > z_lens + 0.2; weights w (1 - chi_l/chi_s)^2 from the same code for both '
                        'sets; mean_eff is sum w eff^3 / sum w eff^2. Bootstrap: Poisson source counts per lens; stack '
                        'jackknife over lenses. STACK:all pools sums over lenses, so it mixes within-lens removal with the '
                        'change in lens weights (S23B_3620 keeps about 0.74 of its weight); STACK:fixed_all_source_lens_weight '
                        'averages the per-lens changes with all-source lens weights held fixed.',
            slope='OLS of winsorised dz on r FWHM over 0.8-1.092 arcsec. The main error is the leave-one-tract-out '
                  'jackknife (tract = object_id >> 42); the source-iid bootstrap is kept for reference only. core_slope '
                  'fits sources with |scaled dz| <= 0.15 without winsorising.',
            missing_groups='r_fwhm_le_threshold: missing sources whose r FWHM is below the threshold; any_band columns '
                           'show another band above the threshold or a band without seeing (no coadd input).'),
        limitations=['Ephor(cmodel) and Mizuki use non-PSF-matched photometry; their own seeing dependence is measured only up to '
                     'the threshold, so the missing sources (above it) are an extrapolation.',
                     'FWHM is mostly a whole-field property; FWHM trends are confounded with field identity.',
                     'Code z_best is a point estimate; differences are in that estimator, not in true redshift.',
                     'Sources are reused across lens and random centres; per-lens bootstraps ignore that and clustering.',
                     'Mixed cells lie mostly in one tract (the S23B_3620 field); the same-cell result is close to a one-field result.',
                     'The per-lens mean z responds weakly to a redshift offset of the missing sources (the z > z_l + 0.2 cut '
                     'damps it); it measures how much the weighted sample changes, not a bound on that offset.',
                     'Nothing here tests whether the COSMOS calibration of the multiplier transfers to near-threshold seeing.'],
        provenance=dict(**provenance, **ab_provenance, **cm_provenance), runner_sha256=sha(Path(__file__)),
        output_sha256={name: sha(output/name) for name in tables})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(objects=len(ids), missing=int(missing.sum()), mixed=int(mixed.sum())), indent=2))


def summarise(sums, boot):
    def values(s):
        return dict(mean_z=s[..., 0]/s[..., 2], mean_eff=s[..., 1]/s[..., 2], weight=s[..., 2])
    a, b = values(sums['all']), values(sums['available'])
    out = dict(mean_z_all=float(a['mean_z']), mean_z_available=float(b['mean_z']),
               mean_z_available_minus_all=float(b['mean_z']-a['mean_z']),
               mean_eff_all=float(a['mean_eff']), mean_eff_available=float(b['mean_eff']),
               mean_eff_available_over_all_minus_1=float(b['mean_eff']/a['mean_eff']-1),
               retained_eff2_weight_fraction=float(b['weight']/a['weight']))
    if boot is not None:
        ba, bb = values(boot['all']), values(boot['available'])
        out['mean_z_available_minus_all_bootstrap_se'] = float(np.std(bb['mean_z']-ba['mean_z'], ddof=1))
        out['mean_eff_available_over_all_minus_1_bootstrap_se'] = float(np.std(bb['mean_eff']/ba['mean_eff']-1, ddof=1))
    return out


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
