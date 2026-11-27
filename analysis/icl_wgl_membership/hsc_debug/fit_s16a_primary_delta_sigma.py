#!/usr/bin/env python3
"""Covariance and stacked NFW mass fit of the S16A primary-definition Delta Sigma (measure_s16a_primary_delta_sigma). Local only.

Fixed before the run (user: "1-->2", 2026-11-27):
1. Covariance per sample (all colours, blue, all without S23B_9307), tangential and cross, 8 proper-radius bins:
   - C_shape: exact linear propagation of per-component shape noise (rms_e^2 + sigma_e^2) through the estimator, with every
     source counted once however many lens/random discs reuse it (shared-source terms included);
   - C_jk: group jackknife (11 groups, lens with its randoms);
   - C_pz: COSMOS object bootstrap of the photo-z multiplier (200, seed 2026112701, as the measurement);
   - adopted C = C_shape + diag(max(0, diag C_jk - diag C_shape)) + C_pz (halo-to-halo/LSS excess on the diagonal only);
     sensitivity: C_jk + C_pz with the Hartlap factor (n - p - 2)/(n - 1) on the inverse.
2. Fit 0.4-3 Mpc (5 bins), one common M200c (flat prior in log10 M, 13-15.5), each lens modelled at its own z_l with
   ogfkit HaloModel (NFW + linear two-halo with Tinker bias), Duffy et al. 2008 c(M, z); model in a bin = w_ls-weighted
   mean of DeltaSigma_l(R) over the lens pairs of the bin. Variants: C_jk-Hartlap, no two-halo, miscentring
   (f_cen 0.7, sigma_off 0.2 Mpc), free c (flat in log c, 1-20), blue, without S23B_9307, and an extra-light upper-bound
   correction (all colours, DeltaSigma x 0.76 at 0.4-0.7 Mpc and x 0.90 at 0.7-1.5 Mpc). Shear treated as unreduced.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import coo_matrix
from scipy.stats import chi2 as chi2_dist

from audit_s23b_s16a_shape_qa_pilot import RADIAL_EDGES
from measure_s16a_primary_delta_sigma import BOOT, OUTLIER, SEED, load_primary
from select_s16a_random_centres import sha

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wl_photoz import stack_photoz_response  # noqa: E402
from ogfkit.wlcluster import Cosmology, HaloModel  # noqa: E402

MEASURED = Path('/home/kjhan/BACKUP/Astrafek/s16a_primary_delta_sigma_20261127_01/delta_sigma.tsv')
MEASURED_SHA = '4c874fe9a6e5e04f22397212055d81957cdca30afbf753ff569f8799a506ffe4'
FIT_BINS = [3, 4, 5, 6, 7]
LOGM = np.arange(13.0, 15.5001, 0.01)
LOGC = np.linspace(0.0, np.log10(20.0), 41)
RGRID = np.logspace(-2, np.log10(4.0), 160)
EXTRA_LIGHT = {3: 0.76, 4: 0.90, 5: 0.90}


def run(output):
    if output.exists():
        raise FileExistsError(output)
    if sha(MEASURED) != MEASURED_SHA:
        raise ValueError('measurement table changed')
    with MEASURED.open() as stream:
        measured = list(csv.DictReader(stream, delimiter='\t'))
    cosmo = Cosmology()
    S = load_primary(cosmo)
    data, groups, lens_ids, new_rand, lens_resp, responses, nb = (S[k] for k in
        ('data', 'groups', 'lens_ids', 'new_rand', 'lens_resp', 'responses', 'nb'))
    rng = np.random.default_rng(SEED)
    boot = [{g: responses(lens_resp[g]['lens_z'], idx)[0] for g in groups} for idx in (rng.integers(0, S['nC'], S['nC']) for _ in range(BOOT))]
    samples = {'all colours': (groups, 'all'), 'blue': (groups, 'blue'), f'all colours, no {OUTLIER}': ([g for g in groups if g != OUTLIER], 'all')}

    def pairs(gs, cname):
        """Per-pair arrays for the lens and random sides of the sample."""
        out = []
        for side, cids in (('lens', [lens_ids[g] for g in gs]), ('random', [r for g in gs for r in new_rand[g]])):
            for c in cids:
                d = data[c]
                sel = (d['bin'] >= 0) & (d['bin'] < nb) & ((d['gi'] <= S['t1']) if cname == 'blue' else True)
                gl = next(g for g in gs if c == lens_ids[g] or c in new_rand[g])
                out.append(dict(side=side, group=gl, cid=c, **{k: d[k][sel] for k in ('ids', 'phi', 'bin', 'et', 'ex', 'sinv', 'wls', 'm', 'rms', 'sig', 'rprop')}))
        return out

    def coefficients(P, gs, resp):
        """a_p with DeltaSigma_k = sum a_p e_t,p (tangential), per side normalisation and stacked multiplier."""
        norm, wl = {}, np.zeros((len(gs), nb))
        for side in ('lens', 'random'):
            w = np.zeros(nb); wm = np.zeros(nb); wr = np.zeros(nb)
            for p in P:
                if p['side'] == side:
                    np.add.at(w, p['bin'], p['wls']); np.add.at(wm, p['bin'], p['wls']*p['m']); np.add.at(wr, p['bin'], p['wls']*p['rms']**2)
                    if side == 'lens':
                        np.add.at(wl[gs.index(p['group'])], p['bin'], p['wls'])
            norm[side] = 2*(1-wr/np.where(w > 0, w, 1))*(1+wm/np.where(w > 0, w, 1))*w
        mult = np.array([1/stack_photoz_response([resp[g] for g in gs], wl[:, k])['response'] if wl[:, k].sum() > 0 else np.nan for k in range(nb)])
        return [(mult*(1 if p['side'] == 'lens' else -1)/norm[p['side']])[p['bin']]*p['wls']/p['sinv'] for p in P], mult

    def vector(P, A):
        t, x = np.zeros(nb), np.zeros(nb)
        for p, a in zip(P, A):
            np.add.at(t, p['bin'], a*p['et']); np.add.at(x, p['bin'], a*p['ex'])
        return t, x
    base_resp = {g: lens_resp[g]['primary_response'] for g in groups}
    results, fits, covs = [], [], {}
    for label, (gs, cname) in samples.items():
        P = pairs(gs, cname)
        A, mult = coefficients(P, gs, base_resp)
        t, x = vector(P, A)
        ref = [r for r in measured if r['sample'] == label]
        if not np.allclose(t, [float(r['delta_sigma_t']) for r in ref], rtol=1e-9, atol=1e-3):
            raise ValueError('linear estimator does not reproduce the measurement: '+label)
        # shape noise, shared sources once
        ids = np.concatenate([p['ids'] for p in P]); _, inv = np.unique(ids, return_inverse=True)
        phi = np.concatenate([p['phi'] for p in P]); b = np.concatenate([p['bin'] for p in P]); a = np.concatenate(A)
        s = np.sqrt(np.concatenate([p['rms']**2+p['sig']**2 for p in P]))
        c2, s2 = np.cos(2*phi), np.sin(2*phi)
        rows_ = np.r_[b, b, b+nb, b+nb]; cols = np.r_[2*inv, 2*inv+1, 2*inv, 2*inv+1]
        vals = np.r_[-a*s*c2, -a*s*s2, a*s*s2, -a*s*c2]
        op = coo_matrix((vals, (rows_, cols)), shape=(2*nb, 2*(inv.max()+1))).tocsr()
        Csh = (op @ op.T).toarray()
        reps = []
        for drop in gs:
            g2 = [g for g in gs if g != drop]
            P2 = [p for p in P if p['group'] != drop]
            A2, _ = coefficients(P2, g2, base_resp)
            reps.append(np.r_[vector(P2, A2)])
        reps = np.array(reps); n = len(gs)
        Cjk = (n-1)/n*(reps-reps.mean(0)).T @ (reps-reps.mean(0))
        raw_t, raw_x = t/mult, x/mult
        breps = []
        for bb in boot:
            _, mb = coefficients(P, gs, bb)
            breps.append(np.r_[raw_t*mb, raw_x*mb])
        breps = np.array(breps)
        Cpz = np.cov(breps.T, ddof=1)
        Cad = Csh+np.diag(np.maximum(0, np.diag(Cjk)-np.diag(Csh)))+Cpz
        covs[label] = dict(shape=Csh, jackknife=Cjk, photoz=Cpz, adopted=Cad, n_groups=n)
        for k in range(nb):
            results.append(dict(sample=label, r_lo=RADIAL_EDGES[k], r_hi=RADIAL_EDGES[k+1], in_fit_range=k in FIT_BINS,
                                delta_sigma_t=t[k], se_shape=np.sqrt(Csh[k, k]), se_jackknife=np.sqrt(Cjk[k, k]), se_photoz=np.sqrt(Cpz[k, k]),
                                se_adopted=np.sqrt(Cad[k, k]), delta_sigma_x=x[k], se_x_adopted=np.sqrt(Cad[nb+k, nb+k]),
                                shared_source_shape_se_ratio=np.sqrt(Csh[k, k]/((a[b == k])**2*s[b == k]**2).sum())))
        covs[label]['t'], covs[label]['x'], covs[label]['P'] = t, x, P
    # model: per lens pairs in fit bins
    models = {}

    def model_vector(label, logm, logc=None, two_halo=True, miscenter=False):
        gs, cname = samples[label]
        key = (label, two_halo, miscenter)
        if key not in models:
            P = [p for p in covs[label]['P'] if p['side'] == 'lens']
            hm = {g: HaloModel(cosmo, lens_resp[g]['lens_z'], RGRID, two_halo=two_halo, miscenter=miscenter, f_cen=0.7, sigma_off=0.2) for g in gs}
            models[key] = (P, hm)
        P, hm = models[key]
        M = 10**logm
        num_, den = np.zeros(nb), np.zeros(nb)
        prof = {}
        for g in gs:
            c = 10**logc if logc is not None else float(cosmo.duffy_c(M, lens_resp[g]['lens_z']))
            prof[g] = hm[g].profiles(M, c)[1]
        for p in P:
            v = np.interp(np.log(p['rprop']), np.log(RGRID), prof[p['group']])
            np.add.at(num_, p['bin'], p['wls']*v); np.add.at(den, p['bin'], p['wls'])
        return num_/np.where(den > 0, den, 1)

    def fit(name, label, cov='adopted', two_halo=True, miscenter=False, free_c=False, scale=None):
        C = covs[label]
        y = C['t'][FIT_BINS].copy()
        if scale:
            y *= np.array([scale.get(k, 1.0) for k in FIT_BINS])
        if cov == 'adopted':
            Cm = C['adopted'][np.ix_(FIT_BINS, FIT_BINS)]; icov = np.linalg.inv(Cm); hart = 1.0
        else:
            Cm = (C['jackknife']+C['photoz'])[np.ix_(FIT_BINS, FIT_BINS)]
            hart = (C['n_groups']-len(FIT_BINS)-2)/(C['n_groups']-1); icov = hart*np.linalg.inv(Cm)
        grid = [None] if not free_c else list(LOGC)
        def c2(lm, lc):
            r = y-model_vector(label, lm, lc, two_halo, miscenter)[FIT_BINS]
            return float(r @ icov @ r)
        chi = np.array([[c2(lm, lc) for lc in grid] for lm in LOGM])
        post = np.exp(-0.5*(chi-chi.min())); pm = post.sum(1); pm /= pm.sum(); cdf = np.cumsum(pm)
        i, j = np.unravel_index(np.argmin(chi), chi.shape)
        q = lambda f: float(np.interp(f, cdf, LOGM))
        dof = len(FIT_BINS)-(2 if free_c else 1)
        null = float(y @ icov @ y)
        out = dict(name=name, sample=label, covariance=cov, hartlap=hart, two_halo=two_halo, miscenter=miscenter, free_c=free_c,
                   extra_light_scaled=bool(scale), best_log10_M200c=float(LOGM[i]), median_log10_M200c=q(0.5),
                   lo68=q(0.16), hi68=q(0.84), best_log10_c=None if not free_c else float(LOGC[j]),
                   chi2_min=float(chi.min()), dof=dof, p_value=float(chi2_dist.sf(chi.min(), dof)),
                   chi2_null=null, detection_sigma=float(np.sqrt(max(0, null-chi.min()))),
                   at_prior_edge=bool(i in (0, len(LOGM)-1)))
        if not free_c:
            out['best_model'] = model_vector(label, LOGM[i], None, two_halo, miscenter).tolist()
        return out
    fits = [fit('baseline', 'all colours'),
            fit('jackknife+Hartlap', 'all colours', cov='jk'),
            fit('no two-halo', 'all colours', two_halo=False),
            fit('miscentring f=0.7 s=0.2', 'all colours', miscenter=True),
            fit('free c', 'all colours', free_c=True),
            fit('blue', 'blue'),
            fit(f'no {OUTLIER}', f'all colours, no {OUTLIER}'),
            fit('extra-light upper bound', 'all colours', scale=EXTRA_LIGHT)]
    chi_x = {lab: float(C['x'][FIT_BINS] @ np.linalg.inv(C['adopted'][np.ix_([nb+k for k in FIT_BINS], [nb+k for k in FIT_BINS])]) @ C['x'][FIT_BINS])
             for lab, C in covs.items()}
    output.mkdir(parents=True)
    with (output/'profile_errors.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(results)
    np.savez(output/'covariances.npz', **{f'{lab}|{k}': v for lab, C in covs.items() for k, v in C.items() if k in ('shape', 'jackknife', 'photoz', 'adopted')})
    summary = dict(schema='astrafek-s16a-primary-delta-sigma-fit-v1', science_measurement=True, new_query=False, measured_sha256=MEASURED_SHA,
                   fit_bins=[(RADIAL_EDGES[k], RADIAL_EDGES[k+1]) for k in FIT_BINS], log10_mass_grid=[LOGM[0], LOGM[-1]],
                   mass_definition='M200c, Msun (no h), flat LCDM H0=70 Om=0.3', fits=fits, cross_chi2_fit_range_adopted=chi_x,
                   correlation_adopted_all={'matrix': (lambda Cm: (Cm/np.sqrt(np.outer(np.diag(Cm), np.diag(Cm)))).tolist())(
                       covs['all colours']['adopted'][np.ix_(FIT_BINS, FIT_BINS)])},
                   limitations=['Halo-to-halo/LSS excess on the diagonal only; jackknife with 11 groups is noisy.',
                                'Unreduced shear; boost 1; crowding loss and extra-light bias not corrected in the baseline.',
                                'Single common mass for a stack spanning z 0.35-0.54; no mass-richness scatter.',
                                'COSMOS seeing transfer untested; COSMOS bootstrap ignores field variance.'],
                   runner_sha256=sha(Path(__file__)), output_sha256={n: sha(output/n) for n in ('profile_errors.tsv', 'covariances.npz')})
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return results, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_primary_delta_sigma_fit_20261127_01'))
    results, summary = run(parser.parse_args().output)
    u = 1e-12
    for r in results:
        print(f"{r['sample']:26s} {r['r_lo']:>4}-{r['r_hi']:<4} DS={r['delta_sigma_t']*u:7.2f} se shape {r['se_shape']*u:6.2f} "
              f"jk {r['se_jackknife']*u:6.2f} pz {r['se_photoz']*u:5.2f} adopted {r['se_adopted']*u:6.2f} "
              f"DSx={r['delta_sigma_x']*u:6.2f}±{r['se_x_adopted']*u:5.2f} shared/indep {r['shared_source_shape_se_ratio']:.3f}")
    for f in summary['fits']:
        print({k: v for k, v in f.items() if k != 'best_model'})
    print('cross chi2', summary['cross_chi2_fit_range_adopted'])
    print(np.round(np.array(summary['correlation_adopted_all']['matrix']), 2))
