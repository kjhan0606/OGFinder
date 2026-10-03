#!/usr/bin/env python3
"""Photo-z distribution tools: calibration (PIT), percentile accuracy and sample representativeness of a catalog.

    photoz_quality.py --catalog CAT --work DIR --mode quality [--zspec-col Z_SPEC | --zspec-file TSV] [--zphot-col PHOTO_Z]
        [--zerr-col PHOTO_Z_ERR | --q16-col EZ_Z16 --q84-col EZ_Z84 | --mixture-file NPZ | --pdf-file NPZ] [--pdf-kind auto|gauss|split|mixture|grid]
        [--zmin 0] [--outlier 0.15] [--mag-col MAG_AUTO] [--z-bins 0,0.1,0.2,..] [--mag-bins ..] [--weights-col PZ_SPECW] [--recalibrate]
    photoz_quality.py --catalog CAT --work DIR --mode repr [--zspec-col ..] [--features MAG_AUTO,MAG_F606W-MAG_F850LP,...] [--knn 5] [--coverage-q 95]

quality  rows with a spectroscopic redshift (the calibration sample): point metrics (bias, sigma_NMAD, outlier fraction, percentiles of |dz|, bootstrap errors,
         binned in z_spec / z_phot / magnitude), and, with a predictive distribution per object (Gaussian zphot+zerr, asymmetric from 16/84 percentile columns,
         Gaussian mixture of the MDN, or tabulated p(z)): PIT histogram + KS / Cramer-von Mises tests, credible-interval coverage, CRPS, log score, z-score,
         optimal width scale and PIT recalibration, stacked p(z) vs spec-z histogram.  Weighted by --weights-col (representativeness weights).
repr     is the spectroscopic sample representative of the target sample (all rows with the features)?  1-D KS / mean difference per feature, k-NN two-sample
         test in feature space (AUC + permutation p), distance of every target object to the spec sample (support), k-NN density-ratio weights (Lima et al. 2008)
         and the effective sample size; the weighted spectroscopic accuracy is the estimate for the target sample.
Catalog columns (add_columns): quality: PZ_PIT, PZ_CRPS, PZ_ZSCORE, PZ_DZ, PZ_INCI68, PZ_OUTLIER; repr: PZ_SPECDIST, PZ_INSPEC, PZ_SPECW.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, photoz_stats as ps, meta as ometa, pz_closure as pcl  # noqa: E402

QCOLS = ['PZ_PIT', 'PZ_CRPS', 'PZ_ZSCORE', 'PZ_DZ', 'PZ_INCI68', 'PZ_OUTLIER']
RCOLS = ['PZ_SPECDIST', 'PZ_INSPEC', 'PZ_SPECW']
CCOLS = ['PZ_ZDIFF', 'PZ_INCONSIST']
ZSPEC_NAMES = ('Z_SPEC', 'ZSPEC', 'SPEC_Z', 'SPECZ', 'Z_SPECTRO', 'REDSHIFT', 'Z')


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return clean(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    return o


def colvals(rows, name):
    return np.array([tsvio.fnum(r.get(name)) for r in rows])


def find_zspec(cols, rows, a, W):
    if a.zspec_file:
        c2, r2 = tsvio.read_catalog(a.zspec_file)
        name = a.zspec_col if a.zspec_col in c2 else next((n for n in ZSPEC_NAMES if n in c2), None)
        if name is None:
            sys.exit('ERROR: no spec-z column in %s (columns %s)' % (a.zspec_file, c2))
        m = {r['NUMBER']: tsvio.fnum(r.get(name)) for r in r2}
        return np.array([m.get(r['NUMBER'], float('nan')) for r in rows]), name
    name = a.zspec_col if a.zspec_col in cols else next((n for n in ZSPEC_NAMES if n in cols), None)
    if name is None:
        sys.exit('ERROR: spec-z column %r not found (set --zspec-col or --zspec-file)' % a.zspec_col)
    return colvals(rows, name), name


def build_predictive(cols, rows, a, nums):
    zp = colvals(rows, a.zphot_col) if a.zphot_col in cols else None
    mixf = a.mixture_file or ''
    if a.pdf_file and a.pdf_kind in ('auto', 'grid') and not mixf:
        d = np.load(a.pdf_file)
        if 'pdf' in d and 'zgrid' in d:
            mixf = ''
            idx = {int(n): i for i, n in enumerate(d['number'])}
            sel = np.array([idx.get(int(float(n)), -1) for n in nums])
            ok = sel >= 0
            pdf = np.zeros((len(nums), len(d['zgrid'])))
            pdf[ok] = d['pdf'][sel[ok]]
            pdf[~ok] = 1.0
            return ps.Predictive('grid', zmin=a.zmin, zgrid=d['zgrid'], pdf=pdf), 'grid', ok, zp
    if mixf:
        d = np.load(mixf)
        idx = {int(n): i for i, n in enumerate(d['number'])}
        sel = np.array([idx.get(int(float(n)), -1) for n in nums])
        ok = sel >= 0
        K = d['pi'].shape[1]
        pi = np.full((len(nums), K), 1.0 / K); mu = np.zeros((len(nums), K)); sg = np.ones((len(nums), K))
        pi[ok] = d['pi'][sel[ok]]; mu[ok] = d['mu'][sel[ok]]; sg[ok] = d['sigma'][sel[ok]]
        return ps.Predictive('mixture', zmin=a.zmin, pi=pi, mu=mu, sigma=sg), 'mixture', ok, zp
    kind = a.pdf_kind
    if kind == 'auto':
        kind = 'split' if (a.q16_col in cols and a.q84_col in cols) else 'gauss'
    if zp is None:
        sys.exit('ERROR: photo-z column %r not found' % a.zphot_col)
    if kind == 'split':
        q16 = colvals(rows, a.q16_col); q84 = colvals(rows, a.q84_col)
        ok = np.isfinite(zp) & np.isfinite(q16) & np.isfinite(q84)
        zq = np.where(ok, zp, 0.5)
        return ps.Predictive('split', zmin=a.zmin, zphot=zq, lo=np.where(ok, zp - q16, 0.1), hi=np.where(ok, q84 - zp, 0.1)), 'split', ok, zp
    err = colvals(rows, a.zerr_col) if a.zerr_col in cols else np.full(len(rows), float('nan'))
    if a.zerr_col not in cols:
        sys.exit('ERROR: photo-z error column %r not found (or give --q16-col/--q84-col or --mixture-file)' % a.zerr_col)
    ok = np.isfinite(zp) & np.isfinite(err) & (err > 0)
    return ps.Predictive('gauss', zmin=a.zmin, zphot=np.where(ok, zp, 0.5), sigma=np.where(ok, err, 0.1)), 'gauss', ok, zp


def plot_quality(path, rep, pit, zs, zp, pred_ok, stack, binned, label):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(10, 8))
    if pit is not None:
        h = rep['pit']['hist']
        n = len(h)
        ax[0, 0].bar((np.arange(n) + 0.5) / n, np.array(h) / max(sum(h), 1) * n, width=1.0 / n, color='C0', alpha=0.8)
        ax[0, 0].axhline(1, color='k', ls='--')
        ax[0, 0].set_title('PIT: KS p=%.3g, CvM p=%.3g\n%s' % (rep['pit']['ks_p'], rep['pit']['cvm_p'], rep['pit']['shape']), fontsize=9)
        ax[0, 0].set_xlabel('PIT'); ax[0, 0].set_ylabel('density')
        cv = rep['pit']['coverage']
        lv = [c['level'] for c in cv]
        ax[0, 1].plot([0, 1], [0, 1], 'k--')
        ax[0, 1].errorbar(lv, [c['observed'] for c in cv], yerr=[[c['observed'] - c['lo'] for c in cv], [c['hi'] - c['observed'] for c in cv]], fmt='o-')
        ax[0, 1].set_xlabel('nominal credible level'); ax[0, 1].set_ylabel('fraction of z_spec inside'); ax[0, 1].set_title('coverage')
        q = np.sort(pit)
        ax[1, 1].plot(np.linspace(0, 1, len(q)), q, 'C1'); ax[1, 1].plot([0, 1], [0, 1], 'k--')
        ax[1, 1].set_xlabel('uniform quantile'); ax[1, 1].set_ylabel('PIT quantile'); ax[1, 1].set_title('Q-Q of the PIT')
    pm = rep['point']
    ax[1, 0].plot(zs, zp, '.', ms=2, alpha=0.4)
    zz = np.array([0, max(np.nanmax(zs), 0.1)])
    ax[1, 0].plot(zz, zz, 'k-', lw=0.8); ax[1, 0].plot(zz, zz + 0.15 * (1 + zz), 'r:', lw=0.8); ax[1, 0].plot(zz, zz - 0.15 * (1 + zz), 'r:', lw=0.8)
    ax[1, 0].set_xlabel('z_spec'); ax[1, 0].set_ylabel('z_phot')
    ax[1, 0].set_title('N=%d  bias %.4f  NMAD %.4f  outl %.1f%%' % (pm['n'], pm['bias_mean'], pm['sigma_nmad'], 100 * pm['outlier_frac']), fontsize=9)
    fig.suptitle(label, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def do_quality(cols, rows, a, W):
    nums = [r['NUMBER'] for r in rows]
    zs, zname = find_zspec(cols, rows, a, W)
    pred, kind, pok, zp = build_predictive(cols, rows, a, nums)
    if zp is None or not np.any(np.isfinite(zp)):
        zp = pred.mean_std()[0] if pred is not None else None
    if a.zphot_col not in cols:                                           # point estimate from the predictive (median)
        zp = pred.quantile(0.5) if kind in ('grid', 'split') else pred.mean_std()[0]
    ok = np.isfinite(zs) & (zs > 0) & np.isfinite(zp) & pok
    w = None
    if a.weights_col and a.weights_col in cols:
        wall = colvals(rows, a.weights_col)
        ok &= np.isfinite(wall) & (wall > 0)
        w = wall[ok]
    if ok.sum() < 20:
        sys.exit('ERROR: only %d objects with a spec-z and a photo-z (need >= 20)' % ok.sum())
    idx = np.where(ok)[0]
    sub = Predictive_subset(pred, idx)
    rep = dict(n_catalog=len(rows), n_spec=int(ok.sum()), zspec_column=zname, pdf_kind=kind, zmin=a.zmin, weighted=w is not None)
    rep['point'] = ps.point_metrics(zp[idx], zs[idx], a.outlier)
    rep['point'].update(ps.bootstrap_point(zp[idx], zs[idx], n=a.nboot, outlier=a.outlier))
    if w is not None:
        rep['point_weighted'] = ps.point_metrics(zp[idx], zs[idx], a.outlier, weights=w)
    # PDF part
    summ, pit, cr, ls, zsc = ps.calibration_summary(sub, zs[idx])
    rep.update(summ)
    rep['recalibration'] = {}
    if kind != 'grid':
        rep['recalibration'] = ps.optimal_scale(sub, zs[idx])
        sc = rep['recalibration']['scale']
        rep['recalibration']['pit_ks_p_scaled'] = float(__import__('scipy.stats', fromlist=['kstest']).kstest(sub.scaled(sc).cdf(zs[idx]), 'uniform').pvalue)
        rep['recalibration']['crps_scaled'] = float(np.mean(ps.crps(sub.scaled(sc), zs[idx])))
    # PIT recalibration with 2-fold cross-fit (an honest estimate of what PIT-recalibration achieves)
    rng = np.random.default_rng(1)
    f = rng.integers(0, 2, len(pit))
    pc = np.where(f == 0, ps.recalibrate_pit(pit[f == 1], pit), ps.recalibrate_pit(pit[f == 0], pit))
    from scipy import stats as sst
    rep['recalibration']['pit_crossfit_ks_p'] = float(sst.kstest(pc, 'uniform').pvalue)
    # closure: K-fold PIT recalibration fed back into the predictive (p'(z) = p(z) g(F(z))), out-of-fold scores before / after; the map learnt on all objects can be saved
    # (--recal-out) and applied to another catalogue of the same photo-z product (--recal-in)
    try:
        cf = pcl.crossfit(sub, zs[idx], k=5)
        rep['closure'] = dict(before=cf['before'], after_oof=cf['after_oof'], note='out-of-fold: each fold recalibrated with a map learnt on the other folds')
        if a.recal_out:
            with open(a.recal_out, 'w') as fh:
                json.dump(dict(map=cf['map'].to_dict(), n=int(len(idx)), pdf_kind=kind, source=os.path.basename(a.catalog)), fh)
    except Exception as e:                                                 # a failing closure must not break the quality report
        rep['closure'] = dict(error=str(e))
    if a.recal_in and os.path.isfile(a.recal_in):
        rcm = pcl.PITRecal.from_dict(json.load(open(a.recal_in))['map'])
        rep['recal_applied'] = pcl.summary(pcl.recalibrated_predictive(sub, rcm), zs[idx])[0]
        rows_q = []
        okall = np.where(pok & np.isfinite(zp))[0] if hasattr(pok, '__len__') else np.arange(len(rows))
        for s0 in range(0, len(okall), 4000):
            ch = okall[s0:s0 + 4000]
            pr = pcl.recalibrated_predictive(Predictive_subset(pred, ch), rcm)
            q16, q50, q84 = pr.quantile(0.16), pr.quantile(0.5), pr.quantile(0.84)
            rows_q += [dict(NUMBER=nums[j], PZ_P16_RC=q16[m], PZ_P50_RC=q50[m], PZ_P84_RC=q84[m]) for m, j in enumerate(ch)]
        os.makedirs(a.work, exist_ok=True)
        tsvio.write_table(os.path.join(a.work, 'pzq_recalibrated.tsv'), ['NUMBER', 'PZ_P16_RC', 'PZ_P50_RC', 'PZ_P84_RC'], rows_q)
    # binned
    zsi, zpi = zs[idx], zp[idx]
    zb = [float(v) for v in a.z_bins.split(',')] if a.z_bins else list(np.round(np.linspace(0, max(zsi.max(), 0.1) * 1.0001, 6), 3))
    rep['binned_zspec'] = ps.binned_metrics(zsi, zpi, zsi, zb, sub, a.outlier)
    rep['binned_zphot'] = ps.binned_metrics(zpi, zpi, zsi, zb, sub, a.outlier)
    if a.mag_col and a.mag_col in cols:
        mg = colvals(rows, a.mag_col)[idx]
        mbins = [float(v) for v in a.mag_bins.split(',')] if a.mag_bins else list(np.round(np.nanpercentile(mg, np.linspace(0, 100, 6)), 2))
        mbins[-1] += 1e-6
        rep['binned_mag'] = ps.binned_metrics(mg, zpi, zsi, mbins, sub, a.outlier)
        rep['mag_column'] = a.mag_col
    zg = np.linspace(0.0, max(zsi.max(), zpi.max()) * 1.3 + 0.1, 800)
    st = ps.stacked_nz(sub, zsi, zg, np.linspace(0, zg[-1], 21))
    rep['stacked_nz'] = dict(ks_distance=st['ks_distance'], mean_stacked=st['mean_stacked'], mean_spec=st['mean_spec'])
    # tables
    os.makedirs(a.work, exist_ok=True)
    pit_rows = [dict(lo=i / len(rep['pit']['hist']), hi=(i + 1) / len(rep['pit']['hist']), count=c, expected=rep['pit']['hist_expected']) for i, c in enumerate(rep['pit']['hist'])]
    tsvio.write_table(os.path.join(a.work, 'pzq_pit.tsv'), ['lo', 'hi', 'count', 'expected'], pit_rows)
    tsvio.write_table(os.path.join(a.work, 'pzq_coverage.tsv'), ['level', 'observed', 'lo', 'hi'], rep['pit']['coverage'])
    brow = []
    for nm in ('binned_zspec', 'binned_zphot', 'binned_mag'):
        for r in rep.get(nm, []):
            brow.append(dict(r, variable=nm[7:]))
    tsvio.write_table(os.path.join(a.work, 'pzq_binned.tsv'), ['variable', 'lo', 'hi', 'n', 'bias_mean', 'sigma_nmad', 'outlier_frac', 'absdz_p68', 'absdz_p95', 'pit_ks_p', 'cover68', 'pit_var'], brow)
    tsvio.write_table(os.path.join(a.work, 'pzq_stacked_nz.tsv'), ['z', 'stacked', 'spec'], [dict(z=z, stacked=s, spec=h) for z, s, h in zip(st['z'], st['stacked'], st['spec'])])
    with open(os.path.join(a.work, 'pzq_report.json'), 'w') as f:
        json.dump(clean(rep), f, indent=1)
    try:
        plot_quality(os.path.join(a.work, 'pzq_plot.png'), rep, pit, zsi, zpi, ok, st, None, 'photo-z quality (%s PDFs, %d spec-z)' % (kind, ok.sum()))
    except Exception as e:
        sys.stderr.write('plot failed: %s\n' % e)
    p = rep['point']
    sys.stderr.write('photo-z quality: N=%d bias %.4f sigma_NMAD %.4f outliers %.2f%% | PIT KS p=%.3g CvM p=%.3g: %s | CRPS %.4f\n' % (
        p['n'], p['bias_mean'], p['sigma_nmad'], 100 * p['outlier_frac'], rep['pit']['ks_p'], rep['pit']['cvm_p'], rep['pit']['shape'], rep['crps_mean']))
    # columns
    dz = ps.delta_z(zp, zs)
    full = {c: np.full(len(rows), np.nan) for c in QCOLS}
    full['PZ_PIT'][idx] = pit; full['PZ_CRPS'][idx] = cr; full['PZ_ZSCORE'][idx] = zsc; full['PZ_DZ'][idx] = dz[idx]
    full['PZ_INCI68'][idx] = (np.abs(pit - 0.5) <= 0.34).astype(float)
    full['PZ_OUTLIER'][idx] = (np.abs(dz[idx]) > a.outlier).astype(float)
    return rep, full


def do_consistency(cols, rows, a, W):
    """Two redshift estimates of the same objects (e.g. a template SED fit and the photo-z of this product): normalised difference, flags, accuracy of the agreeing set."""
    for c in (a.zphot_col, a.zerr_col, a.zalt_col, a.zalt_err_col):
        if c not in cols:
            sys.exit('ERROR: column %s not in the catalogue' % c)
    zA, sA, zB, sB = (colvals(rows, c) for c in (a.zphot_col, a.zerr_col, a.zalt_col, a.zalt_err_col))
    ok = np.isfinite(zA) & np.isfinite(zB) & np.isfinite(sA) & np.isfinite(sB) & (sA > 0) & (sB > 0)
    zs, zname = (find_zspec(cols, rows, a, W) if (a.zspec_col in cols or a.zspec_file) else (np.full(len(rows), np.nan), ''))
    has_spec = np.isfinite(zs) & (zs > 0)
    ref = ok & has_spec if (ok & has_spec).sum() >= 50 else ok
    cal = pcl.calibrate(zA[ref], np.minimum(sA[ref], a.err_cap), zB[ref], np.minimum(sB[ref], a.err_cap))
    sAc, sBc = np.minimum(sA, a.err_cap), np.minimum(sB, a.err_cap)
    idx = np.where(ok)[0]
    c = pcl.consistency(zA[idx], sAc[idx], zB[idx], sBc[idx], c=cal['c'], nsig=a.nsig, zspec=zs[idx] if has_spec[idx].all() else None)
    rep = dict(n_catalog=len(rows), n_compared=int(ok.sum()), calibration=cal, calibrated_on='spec-z objects' if ref is not ok else 'all compared objects', nsig=a.nsig, err_cap=a.err_cap,
               frac_flagged=c['frac_flagged'], A=a.zphot_col, B=a.zalt_col)
    if has_spec.any():
        hs = idx[has_spec[idx]]
        cs = pcl.consistency(zA[hs], sAc[hs], zB[hs], sBc[hs], c=cal['c'], nsig=a.nsig, zspec=zs[hs])
        rep['with_specz'] = {k: v for k, v in cs.items() if k not in ('flag', 'd')}
    os.makedirs(a.work, exist_ok=True)
    with open(os.path.join(a.work, 'pzq_consistency.json'), 'w') as fh:
        json.dump(clean(rep), fh, indent=1)
    full = {k: np.full(len(rows), np.nan) for k in CCOLS}
    full['PZ_ZDIFF'][idx] = (c['d'] - cal['median']) / max(cal['c'], 1e-9); full['PZ_INCONSIST'][idx] = c['flag'].astype(float)
    sys.stderr.write('photo-z consistency %s vs %s: %d compared, width c=%.2f, %.1f %% flagged (> %.1f sigma)\n' % (a.zphot_col, a.zalt_col, ok.sum(), cal['c'], 100 * c['frac_flagged'], a.nsig))
    return rep, full


def Predictive_subset(pred, idx):
    c = ps.Predictive.__new__(ps.Predictive)
    c.__dict__.update(pred.__dict__)
    for k in ('pi', 'mu', 'sig', 'med', 'slo', 'shi', 'p', 'cg'):
        if k in c.__dict__:
            v = c.__dict__[k]
            if isinstance(v, np.ndarray) and v.shape[:1] == (pred.n,):
                c.__dict__[k] = v[idx]
    c.n = len(idx)
    return c


def default_features(cols):
    mags = [c for c in cols if c.startswith('MAG_') and not c.startswith('MAGERR')]
    mags = [c for c in mags if c not in ('MAG_APER',)] if len(mags) > 1 else mags
    out = []
    if mags:
        out.append(mags[0])
        out += ['%s-%s' % (mags[i], mags[i + 1]) for i in range(len(mags) - 1)]
    return out


def feature_matrix(cols, rows, spec):
    names = []
    X = []
    for tok in spec:
        if '-' in tok and tok not in cols:
            a_, b_ = tok.split('-', 1)
            if a_ in cols and b_ in cols:
                X.append(colvals(rows, a_) - colvals(rows, b_)); names.append(tok)
        elif tok in cols:
            X.append(colvals(rows, tok)); names.append(tok)
    if not X:
        sys.exit('ERROR: none of the features %s found' % spec)
    X = np.array(X).T
    X[np.abs(X) > 90] = np.nan
    return X, names


def plot_repr(path, tbl, Xs, Xt, names, ratio, w, rep):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    nf = min(len(names), 4)
    fig, ax = plt.subplots(2, max(nf, 2), figsize=(3.4 * max(nf, 2), 6.2))
    for j in range(nf):
        a_, b_ = Xs[:, j][np.isfinite(Xs[:, j])], Xt[:, j][np.isfinite(Xt[:, j])]
        lo, hi = np.percentile(np.r_[a_, b_], [0.5, 99.5])
        bins = np.linspace(lo, hi, 30)
        ax[0, j].hist(b_, bins=bins, density=True, histtype='step', label='target')
        ax[0, j].hist(a_, bins=bins, density=True, histtype='step', label='spec')
        ax[0, j].hist(a_, bins=bins, density=True, histtype='step', weights=w[np.isfinite(Xs[:, j])] if w is not None else None, label='spec weighted', ls='--')
        ax[0, j].set_title('%s  KS=%.2f' % (names[j], tbl[j]['ks']), fontsize=8)
        if j == 0:
            ax[0, j].legend(fontsize=7)
    ax[1, 0].hist(np.clip(ratio, 0, 5), bins=40)
    ax[1, 0].axvline(1, color='r'); ax[1, 0].set_xlabel('distance to nearest spec / p%g of spec NN' % rep['coverage']['q']); ax[1, 0].set_title('outside support: %.1f%% (expected %.0f%%)' % (
        100 * rep['coverage']['frac_outside'], 100 * rep['coverage']['expected']), fontsize=8)
    ax[1, 1].hist(w, bins=40); ax[1, 1].set_xlabel('kNN weight (mean 1)'); ax[1, 1].set_title('ESS fraction %.2f' % rep['weights']['ess_fraction'], fontsize=8)
    for j in range(2, max(nf, 2)):
        ax[1, j].axis('off')
    for j in range(nf, max(nf, 2)):
        ax[0, j].axis('off')
    fig.suptitle('spec sample (%d) vs target (%d): kNN AUC %.3f (p=%.3g)' % (rep['n_spec'], rep['n_target'], rep['knn_test']['auc'], rep['knn_test']['p']), fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def do_repr(cols, rows, a, W):
    zs, zname = find_zspec(cols, rows, a, W)
    feats = [t.strip() for t in a.features.split(',') if t.strip()] or default_features(cols)
    X, names = feature_matrix(cols, rows, feats)
    okx = np.all(np.isfinite(X), axis=1)
    is_spec = np.isfinite(zs) & (zs > 0) & okx
    tgt = okx & (~is_spec if a.target == 'nospec' else True)
    if is_spec.sum() < 30 or tgt.sum() < 30:
        sys.exit('ERROR: need >= 30 spec objects (%d) and >= 30 target objects (%d) with all features' % (is_spec.sum(), tgt.sum()))
    Z = ps.robust_scale(X, X[tgt])
    Xs, Xt = X[is_spec], X[tgt]
    Zs, Zt = Z[is_spec], Z[tgt]
    rep = dict(n_catalog=len(rows), n_spec=int(is_spec.sum()), n_target=int(tgt.sum()), features=names, target=a.target, zspec_column=zname)
    tbl = ps.feature_table(Xs, Xt, names)
    rep['features_1d'] = tbl
    rep['knn_test'] = ps.knn_classifier_test(Zs, Zt, k=a.knn_test_k, nperm=a.nperm, seed=1)
    ratio_t, fo, ex = ps.coverage_distance(Zs, Zt, a.coverage_q)
    rep['coverage'] = dict(q=a.coverage_q, frac_outside=fo, expected=ex)
    w = ps.knn_weights(Zs, Zt, k=a.knn)
    rep['weights'] = dict(k=a.knn, ess_fraction=ps.effective_fraction(w), w_min=float(w.min()), w_max=float(w.max()), w_p99=float(np.percentile(w, 99)))
    # weighted feature moments (does the reweighted spec sample reproduce the target?)
    wm = []
    for j, nme in enumerate(names):
        wm.append(dict(feature=nme, mean_target=float(Xt[:, j].mean()), mean_spec=float(Xs[:, j].mean()), mean_spec_weighted=float(np.average(Xs[:, j], weights=w)),
                       ks_weighted=weighted_ks(Xs[:, j], w, Xt[:, j])))
    rep['weighted_moments'] = wm
    rep['verdict'] = ('representative' if rep['knn_test']['p'] > 0.01 and abs(rep['knn_test']['auc'] - 0.5) < 0.03 and fo < 2 * ex + 0.02 else
                      'NOT representative: use weights / restrict to the covered region')
    # accuracy check when a photo-z and spec-z exist
    zpc = colvals(rows, a.zphot_col) if a.zphot_col in cols else None
    if zpc is not None:
        s = is_spec & np.isfinite(zpc)
        if s.sum() > 20:
            sub_w = np.interp(np.where(s)[0], np.where(is_spec)[0], w)
            rep['accuracy'] = dict(unweighted=ps.point_metrics(zpc[s], zs[s], a.outlier), weighted=ps.point_metrics(zpc[s], zs[s], a.outlier, weights=sub_w))
    os.makedirs(a.work, exist_ok=True)
    tsvio.write_table(os.path.join(a.work, 'pzr_features.tsv'), ['feature', 'n_spec', 'n_target', 'ks', 'ks_p', 'std_mean_diff', 'median_spec', 'median_target', 'p5_spec', 'p95_spec',
                      'p5_target', 'p95_target', 'frac_target_outside_spec_range'], tbl)
    with open(os.path.join(a.work, 'pzr_report.json'), 'w') as f:
        json.dump(clean(rep), f, indent=1)
    try:
        plot_repr(os.path.join(a.work, 'pzr_plot.png'), tbl, Xs, Xt, names, ratio_t, w, rep)
    except Exception as e:
        sys.stderr.write('plot failed: %s\n' % e)
    sys.stderr.write('representativeness: %d spec vs %d target, kNN AUC %.3f (p=%.3g), outside support %.1f%% (exp %.0f%%), ESS %.2f -> %s\n' % (
        rep['n_spec'], rep['n_target'], rep['knn_test']['auc'], rep['knn_test']['p'], 100 * fo, 100 * ex, rep['weights']['ess_fraction'], rep['verdict']))
    full = {c: np.full(len(rows), np.nan) for c in RCOLS}
    from scipy.spatial import cKDTree
    d_all = cKDTree(Zs).query(Z[okx], k=1)[0]
    ref = np.percentile(cKDTree(Zs).query(Zs, k=2)[0][:, 1], a.coverage_q)
    full['PZ_SPECDIST'][okx] = d_all / ref
    full['PZ_INSPEC'][okx] = (d_all / ref <= 1).astype(float)
    full['PZ_SPECW'][np.where(is_spec)[0]] = w
    return rep, full


def weighted_ks(x, w, y):
    """KS distance between the weighted sample (x, w) and the sample y."""
    x = np.asarray(x); o = np.argsort(x); xs = x[o]; cw = np.cumsum(np.asarray(w)[o]); cw = cw / cw[-1]
    ys = np.sort(y)
    g = np.sort(np.r_[xs, ys])
    fx = np.interp(g, xs, cw, left=0.0)
    fy = np.searchsorted(ys, g, side='right') / len(ys)
    return float(np.max(np.abs(fx - fy)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--work', required=True)
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--mode', choices=['quality', 'repr', 'consistency'], default='quality')
    ap.add_argument('--recal-out', default='', help='quality: save the PIT recalibration map learnt on the spec-z objects (JSON)')
    ap.add_argument('--recal-in', default='', help='quality: apply a saved map; writes pzq_recalibrated.tsv (recalibrated 16/50/84 percentiles) and evaluates it')
    ap.add_argument('--zalt-col', default='EZ_Z', help='consistency: second redshift estimate (e.g. SED-fit z)')
    ap.add_argument('--zalt-err-col', default='EZ_ZERR')
    ap.add_argument('--nsig', type=float, default=3.0)
    ap.add_argument('--err-cap', type=float, default=0.3, help='consistency: cap on the quoted redshift errors (multimodal PDFs give huge widths)')
    ap.add_argument('--zspec-col', default='Z_SPEC')
    ap.add_argument('--zspec-file', default='')
    ap.add_argument('--zphot-col', default='PHOTO_Z')
    ap.add_argument('--zerr-col', default='PHOTO_Z_ERR')
    ap.add_argument('--q16-col', default='')
    ap.add_argument('--q84-col', default='')
    ap.add_argument('--mixture-file', default='')
    ap.add_argument('--pdf-file', default='')
    ap.add_argument('--pdf-kind', choices=['auto', 'gauss', 'split', 'mixture', 'grid'], default='auto')
    ap.add_argument('--zmin', type=float, default=0.0)
    ap.add_argument('--outlier', type=float, default=0.15)
    ap.add_argument('--mag-col', default='MAG_AUTO')
    ap.add_argument('--mag-bins', default='')
    ap.add_argument('--z-bins', default='')
    ap.add_argument('--weights-col', default='')
    ap.add_argument('--nboot', type=int, default=200)
    ap.add_argument('--features', default='')
    ap.add_argument('--target', choices=['all', 'nospec'], default='all')
    ap.add_argument('--knn', type=int, default=5)
    ap.add_argument('--knn-test-k', type=int, default=20)
    ap.add_argument('--nperm', type=int, default=200)
    ap.add_argument('--coverage-q', type=float, default=95.0)
    a = ap.parse_args(argv)
    cols, rows = tsvio.read_catalog(a.catalog)
    if not rows or 'NUMBER' not in cols:
        sys.exit('ERROR: catalog is empty or has no NUMBER column')
    rep, full = {'quality': do_quality, 'repr': do_repr, 'consistency': do_consistency}[a.mode](cols, rows, a, None)
    out_cols = {'quality': QCOLS, 'repr': RCOLS, 'consistency': CCOLS}[a.mode]
    out = [(r['NUMBER'], {c: (float(full[c][i]) if np.isfinite(full[c][i]) else None) for c in out_cols}) for i, r in enumerate(rows)]
    tsvio.write_columns(out_cols, out)
    if a.meta_out:
        if a.mode == 'quality':
            p = rep['point']
            ometa.update(a.meta_out, 'photoz_quality', clean(dict(n_spec=rep['n_spec'], sigma_nmad=p['sigma_nmad'], bias=p['bias_mean'], outlier_frac=p['outlier_frac'],
                         pit_ks_p=rep['pit']['ks_p'], pit_cvm_p=rep['pit']['cvm_p'], pit_shape=rep['pit']['shape'], crps=rep['crps_mean'], pdf_kind=rep['pdf_kind'])), nrows=len(rows))
        elif a.mode == 'consistency':
            ometa.update(a.meta_out, 'photoz_consistency', clean(dict(n=rep['n_compared'], width=rep['calibration']['c'], frac_flagged=rep['frac_flagged'])), nrows=len(rows))
        else:
            ometa.update(a.meta_out, 'photoz_repr', clean(dict(n_spec=rep['n_spec'], n_target=rep['n_target'], knn_auc=rep['knn_test']['auc'], knn_p=rep['knn_test']['p'],
                         frac_outside=rep['coverage']['frac_outside'], ess_fraction=rep['weights']['ess_fraction'], verdict=rep['verdict'])), nrows=len(rows))
    return 0


if __name__ == '__main__':
    sys.exit(main())
