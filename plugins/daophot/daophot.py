#!/usr/bin/env python3
"""DAOPHOT / ALLSTAR / DAOMASTER-style crowded-field photometry (CLI of the `daophot` plugin; engine: ogfkit/daophot.py, ogfkit/psfmodel.py).

    daophot.py IMAGE --mode MODE --work DIR [--catalog CAT] [--mask M] [options]

MODE (each reads/writes fixed file names in --work, so the steps can be run one after the other, in the GUI or in a replayed session):
    find      FIND     -> daophot_find.tsv            (threshold, sharpness / roundness cuts)
    phot      PHOT     -> daophot_phot.tsv            (3 apertures, sky annulus, sky mode)
    pickpsf   PICKPSF  -> daophot_psfstars.tsv        (auto selection + --psf-add / --psf-remove)
    psf       PSF      -> daophot_psf.json/.fits      (analytic + lookup table, order 0-3, neighbour subtraction)
    fit       ALLSTAR  -> daophot_stars.tsv, daophot_resid.fits, daophot_result.json, daophot_diag.png   (needs the PSF: built if missing)
    run       all of the above in one go
    substar   SUBSTAR  -> daophot_sub.fits            (--keep-ids stay in the image)
    addstar   ADDSTAR  -> daophot_add.fits + daophot_add_truth.tsv
    apcorr    growth curve / aperture correction from the PSF stars -> printed + daophot_result.json
    match     DAOMATCH/DAOMASTER with --second IMAGE2 (own pipeline run on it) or --second-stars TSV -> daophot_master.tsv, daophot_cmd.tsv, daophot_cmd.png
    arttest   artificial-star completeness test (plugins/completeness engine) with the DAOPHOT fitter as the detection callback

With --catalog the step prints `NUMBER` + DAO_* columns (nearest fitted star within --match-radius) for the add_columns contract.
Coordinates in all tables are 1-based like the catalog.
"""
import argparse
import json
import math
import os
import sys
import warnings

import numpy as np

_here = os.path.dirname(os.path.abspath(__file__))
_root = os.path.abspath(os.path.join(_here, '..', '..'))
for _p in (_root,):
    if _p not in sys.path:
        sys.path.insert(0, _p)
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio, psfmodel as pm, daophot as dp  # noqa: E402

COLUMNS = ['DAO_FLUX', 'DAO_FLUXERR', 'DAO_MAG', 'DAO_MAGERR', 'DAO_MAG_APC', 'DAO_X', 'DAO_Y', 'DAO_CHI', 'DAO_SHARP', 'DAO_ROUND', 'DAO_SNR',
           'DAO_SKY', 'DAO_NGRP', 'DAO_FLAG', 'DAO_GOOD']
STAR_COLS = ['ID', 'X_IMAGE', 'Y_IMAGE', 'FLUX', 'FLUXERR', 'MAG', 'MAGERR', 'MAG_APC', 'MAG_AP1', 'MAGERR_AP1', 'MAG_AP2', 'MAGERR_AP2', 'MAG_AP3', 'MAGERR_AP3',
             'SKY', 'CHI', 'SHARP', 'ROUND', 'SNR', 'NGRP', 'NIT', 'FOUND_ITER', 'FLAG', 'GOOD', 'XERR', 'YERR']
MATCH_COLUMNS = ['DAO_MATCH', 'DAO_MX', 'DAO_MY', 'DAO_MSEP', 'DAO_MMAG2', 'DAO_COLOR']


def W(work, name):
    return os.path.join(work, 'daophot_' + name)


def log(*a):
    sys.stderr.write(' '.join(str(x) for x in a) + '\n')


def stars_to_rows(stars):
    rows = []
    for s in stars:
        r = dict(ID=s['id'], X_IMAGE=s['x'] + 1, Y_IMAGE=s['y'] + 1, FLUX=s['flux'], FLUXERR=s['fluxerr'], MAG=s['mag'], MAGERR=s['magerr'], MAG_APC=s.get('mag_apc'),
                 SKY=s.get('sky_ap', s.get('sky')), CHI=s['chi'], SHARP=s['sharp'], ROUND=s['round'], SNR=s['snr'], NGRP=s['gsize'], NIT=s['nit'],
                 FOUND_ITER=s['found_iter'], FLAG=s['flag'], GOOD=int(s['good']), XERR=s.get('xerr'), YERR=s.get('yerr'))
        for k in range(3):
            am = s.get('ap_mag', [None] * 3)
            ae = s.get('ap_magerr', [None] * 3)
            r['MAG_AP%d' % (k + 1)] = am[k] if k < len(am) else None
            r['MAGERR_AP%d' % (k + 1)] = ae[k] if k < len(ae) else None
        rows.append(r)
    return rows


def read_stars(path):
    cols, rows = tsvio.read_catalog(path)
    f = tsvio.fnum
    out = []
    for r in rows:
        out.append(dict(id=int(f(r['ID'])), x=f(r['X_IMAGE']) - 1, y=f(r['Y_IMAGE']) - 1, flux=f(r['FLUX']), fluxerr=f(r['FLUXERR']), mag=f(r['MAG']),
                        magerr=f(r['MAGERR']), mag_apc=f(r.get('MAG_APC')), chi=f(r.get('CHI')), sharp=f(r.get('SHARP')), round=f(r.get('ROUND')), snr=f(r.get('SNR')),
                        gsize=int(f(r.get('NGRP'), 1)), nit=int(f(r.get('NIT'), 0)), found_iter=int(f(r.get('FOUND_ITER'), 0)), flag=int(f(r.get('FLAG'), 0)),
                        good=bool(f(r.get('GOOD'), 1)), sky=f(r.get('SKY'), 0.0)))
    return out


def map_catalog(cat, stars, radius, build):
    """add_columns rows: for every catalog object the nearest star within `radius` (1-based coordinates in both)."""
    from scipy.spatial import cKDTree
    cols, rows = tsvio.read_catalog(cat)
    cx = np.array([tsvio.fnum(r.get('X_IMAGE')) for r in rows])
    cy = np.array([tsvio.fnum(r.get('Y_IMAGE')) for r in rows])
    out = []
    if not stars:
        return [(r['NUMBER'], {}) for r in rows]
    tr = cKDTree(np.array([[s['x'] + 1, s['y'] + 1] for s in stars]))
    d, j = tr.query(np.c_[np.nan_to_num(cx, nan=-1e6), np.nan_to_num(cy, nan=-1e6)])
    for r, dd, jj in zip(rows, d, j):
        out.append((r['NUMBER'], build(stars[jj], dd) if dd <= radius else {}))
    return out


def star_cols(s, d):
    return dict(DAO_FLUX=s['flux'], DAO_FLUXERR=s['fluxerr'], DAO_MAG=s['mag'], DAO_MAGERR=s['magerr'], DAO_MAG_APC=s.get('mag_apc'), DAO_X=s['x'] + 1, DAO_Y=s['y'] + 1,
                DAO_CHI=s['chi'], DAO_SHARP=s['sharp'], DAO_ROUND=s['round'], DAO_SNR=s['snr'], DAO_SKY=s.get('sky_ap', s.get('sky')), DAO_NGRP=s['gsize'],
                DAO_FLAG=s['flag'], DAO_GOOD=int(s['good']))


def plot_diag(res, path, stars=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    S = res['stars']
    mag = np.array([s['mag'] for s in S], float)
    good = np.array([s['good'] for s in S], bool)
    fig, ax = plt.subplots(2, 3, figsize=(11, 6))
    for a, key, lab in ((ax[0, 0], 'chi', 'chi'), (ax[0, 1], 'sharp', 'sharp'), (ax[0, 2], 'magerr', 'mag error')):
        v = np.array([s[key] for s in S], float)
        a.plot(mag[~good], v[~good], 'r.', ms=3)
        a.plot(mag[good], v[good], 'k.', ms=3)
        a.set_xlabel('PSF mag')
        a.set_ylabel(lab)
        a.grid(alpha=.3)
    ax[0, 0].set_ylim(0, 5)
    ax[0, 2].set_ylim(0, 0.5)
    ax[0, 1].set_ylim(-1.5, 1.5)
    d = res['psf_diag']
    for a, key, lab in ((ax[1, 0], 'fwhm', 'PSF FWHM map (px)'), (ax[1, 1], 'e', 'PSF ellipticity map')):
        im = a.imshow(np.array(d[key]), origin='lower', extent=[d['x'][0], d['x'][-1], d['y'][0], d['y'][-1]])
        fig.colorbar(im, ax=a, fraction=0.046)
        a.set_title(lab, fontsize=9)
    ac = res.get('apcorr', {})
    if ac.get('radii'):
        ax[1, 2].errorbar(ac['radii'], ac['corr'], yerr=np.nan_to_num(ac.get('err', [0] * len(ac['radii']))), fmt='ko-', label='stars')
        ax[1, 2].plot(ac['radii'], ac.get('model_corr', []), 'b--', label='PSF model')
        ax[1, 2].set_xlabel('aperture radius (px)')
        ax[1, 2].set_ylabel('m_ap - m_psf')
        ax[1, 2].legend(fontsize=8)
        ax[1, 2].grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def plot_cmd(rows, path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    c = np.array([r['COLOR'] for r in rows], float)
    m = np.array([r['MAG1'] for r in rows], float)
    ax.plot(c, m, 'k.', ms=4)
    ax.invert_yaxis()
    ax.set_xlabel('mag1 - mag2')
    ax.set_ylabel('mag1')
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def pipeline_kwargs(a):
    sat = a.saturation if a.saturation > 0 else None
    radii = [float(v) for v in a.apertures.replace(' ', '').split(',') if v][:3]
    while len(radii) < 3:
        radii.append(radii[-1] + 2.0 if radii else 5.0)
    return dict(zp=a.mag_zeropoint, fwhm=a.fwhm, find_thresh=a.find_thresh, find_sharp=(a.sharp_lo, a.sharp_hi), find_round=(a.round_lo, a.round_hi), radii=radii,
                sky_inner=a.sky_inner, sky_outer=a.sky_outer, sky_mode=a.sky_mode, psf_kind=a.psf_kind, psf_order=a.psf_order, psf_lookup=not a.no_lookup,
                psf_size=a.psf_size, psf_nstars=a.psf_nstars, psf_snr=a.psf_snr, saturation=sat, psf_add=dp._parse_xy_list(a.psf_add), psf_remove=dp._parse_xy_list(a.psf_remove),
                neighbour_iter=a.neighbour_iter, fit_radius=a.fit_radius, fit_sky=not a.no_fit_sky, gain=a.gain or None, n_iter=a.n_iter, chi_max=a.chi_max,
                snr_min=a.snr_min, sharp_cut=(a.fit_sharp_lo, a.fit_sharp_hi), round_cut=(a.fit_round_lo, a.fit_round_hi), max_group=a.max_group,
                apcorr_radius=a.apcorr_radius, n_workers=a.n_workers)


def write_find(path, f):
    rows = [dict(ID=i + 1, X_IMAGE=f['x'][i] + 1, Y_IMAGE=f['y'][i] + 1, PEAK=f['peak'][i], FLUX=f['flux'][i], SNR=f.get('snr', np.zeros(len(f['x'])))[i], SHARP=f['sharp'][i],
                 ROUND1=f['round1'][i], ROUND2=f['round2'][i]) for i in range(len(f['x']))]
    tsvio.write_table(path, ['ID', 'X_IMAGE', 'Y_IMAGE', 'PEAK', 'FLUX', 'SNR', 'SHARP', 'ROUND1', 'ROUND2'], rows)


def read_find(path):
    cols, rows = tsvio.read_catalog(path)
    g = lambda k: np.array([tsvio.fnum(r[k]) for r in rows], float)
    return dict(x=g('X_IMAGE') - 1, y=g('Y_IMAGE') - 1, peak=g('PEAK'), flux=g('FLUX'), snr=g('SNR'), sharp=g('SHARP'), round1=g('ROUND1'), round2=g('ROUND2'), npix=np.zeros(len(rows)))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--mode', default='run')
    ap.add_argument('--work', default='.')
    ap.add_argument('--catalog', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--n-workers', type=int, default=1)
    ap.add_argument('--match-radius', type=float, default=2.0)
    ap.add_argument('--fwhm', type=float, default=3.0)
    ap.add_argument('--find-thresh', type=float, default=4.0)
    ap.add_argument('--sharp-lo', type=float, default=0.2)
    ap.add_argument('--sharp-hi', type=float, default=1.0)
    ap.add_argument('--round-lo', type=float, default=-1.0)
    ap.add_argument('--round-hi', type=float, default=1.0)
    ap.add_argument('--apertures', default='3,5,8')
    ap.add_argument('--sky-inner', type=float, default=12.0)
    ap.add_argument('--sky-outer', type=float, default=18.0)
    ap.add_argument('--sky-mode', default='median', choices=list(dp.SKY_MODES))
    ap.add_argument('--psf-kind', default='empirical', choices=['empirical'] + list(pm.ANALYTIC_KINDS))
    ap.add_argument('--psf-order', type=int, default=2)
    ap.add_argument('--no-lookup', action='store_true')
    ap.add_argument('--psf-size', type=int, default=0)
    ap.add_argument('--psf-nstars', type=int, default=60)
    ap.add_argument('--psf-snr', type=float, default=20.0)
    ap.add_argument('--saturation', type=float, default=0.0)
    ap.add_argument('--psf-add', default='')
    ap.add_argument('--psf-remove', default='')
    ap.add_argument('--neighbour-iter', type=int, default=2)
    ap.add_argument('--fit-radius', type=float, default=0.0)
    ap.add_argument('--no-fit-sky', action='store_true')
    ap.add_argument('--gain', type=float, default=0.0)
    ap.add_argument('--n-iter', type=int, default=3)
    ap.add_argument('--chi-max', type=float, default=3.0)
    ap.add_argument('--snr-min', type=float, default=3.0)
    ap.add_argument('--fit-sharp-lo', type=float, default=-1.0)
    ap.add_argument('--fit-sharp-hi', type=float, default=1.0)
    ap.add_argument('--fit-round-lo', type=float, default=-1.0)
    ap.add_argument('--fit-round-hi', type=float, default=1.0)
    ap.add_argument('--max-group', type=int, default=20)
    ap.add_argument('--apcorr-radius', type=float, default=10.0)
    ap.add_argument('--keep-ids', default='')
    ap.add_argument('--add-n', type=int, default=100)
    ap.add_argument('--add-mag-min', type=float, default=20.0)
    ap.add_argument('--add-mag-max', type=float, default=24.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--second', default='')
    ap.add_argument('--second-stars', default='')
    ap.add_argument('--match-model', default='affine', choices=['shift', 'similarity', 'affine', 'poly2', 'poly3'])
    ap.add_argument('--art-mag-min', type=float, default=19.0)
    ap.add_argument('--art-mag-max', type=float, default=24.0)
    ap.add_argument('--art-bins', type=int, default=6)
    ap.add_argument('--art-per-bin', type=int, default=12)
    ap.add_argument('--art-per-image', type=int, default=10)
    ap.add_argument('--crop-size', type=int, default=0)
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    data, hdr = imageio.load_image(a.image)
    mask = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    kw = pipeline_kwargs(a)
    mode = a.mode
    emit = None                     # rows for the add_columns contract
    text = []
    f_find = W(a.work, 'find.tsv')

    def get_find():
        if os.path.isfile(f_find):
            return read_find(f_find)
        bkg, rms = pm.background(data, mask)
        f = dp.find(data, bkg, rms, fwhm=a.fwhm, thresh=a.find_thresh, sharp=kw['find_sharp'], round_=kw['find_round'], mask=mask)
        ix = np.clip(np.round(f['y']).astype(int), 0, data.shape[0] - 1)
        jx = np.clip(np.round(f['x']).astype(int), 0, data.shape[1] - 1)
        f['snr'] = f['peak'] / rms[ix, jx] if len(ix) else np.zeros(0)
        write_find(f_find, f)
        return f

    def load_psf():
        p = W(a.work, 'psf.json')
        if not os.path.isfile(p):
            raise SystemExit('no PSF model yet (%s) - run the PSF step first' % p)
        return pm.load_model(p)

    if mode == 'find':
        bkg, rms = pm.background(data, mask)
        f = dp.find(data, bkg, rms, fwhm=a.fwhm, thresh=a.find_thresh, sharp=kw['find_sharp'], round_=kw['find_round'], mask=mask)
        ix = np.clip(np.round(f['y']).astype(int), 0, data.shape[0] - 1)
        jx = np.clip(np.round(f['x']).astype(int), 0, data.shape[1] - 1)
        f['snr'] = f['peak'] / rms[ix, jx] if len(ix) else np.zeros(0)
        write_find(f_find, f)
        text.append('FIND: %d candidates (threshold %.1f sigma, sharpness %.2f-%.2f, roundness %.2f-%.2f) -> %s' % (len(f['x']), a.find_thresh, a.sharp_lo, a.sharp_hi, a.round_lo, a.round_hi, f_find))
    elif mode == 'phot':
        f = get_find()
        bkg, rms = pm.background(data, mask)
        ph = dp.phot(data, np.c_[f['x'], f['y']], radii=kw['radii'], sky_inner=a.sky_inner, sky_outer=a.sky_outer, sky_mode=a.sky_mode, bkg=bkg, rms=rms, zp=a.mag_zeropoint, gain=a.gain or None)
        rows = []
        for i, p in enumerate(ph):
            r = dict(ID=i + 1, X_IMAGE=f['x'][i] + 1, Y_IMAGE=f['y'][i] + 1, SKY=p['sky'], SKYSIG=p['skysig'], NSKY=p['nsky'])
            for k in range(len(kw['radii'])):
                r['MAG_AP%d' % (k + 1)], r['MAGERR_AP%d' % (k + 1)] = p['mag'][k], p['magerr'][k]
            rows.append(r)
        cols = ['ID', 'X_IMAGE', 'Y_IMAGE', 'SKY', 'SKYSIG', 'NSKY'] + sum([['MAG_AP%d' % (k + 1), 'MAGERR_AP%d' % (k + 1)] for k in range(len(kw['radii']))], [])
        tsvio.write_table(W(a.work, 'phot.tsv'), cols, rows)
        text.append('PHOT: %d stars, apertures %s px, sky annulus %.0f-%.0f px (%s) -> %s' % (len(rows), kw['radii'], a.sky_inner, a.sky_outer, a.sky_mode, W(a.work, 'phot.tsv')))
    elif mode == 'pickpsf':
        f = get_find()
        half = int(max(8, math.ceil(2.6 * a.fwhm)))
        xy, info = dp.pick_psf(f, data.shape, n_max=a.psf_nstars, half=half, snr_min=a.psf_snr, saturation=kw['saturation'], add=kw['psf_add'], remove=kw['psf_remove'])
        tsvio.write_table(W(a.work, 'psfstars.tsv'), ['ID', 'X_IMAGE', 'Y_IMAGE'], [dict(ID=i + 1, X_IMAGE=p[0] + 1, Y_IMAGE=p[1] + 1) for i, p in enumerate(xy)])
        text.append('PICKPSF: %d PSF stars (auto %d, added %d, removed %d); rejected: %s -> %s' % (len(xy), info['n_auto'], info['n_added'], info['n_removed'], info['reasons'], W(a.work, 'psfstars.tsv')))
    elif mode == 'psf':
        f = get_find()
        bkg, rms = pm.background(data, mask)
        half = int(max(8, math.ceil(2.6 * a.fwhm)))
        xy, info = dp.pick_psf(f, data.shape, n_max=a.psf_nstars, half=half, snr_min=a.psf_snr, saturation=kw['saturation'], add=kw['psf_add'], remove=kw['psf_remove'])
        tsvio.write_table(W(a.work, 'psfstars.tsv'), ['ID', 'X_IMAGE', 'Y_IMAGE'], [dict(ID=i + 1, X_IMAGE=p[0] + 1, Y_IMAGE=p[1] + 1) for i, p in enumerate(xy)])

        def nb_fit(mdl, _s):
            r, _ = dp.allstar(data, mdl, bkg=bkg, rms=rms, mask=mask, init_stars=f, fwhm=a.fwhm, thresh=a.find_thresh, n_iter=1, refind=False, fit_sky=False, zp=a.mag_zeropoint)
            s = [q for q in r['stars'] if q['flux'] > 0]
            return [q['x'] for q in s], [q['y'] for q in s], [q['flux'] for q in s]

        model, pinfo = pm.build_psf_model(data, bkg=bkg, rms=rms, mask=mask, xy=xy if len(xy) else None, fwhm_prior=a.fwhm, size=a.psf_size or None, degree=a.psf_order,
                                          saturation=kw['saturation'], snr_min=a.psf_snr, neighbour_fit=nb_fit if (a.neighbour_iter > 0 and len(xy) >= 4) else None,
                                          n_iter=max(1, a.neighbour_iter + 1), kind=a.psf_kind, lookup=kw['psf_lookup'])
        model.save(W(a.work, 'psf.json'))
        model.save(W(a.work, 'psf.fits'))
        dg = pm.diagnostics(model, 5)
        with open(W(a.work, 'psf_diag.json'), 'w') as fh:
            json.dump(dict(info={k: v for k, v in pinfo.items() if k != 'stars'}, diag=dg, stars=pinfo.get('stars', [])), fh, default=lambda o: None)
        text.append('PSF: mode %s, %s stars used (of %d picked), polynomial order %s, oversampling %s, stamp %d px' % (model.meta.get('mode'), pinfo.get('n_used', pinfo.get('n_stars')), len(xy), model.degree, model.oversample, model.size))
        text.append('PSF FWHM %.2f-%.2f px (mean %.2f), ellipticity %.3f-%.3f over the field -> %s' % (dg['fwhm_min'], dg['fwhm_max'], dg['fwhm_mean'], dg['e_min'], dg['e_max'], W(a.work, 'psf.json')))
    elif mode in ('fit', 'run'):
        pm_model = None
        if mode == 'fit':
            pm_model = load_psf()
        res, resid, model = dp.run_pipeline(data, mask=mask, psf_model=pm_model, log=log, **kw)
        model.save(W(a.work, 'psf.json'))
        model.save(W(a.work, 'psf.fits'))
        stars = res['stars']
        tsvio.write_table(W(a.work, 'stars.tsv'), STAR_COLS, stars_to_rows(stars))
        imageio.save_fits(W(a.work, 'resid.fits'), resid, header=hdr)
        out = {k: v for k, v in res.items() if k not in ('stars', 'psf_stars')}
        out['image'] = os.path.basename(a.image)
        with open(W(a.work, 'result.json'), 'w') as fh:
            json.dump(out, fh, default=lambda o: None)
        plot_diag(res, W(a.work, 'diag.png'))
        if mode == 'run':
            f = None
        text.append('ALLSTAR: %d stars fitted (%d pass the cuts chi<%.1f, S/N>%.1f, sharp %.2f..%.2f); new stars per iteration %s; fit radius %.1f px' % (res['n_stars'], res['n_good'], a.chi_max, a.snr_min, a.fit_sharp_lo, a.fit_sharp_hi, res['n_found_per_iter'], res['fit_radius']))
        text.append('aperture correction (PSF mag -> %.0f px aperture) %+.4f mag from %d PSF stars; residual rms %.3f (sky rms %.3f)' % (a.apcorr_radius, res['apcorr_value'], res['apcorr']['n'], float(np.nanstd(resid)), res['rms_median']))
        if a.catalog:
            emit = map_catalog(a.catalog, stars, a.match_radius, star_cols)
    elif mode == 'substar':
        model = load_psf()
        stars = read_stars(W(a.work, 'stars.tsv'))
        keep = tsvio.parse_numbers(a.keep_ids)
        out = dp.substar(data, model, stars, keep=keep)
        imageio.save_fits(W(a.work, 'sub.fits'), out, header=hdr)
        text.append('SUBSTAR: subtracted %d of %d fitted stars (kept ids %s); rms before %.3f after %.3f -> %s' % (len(stars) - len(set(keep) & {s['id'] for s in stars}), len(stars), keep or '-', float(np.nanstd(data)), float(np.nanstd(out)), W(a.work, 'sub.fits')))
    elif mode == 'addstar':
        model = load_psf()
        pos = dp.random_positions(data.shape, a.add_n, margin=model.size // 2 + 2, seed=a.seed)
        mags = np.random.default_rng(a.seed + 7).uniform(a.add_mag_min, a.add_mag_max, len(pos))
        out, truth = dp.addstar(data, model, pos, mags, zp=a.mag_zeropoint, gain=a.gain or None, seed=a.seed)
        imageio.save_fits(W(a.work, 'add.fits'), out, header=hdr)
        tsvio.write_table(W(a.work, 'add_truth.tsv'), ['ID', 'X_IMAGE', 'Y_IMAGE', 'MAG', 'FLUX'], [dict(ID=i + 1, X_IMAGE=t['x'] + 1, Y_IMAGE=t['y'] + 1, MAG=t['mag'], FLUX=t['flux']) for i, t in enumerate(truth)])
        text.append('ADDSTAR: %d artificial stars (mag %.1f-%.1f) added -> %s (truth: %s)' % (len(truth), a.add_mag_min, a.add_mag_max, W(a.work, 'add.fits'), W(a.work, 'add_truth.tsv')))
    elif mode == 'apcorr':
        model = load_psf()
        stars = read_stars(W(a.work, 'stars.tsv'))
        bkg, rms = pm.background(data, mask)
        sub = np.nan_to_num(data - bkg, nan=0.0)
        px = W(a.work, 'psfstars.tsv')
        if os.path.isfile(px):
            _, rr = tsvio.read_catalog(px)
            xy = np.array([[tsvio.fnum(r['X_IMAGE']) - 1, tsvio.fnum(r['Y_IMAGE']) - 1] for r in rr])
        else:
            xy = np.array([(s['x'], s['y']) for s in stars if s['good'] and s['snr'] > 30][:40])
        rad = sorted(set([3, 4, 5, 6, 8, 10, 12, 15, a.apcorr_radius]))
        ac = dp.aperture_correction(sub, model, stars, xy, radii=[r for r in rad if r <= model.size // 2 + 4])
        text.append('APERTURE CORRECTION from %d PSF stars (m_ap - m_psf):' % ac['n'])
        for r, c, e, mc in zip(ac['radii'], ac['corr'], ac['err'], ac['model_corr']):
            text.append('  r=%4.1f px  stars %+.4f +- %.4f   PSF model %+.4f' % (r, c, e if e == e else 0, mc))
        try:
            rj = json.load(open(W(a.work, 'result.json')))
        except (OSError, ValueError):
            rj = {}
        rj['apcorr'] = ac
        json.dump(rj, open(W(a.work, 'result.json'), 'w'), default=lambda o: None)
    elif mode == 'match':
        s1 = read_stars(W(a.work, 'stars.tsv'))
        if a.second_stars:
            s2 = read_stars(a.second_stars)
        elif a.second:
            d2, h2 = imageio.load_image(a.second)
            kw2 = dict(kw)
            kw2.update(n_iter=2, neighbour_iter=1, psf_add=[], psf_remove=[])
            r2, _, _ = dp.run_pipeline(d2, log=log, **kw2)
            s2 = r2['stars']
            tsvio.write_table(W(a.work, 'stars2.tsv'), STAR_COLS, stars_to_rows(s2))
        else:
            raise SystemExit('match needs --second IMAGE or --second-stars TSV')
        g1 = [s for s in s1 if s['good']]
        g2 = [s for s in s2 if s['good']]
        A = dict(x=[s['x'] for s in g1], y=[s['y'] for s in g1], mag=[s['mag'] for s in g1])
        B = dict(x=[s['x'] for s in g2], y=[s['y'] for s in g2], mag=[s['mag'] for s in g2])
        m = dp.match_lists(A, B, model=a.match_model, radius=a.match_radius)
        if not m['ok']:
            raise SystemExit('DAOMATCH failed: %s' % m.get('reason'))
        T = m['transform']
        rows = []
        for ia, ib, sep in m['pairs']:
            q = dp.apply_transform(T, [[A['x'][ia], A['y'][ia]]])[0]
            rows.append(dict(ID1=g1[ia]['id'], ID2=g2[ib]['id'], X1=A['x'][ia] + 1, Y1=A['y'][ia] + 1, X2=B['x'][ib] + 1, Y2=B['y'][ib] + 1, SEP=sep, MAG1=A['mag'][ia],
                             MAG2=B['mag'][ib], COLOR=A['mag'][ia] - B['mag'][ib], MAGERR1=g1[ia]['magerr'], MAGERR2=g2[ib]['magerr']))
        cols = ['ID1', 'ID2', 'X1', 'Y1', 'X2', 'Y2', 'SEP', 'MAG1', 'MAG2', 'COLOR', 'MAGERR1', 'MAGERR2']
        tsvio.write_table(W(a.work, 'cmd.tsv'), cols, rows)
        mas = dp.master_list([A, B], ref=1, model=a.match_model, radius=a.match_radius)
        tsvio.write_table(W(a.work, 'master.tsv'), ['id', 'x', 'y', 'mag', 'sigma', 'nframes'], mas['master'])
        plot_cmd(rows, W(a.work, 'cmd.png'))
        json.dump(dict(transform=T, n_match=m['n_match'], rms=m['rms'], dmag_median=m.get('dmag_median'), dmag_scatter=m.get('dmag_scatter'), n1=len(g1), n2=len(g2)),
                  open(W(a.work, 'match.json'), 'w'), default=lambda o: None)
        text.append('DAOMATCH: %d of %d / %d stars matched, transform %s, rms residual %.3f px, median mag offset %.3f (scatter %.3f) -> %s, %s' % (m['n_match'], len(g1), len(g2), T['kind'] + ('/%d' % T.get('degree', 1) if T['kind'] == 'poly' else ''), m['rms'], m.get('dmag_median', float('nan')), m.get('dmag_scatter', float('nan')), W(a.work, 'cmd.tsv'), W(a.work, 'master.tsv')))
        if a.catalog:
            byid = {r['ID1']: r for r in rows}
            def build(s, d):
                r = byid.get(s['id'])
                if not r:
                    return dict(DAO_MATCH=0)
                return dict(DAO_MATCH=1, DAO_MX=r['X2'], DAO_MY=r['Y2'], DAO_MSEP=r['SEP'], DAO_MMAG2=r['MAG2'], DAO_COLOR=r['COLOR'])
            emit = map_catalog(a.catalog, s1, a.match_radius, build)
    elif mode == 'arttest':
        import importlib.util
        spec = importlib.util.spec_from_file_location('ogf_completeness', os.path.join(_root, 'plugins', 'completeness', 'completeness.py'))
        cm = importlib.util.module_from_spec(spec)
        sys.modules['ogf_completeness'] = cm
        spec.loader.exec_module(cm)
        model = load_psf()
        det = dp.DaophotDetector(model.to_dict(), thresh=a.find_thresh, n_iter=2, zp=a.mag_zeropoint, fwhm=a.fwhm, snr_min=a.snr_min, chi_max=a.chi_max)
        d2, m2 = data, mask
        if a.crop_size and min(data.shape) > a.crop_size:
            cy, cx = data.shape[0] // 2, data.shape[1] // 2
            h = a.crop_size // 2
            d2 = np.ascontiguousarray(data[cy - h:cy - h + a.crop_size, cx - h:cx - h + a.crop_size])
        res = cm.run_completeness(d2, det, mask=m2, kind='star', mag_min=a.art_mag_min, mag_max=a.art_mag_max, n_bins=a.art_bins, per_bin=a.art_per_bin, per_image=a.art_per_image,
                                  zp=a.mag_zeropoint, psf=model, match_radius=max(a.match_radius, 1.5), seed=a.seed, avoid_detected=False, min_sep=2.5 * a.fwhm, n_workers=a.n_workers)
        json.dump(cm.json_clean(res), open(W(a.work, 'art.json'), 'w'), indent=1, allow_nan=False)
        text.append('ARTIFICIAL STARS (DAOPHOT fitter as the detector, injected into the crowded field): %d injected, %d recovered; 50%% limit %s, 90%% limit %s' % (res['n_injected'], res['n_recovered'], tsvio.fmt(res['lim50'], 4), tsvio.fmt(res['lim90'], 4)))
        text.append('mag\tn_inj\tn_rec\tfrac\tdmag_bias\tdmag_scatter')
        for b in res['bins']:
            text.append('%.2f\t%d\t%d\t%.3f\t%s\t%s' % (b['mag'], b['n_inj'], b['n_rec'], b['frac'], tsvio.fmt(b['bias'], 3), tsvio.fmt(b['scatter'], 3)))
        if a.catalog:
            cols_, rows_ = tsvio.read_catalog(a.catalog)
            emit = [(r['NUMBER'], dict(DAO_ART_LIM50=res['lim50'], DAO_ART_LIM90=res['lim90'])) for r in rows_]
    else:
        raise SystemExit('unknown --mode %s' % mode)
    if emit is not None:
        cols = COLUMNS if mode in ('fit', 'run') else (MATCH_COLUMNS if mode == 'match' else ['DAO_ART_LIM50', 'DAO_ART_LIM90'])
        for t in text:
            log(t)
        tsvio.write_columns(cols, emit)
    else:
        for t in text:
            print(t)
    return 0


if __name__ == '__main__':
    sys.exit(main())
