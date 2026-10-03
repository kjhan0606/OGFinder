#!/usr/bin/env python3
"""Background and noise model of the current image, and photometric errors that respect correlated noise.

    noisemodel.py IMAGE --work DIR [--catalog TSV] [--mask FITS] [--model mesh|poly|mesh+poly] [--bw 64] [--poly-order 2] [--radii 1,1.5,2,3,4,6,8,12]
        [--n-aper 500] [--aperture auto|fixed] [--aper-radius 3] [--gain 0] [--correct] [--mag-zeropoint 25] [--seed 1]

With --catalog the add_columns contract is followed (stdout = NUMBER + NM_* columns).  Files in DIR/noisemodel_*:
bkg.fits (background model), rms.fits (local pixel rms map), sub.fits (image - model, with --correct), json (summary), noise_curve.png, curve.tsv, acf.tsv.
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

from ogfkit import tsvio, imageio, noise as nz, meta as ometa, photerr as pe  # noqa: E402

COLUMNS = ['NM_SKY', 'NM_RMS', 'NM_NAPER', 'NM_CORR', 'NM_FLUXERR', 'NM_MAGERR', 'NM_SNR', 'NM_FLUX_AP', 'NM_MAG_AP', 'NM_MAGERR_AP', 'NM_SKYERR', 'NM_FLUXERR_LOC', 'NM_CONTAM', 'NM_FLUXERR_TOT']


def build_model(data, bad, model='mesh', bw=64, fw=3, poly_order=2, grow=3.0, thresh=2.0, minarea=8):
    """Returns dict(mask, bkg, rms, poly (surface or None), sub, summary pieces)."""
    mask = nz.source_mask(data, bad, thresh=thresh, minarea=minarea, grow=grow, bw=bw)
    d = np.where(np.isfinite(data), data, 0.0)
    poly = None
    coeffs = None
    if model in ('poly', 'mesh+poly'):
        poly, coeffs, prms = nz.fit_poly_background(np.where(mask, np.nan, d), mask, poly_order)
        base = d - poly
    else:
        base = d
    if model == 'poly':
        bkg = poly
        _, rms = nz.mesh_background(d - poly, mask, bw, fw)
    else:
        b, rms = nz.mesh_background(base, mask, bw, fw)
        bkg = b + (poly if poly is not None else 0.0)
    return dict(mask=mask, bkg=bkg, rms=rms, poly=poly, coeffs=None if coeffs is None else [float(c) for c in coeffs], sub=d - bkg)


def analyse(data, bad, a, radii):
    """Server-ready entry point: returns (summary dict, arrays dict)."""
    m = build_model(data, bad, a.model, a.bw, a.fw, a.poly_order, getattr(a, 'mask_grow', 3.0), getattr(a, 'mask_thresh', 2.0), getattr(a, 'mask_minarea', 8))
    sub, mask = m['sub'], m['mask']
    flat0 = nz.flatness(data - np.median(data[~mask]), mask, 6)          # before: the raw image with only a constant removed
    flat1 = nz.flatness(sub, mask, 6)
    acf = nz.autocorrelation(sub, mask, 10)
    cs = nz.correlation_summary(acf)
    an = nz.aperture_noise(sub, mask, radii, a.n_aper, a.seed)
    law = nz.fit_noise_law(an['N'], an['sigma'], an['sigma_pix'])
    sp = an['sigma_pix']
    corr = [s / (sp * math.sqrt(n)) for s, n in zip(an['sigma'], an['N'])]
    summary = dict(model=a.model, bw=a.bw, poly_order=a.poly_order if a.model != 'mesh' else None, sky_median=float(np.median(m['bkg'])), sky_rms_variation=float(np.std(m['bkg'])),
                   masked_fraction=float(mask.mean()), sigma_pix=sp, rms_map_median=float(np.median(m['rms'])), flatness_before=flat0, flatness_after=flat1,
                   acf=cs, noise_law=law, aperture_noise=an, correlation_factor=corr, poly_coeffs=m['coeffs'])
    return summary, dict(m, acf=acf, law=law, an=an)


def plot(summary, arrays, path):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        return
    an = arrays['an']
    law = summary['noise_law']
    N = np.array(an['N']); s = np.array(an['sigma']); sp = an['sigma_pix']
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.2))
    ax[0].loglog(N, s, 'ko', label='blank apertures')
    ax[0].loglog(N, sp * np.sqrt(N), 'b--', label='uncorrelated sqrt(N)')
    if law.get('ok'):
        ax[0].loglog(N, nz.noise_law_sigma(N, law, sp), 'r-', label='alpha N^beta (beta %.2f)' % law['beta'])
    ax[0].set_xlabel('N pixels'); ax[0].set_ylabel('sigma of the sum'); ax[0].legend(fontsize=6)
    acf = arrays['acf']; mm = acf.shape[0] // 2
    ax[1].imshow(acf, origin='lower', cmap='viridis', vmin=-0.1, vmax=1, extent=[-mm, mm, -mm, mm]); ax[1].set_title('pixel autocorrelation', fontsize=8)
    fl = np.array(summary['flatness_after']['medians'])
    im = ax[2].imshow(fl - np.nanmedian(fl), origin='lower', cmap='RdBu_r'); ax[2].set_title('sky block medians (after)', fontsize=8)
    fig.colorbar(im, ax=ax[2], shrink=0.8)
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)


def catalog_errors(rows, data, arrays, summary, a):
    sep = nz._sep()
    sub, rms, law, sp = arrays['sub'], arrays['rms'], arrays['law'], summary['sigma_pix']
    ny, nx = data.shape
    out = []
    xs = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows])
    ys = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
    ok = np.isfinite(xs) & np.isfinite(ys)
    loc_rms = np.where(ok, nz.sample_map(rms, np.where(ok, xs, 0), np.where(ok, ys, 0)), np.nan)
    loc_sky = np.where(ok, nz.sample_map(arrays['bkg'], np.where(ok, xs, 0), np.where(ok, ys, 0)), np.nan)
    if a.aperture == 'fixed':
        Nap = np.full(len(rows), math.pi * a.aper_radius ** 2)
        fl, _, _ = sep.sum_circle(np.ascontiguousarray(sub, np.float64), np.where(ok, xs, 0), np.where(ok, ys, 0), a.aper_radius, subpix=5)
        flux = np.where(ok, fl, np.nan)
    else:
        flux = np.array([tsvio.fnum(r.get('FLUX_AUTO')) for r in rows])
        A = np.array([tsvio.fnum(r.get('A_IMAGE'), 2.0) for r in rows]); B = np.array([tsvio.fnum(r.get('B_IMAGE'), 2.0) for r in rows]); K = np.array([tsvio.fnum(r.get('KRON_RADIUS'), 3.5) for r in rows])
        Nap = math.pi * (a.kron_fact * K * A) * (a.kron_fact * K * B)
        npx = np.array([tsvio.fnum(r.get('NPIX_ISO')) for r in rows])
        Nap = np.where(np.isfinite(Nap), Nap, npx)
    model = pe.NoiseModel(sp, law['alpha'], law['beta'], bool(law.get('ok')))
    ann = _annulus(a)
    extra = {}
    if ann:
        radii = [float(v) for v in a.radii.replace(';', ',').split(',') if v.strip() and float(v) <= ann[0] - 2]
        model.add_local(sub, arrays['mask'], radii + [r for r in (2.0, 3.0, 5.0, 8.0) if r <= ann[0] - 2 and r not in radii], ann, 'median', n=a.n_aper, seed=a.seed)
        extra = unified_errors(rows, model, ann, xs, ys, ok, flux, Nap, loc_rms, sp, a)
        arrays['local'] = model.local
    ratio = np.where(np.isfinite(loc_rms), loc_rms / sp, 1.0)
    for i, r in enumerate(rows):
        d = {}
        if ok[i]:
            n_ = Nap[i]
            d['NM_SKY'] = loc_sky[i]; d['NM_RMS'] = loc_rms[i]; d['NM_NAPER'] = n_
            if np.isfinite(n_) and n_ > 0:
                corr = nz.noise_law_sigma(n_, law, 1.0) / math.sqrt(n_) if law.get('ok') else 1.0
                d['NM_CORR'] = corr
                f = flux[i]
                e = float(pe.aperture_error(f if np.isfinite(f) else 0.0, n_, loc_rms[i], model, a.gain or None)['ap'])
                d['NM_FLUXERR'] = e
                if np.isfinite(f) and f > 0:
                    d['NM_MAGERR'] = 1.0857 * e / f
                    d['NM_SNR'] = f / e
                if a.aperture == 'fixed' and np.isfinite(f):
                    d['NM_FLUX_AP'] = f
                    if f > 0:
                        d['NM_MAG_AP'] = a.mag_zeropoint - 2.5 * math.log10(f)
                        d['NM_MAGERR_AP'] = 1.0857 * e / f
        for k, v in extra.items():
            if ok[i] and np.isfinite(v[i]):
                d[k] = float(v[i])
        out.append((r['NUMBER'], d))
    return out


def _annulus(a):
    t = [v for v in (getattr(a, 'sky_annulus', '') or '').replace(';', ',').split(',') if v.strip()]
    return (float(t[0]), float(t[1])) if len(t) == 2 and 0 < float(t[0]) < float(t[1]) else None


def _psf_callable(path):
    """PSFModel json (psfex) -> callable(x, y, dx, dy) -> unit-sum stamp; FITS -> constant stamp."""
    if not path or not os.path.isfile(path):
        return None
    if path.lower().endswith('.json'):
        from ogfkit import psfmodel as pm
        mdl = pm.load_model(path)
        return lambda x, y, dx, dy: mdl.stamp(x, y, dx, dy)
    from astropy.io import fits
    st = np.nan_to_num(np.asarray(fits.getdata(path), float))
    st = np.clip(st, 0, None)
    return st / st.sum()


def unified_errors(rows, model, ann, xs, ys, ok, flux, Nap, loc_rms, sp, a):
    """Columns of the unified error model (ogfkit.photerr): NM_SKYERR (sky-estimation error N * sigma_sky, counts, from the annulus factor), NM_FLUXERR_LOC (error of the
    local-sky-subtracted aperture flux: measured blank-aperture law with the same annulus, + Poisson), NM_CONTAM (flux the PSF wings of the other catalogue sources
    put into the aperture net of the annulus sky, counts; needs --psf), NM_FLUXERR_TOT (NM_FLUXERR_LOC and the contamination uncertainty in quadrature).  Apertures with N outside
    [N_min/2, 2 N_max] of the measured local law are left empty."""
    n = len(rows)
    nan = lambda: np.full(n, np.nan)
    out = dict(NM_SKYERR=nan(), NM_FLUXERR_LOC=nan(), NM_CONTAM=nan(), NM_FLUXERR_TOT=nan())
    if not (model.local and model.local.get('ok')):
        return out
    lo, hi = min(model.local['N']) / 2, max(model.local['N']) * 2
    good = ok & np.isfinite(Nap) & (Nap >= lo) & (Nap <= hi)
    N = np.where(good, Nap, 1.0)
    f = np.where(np.isfinite(flux), flux, 0.0)
    e = pe.aperture_error(f, N, np.where(np.isfinite(loc_rms), loc_rms, sp), model, a.gain or None, local=True)
    out['NM_FLUXERR_LOC'] = np.where(good, e['ap'], np.nan)
    nsky = math.pi * (ann[1] ** 2 - ann[0] ** 2)
    out['NM_SKYERR'] = np.where(good, N * pe.sky_est_sigma(nsky, np.where(np.isfinite(loc_rms), loc_rms, sp), model), np.nan)
    out['NM_FLUXERR_TOT'] = out['NM_FLUXERR_LOC'].copy()
    psf = _psf_callable(getattr(a, 'psf', ''))
    if psf is not None and np.isfinite(flux).any():
        xy = np.c_[np.where(ok, xs, 0), np.where(ok, ys, 0)]
        rr = np.round(np.sqrt(np.where(good, N, 1.0) / math.pi) * 2) / 2
        fe = e['ap'] / 1.0
        for r in np.unique(rr[good]):
            tg = np.where(good & (rr == r))[0]
            r_use = min(float(r), ann[0] - 1.0)
            c = pe.neighbour_contamination(xy, np.where(np.isfinite(flux) & ok, flux, 0.0), r_use, psf, sky_annulus=ann, flux_err=np.where(ok, fe, 0.0), psf_rel_err=0.1, targets=tg)
            out['NM_CONTAM'][tg] = c['contam'][tg]
            out['NM_FLUXERR_TOT'][tg] = np.sqrt(out['NM_FLUXERR_LOC'][tg] ** 2 + c['err'][tg] ** 2)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--work', default='.')
    ap.add_argument('--catalog', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--model', default='mesh', choices=['mesh', 'poly', 'mesh+poly'])
    ap.add_argument('--bw', type=int, default=64)
    ap.add_argument('--fw', type=int, default=3)
    ap.add_argument('--poly-order', type=int, default=2)
    ap.add_argument('--mask-thresh', type=float, default=2.0)
    ap.add_argument('--mask-grow', type=float, default=3.0)
    ap.add_argument('--mask-minarea', type=int, default=8)
    ap.add_argument('--radii', default='1,1.5,2,3,4,6,8,12')
    ap.add_argument('--n-aper', type=int, default=500)
    ap.add_argument('--aperture', default='auto', choices=['auto', 'fixed'])
    ap.add_argument('--aper-radius', type=float, default=3.0)
    ap.add_argument('--kron-fact', type=float, default=2.5)
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--gain', type=float, default=0.0)
    ap.add_argument('--correct', action='store_true')
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--sky-annulus', default='', help='R_IN,R_OUT of the local sky annulus: adds NM_SKYERR / NM_FLUXERR_LOC (unified error model)')
    ap.add_argument('--psf', default='', help='PSF model json (psfex) or FITS stamp: adds NM_CONTAM / NM_FLUXERR_TOT (needs --sky-annulus)')
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, 'noisemodel_' + n)
    data, hdr = imageio.load_image(a.image)
    bad = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    bad = (np.zeros(data.shape, bool) if bad is None else bad) | ~np.isfinite(data) | (data == 0)
    radii = [float(v) for v in a.radii.replace(';', ',').split(',') if v.strip()]
    summary, arr = analyse(np.asarray(data, np.float64), bad, a, radii)
    imageio.save_fits(W('bkg.fits'), arr['bkg'].astype(np.float32))
    imageio.save_fits(W('rms.fits'), arr['rms'].astype(np.float32))
    if a.correct:
        imageio.save_fits(W('sub.fits'), arr['sub'].astype(np.float32))
    with open(W('summary.json'), 'w') as fh:
        json.dump(summary, fh, default=lambda o: None)
    tsvio.write_table(W('curve.tsv'), ['R', 'N', 'SIGMA', 'SIGMA_NAIVE', 'SIGMA_LAW', 'CORR'],
                      [dict(R=r, N=n, SIGMA=s, SIGMA_NAIVE=summary['sigma_pix'] * math.sqrt(n), SIGMA_LAW=float(nz.noise_law_sigma(n, arr['law'], summary['sigma_pix'])), CORR=c)
                       for r, n, s, c in zip(arr['an']['radii'], arr['an']['N'], arr['an']['sigma'], summary['correlation_factor'])])
    plot(summary, arr, W('noise_curve.png'))
    law = arr['law']
    f1 = summary['flatness_after']
    lines = ['noise model of %s (%s background, bw %d): sky %.5g (variation %.3g), pixel sigma %.5g, %.1f %% of the pixels masked' % (
        os.path.basename(a.image), a.model, a.bw, summary['sky_median'], summary['sky_rms_variation'], summary['sigma_pix'], 100 * summary['masked_fraction']),
        'pixel correlation: rho(1,0) %.3f rho(0,1) %.3f rho(1,1) %.3f, kernel FWHM ~ %.2f px' % (summary['acf']['rho1_x'], summary['acf']['rho1_y'], summary['acf']['rho_diag'], summary['acf']['kernel_fwhm_px']),
        'aperture noise law: sigma_N = sigma_pix * %.3f * N^%.3f (uncorrelated: 1 * N^0.5); correction vs sqrt(N): %s' % (law['alpha'], law['beta'], ', '.join('r=%g: x%.2f' % (r, c) for r, c in zip(arr['an']['radii'], summary['correlation_factor']))),
        'flatness (6x6 sky blocks): peak-to-peak %.3g before -> %.3g after, chi2 %.1f -> %.1f, plane %.3g' % (summary['flatness_before']['peak_to_peak'], f1['peak_to_peak'], summary['flatness_before']['chi2_flat'], f1['chi2_flat'], f1['plane_amp'])]
    if a.catalog:
        cols, rows = tsvio.read_catalog(a.catalog)
        out = catalog_errors(rows, data, arr, summary, a)
        if a.meta_out:
            ometa.update(a.meta_out, 'noisemodel', dict(sigma_pix=summary['sigma_pix'], alpha=law['alpha'], beta=law['beta'], rho1=summary['acf']['rho1'], model=a.model,
                         kernel_fwhm_px=summary['acf']['kernel_fwhm_px'], local=arr.get('local')), nrows=len(rows))
        for t in lines:
            sys.stderr.write(t + '\n')
        tsvio.write_columns(COLUMNS, out)
    else:
        print('\n'.join(lines))
    return 0


if __name__ == '__main__':
    sys.exit(main())
