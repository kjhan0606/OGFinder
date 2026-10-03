"""Validation of ogfkit/photerr.py (unified photometric error model).

blank   : aperture photometry (local annulus sky) on blank apertures of a correlated-noise image (synthetic, kernel sigma 1.0 px) and of the real HUDF F160W
          drizzled image; pull = flux / error for the naive sqrt(N) error, the legacy DAOPHOT term and the unified model.  Law and annulus factor are measured
          on one set of blank positions, evaluated on an independent set.
inject  : stars of known flux (Moffat, FWHM 3 px, pixel-integrated, with source Poisson noise) injected into the same images, sparse and crowded;
          pull = (estimated total flux - true) / error with / without neighbour-contamination correction and aperture-correction error.
Usage: photerr_validate.py OUT.json [--quick]
"""
import json
import math
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..'))
from ogfkit import noise as nz, photerr as pe, models as M  # noqa: E402

HUDF = '/workspace/fits/hudf_f160w.fits'
FWHM, BETA, PSIZE = 3.0, 2.5, 61


def synth_image(n=1200, sigma=5.0, sk=1.0, seed=1):
    rng = np.random.default_rng(seed)
    w = rng.normal(size=(n, n))
    c = gaussian_filter(w, sk, mode='wrap')
    return 100.0 + c * (sigma / c.std())


def real_image(box=(1000, 1000, 2400, 2400)):
    from astropy.io import fits
    d = fits.getdata(HUDF).astype(float)
    x0, y0, x1, y1 = box
    return np.nan_to_num(d[y0:y1, x0:x1])


def pulls_stats(p):
    p = np.asarray(p, float)
    p = p[np.isfinite(p)]
    if len(p) < 5:
        return dict(n=int(len(p)))
    nm = 1.4826 * np.median(np.abs(p - np.median(p)))
    return dict(n=int(len(p)), mean=float(p.mean()), std=float(p.std()), nmad=float(nm), within1=float((np.abs(p) < 1).mean()), within2=float((np.abs(p) < 2).mean()))


def blank_test(img, label, radii=(2.0, 3.0, 5.0, 8.0), gap=4.0, width=6.0, n=800):
    mask = nz.source_mask(img)
    bkg = np.median(img[~mask])
    sub = img - bkg
    sp = nz.robust_std(sub[~mask][:2000000])
    model = pe.NoiseModel.measure(sub, mask, seed=1)
    out = dict(label=label, sigma_pix=sp, alpha=model.alpha, beta=model.beta, rows=[])
    ann0 = (12.0, 18.0)
    model.add_local(img, mask, radii + (1.5, 4.0, 6.0, 10.0), ann0, 'median', n=n, seed=2)
    for r in radii:
        ann = ann0
        kappa, nsky, nu = pe.annulus_factor(sub, mask, ann[0], ann[1], 'median', n=n, seed=3, sigma_pix=sp)
        rng = np.random.default_rng(77)                       # independent positions
        pos = nz.blank_positions(mask, ann[1] + 1, n, rng)
        ph = pe.aperture_photometry(img, pos, r, ann, 'median', mask=mask)
        N, ns, f = ph['N'], ph['nsky'], ph['flux']
        e_naive = sp * np.sqrt(N)
        e_leg = np.sqrt(N * sp ** 2 * (1 + N / np.maximum(ns, 1)))
        e_law = pe.aperture_error(f, N, sp, model)['ap']                                     # correlated sky term only
        e_uni = pe.aperture_error(f, N, sp, model, nsky=ns, sky_factor=kappa)['ap']          # + measured sky-estimation error
        e_loc = pe.aperture_error(f, N, sp, model, nsky=ns, local=True)['ap']            # measured local-sky-subtracted law
        row = dict(r=r, N=float(np.nanmedian(N)), nsky=float(np.nanmedian(ns)), kappa=kappa, std_flux=float(np.nanstd(f)),
                   err_naive=float(np.nanmedian(e_naive)), err_legacy=float(np.nanmedian(e_leg)), err_law=float(np.nanmedian(e_law)), err_unified=float(np.nanmedian(e_uni)),
                   pull_naive=pulls_stats(f / e_naive), pull_legacy=pulls_stats(f / e_leg), pull_law=pulls_stats(f / e_law), pull_unified=pulls_stats(f / e_uni),
                   pull_local_law=pulls_stats(f / e_loc))
        out['rows'].append(row)
        print('%s r=%g  std(flux)=%.3g  err: naive %.3g legacy %.3g law %.3g unified %.3g  | pull std: naive %.2f legacy %.2f law %.2f unified %.2f (local-law %.2f)' % (
            label, r, row['std_flux'], row['err_naive'], row['err_legacy'], row['err_law'], row['err_unified'], row['pull_naive']['std'], row['pull_legacy']['std'], row['pull_law']['std'],
            row['pull_unified']['std'], row['pull_local_law']['std']), flush=True)
    return out


def inject_test(img0, label, n_stars, seed, r=4.0, ann=(8.0, 14.0), gain=4.0, snr_lim=3.0, hidden_snr=0.0):
    rng = np.random.default_rng(seed)
    mask0 = nz.source_mask(img0)
    sky = np.median(img0[~mask0])
    sub0 = img0 - sky
    sp = nz.robust_std(sub0[~mask0][:2000000])
    model = pe.NoiseModel.measure(sub0, mask0, seed=1)
    kappa, nsky_t, _ = pe.annulus_factor(sub0, mask0, ann[0], ann[1], 'median', n=600, seed=5, sigma_pix=sp)
    ny, nx = img0.shape
    # injection positions: anywhere in the (unmasked) image so the neighbours are real; stars on masked regions are excluded from the evaluation
    x = rng.uniform(40, nx - 41, n_stars)
    y = rng.uniform(40, ny - 41, n_stars)
    scale = sp / 5.0                      # all fluxes in units where sigma_pix = 5 (synthetic image)
    gain = gain / scale
    flux = scale * 10 ** rng.uniform(2.0, 4.0, n_stars)
    # + 40 bright isolated PSF stars (>= 45 px from everything)
    xs, ys, fs = list(x), list(y), list(flux)
    tries = 0
    while len(xs) < n_stars + 40 and tries < 20000:
        tries += 1
        xn, yn = rng.uniform(40, nx - 41), rng.uniform(40, ny - 41)
        if min(np.hypot(np.array(xs) - xn, np.array(ys) - yn)) > 45:
            xs.append(xn); ys.append(yn); fs.append(scale * 10 ** rng.uniform(4.0, 4.7))
    x, y, flux = np.array(xs), np.array(ys), np.array(fs)
    n_stars = len(x)
    img = img0.copy()
    h = PSIZE // 2
    stars = np.zeros((ny, nx))
    for xi, yi, fi in zip(x, y, flux):
        ix, iy = int(round(xi)), int(round(yi))
        st = M.moffat_psf(PSIZE, FWHM, BETA, dx=xi - ix, dy=yi - iy)
        stars[iy - h:iy + h + 1, ix - h:ix + h + 1] += fi * st
    img = img + stars + rng.normal(size=img.shape) * np.sqrt(np.clip(stars, 0, None) / gain)
    xy = np.c_[x, y]
    # PSF used = unit-sum Moffat of the same shape from a stack of isolated bright stars is not needed: the "measured" PSF is the stamp built from the
    # brightest isolated injected stars (as a user would), then fractions come from it
    N = math.pi * r * r
    model.add_local(img0, mask0, (2.0, 3.0, 4.0, 6.0, 8.0), ann, 'median', n=500, seed=2)
    ph = pe.aperture_photometry(img, xy, r, ann, 'median')
    f_ap = ph['flux']; ns = ph['nsky']
    # empirical aperture fraction from isolated bright stars (reference: flux inside a large aperture)
    from scipy.spatial import cKDTree
    tree = cKDTree(xy)
    d2 = tree.query(xy, k=2)[0][:, 1]
    bright = np.where((flux > 8000 * scale) & (d2 > 40))[0]
    from ogfkit import psfext as px
    psf, nst = px.stack_star_psf(img - sky, xy[bright], size=PSIZE, bkg=None, nmax=40, rmax=28)
    frac = pe.aperture_fraction(psf, r, ann)
    bs = []
    for k in range(20):                                    # bootstrap over the PSF stars -> error of the aperture correction
        sel = rng.choice(bright, len(bright), replace=True)
        pk, _ = px.stack_star_psf(img - sky, xy[sel], size=PSIZE, bkg=None, nmax=40, rmax=28)
        bs.append(pe.aperture_fraction(pk, r, ann))
    frac_sem = float(np.std(bs, ddof=1)); frac_scat = frac_sem
    true_frac = pe.aperture_fraction(M.moffat_psf(PSIZE, FWHM, BETA), r, ann)
    f_est0 = f_ap / frac
    # neighbours: detected sources (aperture S/N > snr_lim)
    err0 = pe.aperture_error(f_ap, N, sp, model, gain=gain, local=True)
    det = f_ap / err0['ap'] > snr_lim
    ids = np.where(det)[0]
    fcur = f_est0.copy()
    ecur = err0['total']
    for it in range(3):
        nb = pe.neighbour_contamination(xy[ids], fcur[ids], r, psf, sky_annulus=ann, flux_err=ecur[ids], psf_rel_err=0.1)
        con = np.zeros(len(xy)); cerr = np.zeros(len(xy))
        con[ids] = nb['contam']; cerr[ids] = nb['err']
        # targets that are not detected but evaluated: compute for all targets against the detected list
        tgt = np.setdiff1d(np.arange(len(xy)), ids)
        if len(tgt):
            full = pe.neighbour_contamination(xy, np.where(det, fcur, 0.0), r, psf, sky_annulus=ann, flux_err=np.where(det, ecur, 0.0), psf_rel_err=0.1, targets=tgt)
            con[tgt] = full['contam'][tgt]; cerr[tgt] = full['err'][tgt]
        e = pe.aperture_error(f_ap - con, N, sp, model, gain=gain, contam_err=cerr, frac=frac, frac_err=frac_sem, local=True)
        fcur = e['flux_total']; ecur = e['total']
    # evaluation sample: stars whose aperture and annulus lie on pixels free of the real (masked) sources of the original image
    ok_pos = np.array([not mask0[max(int(yy) - int(ann[1]), 0):int(yy) + int(ann[1]) + 2, max(int(xx) - int(ann[1]), 0):int(xx) + int(ann[1]) + 2].any() for xx, yy in xy])
    nn_d = tree.query(xy, k=2)[0][:, 1]
    # neighbour sample: stars with a close neighbour (within the annulus)
    close = nn_d < ann[1] + r
    e_naive = sp * np.sqrt(N)
    e_leg = np.sqrt(N * sp ** 2 * (1 + N / np.maximum(ns, 1)) + np.clip(f_ap, 0, None) / gain)
    e_loc = pe.aperture_error(f_ap, N, sp, model, gain=gain, frac=frac, frac_err=frac_sem, local=True)
    e_nc = pe.aperture_error(f_ap, N, sp, model, gain=gain, nsky=ns, sky_factor=kappa, frac=frac, frac_err=frac_sem)    # global law + annulus factor, without contamination
    res = dict(label=label, n_stars=int(n_stars), r=r, ann=list(ann), gain=gain, sigma_pix=sp, alpha=model.alpha, beta=model.beta, kappa=kappa,
               frac_measured=frac, frac_true=true_frac, frac_sem=frac_sem, frac_scatter=frac_scat, n_bright_psf=int(len(bright)), n_detected=int(det.sum()), n_close=int(close.sum()))
    sel_all = ok_pos
    sel_close = ok_pos & close
    sel_iso = ok_pos & ~close
    truth = flux
    for nm, sel in (('all', sel_all), ('isolated', sel_iso), ('close_neighbour', sel_close)):
        for kind, est, er in (('naive', f_est0, e_naive / frac), ('legacy', f_est0, e_leg / frac), ('global_law_no_contam', f_est0, e_nc['total']), ('unified_no_contam', f_est0, e_loc['total']), ('unified_full', fcur, ecur)):
            res['%s/%s' % (nm, kind)] = pulls_stats(((est - truth) / er)[sel])
        res['%s/bias_raw_over_flux' % nm] = float(np.median(((f_est0 - truth) / truth)[sel])) if sel.any() else None
        res['%s/bias_corrected_over_flux' % nm] = float(np.median(((fcur - truth) / truth)[sel])) if sel.any() else None
    # by S/N
    for lo, hi in ((100, 300), (300, 1000), (1000, 1e5)):
        s = sel_all & (truth >= lo) & (truth < hi)
        res['flux_%g_%g/unified_full' % (lo, hi)] = pulls_stats(((fcur - truth) / ecur)[s])
        res['flux_%g_%g/naive' % (lo, hi)] = pulls_stats(((f_est0 - truth) / (e_naive / frac))[s])
    print(label, 'frac meas %.4f true %.4f sem %.4f | n_eval %d (close %d)' % (frac, true_frac, frac_sem, sel_all.sum(), sel_close.sum()))
    for nm in ('all', 'isolated', 'close_neighbour'):
        print('  %-16s' % nm, ' '.join('%s std %.2f mean %+.2f (n %d)' % (k, res['%s/%s' % (nm, k)].get('std', float('nan')), res['%s/%s' % (nm, k)].get('mean', float('nan')), res['%s/%s' % (nm, k)]['n']) for k in ('naive', 'legacy', 'global_law_no_contam', 'unified_no_contam', 'unified_full')))
    return res


if __name__ == '__main__':
    out = sys.argv[1] if len(sys.argv) > 1 else 'photerr_report.json'
    quick = '--quick' in sys.argv
    syn = synth_image(900 if quick else 1200)
    report = dict(blank=[], inject=[])
    report['blank'].append(blank_test(syn, 'synthetic (kernel sigma 1.0 px)'))
    if os.path.exists(HUDF):
        hu = real_image()
        report['blank'].append(blank_test(hu, 'HUDF F160W 1400x1400'))
    report['inject'].append(inject_test(syn, 'synthetic sparse', 150, 11))
    report['inject'].append(inject_test(syn, 'synthetic moderate', 500 if not quick else 300, 15))
    report['inject'].append(inject_test(syn, 'synthetic crowded', 1500 if not quick else 900, 12))
    if os.path.exists(HUDF):
        hu = real_image((1000, 1000, 2000, 2000))
        report['inject'].append(inject_test(hu, 'HUDF F160W sparse', 100, 13))
        report['inject'].append(inject_test(hu, 'HUDF F160W crowded', 600, 14))
    json.dump(report, open(out, 'w'), indent=1)
