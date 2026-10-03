"""Real-data check: SDSS J0946+1006 (SLACS, HST ACS/WFC F814W, program 10886), lens-light subtraction + pixelated source inversion with an SIE + shear
lens whose parameters are fitted by maximising the Bayesian evidence.  Data are public (MAST) and are cached under $OGF_DATA_CACHE or
~/.cache/ogfinder_regression/lens (downloaded on demand, ~215 MB, cut to 500 x 500 pixels)."""
import os, sys, json, math, time
import numpy as np
ROOT = os.environ.get('OGF_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
sys.path.insert(0, ROOT)
from ogfkit import lensmodel as L, lensextra as X

CACHE = os.path.expanduser(os.environ.get('OGF_DATA_CACHE', '~/.cache/ogfinder_regression')) 
CACHE = CACHE if CACHE.endswith('lens') else os.path.join(CACHE, 'lens')
URL = 'https://mast.stsci.edu/api/v0.1/Download/file?uri=mast:HST/product/j9op14010_drc.fits'
RA, DEC = 146.73617, 10.11467


def get_cutout(size=500):
    cut = os.path.join(CACHE, 'j0946_f814w_cut.fits')
    if os.path.isfile(cut):
        return cut
    import requests
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    from astropy.coordinates import SkyCoord
    from scipy.ndimage import gaussian_filter
    os.makedirs(CACHE, exist_ok=True)
    full = os.path.join(CACHE, 'j9op14010_drc_full.fits')
    if not os.path.isfile(full):
        with requests.get(URL, stream=True, timeout=600) as r:
            r.raise_for_status()
            with open(full + '.part', 'wb') as f:
                for ch in r.iter_content(1 << 20):
                    f.write(ch)
        os.replace(full + '.part', full)
    h = fits.open(full)
    sci = h['SCI']
    w = WCS(sci.header)
    c0 = Cutout2D(sci.data, SkyCoord(RA, DEC, unit='deg'), (400, 400), wcs=w)
    s = gaussian_filter(np.nan_to_num(c0.data), 5)
    y, x = np.unravel_index(np.argmax(s), s.shape)
    cx, cy = c0.to_original_position((x, y))
    c2 = Cutout2D(sci.data, (cx, cy), (size, size), wcs=w)
    hd = c2.wcs.to_header()
    hd['EXPTIME'] = float(h[0].header.get('EXPTIME', 2096.0))
    for k in ('PHOTFLAM', 'PHOTZPT', 'BUNIT'):
        if k in sci.header:
            hd[k] = sci.header[k]
    hd['FILTER'] = 'F814W'
    fits.writeto(cut, c2.data.astype('float32'), hd, overwrite=True)
    return cut


def sersic2d(x, y, x0, y0, amp, reff, n, q, phi):
    k = 1.9992 * n - 0.3271
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    dx, dy = x - x0, y - y0
    xp, yp = c * dx + s * dy, -s * dx + c * dy
    r = np.sqrt(xp ** 2 + (yp / q) ** 2 + 1e-4)
    return amp * np.exp(-k * ((r / reff) ** (1.0 / n) - 1.0))


def lens_light(img, mask_fit, pixscale, psf_sigma):
    """two concentric elliptical Sersics (free centre, q, PA each) + constant, convolved with a Gaussian PSF, fitted to the pixels in mask_fit"""
    from scipy.optimize import least_squares
    from scipy.ndimage import gaussian_filter
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    sm = gaussian_filter(np.nan_to_num(img), 4)
    y0, x0 = np.unravel_index(np.argmax(sm), sm.shape)

    def model(v):
        a = sersic2d(xx, yy, v[0], v[1], v[2], v[3], v[4], v[5], v[6]) + sersic2d(xx, yy, v[0] + v[12], v[1] + v[13], v[7], v[8], v[9], v[10], v[11]) + v[14]
        return gaussian_filter(a, psf_sigma)

    def res(v):
        return ((model(v) - img) / sig)[mask_fit]
    sig = max(float(np.nanstd(img[:30, :30])), 1e-4)
    v0 = [x0, y0, 2.0, 8.0, 3.5, 0.8, 90.0, 1.0, 40.0, 1.0, 0.8, 90.0, 0.0, 0.0, 0.0]
    lb = [x0 - 5, y0 - 5, 0, 1, 0.5, 0.3, -360, 0, 3, 0.3, 0.3, -360, -5, -5, -1]
    ub = [x0 + 5, y0 + 5, 1e3, 200, 8, 1.0, 360, 1e3, 400, 6, 1.0, 360, 5, 5, 1]
    r = least_squares(res, v0, bounds=(lb, ub), x_scale='jac', max_nfev=60)
    return model(r.x), r.x, sig


def main(outjson):
    from astropy.io import fits
    from astropy.wcs import WCS
    cut = get_cutout()
    h = fits.open(cut)
    img = h[0].data.astype(float)
    hdr = h[0].header
    w = WCS(hdr)
    cd = w.pixel_scale_matrix
    pixscale = float(math.sqrt(abs(np.linalg.det(cd)))) * 3600.0
    psf_sigma = 0.045 / pixscale * 1.0                   # ACS/WFC F814W FWHM ~0.09 arcsec -> sigma 0.04 arcsec (Gaussian stand-in for the PSF)
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    sm_peak = np.unravel_index(np.argmax(__import__('scipy.ndimage', fromlist=['x']).gaussian_filter(img, 4)), img.shape)
    org = (float(sm_peak[1]), float(sm_peak[0]))
    r_as = np.hypot(xx - org[0], yy - org[1]) * pixscale
    # pixels used for the lens-light fit: outside the ring (0.9 .. 2.0 arcsec) and the outer fainter ring zone (2.0 .. 3.2), inside 8 arcsec
    fitmask = (r_as < 0.45) | ((r_as > 3.2) & (r_as < 8.0))
    t0 = time.time()
    ll, v, sig = lens_light(img, fitmask, pixscale, psf_sigma)
    res = img - ll
    return dict(pixscale=pixscale, org=org, lens_light_params=[float(a) for a in v], noise=sig, seconds=time.time() - t0), img, ll, res, org, pixscale, psf_sigma, w


if __name__ == '__main__' and len(sys.argv) <= 2:
    info, img, ll, res, org, ps, psf_sigma, w = main(sys.argv[1] if len(sys.argv) > 1 else '/tmp/x.json')
    print(info)
    np.save('/workspace/work/stage6/j0946_res.npy', res); np.save('/workspace/work/stage6/j0946_img.npy', img)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(15, 5))
    c = 150
    sl = (slice(int(org[1]) - c, int(org[1]) + c), slice(int(org[0]) - c, int(org[0]) + c))
    for a, d, t in zip(ax, (img, ll, res), ('data', 'lens light', 'residual')):
        a.imshow(d[sl], origin='lower', vmin=-0.02, vmax=0.15 if t != 'data' else 0.3, cmap='gray'); a.set_title(t)
    plt.savefig('/workspace/work/stage6/j0946_stage1.png', dpi=60)


def analyse(outjson):
    from scipy.ndimage import gaussian_filter
    info, img, ll, res, org0, pixscale, psf_sigma, w = main(outjson)
    v = info['lens_light_params']
    org = (v[0], v[1])                                    # lens-light centre (pixel frame, 0-based)
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r_as = np.hypot(xx - org[0], yy - org[1]) * pixscale
    sig = info['noise']
    k = X.gaussian_kernel(psf_sigma)
    inner = (r_as > 0.55) & (r_as < 2.0)
    m1 = X.arc_mask(res, sig, org, 2.0 / pixscale, nsig=4.0, dilate=3) & (r_as > 0.55)
    start = dict(theta_E=1.4, q=0.9, phi=0.0, gamma=0.03, phi_g=0.0, x0=0.0, y0=0.0)
    t0 = time.time()
    best = None
    for q0, ph0 in ((0.9, 0.0), (0.9, 90.0)):                           # two starting position angles (the evidence surface has local maxima)
        st = dict(start, q=q0, phi=ph0)
        f = X.fit_lens_pixels(res, sig, m1, pixscale, org, (0.0, 0.0), st, kernel=k, n=28, maxiter=220)
        if best is None or f['evidence'] > best['evidence']:
            best = f
    p = best['params']
    m = L.build_sie_shear(dict(p, x0=0.0, y0=0.0))
    S = X.SourceInversion(m, res, sig, m1, pixscale, org, (0.0, 0.0), kernel=k, n=40)
    r = S.best()
    neff = S.n_eff(r)
    src, mod = S.images(r)
    chi2red = r['chi2'] / (S.ndata - neff)
    # sky position angle of the major axis (east of north) from the pixel-frame angle
    cd = w.pixel_scale_matrix
    d = cd @ np.array([math.cos(math.radians(p['phi'])), math.sin(math.radians(p['phi']))])
    pa_sky = math.degrees(math.atan2(d[0] * -1 if False else d[0], d[1])) % 180.0      # (east component, north component) -> east of north
    out = dict(pixscale=pixscale, noise=sig, psf_sigma_pix=psf_sigma, lens_centre_pix=org, fit=p, evidence=best['evidence'], nfev=best['nfev'],
               chi2=r['chi2'], ndata=S.ndata, neff=neff, chi2_red=chi2red, lam=r['lam'], pa_sky_deg=pa_sky, seconds=time.time() - t0)
    # outer ring: subtract the lensed inner source and look at the residual radial profile
    full_mod = np.zeros_like(res)
    SIfull = X.SourceInversion(m, res, sig, np.ones_like(res, bool) & (r_as < 3.5), pixscale, org, (0.0, 0.0), kernel=k, grid=S.grid, n=S.grid['n'])
    r2 = SIfull.solve(r['lam'])
    _, mod2 = SIfull.images(r2)
    resid2 = res - mod2
    out['_resid2'] = resid2
    np.save('/workspace/work/stage6/j0946_resid2.npy', resid2)
    np.save('/workspace/work/stage6/j0946_model.npy', mod)
    np.save('/workspace/work/stage6/j0946_src.npy', src)
    return out


def regression_metrics(report=None):
    """Fast regression check (about 1 min): with the lens-light-subtracted residual (cached after the first call) and the stored best-fit lens
    (j0946_report.json, produced by `lens_real_j0946.py out.json analyse`), evaluate chi2_red, the evidence at the stored model and the evidence
    losses for theta_E +3 % and q -0.05 (must stay strongly negative, i.e. the data still prefer the stored lens)."""
    report = report or os.path.join(os.path.dirname(os.path.abspath(__file__)), 'j0946_report.json')
    rep = json.load(open(report))
    cache = os.path.join(CACHE, 'j0946_residual.npz')
    if not os.path.exists(cache):
        info, img, ll, res, org0, pixscale, psf_sigma, w = main(None)
        np.savez(cache, res=res, noise=info['noise'], org=info['lens_light_params'][:2], pixscale=pixscale, psf_sigma=psf_sigma)
    z = np.load(cache)
    res, sig, org, pixscale, psf_sigma = z['res'], float(z['noise']), tuple(float(a) for a in z['org']), float(z['pixscale']), float(z['psf_sigma'])
    ny, nx = res.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r_as = np.hypot(xx - org[0], yy - org[1]) * pixscale
    k = X.gaussian_kernel(psf_sigma)
    m1 = X.arc_mask(res, sig, org, 2.0 / pixscale, nsig=4.0, dilate=3) & (r_as > 0.55)

    def ev(p):
        S = X.SourceInversion(L.build_sie_shear(dict(p, x0=0.0, y0=0.0)), res, sig, m1, pixscale, org, (0.0, 0.0), kernel=k, n=40)
        r = S.best()
        return r['evidence'], r['chi2'] / (S.ndata - S.n_eff(r))
    p0 = rep['fit']
    e0, c0 = ev(p0)
    e1, _ = ev(dict(p0, theta_E=p0['theta_E'] * 1.03))
    e2, _ = ev(dict(p0, q=p0['q'] - 0.05))
    return dict(theta_E=p0['theta_E'], q=p0['q'], gamma=p0['gamma'], chi2_red=c0, dlnE_thetaE_p3pct=e1 - e0, dlnE_q_m005=e2 - e0, evidence=e0)


def PolyArea(x, y):
    x = np.asarray(x); y = np.asarray(y)
    return 0.5 * abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))


if __name__ == '__main__' and len(sys.argv) > 2 and sys.argv[2] == 'analyse':
    o = analyse(sys.argv[1])
    o.pop('_resid2')
    json.dump(o, open(sys.argv[1], 'w'), indent=1, default=float)
    print(json.dumps(o, default=float)[:1500])
