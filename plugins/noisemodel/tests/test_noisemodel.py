"""Background / noise model: truth-based validation on synthetic images with known correlated noise, gradients, scattered light, flat-field residuals."""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy.ndimage import gaussian_filter

from ogfkit import noise as nz, tsvio, imageio
import noisemodel as nm

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
HUDF = '/workspace/fits/hudf_f160w.fits'


def correlated_noise(shape, sigma_pix, sk, rng):
    """Gaussian noise convolved with a Gaussian kernel of sigma sk (pixel rms rescaled to sigma_pix): ACF is Gaussian with sigma sk*sqrt(2)."""
    w = rng.normal(size=shape)
    c = gaussian_filter(w, sk, mode='wrap') if sk > 0 else w
    return c * (sigma_pix / c.std())


def true_acf(sk, size=96):
    """Exact ACF (normalised) of white noise filtered with scipy's discrete gaussian_filter(sk): autocorrelation of its impulse response."""
    d = np.zeros((size, size)); d[size // 2, size // 2] = 1.0
    k = gaussian_filter(d, sk, mode='constant') if sk > 0 else d
    f = np.fft.rfft2(k)
    ac = np.fft.fftshift(np.fft.irfft2(np.abs(f) ** 2, s=k.shape))
    return ac / ac.max()


def aperture_weights(r, size=96, sub=5):
    """Fractional-pixel weights of a circular aperture of radius r centred on pixel (size//2, size//2) (what sep.sum_circle(subpix=5) integrates)."""
    u = (np.arange(size * sub) + 0.5) / sub - 0.5 - size // 2
    inside = np.hypot(u[:, None], u[None, :]) <= r
    return inside.reshape(size, sub, size, sub).mean(axis=(1, 3))


def true_aperture_sigma(r, sigma_pix, sk, size=96):
    """Analytic sigma of the sum in a circular aperture of radius r for that noise: sum_{ij} C(i-j) = sum_d C(d) * Mask-autocorrelation(d)."""
    m = aperture_weights(r, size)
    f = np.fft.rfft2(m)
    mac = np.fft.fftshift(np.fft.irfft2(np.abs(f) ** 2, s=m.shape))
    c = true_acf(sk, size)
    yy, xx = np.mgrid[:size, :size] - size // 2
    return sigma_pix * math.sqrt(float(np.sum(mac * c)))


class A:
    model = 'mesh'; bw = 64; fw = 3; poly_order = 2; n_aper = 700; seed = 3


def run(data, bad=None, radii=(1.5, 3, 5, 8), **kw):
    a = A()
    for k, v in kw.items():
        setattr(a, k, v)
    bad = np.zeros(data.shape, bool) if bad is None else bad
    return nm.analyse(data, bad, a, list(radii))


@pytest.fixture(scope='module')
def corr_img():
    rng = np.random.default_rng(11)
    return 100.0 + correlated_noise((900, 900), 5.0, 1.0, rng)


def test_correlation_recovered():
    for sk in (0.8, 1.2):
        d = 100 + correlated_noise((700, 700), 5.0, sk, np.random.default_rng(5))
        s, _ = run(d)
        ta = true_acf(sk); truth = float(ta[48, 49])
        print('sk %.1f rho1 %.4f truth %.4f' % (sk, s['acf']['rho1_x'], truth))
        assert abs(s['acf']['rho1_x'] - truth) < 0.02 and abs(s['acf']['rho1_y'] - truth) < 0.02
        assert abs(s['sigma_pix'] - 5.0) < 0.15


def test_noise_curve_vs_analytic_truth(corr_img):
    s, arr = run(corr_img)
    an = arr['an']
    errs = []
    for r, sig in zip(an['radii'], an['sigma']):
        t = true_aperture_sigma(r, 5.0, 1.0)
        errs.append(sig / t - 1)
        naive = an['sigma_pix'] * math.sqrt(math.pi * r * r)
    print('measured/true-1 per radius', np.round(errs, 3), 'beta %.3f alpha %.3f' % (s['noise_law']['beta'], s['noise_law']['alpha']))
    assert max(abs(e) for e in errs) < 0.10, errs
    # naive sqrt(N) underestimates strongly at large apertures
    assert s['correlation_factor'][-1] > 2.5
    # fitted law reproduces the truth within 10 % at all radii
    for r in (1.5, 3, 5, 8):
        N = np.pi * r * r
        law = nz.noise_law_sigma(N, s['noise_law'], s['sigma_pix'])
        assert abs(law / true_aperture_sigma(r, 5.0, 1.0) - 1) < 0.12


def test_white_noise_gives_sqrt_n():
    d = 100 + np.random.default_rng(2).normal(size=(700, 700)) * 5.0
    s, arr = run(d)
    print('white: beta %.3f alpha %.3f rho1 %.3f' % (s['noise_law']['beta'], s['noise_law']['alpha'], s['acf']['rho1']))
    assert abs(s['acf']['rho1']) < 0.03
    # fractional-pixel apertures: the truth is sigma * sqrt(sum w^2), not sigma sqrt(N)
    for r, sg in zip(arr['an']['radii'], arr['an']['sigma']):
        t = 5.0 * math.sqrt(float(np.sum(aperture_weights(r) ** 2)))
        assert abs(sg / t - 1) < 0.08, (r, sg, t)
    assert abs(s['noise_law']['beta'] - 0.5) < 0.1 and abs(s['noise_law']['alpha'] - 0.85) < 0.2      # fractional-pixel edges lower N_eff of small apertures


def star_catalog(shape, n, flux, rng, fwhm=3.0, sky_rms=5.0, margin=14):
    """Gaussian stars (sigma fwhm/2.355, integrated over pixels via oversampling) on a grid with margin; returns (image add-on, x, y, flux)."""
    img = np.zeros(shape)
    sg = fwhm / 2.355
    x = rng.uniform(margin, shape[1] - margin, n)
    y = rng.uniform(margin, shape[0] - margin, n)
    h = 12
    ax = np.arange(-h, h + 1)
    for xi, yi, fi in zip(x, y, flux):
        ix, iy = int(round(xi)), int(round(yi))
        gx = np.exp(-(ax[None, :] + ix - xi) ** 2 / (2 * sg * sg)); gy = np.exp(-(ax[:, None] + iy - yi) ** 2 / (2 * sg * sg))
        st = gx * gy
        img[iy - h:iy + h + 1, ix - h:ix + h + 1] += fi * st / st.sum()
    return img, x, y


def test_flux_error_pulls_naive_vs_corrected(tmp_path):
    rng = np.random.default_rng(21)
    shape = (900, 900)
    noise = correlated_noise(shape, 5.0, 1.0, rng)
    flux = np.full(400, 400.0)
    add, x, y = star_catalog(shape, 400, flux, rng)
    img = (100.0 + noise + add).astype(np.float32)
    imageio.save_fits(str(tmp_path / 'im.fits'), img)
    cat = tmp_path / 'cat.tsv'
    with open(cat, 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\n')
        for i, (a, b) in enumerate(zip(x, y)):
            f.write('%d\t%.4f\t%.4f\n' % (i + 1, a + 1, b + 1))
    out = subprocess.run([PY, os.path.join(os.path.dirname(HERE), 'noisemodel.py'), str(tmp_path / 'im.fits'), '--work', str(tmp_path / 'w'), '--catalog', str(cat), '--aperture', 'fixed',
                          '--aper-radius', '6', '--n-aper', '600', '--radii', '2,4,6,8'], capture_output=True, text=True, check=True).stdout
    lines = [l.split('\t') for l in out.strip().split('\n')]
    cols = lines[0]
    rows = [dict(zip(cols, l)) for l in lines[1:]]
    f = np.array([tsvio.fnum(r['NM_FLUX_AP']) for r in rows])
    e = np.array([tsvio.fnum(r['NM_FLUXERR']) for r in rows])
    corr = np.array([tsvio.fnum(r['NM_CORR']) for r in rows])
    # aperture r=6 holds >99.9 % of a sigma 1.27 px Gaussian
    pull = (f - 400.0) / e
    pull_naive = (f - 400.0) / (e / corr)
    sp = (f - 400.0).std()
    print('flux scatter %.2f ; corrected err %.2f ; naive err %.2f ; pull std corrected %.3f naive %.3f ; mean pull %.3f' % (sp, e.mean(), (e / corr).mean(), pull.std(), pull_naive.std(), pull.mean()))
    assert 0.85 < pull.std() < 1.2 and abs(pull.mean()) < 0.2
    assert pull_naive.std() > 2.0
    assert list(cols[:2]) == ['NUMBER', 'NM_SKY']


def test_gradient_and_scattered_light_recovery():
    rng = np.random.default_rng(31)
    n = 900
    yy, xx = np.mgrid[:n, :n].astype(float)
    grad = 0.004 * (xx - n / 2) + 0.0025 * (yy - n / 2)                      # +-2.2, +-1.1
    bump = 6.0 * np.exp(-((xx - 300) ** 2 + (yy - 600) ** 2) / (2 * 150.0 ** 2))   # scattered-light halo
    truth = 100.0 + grad + bump
    add, x, y = star_catalog((n, n), 120, np.full(120, 3000.0), rng)
    d = truth + add + correlated_noise((n, n), 5.0, 1.0, rng)
    res = {}
    for model in ('mesh', 'poly', 'mesh+poly'):
        s, arr = run(d, model=model, poly_order=2, radii=(2, 4))
        e = arr['bkg'] - truth
        res[model] = (float(np.sqrt(np.mean((e - np.median(e)) ** 2))), s['flatness_after']['peak_to_peak'])
    flat_raw = nz.flatness(d - np.median(d), np.zeros(d.shape, bool), 6)['peak_to_peak']
    print('bkg model rms error (counts), flatness p2p after:', {k: (round(v[0], 3), round(v[1], 3)) for k, v in res.items()}, 'raw p2p %.2f' % flat_raw)
    assert res['mesh'][0] < 0.8 and res['mesh+poly'][0] < 0.8
    assert res['poly'][0] > res['mesh'][0]                              # a quadratic cannot follow the halo
    assert res['mesh'][1] < 0.25 * flat_raw
    # a plane is recovered exactly by the polynomial model
    d2 = 100 + grad + correlated_noise((n, n), 5.0, 1.0, rng)
    s2, arr2 = run(d2, model='poly', poly_order=1, radii=(2, 4))
    e2 = arr2['bkg'] - (100 + grad)
    assert np.sqrt(np.mean((e2 - e2.mean()) ** 2)) < 0.05


def test_flatfield_residual_detected():
    rng = np.random.default_rng(41)
    n = 900
    yy, xx = np.mgrid[:n, :n].astype(float)
    flat = 1000.0 + 6.0 * np.sin(2 * np.pi * xx / 600.0) * np.sin(2 * np.pi * yy / 700.0 + 0.4)
    d_bad = flat + rng.normal(size=(n, n)) * 10
    d_ok = 1000.0 + rng.normal(size=(n, n)) * 10
    m = np.zeros((n, n), bool)
    fb, fo = nz.flatness(d_bad, m, 6), nz.flatness(d_ok, m, 6)
    print('flatness chi2 flat-field residual %.1f clean %.2f ; block p2p %.2f vs %.2f' % (fb['chi2_flat'], fo['chi2_flat'], fb['peak_to_peak'], fo['peak_to_peak']))
    assert fo['chi2_flat'] < 2.0 and fb['chi2_flat'] > 20
    s, _ = run(d_bad, radii=(2, 4))
    assert s['flatness_before']['chi2_flat'] > 20


def test_nonuniform_rms_map(tmp_path):
    rng = np.random.default_rng(51)
    n = 800
    xx = np.arange(n)[None, :] * np.ones((n, 1))
    sig = 4.0 + 4.0 * xx / (n - 1)                                  # exposure-time map: 4 -> 8
    d = 100 + rng.normal(size=(n, n)) * sig
    s, arr = run(d, radii=(2, 4))
    rr = arr['rms']
    ratio = rr[:, ::50].mean(axis=0) / sig[0, ::50]
    print('rms map / truth: min %.3f max %.3f' % (ratio.min(), ratio.max()))
    assert 0.93 < ratio.min() and ratio.max() < 1.07


@pytest.mark.skipif(not os.path.exists(HUDF), reason='HUDF image not available')
def test_real_drizzled_image():
    from astropy.io import fits
    d = fits.getdata(HUDF).astype(float)[1200:2000, 1200:2000]
    bad = ~np.isfinite(d) | (d == 0)
    s, arr = run(d, bad, radii=(1, 2, 3, 4, 6, 8))
    print('HUDF F160W crop: sigma_pix %.5f beta %.3f alpha %.3f rho1 %.3f fwhm %.2f ; corr factors %s' % (s['sigma_pix'], s['noise_law']['beta'], s['noise_law']['alpha'], s['acf']['rho1'],
          s['acf']['kernel_fwhm_px'], np.round(s['correlation_factor'], 2)))
    assert 0.55 < s['noise_law']['beta'] < 0.8 and s['acf']['rho1'] > 0.3 and s['correlation_factor'][-1] > 2.0
