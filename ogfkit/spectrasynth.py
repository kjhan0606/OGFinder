"""Synthetic spectra / 2D spectra / IFU cubes with known truth (validation of ogfkit.spectra)."""
import math

import numpy as np

from .spectra import LINES, SQRT2PI

STRONG = {n: (l, w) for n, l, w in LINES}


def spectrum_1d(z=0.8, seed=1, wave=(5000.0, 10000.0), dlam=2.5, cont=1.0, slope=-0.15, noise=0.2, lines=None, fwhm_kms=200.0, res_fwhm_A=6.0, ha_flux=None, scale=1.0):
    """Emission-line galaxy spectrum: power-law-ish continuum plus Gaussian lines at (1+z) rest wavelengths.  lines = dict name -> flux (default: strength x 20 x scale
    x noise x dlam units); line sigma = sqrt(intrinsic^2 + instrument^2).  Returns (wave, flux, err, truth dict)."""
    rng = np.random.default_rng(seed)
    w = np.arange(wave[0], wave[1], dlam)
    c = cont * (w / 7500.0) ** slope
    f = c.copy()
    truth = dict(z=z, lines={})
    if lines is None:
        lines = {n: wt * 12.0 * noise * dlam * 4 * scale for n, l, wt in LINES if wt >= 0.3}
    for n, F in lines.items():
        lam = STRONG[n][0] * (1 + z)
        if not (w[0] + 20 < lam < w[-1] - 20):
            continue
        s = math.hypot(fwhm_kms / 299792.458 * lam, res_fwhm_A) / 2.3548
        f += F / (s * SQRT2PI) * np.exp(-0.5 * ((w - lam) / s) ** 2)
        truth['lines'][n] = dict(flux=F, obs=lam, sigma=s)
    err = np.full_like(w, noise)
    f = f + rng.normal(0, noise, len(w))
    return w, f, err, truth


def spectrum_2d(ny=41, row=20.3, z=0.45, seed=1, dlam=3.0, wave0=6000.0, nw=1200, noise=0.1, psf_sigma=1.6, amp=3.0, ha=True):
    """2D spectrum: continuum trace + [OIII]/Ha lines (point source, Gaussian spatial profile sigma psf_sigma at fractional row), slowly varying sky residual per column.
    Returns (wave, data[ny, nw], err, truth dict with the 1D flux of the source)."""
    rng = np.random.default_rng(seed)
    w = wave0 + dlam * np.arange(nw)
    prof = np.exp(-0.5 * ((np.arange(ny) - row) / psf_sigma) ** 2)
    prof /= prof.sum()
    f1 = amp * (w / 8000.0) ** -0.5
    lines = {'[OIII]5008': 400.0, 'Ha': 500.0, '[OII]': 250.0}
    truth = dict(row=row, z=z, lines={})
    for n, F in lines.items():
        lam = STRONG[n][0] * (1 + z)
        if w[0] < lam < w[-1]:
            s = 8.0 / 2.3548 * 1.0
            f1 = f1 + F / (s * SQRT2PI) * np.exp(-0.5 * ((w - lam) / s) ** 2)
            truth['lines'][n] = dict(flux=F, obs=lam)
    data = prof[:, None] * f1[None, :] + rng.normal(0, noise, (ny, nw))
    sky = 0.5 * np.sin(w / 300.0)
    data = data + sky[None, :]
    err = np.full_like(data, noise)
    truth['flux1d'] = f1
    return w, data, err, truth


def cube(nx=21, ny=21, nw=240, wave0=6500.0, dlam=1.5, seed=1, noise=0.05, vmax=150.0, re=4.0, z=0.0, inc_pa=30.0, line_rest=6564.61, flux0=40.0, extra_lines=True):
    """IFU cube with a rotating exponential disc in one emission line (arctan rotation curve, amplitude vmax km/s) plus a flat continuum source.
    Returns (wave, cube[nw, ny, nx], truth dict with velocity map, flux map)."""
    rng = np.random.default_rng(seed)
    w = wave0 + dlam * np.arange(nw)
    yy, xx = np.mgrid[0:ny, 0:nx]
    x0, y0 = (nx - 1) / 2, (ny - 1) / 2
    ph = math.radians(inc_pa)
    dx, dy = xx - x0, yy - y0
    xr = dx * math.cos(ph) + dy * math.sin(ph)
    r = np.hypot(dx, dy)
    flux = flux0 * np.exp(-r / re)
    cosphi = np.where(r > 0, xr / np.maximum(r, 1e-9), 0.0)
    vel = vmax * (2 / math.pi) * np.arctan(r / 2.0) * cosphi
    sig_v = 40.0
    lam_obs = line_rest * (1 + z)
    c = 299792.458
    data = np.zeros((nw, ny, nx))
    for j in range(ny):
        for i in range(nx):
            mu = lam_obs * (1 + vel[j, i] / c)
            s = math.hypot(sig_v / c * lam_obs, 1.5)
            data[:, j, i] = 0.3 * np.exp(-r[j, i] / (2 * re))
            for rest, ratio in (((line_rest, 1.0), (6585.27, 0.35), (6718.29, 0.2), (6732.67, 0.15)) if extra_lines else ((line_rest, 1.0),)):
                m_ = rest * (1 + z) * (1 + vel[j, i] / c)
                data[:, j, i] += ratio * flux[j, i] / (s * SQRT2PI) * np.exp(-0.5 * ((w - m_) / s) ** 2)
    data += rng.normal(0, noise, data.shape)
    return w, data, dict(velocity=vel, flux=flux, lam_obs=lam_obs, center=(x0, y0), sigma_noise=noise)
