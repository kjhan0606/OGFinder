"""Synthetic trails and test fields for ogfkit.trails (validation, tests, GUI check)."""
import math
import numpy as np


def trail_image(shape, theta, rho, fwhm_box, amp, s0=None, s1=None, blur=1.2, gradient=0.0):
    """additive trail: a box of full width `fwhm_box` (px) convolved with a Gaussian of sigma `blur`, peak amplitude `amp` per pixel, along the line
    (theta deg from +x, rho px from the image centre; same convention as ogfkit.trails).  s0, s1 = along-track extent (None = whole line).
    `gradient` makes the amplitude vary linearly along the track by that fraction (+-gradient/2)."""
    from scipy.special import erf
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float64)
    xr, yr = xx - (nx - 1) / 2.0, yy - (ny - 1) / 2.0
    th = math.radians(theta)
    t = -xr * math.sin(th) + yr * math.cos(th) - rho
    s = xr * math.cos(th) + yr * math.sin(th)
    a = fwhm_box / 2.0
    prof = 0.5 * (erf((t + a) / (math.sqrt(2) * blur)) - erf((t - a) / (math.sqrt(2) * blur)))
    m = np.ones(shape)
    if s0 is not None:
        m *= 0.5 * (1 + erf((s - s0) / (math.sqrt(2) * blur)))
    if s1 is not None:
        m *= 0.5 * (1 + erf((s1 - s) / (math.sqrt(2) * blur)))
    if gradient:
        L = math.hypot(nx, ny)
        m *= 1 + gradient * (s / L)
    return (amp * prof * m).astype(np.float32)


def test_field(n=900, seed=3, sky=100.0, sigma=3.0, n_gal=14, n_star=60, corr=0.0):
    """noise + Gaussian-smoothed noise (corr) + stars + Sersic-like galaxies; returns image, truth catalogue (x, y, flux, kind, size)"""
    from scipy.ndimage import gaussian_filter
    rng = np.random.default_rng(seed)
    img = rng.normal(0, 1, (n, n))
    if corr > 0:
        img = gaussian_filter(img, corr)
        img /= img.std()
    img = sky + sigma * img
    yy, xx = np.mgrid[0:n, 0:n].astype(np.float64)
    truth = []
    for i in range(n_star):
        x, y = rng.uniform(10, n - 10, 2)
        f = 10 ** rng.uniform(2.2, 3.8)
        sg = 1.6
        r2 = (xx - x) ** 2 + (yy - y) ** 2
        img += f / (2 * math.pi * sg ** 2) * np.exp(-0.5 * r2 / sg ** 2)
        truth.append((x, y, f, 'star', sg))
    for i in range(n_gal):
        x, y = rng.uniform(40, n - 40, 2)
        f = 10 ** rng.uniform(3.6, 4.6)
        re = rng.uniform(4, 16)
        q = rng.uniform(0.4, 1.0); ph = rng.uniform(0, math.pi)
        c, s = math.cos(ph), math.sin(ph)
        u = (xx - x) * c + (yy - y) * s; v = -(xx - x) * s + (yy - y) * c
        r = np.sqrt(u ** 2 + (v / q) ** 2)
        prof = np.exp(-1.678 * ((r / re) - 1))               # exponential (n=1)
        prof *= f / prof.sum()
        img += gaussian_filter(prof, 1.6)
        truth.append((x, y, f, 'galaxy', re))
    return img.astype(np.float32), truth


def curved_trail_image(shape, p0, theta0, curv_per_100, length, fwhm_box, amp, blur=1.2, duty=None, period=60.0, phase=0.0):
    """additive curved trail: starts at p0 (x, y 0-based) with heading theta0 (deg from +x toward +y), the heading changes by `curv_per_100` degrees per 100 px;
    `duty` (0-1) makes it flicker: on for duty*period px, off for the rest of each `period`.  Same profile as trail_image."""
    from scipy.special import erf
    from scipy.spatial import cKDTree
    ny, nx = shape
    n = int(length * 4)
    s = np.arange(n) / 4.0
    th = np.radians(theta0 + curv_per_100 * s / 100.0)
    # integrate the heading
    x = p0[0] + np.cumsum(np.cos(th)) / 4.0
    y = p0[1] + np.cumsum(np.sin(th)) / 4.0
    tree = cKDTree(np.c_[x, y])
    yy, xx = np.mgrid[0:ny, 0:nx]
    a = fwhm_box / 2.0
    near = (np.zeros(shape, bool))
    pts = np.c_[xx.ravel(), yy.ravel()].astype(float)
    d, idx = tree.query(pts, distance_upper_bound=a + 8.0)
    ok = np.isfinite(d)
    out = np.zeros(ny * nx)
    prof = 0.5 * (erf((d[ok] + a) / (math.sqrt(2) * blur)) - erf((d[ok] - a) / (math.sqrt(2) * blur)))
    on = np.ones(ok.sum())
    if duty is not None:
        on = (((s[idx[ok]] + phase) % period) < duty * period).astype(float)
    out[ok] = amp * prof * on
    return out.reshape(shape).astype(np.float32), (x, y, s)
