"""TOY galaxy SED library for tests and the mock executables - NOT physical, NOT a model of real galaxies.

Five rest-frame f_nu shapes (power law + 4000 A break + simple IGM), a Calzetti-like power-law attenuation and a made-up mass-to-light normalisation per template.
It exists so that the adapter plumbing (band -> filter mapping, mag -> uJy conversion, native file formats, parsing, units) can be tested end to end with a
known truth.  The recovered numbers say nothing about the accuracy of EAZY / CIGALE / Bagpipes / Prospector.
"""
import math

import numpy as np

from . import filters

LAM = np.linspace(0.05, 8.0, 3000)        # micron (rest)
# name, slope of f_nu ~ lam^slope, 4000 A break factor, log age, log Z, log M/L normalisation
TEMPLATES = [('passive', 1.5, 3.5, 10.0, 0.0, 0.0), ('old', 1.0, 2.4, 9.6, -0.2, -0.3), ('intermediate', 0.5, 1.7, 9.0, -0.4, -0.6),
             ('starforming', 0.0, 1.15, 8.5, -0.5, -0.9), ('starburst', -0.4, 1.0, 7.8, -0.7, -1.2)]
NT = len(TEMPLATES)
AV_GRID = np.arange(0.0, 2.01, 0.25)


def rest_fnu(i):
    _, slope, brk, *_ = TEMPLATES[i]
    f = (LAM / 0.55) ** slope
    f = np.where(LAM < 0.4, f / brk, f) if brk != 1.0 else f
    # made-up spectral structure (two incommensurate ripples in ln lambda, template-specific phase) so that colours trace redshift in the broad bands
    x = np.log(LAM / 0.55)
    return f * (1 + 0.45 * np.sin(2 * np.pi * x / 0.9 + 1.3 * i)) * (1 + 0.3 * np.sin(2 * np.pi * x / 0.55 + 2.1 * i))


def igm_transmission(lam_rest, z):
    tau = 0.0036 * (1 + z) ** 3.46 * 0 + 0.0025 * (1 + z) ** 3.46        # mean Lya forest depth, crude
    t = np.where(lam_rest < 0.1216, np.exp(-min(tau * 3.0, 6.0)), 1.0)
    return np.where(lam_rest < 0.0912, 0.0, t)


def atten(lam_rest, av):
    return 10 ** (-0.4 * av * (lam_rest / 0.55) ** -1.0)


def band_flux(i, z, band, av=0.0):
    """Flux (uJy per unit template normalisation) in a registry band: Gaussian response with the registry FWHM."""
    e = filters.lookup(band) if isinstance(band, str) else band
    lam, dl = e['lam'], e['dl']
    lo, hi = lam - 2.5 * dl, lam + 2.5 * dl
    lo = max(lo, 0.06)
    lo = max(lo, LAM[0] * (1 + z))
    lam_o = np.linspace(lo, hi, 80)
    lr = lam_o / (1 + z)
    ok = (lr >= LAM[0]) & (lr <= LAM[-1])
    if ok.sum() < 3:
        return 0.0
    f = np.interp(lr, LAM, rest_fnu(i)) * igm_transmission(lr, z) * atten(lr, av)
    r = np.exp(-0.5 * ((lam_o - lam) / (dl / 2.355)) ** 2) / lam_o
    return float(np.trapezoid(f * r, lam_o) / np.trapezoid(r, lam_o)) if hasattr(np, 'trapezoid') else float(np.trapz(f * r, lam_o) / np.trapz(r, lam_o))


_cache = {}


def design(z, bands, av=0.0):
    """(nbands, NT) matrix of template fluxes for the band list at redshift z (cached on a rounded z)."""
    key = (round(z, 4), tuple(bands), av)
    if key not in _cache:
        _cache[key] = np.array([[band_flux(i, z, b, av) for i in range(NT)] for b in bands])
        if len(_cache) > 20000:
            _cache.clear()
    return _cache[key]


def dl_mpc(z):
    c = 299792.458
    from math import sqrt
    zz = np.linspace(0, z, 200)
    ez = np.sqrt(0.3 * (1 + zz) ** 3 + 0.7)
    dc = c / 70.0 * (np.trapezoid(1 / ez, zz) if hasattr(np, 'trapezoid') else np.trapz(1 / ez, zz))
    return (1 + z) * dc


def log_mass(i, scale, z):
    """Toy stellar mass: scale (template normalisation, uJy at rest 0.55 um) x d_L^2 x template M/L."""
    d = max(dl_mpc(max(z, 0.01)), 1.0)
    return math.log10(max(scale, 1e-30) * d * d / (1 + z)) + 5.0 + TEMPLATES[i][5]


def simulate(i, z, mag_ref, bands, ref_band, av=0.0, rng=None, noise_frac=0.05, floor_ujy=0.02):
    """Photometry (uJy, errors) of template i at z normalised to magnitude mag_ref in ref_band."""
    f = np.array([band_flux(i, z, b, av) for b in bands])
    fr = band_flux(i, z, ref_band, av)
    scale = 10 ** (-0.4 * (mag_ref - 23.9)) / max(fr, 1e-30)
    fl = f * scale
    err = np.maximum(noise_frac * fl, floor_ujy)
    if rng is not None:
        fl = fl + rng.normal(size=fl.size) * err
    return fl, err, scale


def fit_photoz(flux, err, bands, zgrid, combine=False):
    """chi2(z) on a redshift grid.  combine=False: best single template with a free scale (the mock default - the TOY templates are too featureless for linear
    combinations); combine=True: EAZY-like non-negative linear combination of all templates."""
    from scipy.optimize import nnls
    chi = np.zeros(len(zgrid))
    coefs = np.zeros((len(zgrid), NT))
    w = 1.0 / err
    for k, z in enumerate(zgrid):
        D = design(z, bands)
        if combine:
            c, rn = nnls(D * w[:, None], flux * w)
            coefs[k] = c
            chi[k] = rn * rn
        else:
            best = None
            for i in range(NT):
                a = D[:, i] * w
                s_ = max(float(np.dot(a, flux * w) / max(np.dot(a, a), 1e-30)), 0.0)
                c2 = float(np.sum((flux * w - s_ * a) ** 2))
                if best is None or c2 < best[0]:
                    best = (c2, i, s_)
            chi[k] = best[0]
            coefs[k, best[1]] = best[2]
    return chi, coefs


def fit_sed(flux, err, bands, z, rng=None):
    """Grid over template x Av at fixed z with analytic best scale; returns posterior samples summary."""
    w = 1.0 / err
    rows = []
    for ai, av in enumerate(AV_GRID):
        D = design(z, bands, av)
        for i in range(NT):
            a = D[:, i] * w
            s = float(np.dot(a, flux * w) / max(np.dot(a, a), 1e-30))
            s = max(s, 0.0)
            chi = float(np.sum((flux * w - s * a) ** 2))
            rows.append((chi, i, av, s))
    chi = np.array([r[0] for r in rows])
    p = np.exp(-0.5 * (chi - chi.min()))
    p /= p.sum()
    lm = np.array([log_mass(r[1], r[3], z) for r in rows])
    la = np.array([TEMPLATES[r[1]][3] for r in rows]); lz = np.array([TEMPLATES[r[1]][4] for r in rows]); av = np.array([r[2] for r in rows])
    def mean_sd(x):
        m = float(np.sum(p * x)); return m, float(math.sqrt(max(np.sum(p * (x - m) ** 2), 0)))
    out = {}
    out['LOG_MASS'], out['LOG_MASS_ERR'] = mean_sd(lm)
    out['LOG_AGE'], out['LOG_AGE_ERR'] = mean_sd(la)
    out['LOG_Z'], _ = mean_sd(lz)
    out['AV'], _ = mean_sd(av)
    sfr_t = np.array([10 ** (r_la - 9.0) * 0 + 1.0 for r_la in la])
    mass = 10 ** lm
    out['SFR'] = float(np.sum(p * mass * 10 ** (-0.7 * (la - 7.8) - 9.0)))       # toy: specific SFR falls with age
    out['SED_CHI2'] = float(chi.min() / max(len(flux) - 3, 1))
    return out
