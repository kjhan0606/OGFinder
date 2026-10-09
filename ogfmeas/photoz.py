"""Template redshift by non-negative chi-square.

For each template and redshift the template is shifted to the observed frame
and scaled by one non-negative factor. This is the published template
chi-square idea. It is not a stellar-population code and not a substitute for
a full photometric-redshift pipeline.
"""
import numpy as np


def _scale(flux, model, invvar):
    numer = float(np.sum(flux * model * invvar))
    denom = float(np.sum(model * model * invvar))
    if denom <= 0:
        return 0.0
    return max(numer / denom, 0.0)


def template_chi2(wavelength, flux, fluxerr, templates, redshifts):
    """Fit observed samples against rest-frame templates.

    ``templates`` is a sequence of ``(rest_wavelength, rest_flux)``.
    Returns a dict with the minimum-chi-square redshift, the chi-square grid,
    a normalised pdf, the best template index, and its scale.
    """
    wave = np.asarray(wavelength, dtype=np.float64)
    observed = np.asarray(flux, dtype=np.float64)
    err = np.asarray(fluxerr, dtype=np.float64)
    if wave.shape != observed.shape or wave.shape != err.shape:
        raise ValueError("wavelength, flux, and fluxerr must share a shape")
    good = np.isfinite(wave) & np.isfinite(observed) & np.isfinite(err) & (err > 0)
    if int(good.sum()) < 3:
        raise ValueError("need at least three finite samples with positive errors")
    wave = wave[good]
    observed = observed[good]
    invvar = 1.0 / err[good] ** 2
    zs = np.asarray(redshifts, dtype=np.float64)
    ntemp = len(templates)
    chi2 = np.empty((ntemp, zs.size), dtype=np.float64)
    scales = np.empty_like(chi2)
    prepared = []
    for rest_wave, rest_flux in templates:
        prepared.append((np.asarray(rest_wave, dtype=np.float64),
                         np.asarray(rest_flux, dtype=np.float64)))
    for it, (rest_wave, rest_flux) in enumerate(prepared):
        order = np.argsort(rest_wave)
        rw = rest_wave[order]
        rf = rest_flux[order]
        for iz, z in enumerate(zs):
            shifted = rw * (1.0 + float(z))
            model = np.interp(wave, shifted, rf, left=0.0, right=0.0)
            scale = _scale(observed, model, invvar)
            residual = observed - scale * model
            chi2[it, iz] = float(np.sum(residual * residual * invvar))
            scales[it, iz] = scale
    flat = int(np.argmin(chi2))
    it_best, iz_best = divmod(flat, zs.size)
    best = chi2[it_best]
    delta = best - float(best.min())
    logw = -0.5 * delta
    logw -= float(logw.max())
    weights = np.exp(logw)
    weights_sum = float(weights.sum())
    pdf = weights / weights_sum if weights_sum > 0 else np.full(zs.size, 1.0 / zs.size)
    return {
        "z": float(zs[iz_best]),
        "template": int(it_best),
        "scale": float(scales[it_best, iz_best]),
        "zgrid": zs,
        "chi2": best,
        "chi2_all": chi2,
        "pdf": pdf,
    }
