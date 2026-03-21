"""LSBG-specific Sérsic profile fitting.

1D radial profile extraction → 1D Sérsic fit (initial guess) → 2D Sérsic fit.
Parallelized per-source via SharedArray + parallel_map.
"""

import sys
import numpy as np
from scipy.optimize import least_squares, curve_fit
from scipy.special import gamma as gamma_func

from sersic_fit.utils import sersic_bn, ellipse_radius
from sersic_fit.models import sersic_2d


def sersic_total_flux(Ie, re, n, ellip):
    """Compute analytical total flux of a 2D Sérsic profile.

    F_total = 2*pi * Ie * re^2 * (1 - ellip) * n * exp(bn) * Gamma(2n) / bn^(2n)

    Parameters
    ----------
    Ie : float
        Surface brightness at effective radius.
    re : float
        Effective radius (pixels).
    n : float
        Sérsic index.
    ellip : float
        Ellipticity (1 - b/a).

    Returns
    -------
    float : total flux (same units as Ie * pixels^2).
    """
    if n <= 0 or re <= 0 or Ie <= 0:
        return np.nan
    bn = sersic_bn(n)
    q = 1.0 - ellip  # axis ratio b/a
    try:
        ftot = (2.0 * np.pi * Ie * re**2 * q * n
                * np.exp(bn) * gamma_func(2.0 * n) / bn**(2.0 * n))
    except (OverflowError, ValueError):
        return np.nan
    if not np.isfinite(ftot) or ftot <= 0:
        return np.nan
    return float(ftot)


def extract_radial_profile(cutout, xc, yc, ellip, theta, n_bins=20):
    """Extract azimuthally-averaged radial surface brightness profile.

    Uses elliptical annuli for non-circular sources.

    Parameters
    ----------
    cutout : 2D array
        Source cutout (background-subtracted).
    xc, yc : float
        Center in cutout coordinates.
    ellip : float
        Ellipticity (1 - b/a).
    theta : float
        Position angle (radians).
    n_bins : int
        Number of radial bins.

    Returns
    -------
    r_mid : 1D array
        Midpoint radius of each bin (pixels).
    sb : 1D array
        Mean surface brightness in each bin.
    sb_err : 1D array
        Error in surface brightness (std / sqrt(N)).
    """
    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    r = ellipse_radius(xx.astype(float), yy.astype(float),
                       xc, yc, ellip, theta)

    r_max = min(nx / 2, ny / 2, r.max())
    bin_edges = np.linspace(0, r_max, n_bins + 1)
    r_mid = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    sb = np.zeros(n_bins)
    sb_err = np.zeros(n_bins)

    for i in range(n_bins):
        ring = (r >= bin_edges[i]) & (r < bin_edges[i + 1])
        npix = ring.sum()
        if npix > 0:
            vals = cutout[ring]
            sb[i] = np.mean(vals)
            sb_err[i] = np.std(vals) / max(1, np.sqrt(npix))
        else:
            sb[i] = 0.0
            sb_err[i] = 1e10

    return r_mid, sb, sb_err


def sersic_1d(r, Ie, re, n):
    """1D Sérsic profile: I(r) = Ie * exp(-bn * ((r/re)^(1/n) - 1))."""
    bn = sersic_bn(n)
    exponent = bn * ((r / max(re, 0.1)) ** (1.0 / max(n, 0.1)) - 1.0)
    exponent = np.clip(exponent, -50, 50)
    return Ie * np.exp(-exponent)


def fit_sersic_1d(r_mid, sb, sb_err):
    """Fit 1D Sérsic profile to radial surface brightness profile.

    Returns
    -------
    dict with Ie, re, n, or None on failure.
    """
    # Initial guesses
    Ie0 = max(sb[0], 1e-10)
    # Find half-light radius: where sb drops to ~Ie/e
    half_idx = np.searchsorted(-sb, -Ie0 * 0.5)
    re0 = r_mid[min(half_idx, len(r_mid) - 1)]
    re0 = max(re0, 1.0)

    # Only fit bins with positive signal
    good = (sb > 0) & (sb_err > 0) & (sb_err < 1e9)
    if good.sum() < 4:
        return None

    r_fit = r_mid[good]
    sb_fit = sb[good]
    w = 1.0 / np.maximum(sb_err[good], 1e-10)

    try:
        popt, _ = curve_fit(
            sersic_1d, r_fit, sb_fit,
            p0=[Ie0, re0, 1.0],
            bounds=([0, 0.5, 0.2], [np.inf, r_mid[-1] * 2, 10.0]),
            sigma=1.0 / w, absolute_sigma=True,
            maxfev=200
        )
        return {'Ie': popt[0], 're': popt[1], 'n': popt[2]}
    except Exception:
        return None


def fit_sersic_2d(cutout, xc, yc, a, b, theta, flux, config):
    """Fit 2D Sérsic profile with 1D initial guess.

    Steps:
    1. Extract radial profile in elliptical annuli
    2. Fit 1D Sérsic for initial (Ie, re, n)
    3. Fit full 2D Sérsic model (C+OpenMP fast path when available)

    Parameters
    ----------
    cutout : 2D array
        Source cutout (background-subtracted).
    xc, yc : float
        Center in cutout coordinates.
    a, b : float
        Semi-major/minor axes.
    theta : float
        Position angle (radians).
    flux : float
        Total flux (for initial guess fallback).
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    dict with n, re, Ie, ellip, theta_deg, chi2, grade, or None.
    """
    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]

    ellip = max(0.0, min(0.9, 1.0 - b / max(a, 0.5)))
    bg0 = np.median(np.concatenate([
        cutout[0, :], cutout[-1, :], cutout[:, 0], cutout[:, -1]
    ]))

    # Step 1: 1D radial profile
    cutout_sub = cutout - bg0
    r_mid, sb, sb_err = extract_radial_profile(
        cutout_sub, xc, yc, ellip, theta, n_bins=20
    )

    # Step 2: 1D Sérsic fit for initial guess
    fit1d = fit_sersic_1d(r_mid, sb, sb_err)
    if fit1d is not None:
        Ie0 = fit1d['Ie']
        re0 = fit1d['re']
        n0 = fit1d['n']
    else:
        # Fallback: simple estimates
        re0 = max(config.sersic_re_min, np.sqrt(a * b))
        Ie0 = max(1.0, flux / (2.0 * np.pi * re0 * re0 * (1.0 - ellip)))
        n0 = 1.0

    # Step 3: 2D Sérsic fit
    p0 = [Ie0, re0, n0, xc, yc, ellip, theta, bg0]

    bounds_lower = [
        0, config.sersic_re_min, config.sersic_n_min,
        xc - 5, yc - 5, 0.0, -np.pi, -np.inf
    ]
    bounds_upper = [
        np.inf, max(a * 5, 50), config.sersic_n_max,
        xc + 5, yc + 5, 0.95, np.pi, np.inf
    ]

    # Check for C fast path: full LM in C
    try:
        from sersic_fit.fast import is_available, fit_sersic_2d_c
        use_fast = is_available()
    except ImportError:
        use_fast = False

    if use_fast:
        xx_flat = np.ascontiguousarray(xx.ravel(), dtype=np.float64)
        yy_flat = np.ascontiguousarray(yy.ravel(), dtype=np.float64)
        data_flat = np.ascontiguousarray(cutout.ravel(), dtype=np.float64)

        try:
            params_out, chi2 = fit_sersic_2d_c(
                data_flat, xx_flat, yy_flat, p0,
                bounds_lower, bounds_upper,
                max_iter=config.sersic_max_nfev
            )
        except Exception:
            return None

        Ie, re, n, xc_f, yc_f, ellip_f, theta_f, bg_f = params_out
    else:
        def residuals(params):
            model = sersic_2d(params, xx, yy)
            return (cutout - model).ravel()

        try:
            result = least_squares(
                residuals, p0,
                bounds=(bounds_lower, bounds_upper),
                method='trf',
                max_nfev=config.sersic_max_nfev
            )
        except Exception:
            return None

        Ie, re, n, xc_f, yc_f, ellip_f, theta_f, bg_f = result.x
        npix_total = len(result.fun)
        chi2 = np.sum(result.fun ** 2) / max(1, npix_total - 8)

    # Analytical total flux from Sérsic profile
    flux_total = sersic_total_flux(Ie, re, n, ellip_f)

    return {
        'n': float(n),
        're': float(re),
        'Ie': float(Ie),
        'ellip': float(ellip_f),
        'theta_deg': float(np.degrees(theta_f)),
        'chi2': float(chi2),
        'bg': float(bg_f),
        'flux_total': flux_total,
    }


def _extract_cutout(data, x, y, a, b, scale):
    """Extract a square cutout centered on source."""
    ny, nx = data.shape
    r = max(a, b) * scale
    r = max(r, 50)  # minimum 50 pixels (13" at 0.262"/px)
    r = min(r, 200)  # maximum 200 pixels to keep 2D fit tractable
    half = int(r)

    ix, iy = int(round(x)), int(round(y))
    y0 = max(0, iy - half)
    y1 = min(ny, iy + half + 1)
    x0 = max(0, ix - half)
    x1 = min(nx, ix + half + 1)

    cutout = data[y0:y1, x0:x1].copy()
    xc_cut = x - x0
    yc_cut = y - y0

    return cutout, xc_cut, yc_cut


def _worker_sersic(args):
    """Parallel worker for Sérsic fitting."""
    import os
    from parallel import SharedArraySpec

    # Limit OpenMP threads inside multiprocessing workers
    if 'OMP_NUM_THREADS' not in os.environ:
        os.environ['OMP_NUM_THREADS'] = '2'

    data_spec, obj_dict, config = args

    data_arr, data_shm = data_spec.attach()

    try:
        x = obj_dict['x']
        y = obj_dict['y']
        a = obj_dict['a']
        b = obj_dict['b']
        theta = obj_dict['theta']
        flux = obj_dict['flux']

        cutout, xc, yc = _extract_cutout(
            data_arr, x, y, a, b, config.sersic_cutout_scale
        )

        if cutout.size < 25:
            return None

        return fit_sersic_2d(cutout, xc, yc, a, b, theta, flux, config)
    finally:
        data_shm.close()


def fit_sources_sersic(data, objects, config):
    """Fit 2D Sérsic profiles to all detected LSBG candidates.

    Parallelized using SharedArray + parallel_map.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    objects : structured array
        SEP detection catalog.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    results : list of dict or None
        Sérsic fit results for each source (None if fit failed).
    """
    from parallel import SharedArray, parallel_map, resolve_n_workers

    n_sources = len(objects)
    if n_sources == 0:
        return []

    actual_workers = resolve_n_workers(config.n_workers)

    # Sequential mode
    if actual_workers <= 1 or n_sources < 10:
        results = []
        for i in range(n_sources):
            obj = objects[i]
            x = float(obj['x'])
            y = float(obj['y'])
            a = float(obj['a'])
            b = float(obj['b'])
            theta = float(obj['theta'])
            flux = float(obj['flux'])

            cutout, xc, yc = _extract_cutout(
                data, x, y, a, b, config.sersic_cutout_scale
            )

            if cutout.size < 25:
                results.append(None)
                continue

            fit = fit_sersic_2d(cutout, xc, yc, a, b, theta, flux, config)
            results.append(fit)

            if (i + 1) % 50 == 0:
                n_ok = sum(1 for r in results if r is not None)
                print(f"  Sérsic fit: {i + 1}/{n_sources} "
                      f"({n_ok} succeeded)", file=sys.stderr)

        n_ok = sum(1 for r in results if r is not None)
        print(f"Sérsic fit: {n_ok}/{n_sources} succeeded", file=sys.stderr)
        return results

    # Parallel mode
    print(f"Sérsic fit: {n_sources} sources, {actual_workers} workers",
          file=sys.stderr)

    data_c = np.ascontiguousarray(data, dtype=np.float64)

    with SharedArray(data_c) as data_spec:
        tasks = [
            (data_spec,
             {name: float(objects[i][name]) if objects.dtype[name] != np.int32
              else int(objects[i][name])
              for name in objects.dtype.names},
             config)
            for i in range(n_sources)
        ]
        results = parallel_map(_worker_sersic, tasks,
                               n_workers=config.n_workers,
                               label="Sérsic fit")

    n_ok = sum(1 for r in results if r is not None)
    print(f"Sérsic fit: {n_ok}/{n_sources} succeeded", file=sys.stderr)

    return results
