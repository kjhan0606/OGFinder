"""2-component Bulge+Disk fitting engine."""

import numpy as np
from scipy.optimize import least_squares

from sersic_fit.fitter import fit_sersic
from .config import BulgeDiskConfig
from .models import bulge_disk_2d, compute_bt_ratio
from .psf_convolve import convolve_model


def decompose_source(cutout, x0, y0, a, b, theta, flux, psf=None, cfg=None):
    """Fit a 2-component Bulge+Disk model to a source cutout.

    Parameters
    ----------
    cutout : 2D array
        Image cutout around the source.
    x0, y0 : float
        Source center relative to cutout (pixels).
    a, b : float
        Semi-major and semi-minor axes (pixels).
    theta : float
        Position angle (radians).
    flux : float
        Total flux (initial estimate).
    psf : 2D array or None
        PSF kernel (normalized). If None, no convolution is applied.
    cfg : BulgeDiskConfig or None
        Configuration parameters.

    Returns
    -------
    dict with keys:
        bt_ratio, bulge_re, bulge_mag, disk_rs, disk_mag,
        chi2, flag (0=ok, 1=fit failed)
    """
    if cfg is None:
        cfg = BulgeDiskConfig()

    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]

    fail_result = {
        'bt_ratio': np.nan, 'bulge_re': np.nan, 'bulge_mag': np.nan,
        'disk_rs': np.nan, 'disk_mag': np.nan, 'chi2': np.nan, 'flag': 1,
    }

    # Step 1: single Sérsic fit for initial estimates
    from sersic_fit.config import SersicConfig
    scfg = SersicConfig(re_min=cfg.re_min, mag_zeropoint=cfg.mag_zeropoint)
    sersic_result = fit_sersic(cutout, x0, y0, a, b, theta, flux, scfg)

    if sersic_result is None:
        return fail_result

    n_single = sersic_result['n']
    re_single = sersic_result['re']
    Ie_single = sersic_result['Ie']
    ellip_single = sersic_result['ellip']
    theta_single = sersic_result['theta'] * np.pi / 180.0  # deg -> rad

    # Step 2: set initial params based on single Sérsic index
    bg0 = np.median(cutout[:3, :].ravel())

    if n_single > 2.5:
        # Bulge-dominated initial guess
        Ie_b0 = Ie_single * 0.7
        re_b0 = max(cfg.re_min, re_single * 0.6)
        e_b0 = min(0.5, ellip_single)
        theta_b0 = theta_single

        I0_d0 = Ie_single * 0.3
        rs_d0 = max(cfg.re_min, re_single * 1.5)
        e_d0 = min(0.7, ellip_single + 0.1)
        theta_d0 = theta_single
    else:
        # Disk-dominated initial guess
        Ie_b0 = Ie_single * 0.2
        re_b0 = max(cfg.re_min, re_single * 0.4)
        e_b0 = min(0.3, ellip_single * 0.5)
        theta_b0 = theta_single

        I0_d0 = Ie_single * 0.8
        rs_d0 = max(cfg.re_min, re_single * 1.2)
        e_d0 = min(0.8, ellip_single)
        theta_d0 = theta_single

    # Step 3: param vector [Ie_b, re_b, e_b, theta_b, I0_d, rs_d, e_d, theta_d, xc, yc, sky]
    p0 = [Ie_b0, re_b0, e_b0, theta_b0,
          I0_d0, rs_d0, e_d0, theta_d0,
          x0, y0, bg0]

    # Step 4: bounds
    re_max = max(a * 10, 100.0)
    bounds_lower = [0.0, cfg.re_min, 0.0, -np.pi,
                    0.0, cfg.re_min, 0.0, -np.pi,
                    x0 - 10, y0 - 10, -np.inf]
    bounds_upper = [np.inf, re_max, 0.95, np.pi,
                    np.inf, re_max, 0.95, np.pi,
                    x0 + 10, y0 + 10, np.inf]

    n_bulge = cfg.n_bulge
    n_disk = cfg.n_disk

    # Step 5: residuals function
    def residuals(params):
        model = bulge_disk_2d(params, xx, yy, n_bulge=n_bulge, n_disk=n_disk)
        if psf is not None:
            model = convolve_model(model, psf)
        return (cutout - model).ravel()

    # Step 6: optimize
    try:
        result = least_squares(residuals, p0,
                               bounds=(bounds_lower, bounds_upper),
                               method='trf', max_nfev=1000)
    except Exception:
        return fail_result

    popt = result.x
    npix = len(result.fun)
    npar = len(popt)
    chi2 = np.sum(result.fun ** 2) / max(1, npix - npar)

    # Step 7: compute B/T ratio
    bt = compute_bt_ratio(popt, xx, yy, n_bulge=n_bulge, n_disk=n_disk)

    # Extract component parameters
    Ie_b, re_b, e_b, theta_b, I0_d, rs_d, e_d, theta_d, xc, yc, sky = popt

    # Compute component fluxes for magnitudes
    from sersic_fit.utils import sersic_bn
    bn_b = sersic_bn(n_bulge)
    r_b = np.sqrt((xx - xc) ** 2 + (yy - yc) ** 2)
    exp_b = bn_b * ((r_b / max(re_b, 0.1)) ** (1.0 / max(n_bulge, 0.1)) - 1.0)
    exp_b = np.clip(exp_b, -50, 50)
    flux_bulge = float(np.sum(Ie_b * np.exp(-exp_b)))

    bn_d = sersic_bn(n_disk)
    exp_d = bn_d * ((r_b / max(rs_d, 0.1)) ** (1.0 / max(n_disk, 0.1)) - 1.0)
    exp_d = np.clip(exp_d, -50, 50)
    flux_disk = float(np.sum(I0_d * np.exp(-exp_d)))

    # Convert to magnitudes
    zp = cfg.mag_zeropoint
    bulge_mag = zp - 2.5 * np.log10(max(flux_bulge, 1e-30))
    disk_mag = zp - 2.5 * np.log10(max(flux_disk, 1e-30))

    return {
        'bt_ratio': bt,
        'bulge_re': re_b,
        'bulge_mag': bulge_mag,
        'disk_rs': rs_d,
        'disk_mag': disk_mag,
        'chi2': chi2,
        'flag': 0,
    }


def _fit_one_source(data, src, psf, cfg):
    """Extract cutout and fit Bulge+Disk for a single source.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    src : dict
        Source info with keys: num, x, y, a, b, theta, kron_r, flux.
    psf : 2D array or None
    cfg : BulgeDiskConfig

    Returns
    -------
    dict : decomposition result with 'num' added.
    """
    ny, nx = data.shape
    x, y = src['x'], src['y']
    a = src['a']
    kron_r = src['kron_r']

    halfsize = int(cfg.cutout_scale * kron_r * max(a, 3.0))
    halfsize = max(halfsize, 15)
    halfsize = min(halfsize, 150)

    ix, iy = int(x), int(y)
    x0 = max(0, ix - halfsize)
    y0 = max(0, iy - halfsize)
    x1 = min(nx, ix + halfsize + 1)
    y1 = min(ny, iy + halfsize + 1)

    nan_result = {
        'num': src['num'], 'bt_ratio': np.nan, 'bulge_re': np.nan,
        'bulge_mag': np.nan, 'disk_rs': np.nan, 'disk_mag': np.nan,
        'chi2': np.nan, 'flag': 1,
    }

    if x1 - x0 < 10 or y1 - y0 < 10:
        return nan_result

    cutout = data[y0:y1, x0:x1].copy()
    cx = x - x0
    cy = y - y0

    result = decompose_source(cutout, cx, cy, a, src['b'], src['theta'],
                               src['flux'], psf=psf, cfg=cfg)
    result['num'] = src['num']
    return result


def _worker_decompose(args):
    """Worker for parallel_map: decompose one source from shared data."""
    data_spec, src, psf, cfg = args
    shm = None
    try:
        data, shm = data_spec.attach()
        return _fit_one_source(data, src, psf, cfg)
    except Exception:
        return {
            'num': src.get('num', 0), 'bt_ratio': np.nan,
            'bulge_re': np.nan, 'bulge_mag': np.nan,
            'disk_rs': np.nan, 'disk_mag': np.nan,
            'chi2': np.nan, 'flag': 1,
        }
    finally:
        if shm is not None:
            shm.close()


def decompose_catalog(data, catalog, psf=None, cfg=None, n_workers=0):
    """Decompose all sources in a catalog.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    catalog : list of dict
        Each dict has keys: num, x, y, a, b, theta, kron_r, flux.
    psf : 2D array or None
        PSF kernel (normalized).
    cfg : BulgeDiskConfig or None
    n_workers : int
        0=auto, 1=sequential, N=N workers.

    Returns
    -------
    list of dict : One result per source.
    """
    if cfg is None:
        cfg = BulgeDiskConfig()

    from parallel import parallel_map, resolve_n_workers, SharedArray

    actual = resolve_n_workers(n_workers)

    if actual == 1:
        results = [_fit_one_source(data, src, psf, cfg) for src in catalog]
    else:
        with SharedArray(data) as spec:
            tasks = [(spec, src, psf, cfg) for src in catalog]
            results = parallel_map(_worker_decompose, tasks,
                                   n_workers=n_workers,
                                   label="Bulge+Disk")

    return results
