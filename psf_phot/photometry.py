"""PSF photometry for all sources in a catalog."""

import numpy as np
from .config import PSFPhotConfig
from .fitter import fit_psf_single


def _fit_one_source(data, psf, src, cfg):
    """Fit PSF to a single source (used by both sequential and parallel paths)."""
    x = src['x']
    y = src['y']
    num = src['number']

    ny_img, nx_img = data.shape
    psf_hy, psf_hx = psf.shape[0] // 2, psf.shape[1] // 2
    fit_r = cfg.fit_radius

    half = fit_r + psf_hx
    x0 = int(x) - half
    y0 = int(y) - half
    x1 = int(x) + half + 1
    y1 = int(y) + half + 1

    nan_result = {
        'NUMBER': num, 'FLUX_PSF': np.nan, 'FLUXERR_PSF': np.nan,
        'MAG_PSF': 99.0, 'MAGERR_PSF': 99.0, 'CHI2_PSF': np.nan,
        'X_PSF': x + 1, 'Y_PSF': y + 1,
    }

    if x0 < 0 or y0 < 0 or x1 > nx_img or y1 > ny_img:
        return nan_result

    cutout = data[y0:y1, x0:x1].copy()

    psf_cut = psf.copy()
    if psf_cut.shape[0] > cutout.shape[0] or psf_cut.shape[1] > cutout.shape[1]:
        cy, cx = psf_cut.shape[0] // 2, psf_cut.shape[1] // 2
        hy = min(cy, cutout.shape[0] // 2)
        hx = min(cx, cutout.shape[1] // 2)
        psf_cut = psf_cut[cy-hy:cy+hy+1, cx-hx:cx+hx+1]

    dx0 = x - int(x)
    dy0 = y - int(y)

    result = fit_psf_single(cutout, psf_cut, x0=dx0, y0=dy0,
                             max_shift=cfg.max_shift)

    if result is None:
        return nan_result

    flux = result['flux']
    fluxerr = result['fluxerr']
    mag = -2.5 * np.log10(flux) + cfg.mag_zeropoint if flux > 0 else 99.0
    magerr = 1.0857 * fluxerr / flux if flux > 0 and np.isfinite(fluxerr) else 99.0

    return {
        'NUMBER': num,
        'FLUX_PSF': flux,
        'FLUXERR_PSF': fluxerr,
        'MAG_PSF': mag,
        'MAGERR_PSF': magerr,
        'CHI2_PSF': result['chi2'],
        'X_PSF': x + 1 + result['dx'],
        'Y_PSF': y + 1 + result['dy'],
    }


def _worker_fit_psf(args):
    """Worker for parallel_map: fit PSF to one source from shared data + PSF."""
    data_spec, psf_spec, src, cfg = args
    shm_data = None
    shm_psf = None
    try:
        data, shm_data = data_spec.attach()
        psf, shm_psf = psf_spec.attach()
        return _fit_one_source(data, psf, src, cfg)
    except Exception:
        return {
            'NUMBER': src.get('number', 0),
            'FLUX_PSF': np.nan, 'FLUXERR_PSF': np.nan,
            'MAG_PSF': 99.0, 'MAGERR_PSF': 99.0, 'CHI2_PSF': np.nan,
            'X_PSF': src.get('x', 0) + 1, 'Y_PSF': src.get('y', 0) + 1,
        }
    finally:
        if shm_data is not None:
            shm_data.close()
        if shm_psf is not None:
            shm_psf.close()


def do_psf_photometry(data, psf, sources, cfg=None, n_workers=0):
    """Perform PSF photometry on all sources.

    Parameters
    ----------
    data : 2D array, image data
    psf : 2D array, PSF image
    sources : list of dicts with x, y, number
    cfg : PSFPhotConfig
    n_workers : int
        0 = auto, 1 = sequential, N = N workers

    Returns
    -------
    list of result dicts
    """
    if cfg is None:
        cfg = PSFPhotConfig()

    from parallel import parallel_map, resolve_n_workers, SharedArray

    actual = resolve_n_workers(n_workers)

    if actual == 1:
        return [_fit_one_source(data, psf, src, cfg) for src in sources]

    with SharedArray(data) as data_spec, SharedArray(psf) as psf_spec:
        tasks = [(data_spec, psf_spec, src, cfg) for src in sources]
        results = parallel_map(_worker_fit_psf, tasks,
                               n_workers=n_workers, label="PSF Photometry")

    return results
