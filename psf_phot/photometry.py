"""PSF photometry for all sources in a catalog."""

import numpy as np
from .config import PSFPhotConfig
from .fitter import fit_psf_single


def do_psf_photometry(data, psf, sources, cfg=None):
    """Perform PSF photometry on all sources.

    Parameters
    ----------
    data : 2D array, image data
    psf : 2D array, PSF image
    sources : list of dicts with x, y, number
    cfg : PSFPhotConfig

    Returns
    -------
    list of result dicts
    """
    if cfg is None:
        cfg = PSFPhotConfig()

    ny_img, nx_img = data.shape
    psf_hy, psf_hx = psf.shape[0] // 2, psf.shape[1] // 2
    fit_r = cfg.fit_radius

    results = []
    for src in sources:
        x = src['x']
        y = src['y']
        num = src['number']

        # Extract cutout
        half = fit_r + psf_hx
        x0 = int(x) - half
        y0 = int(y) - half
        x1 = int(x) + half + 1
        y1 = int(y) + half + 1

        if x0 < 0 or y0 < 0 or x1 > nx_img or y1 > ny_img:
            results.append({
                'NUMBER': num, 'FLUX_PSF': np.nan, 'FLUXERR_PSF': np.nan,
                'MAG_PSF': 99.0, 'MAGERR_PSF': 99.0, 'CHI2_PSF': np.nan,
                'X_PSF': x + 1, 'Y_PSF': y + 1,
            })
            continue

        cutout = data[y0:y1, x0:x1].copy()

        # Trim PSF to match cutout if needed
        psf_cut = psf.copy()
        if psf_cut.shape[0] > cutout.shape[0] or psf_cut.shape[1] > cutout.shape[1]:
            cy, cx = psf_cut.shape[0] // 2, psf_cut.shape[1] // 2
            hy = min(cy, cutout.shape[0] // 2)
            hx = min(cx, cutout.shape[1] // 2)
            psf_cut = psf_cut[cy-hy:cy+hy+1, cx-hx:cx+hx+1]

        # Sub-pixel offset
        dx0 = x - int(x)
        dy0 = y - int(y)

        result = fit_psf_single(cutout, psf_cut, x0=dx0, y0=dy0,
                                 max_shift=cfg.max_shift)

        if result is None:
            results.append({
                'NUMBER': num, 'FLUX_PSF': np.nan, 'FLUXERR_PSF': np.nan,
                'MAG_PSF': 99.0, 'MAGERR_PSF': 99.0, 'CHI2_PSF': np.nan,
                'X_PSF': x + 1, 'Y_PSF': y + 1,
            })
            continue

        flux = result['flux']
        fluxerr = result['fluxerr']
        mag = -2.5 * np.log10(flux) + cfg.mag_zeropoint if flux > 0 else 99.0
        magerr = 1.0857 * fluxerr / flux if flux > 0 and np.isfinite(fluxerr) else 99.0

        results.append({
            'NUMBER': num,
            'FLUX_PSF': flux,
            'FLUXERR_PSF': fluxerr,
            'MAG_PSF': mag,
            'MAGERR_PSF': magerr,
            'CHI2_PSF': result['chi2'],
            'X_PSF': x + 1 + result['dx'],
            'Y_PSF': y + 1 + result['dy'],
        })

    return results
