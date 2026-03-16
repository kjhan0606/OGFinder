"""Multi-PSF fitting for source groups (NSTAR-style)."""

import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import shift


def fit_group(data, psf, sources, max_shift=2.0):
    """Simultaneously fit N PSFs to a group of sources.

    Parameters
    ----------
    data : 2D array, image data
    psf : 2D array, PSF image (normalized)
    sources : list of dicts with x, y, number
    max_shift : float, max centroid shift

    Returns
    -------
    list of result dicts
    """
    n = len(sources)
    if n == 0:
        return []

    ny, nx = data.shape
    psf_hy, psf_hx = psf.shape[0] // 2, psf.shape[1] // 2

    # Determine fitting region
    xs = [s['x'] for s in sources]
    ys = [s['y'] for s in sources]
    margin = psf_hx + 5

    x0 = max(0, int(min(xs)) - margin)
    y0 = max(0, int(min(ys)) - margin)
    x1 = min(nx, int(max(xs)) + margin + 1)
    y1 = min(ny, int(max(ys)) + margin + 1)

    cutout = data[y0:y1, x0:x1].copy()
    cut_ny, cut_nx = cutout.shape

    psf_sum = np.sum(psf)
    psf_norm = psf / psf_sum if psf_sum > 0 else psf

    # Initial parameters: [amp1, dx1, dy1, amp2, dx2, dy2, ..., bg]
    p0 = []
    bg0 = np.median(cutout[:3, :].ravel())
    for s in sources:
        amp0 = max(1.0, np.sum(cutout) / n - bg0 * cutout.size / n)
        p0.extend([amp0, 0.0, 0.0])
    p0.append(bg0)

    bounds_lower = []
    bounds_upper = []
    for _ in sources:
        bounds_lower.extend([0, -max_shift, -max_shift])
        bounds_upper.extend([np.inf, max_shift, max_shift])
    bounds_lower.append(-np.inf)
    bounds_upper.append(np.inf)

    def build_model(params):
        model = np.full_like(cutout, params[-1])  # background
        for i, s in enumerate(sources):
            amp = params[3 * i]
            dx = params[3 * i + 1]
            dy = params[3 * i + 2]

            # PSF position in cutout coords
            cx = s['x'] - x0 + dx
            cy = s['y'] - y0 + dy

            px0 = int(cx) - psf_hx
            py0 = int(cy) - psf_hy
            px1 = px0 + psf.shape[1]
            py1 = py0 + psf.shape[0]

            # Clip
            spx0 = max(0, -px0)
            spy0 = max(0, -py0)
            spx1 = psf.shape[1] - max(0, px1 - cut_nx)
            spy1 = psf.shape[0] - max(0, py1 - cut_ny)

            dpx0 = max(0, px0)
            dpy0 = max(0, py0)
            dpx1 = min(cut_nx, px1)
            dpy1 = min(cut_ny, py1)

            if dpx1 > dpx0 and dpy1 > dpy0:
                frac_dx = cx - int(cx)
                frac_dy = cy - int(cy)
                shifted_psf = shift(psf_norm, [frac_dy, frac_dx], order=1,
                                     mode='constant', cval=0)
                model[dpy0:dpy1, dpx0:dpx1] += amp * shifted_psf[spy0:spy1, spx0:spx1]

        return model

    def residuals(params):
        model = build_model(params)
        return (cutout - model).ravel()

    try:
        result = least_squares(residuals, p0,
                                bounds=(bounds_lower, bounds_upper),
                                method='trf', max_nfev=500)
    except Exception:
        return [{'NUMBER': s['number'], 'FLUX_CROWD': np.nan, 'FLUXERR_CROWD': np.nan,
                 'MAG_CROWD': 99.0, 'X_CROWD': s['x'] + 1, 'Y_CROWD': s['y'] + 1,
                 'N_NEIGHBORS': n - 1} for s in sources]

    chi2 = np.sum(result.fun**2) / max(1, len(result.fun) - len(p0))

    # Extract results
    results = []
    for i, s in enumerate(sources):
        amp = result.x[3 * i]
        dx = result.x[3 * i + 1]
        dy = result.x[3 * i + 2]
        flux = amp * psf_sum

        # Error from Jacobian
        try:
            J = result.jac
            cov = np.linalg.inv(J.T @ J) * chi2
            amp_err = np.sqrt(cov[3 * i, 3 * i])
            fluxerr = amp_err * psf_sum
        except (np.linalg.LinAlgError, ValueError):
            fluxerr = np.nan

        results.append({
            'NUMBER': s['number'],
            'FLUX_CROWD': flux,
            'FLUXERR_CROWD': fluxerr,
            'MAG_CROWD': 99.0,  # Computed by caller
            'X_CROWD': s['x'] + 1 + dx,
            'Y_CROWD': s['y'] + 1 + dy,
            'N_NEIGHBORS': n - 1,
        })

    return results
