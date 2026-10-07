"""ICL fraction and isophotal radius measurements."""

import numpy as np


def compute_icl_fraction(profile, mu_threshold=26.5):
    """Compute ICL fraction from surface brightness profile.

    f_ICL = L_ICL / (L_ICL + L_gal)
    where L_ICL = sum of flux in annuli with mu > mu_threshold
    and   L_gal = sum of flux in annuli with mu <= mu_threshold.

    Parameters
    ----------
    profile : list of dict
        SB profile from measure_sb_profile.
    mu_threshold : float
        Surface brightness threshold separating ICL from galaxies
        (mag/arcsec^2).

    Returns
    -------
    f_icl : float
        ICL fraction (0 to 1).
    l_icl : float
        Total ICL flux.
    l_gal : float
        Total galaxy flux.
    """
    l_icl = 0.0
    l_gal = 0.0

    for p in profile:
        mu = p['MU']
        flux = p['FLUX']
        if mu >= 99.0 or flux <= 0:
            continue
        if mu > mu_threshold:
            l_icl += flux
        else:
            l_gal += flux

    l_total = l_icl + l_gal
    f_icl = l_icl / l_total if l_total > 0 else 0.0

    return f_icl, l_icl, l_gal


def compute_icl_fraction_multi(profile, mu_min=25.0, mu_max=28.0,
                                mu_step=0.5, n_bootstrap=200):
    """Compute ICL fraction at multiple SB thresholds with bootstrap errors.

    Parameters
    ----------
    profile : list of dict
        SB profile with keys MU, FLUX.
    mu_min : float
        Minimum SB threshold (mag/arcsec^2).
    mu_max : float
        Maximum SB threshold (mag/arcsec^2).
    mu_step : float
        Step between thresholds.
    n_bootstrap : int
        Number of bootstrap resamples for error estimation.

    Returns
    -------
    results : list of dict
        Each dict has MU_THRESHOLD, F_ICL, F_ICL_ERR, L_ICL, L_GAL.
    """
    thresholds = np.arange(mu_min, mu_max + mu_step / 2.0, mu_step)

    # Filter valid profile entries
    valid = [p for p in profile if p['MU'] < 99.0 and p['FLUX'] > 0]
    if not valid:
        return [{'MU_THRESHOLD': float(t), 'F_ICL': 0.0, 'F_ICL_ERR': 0.0,
                 'L_ICL': 0.0, 'L_GAL': 0.0} for t in thresholds]

    results = []
    for mu_thresh in thresholds:
        f_icl, l_icl, l_gal = compute_icl_fraction(profile, mu_thresh)

        # Bootstrap error
        if n_bootstrap > 0 and len(valid) > 2:
            rng = np.random.default_rng(42)
            f_boot = []
            for _ in range(n_bootstrap):
                idx = rng.choice(len(valid), size=len(valid), replace=True)
                sample = [valid[i] for i in idx]
                fb, _, _ = compute_icl_fraction(sample, mu_thresh)
                f_boot.append(fb)
            f_icl_err = float(np.std(f_boot))
        else:
            f_icl_err = 0.0

        results.append({
            'MU_THRESHOLD': float(mu_thresh),
            'F_ICL': float(f_icl),
            'F_ICL_ERR': float(f_icl_err),
            'L_ICL': float(l_icl),
            'L_GAL': float(l_gal),
        })

    return results


def compute_isophotal_radii(profile, mu_levels):
    """Compute radii at specified surface brightness levels.

    Interpolates the profile to find R where mu(R) = level.

    Parameters
    ----------
    profile : list of dict
        SB profile from measure_sb_profile.
    mu_levels : list of float
        Surface brightness levels (mag/arcsec^2).

    Returns
    -------
    radii : dict
        {mu_level: R_pix} for each level. NaN if level not reached.
    """
    # Extract valid profile points
    r_arr = []
    mu_arr = []
    for p in profile:
        if p['MU'] < 99.0 and p['FLUX'] > 0:
            r_arr.append(p['R_PIX'])
            mu_arr.append(p['MU'])

    if len(r_arr) < 2:
        return {level: float('nan') for level in mu_levels}

    r_arr = np.array(r_arr)
    mu_arr = np.array(mu_arr)

    radii = {}
    for level in mu_levels:
        # Find where mu crosses the level (mu increases outward)
        # Look for consecutive points bracketing the level
        found = False
        for i in range(len(mu_arr) - 1):
            if (mu_arr[i] <= level <= mu_arr[i + 1]) or \
               (mu_arr[i] >= level >= mu_arr[i + 1]):
                # Linear interpolation
                dmu = mu_arr[i + 1] - mu_arr[i]
                if abs(dmu) > 1e-10:
                    frac = (level - mu_arr[i]) / dmu
                    r_interp = r_arr[i] + frac * (r_arr[i + 1] - r_arr[i])
                    radii[level] = float(r_interp)
                else:
                    radii[level] = float(r_arr[i])
                found = True
                break

        if not found:
            radii[level] = float('nan')

    return radii


def summarize_icl(profile, config):
    """Compute comprehensive ICL measurements.

    Parameters
    ----------
    profile : list of dict
        SB profile.
    config : ICLConfig
        Configuration with icl_mu_threshold, icl_mu_levels.

    Returns
    -------
    summary : dict
        Contains F_ICL, L_ICL, L_GAL, L_TOTAL, and R_MU* for
        each isophotal level.
    """
    mu_threshold = config.icl_mu_threshold
    mu_levels = [float(x) for x in config.icl_mu_levels.split(',')]

    f_icl, l_icl, l_gal = compute_icl_fraction(profile, mu_threshold)
    iso_radii = compute_isophotal_radii(profile, mu_levels)

    summary = {
        'F_ICL': f_icl,
        'L_ICL': l_icl,
        'L_GAL': l_gal,
        'L_TOTAL': l_icl + l_gal,
        'MU_THRESHOLD': mu_threshold,
    }

    for level, r in iso_radii.items():
        key = f'R_MU{level:.0f}'
        summary[key] = r
        summary[f'{key}_ARCSEC'] = r * config.pixel_scale if not np.isnan(r) else float('nan')

    return summary


def pixel_sb_map(image, zeropoint, pixel_scale, mu_offset=0.0):
    """Per-pixel surface brightness (mag/arcsec^2); NaN where flux <= 0.

    ``mu_offset`` is subtracted, e.g. 10*log10(1+z) to correct cosmological
    dimming (+ a K-correction) so thresholds can be applied in the rest
    frame: mu_rest = mu_obs - mu_offset.
    """
    img = np.asarray(image, dtype=np.float64)
    with np.errstate(invalid='ignore', divide='ignore'):
        mu = (-2.5 * np.log10(img) + zeropoint
              + 2.5 * np.log10(pixel_scale ** 2) - mu_offset)
    mu[~(img > 0)] = np.nan
    return mu


def icl_fraction_pixels(image, valid, zeropoint, pixel_scale, mu_bright,
                        mu_faint=None, mu_offset=0.0, smooth_sigma=0.0,
                        n_bootstrap=200, block=64, seed=42,
                        include_nonpositive=False):
    """ICL fraction from a per-pixel surface-brightness threshold.

    The isophotal definition used by Rudick et al. (2011), Burke et al.
    (2015) and Montes & Trujillo (2018, Sec. 3.4)::

        f_ICL = sum(F_pix : mu_bright < mu_pix [<= mu_faint]) / sum(F_pix)

    over the ``valid`` pixels (unmasked = cluster galaxies + BCG + ICL).
    The pixel SB is classified on an optionally Gaussian-smoothed copy of
    the image (``smooth_sigma`` px, NaN-aware) to reduce noise scatter
    across the threshold; the summed fluxes are always the unsmoothed ones.
    Pixels with flux <= 0 in the classification image have no defined
    surface brightness.  By default they are not ICL pixels in either mode
    ("all the pixels fainter than the threshold" = pixels with a defined
    mu > mu_bright), so the open-ended set always contains the slice and
    f(mu > mu_bright) >= f(mu_bright < mu <= mu_faint).  Counting them
    (``include_nonpositive=True``, open-ended mode only) adds the negative
    sky-noise tail of the whole faint area, which can make the open-ended
    fraction smaller than the slice.  The denominator is the total flux of
    all valid pixels (positive and negative, so that sky noise cancels).

    Uncertainty: block bootstrap over ``block`` x ``block`` px tiles.

    Returns
    -------
    dict with F_ICL, F_ICL_ERR, L_ICL, L_TOTAL, N_PIX_ICL, N_PIX.
    """
    img = np.asarray(image, dtype=np.float64)
    ok = np.asarray(valid, dtype=bool) & np.isfinite(img)
    cls = img
    if smooth_sigma and smooth_sigma > 0:
        from scipy.ndimage import gaussian_filter
        z = np.where(ok, img, 0.0)
        w = gaussian_filter(ok.astype(np.float64), smooth_sigma)
        with np.errstate(invalid='ignore', divide='ignore'):
            cls = gaussian_filter(z, smooth_sigma) / w
    mu = pixel_sb_map(cls, zeropoint, pixel_scale, mu_offset)
    if mu_faint is not None:
        sel = (mu > mu_bright) & (mu <= mu_faint)
    elif include_nonpositive:
        sel = ~(mu <= mu_bright)               # includes NaN (flux <= 0)
    else:
        sel = mu > mu_bright                   # defined (positive-flux) pixels only
    sel &= ok

    f_all = np.where(ok, img, 0.0)
    f_icl = np.where(sel, img, 0.0)
    l_tot = float(f_all.sum())
    l_icl = float(f_icl.sum())
    frac = l_icl / l_tot if l_tot > 0 else float('nan')

    err = float('nan')
    if n_bootstrap and n_bootstrap > 0:
        ny, nx = img.shape
        by, bx = int(np.ceil(ny / block)), int(np.ceil(nx / block))
        pad = ((0, by * block - ny), (0, bx * block - nx))
        t_all = np.pad(f_all, pad).reshape(by, block, bx, block).sum((1, 3)).ravel()
        t_icl = np.pad(f_icl, pad).reshape(by, block, bx, block).sum((1, 3)).ravel()
        keep = np.pad(ok, pad).reshape(by, block, bx, block).any((1, 3)).ravel()
        t_all, t_icl = t_all[keep], t_icl[keep]
        rng = np.random.default_rng(seed)
        boots = []
        for _ in range(n_bootstrap):
            idx = rng.integers(0, len(t_all), len(t_all))
            den = t_all[idx].sum()
            if den > 0:
                boots.append(t_icl[idx].sum() / den)
        if len(boots) > 1:
            err = float(np.std(boots))
    return {'F_ICL': frac, 'F_ICL_ERR': err, 'L_ICL': l_icl,
            'L_TOTAL': l_tot, 'N_PIX_ICL': int(sel.sum()),
            'N_PIX': int(ok.sum()), 'MU_BRIGHT': float(mu_bright),
            'MU_FAINT': None if mu_faint is None else float(mu_faint)}
