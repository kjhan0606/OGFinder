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
