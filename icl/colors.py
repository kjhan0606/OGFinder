"""Multi-band ICL color profiles (Phase 2).

Measures ICL surface brightness profiles in multiple bands
and computes color index profiles (e.g., B-V, g-r).
"""

import numpy as np


def measure_multiband_icl_profiles(band_files, mask, cx, cy, config):
    """Measure ICL SB profiles in multiple bands.

    Parameters
    ----------
    band_files : dict
        {band_name: fits_path} mapping.
    mask : 2D bool array
        Source mask (shared across bands).
    cx, cy : float
        BCG center (0-indexed).
    config : ICLConfig
        Profile configuration.

    Returns
    -------
    profiles : dict
        {band_name: profile_list} for each band.
    """
    from .profile import measure_sb_profile
    from astropy.io import fits

    profiles = {}
    for band_name, fpath in band_files.items():
        with fits.open(fpath) as hdul:
            # Find image data in primary or first image extension
            if hdul[0].data is not None:
                data = hdul[0].data.astype(np.float64)
            else:
                data = None
                for ext in hdul[1:]:
                    if ext.data is not None and ext.data.ndim >= 2:
                        data = ext.data.astype(np.float64)
                        break
                if data is None:
                    raise ValueError(f"No image data in {fpath}")
            while data.ndim > 2:
                data = data[0]

        # Zero out masked regions (background already subtracted)
        profiles[band_name] = measure_sb_profile(data, cx, cy, config)

    return profiles


def compute_color_profiles(profiles, band1, band2):
    """Compute color index profile from two bands.

    color(r) = mu_band1(r) - mu_band2(r)

    Parameters
    ----------
    profiles : dict
        {band_name: profile_list} from measure_multiband_icl_profiles.
    band1, band2 : str
        Band names to compute color between.

    Returns
    -------
    color_profile : list of dict
        Each dict has: R_PIX, R_ARCSEC, COLOR, COLOR_ERR
    """
    p1 = profiles.get(band1)
    p2 = profiles.get(band2)

    if p1 is None or p2 is None:
        return []

    # Match by radius
    r1_map = {round(p['R_PIX'], 3): p for p in p1}
    r2_map = {round(p['R_PIX'], 3): p for p in p2}

    common_r = sorted(set(r1_map.keys()) & set(r2_map.keys()))

    color_profile = []
    for r in common_r:
        pp1 = r1_map[r]
        pp2 = r2_map[r]

        mu1 = pp1['MU']
        mu2 = pp2['MU']

        if mu1 >= 99.0 or mu2 >= 99.0:
            continue

        color = mu1 - mu2
        color_err = np.sqrt(pp1['MU_ERR'] ** 2 + pp2['MU_ERR'] ** 2)

        color_profile.append({
            'R_PIX': pp1['R_PIX'],
            'R_ARCSEC': pp1['R_ARCSEC'],
            'COLOR': color,
            'COLOR_ERR': color_err,
        })

    return color_profile
