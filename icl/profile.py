"""Surface brightness profile measurement for ICL analysis.

Measures azimuthally-averaged (or sector) surface brightness
profiles in concentric annuli around the BCG center.
"""

import numpy as np


def find_bcg(data, catalog):
    """Find the brightest cluster galaxy (BCG) in the catalog.

    Selects the brightest extended source: highest flux among
    sources with large A_IMAGE and low CLASS_STAR.

    Parameters
    ----------
    data : 2D array
        Image data (used only for shape/bounds check).
    catalog : dict
        Must contain 'x_image', 'y_image', 'flux_auto'.
        Optionally 'a_image', 'class_star'.

    Returns
    -------
    bcg_x, bcg_y : float
        BCG center in 0-indexed pixel coordinates.
    """
    flux = catalog.get('flux_auto')
    if flux is None:
        raise ValueError("Catalog must contain flux_auto")

    x_img = catalog['x_image']
    y_img = catalog['y_image']

    # Prefer extended sources
    class_star = catalog.get('class_star', np.zeros(len(flux)))
    a_image = catalog.get('a_image', np.ones(len(flux)))

    # Score: high flux * large size * low stellarity
    score = flux * a_image * (1.0 - class_star)
    best = np.argmax(score)

    # Convert to 0-based
    bcg_x = x_img[best] - 1.0
    bcg_y = y_img[best] - 1.0

    return bcg_x, bcg_y


def measure_sb_profile(data, cx, cy, config):
    """Measure surface brightness profile in concentric annuli.

    Uses pixel-by-pixel computation with NaN exclusion for robust
    measurement at large radii where edge/masked pixels are common.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    cx, cy : float
        Center coordinates (0-indexed).
    config : ICLConfig
        Configuration with profile_rmin, profile_rmax, profile_nsteps,
        profile_spacing, profile_ellipticity, profile_pa,
        mag_zeropoint, pixel_scale.

    Returns
    -------
    profile : list of dict
        Each dict has: R_PIX, R_ARCSEC, MU, MU_ERR, FLUX, AREA, NPIX
    """
    rmin = config.profile_rmin
    rmax = config.profile_rmax
    nsteps = config.profile_nsteps
    zp = config.mag_zeropoint
    pscale = config.pixel_scale
    ell = config.profile_ellipticity
    pa = config.profile_pa

    # Generate radii
    if config.profile_spacing == 'log':
        radii = np.logspace(np.log10(max(rmin, 0.5)), np.log10(rmax), nsteps)
    else:
        radii = np.linspace(rmin, rmax, nsteps)

    ny, nx = data.shape

    # Build distance map (elliptical if needed)
    yy, xx = np.mgrid[0:ny, 0:nx]
    dx = (xx - cx).astype(np.float64)
    dy = (yy - cy).astype(np.float64)

    if ell > 0:
        theta = np.radians(pa)
        cos_t, sin_t = np.cos(theta), np.sin(theta)
        dx_rot = dx * cos_t + dy * sin_t
        dy_rot = -dx * sin_t + dy * cos_t
        r_map = np.sqrt(dx_rot ** 2 + (dy_rot / (1.0 - ell)) ** 2)
    else:
        r_map = np.sqrt(dx * dx + dy * dy)

    # Valid pixel mask (exclude NaN/inf)
    valid = np.isfinite(data)

    profile = []
    for i in range(len(radii)):
        r_out = radii[i]
        r_in = radii[i - 1] if i > 0 else 0.0

        annulus = (r_map >= r_in) & (r_map < r_out) & valid
        npix = int(annulus.sum())

        if npix < 1:
            continue

        vals = data[annulus]
        flux = float(np.sum(vals))
        area = float(npix)

        # Use median for robust SB (less sensitive to bright source residuals)
        median_val = float(np.median(vals))
        flux_per_pix = median_val

        if flux_per_pix > 0:
            mu = -2.5 * np.log10(flux_per_pix) + zp + 2.5 * np.log10(pscale ** 2)
            # Error from scatter in pixel values
            std_val = float(np.std(vals))
            if std_val > 0 and npix > 1:
                se = std_val / np.sqrt(npix)
                # Propagate to magnitude
                mu_err = 2.5 / np.log(10) * se / flux_per_pix
            else:
                mu_err = 99.0
        else:
            mu = 99.0
            mu_err = 99.0

        r_mid = (r_in + r_out) / 2.0 if r_in > 0 else r_out / 2.0

        profile.append({
            'R_PIX': r_mid,
            'R_ARCSEC': r_mid * pscale,
            'MU': mu,
            'MU_ERR': mu_err,
            'FLUX': flux,
            'AREA': area,
            'NPIX': npix,
        })

    return profile


def measure_sector_profile(data, cx, cy, sector_pa, sector_width, config):
    """Measure SB profile in a specific angular sector.

    Useful for analyzing tidal streams or asymmetric ICL features.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    cx, cy : float
        Center (0-indexed).
    sector_pa : float
        Position angle of sector center (degrees, E of N).
    sector_width : float
        Angular width of sector (degrees).
    config : ICLConfig
        Profile configuration.

    Returns
    -------
    profile : list of dict
        Same format as measure_sb_profile.
    """
    from .config import ICLConfig
    import copy

    # Override sector params in config
    cfg = copy.copy(config)
    cfg.profile_sector_pa = sector_pa
    cfg.profile_sector_width = sector_width

    ny, nx = data.shape
    rmin = cfg.profile_rmin
    rmax = cfg.profile_rmax
    nsteps = cfg.profile_nsteps
    zp = cfg.mag_zeropoint
    pscale = cfg.pixel_scale

    # Generate radii
    if cfg.profile_spacing == 'log':
        radii = np.logspace(np.log10(max(rmin, 0.5)), np.log10(rmax), nsteps)
    else:
        radii = np.linspace(rmin, rmax, nsteps)

    # Build pixel coordinate arrays
    yy, xx = np.mgrid[0:ny, 0:nx]
    dx = xx - cx
    dy = yy - cy
    r_map = np.sqrt(dx * dx + dy * dy)
    angle_map = np.degrees(np.arctan2(dy, dx))

    # Sector angular mask
    angle_diff = (angle_map - sector_pa + 180) % 360 - 180
    in_sector = np.abs(angle_diff) <= sector_width / 2.0

    profile = []
    for i in range(len(radii)):
        r_out = radii[i]
        r_in = radii[i - 1] if i > 0 else 0.0

        valid = np.isfinite(data)
        annulus = (r_map >= r_in) & (r_map < r_out) & in_sector & valid
        npix = annulus.sum()

        if npix < 1:
            continue

        vals = data[annulus]
        flux = float(np.sum(vals))
        area = float(npix)
        flux_per_pix = float(np.median(vals))

        if flux_per_pix > 0:
            mu = -2.5 * np.log10(flux_per_pix) + zp + 2.5 * np.log10(pscale ** 2)
            std_val = float(np.std(vals))
            if std_val > 0 and npix > 1:
                se = std_val / np.sqrt(npix)
                mu_err = 2.5 / np.log(10) * se / flux_per_pix
            else:
                mu_err = 99.0
        else:
            mu = 99.0
            mu_err = 99.0

        r_mid = (r_in + r_out) / 2.0 if r_in > 0 else r_out / 2.0

        profile.append({
            'R_PIX': r_mid,
            'R_ARCSEC': r_mid * pscale,
            'MU': mu,
            'MU_ERR': mu_err,
            'FLUX': flux,
            'AREA': area,
            'NPIX': int(npix),
        })

    return profile
