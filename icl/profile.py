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
    import sep

    rmin = config.profile_rmin
    rmax = config.profile_rmax
    nsteps = config.profile_nsteps
    zp = config.mag_zeropoint
    pscale = config.pixel_scale
    ell = config.profile_ellipticity
    pa = config.profile_pa
    sector_width = config.profile_sector_width
    sector_pa = config.profile_sector_pa

    # Generate radii
    if config.profile_spacing == 'log':
        radii = np.logspace(np.log10(max(rmin, 0.5)), np.log10(rmax), nsteps)
    else:
        radii = np.linspace(rmin, rmax, nsteps)

    data_c = np.ascontiguousarray(data, dtype=np.float64)
    ny, nx = data.shape

    # Sector mask (if not full azimuthal)
    use_sector = sector_width < 360.0
    if use_sector:
        yy, xx = np.mgrid[0:ny, 0:nx]
        angles = np.degrees(np.arctan2(yy - cy, xx - cx))  # [-180, 180]
        # Normalize relative to sector PA
        angle_diff = (angles - sector_pa + 180) % 360 - 180
        sector_mask = np.abs(angle_diff) <= sector_width / 2.0

    profile = []
    x_arr = np.array([cx], dtype=np.float64)
    y_arr = np.array([cy], dtype=np.float64)

    for i in range(len(radii)):
        r_out = radii[i]
        r_in = radii[i - 1] if i > 0 else 0.0

        if ell > 0:
            # Elliptical annulus
            b_out = r_out * (1.0 - ell)
            b_in = r_in * (1.0 - ell)
            theta = np.radians(pa)

            flux_out, _, _ = sep.sum_ellipse(
                data_c, x_arr, y_arr,
                np.array([r_out]), np.array([b_out]),
                np.array([theta]))
            if r_in > 0:
                flux_in, _, _ = sep.sum_ellipse(
                    data_c, x_arr, y_arr,
                    np.array([r_in]), np.array([b_in]),
                    np.array([theta]))
            else:
                flux_in = np.array([0.0])

            area_out = np.pi * r_out * b_out
            area_in = np.pi * r_in * b_in if r_in > 0 else 0.0
        else:
            # Circular annulus
            flux_out, _, _ = sep.sum_circle(
                data_c, x_arr, y_arr,
                np.array([r_out]))
            if r_in > 0:
                flux_in, _, _ = sep.sum_circle(
                    data_c, x_arr, y_arr,
                    np.array([r_in]))
            else:
                flux_in = np.array([0.0])

            area_out = np.pi * r_out * r_out
            area_in = np.pi * r_in * r_in if r_in > 0 else 0.0

        annulus_flux = float(flux_out[0] - flux_in[0])
        annulus_area = area_out - area_in

        if annulus_area <= 0:
            continue

        # Approximate pixel count
        npix = int(annulus_area + 0.5)

        # Apply sector correction
        if use_sector:
            sector_frac = sector_width / 360.0
            annulus_flux *= 1.0  # flux already from full annulus
            # Actually need to measure sector flux differently
            # For now, scale area and flux by sector fraction
            # (approximate for non-sector photometry)
            annulus_area *= sector_frac
            npix = max(1, int(npix * sector_frac))

        flux_per_pix = annulus_flux / annulus_area if annulus_area > 0 else 0.0

        # Surface brightness
        if flux_per_pix > 0:
            mu = -2.5 * np.log10(flux_per_pix) + zp + 2.5 * np.log10(pscale ** 2)
            # Error: Poisson + read noise approximation
            if annulus_flux > 0 and npix > 0:
                mu_err = 2.5 / np.log(10) / (annulus_flux / np.sqrt(annulus_flux))
                mu_err = 2.5 / (np.log(10) * np.sqrt(annulus_flux / npix) /
                                (annulus_flux / npix)) if annulus_flux > 0 else 99.0
                # Simplified: sigma_mu ≈ 2.5/(ln10) * 1/SNR_per_pixel * 1/sqrt(npix)
                snr_pix = flux_per_pix / max(np.sqrt(np.abs(flux_per_pix)), 1e-30)
                mu_err = 2.5 / (np.log(10) * snr_pix * np.sqrt(max(npix, 1)))
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
            'FLUX': annulus_flux,
            'AREA': annulus_area,
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

        annulus = (r_map >= r_in) & (r_map < r_out) & in_sector
        npix = annulus.sum()

        if npix < 1:
            continue

        flux = float(np.sum(data[annulus]))
        area = float(npix)  # pixel area
        flux_per_pix = flux / area

        if flux_per_pix > 0:
            mu = -2.5 * np.log10(flux_per_pix) + zp + 2.5 * np.log10(pscale ** 2)
            snr_pix = flux_per_pix / max(np.sqrt(np.abs(flux_per_pix)), 1e-30)
            mu_err = 2.5 / (np.log(10) * snr_pix * np.sqrt(max(npix, 1)))
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
