"""Sersic galaxy stamps.

The default stamp is drawn in-tree: a Sérsic surface brightness, reduced
shear, pixel integration, and a Gaussian PSF. GalSim runs only when
``backend='galsim'`` and the user-installed package imports.
"""

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.special import gammaincinv, gammaln

# Samples along each pixel edge. The Sérsic is integrated, then the PSF is applied.
_OVERSAMPLE = 5
_FWHM_TO_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


def make_galaxy_stamp(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                      psf_fwhm, x_offset=0.0, y_offset=0.0, *,
                      backend="numpy"):
    """
    Generate a single Sersic galaxy stamp.

    Parameters
    ----------
    stamp_size : int
        Size of the output image (stamp_size x stamp_size).
    flux : float
        Total flux of the galaxy.
    re : float
        Effective (half-light) radius in pixels.
    sersic_n : float
        Sersic index.
    ellip : float
        Ellipticity (1 - b/a), in [0, 1).
    theta_rad : float
        Position angle in radians.
    psf_fwhm : float
        PSF FWHM in pixels.
    x_offset, y_offset : float
        Sub-pixel offset from stamp center.

    Returns
    -------
    image : np.ndarray
        2D stamp image (float64).
    backend : {'numpy', 'galsim'}
        ``numpy`` is the default. ``galsim`` uses the installed package.
    """
    if backend == "galsim":
        return _make_galsim(stamp_size, flux, re, sersic_n, ellip,
                            theta_rad, psf_fwhm, x_offset, y_offset)
    if backend != "numpy":
        raise ValueError("backend must be 'numpy' or 'galsim'")
    return _make_numpy(stamp_size, flux, re, sersic_n, ellip,
                       theta_rad, psf_fwhm, x_offset, y_offset)


def _make_galsim(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                 psf_fwhm, x_offset, y_offset):
    """GalSim backend. Called only for an explicit backend='galsim'."""
    try:
        import galsim
    except Exception as exc:
        raise RuntimeError(
            "galsim was not found. The in-tree stamp is numpy."
        ) from exc
    # Clamp sersic_n to GalSim's valid range
    n = max(0.3, min(sersic_n, 6.2))
    re_safe = max(0.5, re)

    gal = galsim.Sersic(n=n, half_light_radius=re_safe, flux=flux)

    # Apply shear for ellipticity
    if ellip > 0.01:
        e = min(ellip, 0.95)
        # Convert ellipticity + PA to shear components
        g1 = e * np.cos(2 * theta_rad)
        g2 = e * np.sin(2 * theta_rad)
        gal = gal.shear(g1=g1, g2=g2)

    # PSF convolution
    psf = galsim.Gaussian(fwhm=max(0.5, psf_fwhm))
    final = galsim.Convolve([gal, psf])

    # Offset
    final = final.shift(x_offset, y_offset)

    # Draw image
    img = galsim.Image(stamp_size, stamp_size, scale=1.0)
    final.drawImage(image=img, method='auto')

    return img.array.astype(np.float64)


def sersic_b(n):
    """Half-light parameter b_n from the regularised lower incomplete gamma."""
    return float(gammaincinv(2.0 * n, 0.5))


def sersic_surface_brightness(r, re, n, flux):
    """Round Sérsic surface brightness. The integral over the plane equals ``flux``."""
    r = np.asarray(r, dtype=np.float64)
    re = float(re)
    n = float(n)
    flux = float(flux)
    b = sersic_b(n)
    log_ie = (np.log(flux) + 2.0 * n * np.log(b) - np.log(2.0 * np.pi)
              - np.log(n) - b - gammaln(2.0 * n) - 2.0 * np.log(re))
    expo = -b * (np.power(r / re, 1.0 / n) - 1.0)
    return np.exp(log_ie + expo)


def sheared_radius(x, y, ellip, theta_rad):
    """Intrinsic radius after an area-preserving reduced shear.

    ``ellip`` is the shear magnitude |g|, capped at 0.95. ``theta_rad`` is the
    position angle of that shear, the same angle the explicit GalSim path
    passes to ``shear(g1, g2)``.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    g = min(max(float(ellip), 0.0), 0.95)
    if g < 1e-12:
        return np.hypot(x, y)
    g1 = g * np.cos(2.0 * theta_rad)
    g2 = g * np.sin(2.0 * theta_rad)
    norm = np.sqrt(1.0 - g * g)
    xi = ((1.0 - g1) * x - g2 * y) / norm
    yi = (-g2 * x + (1.0 + g1) * y) / norm
    return np.hypot(xi, yi)


def _draw_sersic_fine(stamp_size, flux, re, n, ellip, theta_rad, x_offset, y_offset):
    """Flux in each subpixel of an oversampled stamp. PSF is not applied yet."""
    step = 1.0 / _OVERSAMPLE
    n_fine = int(stamp_size) * _OVERSAMPLE
    pos = (np.arange(n_fine, dtype=np.float64) + 0.5) * step
    half = stamp_size / 2.0
    x_rel = pos - half - float(x_offset)
    y_rel = pos - half - float(y_offset)
    yy, xx = np.meshgrid(y_rel, x_rel, indexing="ij")
    radius = sheared_radius(xx, yy, ellip, theta_rad)
    return sersic_surface_brightness(radius, re, n, flux) * (step * step)


def convolve_gaussian_fwhm(image, fwhm_pixels, pixel_scale=1.0):
    """Convolve with a unit-sum Gaussian. ``fwhm_pixels`` is the full width at half maximum."""
    sigma = max(float(fwhm_pixels), 0.5) / _FWHM_TO_SIGMA / float(pixel_scale)
    return gaussian_filter(np.asarray(image, dtype=np.float64), sigma=sigma, mode="constant", cval=0.0)


def _make_numpy(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                psf_fwhm, x_offset, y_offset):
    """In-tree Sérsic stamp: shear, pixel integral, then a Gaussian PSF."""
    stamp_size = int(stamp_size)
    if stamp_size < 1 or float(flux) <= 0.0:
        return np.zeros((max(stamp_size, 0), max(stamp_size, 0)), dtype=np.float64)
    n = max(0.3, min(float(sersic_n), 6.2))
    re_safe = max(0.5, float(re))
    fine = _draw_sersic_fine(
        stamp_size, float(flux), re_safe, n, ellip, theta_rad, x_offset, y_offset)
    blurred = convolve_gaussian_fwhm(fine, psf_fwhm, pixel_scale=1.0 / _OVERSAMPLE)
    blocks = blurred.reshape(stamp_size, _OVERSAMPLE, stamp_size, _OVERSAMPLE)
    return blocks.sum(axis=(1, 3))


def make_pair_stamp(stamp_size, params1, params2, psf_fwhm, noise_std):
    """
    Generate a stamp containing two galaxy components.

    Parameters
    ----------
    stamp_size : int
    params1, params2 : dict
        Each with keys: flux, re, sersic_n, ellip, theta_rad, x_offset, y_offset
    psf_fwhm : float
    noise_std : float

    Returns
    -------
    image : np.ndarray
        Noisy stamp with both components.
    clean : np.ndarray
        Noiseless stamp.
    """
    img1 = make_galaxy_stamp(stamp_size, psf_fwhm=psf_fwhm, **params1)
    img2 = make_galaxy_stamp(stamp_size, psf_fwhm=psf_fwhm, **params2)
    clean = img1 + img2
    noise = np.random.normal(0, noise_std, clean.shape) if noise_std > 0 else 0.0
    return (clean + noise), clean
