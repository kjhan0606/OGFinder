"""Sersic galaxy profile generation with GalSim backend and numpy fallback."""

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.special import gammaincinv, gamma

try:
    import galsim
    HAS_GALSIM = True
except ImportError:
    HAS_GALSIM = False


def make_galaxy_stamp(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                      psf_fwhm, x_offset=0.0, y_offset=0.0):
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
    """
    if HAS_GALSIM:
        return _make_galsim(stamp_size, flux, re, sersic_n, ellip,
                            theta_rad, psf_fwhm, x_offset, y_offset)
    else:
        return _make_numpy(stamp_size, flux, re, sersic_n, ellip,
                           theta_rad, psf_fwhm, x_offset, y_offset)


def _make_galsim(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                 psf_fwhm, x_offset, y_offset):
    """GalSim backend."""
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


def _make_numpy(stamp_size, flux, re, sersic_n, ellip, theta_rad,
                psf_fwhm, x_offset, y_offset):
    """Numpy fallback: analytic Sersic + Gaussian PSF convolution."""
    n = max(0.3, min(sersic_n, 6.2))
    re_safe = max(0.5, re)

    # Sersic b_n parameter
    b_n = gammaincinv(2.0 * n, 0.5)

    # Coordinate grid
    cy, cx = stamp_size / 2.0 + y_offset, stamp_size / 2.0 + x_offset
    y, x = np.mgrid[0:stamp_size, 0:stamp_size].astype(np.float64)
    x = x - cx + 0.5
    y = y - cy + 0.5

    # Rotation
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)
    xr = x * cos_t + y * sin_t
    yr = -x * sin_t + y * cos_t

    # Elliptical radius
    q = max(1.0 - ellip, 0.05)
    r = np.sqrt(xr**2 + (yr / q)**2)

    # Sersic profile (unnormalized)
    profile = np.exp(-b_n * ((r / re_safe)**(1.0 / n) - 1.0))

    # Normalize to total flux
    total = profile.sum()
    if total > 0:
        profile *= flux / total

    # PSF convolution using Gaussian filter
    psf_sigma = max(0.3, psf_fwhm) / 2.3548  # FWHM to sigma
    image = gaussian_filter(profile, sigma=psf_sigma)

    return image


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
