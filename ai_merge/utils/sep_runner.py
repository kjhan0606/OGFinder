"""SEP extraction wrapper mirroring ds9_sextract.py pipeline."""

import numpy as np
import sep


def _photometry(data_sub, objects, a_safe, b_safe, bkg_rms, mag_zeropoint, gain=None):
    """Common photometry pipeline for both stamp and FITS modes."""
    n = len(objects)

    kronrad, _ = sep.kron_radius(
        data_sub, objects['x'], objects['y'],
        a_safe, b_safe, objects['theta'], 6.0
    )
    kronrad = np.maximum(kronrad, 3.5)

    # Try elliptical photometry, fall back to circular
    try:
        flux_auto, _, _ = sep.sum_ellipse(
            data_sub, objects['x'], objects['y'],
            a_safe, b_safe, objects['theta'],
            2.5 * kronrad, err=bkg_rms,
            gain=gain, subpix=5
        )
    except Exception:
        # Fallback: circular aperture with equivalent radius
        r_circ = 2.5 * kronrad * np.sqrt(a_safe * b_safe)
        r_circ = np.clip(r_circ, 2.0, 500.0)
        flux_auto, _, _ = sep.sum_circle(
            data_sub, objects['x'], objects['y'],
            r_circ, err=bkg_rms, gain=gain, subpix=5
        )

    flux_radius, _ = sep.flux_radius(
        data_sub, objects['x'], objects['y'],
        6.0 * a_safe, 0.5, normflux=flux_auto
    )

    a, b = a_safe, b_safe
    with np.errstate(divide='ignore', invalid='ignore'):
        ellipticity = np.where(a > 0, 1.0 - b / a, 0.0)
        ellipticity = np.clip(ellipticity, 0.0, 1.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        mag_auto = np.where(
            flux_auto > 0,
            -2.5 * np.log10(flux_auto) + mag_zeropoint,
            99.0
        )

    fwhm_image = 2.0 * np.sqrt(np.log(2.0) * (a**2 + b**2))
    with np.errstate(divide='ignore', invalid='ignore'):
        fwhm_ratio = fwhm_image / 3.0
        class_star = 1.0 / (1.0 + np.exp(3.0 * (fwhm_ratio - 1.5)))
        class_star = np.clip(class_star, 0.0, 1.0)

    return kronrad, flux_auto, flux_radius, ellipticity, mag_auto, class_star


def run_sep_on_stamp(image, detect_thresh=1.5, detect_minarea=5,
                     deblend_nthresh=32, deblend_mincont=0.005,
                     mag_zeropoint=25.0):
    """
    Run SEP source extraction on a 2D stamp image.
    Returns dict or None if no sources detected.
    """
    data = np.ascontiguousarray(image, dtype=np.float64)

    try:
        bkg = sep.Background(data, bw=32, bh=32, fw=3, fh=3)
    except Exception:
        bkg = sep.Background(data)
    data_sub = data - bkg.back()
    bkg_rms = bkg.rms()

    try:
        objects, segmap = sep.extract(
            data_sub, detect_thresh, err=bkg_rms,
            minarea=detect_minarea,
            deblend_nthresh=deblend_nthresh,
            deblend_cont=deblend_mincont,
            segmentation_map=True
        )
    except Exception:
        return None

    n = len(objects)
    if n == 0:
        return None

    a_safe = np.nan_to_num(objects['a'], nan=0.5)
    a_safe = np.maximum(a_safe, 0.5)
    b_safe = np.nan_to_num(objects['b'], nan=0.5)
    b_safe = np.clip(b_safe, 0.5, a_safe)

    try:
        kronrad, flux_auto, flux_radius, ellipticity, mag_auto, class_star = \
            _photometry(data_sub, objects, a_safe, b_safe, bkg_rms, mag_zeropoint)
    except Exception:
        return None

    return {
        'objects': objects,
        'segmap': segmap,
        'n': n,
        'flux_auto': flux_auto,
        'kron_radius': kronrad,
        'flux_radius': flux_radius,
        'ellipticity': ellipticity,
        'mag_auto': mag_auto,
        'class_star': class_star,
        'data_sub': data_sub,
        'bkg_rms': bkg_rms,
    }


def run_sep_on_fits(data, back_size=64, back_filtersize=3, **kwargs):
    """
    Run SEP on full FITS image data (larger background mesh).
    Filters edge objects to prevent aperture overflow.
    """
    data = np.ascontiguousarray(data, dtype=np.float64)

    detect_thresh = kwargs.get('detect_thresh', 1.5)
    detect_minarea = kwargs.get('detect_minarea', 5)
    deblend_nthresh = kwargs.get('deblend_nthresh', 32)
    deblend_mincont = kwargs.get('deblend_mincont', 0.005)
    mag_zeropoint = kwargs.get('mag_zeropoint', 25.0)
    gain = kwargs.get('gain', None)

    bkg = sep.Background(data, bw=back_size, bh=back_size,
                         fw=back_filtersize, fh=back_filtersize)
    data_sub = data - bkg.back()
    bkg_rms = bkg.rms()

    objects, segmap = sep.extract(
        data_sub, detect_thresh, err=bkg_rms,
        minarea=detect_minarea,
        deblend_nthresh=deblend_nthresh,
        deblend_cont=deblend_mincont,
        segmentation_map=True
    )

    n = len(objects)
    if n == 0:
        return None

    a_safe = np.nan_to_num(objects['a'], nan=0.5)
    a_safe = np.maximum(a_safe, 0.5)
    b_safe = np.nan_to_num(objects['b'], nan=0.5)
    b_safe = np.clip(b_safe, 0.5, a_safe)

    # Filter objects whose aperture would extend beyond image
    h, w = data_sub.shape
    aperture_r = 2.5 * 3.5 * a_safe + 5.0  # conservative estimate
    keep = ((objects['x'] - aperture_r > 0) &
            (objects['x'] + aperture_r < w) &
            (objects['y'] - aperture_r > 0) &
            (objects['y'] + aperture_r < h))

    if keep.sum() == 0:
        return None

    objects = objects[keep]
    a_safe = a_safe[keep]
    b_safe = b_safe[keep]
    n = int(keep.sum())

    try:
        kronrad, flux_auto, flux_radius, ellipticity, mag_auto, class_star = \
            _photometry(data_sub, objects, a_safe, b_safe, bkg_rms, mag_zeropoint, gain)
    except Exception:
        return None

    return {
        'objects': objects,
        'segmap': segmap,  # full segmap (labels won't match filtered indices)
        'n': n,
        'flux_auto': flux_auto,
        'kron_radius': kronrad,
        'flux_radius': flux_radius,
        'ellipticity': ellipticity,
        'mag_auto': mag_auto,
        'class_star': class_star,
        'data_sub': data_sub,
        'bkg_rms': bkg_rms,
    }
