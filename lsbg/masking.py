"""Magnitude-based selective masking for LSBG detection.

Masks only bright sources (mag < threshold) so that faint LSBG
candidates are preserved while bright PSF wings/halos are removed.
Reuses ICL masking functions for core operations.
"""

import numpy as np
try:
    import sep
except ImportError:
    import sep_pjw as sep
from scipy.ndimage import binary_dilation

from icl.masking import create_source_mask, mask_bright_stars, interpolate_masked


def mask_by_magnitude(data, segmap, catalog, mag_threshold=22.0,
                      expand_factor=3.0):
    """Selectively mask only sources brighter than mag_threshold.

    Parameters
    ----------
    data : 2D array
        Image data (used for shape).
    segmap : 2D int array
        Segmentation map from SEP (0 = background).
    catalog : structured array
        SEP extraction catalog with 'flux' field.
    mag_threshold : float
        Magnitude limit; sources brighter (mag < threshold) are masked.
    expand_factor : float
        Dilation factor relative to segment equivalent radius.

    Returns
    -------
    mask : 2D bool array
        True where bright sources are masked.
    bright_labels : list
        Segment labels that were masked.
    """
    mask = np.zeros(data.shape, dtype=bool)
    bright_labels = []

    # Compute rough magnitudes from flux
    flux = catalog['flux']
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = np.where(flux > 0, -2.5 * np.log10(flux), 99.0)

    labels = np.unique(segmap)
    labels = labels[labels > 0]

    for i, lab in enumerate(labels):
        if i >= len(mag):
            break
        if mag[i] >= mag_threshold:
            continue

        seg_pix = segmap == lab
        npix = seg_pix.sum()
        if npix == 0:
            continue

        bright_labels.append(lab)

        r_eq = np.sqrt(npix / np.pi)
        r_dilate = max(1, int(r_eq * expand_factor))

        y, x = np.ogrid[-r_dilate:r_dilate + 1, -r_dilate:r_dilate + 1]
        struct = (x * x + y * y) <= r_dilate * r_dilate

        dilated = binary_dilation(seg_pix, structure=struct)
        mask |= dilated

    return mask, bright_labels


def iterative_mask_refine(data, mask, bkg_estimate, sigma_thresh=2.0,
                          minarea=10):
    """Re-detect residual sources and expand mask.

    After initial masking and background subtraction, residual
    sources may remain. This detects them and adds to the mask.

    Parameters
    ----------
    data : 2D array
        Original image data.
    mask : 2D bool array
        Current mask (modified in-place and returned).
    bkg_estimate : 2D array
        Current background estimate.
    sigma_thresh : float
        Detection threshold in sigma for residuals.
    minarea : int
        Minimum area for residual detection.

    Returns
    -------
    mask : 2D bool array
        Updated mask with additional detections.
    n_new : int
        Number of new sources detected.
    """
    residual = data - bkg_estimate
    residual_c = np.ascontiguousarray(residual, dtype=np.float64)

    # Mask current sources in residual
    residual_c[mask] = 0.0

    try:
        bkg = sep.Background(residual_c)
        rms = bkg.rms()
        objects, new_segmap = sep.extract(
            residual_c - bkg.back(), thresh=sigma_thresh,
            err=rms, minarea=minarea, segmentation_map=True
        )
    except Exception:
        return mask, 0

    n_new = len(objects)
    if n_new > 0:
        new_mask = new_segmap > 0
        # Dilate new detections slightly
        struct = np.ones((5, 5), dtype=bool)
        new_mask = binary_dilation(new_mask, structure=struct)
        mask = mask | new_mask

    return mask, n_new


def create_lsbg_mask(data, config):
    """Full LSBG masking pipeline.

    1. Initial SEP detection
    2. Magnitude-based selective masking
    3. Bright star masking
    4. Interpolation of masked regions

    Parameters
    ----------
    data : 2D array
        Input image.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    mask : 2D bool array
        Final mask.
    masked_data : 2D array
        Image with masked regions interpolated.
    """
    import sys

    data_c = np.ascontiguousarray(data, dtype=np.float64)

    # Initial SEP detection for segmentation
    bkg = sep.Background(data_c)
    data_sub = data_c - bkg.back()

    objects, segmap = sep.extract(
        data_sub, thresh=config.mask_detect_thresh,
        err=bkg.rms(), minarea=config.mask_detect_minarea,
        segmentation_map=True
    )
    print(f"LSBG mask: detected {len(objects)} initial sources",
          file=sys.stderr)

    # Magnitude-based selective masking
    mask, bright_labels = mask_by_magnitude(
        data, segmap, objects,
        mag_threshold=config.mask_mag_threshold,
        expand_factor=config.mask_expand_factor
    )
    print(f"LSBG mask: masked {len(bright_labels)} bright sources "
          f"(mag < {config.mask_mag_threshold})", file=sys.stderr)

    # Bright star masking
    if len(objects) > 0:
        catalog = {
            'x_image': objects['x'] + 1.0,
            'y_image': objects['y'] + 1.0,
            'mag_auto': -2.5 * np.log10(np.maximum(objects['flux'], 1e-30))
        }
        mask_bright_stars(mask, catalog,
                          mag_limit=config.bright_star_mag_limit,
                          radius_scale=config.bright_star_radius_scale)

    n_masked = mask.sum()
    print(f"LSBG mask: {n_masked} pixels masked "
          f"({100.0 * n_masked / mask.size:.1f}%)", file=sys.stderr)

    # Interpolate masked regions
    masked_data = interpolate_masked(data, mask, method=config.interp_method)

    return mask, masked_data
