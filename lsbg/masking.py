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
                      expand_factor=3.0, lsb_protect=False,
                      lsb_mu_threshold=24.0, pixel_scale=0.06,
                      mag_zeropoint=25.0):
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
    lsb_protect : bool
        If True, protect sources with low mean surface brightness.
    lsb_mu_threshold : float
        Surface brightness threshold (mag/arcsec^2); protect if fainter.
    pixel_scale : float
        Pixel scale in arcsec/pixel.
    mag_zeropoint : float
        Magnitude zeropoint.

    Returns
    -------
    mask : 2D bool array
        True where bright sources are masked.
    bright_labels : list
        Segment labels that were masked.
    """
    import sys
    from collections import defaultdict

    mask = np.zeros(data.shape, dtype=bool)

    # Compute rough magnitudes from flux
    flux = catalog['flux']
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = np.where(flux > 0, -2.5 * np.log10(flux), 99.0)

    labels = np.unique(segmap)
    labels = labels[labels > 0]

    # O(N) pixel count per label via bincount (single image scan)
    label_counts = np.bincount(segmap.ravel())

    # Identify bright sources and compute dilation radii
    bright_labels = []
    radii = []
    n_protected = 0
    n = min(len(labels), len(mag))
    pixel_area_arcsec2 = pixel_scale * pixel_scale
    for i in range(n):
        if mag[i] >= mag_threshold:
            continue
        lab = labels[i]
        npix = label_counts[lab] if lab < len(label_counts) else 0
        if npix == 0:
            continue

        # LSB structure protection (Greco+2018 style)
        if lsb_protect and npix > 0 and flux[i] > 0:
            area_arcsec2 = npix * pixel_area_arcsec2
            mu_mean = (mag_zeropoint - 2.5 * np.log10(flux[i])
                       + 2.5 * np.log10(area_arcsec2))
            if mu_mean >= lsb_mu_threshold:
                n_protected += 1
                continue  # Skip masking — potential LSBG

        bright_labels.append(lab)
        r_eq = np.sqrt(npix / np.pi)
        radii.append(max(1, int(r_eq * expand_factor)))

    if not bright_labels:
        return mask, bright_labels

    # Group labels by dilation radius to batch dilations
    radius_groups = defaultdict(list)
    for lab, r in zip(bright_labels, radii):
        radius_groups[r].append(lab)

    # Merge into at most ~15 groups for efficiency
    sorted_radii = sorted(radius_groups.keys())
    if len(sorted_radii) > 15:
        all_radii = np.array(radii)
        bin_edges = np.percentile(all_radii, np.linspace(0, 100, 16))
        bin_edges = np.unique(np.round(bin_edges).astype(int))
        merged = defaultdict(list)
        merged_r = {}
        for r in sorted_radii:
            b = int(np.searchsorted(bin_edges, r, side='right'))
            merged[b].extend(radius_groups[r])
            merged_r[b] = max(merged_r.get(b, 0), r)
        radius_groups = {merged_r[b]: labs for b, labs in merged.items()}

    # Dilate each group
    from scipy.ndimage import distance_transform_edt
    n_groups = len(radius_groups)
    for gi, (r, labs) in enumerate(sorted(radius_groups.items())):
        submask = np.isin(segmap, labs)
        if r > 30:
            # Large radius: distance transform is O(N) regardless of r
            dist = distance_transform_edt(~submask)
            mask |= (dist <= r)
        else:
            # Small radius: binary_dilation is fine
            y, x = np.ogrid[-r:r + 1, -r:r + 1]
            struct = (x * x + y * y) <= r * r
            mask |= binary_dilation(submask, structure=struct)
        print(f"  dilation group {gi + 1}/{n_groups} "
              f"(r={r}, {len(labs)} sources)", file=sys.stderr)

    if n_protected > 0:
        print(f"  LSB protection: {n_protected} sources preserved "
              f"(mu > {lsb_mu_threshold})", file=sys.stderr)

    return mask, [int(l) for l in bright_labels]


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

    # Magnitude-based selective masking (with LSB protection)
    mask, bright_labels = mask_by_magnitude(
        data, segmap, objects,
        mag_threshold=config.mask_mag_threshold,
        expand_factor=config.mask_expand_factor,
        lsb_protect=config.lsb_protect,
        lsb_mu_threshold=config.lsb_mu_threshold,
        pixel_scale=config.pixel_scale,
        mag_zeropoint=config.mag_zeropoint,
    )
    print(f"LSBG mask: masked {len(bright_labels)} bright sources "
          f"(mag < {config.mask_mag_threshold})", file=sys.stderr)

    # Bright star masking (apply zeropoint for calibrated magnitudes)
    if len(objects) > 0:
        catalog = {
            'x_image': objects['x'] + 1.0,
            'y_image': objects['y'] + 1.0,
            'mag_auto': -2.5 * np.log10(np.maximum(objects['flux'], 1e-30))
                        + config.mag_zeropoint
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
