"""Building blocks for a mass-independent BCG+diffuse pilot.

The primary estimator is an equal-cluster mean of native-pixel annular
profiles. Area-overlap images are a separate diagnostic on the input image
axes; they are not a measurement of aligned cluster morphology. This module
does not certify sky subtraction, PSF wings, selection, or ICL separation.
"""
import math

import numpy as np
from scipy import sparse

from ogfkit.wlcluster import ARCSEC, Cosmology


def _positive(value, name):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError('%s must be finite and positive' % name)
    return value


def calibrate_flux(image, calibration):
    """Convert a verified per-pixel calibration to extinction-corrected nJy.

    The AB zeropoint must refer to the actual stored pixel values, including
    any exposure normalization already performed. There is no exposure-time
    guess and no default zeropoint or extinction.
    """
    if not calibration.get('evidence') or not calibration.get('extinction_evidence'):
        raise ValueError('calibration and extinction evidence are required')
    extinction = float(calibration['a_i_mag'])
    if not np.isfinite(extinction) or extinction < 0:
        raise ValueError('a_i_mag must be finite and non-negative')
    unit = calibration.get('input_unit')
    if unit == 'nJy/pixel':
        factor = 1.0
    elif unit == 'instrumental_flux/pixel':
        zp = float(calibration['ab_zeropoint'])
        if not np.isfinite(zp):
            raise ValueError('ab_zeropoint must be finite')
        factor = 3631e9 * 10.0 ** (-0.4 * zp)
    else:
        raise ValueError('input_unit must explicitly name flux per pixel')
    factor *= 10.0 ** (0.4 * extinction)
    if not np.isfinite(factor) or factor <= 0:
        raise ValueError('invalid calibration factor')
    return np.asarray(image, float) * factor, factor


def _geometry(shape, center, pixel_scale, z, cosmo=None):
    scale = _positive(pixel_scale, 'pixel_scale')
    z = _positive(z, 'redshift')
    distance = float((cosmo or Cosmology()).angular_diameter(z))
    kpc_per_pixel = _positive(distance * ARCSEC * 1000.0 * scale, 'kpc_per_pixel')
    center = np.asarray(center, float)
    if center.shape != (2,) or not np.isfinite(center).all():
        raise ValueError('center must be finite 0-based x,y')
    yy, xx = np.indices(shape)
    return np.hypot(xx - center[0], yy - center[1]) * kpc_per_pixel, kpc_per_pixel


def native_profile(image_njy, mask, center, pixel_scale, z, edges_kpc,
                   cosmo=None, min_coverage=0.5):
    """Open-area annular SB, with cropping included in the coverage denominator.

    Radial membership uses native pixel centres. The denominator for coverage
    is the full geometric annulus, not just the part inside the input image.
    Pixel errors are not treated as independent survey errors.
    """
    image = np.asarray(image_njy, float)
    mask = np.asarray(mask, bool)
    if image.ndim != 2 or mask.shape != image.shape:
        raise ValueError('image and mask must share a 2-D shape')
    edges = np.asarray(edges_kpc, float)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or edges[0] < 0 or np.any(np.diff(edges) <= 0):
        raise ValueError('radial edges must be non-negative and increasing')
    if not 0 < min_coverage <= 1:
        raise ValueError('min_coverage must lie in (0,1]')
    radius, kpp = _geometry(image.shape, center, pixel_scale, z, cosmo)
    valid = ~mask & np.isfinite(image)
    sb = image / float(pixel_scale) ** 2
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        annulus = (radius >= lo) & (radius < hi)
        use = annulus & valid
        n = int(np.count_nonzero(use))
        expected = math.pi * (hi ** 2 - lo ** 2) / kpp ** 2
        coverage = min(1.0, n / expected)
        rows.append({
            'r_kpc_min': float(lo), 'r_kpc_max': float(hi),
            'sb_njy_arcsec2': float(np.mean(sb[use])) if n and coverage >= min_coverage else None,
            'n_valid_pixels': n, 'n_image_pixels': int(np.count_nonzero(annulus)),
            'valid_area_arcsec2': float(n * float(pixel_scale) ** 2),
            'coverage': float(coverage),
        })
    return rows


def _overlaps(n_input, centre, axis, step):
    """Exact 1-D overlap of output cells with native unit-width pixels."""
    rows, cols, vals = [], [], []
    for i, offset in enumerate(axis):
        left = centre + offset - step / 2.0
        right = left + step
        first = max(0, int(math.floor(left + 0.5)))
        last = min(n_input - 1, int(math.ceil(right + 0.5)) - 1)
        for j in range(first, last + 1):
            length = min(right, j + 0.5) - max(left, j - 0.5)
            if length > 0:
                rows.append(i); cols.append(j); vals.append(length)
    return sparse.csr_matrix((vals, (rows, cols)), shape=(len(axis), n_input))


def area_warp_metric(image_njy, mask, center, pixel_scale, z, *,
                     kpc_max=500.0, n_side=201, cosmo=None):
    """Area-average square pixels; masked flux never enters the numerator.

    Input pixels must be locally square with a checked WCS scale. Axes follow
    the input image, not celestial north or a BCG position angle. Coverage is
    the unmasked fraction of the *whole* output cell, including off-image area.
    No resampled variance or independent-pixel uncertainty is returned.
    """
    image = np.asarray(image_njy, float)
    mask = np.asarray(mask, bool)
    if image.ndim != 2 or mask.shape != image.shape:
        raise ValueError('image and mask must share a 2-D shape')
    if int(n_side) != n_side or n_side < 3:
        raise ValueError('n_side must be an integer >= 3')
    limit = _positive(kpc_max, 'kpc_max')
    # Avoid allocating a full native radius plane for this operation.
    _radius, kpp = _geometry((1, 1), center, pixel_scale, z, cosmo)
    axis = np.linspace(-limit, limit, int(n_side))
    offsets = axis / kpp
    step = float((axis[1] - axis[0]) / kpp)
    wx = _overlaps(image.shape[1], float(center[0]), offsets, step)
    wy = _overlaps(image.shape[0], float(center[1]), offsets, step)
    good = (~mask & np.isfinite(image)).astype(float)
    numerator = wx.dot(wy.dot(np.where(good > 0, image / float(pixel_scale) ** 2, 0.0)).T).T
    area = wx.dot(wy.dot(good).T).T
    mean = np.full(area.shape, np.nan)
    np.divide(numerator, area, out=mean, where=area > 0)
    return {'sb': mean, 'coverage': np.clip(area / step ** 2, 0.0, 1.0),
            'kpc': axis, 'output_area_arcsec2': (step * float(pixel_scale)) ** 2,
            'axis_orientation': 'input-image axes; radial diagnostic only'}


def detection_mask(image, center, pixel_scale, bad_mask, bright_mask, *,
                   kron_factor=2.5, bcg_policy='keep', member_xy=(), threshold=1.5):
    """Mask every detected source except the explicitly identified BCG.

    SEP's background is used ONLY for detection and Kron radii. The caller
    measures the original calibrated image. Unmatched member positions get
    fallback circles. Segmentation footprints supplement Kron ellipses.
    A blended neighbour may still obscure BCG pixels; that needs inspection.
    """
    import sep
    image = np.asarray(image, float)
    bad = np.asarray(bad_mask, bool) | np.asarray(bright_mask, bool) | ~np.isfinite(image)
    if image.ndim != 2 or bad.shape != image.shape:
        raise ValueError('mask shapes must match the 2-D image')
    if bcg_policy not in ('keep', 'mask', 'none'):
        raise ValueError('bcg_policy must be keep, mask or none')
    factor = _positive(kron_factor, 'kron_factor')
    if not np.any(~bad):
        raise ValueError('no usable pixels')
    work = np.ascontiguousarray(np.where(bad, np.median(image[~bad]), image), dtype=np.float64)
    background = sep.Background(work, mask=bad)
    residual = np.ascontiguousarray(work - np.asarray(background), dtype=np.float64)
    rms = _positive(background.globalrms, 'SEP background RMS')
    objects, segmentation = sep.extract(residual, threshold, err=rms, mask=bad, segmentation_map=True)
    match = max(1.0 / _positive(pixel_scale, 'pixel_scale'), 1.5)
    distances = np.hypot(objects['x'] - center[0], objects['y'] - center[1])
    bcg = int(np.argmin(distances)) if len(objects) else -1
    if bcg_policy != 'none' and (bcg < 0 or distances[bcg] > match):
        raise ValueError('BCG detection is missing or unmatched; inspect the frame')
    if bcg_policy == 'none':
        bcg = -1
    if len(objects):
        kron, flags = sep.kron_radius(residual, objects['x'], objects['y'], objects['a'],
                                      objects['b'], objects['theta'], 6.0, mask=bad)
    else:
        kron, flags = np.array([]), np.array([], dtype=int)
    mask = bad.copy()
    mask |= (segmentation > 0) & ((segmentation != bcg + 1) if bcg_policy == 'keep' else True)
    yy, xx = np.indices(image.shape)
    table = []
    for i, obj in enumerate(objects):
        keep = i == bcg and bcg_policy == 'keep'
        scale = factor * float(kron[i]) if not flags[i] and np.isfinite(kron[i]) and kron[i] > 0 else factor * 1.2
        a, b = max(1.0, scale * float(obj['a'])), max(1.0, scale * float(obj['b']))
        if not keep:
            # Only allocate a small local patch per source.
            x, y, theta = float(obj['x']), float(obj['y']), float(obj['theta'])
            reach = int(math.ceil(max(a, b))) + 1
            x0, x1 = max(0, int(x) - reach), min(image.shape[1], int(x) + reach + 1)
            y0, y1 = max(0, int(y) - reach), min(image.shape[0], int(y) + reach + 1)
            dx, dy = xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y
            ct, st = math.cos(theta), math.sin(theta)
            mask[y0:y1, x0:x1] |= ((dx * ct + dy * st) / a) ** 2 + ((-dx * st + dy * ct) / b) ** 2 <= 1
        table.append({'x': float(obj['x']), 'y': float(obj['y']), 'a_mask': a,
                      'b_mask': b, 'bcg': bool(i == bcg), 'masked': not keep,
                      'kron_flag': int(flags[i])})
    fallback_count = 0
    for x, y in member_xy:
        if bcg_policy == 'keep' and math.hypot(x - center[0], y - center[1]) <= match:
            continue
        dist = np.hypot(objects['x'] - x, objects['y'] - y)
        if len(dist) and np.min(dist) <= match:
            continue
        radius = max(2.0 / float(pixel_scale), 1.5)
        mask |= (xx - x) ** 2 + (yy - y) ** 2 <= radius ** 2
        fallback_count += 1
    return mask, {'detections': table, 'n_detected': len(objects),
                  'n_masked_sources': len(objects) - int(bcg_policy == 'keep'),
                  'n_member_fallback': fallback_count, 'bcg_policy': bcg_policy,
                  'detection_background_subtracted_from_science': False}


def profile_matrix(profiles):
    return np.array([[np.nan if row['sb_njy_arcsec2'] is None else row['sb_njy_arcsec2']
                      for row in profile] for profile in profiles], dtype=float)


def _mean(matrix):
    finite = np.isfinite(matrix)
    count = np.sum(finite, axis=0)
    result = np.full(matrix.shape[1], np.nan)
    np.divide(np.sum(np.where(finite, matrix, 0.0), axis=0), count, out=result, where=count > 0)
    return result, count


def paired_summary(cluster_profiles, random_profiles, region_ids):
    """Equal-cluster mean with equal weight for each lens's random-field mean.

    random_profiles is a list of lists, one matched random set per cluster.
    Regional delete-one jackknife recomputes BOTH means in each replicate.
    Covariance is restricted to bins supported in every replicate. It is a
    pilot estimate, not proof of sky transfer, PSF, or sample completeness.
    """
    if not cluster_profiles or len(random_profiles) != len(cluster_profiles) or len(region_ids) != len(cluster_profiles):
        raise ValueError('clusters, matched random sets and regions must align')
    cluster = profile_matrix(cluster_profiles)
    unpaired_raw, unpaired_count = _mean(cluster)
    random = np.full(cluster.shape, np.nan)
    for i, profiles in enumerate(random_profiles):
        if profiles:
            values = profile_matrix(profiles)
            if values.shape[1] != cluster.shape[1]:
                raise ValueError('random profile bins do not match')
            # Every random field in the matched set must cover a bin.
            complete = np.all(np.isfinite(values), axis=0)
            random[i, complete] = np.mean(values[:, complete], axis=0)
    paired = np.isfinite(cluster) & np.isfinite(random)
    cmean, count = _mean(np.where(paired, cluster, np.nan))
    rmean, _ = _mean(np.where(paired, random, np.nan))
    corrected = np.where(paired, cluster - random, np.nan)
    corrected_mean, _ = _mean(corrected)
    covariance = np.full((cluster.shape[1], cluster.shape[1]), np.nan)
    labels = [None if r is None or str(r).strip() == '' else str(r) for r in region_ids]
    groups = sorted(set(labels)) if all(r is not None for r in labels) else []
    cov_bins = []
    if len(groups) >= 3:
        samples = np.array([_mean(corrected[np.array(labels) != group])[0] for group in groups])
        good = np.all(np.isfinite(samples), axis=0)
        cov_bins = np.flatnonzero(good).tolist()
        residual = samples[:, good] - np.mean(samples[:, good], axis=0)
        covariance[np.ix_(good, good)] = (len(groups) - 1.0) / len(groups) * residual.T.dot(residual)
    return {'raw_sb': cmean, 'random_sb': rmean, 'corrected_sb': corrected_mean,
            'unpaired_raw_sb': unpaired_raw, 'unpaired_n_clusters': unpaired_count,
            'n_clusters': count, 'paired': paired, 'covariance': covariance,
            'covariance_bins': cov_bins, 'n_regions': len(groups),
            'covariance_method': 'regional delete-one jackknife' if cov_bins else 'not estimated',
            'science_measurement': False}
