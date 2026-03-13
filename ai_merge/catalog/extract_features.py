"""Extract 20 numerical features for merge/separate classification."""

import numpy as np


def extract_pair_features(i, j, sep_result):
    """
    Extract 20 features for a candidate pair (i, j).

    Parameters
    ----------
    i, j : int
        Indices into sep_result arrays.
    sep_result : dict
        Output from sep_runner.run_sep_on_stamp or run_sep_on_fits.

    Returns
    -------
    features : np.ndarray
        20-element feature vector, or None if extraction fails.
    names : list of str
        Feature names (same order every call).
    """
    obj = sep_result['objects']
    segmap = sep_result['segmap']
    data_sub = sep_result['data_sub']
    flux_auto = sep_result['flux_auto']
    kron_radius = sep_result['kron_radius']
    flux_radius = sep_result['flux_radius']
    ellipticity = sep_result['ellipticity']
    class_star = sep_result['class_star']

    try:
        feats = _compute_features(i, j, obj, segmap, data_sub,
                                  flux_auto, kron_radius, flux_radius,
                                  ellipticity, class_star)
        return feats, FEATURE_NAMES
    except Exception:
        return None, FEATURE_NAMES


FEATURE_NAMES = [
    # Spatial (5)
    'sep_ratio',         # separation / (R_i + R_j)
    'overlap_frac',      # fraction of overlapping segmentation pixels
    'angle_to_line',     # angle between PA and connecting line
    'pa_alignment',      # |PA_i - PA_j| normalized
    'centroid_offset',   # centroid offset from geometric center
    # Photometric (4)
    'flux_ratio',        # flux_j / flux_i (fainter / brighter)
    'mag_diff',          # |mag_i - mag_j|
    'sb_contrast',       # surface brightness contrast at midpoint
    'bridge_flux_frac',  # flux in bridge region / total flux
    # Morphological (6)
    'ellip_i',           # ellipticity of source i
    'ellip_j',           # ellipticity of source j
    'size_ratio',        # R_j / R_i (smaller / larger)
    'class_star_i',      # star/galaxy class of i
    'class_star_j',      # star/galaxy class of j
    'kron_ratio',        # kron_j / kron_i
    # Pixel-level (5)
    'saddle_depth',      # dip between peaks along connecting line
    'npix_ratio',        # npix_j / npix_i
    'seg_adjacency',     # are segmentation regions adjacent?
    'asymmetry_ratio',   # asymmetry difference
    'concentration_diff', # concentration index difference
]


def _safe_div(a, b, default=0.0):
    """Safe division avoiding div-by-zero."""
    return a / b if abs(b) > 1e-30 else default


def _compute_features(i, j, obj, segmap, data_sub, flux_auto,
                      kron_radius, flux_radius, ellipticity, class_star):
    """Core feature computation. Uses local cutouts for speed on large images."""
    from scipy.ndimage import binary_dilation

    # Ensure i is the brighter source
    if flux_auto[j] > flux_auto[i]:
        i, j = j, i

    xi, yi = obj['x'][i], obj['y'][i]
    xj, yj = obj['x'][j], obj['y'][j]
    ai, bi = obj['a'][i], obj['b'][i]
    aj, bj = obj['a'][j], obj['b'][j]
    ti, tj = obj['theta'][i], obj['theta'][j]

    # Effective radii for separation ratio
    ri = max(kron_radius[i] * ai, flux_radius[i], 1.0)
    rj = max(kron_radius[j] * aj, flux_radius[j], 1.0)

    dist = np.sqrt((xi - xj)**2 + (yi - yj)**2)

    # === Spatial (5) ===
    sep_ratio = _safe_div(dist, ri + rj, 1.0)

    # Local cutout for segmentation analysis (avoid full-image boolean ops)
    h_full, w_full = segmap.shape
    margin = max(ri, rj, dist / 2) + 10
    x0 = max(0, int(min(xi, xj) - margin))
    y0 = max(0, int(min(yi, yj) - margin))
    x1 = min(w_full, int(max(xi, xj) + margin) + 1)
    y1 = min(h_full, int(max(yi, yj) + margin) + 1)

    seg_cut = segmap[y0:y1, x0:x1]
    label_i = i + 1  # segmap labels are 1-indexed
    label_j = j + 1
    mask_i = (seg_cut == label_i)
    mask_j = (seg_cut == label_j)
    npix_i = max(mask_i.sum(), 1)
    npix_j = max(mask_j.sum(), 1)

    # Expand masks by 1px to check adjacency/overlap (on small cutout)
    dilated_i = binary_dilation(mask_i, iterations=1)
    dilated_j = binary_dilation(mask_j, iterations=1)
    overlap = (dilated_i & dilated_j).sum()
    overlap_frac = _safe_div(overlap, npix_i + npix_j)

    # Angle between PA_i and the line connecting centroids
    line_angle = np.arctan2(yj - yi, xj - xi)
    angle_to_line = abs(ti - line_angle)
    angle_to_line = min(angle_to_line, np.pi - angle_to_line) / (np.pi / 2)

    # PA alignment
    pa_diff = abs(ti - tj)
    pa_alignment = min(pa_diff, np.pi - pa_diff) / (np.pi / 2)

    # Centroid offset using local cutout
    mid_x, mid_y = (xi + xj) / 2, (yi + yj) / 2
    combined = mask_i | mask_j
    if combined.any():
        data_cut = data_sub[y0:y1, x0:x1]
        yy, xx = np.where(combined)
        weights = np.maximum(data_cut[combined], 0)
        w_sum = weights.sum()
        if w_sum > 0:
            # Convert local coords back to global
            cx = (xx * weights).sum() / w_sum + x0
            cy = (yy * weights).sum() / w_sum + y0
            centroid_offset = np.sqrt((cx - mid_x)**2 + (cy - mid_y)**2)
            centroid_offset = _safe_div(centroid_offset, dist, 0.0) if dist > 0 else 0.0
        else:
            centroid_offset = 0.0
    else:
        centroid_offset = 0.0

    # === Photometric (4) ===
    flux_ratio = _safe_div(flux_auto[j], flux_auto[i], 0.0)
    flux_ratio = np.clip(flux_ratio, 0.0, 1.0)

    with np.errstate(divide='ignore', invalid='ignore'):
        mag_i = -2.5 * np.log10(max(flux_auto[i], 1e-30)) + 25.0
        mag_j = -2.5 * np.log10(max(flux_auto[j], 1e-30)) + 25.0
    mag_diff = min(abs(mag_i - mag_j), 10.0)  # clamp to 10 mag

    # Surface brightness contrast at midpoint
    h, w = data_sub.shape
    mx, my = int(np.clip(mid_x, 0, w - 1)), int(np.clip(mid_y, 0, h - 1))
    sb_mid = data_sub[my, mx]
    peak_i = obj['peak'][i]
    peak_j = obj['peak'][j]
    sb_contrast = _safe_div(sb_mid, max(peak_i, peak_j), 0.0)

    # Bridge flux fraction: flux along the connecting line between peaks
    bridge_flux_frac = _compute_bridge_frac(data_sub, xi, yi, xj, yj, ri, rj)
    bridge_flux_frac = np.clip(bridge_flux_frac, 0.0, 1.0)

    # === Morphological (6) ===
    ellip_i = ellipticity[i]
    ellip_j = ellipticity[j]
    size_ratio = _safe_div(rj, ri, 0.0)
    cs_i = class_star[i]
    cs_j = class_star[j]
    kron_ratio = _safe_div(kron_radius[j], kron_radius[i], 0.0)

    # === Pixel-level (5) ===
    saddle_depth = np.clip(_compute_saddle_depth(data_sub, xi, yi, xj, yj, peak_i, peak_j), 0.0, 1.0)
    npix_ratio = np.clip(_safe_div(npix_j, npix_i, 0.0), 0.0, 5.0)

    # Segmentation adjacency: 1 if segments touch, 0 otherwise
    seg_adjacency = 1.0 if (dilated_i & mask_j).any() or (dilated_j & mask_i).any() else 0.0

    # Asymmetry ratio (use local cutouts)
    asym_i = _compute_asymmetry(data_sub, xi, yi, mask_i, x0, y0)
    asym_j = _compute_asymmetry(data_sub, xj, yj, mask_j, x0, y0)
    asymmetry_ratio = abs(asym_i - asym_j)

    # Concentration difference
    conc_i = _safe_div(flux_radius[i], kron_radius[i] * ai, 0.5)
    conc_j = _safe_div(flux_radius[j], kron_radius[j] * aj, 0.5)
    concentration_diff = abs(conc_i - conc_j)

    return np.array([
        sep_ratio, overlap_frac, angle_to_line, pa_alignment, centroid_offset,
        flux_ratio, mag_diff, sb_contrast, bridge_flux_frac,
        ellip_i, ellip_j, size_ratio, cs_i, cs_j, kron_ratio,
        saddle_depth, npix_ratio, seg_adjacency, asymmetry_ratio, concentration_diff,
    ], dtype=np.float32)


def _compute_bridge_frac(data_sub, x1, y1, x2, y2, r1, r2):
    """Compute fraction of flux in the bridge region between two sources."""
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    if dist < 1.0:
        return 0.0

    n_pts = max(int(dist * 2), 10)
    t = np.linspace(0, 1, n_pts)
    xs = (x1 + t * (x2 - x1)).astype(int)
    ys = (y1 + t * (y2 - y1)).astype(int)

    h, w = data_sub.shape
    valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
    xs, ys = xs[valid], ys[valid]

    if len(xs) == 0:
        return 0.0

    line_flux = data_sub[ys, xs]

    # Bridge is the middle region (exclude source cores)
    frac_start = min(0.3, r1 / dist) if dist > 0 else 0.3
    frac_end = max(0.7, 1.0 - r2 / dist) if dist > 0 else 0.7
    n = len(line_flux)
    bridge_mask = np.zeros(n, dtype=bool)
    bridge_mask[int(n * frac_start):int(n * frac_end)] = True

    if bridge_mask.sum() == 0:
        return 0.0

    bridge_sum = float(max(line_flux[bridge_mask].sum(), 0))
    total_sum = float(max(abs(line_flux.sum()), 1e-30))
    return min(bridge_sum / total_sum, 1.0)


def _compute_saddle_depth(data_sub, x1, y1, x2, y2, peak1, peak2):
    """Compute the depth of the saddle point between two peaks."""
    dist = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    if dist < 1.0:
        return 0.0

    n_pts = max(int(dist * 2), 10)
    t = np.linspace(0, 1, n_pts)
    xs = (x1 + t * (x2 - x1)).astype(int)
    ys = (y1 + t * (y2 - y1)).astype(int)

    h, w = data_sub.shape
    valid = (xs >= 0) & (xs < w) & (ys >= 0) & (ys < h)
    xs, ys = xs[valid], ys[valid]

    if len(xs) == 0:
        return 0.0

    profile = data_sub[ys, xs]
    min_val = profile.min()
    max_peak = max(peak1, peak2)
    if max_peak <= 0:
        return 0.0

    return max(0, 1.0 - min_val / max_peak)


def _compute_asymmetry(data_sub, cx, cy, mask, cut_x0=0, cut_y0=0):
    """Compute rotational asymmetry of a source.
    mask is in local cutout coords; cut_x0/cut_y0 are cutout origin in global coords.
    """
    if not mask.any():
        return 0.0

    yy, xx = np.where(mask)
    if len(xx) < 4:
        return 0.0

    # Cutout bounding box (local coords within the cutout)
    y0, y1 = yy.min(), yy.max() + 1
    x0, x1 = xx.min(), xx.max() + 1
    # data_sub is global; convert local cutout coords to global
    gy0, gy1 = y0 + cut_y0, y1 + cut_y0
    gx0, gx1 = x0 + cut_x0, x1 + cut_x0
    h, w = data_sub.shape
    gy0, gy1 = max(0, gy0), min(h, gy1)
    gx0, gx1 = max(0, gx0), min(w, gx1)
    if gy1 <= gy0 or gx1 <= gx0:
        return 0.0
    cutout = data_sub[gy0:gy1, gx0:gx1].copy()

    # 180-degree rotated version
    rotated = cutout[::-1, ::-1]

    # Match shapes
    h, w = cutout.shape
    rh, rw = rotated.shape
    mh, mw = min(h, rh), min(w, rw)
    diff = np.abs(cutout[:mh, :mw] - rotated[:mh, :mw])

    total = np.abs(cutout[:mh, :mw]).sum()
    if total < 1e-30:
        return 0.0

    return diff.sum() / total


def extract_pair_features_from_objects(i, j, objects, segmap, data_sub,
                                       flux_auto, kron_radius, flux_radius,
                                       ellipticity, class_star):
    """
    Extract features using raw arrays (for evaluate.py compatibility).
    Constructs a sep_result-like dict and delegates.
    """
    sep_result = {
        'objects': objects,
        'segmap': segmap,
        'data_sub': data_sub,
        'flux_auto': flux_auto,
        'kron_radius': kron_radius,
        'flux_radius': flux_radius,
        'ellipticity': ellipticity,
        'class_star': class_star,
    }
    return extract_pair_features(i, j, sep_result)
