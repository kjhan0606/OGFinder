"""PSF-matching star/galaxy classification.

Instead of a pre-trained CNN, this approach:
1. Auto-detects PSF stars from the image itself
2. Builds an empirical PSF template (median stack)
3. Cross-correlates every source cutout with the PSF
4. Stars = high correlation (source matches PSF shape)
   Galaxies = low correlation (source differs from PSF)

Image-adaptive: works on any telescope without pre-training.
No model checkpoint needed.
"""

import numpy as np


def _extract_cutout(data, x, y, size):
    """Extract zero-padded cutout centered on (x, y). 0-based coords."""
    h, w = data.shape
    ix, iy = int(round(x)), int(round(y))
    half = size // 2
    x0, y0 = ix - half, iy - half
    x1, y1 = ix + half, iy + half
    cx0, cy0 = max(0, x0), max(0, y0)
    cx1, cy1 = min(w, x1), min(h, y1)
    if cx1 <= cx0 or cy1 <= cy0:
        return None
    raw = np.zeros((y1 - y0, x1 - x0), dtype=np.float64)
    raw[cy0 - y0:cy1 - y0, cx0 - x0:cx1 - x0] = data[cy0:cy1, cx0:cx1]
    return raw


def _normalize_cutout(cutout):
    """Normalize cutout for correlation: zero-mean, unit-variance."""
    c = cutout.astype(np.float64)
    c = c - np.mean(c)
    std = np.std(c)
    if std > 0:
        c /= std
    return c


def _pearson_correlation(a, b):
    """Pearson correlation between two 2D arrays (both already zero-mean)."""
    n = a.size
    dot = np.sum(a * b)
    return dot / n


def _compute_fwhm(cutout):
    """Estimate FWHM from cutout via half-maximum contour area."""
    peak = cutout.max()
    if peak <= 0:
        return 0.0
    half_max = peak / 2.0
    n_above = np.sum(cutout > half_max)
    # Circular area: n = pi * (fwhm/2)^2 -> fwhm = 2 * sqrt(n/pi)
    return 2.0 * np.sqrt(n_above / np.pi)


def auto_detect_psf_stars(data, catalog, n_stars=25, size=32):
    """Auto-detect likely PSF stars from the image.

    Strategy:
    1. If CLASS_STAR available: use CLASS_STAR > 0.7 as initial filter
    2. Otherwise: use compact + bright criteria
    3. Check FWHM consistency (stars have similar FWHM)
    4. Return top N most reliable star positions

    Args:
        data: 2D background-subtracted image
        catalog: dict with 'x', 'y', and optionally 'flux', 'a_image',
                 'b_image', 'class_star' arrays (0-based coords)
        n_stars: max number of PSF stars to return
        size: cutout size for FWHM measurement

    Returns:
        star_indices: array of indices into catalog
    """
    n = len(catalog['x'])
    if n == 0:
        return np.array([], dtype=int)

    # Extract basic measurements for all sources
    fluxes = catalog.get('flux', np.ones(n))
    a_image = catalog.get('a_image', np.full(n, 3.0))
    b_image = catalog.get('b_image', np.full(n, 3.0))
    class_star = catalog.get('class_star', None)

    # Step 1: Initial candidates
    if class_star is not None and np.any(class_star > 0.5):
        # Use CLASS_STAR as initial filter
        candidates = np.where(class_star > 0.7)[0]
        if len(candidates) < 5:
            candidates = np.where(class_star > 0.5)[0]
    else:
        # Fallback: compact + bright
        compactness = a_image * b_image
        bright = fluxes > np.percentile(fluxes, 50)
        compact = compactness < np.percentile(compactness, 30)
        circularity = b_image / np.maximum(a_image, 0.1)
        circular = circularity > 0.6
        candidates = np.where(bright & compact & circular)[0]

    if len(candidates) < 3:
        # Very few candidates, use top bright compact sources
        compactness = a_image * b_image
        order = np.argsort(compactness)
        bright_mask = fluxes > np.percentile(fluxes, 30)
        candidates = order[bright_mask[order]][:max(10, n_stars)]

    # Step 2: Measure FWHM for each candidate
    fwhms = []
    valid_cands = []
    h, w = data.shape

    for idx in candidates:
        x, y = catalog['x'][idx], catalog['y'][idx]
        # Skip edge sources
        if x < size or x > w - size or y < size or y > h - size:
            continue
        cut = _extract_cutout(data, x, y, size)
        if cut is None:
            continue
        fwhm = _compute_fwhm(cut)
        if fwhm > 0.5:
            fwhms.append(fwhm)
            valid_cands.append(idx)

    if len(valid_cands) < 3:
        return np.array(valid_cands, dtype=int)

    fwhms = np.array(fwhms)
    valid_cands = np.array(valid_cands)

    # Step 3: Find sources with consistent FWHM (the "star locus")
    median_fwhm = np.median(fwhms)
    fwhm_tolerance = max(0.5, median_fwhm * 0.3)
    consistent = np.abs(fwhms - median_fwhm) < fwhm_tolerance
    star_cands = valid_cands[consistent]

    if len(star_cands) < 3:
        # Relax tolerance
        consistent = np.abs(fwhms - median_fwhm) < median_fwhm * 0.5
        star_cands = valid_cands[consistent]

    # Step 4: Sort by flux (prefer bright stars) and take top N
    star_fluxes = fluxes[star_cands]
    order = np.argsort(star_fluxes)[::-1]
    star_cands = star_cands[order[:n_stars]]

    return star_cands


def build_empirical_psf(data, catalog, star_indices, size=64):
    """Build PSF template by median-stacking star cutouts.

    Each star cutout is normalized (divided by its total positive flux)
    before stacking, so bright and faint stars contribute equally.

    Args:
        data: 2D image
        catalog: dict with 'x', 'y'
        star_indices: indices of PSF star candidates
        size: cutout size

    Returns:
        psf: (size, size) float64 normalized PSF template
    """
    cutouts = []
    for idx in star_indices:
        x, y = catalog['x'][idx], catalog['y'][idx]
        cut = _extract_cutout(data, x, y, size)
        if cut is None:
            continue
        # Normalize by positive flux
        pos_sum = np.sum(np.maximum(cut, 0))
        if pos_sum > 0:
            cut = cut / pos_sum
        cutouts.append(cut)

    if len(cutouts) == 0:
        return None

    # Median stack (robust against outliers: cosmic rays, neighbors)
    stack = np.median(cutouts, axis=0)

    # Re-normalize
    pos_sum = np.sum(np.maximum(stack, 0))
    if pos_sum > 0:
        stack /= pos_sum

    return stack


def classify_sources(data, catalog, size=64, psf_template=None,
                     star_threshold=0.6):
    """Classify all sources by PSF cross-correlation.

    Args:
        data: 2D background-subtracted image
        catalog: dict with 'number', 'x', 'y', and optionally
                 'flux', 'a_image', 'b_image', 'class_star'
        size: cutout size for correlation
        psf_template: pre-built PSF template (if None, auto-detect)
        star_threshold: correlation threshold for star classification

    Returns:
        results: list of dicts with:
            'number': source number
            'correlation': PSF correlation score (0-1)
            'concentration': concentration index
            'label': 'star' or 'galaxy'
            'color': marker color
    """
    import sys

    n = len(catalog['number'])
    if n == 0:
        return []

    # Build PSF template if not provided
    if psf_template is None:
        print("Auto-detecting PSF stars ...", file=sys.stderr)
        star_indices = auto_detect_psf_stars(data, catalog, n_stars=25,
                                              size=min(size, 32))
        print(f"  Found {len(star_indices)} PSF star candidates",
              file=sys.stderr)

        if len(star_indices) < 3:
            print("WARNING: Too few PSF stars found, classification "
                  "may be unreliable", file=sys.stderr)
            if len(star_indices) == 0:
                return [{'number': int(catalog['number'][i]),
                         'correlation': 0.0, 'concentration': 0.0,
                         'label': 'galaxy', 'color': 'yellow'}
                        for i in range(n)]

        psf_template = build_empirical_psf(data, catalog, star_indices,
                                            size=size)
        if psf_template is None:
            print("ERROR: Could not build PSF template", file=sys.stderr)
            return []

    # Normalize PSF template for correlation
    psf_norm = _normalize_cutout(psf_template)

    # Classify each source
    results = []
    for i in range(n):
        x, y = catalog['x'][i], catalog['y'][i]
        src_num = int(catalog['number'][i])

        cut = _extract_cutout(data, x, y, size)
        if cut is None:
            results.append({
                'number': src_num,
                'correlation': 0.0,
                'concentration': 0.0,
                'label': 'galaxy',
                'color': 'yellow',
            })
            continue

        # PSF correlation
        cut_norm = _normalize_cutout(cut)
        corr = _pearson_correlation(cut_norm, psf_norm)
        corr = max(0.0, min(1.0, corr))  # clamp to [0, 1]

        # Concentration index
        conc = _concentration_index(cut)

        # Classification
        if corr >= star_threshold:
            label = 'star'
            color = 'cyan'
        else:
            label = 'galaxy'
            color = 'yellow'

        results.append({
            'number': src_num,
            'correlation': float(corr),
            'concentration': float(conc),
            'label': label,
            'color': color,
        })

    return results


def _concentration_index(cutout):
    """Concentration index C = 5 * log10(r80 / r20).

    r_N = radius containing N% of total flux.
    Stars: C ~ 3-4 (concentrated), Galaxies: C ~ 2-3 (extended).
    """
    h, w = cutout.shape
    cy, cx = h / 2.0, w / 2.0
    y, x = np.mgrid[:h, :w].astype(np.float64)
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)

    # Sort pixels by radius
    r_flat = r.ravel()
    flux_flat = cutout.ravel()
    order = np.argsort(r_flat)
    r_sorted = r_flat[order]
    cumflux = np.cumsum(flux_flat[order])

    total_flux = cumflux[-1]
    if total_flux <= 0:
        return 0.0

    cumflux_norm = cumflux / total_flux

    # Find r20 and r80
    r20_idx = np.searchsorted(cumflux_norm, 0.20)
    r80_idx = np.searchsorted(cumflux_norm, 0.80)

    r20 = r_sorted[min(r20_idx, len(r_sorted) - 1)]
    r80 = r_sorted[min(r80_idx, len(r_sorted) - 1)]

    if r20 <= 0 or r80 <= 0:
        return 0.0

    return 5.0 * np.log10(max(r80 / r20, 1.01))
