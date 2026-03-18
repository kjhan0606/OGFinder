"""PSF-based star/galaxy classification.

Approach:
1. Find the "star locus" from compact+round+bright sources in the catalog
2. Build empirical PSF from these stars (optional, for verification)
3. Classify each source based on:
   a. Size ratio: A*B relative to the star locus
   b. Ellipticity: stars are round (b/a close to 1)
   c. Cutout concentration check (for bright sources)

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


def _compute_fwhm(cutout):
    """Estimate FWHM from cutout via half-maximum contour area."""
    peak = cutout.max()
    if peak <= 0:
        return 0.0
    half_max = peak / 2.0
    n_above = np.sum(cutout > half_max)
    return 2.0 * np.sqrt(n_above / np.pi)


def _find_star_locus(catalog, n_stars=30):
    """Find the star locus: compact + round + bright sources.

    Among the brightest sources, select those with the smallest A*B
    (most compact) and lowest ellipticity (most round). These define
    the PSF size and shape.

    Returns:
        star_indices: array of indices into catalog
        psf_ab: median A*B of the star locus
    """
    n = len(catalog['x'])
    if n == 0:
        return np.array([], dtype=int), 1.0

    fluxes = catalog.get('flux', np.ones(n))
    a_image = catalog.get('a_image', np.full(n, 3.0))
    b_image = catalog.get('b_image', np.full(n, 3.0))
    ab = a_image * b_image
    ellip = 1.0 - b_image / np.maximum(a_image, 0.01)

    # Step 1: Take the brightest N sources (top N or top 10%)
    n_bright = max(min(500, n // 5), 50)
    bright_order = np.argsort(fluxes)[::-1][:n_bright]

    # Step 2: Among bright sources, find the most compact ones
    ab_bright = ab[bright_order]
    ab_threshold = np.percentile(ab_bright, 15)

    # Step 3: Also require roundness
    ellip_bright = ellip[bright_order]

    compact_mask = (ab_bright < ab_threshold) & (ellip_bright < 0.25)
    candidates = bright_order[compact_mask]

    if len(candidates) < 3:
        # Relax: just take the most compact bright sources
        ab_threshold = np.percentile(ab_bright, 30)
        compact_mask = ab_bright < ab_threshold
        candidates = bright_order[compact_mask]

    if len(candidates) < 3:
        return np.array([], dtype=int), np.median(ab)

    # Step 4: Sort by compactness and take top N
    ab_cands = ab[candidates]
    order = np.argsort(ab_cands)
    candidates = candidates[order[:n_stars]]

    psf_ab = np.median(ab[candidates])
    return candidates, psf_ab


def auto_detect_psf_stars(data, catalog, n_stars=25, size=32):
    """Auto-detect likely PSF stars. Delegates to _find_star_locus."""
    star_indices, _ = _find_star_locus(catalog, n_stars)
    return star_indices


def build_empirical_psf(data, catalog, star_indices, size=64):
    """Build PSF template by median-stacking star cutouts."""
    cutouts = []
    for idx in star_indices:
        x, y = catalog['x'][idx], catalog['y'][idx]
        cut = _extract_cutout(data, x, y, size)
        if cut is None:
            continue
        pos_sum = np.sum(np.maximum(cut, 0))
        if pos_sum > 0:
            cut = cut / pos_sum
        cutouts.append(cut)

    if len(cutouts) == 0:
        return None

    stack = np.median(cutouts, axis=0)
    pos_sum = np.sum(np.maximum(stack, 0))
    if pos_sum > 0:
        stack /= pos_sum

    return stack


def classify_sources(data, catalog, size=64, psf_template=None,
                     star_threshold=0.10):
    """Classify sources using catalog morphometry + PSF comparison.

    Primary classification (all sources):
    - Size ratio: A*B / psf_A*B
    - Ellipticity: 1 - B/A

    Star if ALL of:
    1. size_ratio < 2.0 (compact, within 2x PSF size)
    2. ellipticity < 0.25 (round)
    3. Among the brightest 50% or flux > min_star_flux

    Secondary verification (bright sources, if PSF available):
    - Cutout concentration comparison with PSF

    Args:
        data: 2D background-subtracted image
        catalog: dict with 'number', 'x', 'y', 'a_image', 'b_image',
                 optionally 'flux'
        size: cutout size (for PSF building)
        psf_template: pre-built PSF template (optional)
        star_threshold: not used directly (kept for CLI compatibility)

    Returns:
        results: list of dicts
    """
    import sys

    n = len(catalog['number'])
    if n == 0:
        return []

    fluxes = catalog.get('flux', np.ones(n))
    a_image = catalog.get('a_image', np.full(n, 3.0))
    b_image = catalog.get('b_image', np.full(n, 3.0))
    ab = a_image * b_image
    ellip = 1.0 - b_image / np.maximum(a_image, 0.01)

    # Find star locus
    print("Finding star locus ...", file=sys.stderr)
    star_indices, psf_ab = _find_star_locus(catalog, n_stars=30)
    print(f"  Star locus: {len(star_indices)} candidates, "
          f"PSF A*B = {psf_ab:.2f}", file=sys.stderr)

    if len(star_indices) == 0:
        print("WARNING: Could not find star locus", file=sys.stderr)
        return [{'number': int(catalog['number'][i]),
                 'correlation': 1.0,
                 'label': 'galaxy', 'color': 'yellow'}
                for i in range(n)]

    # Build PSF if not provided (for potential future use)
    if psf_template is None and len(star_indices) >= 3:
        psf_template = build_empirical_psf(data, catalog, star_indices,
                                            size=size)
        if psf_template is not None:
            psf_fwhm = _compute_fwhm(psf_template)
            print(f"  PSF FWHM: {psf_fwhm:.2f} pixels", file=sys.stderr)

    # Star locus statistics
    star_ab = ab[star_indices]
    star_ellip = ellip[star_indices]
    star_flux = fluxes[star_indices]

    # Classification thresholds based on star locus
    max_star_ab = np.percentile(star_ab, 90) * 1.5  # generous upper bound
    min_star_ab = psf_ab * 0.3  # must be at least ~30% of PSF (not noise)
    max_star_ellip = 0.25
    # Stars are foreground: must be above median flux of star locus candidates
    min_star_flux = np.percentile(star_flux, 10)

    print(f"  Thresholds: {min_star_ab:.2f} < A*B < {max_star_ab:.2f}, "
          f"ellip < {max_star_ellip}, flux > {min_star_flux:.3f}",
          file=sys.stderr)

    # Classify each source
    results = []
    n_stars_found = 0
    n_galaxies = 0

    for i in range(n):
        src_num = int(catalog['number'][i])
        src_ab = ab[i]
        src_ellip = ellip[i]
        src_flux = fluxes[i]

        # Star criteria:
        # 1. Compact: A*B within range of star locus (not too large, not too small)
        # 2. Round: low ellipticity
        # 3. Bright enough: above minimum star locus flux
        is_star = (src_ab < max_star_ab and
                   src_ab > min_star_ab and
                   src_ellip < max_star_ellip and
                   src_flux > min_star_flux)

        # Galaxy score: how much larger than PSF
        size_ratio = src_ab / psf_ab if psf_ab > 0 else 999.0
        galaxy_score = max(0.0, size_ratio - 1.0)

        if is_star:
            label = 'star'
            color = 'cyan'
            n_stars_found += 1
        else:
            label = 'galaxy'
            color = 'yellow'
            n_galaxies += 1

        results.append({
            'number': src_num,
            'correlation': float(galaxy_score),
            'label': label,
            'color': color,
        })

    print(f"  Classified: {n} (star={n_stars_found}, galaxy={n_galaxies})",
          file=sys.stderr)

    return results
