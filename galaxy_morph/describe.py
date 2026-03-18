"""Compact morphological notation from CNN probabilities + image analysis.

Generates de Vaucouleurs-style descriptors like:
    E0, E/S0, S0, Sa, SBb/2/r, Sc/3, Sd(e)/d, M, Irr/t

Notation:
    TYPE[subtype][(e)][/N][/features]

    TYPE:  E (elliptical), S (spiral), SB (barred spiral),
           S0 (lenticular), M (merger), Irr (irregular)
    subtype: 0-7 for E, a/b/c/d for spirals
    (e):   edge-on
    /N:    N spiral arms (1-4+)
    /features: r=ring, R=outer ring, b=bar(weak), d=dust lane,
               t=tidal, c=compact, s=shells
"""

import numpy as np


# Galaxy10 class -> base Hubble type
_BASE_TYPE = {
    0: 'Irr',     # Disturbed
    1: 'M',       # Merging
    2: 'E',       # Round Smooth
    3: 'E/S0',    # In-between Smooth
    4: 'S0',      # Cigar Smooth
    5: 'SB',      # Barred Spiral
    6: 'S',       # Unbarred Tight Spiral
    7: 'S',       # Unbarred Loose Spiral
    8: 'S(e)',     # Edge-on w/o Bulge
    9: 'S(e)',     # Edge-on w/ Bulge
}


def _ellipticity_subtype(cutout):
    """Estimate elliptical subtype E0-E7 from axis ratio."""
    if cutout is None:
        return '0'
    h, w = cutout.shape[-2:]
    cy, cx = h // 2, w // 2
    r = min(h, w) // 3

    # Measure intensity-weighted second moments
    yy, xx = np.mgrid[-r:r+1, -r:r+1]
    stamp = cutout[cy-r:cy+r+1, cx-r:cx+r+1] if cutout.ndim == 2 else \
            cutout[0, cy-r:cy+r+1, cx-r:cx+r+1] if cutout.ndim == 3 else None
    if stamp is None or stamp.size == 0:
        return '0'

    w_map = np.maximum(stamp - np.median(stamp), 0)
    total = w_map.sum()
    if total <= 0:
        return '0'

    Ixx = (w_map * xx**2).sum() / total
    Iyy = (w_map * yy**2).sum() / total
    Ixy = (w_map * xx * yy).sum() / total

    # Eigenvalues of inertia tensor -> axis ratio
    trace = Ixx + Iyy
    det = Ixx * Iyy - Ixy**2
    disc = max(0, trace**2 / 4.0 - det)
    lam1 = trace / 2.0 + np.sqrt(disc)
    lam2 = max(trace / 2.0 - np.sqrt(disc), 0.01)

    q = np.sqrt(lam2 / lam1) if lam1 > 0 else 1.0
    e = 1.0 - q  # ellipticity
    n = min(7, max(0, int(round(10 * e))))
    return str(n)


def _spiral_subtype(class_id, proba):
    """Determine spiral subtype a/b/c/d from class probabilities."""
    # Tight spiral (class 6) -> a or b
    # Loose spiral (class 7) -> c or d
    # Barred spiral (class 5) -> a/b/c based on tight/loose ratio
    p_tight = proba[6] if len(proba) > 6 else 0
    p_loose = proba[7] if len(proba) > 7 else 0
    p_edge_nb = proba[8] if len(proba) > 8 else 0  # edge-on no bulge -> late type
    p_edge_b = proba[9] if len(proba) > 9 else 0   # edge-on with bulge -> early type

    # Tightness score: higher = tighter winding = earlier type
    tight_score = p_tight + p_edge_b * 0.5
    loose_score = p_loose + p_edge_nb * 0.5

    total = tight_score + loose_score + 0.001

    if class_id == 6:  # Unbarred Tight
        if tight_score / total > 0.7:
            return 'a'
        return 'b'
    elif class_id == 7:  # Unbarred Loose
        if loose_score / total > 0.7:
            return 'c'
        return 'b'
    elif class_id == 5:  # Barred Spiral
        ratio = tight_score / total
        if ratio > 0.6:
            return 'a'
        elif ratio > 0.35:
            return 'b'
        return 'c'
    elif class_id == 8:  # Edge-on no bulge -> late
        return 'd'
    elif class_id == 9:  # Edge-on with bulge -> early/mid
        return 'b'
    return ''


def count_spiral_arms(cutout, r_inner_frac=0.25, r_outer_frac=0.75):
    """Count spiral arms via azimuthal Fourier decomposition.

    Args:
        cutout: 2D float array (H, W), asinh-normalized
        r_inner_frac: inner radius fraction of half-size
        r_outer_frac: outer radius fraction of half-size

    Returns:
        n_arms: int (0 if no clear spiral, 1-4 otherwise)
        arm_strength: float (0-1, strength of dominant mode)
    """
    if cutout is None:
        return 0, 0.0

    if cutout.ndim == 3:
        cutout = cutout[0]

    h, w = cutout.shape
    cy, cx = h / 2.0, w / 2.0
    half = min(h, w) / 2.0

    r_inner = half * r_inner_frac
    r_outer = half * r_outer_frac

    # Subtract median background
    img = cutout - np.median(cutout)

    # Sample azimuthal profile at multiple radii
    n_theta = 64
    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)

    # Accumulate Fourier power over radii
    mode_power = np.zeros(5)  # modes m=0,1,2,3,4
    n_radii = 0

    for r in np.linspace(r_inner, r_outer, 12):
        if r < 2:
            continue
        # Sample intensity around annulus
        profile = np.zeros(n_theta)
        for k, theta in enumerate(thetas):
            xi = cx + r * np.cos(theta)
            yi = cy + r * np.sin(theta)
            ix, iy = int(round(xi)), int(round(yi))
            if 0 <= ix < w and 0 <= iy < h:
                profile[k] = img[iy, ix]

        if profile.std() < 1e-8:
            continue

        # FFT
        fft = np.fft.rfft(profile)
        power = np.abs(fft) ** 2
        # Normalize by m=0 power
        if power[0] > 0:
            norm_power = power / power[0]
            for m in range(1, min(5, len(norm_power))):
                mode_power[m] += norm_power[m]
            n_radii += 1

    if n_radii == 0:
        return 0, 0.0

    mode_power /= n_radii

    # Find dominant mode (m=1,2,3,4)
    dominant_m = np.argmax(mode_power[1:]) + 1
    strength = mode_power[dominant_m]

    # Threshold: need significant power above noise
    if strength < 0.05:
        return 0, 0.0

    return int(dominant_m), float(strength)


def detect_ring(cutout, n_bins=20):
    """Detect ring feature from radial intensity profile.

    Returns:
        has_ring: bool
        ring_type: 'r' (inner ring) or 'R' (outer ring) or ''
    """
    if cutout is None:
        return False, ''

    if cutout.ndim == 3:
        cutout = cutout[0]

    h, w = cutout.shape
    cy, cx = h / 2.0, w / 2.0
    half = min(h, w) / 2.0

    # Build radial profile
    yy, xx = np.mgrid[0:h, 0:w]
    rr = np.sqrt((xx - cx)**2 + (yy - cy)**2)

    radii = np.linspace(0, half * 0.9, n_bins + 1)
    profile = np.zeros(n_bins)

    for i in range(n_bins):
        mask = (rr >= radii[i]) & (rr < radii[i+1])
        if mask.sum() > 0:
            profile[i] = np.mean(cutout[mask])

    if profile.max() <= 0:
        return False, ''

    # Normalize
    profile = profile / profile.max()

    # Look for local maxima that deviate from monotonic decline
    # A ring shows as a bump in the radial profile
    for i in range(2, n_bins - 2):
        # Is this a local maximum?
        if profile[i] > profile[i-1] and profile[i] > profile[i+1]:
            # Is it significantly above neighbors?
            baseline = (profile[i-2] + profile[i+2]) / 2.0
            bump = profile[i] - baseline
            if bump > 0.08:  # significant bump
                r_frac = (radii[i] + radii[i+1]) / 2.0 / half
                if r_frac < 0.5:
                    return True, 'r'  # inner ring
                else:
                    return True, 'R'  # outer ring

    return False, ''


def detect_bar_strength(proba):
    """Estimate bar presence from CNN probabilities.

    Returns:
        has_bar: bool
        bar_code: 'B' (strong bar in type) or 'b' (weak bar as feature) or ''
    """
    p_barred = proba[5] if len(proba) > 5 else 0
    p_unbarred = (proba[6] + proba[7]) if len(proba) > 7 else 0

    if p_barred > 0.4:
        return True, 'B'  # strong bar -> SB type
    elif p_barred > 0.15:
        return True, 'b'  # weak bar -> feature
    return False, ''


def detect_tidal(proba):
    """Detect tidal/disturbed features from probabilities."""
    p_disturbed = proba[0] if len(proba) > 0 else 0
    p_merging = proba[1] if len(proba) > 1 else 0
    return (p_disturbed + p_merging) > 0.25


def measure_concentration(cutout):
    """Measure concentration index C = 5 * log10(r80/r20)."""
    if cutout is None:
        return 0.0

    if cutout.ndim == 3:
        cutout = cutout[0]

    h, w = cutout.shape
    cy, cx = h / 2.0, w / 2.0

    yy, xx = np.mgrid[0:h, 0:w]
    rr = np.sqrt((xx - cx)**2 + (yy - cy)**2)

    # Sort pixels by radius
    flat_r = rr.ravel()
    flat_v = cutout.ravel()
    sort_idx = np.argsort(flat_r)
    sorted_r = flat_r[sort_idx]
    sorted_v = flat_v[sort_idx]

    cumsum = np.cumsum(np.maximum(sorted_v - np.median(sorted_v), 0))
    if cumsum[-1] <= 0:
        return 0.0

    cumsum /= cumsum[-1]

    r20 = sorted_r[np.searchsorted(cumsum, 0.2)]
    r80 = sorted_r[np.searchsorted(cumsum, 0.8)]

    if r20 <= 0:
        r20 = 0.5

    return 5.0 * np.log10(r80 / r20)


def build_descriptor(class_id, proba, cutout=None):
    """Build compact morphological descriptor.

    Args:
        class_id: int, Galaxy10 class index (0-9)
        proba: array of 10 class probabilities
        cutout: 2D or 3D float array (optional, for image analysis)

    Returns:
        descriptor: str like 'SBb/2/r', 'E3', 'Sc/3', 'M/t'
    """
    proba = np.asarray(proba)
    parts = []

    # --- Primary type ---
    base = _BASE_TYPE.get(class_id, '?')
    is_spiral = class_id in (5, 6, 7, 8, 9)
    is_elliptical = class_id in (2, 3)
    is_edge_on = class_id in (8, 9)

    if is_elliptical:
        if class_id == 2:  # Round Smooth
            subtype = _ellipticity_subtype(cutout)
            parts.append(f'E{subtype}')
        else:  # In-between Smooth
            parts.append('E/S0')

    elif class_id == 4:  # Cigar Smooth -> S0
        parts.append('S0')

    elif is_spiral:
        subtype = _spiral_subtype(class_id, proba)

        # Check bar
        has_bar, bar_code = detect_bar_strength(proba)

        if class_id == 5 or bar_code == 'B':
            # Barred spiral
            parts.append(f'SB{subtype}')
        elif is_edge_on:
            parts.append(f'S{subtype}(e)')
        else:
            # Check for weak bar
            if bar_code == 'b':
                parts.append(f'SAB{subtype}')
            else:
                parts.append(f'S{subtype}')

    elif class_id == 1:  # Merging
        parts.append('M')

    elif class_id == 0:  # Disturbed
        parts.append('Irr')

    else:
        parts.append(base)

    # --- Features (only for spirals with cutout) ---
    features = []

    if is_spiral and cutout is not None:
        # Spiral arm count
        n_arms, arm_str = count_spiral_arms(cutout)
        if n_arms >= 1 and arm_str > 0.03:
            features.append(str(n_arms))

        # Ring detection
        has_ring, ring_code = detect_ring(cutout)
        if has_ring:
            features.append(ring_code)

    # Edge-on dust lane proxy: if edge-on, check asymmetry
    if is_edge_on and cutout is not None:
        c = cutout[0] if cutout.ndim == 3 else cutout
        h, w = c.shape
        top_half = c[:h//2, :].mean()
        bot_half = c[h//2:, :].mean()
        if abs(top_half - bot_half) / max(top_half, bot_half, 0.01) > 0.15:
            features.append('d')  # dust lane

    # Tidal features
    if detect_tidal(proba) and class_id not in (0, 1):
        features.append('t')

    # Compact indicator
    if cutout is not None:
        conc = measure_concentration(cutout)
        if conc > 4.0 and not is_elliptical:
            features.append('c')

    # Weak bar as feature (if not already SB)
    if is_spiral and class_id != 5:
        _, bar_code = detect_bar_strength(proba)
        if bar_code == 'b' and 'SAB' not in parts[0]:
            features.append('b')

    # Build final descriptor
    descriptor = parts[0]
    if features:
        descriptor += '/' + '/'.join(features)

    return descriptor
