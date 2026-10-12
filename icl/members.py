"""Separate cluster satellites from background galaxies.

The spectroscopic window is the caller's. Nothing here adopts the Abell S1063
reproduction cuts. A finite redshift inside |z - z_cl| < dz is a satellite. A
finite redshift outside that window is a background galaxy. A missing redshift
stays unknown until photometry is supplied.

Photometry runs only on rows that still have no finite spectroscopic redshift.
A point photo-z or a P(z) uses the caller's cluster window and background window.
A row in one window takes that label. A row in both, or in neither, stays unknown.
The red sequence is a straight line fit to spectroscopic satellites only. The
residual width and the magnitude limit are the caller's. A row still unknown can
become a satellite on that line. It does not become a background galaxy by missing
the line. A finite spectroscopic redshift is left as it is.

A chosen BCG index stays the BCG, and a star flag stays a star, ahead of every
redshift and colour label.

ICL surface brightness hides every resolved source except the BCG. An ICL light
fraction keeps the BCG and the satellites. Weak lensing keeps a row only when this
table says background and the source redshift is finite and behind the lens by the
usual redshift gap. A photometric label does not replace that redshift.
"""

import numpy as np

BCG = 'bcg'
SATELLITE = 'satellite'
BACKGROUND = 'background'
STAR = 'star'
UNKNOWN = 'unknown'

# Resolved light removed before an ICL surface-brightness stack. The BCG stays.
DIFFUSE_MASKED = (SATELLITE, BACKGROUND, STAR, UNKNOWN)
# Light kept in an ICL fraction whose denominator is BCG + satellites + ICL.
FRACTION_KEPT = (BCG, SATELLITE)


def classify_spectroscopic(z, z_cluster, dz, bcg_index=None, star=None):
    """Label each row from a spectroscopic redshift and the caller's window.

    z is one row per object. dz is half the width of the cluster window and must
    be positive. The boundary |z - z_cl| = dz is outside the window.
    """
    redshift = np.asarray(z, float)
    if redshift.ndim != 1:
        raise ValueError('redshift must be one-dimensional')
    dz = float(dz)
    z_cluster = float(z_cluster)
    if not np.isfinite(dz) or dz <= 0:
        raise ValueError('dz must be positive')
    if not np.isfinite(z_cluster) or z_cluster < 0:
        raise ValueError('cluster redshift must be non-negative')
    labels = np.full(redshift.shape, UNKNOWN, dtype='U16')
    finite = np.isfinite(redshift)
    inside = finite & (np.abs(redshift - z_cluster) < dz)
    labels[finite & ~inside] = BACKGROUND
    labels[inside] = SATELLITE
    return _apply_star_and_bcg(labels, star, bcg_index)


def _apply_star_and_bcg(labels, star, bcg_index):
    labels = np.array(labels, dtype='U16', copy=True)
    if star is not None:
        star = np.asarray(star, bool)
        if star.shape != labels.shape:
            raise ValueError('star flag length does not match the catalogue')
        labels[star] = STAR
    if bcg_index is not None:
        index = int(bcg_index)
        if index < 0 or index >= labels.size:
            raise ValueError('bcg index is outside the catalogue')
        labels[index] = BCG
    return labels


def _interval(window, name):
    if window is None or len(window) != 2:
        raise ValueError('%s must be (low, high)' % name)
    low, high = float(window[0]), float(window[1])
    if not np.isfinite(low) or not np.isfinite(high) or not high > low:
        raise ValueError('%s must be a finite increasing interval' % name)
    return low, high


def _in_half_open(z, window):
    low, high = window
    redshift = np.asarray(z, float)
    return np.isfinite(redshift) & (redshift >= low) & (redshift < high)


def _point_labels(z_phot, cluster, background):
    in_cluster = _in_half_open(z_phot, cluster)
    in_background = _in_half_open(z_phot, background)
    labels = np.full(np.shape(in_cluster), UNKNOWN, dtype='U16')
    labels[in_cluster & ~in_background] = SATELLITE
    labels[in_background & ~in_cluster] = BACKGROUND
    return labels


def _interval_integral(z_grid, pz, low, high):
    """Integrate each P(z) row on [low, high) by clipping every grid segment."""
    z = np.asarray(z_grid, float)
    p = np.asarray(pz, float)
    acc = np.zeros(p.shape[0], float)
    for i in range(z.size - 1):
        z0, z1 = float(z[i]), float(z[i + 1])
        if not np.isfinite(z0) or not np.isfinite(z1) or z0 == z1:
            continue
        left = max(min(z0, z1), low)
        right = min(max(z0, z1), high)
        if not np.isfinite(left) or not np.isfinite(right) or right <= left:
            continue
        span = z1 - z0
        t0 = (left - z0) / span
        t1 = (right - z0) / span
        if t1 < t0:
            t0, t1 = t1, t0
        p0, p1 = p[:, i], p[:, i + 1]
        f0 = p0 + (p1 - p0) * t0
        f1 = p0 + (p1 - p0) * t1
        good = np.isfinite(f0) & np.isfinite(f1)
        acc[good] += 0.5 * (f0[good] + f1[good]) * (right - left)
    return acc


def photometric_redshift_labels(z_phot=None, cluster_window=None, background_window=None,
                                pz=None, z_grid=None, p_min=None):
    """Label rows from a point photo-z, a P(z), or both.

    cluster_window and background_window are half-open intervals [low, high).
    A row in exactly one window takes that label. Both and neither stay unknown.
    When P(z) is usable, the row is decided by P(z): one window must hold at least
    p_min of that row and the other must hold less than p_min. An unusable P(z)
    row falls through to the point photo-z.
    """
    cluster = _interval(cluster_window, 'cluster window')
    background = _interval(background_window, 'background window')
    if z_phot is None and pz is None:
        raise ValueError('a point photo-z or P(z) is required')
    if pz is not None:
        if p_min is None:
            raise ValueError('p_min is required with P(z)')
        p_min = float(p_min)
        if not np.isfinite(p_min) or p_min <= 0 or p_min > 1:
            raise ValueError('p_min must be in (0, 1]')
        probability = np.asarray(pz, float)
        grid = np.asarray(z_grid, float)
        if probability.ndim != 2 or grid.ndim != 1 or probability.shape[1] != grid.size:
            raise ValueError('P(z) shape must be (n, len(z_grid))')
        if grid.size < 2:
            raise ValueError('redshift grid needs at least two points')
        n = probability.shape[0]
    else:
        probability = None
        n = np.asarray(z_phot, float).size
    if z_phot is not None:
        z_phot = np.asarray(z_phot, float)
        if z_phot.ndim != 1 or z_phot.size != n:
            raise ValueError('point photo-z length does not match P(z)')
    labels = np.full(n, UNKNOWN, dtype='U16')
    used = np.zeros(n, bool)
    if probability is not None:
        p_cluster = _interval_integral(grid, probability, cluster[0], cluster[1])
        p_background = _interval_integral(grid, probability, background[0], background[1])
        p_total = _interval_integral(grid, probability, -np.inf, np.inf)
        usable = np.isfinite(p_total) & (p_total > 0) & np.isfinite(p_cluster) & np.isfinite(p_background)
        fraction_cluster = np.zeros(n, float)
        fraction_background = np.zeros(n, float)
        fraction_cluster[usable] = p_cluster[usable] / p_total[usable]
        fraction_background[usable] = p_background[usable] / p_total[usable]
        satellite = usable & (fraction_cluster >= p_min) & (fraction_background < p_min)
        background_row = usable & (fraction_background >= p_min) & (fraction_cluster < p_min)
        labels[satellite] = SATELLITE
        labels[background_row] = BACKGROUND
        used[usable] = True
    if z_phot is not None:
        point = _point_labels(z_phot, cluster, background)
        labels[~used] = point[~used]
    return labels


def fit_member_red_sequence(color, mag, labels, mag_limit, min_satellites, width):
    """Fit colour = a + b (mag - m0) to spectroscopic satellites.

    Training rows are labelled satellite, at or brighter than mag_limit, with a
    finite colour. Stars and the BCG are out when their labels say so. Clipping
    uses the caller width: a residual at that width leaves the fit. Too few
    remaining satellites returns ok False and adds no photometric satellites.
    """
    colour = np.asarray(color, float)
    magnitude = np.asarray(mag, float)
    label = np.asarray(labels)
    if colour.ndim != 1 or colour.shape != magnitude.shape or colour.shape != label.shape:
        raise ValueError('colour, magnitude and labels must be one catalogue')
    mag_limit = float(mag_limit)
    width = float(width)
    min_satellites = int(min_satellites)
    if not np.isfinite(mag_limit):
        raise ValueError('magnitude limit must be finite')
    if not np.isfinite(width) or width <= 0:
        raise ValueError('red-sequence width must be positive')
    if min_satellites < 2:
        raise ValueError('min_satellites must be at least 2')
    train = ((label == SATELLITE) & np.isfinite(colour) & np.isfinite(magnitude)
             & (magnitude <= mag_limit))
    use = np.flatnonzero(train)
    empty = dict(ok=False, n=int(use.size), a=float('nan'), b=float('nan'), m0=float('nan'),
                 rms=float('nan'), mag_limit=mag_limit, width=width,
                 min_satellites=min_satellites)
    if use.size < min_satellites:
        return empty
    a = b = m0 = float('nan')
    stable = False
    for _ in range(20):
        dm = magnitude[use] - float(np.median(magnitude[use]))
        design = np.column_stack([np.ones(use.size), dm])
        sol, *_ = np.linalg.lstsq(design, colour[use], rcond=None)
        a, b = float(sol[0]), float(sol[1])
        m0 = float(np.median(magnitude[use]))
        resid = colour[use] - (a + b * (magnitude[use] - m0))
        keep = np.abs(resid) < width
        if int(keep.sum()) < min_satellites:
            empty['n'] = int(keep.sum())
            return empty
        nxt = use[keep]
        if nxt.size == use.size:
            stable = True
            break
        use = nxt
    if not stable:
        dm = magnitude[use] - float(np.median(magnitude[use]))
        design = np.column_stack([np.ones(use.size), dm])
        sol, *_ = np.linalg.lstsq(design, colour[use], rcond=None)
        a, b = float(sol[0]), float(sol[1])
        m0 = float(np.median(magnitude[use]))
        resid = colour[use] - (a + b * (magnitude[use] - m0))
    rms = float(np.sqrt(np.mean(resid ** 2)))
    return dict(ok=True, n=int(use.size), a=a, b=b, m0=m0, rms=rms,
                mag_limit=mag_limit, width=width, min_satellites=min_satellites)


def red_sequence_match(color, mag, fit):
    """True on the fitted line, within the caller width and magnitude limit.

    A false value is not a background label. The boundary |residual| = width is
    outside the sequence.
    """
    colour = np.asarray(color, float)
    magnitude = np.asarray(mag, float)
    if not fit or not fit.get('ok'):
        return np.zeros(colour.shape, bool)
    pred = fit['a'] + fit['b'] * (magnitude - fit['m0'])
    return (np.isfinite(colour) & np.isfinite(magnitude)
            & (magnitude <= fit['mag_limit'])
            & (np.abs(colour - pred) < fit['width']))


def classify_members(z, z_cluster, dz, bcg_index=None, star=None, z_phot=None,
                     cluster_window=None, background_window=None, pz=None, z_grid=None,
                     p_min=None, color=None, mag=None, rs_width=None, mag_limit=None,
                     min_satellites=None):
    """Spectroscopic window, then photo-z, then the red sequence.

    Returns a dict with `labels` and `red_sequence`. The fit is None when no
    colour-magnitude arguments are passed. A finite spectroscopic redshift is
    written first and is not replaced by photometry.
    """
    spec = classify_spectroscopic(z, z_cluster, dz)
    labels = np.array(spec, dtype='U16', copy=True)
    photo_requested = (z_phot is not None or pz is not None or cluster_window is not None
                       or background_window is not None or p_min is not None)
    if photo_requested:
        photo = photometric_redshift_labels(
            z_phot=z_phot, cluster_window=cluster_window, background_window=background_window,
            pz=pz, z_grid=z_grid, p_min=p_min)
        if photo.shape != labels.shape:
            raise ValueError('photometric redshift length does not match the catalogue')
        fill = (labels == UNKNOWN) & (photo != UNKNOWN)
        labels[fill] = photo[fill]
    colour_args = (color, mag, rs_width, mag_limit, min_satellites)
    fit = None
    if any(item is not None for item in colour_args):
        if any(item is None for item in colour_args):
            raise ValueError('red sequence needs colour, magnitude, width, magnitude limit and min_satellites')
        train = np.array(spec, dtype='U16', copy=True)
        train = _apply_star_and_bcg(train, star, bcg_index)
        fit = fit_member_red_sequence(color, mag, train, mag_limit, min_satellites, rs_width)
        on_sequence = red_sequence_match(color, mag, fit)
        fill = (labels == UNKNOWN) & on_sequence
        labels[fill] = SATELLITE
    labels = _apply_star_and_bcg(labels, star, bcg_index)
    return {'labels': labels, 'red_sequence': fit}


def diffuse_stack_masked(labels):
    """True for rows whose pixels are removed from an ICL surface-brightness stack."""
    return np.isin(np.asarray(labels), DIFFUSE_MASKED)


def fraction_kept(labels):
    """True for rows kept in the BCG + satellite part of an ICL light fraction."""
    return np.isin(np.asarray(labels), FRACTION_KEPT)


def background_for_lensing(labels, z, z_lens, z_gap):
    """True for background rows whose redshift is finite and behind the lens.

    The redshift is the source redshift used for Sigma_crit. A photometric
    background label with a missing redshift stays out.
    """
    redshift = np.asarray(z, float)
    if redshift.shape != np.asarray(labels).shape:
        raise ValueError('redshift length does not match the labels')
    return ((np.asarray(labels) == BACKGROUND)
            & np.isfinite(redshift)
            & (redshift > float(z_lens) + float(z_gap)))


def object_mask(segmap, labels, mode):
    """Pixel mask from a segmentation map. Segment 0 is empty sky and stays open.

    mode 'diffuse' masks satellites, background, stars and unknowns.
    mode 'nonmember' masks every row that is not the BCG or a satellite.
    Segment ids are 1-based positions in `labels`.
    """
    if mode == 'diffuse':
        reject = diffuse_stack_masked(labels)
    elif mode == 'nonmember':
        reject = ~fraction_kept(labels)
    else:
        raise ValueError('mode must be diffuse or nonmember')
    segments = np.asarray(segmap, int)
    by_id = np.concatenate([[False], np.asarray(reject, bool)])
    if segments.size and int(segments.max()) >= by_id.size:
        raise ValueError('segmentation id has no label')
    if segments.size and int(segments.min()) < 0:
        raise ValueError('segmentation id must be non-negative')
    return by_id[segments]
