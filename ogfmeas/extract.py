"""Threshold detection, 8-connected segmentation, and multi-threshold deblend.

Detection follows the published outline in Bertin & Arnouts 1996: an optional
convolution decides which pixels belong to an object, and the flux and the
second moments are measured on the unfiltered image. Coordinates are 0-based.
The centre of the first pixel is (0, 0). ``a`` and ``b`` are the semi-axes from
the second-moment eigenvalues, with the usual 1/12 pixel variance added to x2
and y2. ``theta`` is radians in [-pi/2, pi/2].

``set_extract_pixstack`` and ``set_sub_object_limit`` record memory limits. Exceeding
them raises RuntimeError with the text ``pixel buffer full`` or ``deblending overflow``.
This module does not translate another extractor's source.
"""
import math

import numpy as np
from scipy import ndimage as ndi

_PIXSTACK = 30_000_000
_SUB_LIMIT = 4096
_STRUCT = np.ones((3, 3), dtype=np.int32)
_DEFAULT_KERNEL = np.array([[1.0, 2.0, 1.0], [2.0, 4.0, 2.0], [1.0, 2.0, 1.0]]) / 16.0
_MISSING = object()

CAT_DTYPE = np.dtype([
    ("thresh", "f8"),
    ("npix", "i4"),
    ("tnpix", "i4"),
    ("xmin", "i4"), ("xmax", "i4"), ("ymin", "i4"), ("ymax", "i4"),
    ("x", "f8"), ("y", "f8"),
    ("x2", "f8"), ("y2", "f8"), ("xy", "f8"),
    ("errx2", "f8"), ("erry2", "f8"), ("errxy", "f8"),
    ("a", "f8"), ("b", "f8"), ("theta", "f8"),
    ("cxx", "f8"), ("cyy", "f8"), ("cxy", "f8"),
    ("cflux", "f8"), ("flux", "f8"),
    ("cpeak", "f8"), ("peak", "f8"),
    ("xcpeak", "f8"), ("ycpeak", "f8"),
    ("xpeak", "i4"), ("ypeak", "i4"),
    ("flag", "i4"),
])


def set_extract_pixstack(n):
    global _PIXSTACK
    _PIXSTACK = max(int(n), 1)


def set_sub_object_limit(n):
    global _SUB_LIMIT
    _SUB_LIMIT = max(int(n), 1)


def _positive_sum(image, footprint):
    return float(np.sum(np.clip(image[footprint], 0, None)))


def _peaks_of(det, cores):
    peaks = []
    for i, comp in enumerate(cores, start=1):
        yy, xx = np.nonzero(comp)
        j = int(np.argmax(det[yy, xx]))
        peaks.append((i, int(yy[j]), int(xx[j])))
    return peaks


def _assign(det, footprint, cores):
    assigned = np.zeros(footprint.shape, dtype=np.int16)
    for i, comp in enumerate(cores, start=1):
        assigned[comp] = i
    unknown = footprint & (assigned == 0)
    peaks = _peaks_of(det, cores)
    guard = 0
    limit = int(footprint.size) + 2
    while unknown.any() and guard < limit:
        guard += 1
        claim = np.zeros(footprint.shape, dtype=np.int16)
        conflict = np.zeros(footprint.shape, dtype=bool)
        grew = False
        for i in range(1, len(cores) + 1):
            border = ndi.binary_dilation(assigned == i, structure=_STRUCT) & unknown
            if not border.any():
                continue
            grew = True
            both = border & (claim > 0)
            conflict |= both
            claim[border & (claim == 0)] = i
        if not grew:
            break
        claim[conflict] = 0
        assigned[claim > 0] = claim[claim > 0]
        if conflict.any():
            cy, cx = np.nonzero(conflict)
            for y, x in zip(cy, cx):
                best, bestd = 1, 1e300
                for i, py, px in peaks:
                    d2 = (int(y) - py) ** 2 + (int(x) - px) ** 2
                    if d2 < bestd:
                        best, bestd = i, d2
                assigned[int(y), int(x)] = best
        unknown = footprint & (assigned == 0)
    if unknown.any():
        cy, cx = np.nonzero(unknown)
        for y, x in zip(cy, cx):
            best, bestd = 1, 1e300
            for i, py, px in peaks:
                d2 = (int(y) - py) ** 2 + (int(x) - px) ** 2
                if d2 < bestd:
                    best, bestd = i, d2
            assigned[int(y), int(x)] = best
    return assigned


def _split(det, orig, footprint, level, minarea, nthresh, mincont, depth, counter):
    area = int(footprint.sum())
    if area < minarea:
        return []
    if counter[0] > _SUB_LIMIT:
        raise RuntimeError("deblending overflow")
    peak = float(np.max(det[footprint]))
    if nthresh <= 1 or depth >= int(nthresh) or not np.isfinite(peak):
        return [footprint]
    base = max(float(level), peak * 1e-6, 1e-12)
    if peak <= base * 1.02:
        return [footprint]
    steps = max(int(nthresh) - int(depth), 1)
    next_level = base * (peak / base) ** (1.0 / steps)
    if next_level <= base:
        next_level = base * 1.05
    if next_level >= peak * 0.999:
        return [footprint]
    high = footprint & (det >= next_level)
    labeled, nlab = ndi.label(high, structure=_STRUCT)
    cores = []
    for k in range(1, int(nlab) + 1):
        comp = labeled == k
        if int(comp.sum()) >= minarea and _positive_sum(orig, comp) > 0:
            cores.append(comp)
    if len(cores) < 2:
        return _split(det, orig, footprint, next_level, minarea, nthresh, mincont, depth + 1, counter)
    assigned = _assign(det, footprint, cores)
    parent_flux = _positive_sum(orig, footprint)
    children = []
    rejected = np.zeros(footprint.shape, dtype=bool)
    for i in range(1, len(cores) + 1):
        child = assigned == i
        flux = _positive_sum(orig, child)
        if int(child.sum()) >= minarea and parent_flux > 0 and flux >= float(mincont) * parent_flux:
            children.append(child)
        else:
            rejected |= child
    if len(children) < 2:
        return _split(det, orig, footprint, next_level, minarea, nthresh, mincont, depth + 1, counter)
    counter[0] += len(children)
    if counter[0] > _SUB_LIMIT:
        raise RuntimeError("deblending overflow")
    brightest = int(np.argmax([_positive_sum(orig, c) for c in children]))
    children[brightest] = children[brightest] | rejected
    leaves = []
    for child in children:
        leaves.extend(_split(det, orig, child, next_level, minarea, nthresh, mincont, depth + 1, counter))
    if len(leaves) < 2:
        return [footprint]
    return leaves


def _moments(ys, xs, values):
    w = np.clip(values, 0, None)
    total = float(w.sum())
    if total <= 0:
        w = np.ones(ys.shape, dtype=np.float64)
        total = float(w.sum())
    x = float(np.sum(xs * w) / total)
    y = float(np.sum(ys * w) / total)
    dx = xs - x
    dy = ys - y
    x2 = float(np.sum(dx * dx * w) / total) + 1.0 / 12.0
    y2 = float(np.sum(dy * dy * w) / total) + 1.0 / 12.0
    xy = float(np.sum(dx * dy * w) / total)
    return x, y, x2, y2, xy


def _ellipse(x2, y2, xy):
    det = x2 * y2 - xy * xy
    if det < 1e-8:
        x2 = max(x2, 1.0 / 12.0)
        y2 = max(y2, 1.0 / 12.0)
        xy = 0.0
        det = x2 * y2
    diff = x2 - y2
    tmp = math.sqrt(max(diff * diff + 4.0 * xy * xy, 0.0))
    a2 = 0.5 * (x2 + y2 + tmp)
    b2 = 0.5 * (x2 + y2 - tmp)
    a = math.sqrt(max(a2, 0.0))
    b = math.sqrt(max(b2, 0.0))
    if b > a:
        a, b = b, a
    theta = 0.5 * math.atan2(2.0 * xy, diff)
    if theta > math.pi / 2:
        theta -= math.pi
    elif theta < -math.pi / 2:
        theta += math.pi
    cxx = y2 / det
    cyy = x2 / det
    cxy = -2.0 * xy / det
    return a, b, theta, cxx, cyy, cxy


def _record(ys, xs, data, filtered, threshold, ny, nx):
    values = data[ys, xs]
    x, y, x2, y2, xy = _moments(ys, xs, values)
    a, b, theta, cxx, cyy, cxy = _ellipse(x2, y2, xy)
    j = int(np.argmax(values))
    jf = int(np.argmax(filtered[ys, xs]))
    if np.ndim(threshold) == 0:
        thr = float(threshold)
    else:
        thr = float(np.median(threshold[ys, xs]))
    flag = 0
    if int(ys.min()) == 0 or int(xs.min()) == 0 or int(ys.max()) == ny - 1 or int(xs.max()) == nx - 1:
        flag = 1
    npix = int(ys.size)
    rec = np.zeros(1, dtype=CAT_DTYPE)[0]
    rec["thresh"] = thr
    rec["npix"] = npix
    rec["tnpix"] = npix
    rec["xmin"] = int(xs.min())
    rec["xmax"] = int(xs.max())
    rec["ymin"] = int(ys.min())
    rec["ymax"] = int(ys.max())
    rec["x"] = x
    rec["y"] = y
    rec["x2"] = x2
    rec["y2"] = y2
    rec["xy"] = xy
    rec["a"] = a
    rec["b"] = b
    rec["theta"] = theta
    rec["cxx"] = cxx
    rec["cyy"] = cyy
    rec["cxy"] = cxy
    rec["flux"] = float(np.nansum(values))
    rec["cflux"] = float(np.nansum(filtered[ys, xs]))
    rec["peak"] = float(values[j])
    rec["cpeak"] = float(filtered[ys[jf], xs[jf]])
    rec["xpeak"] = int(xs[j])
    rec["ypeak"] = int(ys[j])
    rec["xcpeak"] = float(xs[jf])
    rec["ycpeak"] = float(ys[jf])
    rec["flag"] = flag
    return rec


def extract(data, thresh, err=None, mask=None, minarea=5, filter_kernel=_MISSING,
            deblend_nthresh=32, deblend_cont=0.005, clean=True, segmentation_map=False,
            maskthresh=0.0, **kwargs):
    """Detect objects on a background-subtracted image.

    ``thresh`` is in units of ``err`` when ``err`` is given, and an absolute level
    when ``err`` is None. The omitted filter is the 3x3 kernel [[1,2,1],[2,4,2],[1,2,1]]/16.
    ``filter_kernel=None`` or a 1x1 kernel does not filter. ``clean`` is accepted
    and does not delete deblended components. Unknown keywords are ignored so
    callers can pass extractor options this tree does not use.
    """
    del clean, kwargs
    image = np.ascontiguousarray(data, dtype=np.float64)
    if image.ndim != 2:
        raise ValueError("extract data must be 2-D")
    ny, nx = image.shape
    bad = ~np.isfinite(image)
    if mask is not None:
        bad = bad | (np.asarray(mask) > maskthresh)
    work = image.copy()
    work[bad] = 0.0
    if filter_kernel is _MISSING:
        kernel = _DEFAULT_KERNEL
    elif filter_kernel is None:
        kernel = None
    else:
        kernel = np.asarray(filter_kernel, dtype=np.float64)
        if kernel.ndim != 2 or kernel.size <= 1:
            kernel = None
    if kernel is None:
        filtered = work
    else:
        filtered = ndi.convolve(work, kernel, mode="constant", cval=0.0)
    if thresh is None:
        thresh = 0.0
    if err is None:
        level = float(thresh)
    else:
        err_arr = np.asarray(err, dtype=np.float64)
        level = float(thresh) * err_arr
    above = (filtered > level) & ~bad
    n_above = int(np.count_nonzero(above))
    if n_above > _PIXSTACK:
        raise RuntimeError("pixel buffer full")
    minarea = max(int(minarea), 1)
    catalog = np.zeros(0, dtype=CAT_DTYPE)
    segmap = np.zeros((ny, nx), dtype=np.int32)
    if n_above == 0:
        return (catalog, segmap) if segmentation_map else catalog
    labeled, nlab = ndi.label(above, structure=_STRUCT)
    pieces = []
    nthresh = int(deblend_nthresh)
    for k in range(1, int(nlab) + 1):
        ys, xs = np.nonzero(labeled == k)
        if ys.size < minarea:
            continue
        y0, y1 = int(ys.min()), int(ys.max()) + 1
        x0, x1 = int(xs.min()), int(xs.max()) + 1
        foot = np.zeros((y1 - y0, x1 - x0), dtype=bool)
        foot[ys - y0, xs - x0] = True
        if np.ndim(level) == 0:
            base = float(level)
        else:
            base = float(np.median(level[ys, xs]))
            if not np.isfinite(base) or base <= 0:
                base = float(np.max(filtered[ys, xs])) * 1e-3
        counter = [0]
        leaves = _split(
            filtered[y0:y1, x0:x1], work[y0:y1, x0:x1], foot, base,
            minarea, nthresh, float(deblend_cont), 0, counter)
        for leaf in leaves:
            ly, lx = np.nonzero(leaf)
            if ly.size < minarea:
                continue
            pieces.append((ly + y0, lx + x0))
    if not pieces:
        return (catalog, segmap) if segmentation_map else catalog
    pieces.sort(key=lambda p: (int(p[0].min()), int(p[1].min()), int(p[0][0]), int(p[1][0])))
    rows = [_record(ys, xs, work, filtered, level, ny, nx) for ys, xs in pieces]
    catalog = np.array(rows, dtype=CAT_DTYPE)
    if segmentation_map:
        for i, (ys, xs) in enumerate(pieces, start=1):
            segmap[ys, xs] = i
        return catalog, segmap
    return catalog
