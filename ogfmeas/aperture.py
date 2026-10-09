"""Circular and elliptical aperture sums, Kron radius, flux radius, windowed centroid.

Pixel (row i, column j) is centred at (x=j, y=i) and covers [j-0.5, j+0.5] by [i-0.5, i+0.5].
Sub-pixel samples sit at (k+0.5)/subpix - 0.5. Flag value 1: a sample hit a masked or
non-finite pixel. Flag value 2: the aperture extends outside the image.

``kron_radius`` returns the first moment sum(rho * I) / sum(I) inside the elliptical
limit. It is not multiplied by 2.5. ``theta`` is radians, counter-clockwise from +x.
"""
import math

import numpy as np


def _broadcast(arrays):
    cols = [np.atleast_1d(np.asarray(a, dtype=np.float64)).ravel() for a in arrays]
    n = max(c.size for c in cols)
    out = []
    for c in cols:
        if c.size == n:
            out.append(c)
        elif c.size == 1:
            out.append(np.full(n, c[0]))
        else:
            raise ValueError("aperture coordinate arrays have different lengths")
    return n, out


def _empty(kind):
    z = np.zeros(0, dtype=np.float64)
    f = np.zeros(0, dtype=np.int32)
    if kind == "sum":
        return z, z, f
    if kind == "kron":
        return z, f
    if kind == "win":
        return z, z, f
    raise ValueError(kind)


def _bad_map(data, mask, maskthresh):
    bad = ~np.isfinite(data)
    if mask is not None:
        bad = bad | (np.asarray(mask) > maskthresh)
    return bad


def _sample_offsets(subpix):
    subpix = max(int(subpix), 1)
    return (np.arange(subpix, dtype=np.float64) + 0.5) / subpix - 0.5


def _flux_err(var, flux, gain):
    if gain is not None and float(gain) > 0:
        var = var + max(float(flux), 0.0) / float(gain)
    if var > 0:
        return math.sqrt(var)
    return 0.0


def _sum_one(data, bad, x, y, inside_of, err, gain, subpix, radius_pad):
    """``inside_of(dx, dy)`` is true for a sample inside the aperture."""
    ny, nx = data.shape
    pad = float(radius_pad)
    x0 = int(math.floor(x - pad))
    x1 = int(math.ceil(x + pad))
    y0 = int(math.floor(y - pad))
    y1 = int(math.ceil(y + pad))
    flag = 0
    if x0 < 0 or y0 < 0 or x1 >= nx or y1 >= ny:
        flag |= 2
    xa0, xa1 = max(0, x0), min(nx - 1, x1)
    ya0, ya1 = max(0, y0), min(ny - 1, y1)
    if xa1 < xa0 or ya1 < ya0:
        return 0.0, 0.0, flag
    yy, xx = np.mgrid[ya0:ya1 + 1, xa0:xa1 + 1]
    offs = _sample_offsets(subpix)
    pix = data[ya0:ya1 + 1, xa0:xa1 + 1]
    bad_box = bad[ya0:ya1 + 1, xa0:xa1 + 1]
    # Mean of the sub-pixel indicators. A bad pixel contributes nothing.
    acc = np.zeros(pix.shape, dtype=np.float64)
    n_sub = 0
    for oy in offs:
        for ox in offs:
            n_sub += 1
            acc += inside_of(xx + ox - x, yy + oy - y)
    frac = acc / n_sub
    if np.any((frac > 0) & bad_box):
        flag |= 1
    weight = np.where(bad_box, 0.0, frac)
    values = np.where(np.isfinite(pix), pix, 0.0)
    flux = float(np.sum(values * weight))
    var = 0.0
    if err is not None:
        err_arr = np.asarray(err, dtype=np.float64)
        if err_arr.ndim == 0:
            var = float(err_arr) ** 2 * float(np.sum(weight ** 2))
        else:
            em = err_arr[ya0:ya1 + 1, xa0:xa1 + 1]
            em = np.where(np.isfinite(em), em, 0.0)
            var = float(np.sum((em * weight) ** 2))
    return flux, _flux_err(var, flux, gain), flag


def sum_circle(data, x, y, r, err=None, var=None, gain=None, mask=None,
               maskthresh=0.0, subpix=5):
    data = np.asarray(data, dtype=np.float64)
    if err is None and var is not None:
        err = np.sqrt(np.asarray(var, dtype=np.float64))
    n, (xs, ys, rs) = _broadcast((x, y, r))
    if n == 0:
        return _empty("sum")
    bad = _bad_map(data, mask, maskthresh)
    flux = np.empty(n)
    fluxerr = np.empty(n)
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        radius = float(rs[i])
        if not np.isfinite(radius) or radius <= 0:
            flux[i], fluxerr[i], flag[i] = 0.0, 0.0, 0
            continue

        def inside(dx, dy, radius=radius):
            return (dx * dx + dy * dy) <= radius * radius

        flux[i], fluxerr[i], flag[i] = _sum_one(
            data, bad, float(xs[i]), float(ys[i]), inside, err, gain, subpix, radius)
    return flux, fluxerr, flag


def sum_circann(data, x, y, rin, rout, err=None, var=None, gain=None, mask=None,
                maskthresh=0.0, subpix=5):
    data = np.asarray(data, dtype=np.float64)
    if err is None and var is not None:
        err = np.sqrt(np.asarray(var, dtype=np.float64))
    n, (xs, ys, rins, routs) = _broadcast((x, y, rin, rout))
    if n == 0:
        return _empty("sum")
    bad = _bad_map(data, mask, maskthresh)
    flux = np.empty(n)
    fluxerr = np.empty(n)
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        r0 = max(float(rins[i]), 0.0)
        r1 = float(routs[i])
        if not np.isfinite(r1) or r1 <= r0:
            flux[i], fluxerr[i], flag[i] = 0.0, 0.0, 0
            continue

        def inside(dx, dy, r0=r0, r1=r1):
            d2 = dx * dx + dy * dy
            return (d2 <= r1 * r1) & (d2 >= r0 * r0)

        flux[i], fluxerr[i], flag[i] = _sum_one(
            data, bad, float(xs[i]), float(ys[i]), inside, err, gain, subpix, r1)
    return flux, fluxerr, flag


def sum_ellipse(data, x, y, a, b, theta, r, err=None, var=None, gain=None,
                mask=None, maskthresh=0.0, subpix=5):
    """Elliptical sum. The 7th positional argument ``r`` scales the semi-axes ``a``, ``b``."""
    data = np.asarray(data, dtype=np.float64)
    if err is None and var is not None:
        err = np.sqrt(np.asarray(var, dtype=np.float64))
    n, cols = _broadcast((x, y, a, b, theta, r))
    if n == 0:
        return _empty("sum")
    xs, ys, aa, bb, th, rr = cols
    bad = _bad_map(data, mask, maskthresh)
    flux = np.empty(n)
    fluxerr = np.empty(n)
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        ai = max(float(aa[i]), 1e-6)
        bi = max(float(bb[i]), 1e-6)
        scale = float(rr[i])
        if not np.isfinite(scale) or scale <= 0:
            flux[i], fluxerr[i], flag[i] = 0.0, 0.0, 0
            continue
        ct = math.cos(float(th[i]))
        st = math.sin(float(th[i]))
        # Bounding radius of the scaled ellipse.
        pad = scale * max(ai, bi)

        def inside(dx, dy, ai=ai, bi=bi, scale=scale, ct=ct, st=st):
            xr = dx * ct + dy * st
            yr = -dx * st + dy * ct
            return (xr / ai) ** 2 + (yr / bi) ** 2 <= scale * scale

        flux[i], fluxerr[i], flag[i] = _sum_one(
            data, bad, float(xs[i]), float(ys[i]), inside, err, gain, subpix, pad)
    return flux, fluxerr, flag


def _ellipse_rho(dx, dy, a, b, theta):
    ct = math.cos(float(theta))
    st = math.sin(float(theta))
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    return np.sqrt((xr / a) ** 2 + (yr / b) ** 2)


def kron_radius(data, x, y, a, b, theta, r, mask=None, maskthresh=0.0):
    """First moment of the light inside elliptical radius ``r`` (units of a, b)."""
    data = np.asarray(data, dtype=np.float64)
    n, cols = _broadcast((x, y, a, b, theta, r))
    if n == 0:
        return _empty("kron")
    xs, ys, aa, bb, th, rr = cols
    bad = _bad_map(data, mask, maskthresh)
    ny, nx = data.shape
    out = np.empty(n)
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        ai = max(float(aa[i]), 1e-6)
        bi = max(float(bb[i]), 1e-6)
        limit = float(rr[i])
        xi, yi = float(xs[i]), float(ys[i])
        if not np.isfinite(limit) or limit <= 0:
            out[i], flag[i] = np.nan, 0
            continue
        pad = limit * max(ai, bi)
        x0 = int(math.floor(xi - pad))
        x1 = int(math.ceil(xi + pad))
        y0 = int(math.floor(yi - pad))
        y1 = int(math.ceil(yi + pad))
        fl = 0
        if x0 < 0 or y0 < 0 or x1 >= nx or y1 >= ny:
            fl |= 2
        xa0, xa1 = max(0, x0), min(nx - 1, x1)
        ya0, ya1 = max(0, y0), min(ny - 1, y1)
        if xa1 < xa0 or ya1 < ya0:
            out[i], flag[i] = np.nan, fl
            continue
        yy, xx = np.mgrid[ya0:ya1 + 1, xa0:xa1 + 1]
        pix = data[ya0:ya1 + 1, xa0:xa1 + 1]
        good = np.isfinite(pix) & ~bad[ya0:ya1 + 1, xa0:xa1 + 1]
        if np.any(~good):
            # Only flag when a pixel inside the ellipse is bad.
            rho_all = _ellipse_rho(xx - xi, yy - yi, ai, bi, th[i])
            if np.any((rho_all <= limit) & ~good):
                fl |= 1
        rho = _ellipse_rho(xx - xi, yy - yi, ai, bi, th[i])
        use = good & (rho <= limit)
        weight = np.where(use, np.clip(pix, 0, None), 0.0)
        total = float(weight.sum())
        if total <= 0:
            out[i] = np.nan
        else:
            out[i] = float(np.sum(rho * weight) / total)
        flag[i] = fl
    return out, flag


def flux_radius(data, x, y, rmax, frac, normflux=None, mask=None,
                maskthresh=0.0, subpix=5):
    """Circular radius enclosing ``frac`` of the flux inside ``rmax`` (or of ``normflux``)."""
    data = np.asarray(data, dtype=np.float64)
    frac_in = np.asarray(frac, dtype=np.float64)
    scalar = frac_in.ndim == 0
    fracs = np.atleast_1d(frac_in).ravel()
    k = fracs.size
    n, (xs, ys, rmaxs) = _broadcast((x, y, rmax))
    if normflux is None:
        norms = None
    else:
        norms = np.atleast_1d(np.asarray(normflux, dtype=np.float64)).ravel()
        if norms.size == 1 and n > 1:
            norms = np.full(n, float(norms[0]))
        elif norms.size != n:
            raise ValueError("normflux length does not match the number of positions")
    if n == 0:
        shape = (0,) if scalar else (0, k)
        return np.zeros(shape), np.zeros(0, dtype=np.int32)
    bad = _bad_map(data, mask, maskthresh)
    ny, nx = data.shape
    offs = _sample_offsets(subpix)
    area = 1.0 / offs.size ** 2
    out = np.empty((n, k))
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        xi, yi = float(xs[i]), float(ys[i])
        radius = float(rmaxs[i])
        fl = 0
        if not np.isfinite(radius) or radius <= 0:
            out[i] = np.nan
            flag[i] = 0
            continue
        x0 = int(math.floor(xi - radius))
        x1 = int(math.ceil(xi + radius))
        y0 = int(math.floor(yi - radius))
        y1 = int(math.ceil(yi + radius))
        if x0 < 0 or y0 < 0 or x1 >= nx or y1 >= ny:
            fl |= 2
        xa0, xa1 = max(0, x0), min(nx - 1, x1)
        ya0, ya1 = max(0, y0), min(ny - 1, y1)
        if xa1 < xa0 or ya1 < ya0:
            out[i] = np.nan
            flag[i] = fl
            continue
        yy, xx = np.mgrid[ya0:ya1 + 1, xa0:xa1 + 1]
        pix = data[ya0:ya1 + 1, xa0:xa1 + 1]
        good = np.isfinite(pix) & ~bad[ya0:ya1 + 1, xa0:xa1 + 1]
        samples_r = []
        samples_f = []
        hit_bad = False
        for oy in offs:
            for ox in offs:
                dx = xx + ox - xi
                dy = yy + oy - yi
                rad = np.sqrt(dx * dx + dy * dy)
                inside = rad <= radius
                if np.any(inside & ~good):
                    hit_bad = True
                use = inside & good
                if not use.any():
                    continue
                samples_r.append(rad[use])
                samples_f.append(np.clip(pix[use], 0, None) * area)
        if hit_bad:
            fl |= 1
        if not samples_r:
            out[i] = np.nan
            flag[i] = fl
            continue
        rr = np.concatenate(samples_r)
        ff = np.concatenate(samples_f)
        order = np.argsort(rr)
        rr = rr[order]
        csum = np.cumsum(ff[order])
        total = float(csum[-1]) if norms is None else float(norms[i])
        if not np.isfinite(total) or total <= 0:
            out[i] = np.nan
            flag[i] = fl
            continue
        for j, frac_j in enumerate(fracs):
            target = float(frac_j) * total
            if target <= 0:
                out[i, j] = 0.0
                continue
            pos = int(np.searchsorted(csum, target, side="left"))
            if pos >= rr.size:
                out[i, j] = radius
            elif pos == 0:
                out[i, j] = float(rr[0])
            else:
                f0, f1 = float(csum[pos - 1]), float(csum[pos])
                r0, r1 = float(rr[pos - 1]), float(rr[pos])
                span = f1 - f0
                t = 0.0 if span <= 0 else (target - f0) / span
                out[i, j] = r0 + t * (r1 - r0)
        flag[i] = fl
    if scalar:
        out = out[:, 0]
    return out, flag


def winpos(data, x, y, sig, mask=None, maskthresh=0.0, subpix=5):
    """Gaussian-windowed centroid. ``sig`` is the window sigma in pixels. A few iterations."""
    data = np.asarray(data, dtype=np.float64)
    n, (xs, ys, sigs) = _broadcast((x, y, sig))
    if n == 0:
        return _empty("win")
    bad = _bad_map(data, mask, maskthresh)
    ny, nx = data.shape
    offs = _sample_offsets(subpix)
    xo = np.empty(n)
    yo = np.empty(n)
    flag = np.empty(n, dtype=np.int32)
    for i in range(n):
        xi, yi = float(xs[i]), float(ys[i])
        sigma = max(float(sigs[i]), 0.3)
        fl = 0
        for _ in range(4):
            pad = 4.0 * sigma
            x0 = int(math.floor(xi - pad))
            x1 = int(math.ceil(xi + pad))
            y0 = int(math.floor(yi - pad))
            y1 = int(math.ceil(yi + pad))
            if x0 < 0 or y0 < 0 or x1 >= nx or y1 >= ny:
                fl |= 2
            xa0, xa1 = max(0, x0), min(nx - 1, x1)
            ya0, ya1 = max(0, y0), min(ny - 1, y1)
            if xa1 < xa0 or ya1 < ya0:
                fl |= 1
                break
            yy, xx = np.mgrid[ya0:ya1 + 1, xa0:xa1 + 1]
            pix = data[ya0:ya1 + 1, xa0:xa1 + 1]
            good = np.isfinite(pix) & ~bad[ya0:ya1 + 1, xa0:xa1 + 1]
            acc_w = np.zeros(pix.shape, dtype=np.float64)
            for oy in offs:
                for ox in offs:
                    dx = xx + ox - xi
                    dy = yy + oy - yi
                    acc_w += np.exp(-0.5 * (dx * dx + dy * dy) / (sigma * sigma))
            weight = acc_w / offs.size ** 2
            if np.any((weight > 1e-3) & ~good):
                fl |= 1
            w = np.where(good, np.clip(pix, 0, None) * weight, 0.0)
            total = float(w.sum())
            if total <= 0:
                fl |= 1
                break
            xi = float(np.sum(xx * w) / total)
            yi = float(np.sum(yy * w) / total)
        xo[i], yo[i], flag[i] = xi, yi, fl
    return xo, yo, flag
