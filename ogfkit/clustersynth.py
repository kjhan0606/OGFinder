"""Synthetic clusters and lensed-arc images with known truth (validation of ogfkit.cluster; every function takes a seed)."""
import math

import numpy as np
from scipy import ndimage as ndi


def cluster_catalog(seed=1, shape=(4000, 4000), centre=(2200.0, 1800.0), n_field=3000, n_mem=150, rcore=150.0, rmax=1500.0, zc=0.4,
                    rs_a=1.20, rs_b=-0.04, rs_m0=21.0, rs_scatter=0.05, mag_lim=(19.0, 25.0), lf_slope=0.35):
    """Field + cluster catalog.  Returns dict of arrays: x, y (0-based), mag, color, color_err, z, z_err, member (bool), color_true.

    Magnitudes follow dN/dm ~ 10^(lf_slope m) for both populations; cluster members sit on a red sequence colour = rs_a + rs_b (m - rs_m0)
    with intrinsic scatter rs_scatter and a Plummer profile; field galaxies are uniform with a broad colour distribution (40 % blue, 40 % red,
    20 % uniform) that overlaps the sequence; photometric errors grow with magnitude."""
    rng = np.random.default_rng(seed)
    a = lf_slope * math.log(10)
    lo, hi = mag_lim

    def mags(n):
        u = rng.random(n)
        return lo + np.log(1 + u * (math.exp(a * (hi - lo)) - 1)) / a

    mf = mags(n_field)
    xf = rng.uniform(0, shape[1], n_field); yf = rng.uniform(0, shape[0], n_field)
    kind = rng.random(n_field)
    cf = np.where(kind < 0.4, rng.normal(0.5, 0.25, n_field), np.where(kind < 0.8, rng.normal(1.0, 0.30, n_field), rng.uniform(-0.2, 2.0, n_field)))
    zf = rng.uniform(0.05, 1.6, n_field)
    xm, ym = [], []
    while len(xm) < n_mem:
        u = rng.random()
        r = rcore * math.sqrt(u / (1 - u)) if u < 0.999 else rmax
        if r > rmax:
            continue
        ph = rng.uniform(0, 2 * math.pi)
        px, py = centre[0] + r * math.cos(ph), centre[1] + r * math.sin(ph)
        if 0 <= px < shape[1] and 0 <= py < shape[0]:
            xm.append(px); ym.append(py)
    xm, ym = np.array(xm), np.array(ym)
    mm = mags(n_mem)
    cm = rs_a + rs_b * (mm - rs_m0) + rng.normal(0, rs_scatter, n_mem)
    zm = rng.normal(zc, 0.003 * (1 + zc), n_mem)
    x = np.concatenate([xf, xm]); y = np.concatenate([yf, ym]); m = np.concatenate([mf, mm]); ct = np.concatenate([cf, cm]); zt = np.concatenate([zf, zm])
    err = 0.01 + 0.012 * 10 ** (0.4 * (m - 22.0))
    err = np.minimum(err, 0.3)
    c = ct + rng.normal(0, 1, len(x)) * err
    zerr = 0.03 * (1 + zt)
    z = zt + rng.normal(0, 1, len(x)) * zerr
    return dict(x=x, y=y, mag=m, color=c, color_err=err, z=z, z_err=zerr, member=np.concatenate([np.zeros(n_field, bool), np.ones(n_mem, bool)]), color_true=ct,
                centre=centre, shape=shape, truth=dict(a=rs_a, b=rs_b, m0=rs_m0, scatter=rs_scatter, zc=zc))


def _taper(u, half, soft=0.6):
    return 0.5 * (1 + np.tanh((half - np.abs(u)) / soft))


def _render(shape, x0, y0, half_box, fn, os_=3):
    """Evaluate fn(X, Y) (image pixel coordinates) on an os_-times oversampled grid inside a box and bin down."""
    xa, xb = int(max(0, math.floor(x0 - half_box))), int(min(shape[1], math.ceil(x0 + half_box)))
    ya, yb = int(max(0, math.floor(y0 - half_box))), int(min(shape[0], math.ceil(y0 + half_box)))
    if xb <= xa or yb <= ya:
        return None
    gx = xa - 0.5 + (np.arange((xb - xa) * os_) + 0.5) / os_
    gy = ya - 0.5 + (np.arange((yb - ya) * os_) + 0.5) / os_
    X, Y = np.meshgrid(gx, gy)
    v = fn(X, Y)
    v = v.reshape(yb - ya, os_, xb - xa, os_).mean(axis=(1, 3))
    return xa, ya, v


def add_arc(img, cx, cy, R, phi0, span, fwhm_w, amp):
    """Circular arc of radius R centred on (cx, cy), centre angle phi0 (rad), angular span `span` (rad), Gaussian cross-section FWHM fwhm_w,
    peak surface brightness amp (before the PSF).  Returns (x, y) of the arc's centre-point."""
    sw = fwhm_w / 2.3548
    def fn(X, Y):
        rho = np.hypot(X - cx, Y - cy)
        dphi = np.angle(np.exp(1j * (np.arctan2(Y - cy, X - cx) - phi0)))
        return amp * np.exp(-0.5 * ((rho - R) / sw) ** 2) * _taper(dphi * R, 0.5 * span * R)
    mx, my = cx + R * math.cos(phi0), cy + R * math.sin(phi0)
    box = R + 4 * sw + 3 if R * span < 400 else R + 4 * sw + 3
    # restrict the box to the arc's neighbourhood
    half = min(R * 2 * math.sin(min(span / 2, math.pi / 2)) / 2 + 4 * sw + 4, R + 4 * sw + 3)
    r = _render(img.shape, mx, my, half + 2, fn)
    if r:
        xa, ya, v = r
        img[ya:ya + v.shape[0], xa:xa + v.shape[1]] += v
    return mx, my


def add_line(img, x0, y0, length, theta, fwhm_w, amp):
    """Straight edge-on-like feature (uniform ridge of the given length, Gaussian cross-section)."""
    sw = fwhm_w / 2.3548
    ct, st = math.cos(theta), math.sin(theta)
    def fn(X, Y):
        u = (X - x0) * ct + (Y - y0) * st; t = -(X - x0) * st + (Y - y0) * ct
        return amp * np.exp(-0.5 * (t / sw) ** 2) * _taper(u, 0.5 * length)
    r = _render(img.shape, x0, y0, 0.5 * length + 4 * sw + 4, fn)
    if r:
        xa, ya, v = r
        img[ya:ya + v.shape[0], xa:xa + v.shape[1]] += v


def add_blob(img, x0, y0, sigma, q, theta, amp):
    ct, st = math.cos(theta), math.sin(theta)
    def fn(X, Y):
        u = (X - x0) * ct + (Y - y0) * st; t = -(X - x0) * st + (Y - y0) * ct
        return amp * np.exp(-0.5 * ((u / sigma) ** 2 + (t / (sigma * q)) ** 2))
    r = _render(img.shape, x0, y0, 5 * sigma + 3, fn)
    if r:
        xa, ya, v = r
        img[ya:ya + v.shape[0], xa:xa + v.shape[1]] += v


def truth_ridge(t, centre, n=60):
    """Sample points (x, y) along the centre-line of a truth object of arc_field."""
    if t['kind'] == 'arc':
        phi0 = math.atan2(t['y'] - centre[1], t['x'] - centre[0])
        ph = phi0 + np.linspace(-0.5, 0.5, n) * math.radians(t['span_deg'])
        return centre[0] + t['R'] * np.cos(ph), centre[1] + t['R'] * np.sin(ph)
    u = np.linspace(-0.5, 0.5, n) * (t['length'] if t['kind'] != 'blob' else 0.0)
    return t['x'] + u * math.cos(t['angle']), t['y'] + u * math.sin(t['angle'])


def arc_field(seed=1, shape=(700, 700), centre=(350.0, 350.0), n_arcs=10, n_radial=6, n_straight_tang=6, n_blobs=30, fwhm_psf=3.0, noise=1.0,
              amp_range=(3.0, 8.0), width_range=(3.0, 6.0), r_range=(60.0, 220.0), span_range_deg=(35.0, 110.0), image=True, min_sep=30.0):
    """Image with tangential arcs about `centre`, radially oriented lines, tangentially oriented straight lines and elliptical blobs.

    amp (peak surface brightness before the PSF) is in units of the pixel noise.  Returns (img float64, truth list of dicts with kind, x, y,
    R, span_deg, length, width (FWHM incl. PSF), amp)."""
    rng = np.random.default_rng(seed)
    clean = np.zeros(shape)
    truth = []
    ridges = []
    def free(x, y, d=25):
        return all(math.hypot(x - t['x'], y - t['y']) > d for t in truth) and 20 < x < shape[1] - 20 and 20 < y < shape[0] - 20
    def clear(rx, ry, d=min_sep):
        return all(np.min(np.hypot(rx[:, None] - ox[None, :], ry[:, None] - oy[None, :])) > d for ox, oy in ridges)
    for kind, n in (('arc', n_arcs), ('radial', n_radial), ('straight_tangential', n_straight_tang), ('blob', n_blobs)):
        for _ in range(n):
            for _try in range(200):
                R = rng.uniform(*r_range); phi = rng.uniform(0, 2 * math.pi)
                wid = rng.uniform(*width_range); amp = rng.uniform(*amp_range) * noise
                span = math.radians(rng.uniform(*span_range_deg))
                if kind == 'arc':
                    x, y = centre[0] + R * math.cos(phi), centre[1] + R * math.sin(phi)
                    if not free(x, y, 40) or not (40 < x < shape[1] - 40 and 40 < y < shape[0] - 40):
                        continue
                    # the whole arc must stay inside the image
                    ends = [(centre[0] + R * math.cos(phi + s * span / 2), centre[1] + R * math.sin(phi + s * span / 2)) for s in (-1, 1)]
                    if any(not (10 < e[0] < shape[1] - 10 and 10 < e[1] < shape[0] - 10) for e in ends):
                        continue
                    t = dict(kind=kind, x=x, y=y, R=R, span_deg=math.degrees(span), length=R * span, width=math.hypot(wid, fwhm_psf), amp=amp, angle=(phi + math.pi / 2) % math.pi)
                    rx, ry = truth_ridge(t, centre)
                    if not clear(rx, ry, min_sep + 2 * wid):
                        continue
                    add_arc(clean, centre[0], centre[1], R, phi, span, wid, amp)
                    ridges.append((rx, ry)); truth.append(t)
                    break
                else:
                    x, y = rng.uniform(40, shape[1] - 40), rng.uniform(40, shape[0] - 40)
                    if not free(x, y, 40) or math.hypot(x - centre[0], y - centre[1]) < 30:
                        continue
                    rad = math.atan2(y - centre[1], x - centre[0])
                    length = rng.uniform(25, 70)
                    tmp = dict(kind=kind, x=x, y=y, length=length if kind != 'blob' else 0.0, angle=rad % math.pi if kind == 'radial' else (rad + math.pi / 2) % math.pi)
                    rx, ry = truth_ridge(tmp, centre)
                    if not clear(rx, ry, min_sep + 8):
                        continue
                    if kind == 'radial':
                        add_line(clean, x, y, length, rad, wid, amp); ang = rad % math.pi
                    elif kind == 'straight_tangential':
                        add_line(clean, x, y, length, rad + math.pi / 2, wid, amp); ang = (rad + math.pi / 2) % math.pi
                    else:
                        sg = rng.uniform(2.5, 5.0); q = rng.uniform(0.5, 1.0); ang = rng.uniform(0, math.pi)
                        add_blob(clean, x, y, sg, q, ang, amp * 2)
                        length, wid = 2.3548 * sg, 2.3548 * sg * q
                    t = dict(kind=kind, x=x, y=y, R=float('inf'), span_deg=0.0, length=length, width=math.hypot(wid, fwhm_psf) if kind != 'blob' else wid, amp=amp, angle=ang)
                    ridges.append(truth_ridge(t, centre) if kind != 'blob' else (np.array([x]), np.array([y])))
                    truth.append(t)
                    break
    sig = fwhm_psf / 2.3548
    img = ndi.gaussian_filter(clean, sig) if sig > 0 else clean
    if image:
        img = img + rng.normal(0, noise, shape)
    return img, truth
