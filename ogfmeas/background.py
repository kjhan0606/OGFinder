"""Mesh background. Bertin & Arnouts 1996: clipped mode per tile, median filter, bilinear zoom.

Original implementation. It follows that published estimator and the call shape used in this
tree. It is not a numerical copy of any existing library.
"""
import numpy as np
from scipy import ndimage as ndi
from scipy.interpolate import RegularGridInterpolator


def _odd(k):
    k = max(int(k), 1)
    if k % 2 == 0:
        k += 1
    return k


def _mode_rms(values):
    """Clipped mode and rms. Mode is 2.5*median - 1.5*mean when that lies between them."""
    v = np.asarray(values, dtype=np.float64)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return np.nan, np.nan
    for _ in range(5):
        if v.size < 4:
            break
        med = np.median(v)
        mad = np.median(np.abs(v - med))
        sig = 1.4826 * mad
        if not np.isfinite(sig) or sig <= 0:
            break
        keep = np.abs(v - med) <= 3.0 * sig
        if keep.all() or int(keep.sum()) < 4:
            v = v[keep] if int(keep.sum()) >= 4 else v
            break
        v = v[keep]
    if v.size == 0:
        return np.nan, np.nan
    mean = float(np.mean(v))
    med = float(np.median(v))
    rms = float(np.std(v))
    if not np.isfinite(rms) or rms <= 0:
        return med, 0.0
    candidate = 2.5 * med - 1.5 * mean
    lo, hi = (mean, med) if mean <= med else (med, mean)
    if lo <= candidate <= hi:
        mode = candidate
    else:
        mode = med
    return float(mode), rms


def _tiles(n, step):
    step = max(int(step), 1)
    starts = list(range(0, int(n), step))
    if not starts:
        starts = [0]
    return starts, step


def _fill(mesh):
    finite = np.isfinite(mesh)
    if not finite.any():
        return np.zeros_like(mesh)
    filled = mesh.copy()
    filled[~finite] = np.median(mesh[finite])
    return filled


def _zoom(mesh, cy, cx, shape):
    ny, nx = shape
    mesh = np.asarray(mesh, dtype=np.float64)
    nby, nbx = mesh.shape
    if nby == 1 and nbx == 1:
        return np.full((ny, nx), float(mesh[0, 0]), dtype=np.float64)
    if nby == 1:
        row = mesh[0]
        x = np.arange(nx, dtype=np.float64)
        return np.tile(np.interp(x, cx, row), (ny, 1))
    if nbx == 1:
        col = mesh[:, 0]
        y = np.arange(ny, dtype=np.float64)
        return np.tile(np.interp(y, cy, col), (nx, 1)).T
    interp = RegularGridInterpolator(
        (np.asarray(cy, dtype=np.float64), np.asarray(cx, dtype=np.float64)),
        mesh, method="linear", bounds_error=False, fill_value=None)
    yy, xx = np.mgrid[0:ny, 0:nx]
    pts = np.column_stack([yy.ravel(), xx.ravel()])
    return interp(pts).reshape(ny, nx)


class Background:
    """Spatially varying background and rms.

    Mask values above ``maskthresh`` are ignored (a boolean mask uses True as bad
    when ``maskthresh`` is 0). ``back()`` is the bilinear full-resolution map.
    ``np.array(self)`` returns that map. ``globalback`` and ``globalrms`` are the
    means of the filtered mesh.
    """

    def __init__(self, data, mask=None, maskthresh=0.0, bw=64, bh=64,
                 fw=3, fh=3, fthresh=0.0):
        src = np.asarray(data)
        if src.ndim != 2:
            raise ValueError("background data must be 2-D")
        self._shape = src.shape
        self._dtype = src.dtype if np.issubdtype(src.dtype, np.floating) else np.float64
        image = np.asarray(src, dtype=np.float64)
        bad = ~np.isfinite(image)
        if mask is not None:
            bad = bad | (np.asarray(mask) > maskthresh)
        ny, nx = image.shape
        y_starts, bh = _tiles(ny, bh)
        x_starts, bw = _tiles(nx, bw)
        nby, nbx = len(y_starts), len(x_starts)
        back_m = np.empty((nby, nbx), dtype=np.float64)
        rms_m = np.empty((nby, nbx), dtype=np.float64)
        cy = np.empty(nby, dtype=np.float64)
        cx = np.empty(nbx, dtype=np.float64)
        for iy, y0 in enumerate(y_starts):
            y1 = min(ny, y0 + bh)
            cy[iy] = 0.5 * (y0 + y1 - 1)
            rows = image[y0:y1]
            bad_rows = bad[y0:y1]
            for ix, x0 in enumerate(x_starts):
                x1 = min(nx, x0 + bw)
                if iy == 0:
                    cx[ix] = 0.5 * (x0 + x1 - 1)
                tile = rows[:, x0:x1]
                keep = ~bad_rows[:, x0:x1]
                back_m[iy, ix], rms_m[iy, ix] = _mode_rms(tile[keep])
        back_m = _fill(back_m)
        rms_m = _fill(rms_m)
        fh_f, fw_f = _odd(fh), _odd(fw)
        if fh_f > 1 or fw_f > 1:
            filt_b = ndi.median_filter(back_m, size=(fh_f, fw_f), mode="nearest")
            filt_r = ndi.median_filter(rms_m, size=(fh_f, fw_f), mode="nearest")
            # fthresh 0 replaces every tile that differs, which is the usual default.
            if fthresh <= 0:
                back_m, rms_m = filt_b, filt_r
            else:
                use = np.abs(back_m - filt_b) > float(fthresh)
                back_m = np.where(use, filt_b, back_m)
                rms_m = np.where(use, filt_r, rms_m)
        self._mesh_back = back_m
        self._mesh_rms = rms_m
        self._cy = cy
        self._cx = cx
        self._back = None
        self._rms = None
        self.globalback = float(np.mean(back_m))
        self.globalrms = float(np.mean(rms_m))

    def back(self, dtype=None):
        if self._back is None:
            self._back = _zoom(self._mesh_back, self._cy, self._cx, self._shape)
        out = self._back if dtype is None else self._back.astype(dtype, copy=False)
        if dtype is None and out.dtype != self._dtype:
            out = out.astype(self._dtype, copy=False)
        return out

    def rms(self, dtype=None):
        if self._rms is None:
            self._rms = _zoom(self._mesh_rms, self._cy, self._cx, self._shape)
        out = self._rms if dtype is None else self._rms.astype(dtype, copy=False)
        if dtype is None and out.dtype != self._dtype:
            out = out.astype(self._dtype, copy=False)
        return out

    def subfrom(self, data):
        """Subtract the background from ``data`` in place."""
        np.subtract(data, self.back(dtype=getattr(data, "dtype", None)), out=data)

    def __array__(self, dtype=None, copy=None):
        arr = np.asarray(self.back(), dtype=dtype)
        if copy:
            return arr.copy()
        return arr
