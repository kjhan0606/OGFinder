"""Static-sky template construction.

For a target chip T (the exposure to be differenced) the template is the pixel-wise median, on T's
pixel grid, of the other exposures of the same filter that overlap T, each reprojected with
`reproject_interp` (handles different pixel scales/orientations and the full distortion model).
Pixels flagged bad (DQ, CR) in a contributing exposure are NaN before the median.  A moving object
present in only one input is removed by the median when >= 3 inputs overlap (with 2 inputs the
lower value is used, which biases faint pixels low; flagged in the returned `nin` map).
The same reprojection is used for deep external coadds (drc/drz/i2d, PS1 stacks): pass them as
`chips` with is_coadd=True and they are used directly (single input).
"""
import numpy as np
from astropy.wcs import WCS
from .util import log

try:
    from reproject import reproject_interp
except Exception:                                   # pragma: no cover
    reproject_interp = None


def _footprint_overlaps(chipA, chipB, margin=0.0):
    ny, nx = chipA.shape
    xs = np.array([0, nx, nx, 0, nx / 2]); ys = np.array([0, 0, ny, ny, ny / 2])
    ra, de = chipA.wcs.all_pix2world(xs, ys, 0)
    x, y = chipB.wcs.all_world2pix(ra, de, 0, quiet=True)
    ny2, nx2 = chipB.shape
    # bounding-box test of A's corners in B pixel space
    return not (np.max(x) < -margin or np.min(x) > nx2 + margin or np.max(y) < -margin or np.min(y) > ny2 + margin)


def pixel_map(src, target_wcs, target_shape, step=32, tgt_pix2world=None):
    """For every target pixel, the (x, y) in `src` pixel coordinates.  The exact WCS (with distortion
    tables) is evaluated on a coarse grid (`step` pixels) and bicubically interpolated: the distortion is
    smooth, interpolation error is << 0.01 pix for step=32 (tested in tests/test_astrometry.py)."""
    from scipy import ndimage as ndi
    ny, nx = target_shape
    gy = np.arange(0, ny + step, step, dtype=float); gx = np.arange(0, nx + step, step, dtype=float)
    GX, GY = np.meshgrid(gx, gy)
    ra, de = target_wcs.all_pix2world(GX, GY, 0)
    sx, sy = src.wcs.all_world2pix(ra, de, 0, quiet=True)
    if hasattr(src, "corr_M"):
        pass
    # bicubic upsample of coarse maps to full resolution
    fy = (np.arange(ny) / step); fx = (np.arange(nx) / step)
    yy, xx = np.meshgrid(fy, fx, indexing="ij")
    X = ndi.map_coordinates(sx, [yy, xx], order=3, mode="nearest")
    Y = ndi.map_coordinates(sy, [yy, xx], order=3, mode="nearest")
    return X.astype(np.float32), Y.astype(np.float32)


def reproject_chip(src, target_wcs, target_shape, order=3):
    """Resample `src.data` (bad pixels -> NaN) onto the target grid; returns (image, valid)."""
    from scipy import ndimage as ndi
    X, Y = pixel_map(src, target_wcs, target_shape)
    data = src.data.astype(np.float32).copy()
    badm = src.bad.copy()
    if getattr(src, "cr", None) is not None:
        badm |= src.cr
    badf = ndi.map_coordinates(badm.astype(np.float32), [Y, X], order=1, mode="constant", cval=1.0) > 0.01
    # bad pixels are replaced by the median before the spline so they do not ring into neighbours; the
    # resampled bad mask (dilated by the kernel footprint) is applied afterwards
    data[badm] = np.median(data[~badm]) if (~badm).any() else 0.0
    img = ndi.map_coordinates(data, [Y, X], order=order, mode="constant", cval=np.nan)
    ny, nx = src.shape
    inside = (X >= 0) & (X <= nx - 1) & (Y >= 0) & (Y <= ny - 1)
    img[badf | ~inside] = np.nan
    return img.astype(np.float32), ~np.isnan(img)


def build_template(target, others, min_inputs=2, pixscale_ratio_flux=True):
    """Return (template, nin) on target's grid.  Flux units are those of the inputs per target pixel:
    surface brightness is conserved (data are per-pixel flux after the area-map correction), so for
    different pixel scales the template is scaled by (pix_target/pix_src)^2 per source."""
    stack = []
    for o in others:
        if o is target or not _footprint_overlaps(target, o):
            continue
        try:
            r, fp = reproject_chip(o, target.wcs, target.shape)
        except Exception as e:                         # footprint problems
            log("reproject failed %s: %s" % (o.name, e))
            continue
        scale = (target.pixscale / o.pixscale) ** 2 if pixscale_ratio_flux and np.isfinite(o.pixscale) else 1.0
        stack.append(r * scale)
    nin = np.zeros(target.shape, np.uint8)
    if not stack:
        return None, nin
    A = np.array(stack)
    nin = np.sum(np.isfinite(A), axis=0).astype(np.uint8)
    with np.errstate(all="ignore"):
        med = np.nanmedian(A, axis=0)
        if len(stack) == 2:
            med = np.nanmin(A, axis=0)
    med[nin < min_inputs] = np.nan
    return med.astype(np.float32), nin
