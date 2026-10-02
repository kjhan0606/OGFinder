"""Image I/O and basic processing for calibrated HST/JWST/PS1/generic FITS exposures.

A `Chip` is one science array with its WCS (including distortion), error and bad-pixel masks, in
units of electrons/s (HST) or the native unit (others) and corrected for the geometric pixel-area
variation (so that the surface brightness is constant across the detector).
"""
import os, json
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.time import Time
from scipy import ndimage as ndi
from .util import log

# DQ bits considered bad (HST ACS/WFC3): 1 bad pixel, 2 data lost, 4 hot? ... CR=4096, WARM=16(not bad), HOT=16
HST_BAD_DQ = 1 | 2 | 4 | 8 | 32 | 256 | 512 | 1024 | 8192 | 16384
CR_DQ = 4096                       # cosmic ray flagged by the archive pipeline (kept separate: fast movers can be mis-flagged)
SATURATED = 256


class Chip:
    def crop(self, x0, x1, y0, y1):
        """New Chip restricted to pixels [x0:x1, y0:y1] (WCS sliced, distortion tables preserved)."""
        import copy
        c = copy.copy(self)
        sl = (slice(y0, y1), slice(x0, x1))
        c.wcs = self.wcs[sl]
        for k in ("data", "err", "bad", "saturated", "cr", "pam"):
            v = getattr(self, k, None)
            if v is not None:
                setattr(c, k, v[sl].copy())
        c.shape = c.data.shape
        c.origin = (x0, y0)
        c.name = self.name + "[%d:%d,%d:%d]" % (x0, x1, y0, y1)
        return c


def mid_mjd(h0, h1):
    for k in ("EXPSTART", "EXPEND"):
        if k not in h1 and k not in h0:
            break
    else:
        es = h1.get("EXPSTART", h0.get("EXPSTART")); ee = h1.get("EXPEND", h0.get("EXPEND"))
        return 0.5 * (es + ee), es, ee
    t = Time(h0["DATE-OBS"] + "T" + h0["TIME-OBS"], scale="utc").mjd
    ex = float(h0.get("EXPTIME", 0))
    return t + 0.5 * ex / 86400.0, t, t + ex / 86400.0


def pixel_area_map(wcs, shape, step=64):
    """Relative pixel area on the sky (normalised to the median) from the WCS Jacobian,
    evaluated on a coarse grid and bilinearly interpolated.  Stand-in for the pixel area map file."""
    ny, nx = shape
    xs = np.linspace(0, nx - 1, max(2, nx // step + 1)); ys = np.linspace(0, ny - 1, max(2, ny // step + 1))
    X, Y = np.meshgrid(xs, ys)
    def sky(x, y):
        r, d = wcs.all_pix2world(x, y, 0)
        return r, d
    h = 0.5
    r1, d1 = sky(X + h, Y); r0, d0 = sky(X - h, Y)
    r3, d3 = sky(X, Y + h); r2, d2 = sky(X, Y - h)
    cd = np.cos(np.radians(d1))
    dxr = ((r1 - r0 + 180) % 360 - 180) * cd; dxd = d1 - d0
    dyr = ((r3 - r2 + 180) % 360 - 180) * cd; dyd = d3 - d2
    area = np.abs(dxr * dyd - dxd * dyr)
    area /= np.median(area)
    f = ndi.zoom(area, (ny / (area.shape[0] - 1 + 1e-9) * 0 + 1,) * 2, order=1)  # placeholder (replaced below)
    from scipy.interpolate import RegularGridInterpolator
    rgi = RegularGridInterpolator((ys, xs), area, bounds_error=False, fill_value=None)
    yy, xx = np.mgrid[0:ny, 0:nx]
    # evaluate in strips to bound memory
    out = np.empty((ny, nx), np.float32)
    for y0 in range(0, ny, 512):
        y1 = min(ny, y0 + 512)
        out[y0:y1] = rgi(np.stack([yy[y0:y1].ravel(), xx[y0:y1].ravel()], 1)).reshape(y1 - y0, nx)
    return out


def load_chips(path, align_dir=None, apply_pam=True, max_chips=None):
    """Load all science chips of a calibrated file.  Returns list of Chip."""
    chips = []
    hdul = fits.open(path, memmap=True)
    h0 = hdul[0].header
    inst = str(h0.get("INSTRUME", "")).strip()
    tel = str(h0.get("TELESCOP", "")).strip()
    sci = [i for i, h in enumerate(hdul) if h.name == "SCI" and h.data is not None]
    if not sci:                                   # plain image (PS1 cutout, drizzle product, ...)
        sci = [i for i, h in enumerate(hdul) if h.data is not None and h.data.ndim == 2]
    for n, i in enumerate(sci):
        if max_chips and n >= max_chips:
            break
        h = hdul[i]; hd = h.header
        c = Chip()
        c.path = path; c.ext = i; c.extver = hd.get("EXTVER", 1); c.name = "%s[%s,%s]" % (os.path.basename(path), h.name, c.extver)
        c.tel = tel; c.inst = inst
        try:
            c.wcs = WCS(hd, fobj=hdul)
        except Exception:
            c.wcs = WCS(hd)
        c.shape = h.data.shape
        data = np.array(h.data, dtype=np.float32)
        bunit = str(hd.get("BUNIT", h0.get("BUNIT", ""))).upper()
        texp = float(h0.get("EXPTIME", hd.get("EXPTIME", 1.0)) or 1.0)
        c.texp = texp
        c.filter = next((str(h0[k]).strip() for k in ("FILTER", "FILTER1", "FILTER2", "PUPIL") if k in h0
                         and str(h0[k]).strip() and not str(h0[k]).startswith(("CLEAR", "N/A"))), "")
        if tel == "HST" and "ELECTRON" in bunit and "/S" not in bunit:
            fac = 1.0 / texp
        else:
            fac = 1.0
        c.counts_per_s_factor = fac
        err = None; dq = None
        for j, hh in enumerate(hdul):
            if hh.name == "ERR" and hh.header.get("EXTVER", 1) == c.extver and hh.data is not None:
                err = np.array(hh.data, np.float32)
            if hh.name == "DQ" and hh.header.get("EXTVER", 1) == c.extver and hh.data is not None:
                dq = np.array(hh.data)
        data *= fac
        err = err * fac if err is not None else None
        bad = ~np.isfinite(data)
        if dq is not None:
            bad |= (dq.astype(np.int64) & HST_BAD_DQ) != 0
            c.saturated = (dq.astype(np.int64) & SATURATED) != 0
            c.cr = (dq.astype(np.int64) & CR_DQ) != 0
        else:
            c.saturated = np.zeros(data.shape, bool)
            c.cr = np.zeros(data.shape, bool)
        if apply_pam and tel == "HST":
            pam = pixel_area_map(c.wcs, data.shape)
            data /= pam
            if err is not None:
                err /= pam
            c.pam = pam
        c.data = np.where(bad, 0, data).astype(np.float32)
        c.err = err
        c.bad = bad
        c.mjd, c.mjd_start, c.mjd_end = mid_mjd(h0, hd)
        c.photflam = hd.get("PHOTFLAM", h0.get("PHOTFLAM"))
        c.photplam = hd.get("PHOTPLAM", h0.get("PHOTPLAM"))
        c.photzpt = hd.get("PHOTZPT", h0.get("PHOTZPT", -21.10))
        c.pixscale = float(np.sqrt(abs(np.linalg.det(c.wcs.pixel_scale_matrix)))) * 3600.0 if c.wcs.has_celestial else np.nan
        c.header0 = h0; c.header = hd
        c.zp_ab = None
        if c.photflam is not None and c.photplam is not None:
            c.zp_ab = -2.5 * np.log10(float(c.photflam)) - 5 * np.log10(float(c.photplam)) - 2.408
        if tel == "HST" and "ELECTRON" in bunit:
            c.gain = float(h0.get("CCDGAIN", 1.0) or 1.0)
        chips.append(c)
    return chips


def align_sidecar(path, align_dir):
    return os.path.join(align_dir, os.path.basename(path).replace(".fits", "") + ".align.json")


def background(img, mask=None, box=256):
    """Robust mesh background and rms using sep."""
    import sep
    d = np.ascontiguousarray(img.astype(np.float32))
    m = None if mask is None else np.ascontiguousarray(mask)
    bkg = sep.Background(d, mask=m, bw=box, bh=box, fw=3, fh=3)
    return bkg.back().astype(np.float32), bkg.rms().astype(np.float32)


def lacosmic_mask(img, sigclip=5.0, sigfrac=0.3, objlim=5.0, niter=3, gain=1.0, readnoise=4.0, bad=None):
    """Simplified L.A.Cosmic (van Dokkum 2001): Laplacian-edge detection on the 2x subsampled
    image, normalised by a noise model, with fine-structure rejection to protect stars.
    Operates on e-/s data scaled by `gain` internally.  Returns boolean CR mask."""
    d = np.array(img, np.float64)
    crmask = np.zeros(d.shape, bool)
    kern = np.array([[0, -1, 0], [-1, 4, -1], [0, -1, 0]], float)
    for _ in range(niter):
        work = d.copy()
        work[crmask] = np.nan
        med5 = ndi.median_filter(np.where(np.isfinite(work), work, np.nanmedian(work)), size=5)
        work = np.where(np.isfinite(work), work, med5)
        sub = np.kron(work, np.ones((2, 2)))
        lap = ndi.convolve(sub, kern, mode="mirror")
        lap[lap < 0] = 0
        lap = lap.reshape(d.shape[0], 2, d.shape[1], 2).mean(axis=(1, 3)) * 1.0  # block-average
        med5 = ndi.median_filter(work, size=5)
        noise = np.sqrt(np.clip(med5 * gain, 0, None) + readnoise ** 2) / gain
        S = lap / (2.0 * noise + 1e-12)
        S = S - ndi.median_filter(S, size=5)
        med3 = ndi.median_filter(work, size=3)
        med37 = ndi.median_filter(med3, size=7)
        F = (med3 - med37) / (noise + 1e-12)
        F = np.clip(F, 0.01, None)
        cand = (S > sigclip) & (lap / (F + 1e-12) / (2 * noise + 1e-12) > objlim)
        cand = (S > sigclip) & (S / F > objlim)
        # grow into neighbours with lower threshold
        grown = ndi.binary_dilation(cand, structure=np.ones((3, 3))) & (S > sigclip * sigfrac)
        new = cand | grown
        if bad is not None:
            new &= ~bad
        if not (new & ~crmask).any():
            break
        crmask |= new
    return crmask


def psf_star_list(chip, snr_min=10.0):
    """Isolated, unsaturated, compact stars of one chip as (x, y) using the half-light-radius locus."""
    from . import align as A
    st = A.detect_stars(chip, snr_min=snr_min)
    return st


def collect_star_cutouts(chips, size=25, nmax=200, snr_min=10.0, confirm=True):
    """Sub-pixel re-centred, flux-normalised cutouts of isolated compact unsaturated stars confirmed in >= 1 other exposure.
    Returns (stack list, fwhm list, x list, y list, chip-name list); x, y are the star positions in the pixel frame of the chip
    they were measured on (full-chip frame; `Chip.origin` is NOT added).  `confirm=False` skips the other-exposure confirmation (for a template that is
    already a median of several exposures, where cosmic rays are gone)."""
    from scipy.spatial import cKDTree
    chips = [chips] if isinstance(chips, Chip) else list(chips)
    r = size // 2
    stack = []; fw = []; sx = []; sy = []; sc = []
    # cosmic rays are random: keep only stars whose sky position is also found (compact source) in >= 1 other exposure
    allst = {}
    for c in chips:
        st = psf_star_list(c, snr_min)
        if st is None:
            continue
        ra_, de_ = c.wcs.all_pix2world(st["x"], st["y"], 0)
        allst[c.name] = (c, st, np.c_[np.cos(np.radians(de_)) * ra_, de_])
    trees = {k: cKDTree(v[2]) for k, v in allst.items()}
    for c in chips:
        if c.name not in allst:
            continue
        _, st, xy = allst[c.name]
        nconf = np.zeros(len(xy), int)
        for k2, (c2, st2, xy2) in allst.items():
            if c2.path == c.path:
                continue
            d_, _j = trees[k2].query(xy, distance_upper_bound=1.0 / 3600.0)
            nconf += np.isfinite(d_)
        keepm = (nconf >= 1) if confirm else np.ones(len(nconf), bool)
        if keepm.sum() < 1:
            continue
        st = dict(st); 
        for kk in ("x", "y", "flux", "snr"):
            st[kk] = st[kk][keepm]
        bkg, rms = background(c.data, c.bad)
        sub = (c.data - bkg).astype(np.float32)
        ny, nx = c.shape
        # all significant sources (for crowding test)
        import sep
        o = sep.extract(np.ascontiguousarray(sub), 4.0, err=np.ascontiguousarray(rms), mask=np.ascontiguousarray(c.bad), minarea=4)
        tree = cKDTree(np.c_[o["x"], o["y"]]) if len(o) else None
        for x0, y0 in zip(st["x"], st["y"]):
            xi, yi = int(round(x0)), int(round(y0))
            if xi < r + 3 or yi < r + 3 or xi > nx - r - 4 or yi > ny - r - 4:
                continue
            if tree is not None:
                nb = [j for j in tree.query_ball_point([x0, y0], 9.0)]
                # crowding: any other source within 9 px with > 5% of this star's flux
                f0 = float(st["flux"][list(st["x"]).index(x0)]) if "flux" in st else 1.0
                if sum(1 for j in nb if abs(o["x"][j] - x0) + abs(o["y"][j] - y0) > 1.0 and o["flux"][j] > 0.05 * f0) > 0:
                    continue
            if c.bad[yi - 6:yi + 7, xi - 6:xi + 7].any() or c.saturated[yi - 6:yi + 7, xi - 6:xi + 7].any():
                continue
            # the archive CR flag is dense in these data and also marks real stellar cores: not used here; instead
            # cutouts are median-combined (robust to CR hits in individual cutouts)
            cut = sub[yi - r - 1:yi + r + 2, xi - r - 1:xi + r + 2].astype(np.float64)
            sh = ndi.shift(cut, (yi - y0, xi - x0), order=3, mode="nearest")[1:-1, 1:-1]
            ssum = sh.sum()
            if ssum <= 0:
                continue
            stack.append(sh / ssum)
            fw.append(st["fwhm"]); sx.append(float(x0)); sy.append(float(y0)); sc.append(c.name)
            if len(stack) >= nmax:
                break
    return stack, fw, sx, sy, sc


def _combine_cutouts(stack, size):
    """Median of normalised cutouts, apodised at the edge, unit sum.  Returns (psf, fwhm_pix)."""
    r = size // 2
    psf = np.median(np.array(stack), axis=0)
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    rr = np.hypot(xx, yy)
    psf = np.clip(psf, 0, None) * np.clip((r - rr) / 3.0, 0, 1)
    psf /= psf.sum()
    hm = psf > psf.max() / 2
    return psf, float(2.0 * np.sqrt(hm.sum() / np.pi))


def estimate_psf(chips, size=25, nmax=200, snr_min=10.0):
    """Empirical PSF pooled over one or several chips (same detector/filter): median of sub-pixel
    re-centred, flux-normalised cutouts of isolated compact unsaturated stars.  Returns
    (psf normalised to 1, fwhm_pix, n_stars) or (None, nan, n) if fewer than 3 stars."""
    stack = collect_star_cutouts(chips, size, nmax, snr_min)[0]
    if len(stack) < 3:
        return None, np.nan, len(stack)
    psf, fw = _combine_cutouts(stack, size)
    return psf, fw, len(stack)


class PSFField:
    """Spatially varying PSF on a grid of tiles of one chip.  `at(x, y)` returns the PSF (unit sum, odd size) of the tile containing the pixel
    position; tiles with fewer than `min_stars` stars (own tile + `reach` tile widths around it) use the constant fallback PSF.
    Attributes: `tiles` (ny_t, nx_t) tile size, `n_stars` (grid of counts), `from_stars` (bool grid), `fwhm` (grid, pixels), `constant` (fallback)."""

    def __init__(self, shape, tile, psfs, n_stars, from_stars, fwhm, constant, constant_fwhm):
        self.shape = shape; self.tile = tile; self.psfs = psfs; self.n_stars = n_stars; self.from_stars = from_stars
        self.fwhm = fwhm; self.constant = constant; self.constant_fwhm = constant_fwhm

    def at(self, x, y):
        j = int(min(max(x // self.tile[1], 0), self.psfs.shape[1] - 1)); i = int(min(max(y // self.tile[0], 0), self.psfs.shape[0] - 1))
        return self.psfs[i, j]

    @property
    def varies(self):
        return bool(self.from_stars.any())

    def summary(self):
        return dict(tile=list(self.tile), n_tiles=int(self.psfs.size), tiles_from_stars=int(self.from_stars.sum()), stars=int(self.n_stars.sum()),
                    fwhm_pix_min=float(np.nanmin(self.fwhm)), fwhm_pix_max=float(np.nanmax(self.fwhm)), constant_fwhm_pix=float(self.constant_fwhm))


def psf_model_field(chip, psf, fwhm, fwhm_floor=0.0, snr_min=15.0):
    """Spatially varying PSF of one chip as a `ogfkit.psfmodel.PSFModel` (has `at(x, y)` and `summary()` like `PSFField`).  The chip's bad / saturated /
    cosmic-ray pixels are masked.  A model narrower than `fwhm_floor` (hot-pixel contamination) or built from < 4 stars is replaced by the
    constant fallback PSF `psf` (so the tiled ZOGY never sees a worse PSF than the default one)."""
    import os, sys
    root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    if root not in sys.path:
        sys.path.insert(0, root)
    from ogfkit import psfmodel as PM
    mask = np.zeros(chip.data.shape, bool)
    for k in ("bad", "saturated", "cr"):
        v = getattr(chip, k, None)
        if v is not None:
            mask |= np.asarray(v).astype(bool)
    mdl, info = PM.build_psf_model(np.asarray(chip.data, np.float32), mask=mask, fwhm_prior=max(float(fwhm), 1.5), snr_min=snr_min, n_iter=2)
    if info.get("mode") not in ("empirical",) or mdl.fwhm_estimate() < fwhm_floor:
        return PM.from_image(psf)
    return mdl


def constant_psf_field(shape, psf, fwhm, tile=(512, 512)):
    ny, nx = shape; nty = int(np.ceil(ny / tile[0])); ntx = int(np.ceil(nx / tile[1]))
    arr = np.empty((nty, ntx), object)
    for i in range(nty):
        for j in range(ntx):
            arr[i, j] = psf
    return PSFField(shape, tile, arr, np.zeros((nty, ntx), int), np.zeros((nty, ntx), bool), np.full((nty, ntx), float(fwhm)), psf, float(fwhm))


def measure_psf_field(chips, shape, tile=(512, 512), size=25, min_stars=8, reach=0.5, snr_min=6.0, constant=None, origin=(0, 0), confirm=True,
                      min_fwhm_pix=0.0):
    """Measure the PSF per tile from field stars (`collect_star_cutouts`; stars of all `chips` given, which must share the detector pixel frame,
    e.g. the chip and the other exposures' same chip) and return a `PSFField` for an image of `shape` (ny, nx) whose pixel (0, 0) is chip pixel
    `origin` = (x0, y0).  A tile uses its own stars plus those within `reach` tile widths beyond its border; with fewer than `min_stars` it
    falls back to `constant` = (psf, fwhm) (default: pooled over all stars, else a 2.2 px Gaussian).  A tile PSF narrower than `min_fwhm_pix` (cosmic-ray /
    hot-pixel contamination of the star list) is rejected and the tile uses the fallback as well."""
    stack, fw, sx, sy, sc = collect_star_cutouts(chips, size, 10 ** 6, snr_min, confirm)
    ny, nx = shape
    nty = int(np.ceil(ny / tile[0])); ntx = int(np.ceil(nx / tile[1]))
    sx = np.asarray(sx) - origin[0]; sy = np.asarray(sy) - origin[1]
    if constant is None and len(stack) >= 3:
        constant = _combine_cutouts(stack, size)
    if constant is None:
        constant = (gaussian_psf(2.2, size), 2.2)
    arr = np.empty((nty, ntx), object); ns = np.zeros((nty, ntx), int); fs = np.zeros((nty, ntx), bool); fwg = np.zeros((nty, ntx))
    for i in range(nty):
        for j in range(ntx):
            y0, y1, x0, x1 = i * tile[0], (i + 1) * tile[0], j * tile[1], (j + 1) * tile[1]
            m = (sx >= x0 - reach * tile[1]) & (sx < x1 + reach * tile[1]) & (sy >= y0 - reach * tile[0]) & (sy < y1 + reach * tile[0])
            ns[i, j] = int(m.sum())
            ok = False
            if m.sum() >= min_stars:
                p, f = _combine_cutouts([stack[k] for k in np.nonzero(m)[0]], size)
                if f >= min_fwhm_pix:
                    arr[i, j] = p; fs[i, j] = True; fwg[i, j] = f; ok = True
            if not ok:
                arr[i, j] = constant[0]; fwg[i, j] = constant[1]
    return PSFField((ny, nx), tile, arr, ns, fs, fwg, constant[0], constant[1])




def gaussian_psf(fwhm, size=25):
    r = size // 2
    yy, xx = np.mgrid[-r:r + 1, -r:r + 1]
    s = fwhm / 2.3548
    p = np.exp(-(xx ** 2 + yy ** 2) / (2 * s * s))
    return p / p.sum()
