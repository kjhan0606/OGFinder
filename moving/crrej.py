"""Cosmic-ray handling for the detection stage (L.A.Cosmic) and CR-aware per-pixel utilities.

* `lac_mask(chip)`: van Dokkum (2001) Laplacian-edge detection through `astroscrappy` (a C implementation of the same algorithm; used when
  importable, `pip install astroscrappy` in the venv), else the simplified numpy version `imaging.lacosmic_mask`.  The mask is independent of
  the archive pipeline flag (`chip.cr`, DQ bit 4096), which is dense on the BB89 chips and also flags real movers.
* `inpaint(data, mask)`: replace masked pixels by the median of the unmasked pixels of a 5x5 neighbourhood (iterated).
* `consistent_cr(...)`: multi-exposure consistency - a pixel excess of a target exposure that is absent in the other exposures cannot be
  told from a mover by position alone, so this is only used as a *feature* (`cr_novel`) of the real/bogus classifier, not as a veto.

A CR mask is NOT a veto for detections: with a Gaussian-PSF mover of FWHM 1.7 px, 65 % of the injected non-trail detections sit on L.A.Cosmic
pixels (124 of 191 in inj4), so a hard veto would remove most faint real movers (measured, docs/moving_objects.md).  The mask enters the
real/bogus classifier (`realbogus.py`) and the centroid refinement (`centroid.py`, CR pixels are down-weighted)."""
import numpy as np
from scipy import ndimage as ndi
from .util import log

try:                                                    # pragma: no cover - availability depends on the venv
    import astroscrappy
except Exception:
    astroscrappy = None


def backend():
    return "astroscrappy" if astroscrappy is not None else "numpy"


def lac_mask(data_eps, bad, texp, gain=1.0, readnoise=5.0, sigclip=4.5, sigfrac=0.3, objlim=5.0, niter=4, grow=False):
    """L.A.Cosmic mask of an image in e-/s (`data_eps`), `bad` = pixels excluded (boolean), exposure time `texp` [s]."""
    img = (np.asarray(data_eps, np.float32) * float(texp)).astype(np.float32)
    bad = np.asarray(bad, bool)
    if astroscrappy is not None:
        m, _ = astroscrappy.detect_cosmics(img, inmask=bad, gain=float(gain), readnoise=float(readnoise), satlevel=1e12, sigclip=sigclip,
                                           sigfrac=sigfrac, objlim=objlim, niter=niter, verbose=False)
        m = np.asarray(m, bool)
    else:
        from .imaging import lacosmic_mask
        m = lacosmic_mask(img, sigclip=sigclip, sigfrac=sigfrac, objlim=objlim, niter=min(niter, 3), gain=gain, readnoise=readnoise, bad=bad)
    if grow:
        m = ndi.binary_dilation(m, iterations=1)
    return m


def chip_lac(chip, readnoise=None, **kw):
    """LAC mask of a `Chip` (ACS: read noise from READNSEA..D header keywords, else 5 e-)."""
    rn = readnoise
    if rn is None:
        h0 = getattr(chip, "header0", None) or {}
        try:
            rn = float(np.mean([float(h0["READNSE" + k]) for k in "ABCD" if "READNSE" + k in h0]))
        except Exception:
            rn = 5.0
        if not np.isfinite(rn) or rn <= 0:
            rn = 5.0
    gain = float(getattr(chip, "gain", 1.0) or 1.0)
    gain = 1.0 if gain <= 0 or gain > 10 else gain          # data of HST/ACS flc are already in electrons
    return lac_mask(chip.data, chip.bad | getattr(chip, "saturated", np.zeros(chip.shape, bool)), chip.texp, gain=1.0, readnoise=rn, **kw)


def inpaint(data, mask, size=5, iters=3):
    """Median of the unmasked neighbours (size x size) in place of masked pixels; returns a new array."""
    out = np.array(data, np.float32)
    m = np.asarray(mask, bool).copy()
    if not m.any():
        return out
    for _ in range(iters):
        if not m.any():
            break
        w = np.where(m, np.nan, out)
        # nan-aware median via a generic filter on the (small) set of masked pixels only
        ys, xs = np.nonzero(m)
        h = size // 2
        pad = np.pad(w, h, constant_values=np.nan)
        stack = np.stack([pad[h + dy + ys, h + dx + xs] for dy in range(-h, h + 1) for dx in range(-h, h + 1)], 0)
        with np.errstate(all="ignore"):
            med = np.nanmedian(stack, axis=0)
        ok = np.isfinite(med)
        out[ys[ok], xs[ok]] = med[ok]
        m[ys[ok], xs[ok]] = False
    return out


def novelty(target_data, others_stack, sig):
    """Multi-exposure consistency map: (target - median(others)) / sig where the other exposures (reprojected stack, NaN = no data) agree
    with each other (their scatter is small); a CR hit is a single-exposure excess, a PSF-shaped static source is present in all.
    Used as a feature only."""
    with np.errstate(all="ignore"):
        med = np.nanmedian(others_stack, axis=0)
    return (target_data - med) / np.maximum(sig, 1e-6)
