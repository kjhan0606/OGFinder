#!/usr/bin/env python3
"""Per-pixel ridge of each SDSS detection in a heldout.py JSON.

heldout.py's ``amp`` column is ``profile_snr`` (fitted cross-track amplitude
divided by the error of the along-track median).  This script measures the
ridge itself, in units of the local pixel noise:

* sample the frame along the reported endpoint chord (rho is not in the JSON)
* background at each along-track sample: median of the wings 28 <= |t| <= 46 px
* sigma: 1.4826 * MAD of those wing pixels
* drop the brightest 3 percent of along-track columns (a star on the line)
* ridge = peak of the remaining cross-track median, minus the median at |t| > 25,
  divided by sigma
* FWHM = full width of that profile at half the ridge
* core/amp = brightest excess inside |t| <= 3, searched 1.5 L past each end,
  divided by the stored profile amplitude in image units.  sep is the distance
  of that peak from the segment, in pixels (0 if the peak lies on it).

usage: TRAILS_SDSS=<frame dir> python3 sdss_ridge.py --json heldout.json
"""
import argparse, json, math, os
import numpy as np
from astropy.io import fits
from scipy.ndimage import map_coordinates

SDSS = os.environ.get("TRAILS_SDSS", "/workspace/work/item4/sdss")


def load(band, field, cache):
    key = (band, field)
    if key not in cache:
        path = os.path.join(SDSS, "frame-%s-000094-3-%04d.fits" % (band, field))
        with fits.open(path) as hdul:
            cache[key] = np.asarray(hdul[0].data, np.float32)
    return cache[key]


def sample(img, x1, y1, x2, y2, s, t):
    length = math.hypot(x2 - x1, y2 - y1)
    ux, uy = (x2 - x1) / length, (y2 - y1) / length
    px, py = -uy, ux
    ss, tt = np.meshgrid(s, t)
    return map_coordinates(img, [y1 + ss * uy + tt * py, x1 + ss * ux + tt * px], order=1, mode="constant", cval=np.nan)


def ridge(img, det):
    x1, y1 = det["x1"] - 1.0, det["y1"] - 1.0
    x2, y2 = det["x2"] - 1.0, det["y2"] - 1.0
    length = math.hypot(x2 - x1, y2 - y1)
    s_seg = np.arange(0.0, length + 1.0, 1.0)
    s_ext = np.arange(-1.5 * length, length * 2.5 + 1.0, 1.0)
    t = np.arange(-48.0, 49.0, 1.0)
    with np.errstate(all="ignore"):
        seg = sample(img, x1, y1, x2, y2, s_seg, t)
        ext = sample(img, x1, y1, x2, y2, s_ext, t)

    def sl(a, b, arr=seg):
        i0 = int(round(a - t[0]))
        i1 = int(round(b - t[0])) + 1
        return arr[max(i0, 0):i1]

    wing = np.concatenate([sl(28, 46), sl(-46, -28)], axis=0)
    finite = wing[np.isfinite(wing)]
    sigma = 1.4826 * np.median(np.abs(finite - np.median(finite))) if finite.size else float("nan")
    bkg = np.nanmedian(wing, axis=0)
    centre = np.nanmedian(sl(-2, 2), axis=0) - bkg
    keep = np.isfinite(centre)
    if keep.sum() > 20:
        keep &= centre < np.nanpercentile(centre[keep], 97)
    cross = np.nanmedian((seg - bkg)[:, keep], axis=1) if keep.any() else np.full(t.shape, np.nan)
    outer = np.abs(t) > 25
    base = np.nanmedian(cross[outer]) if np.isfinite(cross[outer]).any() else 0.0
    peak = np.nanmax(cross) if np.isfinite(cross).any() else float("nan")
    amp = peak - base
    half = base + 0.5 * amp
    above = t[np.isfinite(cross) & (cross >= half)] if np.isfinite(amp) else np.array([])
    fwhm = float(above.max() - above.min()) if above.size else float("nan")
    # fraction of 20 px bins whose centre excess exceeds 0.5 sigma
    c = centre / sigma
    nb = int(len(c) // 20)
    bins = np.array([np.nanmedian(c[i * 20:(i + 1) * 20]) for i in range(nb)]) if nb else np.array([np.nan])
    # core search on the extended strip, image units, same 1.5 L window as bleed_like
    def rows_of(arr, a, b):
        i0 = int(round(a - t[0]))
        i1 = int(round(b - t[0])) + 1
        return arr[max(i0, 0):i1]
    wing_e = np.concatenate([rows_of(ext, 28, 46), rows_of(ext, -46, -28)], axis=0)
    with np.errstate(all="ignore"):
        bkg_e = np.nanmedian(wing_e, axis=0) if np.isfinite(wing_e).any() else np.zeros(ext.shape[1])
        core_col = np.nanmax(rows_of(ext, -3, 3), axis=0) - bkg_e
    ok = np.isfinite(core_col)
    stored = float(det["amp"]) if det.get("amp") else float("nan")
    if ok.any() and stored > 0:
        k = int(np.nanargmax(np.where(ok, core_col, -np.inf)))
        core_s = float(s_ext[k])
        core_ratio = float(core_col[k]) / stored
        core_sep = 0.0 if 0.0 <= core_s <= length else float(-core_s if core_s < 0 else core_s - length)
    else:
        core_ratio = core_sep = float("nan")
    return dict(ridge=float(amp / sigma), fwhm=fwhm, frac=float(np.mean(bins > 0.5)),
                core=core_ratio, sep=core_sep, sigma=float(sigma))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--json", required=True)
    a = ap.parse_args()
    raw = json.load(open(a.json))
    rows = []
    for fr in raw["sdss"]:
        for d in fr["det"]:
            rows.append(dict(d, band=fr["band"], field=fr["field"]))
    rows.sort(key=lambda d: (d["field"], d["band"], d["theta"]))
    cache, seen = {}, {}
    print("| band | field | id | theta | length | auto | ridge [sigma_pix] | FWHM [px] | bins>0.5σ | core/amp | sep [px] |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for d in rows:
        key = "%s_%04d_%d%s" % (d["band"], d["field"], d["id"], "c" if d["curved"] else "")
        seen[key] = seen.get(key, 0) + 1
        m = ridge(load(d["band"], d["field"], cache), d)
        auto = "static %d" % d["n_present"] if d["n_present"] else "one-band"
        print("| %s | %d | %d%s | %.1f | %.0f | %s | %.2f | %.0f | %.2f | %.1f | %.0f |" % (
            d["band"], d["field"], d["id"], "c" if d["curved"] else "", d["theta"], d["length"], auto,
            m["ridge"], m["fwhm"], m["frac"], m["core"], m["sep"]))


if __name__ == "__main__":
    main()
