#!/usr/bin/env python3
"""Fit one de Vaucouleurs or exponential disc on catalog sources.

The menu names map onto presets already in ogfkit.multifit:

- elliptical, dev -> preset dev, Sérsic index fixed at 4
- spiral, exp -> preset exp, Sérsic index fixed at 1

The arm preset named spiral in that module is not used. Bulge+disk is not used.
Magnitudes use the given zeropoint, the same default 25 as plugins/multifit/multifit.py.
"""
import argparse
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit import imageio, multifit as mf, tsvio  # noqa: E402

PRESETS = {"elliptical": "dev", "dev": "dev", "spiral": "exp", "exp": "exp"}
FIXED_N = {"dev": 4.0, "exp": 1.0}
RESULT_COLUMNS = ["NUMBER", "GM_MODEL", "GM_MAG", "GM_RE", "GM_N", "GM_Q", "GM_PA", "GM_CHI2", "GM_FLAG"]
CATALOG_COLUMNS = ["GM_MODEL", "GM_MAG", "GM_RE", "GM_N", "GM_Q", "GM_PA", "GM_CHI2"]


def _num(row, key, default=float("nan")):
    return tsvio.fnum(row.get(key), default)


def start_values(row, zp):
    """Same catalog columns as plugins/multifit/multifit.py start_values."""
    x = _num(row, "X_IMAGE") - 1.0
    y = _num(row, "Y_IMAGE") - 1.0
    flux = _num(row, "FLUX_AUTO")
    if not (flux > 0):
        mag = _num(row, "MAG_AUTO")
        flux = 10 ** (-0.4 * (mag - zp)) if math.isfinite(mag) and mag < 90 else float("nan")
    a_image = _num(row, "A_IMAGE", 2.0)
    b_image = _num(row, "B_IMAGE", 2.0)
    if not (a_image > 0):
        a_image = 2.0
    if not (b_image > 0):
        b_image = a_image
    theta = _num(row, "THETA_IMAGE", 0.0)
    flux_radius = _num(row, "FLUX_RADIUS")
    re = max(flux_radius * 1.1 if flux_radius > 0 else 1.5 * a_image, 1.0)
    return dict(
        x=x, y=y, flux=flux, A=a_image, B=min(b_image, a_image),
        pa=theta if math.isfinite(theta) else 0.0, re=re,
        q=min(max(b_image / a_image, 0.15), 1.0),
        cs=_num(row, "CLASS_STAR", 0.0),
    )


def ellipse_mask(shape, x, y, a, b, pa_deg):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    theta = math.radians(pa_deg)
    dx, dy = xx - x, yy - y
    u = dx * math.cos(theta) + dy * math.sin(theta)
    v = -dx * math.sin(theta) + dy * math.cos(theta)
    return (u / max(a, 0.5)) ** 2 + (v / max(b, 0.5)) ** 2 <= 1.0


def _load_psf(path):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None and np.ndim(hdu.data) == 2:
                return np.asarray(hdu.data, dtype=np.float64)
    raise ValueError("PSF FITS has no 2-D image")


def _neighbour_mask(shape, starts, target, x0, y0):
    mask = np.zeros(shape, dtype=bool)
    tx = starts[target]["x"] - x0
    ty = starts[target]["y"] - y0
    for j, src in enumerate(starts):
        if j == target or not (src["flux"] > 0):
            continue
        x = src["x"] - x0
        y = src["y"] - y0
        if math.hypot(x - tx, y - ty) < 4.0:
            continue
        mask |= ellipse_mask(shape, x, y, 2.5 * src["A"], 2.5 * src["B"], src["pa"])
    mask &= ~ellipse_mask(shape, tx, ty, 1.5, 1.5, 0.0)
    return mask


def fit_sources(image, rows, preset, zp, psf=None, rms=0.0, max_objects=100, max_nfev=150):
    """Fit the brightest extended sources. Returns (records, n_fit)."""
    if preset not in FIXED_N:
        raise ValueError("model must be elliptical or spiral")
    starts = [start_values(row, zp) for row in rows]
    chosen = []
    for i, src in enumerate(starts):
        if not (src["flux"] > 0 and math.isfinite(src["x"]) and math.isfinite(src["y"])):
            continue
        if src["cs"] > 0.8:
            continue
        chosen.append(i)
    chosen.sort(key=lambda i: -starts[i]["flux"])
    if max_objects:
        chosen = chosen[:max_objects]
    sigma = float(rms) if rms and rms > 0 else imageio.robust_sigma(image)
    if not (sigma > 0) or not math.isfinite(sigma):
        sigma = 1.0
    ny, nx = image.shape
    records = []
    for i in chosen:
        src = starts[i]
        half = int(min(max(6.0 * src["A"], 15.5), 60.5))
        x_lo = int(round(src["x"])) - half
        x_hi = int(round(src["x"])) + half + 1
        y_lo = int(round(src["y"])) - half
        y_hi = int(round(src["y"])) + half + 1
        cx0, cx1 = max(0, x_lo), min(nx, x_hi)
        cy0, cy1 = max(0, y_lo), min(ny, y_hi)
        cut = np.asarray(image[cy0:cy1, cx0:cx1], dtype=np.float64)
        number = rows[i].get("NUMBER", "")
        base = {"NUMBER": number, "GM_MODEL": preset, "GM_N": FIXED_N[preset]}
        if cut.size < 9 or not np.isfinite(cut).any():
            records.append(dict(base, GM_MAG=None, GM_RE=None, GM_Q=None, GM_PA=None, GM_CHI2=None, GM_FLAG=128))
            continue
        mask = _neighbour_mask(cut.shape, starts, i, cx0, cy0)
        mask |= ~np.isfinite(cut)
        try:
            result, _comps = mf.fit_preset(
                cut, preset, src["x"] - cx0, src["y"] - cy0, src["flux"], src["re"], src["q"], src["pa"],
                psf=psf, rms=sigma, mask=mask, sky="const", zp=zp, max_nfev=max_nfev,
            )
            comp = result["components"][0]
            records.append(dict(
                base,
                GM_MAG=comp.get("mag"),
                GM_RE=comp.get("re"),
                GM_N=comp.get("n", FIXED_N[preset]),
                GM_Q=comp.get("q"),
                GM_PA=comp.get("pa"),
                GM_CHI2=result.get("chi2_red"),
                GM_FLAG=int(result.get("flags", 0)),
            ))
        except Exception as exc:
            print("galaxy fit %s: %s" % (number, exc), file=sys.stderr)
            records.append(dict(base, GM_MAG=None, GM_RE=None, GM_Q=None, GM_PA=None, GM_CHI2=None, GM_FLAG=128))
    n_fit = sum(1 for rec in records if rec.get("GM_MAG") is not None)
    return records, n_fit


def apply_catalog(columns, rows, results):
    """Copy GM_* cells onto catalog rows with the same NUMBER."""
    by_number = {str(rec.get("NUMBER")): rec for rec in results}
    out_columns = list(columns)
    for name in CATALOG_COLUMNS:
        if name not in out_columns:
            out_columns.append(name)
    merged = []
    for row in rows:
        new = {col: row.get(col, "") for col in out_columns}
        src = by_number.get(str(row.get("NUMBER")))
        if src is not None:
            for name in CATALOG_COLUMNS:
                new[name] = src.get(name, "")
        merged.append(new)
    return out_columns, merged


def _emit_table(columns, records, fh):
    fh.write("\t".join(columns) + "\n")
    for record in records:
        fh.write("\t".join(tsvio.fmt(record.get(column)) for column in columns) + "\n")


def _write_table(path, columns, records):
    parent = os.path.dirname(os.path.abspath(path))
    if parent:
        os.makedirs(parent, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        _emit_table(columns, records, handle)


def _read_results(path):
    columns, rows = tsvio.read_catalog(path)
    if "NUMBER" not in columns:
        raise ValueError("results table has no NUMBER column")
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description="De Vaucouleurs or exponential disc fit for the galaxy menu.")
    parser.add_argument("image", nargs="?", default="")
    parser.add_argument("--catalog", default="")
    parser.add_argument("--model", default="")
    parser.add_argument("--output", default="")
    parser.add_argument("--psf", default="")
    parser.add_argument("--mag-zeropoint", type=float, default=25.0)
    parser.add_argument("--max-objects", type=int, default=100)
    parser.add_argument("--max-nfev", type=int, default=150)
    parser.add_argument("--rms", type=float, default=0.0)
    parser.add_argument("--apply-catalog", default="")
    parser.add_argument("--results", default="")
    args = parser.parse_args(argv)

    if args.apply_catalog:
        if not args.results:
            print("ERROR: --results is required with --apply-catalog", file=sys.stderr)
            return 2
        if not os.path.isfile(args.apply_catalog) or not os.path.isfile(args.results):
            print("ERROR: catalog or results file was not found", file=sys.stderr)
            return 1
        columns, rows = tsvio.read_catalog(args.apply_catalog)
        merged_columns, merged_rows = apply_catalog(columns, rows, _read_results(args.results))
        if args.output:
            _write_table(args.output, merged_columns, merged_rows)
        _emit_table(merged_columns, merged_rows, sys.stdout)
        return 0

    preset = PRESETS.get(str(args.model).strip().lower())
    if preset is None:
        print("ERROR: model must be elliptical or spiral", file=sys.stderr)
        return 2
    if not args.image or not args.catalog:
        print("ERROR: an image and a catalog are required", file=sys.stderr)
        return 2
    if not os.path.isfile(args.image) or not os.path.isfile(args.catalog):
        print("ERROR: image or catalog was not found", file=sys.stderr)
        return 1
    if args.psf and not os.path.isfile(args.psf):
        print("ERROR: PSF file was not found", file=sys.stderr)
        return 1
    image, _header = imageio.load_image(args.image)
    _columns, rows = tsvio.read_catalog(args.catalog)
    if not rows or "NUMBER" not in _columns:
        print("ERROR: catalog has no NUMBER rows", file=sys.stderr)
        return 1
    psf = _load_psf(args.psf) if args.psf else None
    records, n_fit = fit_sources(
        image, rows, preset, args.mag_zeropoint, psf=psf, rms=args.rms,
        max_objects=args.max_objects, max_nfev=args.max_nfev,
    )
    if args.output:
        parent = os.path.dirname(os.path.abspath(args.output))
        if parent:
            os.makedirs(parent, exist_ok=True)
        _write_table(args.output, RESULT_COLUMNS, records)
    summary = (
        "#GALAXY_FIT\tMODEL=%s\tN=%s\tN_FIT=%d\tN_SOURCES=%d\tZP=%s"
        % (preset, FIXED_N[preset], n_fit, len(rows), args.mag_zeropoint)
    )
    print(summary)
    print("Saved galaxy fit: %s" % (args.output or "stdout"), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
