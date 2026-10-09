#!/usr/bin/env python3
"""Fit caller-supplied rest-frame templates with ogfmeas.photoz.template_chi2.

One non-negative scale is fit for each template and redshift. This is not a
stellar-population code. It does not read a template library or a filter
response file.

stdin JSON
----------
templates:
    list of [rest_wavelength, rest_flux]. Required. A file path is refused.
redshifts:
    optional list. Otherwise z_min, z_max, z_step (defaults 0.01, 6, 0.01).
ids:
    one id per object.
wavelength, flux, err:
    spectrum samples. flux and err have one row per id, or one 1-d spectrum
    when there is one id. The wavelength unit must match the templates.
bands, flux, err:
    used when wavelength is absent. Each known band is sampled at the
    registry pivot (micron converted to Angstrom). That is a coarse point
    sample, not a filter integral.

Prints {"version": "ogfmeas.photoz", "rows": [...]} on success.
Without numeric templates, exits 2 and writes the reason to stderr.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

_REFUSAL = (
    "the external package is not called; numeric rest-frame templates "
    "were not supplied. The in-tree fitter is ogfmeas.photoz.template_chi2.\n"
)


def _templates(req):
    raw = req.get("templates")
    if not isinstance(raw, list) or not raw:
        return None
    out = []
    for item in raw:
        if not isinstance(item, (list, tuple)) or len(item) != 2:
            return None
        wave, flux = item
        if isinstance(wave, str) or isinstance(flux, str):
            return None
        w = np.asarray(wave, dtype=np.float64)
        f = np.asarray(flux, dtype=np.float64)
        if w.ndim != 1 or w.shape != f.shape or w.size < 2:
            return None
        out.append((w, f))
    return out


def _redshifts(req):
    if req.get("redshifts") is not None:
        zs = np.asarray(req["redshifts"], dtype=np.float64)
        if zs.ndim != 1 or zs.size < 1:
            raise ValueError("redshifts must be a non-empty list")
        return zs
    z_min = float(req.get("z_min", 0.01))
    z_max = float(req.get("z_max", 6.0))
    z_step = float(req.get("z_step", 0.01))
    if not (z_step > 0) or not (z_max >= z_min):
        raise ValueError("redshift grid is empty")
    n = int(np.floor((z_max - z_min) / z_step + 0.5)) + 1
    return z_min + np.arange(n, dtype=np.float64) * z_step


def _matrix(value, n_id):
    a = np.asarray(value, dtype=np.float64)
    if a.ndim == 1:
        if n_id != 1:
            raise ValueError("flux and err must have one row per id")
        a = a.reshape(1, -1)
    if a.ndim != 2 or a.shape[0] != n_id:
        raise ValueError("flux and err must have one row per id")
    return a


def _observed(req, n_id):
    if req.get("wavelength") is not None:
        wave = np.asarray(req["wavelength"], dtype=np.float64)
        flux = _matrix(req["flux"], n_id)
        err = _matrix(req["err"], n_id)
        if wave.ndim != 1 or flux.shape[1] != wave.size or err.shape != flux.shape:
            raise ValueError("wavelength, flux, and err must share the sample axis")
        return wave, flux, err
    bands = req.get("bands")
    if not bands:
        raise ValueError("wavelength or bands are required")
    from sed_adapters import filters as reg
    wave_list = []
    keep = []
    for j, band in enumerate(bands):
        entry = reg.B.get(str(band).upper())
        if entry is None:
            continue
        wave_list.append(float(entry["lam"]) * 10000.0)
        keep.append(j)
    if len(keep) < 3:
        raise ValueError("need at least three registry bands for a pivot sample")
    flux = _matrix(req["flux"], n_id)[:, keep]
    err = _matrix(req["err"], n_id)[:, keep]
    return np.asarray(wave_list, dtype=np.float64), flux, err


def _summary(zgrid, pdf):
    c = np.cumsum(np.asarray(pdf, dtype=np.float64))
    total = float(c[-1]) if c.size else 0.0
    if not (total > 0):
        return None, None, None, None
    c = c / total
    def q(frac):
        return float(np.interp(frac, c, zgrid))
    p16, p50, p84 = q(0.1587), q(0.5), q(0.8413)
    return p16, p50, p84, 0.5 * (p84 - p16)


def main():
    try:
        req = json.load(sys.stdin)
    except Exception:
        sys.stderr.write(_REFUSAL)
        return 2
    if not isinstance(req, dict):
        sys.stderr.write(_REFUSAL)
        return 2
    templates = _templates(req)
    if templates is None:
        sys.stderr.write(_REFUSAL)
        return 2
    from ogfmeas.photoz import template_chi2
    ids = list(req.get("ids") or ["1"])
    try:
        zs = _redshifts(req)
        wave, flux, err = _observed(req, len(ids))
    except (KeyError, TypeError, ValueError) as exc:
        sys.stderr.write("the external package is not called; %s\n" % exc)
        return 2
    rows = []
    for i in range(len(ids)):
        try:
            fit = template_chi2(wave, flux[i], err[i], templates, zs)
        except ValueError as exc:
            rows.append({"z": None, "error": str(exc)})
            continue
        p16, p50, p84, sigma = _summary(fit["zgrid"], fit["pdf"])
        n_good = int(np.sum(np.isfinite(flux[i]) & np.isfinite(err[i]) & (err[i] > 0)))
        rows.append({
            "z": float(fit["z"]),
            "p16": p16,
            "p50": p50,
            "p84": p84,
            "sigma": sigma,
            "chi2": float(np.min(fit["chi2"])) / max(n_good - 1, 1),
        })
    json.dump({"version": "ogfmeas.photoz", "rows": rows}, sys.stdout)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
