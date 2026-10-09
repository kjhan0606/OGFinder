#!/usr/bin/env python3
"""Fit one elliptical Sérsic on case-A cutouts with ogfmeas.

    python pysersic_map.py WORKDIR OUT.json [--n 30]

The filename is historical. This script fits with ogfmeas.sersic only. It does
not load an external profile package. Fourier modes, a sky gradient, and more
than one component are refused.

Each output row uses 1-based x and y. ``pa`` is the feedme position angle
(degrees, +y toward -x) folded into [-90, 90). ``mag`` is
zp - 2.5 log10(total flux) with the flux integral in ogfmeas.sersic.
"""
import json
import math
import os
import sys
import time

import numpy as np
from astropy.io import fits

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from ogfkit import galfitio as GI  # noqa: E402
from ogfmeas.sersic import (  # noqa: E402
    fit_sersic, galfit_pa_from_theta, ie_from_flux, render, total_flux,
)


class SingleSersicOnly(ValueError):
    pass


def _shape_keys(component):
    return [key for key in component if len(key) > 1 and key[0] in "fb" and key[1].isdigit()]


def _read_start(wd):
    path = os.path.join(wd, "start.feedme")
    if not os.path.isfile(path):
        raise SingleSersicOnly("start.feedme is missing")
    cfg = GI.parse_feedme(path, strict=False, base_dir=wd)
    comps = cfg.get("components") or []
    if len(comps) != 1 or comps[0].get("galfit_type") != "sersic":
        raise SingleSersicOnly("single elliptical Sérsic only")
    component = comps[0]
    if "c0" in component or _shape_keys(component):
        raise SingleSersicOnly("Fourier, boxiness, and bending modes are out of scope")
    grad = cfg.get("sky_grad") or [0.0, 0.0]
    if abs(float(grad[0])) > 1e-8 or abs(float(grad[1])) > 1e-8:
        raise SingleSersicOnly("a sky gradient is out of scope")
    return cfg, component


def _image(wd, name):
    path = os.path.join(wd, name)
    if not os.path.isfile(path):
        return None
    return np.asarray(fits.getdata(path), dtype=np.float64)


def fit_directory(wd, sigma=None):
    """Fit one cutout. x and y in the result are 1-based."""
    data = _image(wd, "data.fits")
    if data is None:
        raise FileNotFoundError(os.path.join(wd, "data.fits"))
    if data.ndim == 3:
        data = data[0]
    cfg, component = _read_start(wd)
    zp = float(cfg.get("zp", 25.0))
    psf = _image(wd, "psf.fits")
    mask_image = _image(wd, "mask.fits")
    mask = None if mask_image is None else mask_image > 0
    if sigma is None:
        sigma = _image(wd, "sigma.fits")
    flux = 10.0 ** (-0.4 * (float(component["mag"]) - zp))
    sky0 = cfg.get("sky_value_galfit")
    fit = fit_sersic(
        data,
        x0=float(component["x"]) - 1.0,
        y0=float(component["y"]) - 1.0,
        mask=mask,
        psf=psf,
        sigma=sigma,
        re0=float(component["re"]),
        n0=float(component["n"]),
        q0=float(component["q"]),
        theta0=math.radians(float(component["pa"])),
        ie0=ie_from_flux(flux, component["re"], component["n"], component["q"]),
        background0=0.0 if sky0 is None else float(sky0),
    )
    integrated = total_flux(fit["ie"], fit["re"], fit["n"], fit["q"])
    mag = 99.0 if integrated <= 0.0 else zp - 2.5 * math.log10(integrated)
    fit["chi2_fit"] = fit.pop("chi2")
    fit.update(x=fit["x"] + 1.0, y=fit["y"] + 1.0, mag=mag, pa=galfit_pa_from_theta(fit["theta"]),
               sky=fit["background"], zp=zp, fitter="ogfmeas.sersic")
    return fit


def pixel_chi2(wd, reported, sigma=None):
    """Sum of squared residuals of one reported single Sérsic, and the pixel count.

    The model is ogfmeas. Reported x and y are 1-based. When ``ie`` and
    ``theta`` are present they are used; otherwise magnitude and ``pa`` are
    converted with the ogfmeas flux integral and the feedme angle.
    """
    data = _image(wd, "data.fits")
    if data is None:
        raise FileNotFoundError(os.path.join(wd, "data.fits"))
    if data.ndim == 3:
        data = data[0]
    if sigma is None:
        sigma = _image(wd, "sigma.fits")
    if sigma is None:
        sigma = 1.0
    zp = float(reported.get("zp", 25.0))
    if "ie" in reported and "theta" in reported:
        ie = float(reported["ie"])
        theta = float(reported["theta"])
    else:
        flux = 10.0 ** (-0.4 * (float(reported["mag"]) - zp))
        ie = ie_from_flux(flux, reported["re"], reported["n"], reported["q"])
        theta = math.radians(float(reported["pa"]) + 90.0)
    model = render(
        data.shape, float(reported["x"]) - 1.0, float(reported["y"]) - 1.0,
        float(reported["re"]), float(reported["n"]), ie, float(reported["q"]), theta,
        float(reported.get("sky", reported.get("background", 0.0))), psf=_image(wd, "psf.fits"),
    )
    mask_image = _image(wd, "mask.fits")
    bad = ~np.isfinite(data) | ~np.isfinite(model)
    if mask_image is not None:
        bad |= mask_image > 0
    resid = (data - model) / np.maximum(np.asarray(sigma, dtype=np.float64), 1e-8)
    use = (~bad) & np.isfinite(resid)
    return float(np.sum(resid[use] ** 2)), int(np.count_nonzero(use))


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if "--n" in argv:
        n = int(argv[argv.index("--n") + 1])
        del argv[argv.index("--n") + 1]
        argv.remove("--n")
    else:
        n = 30
    if len(argv) != 2:
        sys.exit("usage: pysersic_map.py WORKDIR OUT.json [--n 30]")
    work, out = argv
    rows = []
    for seed in range(n):
        wd = os.path.join(work, "A%02d" % seed)
        started = time.time()
        try:
            row = fit_directory(wd)
            row.update(seed=seed, time=time.time() - started)
        except Exception as exc:
            row = dict(seed=seed, error=str(exc)[:200], time=time.time() - started, fitter="ogfmeas.sersic")
        rows.append(row)
        json.dump(rows, open(out, "w"))
        print(seed, row.get("mag"), "%.1fs" % row["time"], flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
