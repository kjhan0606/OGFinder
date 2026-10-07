"""Desktop driver of the shared-core aperture photometry (plugins/ds10core/run_aper_phot.py): clicked ds9 regions -> star list, single image
and time series through the driver as the step runs it (separate process; the core is optional and the run tests skip without it)."""
import csv
import math
import os
import subprocess
import sys

import numpy as np
import pytest

import ds10 as launcher
import run_aper_phot as drv

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
CORE, WHY = launcher.find_core()
needs_core = pytest.mark.skipif(CORE is None, reason="ds10core (ds10-web/server) not found: " + str(WHY))
XY = [(40.0, 50.0), (120.0, 60.0), (90.0, 130.0), (150.0, 150.0)]
FL = [5e4, 7e4, 6e4, 4e4]


def star_image(path, xy, flux, fwhm=3.0, seed=0, mjd=60000.0):
    from astropy.io import fits
    rng = np.random.default_rng(seed)
    n = 200
    yy, xx = np.mgrid[0:n, 0:n]
    s = fwhm / 2.3548
    img = np.full((n, n), 100.0)
    for (x, y), f in zip(xy, flux):
        img += f / (2 * math.pi * s * s) * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / (2 * s * s))
    img = rng.poisson(img * 2.0) / 2.0
    h = fits.PrimaryHDU(img.astype(np.float32))
    h.header.update({"GAIN": 2.0, "EXPTIME": 30.0, "MJD-OBS": mjd, "FILTER": "V", "MAGZERO": 25.0, "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CRPIX1": 100.5,
                     "CRPIX2": 100.5, "CRVAL1": 200.0, "CRVAL2": 30.0, "CD1_1": -1e-4, "CD1_2": 0.0, "CD2_1": 0.0, "CD2_2": 1e-4})
    h.writeto(path, overwrite=True)


def test_picks_from_ds9_regions():
    txt = ("# Region file format: DS9 version 4.1\nglobal color=green\nimage\ncircle(41.2,51.3,5) # text={T}\npoint(121,61) # point=circle\n"
           "box(10,10,4,4,0)\nfk5\ncircle(200.0,30.0,2\")\n")
    p = drv.picks_from_regions(txt)
    assert p == [(41.2, 51.3, "T"), (121.0, 61.0, "")]


def _rows(p):
    with open(p, encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@needs_core
def test_driver_single_image_from_clicked_regions(tmp_path):
    img = tmp_path / "s.fits"
    star_image(img, XY, FL)
    picks = tmp_path / "picks.reg"
    picks.write_text("image\n" + "".join(f"point({x + 1.8:.1f},{y + 0.3:.1f})\n" for x, y in XY))     # clicks ~1.8 px off (1-based)
    out = tmp_path / "cat.tsv"
    r = subprocess.run([sys.executable, os.path.join(PLUGIN, "run_aper_phot.py"), "--stars", "regions", "--picks", str(picks), str(img), "--out", str(out),
                        "--regions", str(tmp_path / "ap.reg"), "--radii", "2", "--unit", "fwhm", "--apcor", "infinity"], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    assert r.stdout.startswith("NUMBER\t")                        # catalogue echoed for the table (output mode set)
    rows = _rows(out)
    assert len(rows) == 4
    for row, (x, y), f in zip(rows, XY, FL):
        assert abs(float(row["X_IMAGE"]) - 1 - x) < 0.1 and abs(float(row["Y_IMAGE"]) - 1 - y) < 0.1
        assert abs(float(row["FLUX"]) / f - 1) < 0.02
        assert abs(float(row["MAG"]) - (25 - 2.5 * math.log10(f))) < 0.03
    assert (tmp_path / "ap.reg").read_text().count("annulus(") == 4


@needs_core
def test_driver_time_series_and_plot(tmp_path):
    paths = []
    for k in range(4):
        f = list(FL)
        f[0] *= 1 + 0.1 * k
        p = tmp_path / f"f{k}.fits"
        star_image(p, [(x + 0.5 * k, y) for x, y in XY], f, seed=k, mjd=60000.0 + 0.01 * k)
        paths.append(str(p))
    picks = tmp_path / "picks.reg"
    picks.write_text("image\n" + "".join(f"circle({x + 1},{y + 1},5)\n" for x, y in XY))
    r = subprocess.run([sys.executable, os.path.join(PLUGIN, "run_aper_phot.py"), "--stars", "regions", "--picks", str(picks), "--frames", ",".join(paths[1:]),
                        "--lc-png", str(tmp_path / "lc.png"), paths[0], "--out", str(tmp_path / "c.tsv"), "--lc-out", str(tmp_path / "lc.tsv"),
                        "--aavso-out", str(tmp_path / "a.txt"), "--target", "1", "--comps", "2,3", "--check", "4"], capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-2000:]
    lc = _rows(tmp_path / "lc.tsv")
    rel = np.array([float(x["REL_FLUX"]) for x in lc])
    assert np.allclose(rel / rel[0], [1.0, 1.1, 1.2, 1.3], atol=0.02)
    try:
        import matplotlib  # noqa: F401
        assert (tmp_path / "lc.png").stat().st_size > 5000
    except ImportError:
        pass
