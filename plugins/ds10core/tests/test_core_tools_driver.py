"""Desktop drivers for the licence-clean shared-core tools that were missing from the menu.

Each driver starts ``python -m ds10core.tools.<name>`` in a separate process. These tests never import ds10core.
Runs that need the core are skipped when the Astrafex Web checkout is not next to OGFinder.
"""
import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

import ds10 as launcher

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
CORE, WHY = launcher.find_core()
needs_core = pytest.mark.skipif(CORE is None, reason="ds10core (ds10-web/server) not found: " + str(WHY))
_FWHM = 2.3548200450309493
_D2R = math.pi / 180.0
_R2D = 180.0 / math.pi

_NEW = {
    "wcs": ("run_wcs.py", ["image"], "text", "measure"),
    "dither": ("run_dither.py", [], "text", "measure"),
    "regrid": ("run_regrid.py", ["image"], "text", "measure"),
    "crmask": ("run_crmask.py", ["image"], "text", "measure"),
    "detect": ("run_detect.py", ["image"], "set", "detect"),
    "extinct": ("run_extinct.py", ["catalog"], "set", "measure"),
    "period": ("run_period.py", ["catalog"], "set", "measure"),
    "dropout": ("run_dropout.py", ["catalog"], "set", "measure"),
    "counts": ("run_counts.py", ["image"], "set", "measure"),
    "lf": ("run_lf.py", ["catalog"], "set", "measure"),
}


def _env():
    e = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        e[name] = "2"
    return e


def _run(script, args):
    return subprocess.run(
        [sys.executable, os.path.join(PLUGIN, script), "--core-dir", CORE or "", *args],
        capture_output=True, text=True, env=_env(), timeout=180,
    )


def _write_fits(path, data, **cards):
    from astropy.io import fits
    hdu = fits.PrimaryHDU(np.asarray(data, np.float32))
    for key, value in cards.items():
        hdu.header[key] = value
    hdu.writeto(path, overwrite=True)


def _plant(img, x, y, flux, fwhm):
    sig = fwhm / _FWHM
    r = int(math.ceil(5.0 * sig)) + 2
    ny, nx = img.shape
    x0, x1 = max(0, int(math.floor(x)) - r), min(nx, int(math.floor(x)) + r + 1)
    y0, y1 = max(0, int(math.floor(y)) - r), min(ny, int(math.floor(y)) + r + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    img[y0:y1, x0:x1] += flux * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / (sig * sig))


def _tan_cards(nx, ny, ra, dec, scale_arcsec):
    s = float(scale_arcsec) / 3600.0
    return {
        "CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN",
        "CRPIX1": (nx + 1) / 2.0, "CRPIX2": (ny + 1) / 2.0,
        "CRVAL1": float(ra), "CRVAL2": float(dec),
        "CD1_1": -s, "CD1_2": 0.0, "CD2_1": 0.0, "CD2_2": s,
        "RADESYS": "ICRS", "EQUINOX": 2000.0,
    }


def _pix2world(x, y, nx, ny, ra0, dec0, scale_arcsec):
    """Same TAN formulas as the core, so the test catalogue lands on the planted pixels."""
    s = float(scale_arcsec) / 3600.0
    crpix = ((nx + 1) / 2.0, (ny + 1) / 2.0)
    u, v = (x + 1.0) - crpix[0], (y + 1.0) - crpix[1]
    xi, eta = -s * u, s * v
    radius = math.hypot(xi, eta)
    phi = math.atan2(xi, -eta)
    theta = math.atan2(_R2D, radius)
    a0, d0 = ra0 * _D2R, dec0 * _D2R
    dphi = phi - math.pi
    st, ct = math.sin(theta), math.cos(theta)
    sd = st * math.sin(d0) + ct * math.cos(d0) * math.cos(dphi)
    dec = math.asin(max(-1.0, min(1.0, sd)))
    ra = a0 + math.atan2(-ct * math.sin(dphi), st * math.cos(d0) - ct * math.sin(d0) * math.cos(dphi))
    return (ra * _R2D) % 360.0, dec * _R2D


def _flat_flags(step):
    flags = []
    for item in step["cli"]:
        chunk = [item] if isinstance(item, str) else item["argv"]
        flags.extend(flag for flag in chunk if isinstance(flag, str) and flag.startswith("--"))
    return flags


def test_manifest_keeps_the_clean_steps_on_the_process_boundary():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json"), encoding="utf-8"))
    assert m["web"] is False and "separate process" in m["license"]["boundary"]
    steps = {s["id"]: s for s in m["steps"]}
    assert "mask" not in steps and "completeness" not in steps
    params = {p["name"]: p for p in m["params"]}
    for sid, (drv, needs, mode, stage) in _NEW.items():
        step = steps[sid]
        visible = step["label"] + step["title"] + step["done_status"]
        assert "ds10" not in visible, sid
        assert step["record"] == "ds10core." + sid
        assert step["needs"] == needs and step["output"]["mode"] == mode and step["stage"] == stage
        src = open(os.path.join(PLUGIN, drv), encoding="utf-8").read()
        assert "import ds10core" not in src and "from ds10core" not in src
        for flag in _flat_flags(step):
            assert flag in src, (sid, flag)
    launch = open(os.path.join(PLUGIN, "launch.py"), encoding="utf-8").read()
    assert "import ds10core" not in launch and "from ds10core" not in launch
    detect = steps["detect"]["cli"]
    assert "--deblend" not in [item for item in detect if isinstance(item, str)]
    assert params["det-deblend"]["default"] is False
    assert params["counts-deblend"]["default"] is False
    assert params["counts-n-bins"]["default"] == 4 and params["counts-per-bin"]["default"] == 6
    assert params["wcs-survey"]["default"] == "none"
    assert steps["wcs"].get("network") is True
    assert all(not s.get("network") for s in m["steps"] if s["id"] != "wcs")
    assert "--survey" not in [item for item in steps["wcs"]["cli"] if isinstance(item, str)]


@needs_core
def test_wcs_driver_matches_a_catalogue_and_refuses_a_bare_image(tmp_path):
    nx = ny = 64
    ra0, dec0, scale = 150.0, 20.0, 1.0
    rng = np.random.default_rng(7)
    img = np.full((ny, nx), 100.0) + rng.normal(0.0, 0.8, (ny, nx))
    pix = [(16, 18), (40, 20), (22, 44), (48, 46), (30, 32), (12, 40),
           (50, 14), (36, 50), (18, 28), (44, 36), (26, 14), (52, 28)]
    for x, y in pix:
        _plant(img, x, y, 400.0, 2.8)
    path = tmp_path / "field.fits"
    _write_fits(path, img, **_tan_cards(nx, ny, ra0, dec0, scale))
    lines = ["ra\tdec"]
    for x, y in pix:
        ra, dec = _pix2world(x, y, nx, ny, ra0, dec0, scale)
        lines.append(f"{ra:.8f}\t{dec:.8f}")
    cat = tmp_path / "cat.tsv"
    cat.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out, summary = tmp_path / "wcs_solved.fits", tmp_path / "wcs_summary.json"
    bare = _run("run_wcs.py", [str(path), "--out", str(tmp_path / "no.fits"), "--summary", str(tmp_path / "no.json")])
    assert bare.returncode == 2, bare.stderr[-800:]
    assert "catalogue" in bare.stderr.lower() or "survey" in bare.stderr.lower()
    got = _run("run_wcs.py", [str(path), "--out", str(out), "--summary", str(summary), "--cat", str(cat)])
    assert got.returncode == 0, got.stderr[-1500:]
    assert got.stdout.startswith("solved ")
    body = json.loads(summary.read_text())
    assert body["matched"] >= 8 and body["guess"] == "header"


@needs_core
def test_dither_regrid_and_crmask_drivers(tmp_path):
    a = np.full((32, 32), 40.0, np.float32)
    _plant(a, 16, 14, 900.0, 2.5)
    b = np.roll(a, 2, axis=1)
    pa, pb = tmp_path / "a.fits", tmp_path / "b.fits"
    _write_fits(pa, a)
    _write_fits(pb, b)
    stacked = _run("run_dither.py", ["--frames", f"{pa}|{pb}", "--combine", "median", "--nsig", "3",
                                     "--out", str(tmp_path / "stack.fits"), "--summary", str(tmp_path / "dither.json")])
    assert stacked.returncode == 0, stacked.stderr[-1500:]
    assert stacked.stdout.startswith("stacked 2 ")
    src = np.zeros((32, 32), np.float32)
    src[16, 18] = 25.0
    ref = np.zeros((20, 28), np.float32)
    _write_fits(tmp_path / "src.fits", src, **_tan_cards(32, 32, 10.0, 20.0, 1.0))
    _write_fits(tmp_path / "ref.fits", ref, **_tan_cards(28, 20, 10.0, 20.0, 1.0))
    resampled = _run("run_regrid.py", [str(tmp_path / "src.fits"), "--ref", str(tmp_path / "ref.fits"),
                                       "--out", str(tmp_path / "resampled.fits"),
                                       "--summary", str(tmp_path / "regrid.json")])
    assert resampled.returncode == 0, resampled.stderr[-1500:]
    assert resampled.stdout.startswith("resampled ")
    body = json.loads((tmp_path / "regrid.json").read_text())
    assert body["shape"] == [20, 28]
    rng = np.random.default_rng(3)
    sky = (50.0 + rng.normal(0.0, 1.0, (40, 36))).astype(np.float32)
    sky[12, 18] = 40000
    _write_fits(tmp_path / "cr.fits", sky)
    masked = _run("run_crmask.py", [str(tmp_path / "cr.fits"), "--out", str(tmp_path / "crmask.fits"),
                                    "--summary", str(tmp_path / "crmask.json"),
                                    "--cr", "--streak", "--nsig", "5", "--grow", "1", "--max-area", "8",
                                    "--star-fwhm", "3", "--min-length", "40", "--clip-percent", "2"])
    assert masked.returncode == 0, masked.stderr[-1500:]
    assert masked.stdout.startswith("mask ")


@needs_core
def test_detect_driver_echoes_a_source_without_deblend(tmp_path):
    rng = np.random.default_rng(4)
    img = (100.0 + rng.normal(0.0, 1.0, (64, 64))).astype(np.float32)
    _plant(img, 30, 28, 8000.0, 2.5)
    path = tmp_path / "sci.fits"
    _write_fits(path, img)
    out = tmp_path / "detect_catalog.tsv"
    got = _run("run_detect.py", [str(path), "--out", str(out), "--segmap", str(tmp_path / "seg.fits"),
                                 "--summary", str(tmp_path / "detect.json"),
                                 "--thresh", "3", "--minarea", "5", "--smooth-fwhm", "2.5",
                                 "--back-size", "32", "--mag-zeropoint", "25", "--max-sources", "100"])
    assert got.returncode == 0, got.stderr[-1500:]
    rows = [ln for ln in got.stdout.splitlines() if ln.strip()]
    assert len(rows) >= 2 and rows[0].startswith("NUMBER") and "X_IMAGE" in rows[0]
    assert json.loads((tmp_path / "detect.json").read_text())["n_sources"] >= 1


@needs_core
def test_catalog_drivers_extinct_period_dropout_and_lf(tmp_path):
    cat = tmp_path / "mags.tsv"
    cat.write_text("NUMBER\tMAG_ISO\n1\t20\n", encoding="utf-8")
    ext = _run("run_extinct.py", ["--ebv", "0.05", "--rv", "3.1", "--filter", "F160W",
                                  "--catalog", str(cat), "--mag-column", "MAG_ISO",
                                  "--out", str(tmp_path / "extinct.tsv"),
                                  "--summary", str(tmp_path / "extinct.json")])
    assert ext.returncode == 0, ext.stderr[-1500:]
    header = ext.stdout.splitlines()[0].split("\t")
    assert header[-2:] == ["MAG_ISO_DERED", "A_LAMBDA"]
    law = json.loads((tmp_path / "extinct.json").read_text())
    assert abs(law["a_lambda"] - 0.0313) < 0.0002

    rng = np.random.default_rng(5)
    t = np.sort(rng.uniform(0.0, 30.0, 80))
    y = 100.0 + 8.0 * np.sin(2.0 * math.pi * t / 1.7) + rng.normal(0.0, 0.2, t.size)
    series = tmp_path / "series.tsv"
    lines = ["JD_UTC\tFLUX"]
    lines.extend(f"{tt:.6f}\t{yy:.6f}" for tt, yy in zip(t, y))
    series.write_text("\n".join(lines) + "\n", encoding="utf-8")
    per = _run("run_period.py", [str(series), "--out", str(tmp_path / "periodogram.tsv"),
                                 "--summary", str(tmp_path / "period.json")])
    assert per.returncode == 0, per.stderr[-1500:]
    assert per.stdout.splitlines()[0] == "FREQUENCY_PER_DAY\tPERIOD_DAY\tPOWER"
    best = json.loads((tmp_path / "period.json").read_text())["best_period_day"]
    assert abs(best - 1.7) < 0.05

    faint = tmp_path / "drop.tsv"
    faint.write_text(
        "MAG_F435W\tMAGERR_F435W\tMAG_F606W\tMAGERR_F606W\tMAG_F775W\tMAGERR_F775W\tMAG_F850LP\tMAGERR_F850LP\n"
        "29.5\t2.0\t26.5\t0.05\t25.0\t0.05\t25.1\t0.05\n",
        encoding="utf-8",
    )
    drop = _run("run_dropout.py", [str(faint), "--out", str(tmp_path / "dropout.tsv"),
                                   "--summary", str(tmp_path / "dropout.json"), "--selection", "B"])
    assert drop.returncode == 0, drop.stderr[-1500:]
    assert json.loads((tmp_path / "dropout.json").read_text())["n_out"] == 0

    lf_cat = tmp_path / "lfcat.tsv"
    lf_cat.write_text("MAG_ABS\tVMAX\n-22\t1000\n-21\t1200\n-20\t1400\n-19\t1600\n", encoding="utf-8")
    lf = _run("run_lf.py", [str(lf_cat), "--out", str(tmp_path / "lf.tsv"), "--summary", str(tmp_path / "lf.json"),
                            "--z-min", "0", "--n-bins", "8", "--completeness", "1", "--h0", "70", "--omega-m", "0.3"])
    assert lf.returncode == 0, lf.stderr[-1500:]
    assert "sty skipped" in lf.stderr
    assert json.loads((tmp_path / "lf.json").read_text())["sty"]["status"] == "skipped"
    assert lf.stdout.splitlines()[0]


@needs_core
def test_counts_driver_runs_a_short_grid_and_drops_an_empty_catalog(tmp_path):
    rng = np.random.default_rng(2)
    img = (100.0 + rng.normal(0.0, 1.0, (48, 48))).astype(np.float32)
    path = tmp_path / "field.fits"
    _write_fits(path, img)
    common = [str(path), "--out", str(tmp_path / "counts.tsv"), "--summary", str(tmp_path / "counts.json"),
              "--kind", "star", "--n-bins", "2", "--per-bin", "2", "--per-image", "2",
              "--mag-zeropoint", "25", "--psf-fwhm", "2.5", "--match-radius", "3",
              "--thresh", "5", "--minarea", "5", "--smooth-fwhm", "2.5", "--back-size", "16", "--seed", "1"]
    got = _run("run_counts.py", common)
    assert got.returncode == 0, got.stderr[-1500:]
    body = json.loads((tmp_path / "counts.json").read_text())
    assert body["n_injected"] == 4 and body["deblend"] is False
    assert got.stdout.splitlines()[0]
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    again = _run("run_counts.py", [*common, "--catalog", str(empty), "--catalog-out", str(tmp_path / "ann.tsv")])
    assert again.returncode == 0, again.stderr[-1500:]
    assert not (tmp_path / "ann.tsv").exists()
