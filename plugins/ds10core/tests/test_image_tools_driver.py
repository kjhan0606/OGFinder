"""Desktop drivers for the shared-core saturation mask, varying PSF, photometric solution and star/galaxy step.

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


def _env():
    e = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        e[name] = "2"
    return e


def _run(script, args):
    return subprocess.run([sys.executable, os.path.join(PLUGIN, script), "--core-dir", CORE or "", *args],
                          capture_output=True, text=True, env=_env(), timeout=180)


def _write_fits(path, data, **cards):
    from astropy.io import fits
    hdu = fits.PrimaryHDU(np.asarray(data, np.float32))
    hdu.header.update(cards)
    hdu.writeto(path, overwrite=True)


def _plant(img, x, y, flux, fwhm):
    sig = fwhm / _FWHM
    r = int(math.ceil(5.0 * sig)) + 2
    ny, nx = img.shape
    x0, x1 = max(0, int(math.floor(x)) - r), min(nx, int(math.floor(x)) + r + 1)
    y0, y1 = max(0, int(math.floor(y)) - r), min(ny, int(math.floor(y)) + r + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    img[y0:y1, x0:x1] += flux * np.exp(-0.5 * ((xx - x) ** 2 + (yy - y) ** 2) / (sig * sig))


def _tsv(text):
    lines = [ln for ln in text.splitlines() if ln.strip()]
    names = lines[0].split("\t")
    return [dict(zip(names, ln.split("\t"))) for ln in lines[1:]]


def test_manifest_steps_stay_on_the_process_boundary():
    m = json.load(open(os.path.join(PLUGIN, "plugin.json"), encoding="utf-8"))
    assert m["web"] is False and "separate process" in m["license"]["boundary"]
    steps = {s["id"]: s for s in m["steps"]}
    assert "mask" not in steps
    params = {p["name"]: p for p in m["params"]}
    expect = {
        "pixmask": ("run_pixmask.py", ["image"], "text"),
        "psfvar": ("run_psfvar.py", ["image"], "text"),
        "photcal": ("run_photcal.py", ["catalog"], "set"),
        "stargal": ("run_stargal.py", ["image", "catalog"], "set"),
    }
    for sid, (drv, needs, mode) in expect.items():
        step = steps[sid]
        visible = step["label"] + step["title"] + step["done_status"]
        assert "ds10" not in visible, sid
        assert step["record"] == "ds10core." + sid
        assert step["needs"] == needs and step["output"]["mode"] == mode
        src = open(os.path.join(PLUGIN, drv), encoding="utf-8").read()
        assert "import ds10core" not in src and "from ds10core" not in src
        flat = []
        for item in step["cli"]:
            if isinstance(item, str):
                flat.append(item)
            else:
                flat.extend(item["argv"])
                for flag in item["argv"]:
                    if flag.startswith("--"):
                        assert flag in src, (sid, flag)
        for flag in flat:
            if isinstance(flag, str) and flag.startswith("--"):
                assert flag in src, (sid, flag)
        assert drv in " ".join(item if isinstance(item, str) else "" for item in step["cli"])
    pix = steps["pixmask"]["cli"]
    assert any(isinstance(e, dict) and e.get("if") == "pixmask-bleed" for e in pix)
    assert any(isinstance(e, dict) and e.get("if") == "pixmask-saturate" for e in pix)
    assert "--saturate" not in [e for e in pix if isinstance(e, str)]
    assert params["pixmask-bleed"]["default"] is True and params["pixmask-dq-bits"]["default"] == 0
    assert params["psfvar-degree"]["default"] == "1"
    assert params["photcal-std-column"]["default"] == "MAG_STD" and params["photcal-airmass"]["default"] == ""
    assert params["stargal-psf-fwhm"]["default"] == "" and params["stargal-snr-min"]["default"] == 8


@needs_core
def test_pixmask_driver_marks_saturation_and_bleed(tmp_path):
    rng = np.random.default_rng(1)
    img = (100.0 + rng.normal(0.0, 1.0, (48, 40))).astype(np.float32)
    img[12:29, 20] = 40000
    img[29:37, 20] = 800
    path = tmp_path / "sci.fits"
    _write_fits(path, img, SATURATE=30000.0)
    out = tmp_path / "pixmask.fits"
    summary = tmp_path / "pixmask_summary.json"
    r = _run("run_pixmask.py", [str(path), "--out", str(out), "--summary", str(summary),
                                "--bleed", "--bleed-axis", "columns", "--bleed-nsig", "10", "--bleed-max", "40",
                                "--dq-bits", "0", "--saturate", "30000"])
    assert r.returncode == 0, r.stderr[-1500:]
    assert r.stdout.startswith("pixmask sat ")
    body = json.loads(summary.read_text())
    assert body["n_saturation"] >= 6 and body["n_bleed"] >= 1 and body["n_persistence"] == 0
    from astropy.io import fits
    mask = fits.getdata(out)
    assert int(mask[20, 20]) & 1 and int(mask[32, 20]) & 2
    missing = _run("run_pixmask.py", [str(tmp_path / "absent.fits"), "--out", str(out), "--summary", str(summary)])
    assert missing.returncode == 2


@needs_core
def test_psfvar_driver_recovers_the_gradient_and_ignores_an_empty_catalog(tmp_path):
    ny = nx = 200
    rng = np.random.default_rng(7)
    img = 100.0 + rng.normal(0.0, 0.5, (ny, nx))
    lines = ["X_IMAGE\tY_IMAGE"]
    xs = np.linspace(28, 172, 5)
    ys = np.linspace(28, 172, 5)
    for y in ys:
        for x in xs:
            _plant(img, float(x), float(y), 8000.0, 2.2 + 1.6 * (float(x) / (nx - 1)))
            lines.append(f"{float(x) + 1:.3f}\t{float(y) + 1:.3f}")
    path = tmp_path / "field.fits"
    _write_fits(path, img.astype(np.float32))
    cat = tmp_path / "stars.tsv"
    cat.write_text("\n".join(lines) + "\n", encoding="utf-8")
    model = tmp_path / "psf_model.json"
    stamp = tmp_path / "psf_center.fits"
    summary = tmp_path / "psfvar_summary.json"
    common = [str(path), "--out-model", str(model), "--out-stamp", str(stamp), "--summary", str(summary),
              "--degree", "1", "--stamp", "21"]
    r = _run("run_psfvar.py", [*common, "--catalog", str(cat)])
    assert r.returncode == 0, r.stderr[-1500:]
    assert r.stdout.startswith("psfvar ") and "degree 1" in r.stdout
    body = json.loads(model.read_text())
    assert body["n_stars"] >= 20 and body["coefficients"]["fwhm"][1] > 0.4
    assert 2.4 < body["center"]["fwhm"] < 3.6
    from astropy.io import fits
    assert abs(float(fits.getdata(stamp).sum()) - 1.0) < 1e-4
    empty = tmp_path / "empty.tsv"
    empty.write_text("", encoding="utf-8")
    r2 = _run("run_psfvar.py", [*common, "--catalog", str(empty)])
    assert r2.returncode == 0, r2.stderr[-1500:]
    detected = json.loads(model.read_text())
    assert detected["degree"] == 1 and detected["n_stars"] >= 6
    assert r2.stdout.startswith("psfvar ")


@needs_core
def test_photcal_driver_replaces_the_table(tmp_path):
    rng = np.random.default_rng(4)
    n = 24
    inst = rng.uniform(12.0, 16.0, n)
    air = rng.uniform(1.05, 1.9, n)
    color = rng.uniform(-0.2, 1.4, n)
    std = inst + 3.25 + 0.18 * air - 0.07 * color + rng.normal(0.0, 0.005, n)
    lines = ["MAG_ISO\tMAG_STD\tAIRMASS\tCOLOR"]
    for i in range(n):
        lines.append(f"{inst[i]:.5f}\t{std[i]:.5f}\t{air[i]:.5f}\t{color[i]:.5f}")
    lines.append("14.00000\t\t1.40000\t0.30000")
    cat = tmp_path / "std.tsv"
    cat.write_text("\n".join(lines) + "\n", encoding="utf-8")
    out = tmp_path / "photcal.tsv"
    summary = tmp_path / "photcal_summary.json"
    r = _run("run_photcal.py", [str(cat), "--out", str(out), "--summary", str(summary),
                                "--std-column", "MAG_STD", "--clip", "3"])
    assert r.returncode == 0, r.stderr[-1500:]
    assert r.stdout.splitlines()[0].split("\t")[0] == "MAG_ISO"
    assert "MAG_CAL" in r.stdout.splitlines()[0]
    assert "photcal zp" in r.stderr
    body = json.loads(summary.read_text())
    assert abs(body["zp"] - 3.25) < 0.03 and abs(body["k"] - 0.18) < 0.03 and abs(body["c"] - (-0.07)) < 0.03
    bad = tmp_path / "bad.tsv"
    bad.write_text("MAG_ISO\n18\n", encoding="utf-8")
    fail = _run("run_photcal.py", [str(bad), "--out", str(out), "--summary", str(summary)])
    assert fail.returncode == 2 and fail.stdout == ""


@needs_core
def test_stargal_driver_separates_stars_and_galaxies(tmp_path):
    rng = np.random.default_rng(9)
    img = 40.0 + rng.normal(0.0, 1.0, (240, 240))
    rows = ["NUMBER\tX_IMAGE\tY_IMAGE\tSNR\tCLASS"]
    n = 1
    for x in (40.0, 90.0, 140.0):
        for y in (40.0, 90.0, 140.0):
            _plant(img, x, y, 5000.0, 3.0)
            rows.append(f"{n}\t{x + 1:.3f}\t{y + 1:.3f}\t40.0\tcompact")
            n += 1
    for x in (50.0, 120.0, 190.0):
        _plant(img, x, 200.0, 9000.0, 10.0)
        rows.append(f"{n}\t{x + 1:.3f}\t201.000\t30.0\tcompact")
        n += 1
    rows.append(f"{n}\t40.000\t80.000\t3.0\tcompact")
    path = tmp_path / "field.fits"
    _write_fits(path, img.astype(np.float32))
    cat = tmp_path / "cat.tsv"
    cat.write_text("\n".join(rows) + "\n", encoding="utf-8")
    out = tmp_path / "stargal.tsv"
    summary = tmp_path / "stargal_summary.json"
    r = _run("run_stargal.py", [str(path), "--catalog", str(cat), "--out", str(out), "--summary", str(summary),
                                "--snr-min", "8", "--k-sigma", "3", "--elong-max", "1.6", "--psf-fwhm", "3"])
    assert r.returncode == 0, r.stderr[-1500:]
    parsed = _tsv(r.stdout)
    kinds = [row["CLASS_SG"] for row in parsed]
    assert kinds[:9] == ["star"] * 9
    assert kinds[9:12] == ["galaxy"] * 3
    assert kinds[12] == "unknown"
    assert all(row["CLASS"] == "compact" for row in parsed)
    body = json.loads(summary.read_text())
    assert body["n_star"] == 9 and body["n_galaxy"] == 3 and body["psf_source"] == "user"
