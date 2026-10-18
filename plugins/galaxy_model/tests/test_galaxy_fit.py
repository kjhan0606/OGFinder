"""Galaxy menu fits: elliptical is de Vaucouleurs, spiral is an exponential disc."""
import os
import subprocess
import sys

import numpy as np
from astropy.io import fits

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
SCRIPT = os.path.join(ROOT, "plugins", "galaxy_model", "galaxy_fit.py")
TCL = os.path.join(ROOT, "plugins", "galaxy_model", "galaxy_model.tcl")
sys.path.insert(0, ROOT)

from ogfkit import multifit as MF  # noqa: E402
from ogfkit import tsvio  # noqa: E402


def _env():
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    env["PYTHONPATH"] = ROOT
    return env


def _run(argv):
    return subprocess.run(argv, capture_output=True, text=True, cwd=ROOT, env=_env())


def _scene(tmp_path, kind, mag, re, q, pa, noise=1.0, seed=1):
    zp = 25.0
    truth = dict(kind=kind, x=40.2, y=39.7, mag=mag, re=re, q=q, pa=pa)
    image = MF.render_model([MF._normalise(truth, zp)], (81, 81), None, 40.0)
    image = image + np.random.default_rng(seed).normal(0.0, noise, image.shape)
    fits_path = tmp_path / ("%s.fits" % kind)
    fits.PrimaryHDU(image.astype(np.float32)).writeto(fits_path, overwrite=True)
    flux = 10 ** (-0.4 * (mag - zp))
    catalog = tmp_path / ("%s.tsv" % kind)
    tsvio.write_table(str(catalog), [
        "NUMBER", "X_IMAGE", "Y_IMAGE", "FLUX_AUTO", "A_IMAGE", "B_IMAGE",
        "THETA_IMAGE", "FLUX_RADIUS", "CLASS_STAR",
    ], [dict(
        NUMBER=7, X_IMAGE=41.2, Y_IMAGE=40.7, FLUX_AUTO=flux, A_IMAGE=8.0,
        B_IMAGE=8.0 * q, THETA_IMAGE=pa, FLUX_RADIUS=re, CLASS_STAR=0.05,
    )])
    return fits_path, catalog, truth


def _fit(tmp_path, kind, menu, mag=16.5, re=5.0, q=0.65, pa=30.0, seed=1):
    image, catalog, truth = _scene(tmp_path, kind, mag, re, q, pa, seed=seed)
    out = tmp_path / "fit.tsv"
    proc = _run([
        sys.executable, SCRIPT, str(image),
        "--catalog", str(catalog),
        "--model", menu,
        "--output", str(out),
        "--rms", "1.0",
        "--mag-zeropoint", "25",
        "--max-nfev", "120",
    ])
    assert proc.returncode == 0, proc.stderr
    _columns, rows = tsvio.read_catalog(str(out))
    assert len(rows) == 1
    return rows[0], truth, proc.stdout


def _delta_pa(got, truth):
    return ((float(got) - truth + 90.0) % 180.0) - 90.0


def test_elliptical_recovers_a_de_vaucouleurs_profile(tmp_path):
    row, truth, stdout = _fit(tmp_path, "dev", "elliptical", seed=3)
    assert "MODEL=dev" in stdout and "N=4.0" in stdout
    assert row["GM_MODEL"] == "dev"
    assert float(row["GM_N"]) == 4.0
    assert int(float(row["GM_FLAG"])) == 0
    dmag = float(row["GM_MAG"]) - truth["mag"]
    dre = float(row["GM_RE"]) / truth["re"] - 1.0
    dq = float(row["GM_Q"]) - truth["q"]
    dpa = _delta_pa(row["GM_PA"], truth["pa"])
    assert abs(dmag) < 0.10 and abs(dre) < 0.16 and abs(dq) < 0.05 and abs(dpa) < 8.0


def test_spiral_recovers_an_exponential_disc(tmp_path):
    row, truth, stdout = _fit(tmp_path, "exp", "spiral", mag=17.0, re=4.0, q=0.7, pa=20.0, seed=5)
    assert "MODEL=exp" in stdout and "N=1.0" in stdout
    assert row["GM_MODEL"] == "exp"
    assert float(row["GM_N"]) == 1.0
    assert int(float(row["GM_FLAG"])) == 0
    dmag = float(row["GM_MAG"]) - truth["mag"]
    dre = float(row["GM_RE"]) / truth["re"] - 1.0
    dq = float(row["GM_Q"]) - truth["q"]
    dpa = _delta_pa(row["GM_PA"], truth["pa"])
    assert abs(dmag) < 0.08 and abs(dre) < 0.10 and abs(dq) < 0.05 and abs(dpa) < 8.0


def test_stars_are_skipped_and_parameters_merge_into_the_catalog(tmp_path):
    image, catalog, _truth = _scene(tmp_path, "exp", 17.0, 4.0, 0.8, 10.0, seed=2)
    _columns, rows = tsvio.read_catalog(str(catalog))
    rows.append(dict(
        NUMBER=9, X_IMAGE=20.0, Y_IMAGE=20.0, FLUX_AUTO=1.0e6, A_IMAGE=2.0,
        B_IMAGE=2.0, THETA_IMAGE=0.0, FLUX_RADIUS=1.5, CLASS_STAR=0.95,
    ))
    tsvio.write_table(str(catalog), _columns, rows)
    out = tmp_path / "fit.tsv"
    proc = _run([
        sys.executable, SCRIPT, str(image), "--catalog", str(catalog),
        "--model", "spiral", "--output", str(out), "--rms", "1.0", "--max-nfev", "80",
    ])
    assert proc.returncode == 0, proc.stderr
    assert "N_FIT=1" in proc.stdout
    assert "N_SOURCES=2" in proc.stdout
    _cols, fitted = tsvio.read_catalog(str(out))
    assert [row["NUMBER"] for row in fitted] == ["7"]

    merged = tmp_path / "catalog.tsv"
    applied = _run([
        sys.executable, SCRIPT, "--apply-catalog", str(catalog),
        "--results", str(out), "--output", str(merged),
    ])
    assert applied.returncode == 0, applied.stderr
    columns, body = tsvio.read_catalog(str(merged))
    for name in ("GM_MODEL", "GM_MAG", "GM_RE", "GM_N", "GM_Q", "GM_PA", "GM_CHI2"):
        assert name in columns
    by_number = {row["NUMBER"]: row for row in body}
    assert by_number["7"]["GM_MODEL"] == "exp"
    assert float(by_number["7"]["GM_N"]) == 1.0
    assert by_number["7"]["X_IMAGE"] == "41.2"
    assert by_number["9"]["GM_MAG"] == ""
    assert by_number["9"]["X_IMAGE"] == "20"


def test_other_models_write_no_file(tmp_path):
    image, catalog, _truth = _scene(tmp_path, "exp", 17.0, 4.0, 0.8, 10.0, seed=1)
    out = tmp_path / "nope.tsv"
    proc = _run([
        sys.executable, SCRIPT, str(image), "--catalog", str(catalog),
        "--model", "bulge+disk", "--output", str(out),
    ])
    assert proc.returncode == 2
    assert "elliptical or spiral" in proc.stderr
    assert not out.exists()
    missing = _run([
        sys.executable, SCRIPT, str(image), "--catalog", str(catalog),
        "--model", "elliptical", "--psf", str(tmp_path / "absent.fits"),
        "--output", str(out),
    ])
    assert missing.returncode == 1
    assert "not found" in missing.stderr
    assert not out.exists()


def test_menu_calls_the_two_presets_and_extracts_columns():
    import importlib.util
    spec = importlib.util.spec_from_file_location("galaxy_fit_menu", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.PRESETS == {"elliptical": "dev", "dev": "dev", "spiral": "exp", "exp": "exp"}
    assert set(module.PRESETS.values()) == {"dev", "exp"}
    text = open(TCL, encoding="utf-8").read()
    assert "not yet implemented" not in text
    assert "galaxy_fit.py" in text
    assert "--model $model" in text
    assert "de Vaucouleurs n=4" in text
    assert "exponential disc n=1" in text
    assert "--apply-catalog" in text
    assert "Extract Parameters writes the catalog columns." in text
