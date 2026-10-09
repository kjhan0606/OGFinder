"""Single elliptical Sérsic, with an optional PSF. Gates were fixed before the run."""
import math
import os
import re
import subprocess
import sys

import numpy as np

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

from astropy.io import fits

from ogfmeas.sersic import (
    convolve, fit_sersic, galfit_pa_from_theta, render, theta_from_galfit_pa, total_flux,
)

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
VAL = os.path.join(ROOT, "plugins", "multifit", "validation")
if VAL not in sys.path:
    sys.path.insert(0, VAL)

# Planted noiseless galaxy. These limits were chosen before fitting it.
X_TRUE, Y_TRUE = 30.0, 28.0
RE_TRUE, N_TRUE, IE_TRUE = 5.5, 2.0, 50.0
Q_TRUE, TH_TRUE, BG_TRUE = 0.6, 0.8, 2.0
SHAPE = (71, 71)


def _theta_delta(a, b):
    return abs((a - b + math.pi / 2.0) % math.pi - math.pi / 2.0)


def _gaussian_psf(fwhm=2.5, size=15):
    sigma = fwhm / (2.0 * math.sqrt(2.0 * math.log(2.0)))
    yy, xx = np.mgrid[:size, :size] - size // 2
    kernel = np.exp(-0.5 * (xx * xx + yy * yy) / sigma ** 2)
    return kernel / kernel.sum()


def _mag(ie, re, n, q, zp=25.0):
    return zp - 2.5 * math.log10(total_flux(ie, re, n, q))


def _feedme(x, y, mag, re, n, q, pa, sky, shape):
    ny, nx = shape
    lines = [
        "A) data.fits", "B) out.fits", "C) none", "D) psf.fits", "E) 1", "F) none", "G) none",
        "H) 1 %d 1 %d" % (nx, ny), "I) 31 31", "J) 25.000", "K) 0.1 0.1", "O) regular", "P) 0", "",
        " 0) sersic",
        " 1) %.4f %.4f 1 1" % (x, y),
        " 3) %.4f 1" % mag,
        " 4) %.4f 1" % re,
        " 5) %.4f 1" % n,
        " 9) %.4f 1" % q,
        "10) %.4f 1" % pa,
        " Z) 0",
        " 0) sky",
        " 1) %.4f 1" % sky,
        " 2) 0.000000 0",
        " 3) 0.000000 0",
        " Z) 0", "",
    ]
    return "\n".join(lines)


def _assert_recovered(fit):
    assert abs(fit["n"] - N_TRUE) < 0.35, fit
    assert abs(fit["re"] - RE_TRUE) < 0.8, fit
    assert abs(fit["q"] - Q_TRUE) < 0.08, fit
    assert abs(fit["x"] - X_TRUE) < 0.4, fit
    assert abs(fit["y"] - Y_TRUE) < 0.4, fit
    assert _theta_delta(fit["theta"], TH_TRUE) < 0.20, fit
    assert abs(fit["background"] - BG_TRUE) < 0.5, fit
    assert abs(_mag(fit["ie"], fit["re"], fit["n"], fit["q"]) - _mag(IE_TRUE, RE_TRUE, N_TRUE, Q_TRUE)) < 0.20, fit


def test_total_flux_matches_the_pixel_sum():
    image = render((201, 201), 100.0, 100.0, re=4.0, n=1.0, ie=10.0, q=0.5, theta=0.3, background=0.0)
    analytic = total_flux(10.0, 4.0, 1.0, 0.5)
    assert abs(float(image.sum()) / analytic - 1.0) < 0.01


def test_a_delta_psf_matches_the_unconvolved_image():
    psf = np.zeros((5, 5))
    psf[2, 2] = 4.0
    plain = render(SHAPE, X_TRUE, Y_TRUE, RE_TRUE, N_TRUE, IE_TRUE, Q_TRUE, TH_TRUE, BG_TRUE)
    convolved = render(SHAPE, X_TRUE, Y_TRUE, RE_TRUE, N_TRUE, IE_TRUE, Q_TRUE, TH_TRUE, BG_TRUE, psf=psf)
    assert np.max(np.abs(plain - convolved)) < 1e-8


def test_position_angle_round_trip():
    for pa in (-89.0, -80.0, -10.0, 0.0, 25.0, 89.0):
        back = galfit_pa_from_theta(theta_from_galfit_pa(pa))
        assert abs((back - pa + 90.0) % 180.0 - 90.0) < 1e-6


def test_convolved_sersic_is_recovered_from_a_perturbed_start():
    psf = _gaussian_psf()
    image = render(SHAPE, X_TRUE, Y_TRUE, RE_TRUE, N_TRUE, IE_TRUE, Q_TRUE, TH_TRUE, BG_TRUE, psf=psf)
    fit = fit_sersic(
        image, x0=X_TRUE + 0.8, y0=Y_TRUE - 0.6, psf=psf,
        re0=4.4, n0=1.4, q0=0.75, theta0=TH_TRUE - 0.25,
        ie0=0.7 * IE_TRUE, background0=BG_TRUE + 0.5,
    )
    _assert_recovered(fit)
    again = render(SHAPE, fit["x"], fit["y"], fit["re"], fit["n"], fit["ie"], fit["q"], fit["theta"], fit["background"], psf=psf)
    assert float(np.sqrt(np.mean((again - image) ** 2))) < 0.30


def _plant(directory):
    os.makedirs(directory, exist_ok=True)
    psf = _gaussian_psf()
    image = render(SHAPE, X_TRUE, Y_TRUE, RE_TRUE, N_TRUE, IE_TRUE, Q_TRUE, TH_TRUE, BG_TRUE, psf=psf)
    fits.PrimaryHDU(image.astype(np.float32)).writeto(os.path.join(directory, "data.fits"), overwrite=True)
    fits.PrimaryHDU(psf.astype(np.float32)).writeto(os.path.join(directory, "psf.fits"), overwrite=True)
    fits.PrimaryHDU(np.ones(SHAPE, np.float32)).writeto(os.path.join(directory, "sigma.fits"), overwrite=True)
    start = _feedme(
        X_TRUE + 1.8, Y_TRUE + 0.4, _mag(0.7 * IE_TRUE, 4.4, 1.4, 0.75),
        4.4, 1.4, 0.75, galfit_pa_from_theta(TH_TRUE - 0.25), BG_TRUE + 0.5, SHAPE,
    )
    truth = _feedme(
        X_TRUE + 1.0, Y_TRUE + 1.0, _mag(IE_TRUE, RE_TRUE, N_TRUE, Q_TRUE),
        RE_TRUE, N_TRUE, Q_TRUE, galfit_pa_from_theta(TH_TRUE), BG_TRUE, SHAPE,
    )
    open(os.path.join(directory, "start.feedme"), "w").write(start)
    open(os.path.join(directory, "truth.feedme"), "w").write(truth)
    return image


def test_directory_fit_and_compare_do_not_need_an_external_solution(tmp_path):
    import compare_pysersic
    import pysersic_map as PM

    work = tmp_path / "work"
    _plant(str(work / "A00"))
    out = tmp_path / "fit.json"
    assert PM.main([str(work), str(out), "--n", "1"]) == 0
    rows = __import__("json").load(open(out))
    assert rows[0]["fitter"] == "ogfmeas.sersic"
    assert "error" not in rows[0]
    assert abs(rows[0]["x"] - (X_TRUE + 1.0)) < 0.4
    assert abs(rows[0]["y"] - (Y_TRUE + 1.0)) < 0.4
    assert abs(rows[0]["n"] - N_TRUE) < 0.35
    assert abs(rows[0]["re"] - RE_TRUE) < 0.8
    assert abs(rows[0]["q"] - Q_TRUE) < 0.08
    cmp_path = tmp_path / "cmp.json"
    __import__("json").dump({"results": [{
        "kind": "A", "seed": 0,
        "truth": [{"type": "sersic", "x": X_TRUE + 1.0, "y": Y_TRUE + 1.0,
                   "mag": _mag(IE_TRUE, RE_TRUE, N_TRUE, Q_TRUE), "re": RE_TRUE,
                   "n": N_TRUE, "q": Q_TRUE, "pa": galfit_pa_from_theta(TH_TRUE)}],
        "truth_sky": BG_TRUE,
    }]}, open(cmp_path, "w"))
    scored = tmp_path / "scored.json"
    compare_pysersic.compare(str(cmp_path), str(out), str(work), str(scored))
    payload = __import__("json").load(open(scored))
    assert payload["summary"]["n"] == 1
    assert payload["summary"]["fitter"] == "ogfmeas.sersic"
    assert "galfit" not in payload["rows"][0]
    assert payload["rows"][0]["ogfmeas_chi2"] < 500.0
    assert abs(payload["rows"][0]["ogfmeas"]["n_rel"]) < 0.35 / N_TRUE


def test_run_one_fits_sersic_without_the_external_binary(tmp_path, monkeypatch):
    import galfit_real_compare as GRC

    def boom(*_args, **_kwargs):
        raise AssertionError("external binary was called")

    monkeypatch.setattr(GRC.GC, "run_galfit", boom)

    def fake_run(argv, **_kwargs):
        assert os.path.basename(str(argv[0])) != "galfit"
        class Result:
            returncode = 1
            stdout = ""
            stderr = "multifit not run in this test"
        return Result()

    monkeypatch.setattr(GRC.subprocess, "run", fake_run)
    root = tmp_path / "real"
    _plant(str(root / "cut"))
    task = (dict(model="sersic", dir="cut", sigma=1.0, band="t", k=0), __import__("types").SimpleNamespace(
        work=str(root), timeout=30, max_nfev=5, galfit="/does/not/exist", ld=""))
    result = GRC.run_one(task)
    assert result["galfit"] is None
    assert "ogfmeas.sersic" in result["galfit_err"]
    fit = result["ogfmeas"]
    assert fit["fitter"] == "ogfmeas.sersic"
    assert abs(fit["n"] - N_TRUE) < 0.35, fit
    assert abs(fit["re"] - RE_TRUE) < 0.8, fit
    assert abs(fit["q"] - Q_TRUE) < 0.08, fit
    assert isinstance(fit["chi2"], tuple)
    assert fit["chi2"][0] < 500.0
    assert result["multifit"] is None


def test_devexp_is_not_sent_to_the_external_binary(tmp_path, monkeypatch):
    import galfit_real_compare as GRC

    def boom(*_args, **_kwargs):
        raise AssertionError("external binary was called")

    monkeypatch.setattr(GRC.GC, "run_galfit", boom)
    called = []

    def fake_run(argv, **_kwargs):
        called.append(argv)
        assert os.path.basename(str(argv[0])) != "galfit"
        class Result:
            returncode = 1
            stdout = ""
            stderr = "multifit not run in this test"
        return Result()

    monkeypatch.setattr(GRC.subprocess, "run", fake_run)
    root = tmp_path / "real"
    _plant(str(root / "cut"))
    task = (dict(model="devexp", dir="cut", sigma=1.0, band="t", k=0), __import__("types").SimpleNamespace(
        work=str(root), timeout=30, max_nfev=5, galfit="galfit", ld=""))
    result = GRC.run_one(task)
    assert result["galfit"] is None
    assert result["ogfmeas"] is None
    assert "out of scope" in result["galfit_err"]
    assert called


def test_render_chi2_uses_the_in_tree_model(tmp_path, monkeypatch):
    import galfit_real_compare as GRC

    def boom(*_args, **_kwargs):
        raise AssertionError("external binary was called")

    monkeypatch.setattr(GRC.GC, "run_galfit", boom)
    directory = tmp_path / "cut"
    _plant(str(directory))
    text = open(directory / "truth.feedme").read()
    scored = GRC.render_chi2(None, str(directory), text, "t", 1.0)
    assert scored is not None
    chi2, npix = scored
    assert npix == SHAPE[0] * SHAPE[1]
    assert math.isfinite(chi2)


def test_validation_modules_do_not_load_pysersic_or_jax():
    for name in ("pysersic_map.py", "compare_pysersic.py"):
        text = open(os.path.join(VAL, name)).read()
        assert re.search(r"(?m)^\s*(?:import|from)\s+pysersic(?:\s|\.|$)", text) is None
        assert re.search(r"(?m)^\s*(?:import|from)\s+jax(?:\s|\.|$)", text) is None
        assert "subprocess" not in text
        assert "run_galfit" not in text
    text = open(os.path.join(VAL, "galfit_real_compare.py")).read()
    assert "run_galfit(" not in text
    script = (
        "import importlib.util, sys\n"
        "names = [%r, %r]\n"
        "for path, modname in zip(names, ['pysersic_map_check', 'compare_pysersic_check']):\n"
        "    spec = importlib.util.spec_from_file_location(modname, path)\n"
        "    mod = importlib.util.module_from_spec(spec)\n"
        "    spec.loader.exec_module(mod)\n"
        "bad = [name for name in sys.modules if name == 'pysersic' or name.startswith('pysersic.') or name == 'jax' or name.startswith('jax.')]\n"
        "if bad:\n"
        "    raise SystemExit('loaded ' + ','.join(bad))\n"
        "print('clean')\n"
    ) % (os.path.join(VAL, "pysersic_map.py"), os.path.join(VAL, "compare_pysersic.py"))
    done = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, cwd=ROOT)
    assert done.returncode == 0, done.stderr
    assert "clean" in done.stdout


def test_convolve_rejects_an_empty_kernel():
    try:
        convolve(np.zeros((4, 4)), np.zeros((3, 3)))
    except ValueError as exc:
        assert "positive" in str(exc)
    else:
        raise AssertionError("empty PSF was accepted")
