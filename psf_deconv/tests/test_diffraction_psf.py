"""In-tree circular diffraction PSF.

The model is an annular Airy pattern plus Gaussian jitter. Its method name
is ogfmeas. WebbPSF and TinyTim stay on their own explicit paths.
"""
import json
import os
import subprocess
import sys

import numpy as np
import pytest
from astropy.io import fits
from scipy.special import jn_zeros

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
CLI = os.path.join(ROOT, "ds9", "library", "ds9_psf_deconv.py")
CLEAN_PATH = "/usr/bin:/bin"

sys.path.insert(0, ROOT)
from psf_deconv.psf.sim_psf import (  # noqa: E402
    _ARCSEC_PER_RAD,
    generate_diffraction,
    resolve_diffraction,
)


def _env():
    env = dict(os.environ)
    env["PATH"] = CLEAN_PATH
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    return env


def _run(argv):
    return subprocess.run(
        argv, capture_output=True, text=True, cwd=ROOT, env=_env(),
    )


def _airy_radius_px(diameter_m, wavelength_um, pixel_scale):
    rho = float(jn_zeros(1, 1)[0])
    theta_arcsec = (rho / np.pi) * (wavelength_um * 1e-6 / diameter_m) * _ARCSEC_PER_RAD
    return theta_arcsec / pixel_scale


def test_unobscured_airy_zero_is_the_first_dark_pixel():
    diameter_m = 2.4
    wavelength_um = 0.5
    pixel_scale = 0.01
    radius = _airy_radius_px(diameter_m, wavelength_um, pixel_scale)
    assert 5.2 < radius < 5.3

    image = generate_diffraction(
        psf_size=65,
        diameter_m=diameter_m,
        obscuration=0.0,
        wavelength_um=wavelength_um,
        pixel_scale=pixel_scale,
    )
    assert image.shape == (65, 65)
    assert abs(float(image.sum()) - 1.0) < 1e-12
    assert np.isfinite(image).all()
    assert (image >= 0).all()
    centre = 32
    assert np.unravel_index(int(np.argmax(image)), image.shape) == (centre, centre)
    cut = image[centre]
    dark = min(range(3, 8), key=lambda k: cut[centre + k])
    assert dark == 5
    assert cut[centre + 5] < cut[centre + 4]
    assert cut[centre + 5] < cut[centre + 6]
    assert abs(float(cut[centre - 5] - cut[centre + 5])) < 1e-15


def test_obscuration_and_jitter_lower_the_peak_and_oversample_keeps_the_shape():
    common = dict(psf_size=65, diameter_m=2.4, wavelength_um=0.5, pixel_scale=0.01)
    clear = generate_diffraction(obscuration=0.0, jitter_sigma=0.0, **common)
    obscured = generate_diffraction(obscuration=0.33, jitter_sigma=0.0, **common)
    blurred = generate_diffraction(obscuration=0.0, jitter_sigma=0.03, **common)
    assert float(obscured.max()) < float(clear.max())
    assert float(blurred.max()) < float(clear.max())

    binned = generate_diffraction(psf_size=32, oversample=2, obscuration=0.0, **{
        "diameter_m": 2.4, "wavelength_um": 0.5, "pixel_scale": 0.01,
    })
    assert binned.shape == (32, 32)
    assert abs(float(binned.sum()) - 1.0) < 1e-12
    centre = 15
    assert np.unravel_index(int(np.argmax(binned)), binned.shape) == (centre, centre)


def test_known_instruments_fill_the_pupil_table():
    nircam = resolve_diffraction("NIRCAM", "F150W")
    assert nircam["telescope"] == "JWST"
    assert nircam["diameter_m"] == 6.5
    assert nircam["obscuration"] == pytest.approx(0.114)
    assert nircam["wavelength_um"] == pytest.approx(1.50)
    assert nircam["wavekind"] == "pivot"
    assert nircam["pixel_scale"] == pytest.approx(0.031)

    long = resolve_diffraction("nircam", "F444W")
    assert long["instrument"] == "NIRCAM"
    assert long["wavelength_um"] == pytest.approx(4.44)
    assert long["pixel_scale"] == pytest.approx(0.063)

    fgs = resolve_diffraction("FGS", None)
    assert fgs["wavelength_um"] == pytest.approx(0.65)
    assert fgs["pixel_scale"] == pytest.approx(0.069)
    assert fgs["wavekind"] == "pivot"

    acs = resolve_diffraction("ACS", "F814W")
    assert acs["instrument"] == "ACS_WFC"
    assert acs["telescope"] == "HST"
    assert acs["diameter_m"] == 2.4
    assert acs["obscuration"] == pytest.approx(0.33)
    assert acs["wavelength_um"] == pytest.approx(0.814)
    assert acs["pixel_scale"] == pytest.approx(0.05)
    assert acs["wavekind"] == "pivot"

    lp = resolve_diffraction("ACS_WFC", "F850LP")
    assert lp["wavelength_um"] == pytest.approx(0.850)
    assert lp["wavekind"] == "cutoff"

    given = resolve_diffraction("ACS_WFC", "F814W", wavelength_um=0.6, pixel_scale=0.04)
    assert given["wavelength_um"] == pytest.approx(0.6)
    assert given["pixel_scale"] == pytest.approx(0.04)
    assert given["wavekind"] == "given"


def test_missing_numbers_raise_before_a_file_is_written(tmp_path):
    out = tmp_path / "missing.fits"
    with pytest.raises(ValueError, match="wavelength was not set"):
        generate_diffraction(instrument="WFPC2", filter_name="F555W", output=str(out))
    assert not out.exists()
    with pytest.raises(ValueError, match="diameter was not set"):
        generate_diffraction(wavelength_um=0.5, pixel_scale=0.01, output=str(out))
    with pytest.raises(ValueError, match="pixel scale was not set"):
        generate_diffraction(diameter_m=2.4, wavelength_um=0.5, output=str(out))
    with pytest.raises(ValueError, match=r"obscuration must be in \[0, 1\)"):
        generate_diffraction(
            diameter_m=2.4, wavelength_um=0.5, pixel_scale=0.01, obscuration=1.0,
            output=str(out))
    with pytest.raises(ValueError, match="must be positive"):
        generate_diffraction(
            diameter_m=-1.0, wavelength_um=0.5, pixel_scale=0.01, output=str(out))
    assert not out.exists()


def test_cli_writes_an_ogfmeas_header_and_refuses_a_focus_offset(tmp_path):
    status = _run([sys.executable, CLI, "dummy.fits", "--mode", "check_sim"])
    assert status.returncode == 0, status.stderr
    assert "#SIM_STATUS\tWEBBPSF=0\tTINYTIM=0\tDIFFRACT=1" in status.stdout
    assert "Diffraction: available" in status.stderr
    assert "WebbPSF: Not found" in status.stderr
    assert "not on PATH" in status.stderr

    out = tmp_path / "acs.fits"
    made = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "ogfmeas",
        "--sim-instrument", "ACS_WFC",
        "--sim-filter", "F814W",
        "--sim-psf-size", "31",
        "--psf-output", str(out),
    ])
    assert made.returncode == 0, made.stderr
    assert f"#PSF_SIM\tMODEL=ogfmeas\tSIZE=31\tOUTPUT={out}" in made.stdout
    assert "webbpsf" not in made.stdout.lower()
    assert "tinytim" not in made.stdout.lower()
    header = fits.getheader(str(out))
    assert header["PSFMETH"] == "ogfmeas"
    assert header["PUPIL"] == "circular"
    assert header["TELESCOP"] == "HST"
    assert header["INSTRUME"] == "ACS_WFC"
    assert header["FILTER"] == "F814W"
    assert header["WAVE_UM"] == pytest.approx(0.814)
    assert header["WAVEKIND"] == "pivot"
    assert header["DIAM_M"] == pytest.approx(2.4)
    assert header["OBSCFRAC"] == pytest.approx(0.33)
    assert header["PIXSCALE"] == pytest.approx(0.05)
    data = fits.getdata(str(out))
    assert data.shape == (31, 31)
    assert abs(float(data.sum()) - 1.0) < 1e-5

    bare = tmp_path / "bare.fits"
    explicit = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "ogfmeas",
        "--sim-instrument", "CUSTOM",
        "--sim-filter", "NONE",
        "--sim-diameter", "1.5",
        "--sim-obscuration", "0.2",
        "--sim-wavelength", "0.7",
        "--sim-pixel-scale", "0.02",
        "--sim-psf-size", "15",
        "--psf-output", str(bare),
    ])
    assert explicit.returncode == 0, explicit.stderr
    bare_header = fits.getheader(str(bare))
    assert bare_header["PSFMETH"] == "ogfmeas"
    assert bare_header["WAVEKIND"] == "given"
    assert bare_header["WAVE_UM"] == pytest.approx(0.7)
    assert "TELESCOP" not in bare_header

    missing = tmp_path / "missing.fits"
    refused = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "ogfmeas",
        "--sim-instrument", "WFPC2",
        "--sim-filter", "F555W",
        "--psf-output", str(missing),
    ])
    assert refused.returncode == 1
    assert "wavelength was not set" in refused.stderr
    assert not missing.exists()

    focused = tmp_path / "focus.fits"
    focus = _run([
        sys.executable, CLI, "dummy.fits",
        "--mode", "sim_psf",
        "--sim-telescope", "ogfmeas",
        "--sim-instrument", "ACS_WFC",
        "--sim-filter", "F814W",
        "--sim-focus-offset", "0.2",
        "--psf-output", str(focused),
    ])
    assert focus.returncode == 1
    assert "focus offset is not part of the circular diffraction model." in focus.stderr
    assert not focused.exists()


def test_menu_points_at_the_diffraction_dialog():
    tcl_path = os.path.join(ROOT, "plugins", "star_psf", "star_psf.tcl")
    tcl = open(tcl_path, encoding="utf-8").read()
    start = tcl.index("proc CatalogPanelSimPSFDiffraction ")
    end = tcl.index("proc CatalogPanelSimPSFUpdateFilters ")
    block = tcl[start:end]
    assert 'wm title $w "Diffraction PSF"' in block
    assert "Circular pupil: in-tree" in block
    assert "--sim-telescope\" \"ogfmeas\"" in block or '"--sim-telescope" "ogfmeas"' in block
    assert "sim-focus-offset" not in block
    assert "NIRSPEC" not in block
    assert "WFPC2" not in block
    assert 'CatalogPanelSimPSFUpdateFilters $f.cfilt ogfmeas' in tcl

    manifest = json.load(open(
        os.path.join(ROOT, "plugins", "star_psf", "plugin.json"), encoding="utf-8"))
    steps = {step["id"]: step for step in manifest["steps"]}
    item = steps["diffraction"]
    assert item["proc"] == "CatalogPanelSimPSFDiffraction"
    assert item["label"] == "Diffraction PSF..."
    assert item["records"] == ["psf.diffraction"]
    telescope = next(p for p in manifest["params"] if p["name"] == "sim-telescope")
    assert telescope["choices"] == ["auto", "ogfmeas", "jwst", "hst"]
