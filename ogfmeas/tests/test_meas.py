"""Synthetic checks for the in-tree measurement library. Expectations are analytic, not another code's numbers."""
import math
import os

import numpy as np
import pytest

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "2")

from ogfmeas import (
    Background, extract, flux_radius, kron_radius, set_extract_pixstack,
    set_sub_object_limit, sum_circle, sum_circann, sum_ellipse, winpos,
)
from ogfmeas.nbody import propagate
from ogfmeas.photoz import template_chi2
from ogfmeas.sersic import fit_sersic, render
from ogfmeas.sextract import catalog_table, main


def _gauss(shape, x0, y0, amp, sig):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return amp * np.exp(-0.5 * ((xx - x0) ** 2 + (yy - y0) ** 2) / sig ** 2)


def test_constant_background_and_mask():
    rng = np.random.default_rng(4)
    image = 7.5 + rng.normal(0.0, 0.15, size=(64, 64))
    bkg = Background(image, bw=32, bh=32)
    assert abs(bkg.globalback - 7.5) < 0.08
    assert 0.05 < bkg.globalrms < 0.3
    assert np.allclose(np.array(bkg), bkg.back())
    sky = np.full((64, 64), 3.0)
    sky[10:30, 10:30] = 800.0
    mask = np.zeros_like(sky, dtype=bool)
    mask[10:30, 10:30] = True
    masked = Background(sky, mask=mask, bw=16, bh=16)
    assert abs(masked.globalback - 3.0) < 0.05
    data = np.full((8, 8), 5.0)
    Background(np.full((8, 8), 3.0)).subfrom(data)
    assert np.allclose(data, 2.0)


def test_one_gaussian_is_found_and_apertures_match():
    amp, sig = 80.0, 2.0
    image = _gauss((64, 64), 32.0, 30.0, amp, sig)
    objects = extract(image, 5.0, err=1.0, minarea=5, deblend_nthresh=16)
    assert len(objects) == 1
    obj = objects[0]
    assert abs(obj["x"] - 32.0) < 0.4
    assert abs(obj["y"] - 30.0) < 0.4
    assert obj["flag"] == 0
    assert obj["a"] >= obj["b"]
    assert -math.pi / 2 <= obj["theta"] <= math.pi / 2
    analytic = amp * 2.0 * math.pi * sig ** 2
    flux, fluxerr, flag = sum_circle(image, 32.0, 30.0, 12.0, err=1.0, subpix=5)
    assert flag[0] == 0
    assert abs(flux[0] - analytic) / analytic < 0.05
    assert fluxerr[0] > 0
    ell, _, ell_flag = sum_ellipse(image, [32.0], [30.0], [12.0], [12.0], [0.4], 1.0, subpix=5)
    assert ell_flag[0] == 0
    assert abs(ell[0] - flux[0]) / flux[0] < 1e-3
    inner, _, _ = sum_circann(image, 32.0, 30.0, 0.0, 4.0, subpix=5)
    outer, _, _ = sum_circann(image, 32.0, 30.0, 4.0, 12.0, subpix=5)
    assert abs((inner[0] + outer[0]) - flux[0]) / flux[0] < 0.02


def test_deblend_two_peaks_and_mask_and_border():
    image = _gauss((48, 48), 16.0, 24.0, 300.0, 2.4) + _gauss((48, 48), 30.0, 24.0, 250.0, 2.4)
    objects, seg = extract(image, 8.0, minarea=8, deblend_nthresh=32, deblend_cont=0.01,
                           segmentation_map=True)
    assert len(objects) == 2
    centers = sorted((float(o["x"]), float(o["y"])) for o in objects)
    assert abs(centers[0][0] - 16.0) < 1.2
    assert abs(centers[1][0] - 30.0) < 1.2
    assert set(np.unique(seg)) >= {0, 1, 2}
    assert int(seg.max()) == 2
    masked = np.zeros(image.shape, dtype=bool)
    masked[:, :24] = True
    kept = extract(image, 8.0, mask=masked, minarea=8, deblend_nthresh=16)
    assert len(kept) == 1
    assert kept[0]["x"] > 24
    edge = _gauss((40, 40), 0.2, 20.0, 100.0, 1.3)
    found = extract(edge, 5.0, minarea=3, deblend_nthresh=1, filter_kernel=None)
    assert len(found) >= 1
    assert int(found[0]["flag"]) & 1


def test_kron_flux_radius_and_winpos():
    yy, xx = np.mgrid[0:80, 0:80]
    radius = np.sqrt((xx - 40.0) ** 2 + (yy - 40.0) ** 2)
    scale = 3.0
    image = np.exp(-radius / scale)
    kron, kflag = kron_radius(image, 40.0, 40.0, scale, scale, 0.0, 6.0)
    assert kflag[0] == 0
    assert 0.5 < kron[0] < 4.0
    half, hflag = flux_radius(_gauss((64, 64), 32.0, 32.0, 50.0, 2.0), 32.0, 32.0, 15.0, 0.5)
    assert hflag[0] == 0
    assert 1.5 < half[0] < 3.3
    shifted = _gauss((48, 48), 20.35, 22.7, 120.0, 1.6)
    xw, yw, wflag = winpos(shifted, 20.0, 23.0, 1.8, subpix=5)
    assert wflag[0] == 0
    assert abs(xw[0] - 20.35) < 0.15
    assert abs(yw[0] - 22.7) < 0.15


def test_limits_raise_the_retry_strings():
    import sys
    meas = sys.modules["ogfmeas.extract"]
    bright = np.full((24, 24), 50.0)
    old_stack = meas._PIXSTACK
    old_sub = meas._SUB_LIMIT
    try:
        set_extract_pixstack(10)
        with pytest.raises(RuntimeError, match="pixel buffer full"):
            extract(bright, 1.0, minarea=1, filter_kernel=None)
        set_extract_pixstack(old_stack)
        image = _gauss((48, 48), 16.0, 24.0, 300.0, 2.4) + _gauss((48, 48), 30.0, 24.0, 250.0, 2.4)
        set_sub_object_limit(1)
        with pytest.raises(RuntimeError, match="deblending overflow"):
            extract(image, 8.0, minarea=8, deblend_nthresh=32, deblend_cont=0.01)
    finally:
        set_extract_pixstack(old_stack)
        set_sub_object_limit(old_sub)


def test_sersic_recovers_a_planted_exponential():
    image = render((51, 51), 25.0, 24.0, re=4.0, n=1.0, ie=40.0, q=0.7, theta=0.35, background=1.0)
    fit = fit_sersic(image, x0=25.0, y0=24.0)
    assert abs(fit["n"] - 1.0) < 0.25
    assert abs(fit["re"] - 4.0) < 0.6
    assert abs(fit["q"] - 0.7) < 0.08


def test_template_redshift_and_circular_orbit():
    wave = np.linspace(3000.0, 9000.0, 240)
    rest = np.exp(-0.5 * ((wave - 5000.0) / 80.0) ** 2)
    z_true = 0.4
    observed = np.interp(wave, wave * (1.0 + z_true), rest, left=0.0, right=0.0)
    fit = template_chi2(wave, observed, np.full_like(observed, 0.02), [(wave, rest)],
                        np.linspace(0.0, 1.0, 51))
    assert abs(fit["z"] - z_true) <= 0.02
    assert abs(float(np.sum(fit["pdf"])) - 1.0) < 1e-6
    separation = 1.0
    mass = 1.0
    omega = math.sqrt(2.0 * mass / separation ** 3)
    period = 2.0 * math.pi / omega
    speed = omega * separation / 2.0
    positions = np.array([[-0.5, 0.0, 0.0], [0.5, 0.0, 0.0]])
    velocities = np.array([[0.0, speed, 0.0], [0.0, -speed, 0.0]])
    steps = 1600
    final, _ = propagate(positions, velocities, [mass, mass], period / steps, steps, G=1.0)
    assert np.allclose(final, positions, atol=0.02)


def test_noise_loader_does_not_import_the_external_package():
    import sys
    sys.modules.pop("sep", None)
    sys.modules.pop("sep_pjw", None)
    from ogfkit.noise import _sep, source_mask
    mod = _sep()
    assert mod.__name__ == "ogfmeas"
    assert "sep" not in sys.modules
    assert "sep_pjw" not in sys.modules
    image = np.zeros((32, 32))
    image[12:18, 14:20] = 30.0
    mask = source_mask(image, thresh=3.0, minarea=4, bw=16)
    assert bool(mask[15, 17])


def test_sextract_writes_a_one_based_catalog(tmp_path):
    image = _gauss((48, 48), 20.0, 22.0, 90.0, 1.8) + 2.0
    text = catalog_table(image, thresh=8.0, minarea=5, bw=48)
    rows = [line.split("\t") for line in text.strip().splitlines()]
    assert rows[0][1] == "X_IMAGE"
    assert len(rows) == 2
    assert abs(float(rows[1][1]) - 21.0) < 0.6
    assert abs(float(rows[1][2]) - 23.0) < 0.6
    from astropy.io import fits
    path = tmp_path / "one.fits"
    out = tmp_path / "cat.tsv"
    fits.PrimaryHDU(image.astype(np.float32)).writeto(path, overwrite=True)
    assert main([str(path), "--out", str(out), "--thresh", "8", "--bw", "48"]) == 0
    written = out.read_text(encoding="utf-8").strip().splitlines()
    assert len(written) == 2
