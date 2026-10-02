"""Registration-error estimate and ZOGY option wiring (item 2).  Synthetic data only."""
import os, sys, types
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pytest
from scipy import ndimage as ndi
from moving import regerr, imaging as I


def _pair(shift, n=400, nstar=40, seed=0, sky_sig=0.02):
    """target and template with the same Gaussian stars; the template is displaced by `shift` = (dx, dy) pixels (+ independent noise)."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:n, 0:n]
    xs = rng.uniform(30, n - 30, nstar); ys = rng.uniform(30, n - 30, nstar); fl = rng.uniform(40, 200, nstar)
    def img(dx, dy):
        a = np.zeros((n, n))
        for x, y, f in zip(xs, ys, fl):
            a += f / (2 * np.pi * 1.7 ** 2) * np.exp(-((xx - x - dx) ** 2 + (yy - y - dy) ** 2) / (2 * 1.7 ** 2))
        return a + rng.normal(0, sky_sig, (n, n))
    T = img(0, 0); R = img(*shift)
    ch = types.SimpleNamespace(data=T.astype(np.float32), bad=np.zeros((n, n), bool), saturated=np.zeros((n, n), bool), cr=np.zeros((n, n), bool), shape=(n, n))
    return ch, R.astype(np.float32)


def _measure(ch, R):
    bt, rt = I.background(ch.data, ch.bad); br, _ = I.background(R, np.zeros(R.shape, bool))
    return regerr.measure_registration(ch, R, np.zeros(R.shape, bool), bt, br, rms_t=rt, nin_mean=3.0)


def test_registration_recovers_known_offsets_per_axis():
    ch, R = _pair((0.30, 0.0), seed=1)
    r = _measure(ch, R)
    assert r is not None and r["n_stars"] >= 20
    assert r["sigma_x"] == pytest.approx(0.30, abs=0.06)       # the x registration error is found ...
    assert r["sigma_y"] < 0.08                                  # ... and y (no error) stays small


def test_registration_zero_for_perfect_pair_and_centroid_noise_subtracted():
    ch, R = _pair((0.0, 0.0), seed=2)
    r = _measure(ch, R)
    assert r["sigma_x"] < 0.06 and r["sigma_y"] < 0.06
    assert r["raw_x"] >= r["sigma_x"]                           # raw rms includes the centroid noise, sigma is the corrected value


def test_registration_none_when_too_few_sources():
    ch, R = _pair((0.2, 0.2), nstar=3, seed=3)
    assert _measure(ch, R) is None


def test_registration_ignores_gross_outliers():
    ch, R = _pair((0.1, 0.1), seed=4)
    R2 = ndi.shift(R, (2.5, 2.5), order=1)                      # a 2.5 px offset is not a "registration error" of the stars (> max_off_pix): blends/wrong match
    r = _measure(ch, R2)
    assert r is None or r["sigma_x"] < 0.1 or r["n_stars"] < 8


def test_cli_options_default_off_and_map_to_difference_kwargs():
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "ds9", "library"))
    import ds9_moving as M
    ns = types.SimpleNamespace
    assert M._detect_opts(ns()) == {}                          # nothing set -> defaults of difference_chip
    assert M._detect_opts(ns(source_noise=False, astrom="off", psf_tile=0, template_psf="target")) == {}
    o = M._detect_opts(ns(source_noise=True, astrom="measure", psf_tile=256, template_psf="measure"))
    assert o == dict(source_noise=True, astrom_sigma="measure", psf_tile=256, template_psf="measure")
    assert M._detect_opts(ns(astrom="sidecar"))["astrom_sigma"] is True
    assert M._detect_opts(ns(astrom="0.3"))["astrom_sigma"] == pytest.approx(0.3)


def test_measure_option_flows_through_difference_chip():
    """difference_chip(astrom_sigma='measure') records the measured registration and the astrometric term is active; falls back cleanly."""
    from moving import detect as D
    ch, R = _pair((0.25, 0.0), n=300, nstar=30, seed=5)
    from astropy.wcs import WCS
    w = WCS(naxis=2); w.wcs.crpix = [150, 150]; w.wcs.crval = [10.0, 10.0]; w.wcs.cdelt = [-0.05 / 3600, 0.05 / 3600]; w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    ch.wcs = w; ch.texp = 100.0; ch.pixscale = 0.05; ch.align = None; ch.name = "t"; ch.path = "t.fits"
    psf = I.gaussian_psf(2.0, 25)
    out = D.difference_chip(ch, R.copy(), nin=np.full(ch.shape, 3.0), psf_t=psf, fw_t=2.0, n_t=0, astrom_sigma="measure", cr_reject=False, trail_fit=False, realbogus=False)
    reg = out["info"]["registration"]
    assert isinstance(reg, dict) and reg["sigma_x"] > 0.1
    assert "astrom_sigma_pix" in out["info"]
