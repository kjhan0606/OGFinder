"""Unit tests for the ICL package (icl/): masking, sky, profiles, colours, f_ICL."""
import os
import sys

import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from icl.config import ICLConfig  # noqa: E402
from icl import colors, measure, background, masking, profile  # noqa: E402


def _write(path, data, **cards):
    from astropy.io import fits
    h = fits.Header()
    for k, v in cards.items():
        h[k] = v
    fits.PrimaryHDU(data.astype(np.float32), h).writeto(path, overwrite=True)


def test_header_ab_zeropoint_hst_and_jwst():
    from astropy.io import fits
    h = fits.Header()
    # HST ACS/WFC F814W (PHOTFLAM, PHOTPLAM) -> 25.94-25.95
    h["PHOTFLAM"] = 7.0723e-20
    h["PHOTPLAM"] = 8045.5
    assert colors.header_ab_zeropoint(h) == pytest.approx(25.947, abs=0.01)
    j = fits.Header()
    j["PIXAR_SR"] = 9.31e-14
    assert colors.header_ab_zeropoint(j) == pytest.approx(-6.10 - 2.5 * np.log10(9.31e-14))
    assert colors.header_ab_zeropoint(fits.Header()) is None


def test_multiband_profiles_apply_mask_and_per_band_zeropoints(tmp_path):
    ny = nx = 101
    flat1 = np.full((ny, nx), 1.0)
    flat2 = np.full((ny, nx), 1.0)
    # a bright contaminant that the mask must remove from both bands
    flat1[60:70, 60:70] = 1e4
    flat2[60:70, 60:70] = 1e4
    mask = np.zeros((ny, nx), bool)
    mask[58:72, 58:72] = True
    p1, p2 = tmp_path / "b1.fits", tmp_path / "b2.fits"
    _write(p1, flat1)
    _write(p2, flat2)
    cfg = ICLConfig(profile_rmin=5, profile_rmax=45, profile_nsteps=6, profile_spacing="linear",
                    mag_zeropoint=25.0, pixel_scale=0.06)
    prof = colors.measure_multiband_icl_profiles({"A": str(p1), "B": str(p2)}, mask, 50, 50, cfg,
                                                 zeropoints={"A": 26.0})
    # the masked contaminant is gone: every annulus has the flat value -> mu = zp + 2.5 log10(0.06^2)
    for p in prof["B"]:
        assert p["MU"] == pytest.approx(25.0 + 2.5 * np.log10(0.06 ** 2), abs=1e-6)
    col = colors.compute_color_profiles(prof, "A", "B")
    assert len(col) == len(prof["A"])
    for c in col:
        assert c["COLOR"] == pytest.approx(1.0, abs=1e-6)   # = ZP difference, no contamination
    # without the mask the contaminated annulus is brighter than the flat level
    prof_nomask = colors.measure_multiband_icl_profiles({"B": str(p2)}, None, 50, 50, cfg)
    assert min(p["FLUX"] / p["NPIX"] for p in prof_nomask["B"]) == pytest.approx(1.0)
    assert max(p["FLUX"] / p["NPIX"] for p in prof_nomask["B"]) > 10


def test_multiband_profiles_mask_shape_mismatch(tmp_path):
    p = tmp_path / "b.fits"
    _write(p, np.ones((20, 20)))
    with pytest.raises(ValueError):
        colors.measure_multiband_icl_profiles({"A": str(p)}, np.zeros((10, 10), bool), 10, 10,
                                              ICLConfig(profile_rmax=8, profile_nsteps=3))


def test_constant_sky_apertures_recovers_offset():
    rng = np.random.default_rng(3)
    img = rng.normal(-0.02, 0.1, (600, 600))
    yy, xx = np.mgrid[0:600, 0:600]
    img += 50.0 * np.exp(-np.hypot(xx - 300, yy - 300) / 30.0)      # "cluster" in the centre
    mask = np.hypot(xx - 300, yy - 300) < 120
    sky, err, cen, means = background.constant_sky_apertures(img, mask, 300, 300, 200, n_apertures=20,
                                                             radius=15, seed=0)
    assert len(means) == 20
    assert np.all(np.hypot(cen[:, 0] - 300, cen[:, 1] - 300) >= 200)
    assert sky == pytest.approx(-0.02, abs=4 * err + 0.003)


def test_hot_cold_mask_finds_big_and_small_sources():
    rng = np.random.default_rng(1)
    img = rng.normal(0, 1.0, (400, 400))
    yy, xx = np.mgrid[0:400, 0:400]
    img += 40 * np.exp(-0.5 * ((xx - 200) ** 2 + (yy - 200) ** 2) / 15.0 ** 2)   # big galaxy
    img += 30 * np.exp(-0.5 * ((xx - 215) ** 2 + (yy - 190) ** 2) / 1.2 ** 2)    # small source on it
    img += 30 * np.exp(-0.5 * ((xx - 80) ** 2 + (yy - 320) ** 2) / 1.2 ** 2)     # isolated small source
    m, seg_c, seg_h, obj = masking.hot_cold_mask(img, cold_thresh=3, cold_minarea=100, hot_thresh=3,
                                                 hot_minarea=5)
    assert m[200, 200] and m[190, 215] and m[320, 80]
    assert seg_c[200, 200] > 0
    assert seg_h[190, 215] > 0            # found on top of the big galaxy by the unsharp-masked pass
    assert not m[50, 350]
    assert m.mean() < 0.2


def test_icl_fraction_pixels_analytic():
    # two-level image: 100 px at mu = 25 (bright) and 300 px at mu = 26.5 (faint), zp 25, scale 1"/px
    img = np.zeros((20, 20))
    img.flat[:100] = 10 ** (-0.4 * (25 - 25))       # flux 1.0 per px -> mu 25
    img.flat[100:400] = 10 ** (-0.4 * (26.5 - 25))  # mu 26.5
    valid = np.ones_like(img, bool)
    lo = 1.0 * 100
    hi = 300 * 10 ** (-0.4 * 1.5)
    r = measure.icl_fraction_pixels(img, valid, 25.0, 1.0, 26.0, 27.0, n_bootstrap=0)
    assert r["F_ICL"] == pytest.approx(hi / (lo + hi))
    assert r["N_PIX_ICL"] == 300
    # a dimming offset of 1 mag moves the faint pixels to mu_rest = 25.5: no ICL above 26
    r2 = measure.icl_fraction_pixels(img, valid, 25.0, 1.0, 26.0, 27.0, mu_offset=1.0, n_bootstrap=0)
    assert r2["F_ICL"] == 0.0
    # masked (invalid) faint pixels leave numerator and denominator
    valid[0, 10:] = False
    r3 = measure.icl_fraction_pixels(img, valid, 25.0, 1.0, 26.0, None, n_bootstrap=50, block=5)
    assert r3["N_PIX"] == 390
    assert 0 < r3["F_ICL"] < 1 and r3["F_ICL_ERR"] >= 0


def test_pixel_sb_map_nan_for_nonpositive():
    mu = measure.pixel_sb_map(np.array([[1.0, 0.0, -1.0]]), 25.0, 0.06)
    assert mu[0, 0] == pytest.approx(25.0 + 2.5 * np.log10(0.06 ** 2))
    assert np.isnan(mu[0, 1]) and np.isnan(mu[0, 2])


def test_sb_profile_flat_image_and_isophotal_radius():
    img = np.full((201, 201), 0.01)
    cfg = ICLConfig(profile_rmin=5, profile_rmax=90, profile_nsteps=10, mag_zeropoint=25, pixel_scale=0.06)
    prof = profile.measure_sb_profile(img, 100, 100, cfg)
    mu0 = 25 - 2.5 * np.log10(0.01) + 2.5 * np.log10(0.06 ** 2)
    assert all(p["MU"] == pytest.approx(mu0) for p in prof)
    # exponential profile: R(mu) interpolation lands between the bracketing annuli
    yy, xx = np.mgrid[0:201, 0:201]
    expo = np.exp(-np.hypot(xx - 100, yy - 100) / 20.0)
    prof2 = profile.measure_sb_profile(expo, 100, 100, cfg)
    target = prof2[4]["MU"]
    r = measure.compute_isophotal_radii(prof2, [target])[target]
    assert r == pytest.approx(prof2[4]["R_PIX"], rel=1e-6)
