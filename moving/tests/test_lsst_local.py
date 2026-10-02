"""LSST/Rubin local-directory path on a SYNTHETIC Rubin-style directory (no real DP1 file was available: DP1 needs a Rubin Science Platform login).
The files imitate the butler FITS layout of a calibrated visit image: PRIMARY (EXPTIME, FILTER/BAND, MJD-OBS), IMAGE (nJy, TAN WCS, 0.2"/px),
MASK (int32, bit numbers in MP_* header keys), VARIANCE.  What is verified: footprint search, candidate fields, the IMAGE-only chip loader with
variance -> err and mask planes -> bad/saturated/cr, the CLI fetch mode.  NOT verified: the real DP1 header keywords (assumed from the LSST stack
documentation)."""
import json, os, subprocess, sys
import numpy as np
import pytest
from astropy.io import fits
from astropy.wcs import WCS

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, ROOT)
from moving import lsst, imaging

RA0, DEC0 = 53.13, -28.10        # ECDFS, one of the DP1 fields


def make_visit(path, mjd, band, shift=0.0, ra0=RA0, dec0=DEC0, n=300, seed=1):
    rng = np.random.default_rng(seed)
    w = WCS(naxis=2)
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]; w.wcs.crval = [ra0, dec0]; w.wcs.crpix = [n / 2 + 0.5, n / 2 + 0.5]
    w.wcs.cdelt = [-0.2 / 3600, 0.2 / 3600]
    img = rng.normal(0, 30.0, (n, n)).astype("f4") + 100.0        # nJy
    yy, xx = np.mgrid[0:n, 0:n]
    for sx, sy in [(60, 70), (200, 90), (150, 220), (240, 240)]:
        img += 5000.0 * np.exp(-((xx - sx) ** 2 + (yy - sy) ** 2) / (2 * 1.4 ** 2))
    img += 3000.0 * np.exp(-((xx - 100 - shift) ** 2 + (yy - 150) ** 2) / (2 * 1.4 ** 2))          # mover
    var = np.full((n, n), 900.0, "f4")
    mask = np.zeros((n, n), "i4")
    mask[:, :4] |= 1 << 0          # BAD
    mask[10:12, 10:12] |= 1 << 1   # SAT
    mask[50, 50] |= 1 << 2         # CR
    hp = fits.PrimaryHDU(); hp.header.update(EXPTIME=30.0, FILTER=band, **{"MJD-OBS": mjd}, TELESCOP="Rubin", INSTRUME="LSSTComCam")
    hi = fits.ImageHDU(img, name="IMAGE", header=w.to_header()); hi.header["BUNIT"] = "nJy"; hi.header["MJD-OBS"] = mjd
    hm = fits.ImageHDU(mask, name="MASK", header=w.to_header()); hm.header.update(MP_BAD=0, MP_SAT=1, MP_CR=2)
    hv = fits.ImageHDU(var, name="VARIANCE", header=w.to_header())
    fits.HDUList([hp, hi, hm, hv]).writeto(path, overwrite=True)


@pytest.fixture(scope="module")
def lsst_dir(tmp_path_factory):
    d = tmp_path_factory.mktemp("rubin")
    for k, (mjd, band) in enumerate([(60600.1, "r"), (60601.1, "r"), (60602.1, "i")]):
        sub = d / "dp1" / "visit_image" / band
        sub.mkdir(parents=True, exist_ok=True)
        make_visit(str(sub / ("visit_%d.fits" % (k + 1))), mjd, band, shift=3.0 * k, seed=k)
    make_visit(str(d / "elsewhere.fits"), 60600.1, "r", ra0=150.0, dec0=2.0)           # different sky position
    (d / "notes.txt").write_text("not a fits file")
    return str(d)


def test_find_local_footprint(lsst_dir):
    c = lsst.find_local(RA0, DEC0, lsst_dir, radius_arcmin=0.5)
    assert sorted(x["obs_id"] for x in c) == ["visit_1.fits", "visit_2.fits", "visit_3.fits"]
    assert {x["filter"] for x in c} == {"r", "i"} and abs(sorted(x["mjd"] for x in c)[0] - 60600.1) < 1e-6 and c[0]["texp"] == 30.0
    assert lsst.find_local(10.0, -60.0, lsst_dir) == []


def test_query_candidates_local_and_field(lsst_dir, monkeypatch):
    monkeypatch.setattr(lsst, "probe", lambda *a, **k: (_ for _ in ()).throw(AssertionError("no network probe expected when local files exist")))
    c, st = lsst.query_candidates(RA0, DEC0, 0.5, local_dir=lsst_dir)
    assert len(c) == 3 and st["available"] and st["dp1_field"] == "ECDFS"
    monkeypatch.setattr(lsst, "probe", lambda *a, **k: {"sia_dp1": (401, 0)})
    c, st = lsst.query_candidates(RA0, DEC0, 0.5, local_dir=os.path.join(lsst_dir, "nothing"))
    assert c == [] and not st["available"] and "login" in st["reason"]


def test_loader_reads_image_variance_mask(lsst_dir):
    ch = imaging.load_chips(os.path.join(lsst_dir, "dp1", "visit_image", "r", "visit_1.fits"))
    assert len(ch) == 1                                                  # MASK/VARIANCE are not science chips
    c = ch[0]
    assert c.shape == (300, 300) and abs(c.pixscale - 0.2) < 1e-3 and c.zp_ab == 31.4 and abs(c.mjd - (60600.1 + 15 / 86400)) < 1e-6
    assert c.bad[:, :4].all() and not c.bad[100:110, 100:110].any()
    assert c.saturated[10:12, 10:12].all() and c.cr[50, 50] and not c.bad[50, 50]
    assert np.allclose(c.err, 30.0)
    assert c.filter == "r"


def test_cli_fetch_local(lsst_dir, tmp_path):
    env = dict(os.environ, PYTHONPATH=ROOT)
    p = subprocess.run([sys.executable, os.path.join(ROOT, "ds9", "library", "ds9_moving.py"), "--mode", "fetch", "--provider", "lsst", "--ra", str(RA0), "--dec", str(DEC0),
                        "--radius", "0.5", "--lsst-local-dir", lsst_dir, "--download", "--workdir", str(tmp_path)], capture_output=True, text=True, env=env, timeout=120)
    assert p.returncode == 0, p.stderr[-600:]
    cj = json.load(open(tmp_path / "candidates.json"))
    assert len(cj["candidates"]) == 3 and cj["status"]["lsst"]["available"]
    fj = json.load(open(tmp_path / "files.json"))
    assert len(fj["files"]) == 3 and all(f["bytes"] > 0 for f in fj["files"])
