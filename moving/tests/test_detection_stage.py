"""Detection-stage tests (item 1): cosmic-ray rejection, trailed-PSF centroiding, real/bogus classifier, linker use of `rb`.  Synthetic data only."""
import os, sys, json
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pytest
from moving import centroid as CE, crrej, realbogus as RB, pipeline as P


def _field(n=200, sigma=1.2, seed=0):
    rng = np.random.default_rng(seed)
    return rng.normal(0, 1.0, (n, n)).astype(np.float32)


def test_segment_model_flux_and_point_limit():
    yy, xx = np.mgrid[0:61, 0:61]
    m = CE.segment_model((61, 61), 30.0, 30.0, 20.0, 0.3, 1.5, 100.0, yy, xx)
    assert abs(m.sum() - 100.0) < 0.5                      # total flux = amp
    p = CE.segment_model((61, 61), 30.0, 30.0, 0.0, 0.3, 1.5, 100.0, yy, xx)
    g = np.exp(-((xx - 30) ** 2 + (yy - 30) ** 2) / (2 * 1.5 ** 2)); g = 100 * g / g.sum()
    assert np.abs(p - g).max() < 0.02 * g.max()            # L -> 0 is a point source


@pytest.mark.parametrize("L,th", [(30.0, 0.4), (55.0, 2.0), (12.0, -0.8)])
def test_fit_segment_recovers_trail_with_cosmic_rays(L, th):
    img = _field(260, seed=3)
    yy, xx = np.mgrid[0:260, 0:260]
    x0, y0 = 128.3, 131.7
    img += CE.segment_model(img.shape, x0, y0, L, th, 1.5, 700.0, yy, xx).astype(np.float32)
    mask = np.zeros(img.shape, bool)
    rng = np.random.default_rng(1)
    for _ in range(6):                                    # CR hits next to / on the trail, flagged in the mask
        cx = int(x0 + rng.uniform(-L / 2, L / 2) * np.cos(th)); cy = int(y0 + rng.uniform(-L / 2, L / 2) * np.sin(th))
        img[cy, cx] += 3000.0; mask[cy - 1:cy + 2, cx - 1:cx + 2] = True
    r = CE.fit_segment(img, x0 + 4.0, y0 - 3.0, L * 0.7, th + 0.1, 1.5, 1.0, mask=mask)
    assert r["ok"]
    assert np.hypot(r["x"] - x0, r["y"] - y0) < 0.6
    assert abs(r["L"] - L) < 0.12 * L + 1.0


def test_inpaint_replaces_masked_pixels():
    a = np.ones((20, 20), np.float32) * 5.0
    m = np.zeros(a.shape, bool); m[10, 10] = True; a[10, 10] = 1e5
    out = crrej.inpaint(a, m)
    assert out[10, 10] == pytest.approx(5.0)


def test_lac_mask_finds_cosmic_ray_not_star():
    rng = np.random.default_rng(5)
    img = rng.normal(0.1, 0.02, (300, 300)).astype(np.float32)          # e-/s
    yy, xx = np.mgrid[0:300, 0:300]
    img += (30.0 * np.exp(-((xx - 150) ** 2 + (yy - 100) ** 2) / (2 * 1.6 ** 2))).astype(np.float32)       # star
    img[200, 60] += 50.0; img[200, 61] += 20.0; img[201, 60] += 15.0                                       # cosmic ray
    m = crrej.lac_mask(img, np.zeros(img.shape, bool), texp=500.0, readnoise=5.0)
    assert m[200, 60]
    assert not m[100, 150]                                                 # PSF-shaped star is protected


def test_realbogus_model_file_is_small_and_valid():
    assert os.path.isfile(RB.MODEL_PATH)
    assert os.path.getsize(RB.MODEL_PATH) < 120_000
    m = RB.load(force=True)
    assert m is not None and len(m["trees"]) >= 10
    assert tuple(json.load(open(RB.MODEL_PATH))["features"]) == RB.FEATURES


def test_realbogus_separates_cr_like_from_psf_like():
    base = dict(sign=1, snr=12.0, channel="point", elong=1.2, a_pix=1.5, b_pix=1.3, npix=14, tpl_snr=0.5, neg_frac=0.0, pos_frac=0.1, sharp=0.5,
                sharp_psf=0.45, on_cr=False, f_r1=0.5, f_r2=0.4, f_r3=0.1, pk_nb=0.6, psf_chi2=1.2, psf_amp=1.0, fine=1.0, lap=3.0, n_hi=6, asym=0.1,
                lac3=0, lac_n7=0, arch3=0, n_near=0, near_bad=False)
    cr = dict(base, sharp=0.95, f_r1=0.93, f_r2=0.07, f_r3=0.0, pk_nb=0.1, psf_chi2=40.0, psf_amp=1.8, fine=25.0, lap=40.0, n_hi=2, npix=3, a_pix=0.7,
              lac3=1, lac_n7=5, arch3=1)
    dets = [dict(base), dict(cr)]
    RB.score_dets(dets)
    assert dets[0]["rb"] > dets[1]["rb"]
    # (only the ordering is asserted: hand-made feature dicts are not the training distribution; real separation is measured on the benchmark)
    more_lac = dict(base, lac3=1, lac_n7=9, arch3=1)
    d2 = [dict(base), more_lac]; RB.score_dets(d2)
    assert d2[1]["rb"] < d2[0]["rb"]


def test_realbogus_tree_walk_matches_sklearn():
    from sklearn.ensemble import GradientBoostingClassifier
    rng = np.random.default_rng(0)
    X = rng.normal(size=(600, len(RB.FEATURES))); y = (X[:, 3] + 0.5 * X[:, 7] + 0.3 * rng.normal(size=600) > 0.8).astype(int)
    gb = GradientBoostingClassifier(n_estimators=20, max_depth=3, random_state=0).fit(X, y)
    m = RB.load_from_dict(RB.export_sklearn(gb, np.zeros(len(RB.FEATURES))))
    z, p = RB.predict(X, m)
    assert np.abs(p - gb.predict_proba(X)[:, 1]).max() < 1e-4


def _det(ex, t, ra, dec, snr, rb=None, cls="point"):
    d = dict(ex=ex, t=t, ra=ra, dec=dec, snr=snr, sign=1, cls=cls, channel="point", flux_e_s=1.0, sig_pos_arcsec=0.05, sharp=0.4, sharp_psf=0.4, on_cr=False,
             texp=500.0, tpl_snr=0.0)
    if rb is not None:
        d["rb"] = rb
    return d


def test_link_pool_orders_by_rb_not_snr():
    """With a per-exposure cap of 1 the real/bogus ordering keeps the faint real detection; without `rb` the brightest one wins (old behaviour)."""
    dets = []
    for e in range(3):
        dets.append(_det(e, 53000.0 + 0.01 * e, 150.0 + 1e-4 * e, 2.0, 50.0, rb=0.01 if True else None))      # bright bogus
        dets.append(_det(e, 53000.0 + 0.01 * e, 150.01 + 1e-4 * e, 2.0, 10.0, rb=0.9))                      # faint real
    seen = {}
    from moving import tracklet as T
    sv = T.link_exposures
    T.link_exposures = lambda sel, **kw: seen.setdefault("sel", sel) and []
    try:
        P.link_detections(dets, None, snr_min=8.0, max_per_exposure=1, veto=False)
        assert all(s["snr"] == 10.0 for s in seen["sel"]) and len(seen["sel"]) == 3
        seen.clear()
        P.link_detections(dets, None, snr_min=8.0, max_per_exposure=1, veto=False, use_rb=False)
        assert all(s["snr"] == 50.0 for s in seen["sel"])
        seen.clear()
        for d in dets:
            d.pop("rb")
        P.link_detections(dets, None, snr_min=8.0, max_per_exposure=1, veto=False)
        assert all(s["snr"] == 50.0 for s in seen["sel"])                    # no rb in the data -> S/N order as before
    finally:
        T.link_exposures = sv
