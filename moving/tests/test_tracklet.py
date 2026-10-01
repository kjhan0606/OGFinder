"""Tracklet fit and linking with and without observer parallax (synthetic, no network)."""
import numpy as np
from moving import tracklet as T


def _obs_offsets(n, rng):
    # HST-like orbital offsets (AU, equatorial): 7000 km radius circle in the x-y plane over ~40 min
    ph = np.linspace(0, 0.45, n)
    return np.c_[7000 * np.cos(ph), 7000 * np.sin(ph), 1500 * np.sin(ph)] / 149597870.7


def test_fit_tracklet_linear_motion():
    t = np.array([0.0, 0.01, 0.02, 0.03]) + 58000.0
    rate = 10.0                                    # arcsec/hour
    dec0 = 2.0
    ra = 150.0 + np.cumsum(np.r_[0, np.diff(t)]) * 24 * rate / 3600 / np.cos(np.radians(dec0))
    dec = np.full(4, dec0)
    f = T.fit_tracklet(t, ra, dec, [0.02] * 4)
    assert abs(f["rate_ash"] - rate) < 0.1
    assert abs(f["pa_deg"] - 90.0) < 0.5            # moving east


def test_parallax_link_recovers_track():
    rng = np.random.default_rng(2)
    n = 4; t = 58000.0 + np.array([0, 0.0083, 0.0166, 0.0250])
    off = _obs_offsets(n, rng)
    k_true = 1 / 2.5                                # Delta = 2.5 AU
    ra0, de0 = 150.0, 2.0
    mu = np.array([6.0, -3.0])                      # arcsec/hour (geocentric motion), x=E, y=N
    pos = []
    for i in range(n):
        dt_h = (t[i] - t[0]) * 24
        x = mu[0] * dt_h - k_true * off[i, 1] * 206264.806 * 0 - k_true * (off[i] @ np.array([-np.sin(np.radians(ra0)), np.cos(np.radians(ra0)), 0])) * 206264.806
        a0, d0 = np.radians(ra0), np.radians(de0)
        ey = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
        y = mu[1] * dt_h - k_true * (off[i] @ ey) * 206264.806
        pos.append((ra0 + x / 3600 / np.cos(np.radians(de0)), de0 + y / 3600))
    pos = np.array(pos)
    dets = [dict(ex=i, t=t[i], ra=pos[i, 0], dec=pos[i, 1], sig=0.05, id=i, flux=1.0) for i in range(n)]
    # add distractors that do not move consistently
    for j in range(40):
        for i in range(n):
            dets.append(dict(ex=i, t=t[i], ra=ra0 + rng.uniform(-0.01, 0.01), dec=de0 + rng.uniform(-0.01, 0.01), sig=0.05, id=100 + j * 10 + i, flux=1.0))
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=4, obs_off_au={i: off[i] for i in range(n)}, min_disp_arcsec=0.5)
    best = trs[0]
    assert set(best["members"]) == {0, 1, 2, 3}
    assert abs(best["rate_ash"] - np.hypot(*mu)) < 0.5


def _synth(n=5, k_true=0.4, mu=(6.0, -3.0), n_dist=60, seed=3, dt_days=0.0083, sig=0.05, ra0=150.0, de0=2.0, extent=0.01, drop=()):
    """One mover with parallax (+ uniform distractors).  Returns dets, obs-offset dict, member ids."""
    rng = np.random.default_rng(seed)
    t = 58000.0 + dt_days * np.arange(n); off = _obs_offsets(n, rng)
    a0, d0 = np.radians(ra0), np.radians(de0)
    ex_ = np.array([-np.sin(a0), np.cos(a0), 0.0]); ey_ = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
    dets = []; ids = []
    for i in range(n):
        if i in drop:
            continue
        h = (t[i] - t[0]) * 24
        x = mu[0] * h - k_true * (off[i] @ ex_) * 206264.806 + rng.normal(0, sig); y = mu[1] * h - k_true * (off[i] @ ey_) * 206264.806 + rng.normal(0, sig)
        dets.append(dict(ex=i, t=t[i], ra=ra0 + x / 3600 / np.cos(d0), dec=de0 + y / 3600, sig=sig, id=i, flux=1.0)); ids.append(i)
    for j in range(n_dist):
        for i in range(n):
            dets.append(dict(ex=i, t=t[i], ra=ra0 + rng.uniform(-extent, extent), dec=de0 + rng.uniform(-extent, extent), sig=sig, id=1000 + j * 10 + i, flux=1.0))
    return dets, {i: off[i] for i in range(n)}, ids


def test_three_detections_link():
    dets, off, ids = _synth(n=3, n_dist=40)
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=3, obs_off_au=off, min_disp_arcsec=0.5)
    assert trs and set(trs[0]["members"]) == set(ids)
    assert trs[0]["n"] == 3


def test_five_detections_one_missing_still_links_as_four():
    dets, off, ids = _synth(n=5, drop=(2,), n_dist=40)
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=3, obs_off_au=off, min_disp_arcsec=0.5)
    assert set(trs[0]["members"]) == set(ids) and trs[0]["n_missing"] == 1
    assert abs(trs[0]["inv_delta"] - 0.4) < 0.15


def test_no_duplicates_and_no_subsets():
    dets, off, ids = _synth(n=5, n_dist=40)
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=3, obs_off_au=off, min_disp_arcsec=0.5)
    seen = set()
    for t in trs:
        pairs = {(a, b) for a in t["members"] for b in t["members"] if a < b}
        assert not (pairs & seen)               # no two tracklets share two detections
        seen |= pairs
    assert sum(1 for t in trs if set(ids) <= set(t["members"])) == 1


def test_unphysical_rate_is_rejected():
    # 5000 arcsec/h with k = 0.1 (Delta = 10 AU) is faster than any bound orbit allows
    dets, off, ids = _synth(n=4, k_true=0.1, mu=(5000.0, 0.0), n_dist=10, extent=0.5)
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=4, obs_off_au=off, min_disp_arcsec=0.5)
    assert not any(set(ids) <= set(t["members"]) for t in trs)


def test_chance_alignments_rank_below_a_real_mover():
    """In a dense random field chance alignments exist (k is free), but the real mover must out-rank all of them."""
    dets, off, ids = _synth(n=4, n_dist=150, extent=0.01, k_true=0.35, mu=(8.0, 2.0))
    trs = T.link_exposures(dets, tol_arcsec=0.3, min_exposures=4, obs_off_au=off, min_disp_arcsec=0.5)
    assert set(trs[0]["members"]) == set(ids)
    assert all(t["score"] > trs[0]["score"] for t in trs[1:])
    bg = T.link_exposures([d for d in dets if d["id"] >= 1000], tol_arcsec=0.3, min_exposures=4, obs_off_au=off, min_disp_arcsec=0.5)
    assert all(t["score"] > trs[0]["score"] for t in bg)


def test_mini_injection_recovery_rate():
    """10 movers (6-40 arcsec/h, random direction / Delta) in a field of 100 clutter detections per exposure, 4 exposures."""
    rng = np.random.default_rng(11)
    n = 4; t = 58000.0 + 0.0083 * np.arange(n); off = _obs_offsets(n, rng)
    ra0, de0 = 150.0, 2.0; a0, d0 = np.radians(ra0), np.radians(de0)
    ex_ = np.array([-np.sin(a0), np.cos(a0), 0.0]); ey_ = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
    dets = []; truth = []
    for m in range(10):
        rate = np.exp(rng.uniform(np.log(6), np.log(40))); pa = rng.uniform(0, 2 * np.pi); k = 1 / rng.uniform(1.5, 4.5)
        x0, y0 = rng.uniform(-60, 60, 2); ids = []
        for i in range(n):
            h = (t[i] - t[0]) * 24
            x = x0 + rate * np.sin(pa) * h - k * (off[i] @ ex_) * 206264.806 + rng.normal(0, 0.15)
            y = y0 + rate * np.cos(pa) * h - k * (off[i] @ ey_) * 206264.806 + rng.normal(0, 0.15)
            dets.append(dict(ex=i, t=t[i], ra=ra0 + x / 3600 / np.cos(d0), dec=de0 + y / 3600, sig=0.1, id=10000 * (m + 1) + i, flux=1.0)); ids.append(10000 * (m + 1) + i)
        truth.append(set(ids))
    for j in range(100):
        for i in range(n):
            dets.append(dict(ex=i, t=t[i], ra=ra0 + rng.uniform(-70, 70) / 3600, dec=de0 + rng.uniform(-70, 70) / 3600, sig=0.1, id=j * 10 + i, flux=1.0))
    trs = T.link_exposures(dets, tol_arcsec=0.6, min_exposures=3, obs_off_au={i: off[i] for i in range(n)}, max_tracklets=60)
    rec = sum(any(len(tr & set(t_["members"])) >= 3 for t_ in trs) for tr in truth)
    assert rec >= 8, rec
    # injected objects rank above the chance alignments
    top = trs[:10]
    assert sum(any(len(tr & set(t_["members"])) >= 3 for tr in truth) for t_ in top) >= 7


def test_legacy_linker_still_available():
    dets, off, ids = _synth(n=4, n_dist=10)
    trs = T.link_exposures_legacy(dets, tol_arcsec=0.4, min_exposures=4, obs_off_au=off, min_disp_arcsec=0.5)
    assert isinstance(trs, list)
