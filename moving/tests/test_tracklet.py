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
    trs = T.link_exposures(dets, tol_arcsec=0.4, min_exposures=4, obs_off_au={i: off[i] for i in range(n)}, min_disp_arcsec=0.5, flux_tol_dex=None)
    best = trs[0]
    assert set(best["members"]) == {0, 1, 2, 3}
    assert abs(best["rate_ash"] - np.hypot(*mu)) < 0.5
