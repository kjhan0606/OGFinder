"""Element distribution for observations the night linker can already keep.

The short arc uses simplified statistical ranging.  The three-night arc uses
the local Laplace sample.  CODES propagation runs only when that integrator
is on disk.  Importing this module must not load rebound, assist, spiceypy,
or neo_orbit_calculator.
"""
import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"


def test_import_does_not_load_nbody_libraries():
    script = (
        "import sys\n"
        "before = set(sys.modules)\n"
        "import moving.nightlink, moving.arc_ranging\n"
        "added = set(sys.modules) - before\n"
        "banned = {'assist', 'rebound', 'spiceypy', 'neo_orbit_calculator'}\n"
        "hit = sorted(n for n in added if n.split('.')[0] in banned)\n"
        "print('\\n'.join(hit))\n"
    )
    proc = subprocess.run([sys.executable, "-c", script], cwd=_ROOT, capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr[-500:]
    assert proc.stdout.strip() == ""


def _truth():
    from moving import kepler as K
    r, v = K.elements_to_state(2.7, 0.15, 8.0, 40.0, 20.0, 30.0, ecliptic=False)
    return np.r_[r, v]


def _observe(times, seed, sig=0.2):
    from moving import nightlink as N
    state = _truth()
    mjd0 = 61000.0
    op, _ = N.observer_helio(times, "500")
    ra, dec, _ = N.predict(state[:3], state[3:], mjd0, times, op)
    rng = np.random.default_rng(seed)
    ra = ra + rng.normal(0.0, sig / 3600.0, len(times)) / np.cos(np.radians(dec))
    dec = dec + rng.normal(0.0, sig / 3600.0, len(times))
    return ra, dec


def _short_tracklet():
    from moving import nightlink as N
    times = 61000.0 + np.arange(3) * 0.5 / 24.0
    ra, dec = _observe(times, 3)
    return N.Tracklet(times, ra, dec, 0.2, code="500", tid=0, night=0)


def _long_tracklets():
    from moving import nightlink as N
    nights = (0, 7, 21)
    times = np.concatenate([61000.0 + nt + np.arange(3) * 0.5 / 24.0 for nt in nights])
    ra, dec = _observe(times, 3)
    trks = []
    for i, nt in enumerate(nights):
        sl = slice(i * 3, i * 3 + 3)
        trks.append(N.Tracklet(times[sl], ra[sl], dec[sl], 0.2, code="500", tid=int(nt), night=int(nt), truth=1))
    return trks


def _width(summary, key):
    q = summary["quantiles_16_50_84"][key]
    return q[2] - q[0], q


def test_one_hour_arc_is_a_broad_ranging_distribution():
    """Three observations in one hour.  The truth semimajor axis sits inside a wide cloud."""
    from moving import arc_ranging as AR
    summary, _pack = AR.distribution([_short_tracklet()], n_samples=80, seed=3)
    width, q = _width(summary, "a")
    print("RANGING_MEASURE short sampler=%s neff=%.2f rms=%.4f a=%.4f/%.4f/%.4f width=%.4f"
          % (summary["sampler"], summary["neff"], summary["best_rms_arcsec"], q[0], q[1], q[2], width))
    assert summary["status"] == "ok"
    assert summary["sampler"] == "iod.ranging"
    assert summary["neff"] > 10.0
    assert summary["best_rms_arcsec"] < 1.0
    assert width > 1.0
    assert q[0] < 2.7 < q[2]


def test_three_night_arc_is_a_narrow_laplace_distribution():
    """Nights 0, 7 and 21.  The local sample is much narrower than the one-hour cloud."""
    from moving import arc_ranging as AR
    short, _ = AR.distribution([_short_tracklet()], n_samples=80, seed=3)
    summary, _pack = AR.distribution(_long_tracklets(), n_samples=80, seed=3)
    width, q = _width(summary, "a")
    we, qe = _width(summary, "e")
    wi, qi = _width(summary, "inc_deg")
    print("RANGING_MEASURE long sampler=%s neff=%.2f rms=%.4f a=%.4f/%.4f/%.4f e50=%.4f i50=%.4f "
          "best_a=%.4f class=%s"
          % (summary["sampler"], summary["neff"], summary["best_rms_arcsec"], q[0], q[1], q[2],
             qe[1], qi[1], summary["best"]["a"], summary["class_probs"]))
    assert summary["status"] == "ok"
    assert summary["sampler"] == "laplace"
    assert summary["best_rms_arcsec"] < 0.5
    assert abs(summary["best"]["a"] - 2.7) < 0.02
    assert abs(q[1] - 2.7) < 0.02
    assert abs(qe[1] - 0.15) < 0.01
    assert abs(qi[1] - 8.0) < 0.05
    assert width < 0.05
    assert _width(short, "a")[0] > 10.0 * width
    assert summary["class_probs"].get("MBA", 0.0) > 0.9


def test_attach_leaves_the_two_body_point_fit():
    from moving import arc_ranging as AR
    from moving import nightlink as N
    trks = _long_tracklets()
    res = N.link_nights(trks, max_gap_days=30)
    assert len(res["groups"]) == 1
    before = (res["groups"][0]["elements"]["a"], res["groups"][0]["chi2_red"], res["groups"][0]["rms_arcsec"])
    out = AR.attach(trks, res, n_samples=40, seed=3, n_propagate=0)
    g = out["groups"][0]
    assert (g["elements"]["a"], g["chi2_red"], g["rms_arcsec"]) == before
    assert g["ranging"]["status"] == "ok"
    assert "quantiles_16_50_84" in g["ranging"]
    assert g["ranging"]["codes"]["status"] == "skipped"
    assert out["ranging"]["n"] == 1


def test_missing_backend_keeps_the_distribution():
    from moving import arc_ranging as AR
    from moving import nightlink as N
    trks = _long_tracklets()
    res = N.link_nights(trks, max_gap_days=30)
    before = res["groups"][0]["elements"]["a"]
    out = AR.attach(trks, res, n_samples=40, seed=3, n_propagate=4, codes_root="", use_default=False, timeout=60)
    g = out["groups"][0]
    assert g["elements"]["a"] == before
    assert g["ranging"]["status"] == "ok"
    assert g["ranging"]["codes"]["status"] == "unavailable"
    assert g["ranging"]["quantiles_16_50_84"]["a"][0] < g["ranging"]["quantiles_16_50_84"]["a"][2]


def _codes_ready():
    from moving import nbody_refine as NR
    info = NR.probe(allow_assist=False)
    return info.get("status") == "ok" and info.get("backend") == "codes"


def test_codes_propagates_accepted_samples():
    if not _codes_ready():
        pytest.skip("CODES integrator not available (set OGF_CODES_ROOT)")
    from moving import arc_ranging as AR
    from moving import nightlink as N
    trks = _long_tracklets()
    res = N.link_nights(trks, max_gap_days=30)
    out = AR.attach(trks, res, n_samples=40, seed=3, n_propagate=6, timeout=180)
    g = out["groups"][0]
    codes = g["ranging"]["codes"]
    print("RANGING_MEASURE codes status=%s backend=%s n=%s epoch=%s a=%s"
          % (codes.get("status"), codes.get("backend"), codes.get("n_propagated"),
             codes.get("mjd_epoch"), (codes.get("quantiles_16_50_84") or {}).get("a")))
    assert codes["status"] == "ok"
    assert codes["backend"] == "codes"
    assert codes["n_propagated"] == 6
    qa = g["ranging"]["quantiles_16_50_84"]["a"]
    pa = codes["quantiles_16_50_84"]["a"]
    assert abs(pa[1] - qa[1]) < 0.05
