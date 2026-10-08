"""N-body refinement of nightlink groups.

The import test always runs: loading the linker must not load rebound, assist,
spiceypy, or the external CODES package.  The orbit-recovery tests run only
when the worker reports a backend (ASSIST, or CODES under OGF_CODES_ROOT /
~/BACKUP/3.5ST).  They generate the observations with that same N-body model.
"""
import os
import subprocess
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _codes_ready():
    if getattr(_codes_ready, "cached", None) is None:
        from moving import nbody_refine as NR
        info = NR.probe(allow_assist=False)
        _codes_ready.cached = info.get("status") == "ok" and info.get("backend") == "codes"
        _codes_ready.info = info
    return _codes_ready.cached


def _need_codes():
    if not _codes_ready():
        pytest.skip("CODES integrator not available (set OGF_CODES_ROOT). ASSIST recovery stays in test_orbit.py")


def test_linker_import_does_not_load_nbody_libraries():
    """A fresh interpreter: importing the linker and the refiner must not load the N-body libraries."""
    script = (
        "import sys\n"
        "before = set(sys.modules)\n"
        "import moving.nightlink, moving.nbody_refine\n"
        "added = set(sys.modules) - before\n"
        "banned = {'assist', 'rebound', 'spiceypy', 'neo_orbit_calculator'}\n"
        "hit = sorted(n for n in added if n.split('.')[0] in banned)\n"
        "print('\\n'.join(hit))\n"
    )
    proc = subprocess.run([sys.executable, "-c", script], cwd=_ROOT, capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr[-500:]
    assert proc.stdout.strip() == ""


def test_missing_backend_keeps_the_two_body_groups():
    from moving import nbody_refine as NR
    from moving import nightlink as N
    rng = np.random.default_rng(1)
    from moving import nightlink_sim as S
    trks = S.make_field(2, [0, 3, 7], rng, sig=0.2, ra_window=(100, 101), dec_window=(0, 1))
    res = N.link_nights(trks, max_gap_days=20)
    assert len(res["groups"]) >= 1
    ids_before = [list(g["ids"]) for g in res["groups"]]
    chi_before = [g["chi2_red"] for g in res["groups"]]
    out = NR.refine_result(trks, res, codes_root="", allow_assist=False, use_default=False, timeout=60)
    assert out["nbody"]["status"] == "unavailable"
    assert [list(g["ids"]) for g in out["groups"]] == ids_before
    assert [g["chi2_red"] for g in out["groups"]] == chi_before


def _mba_state():
    from moving import kepler as K
    # Ecliptic elements; equatorial vectors for the CODES ICRF frame.
    r, v = K.elements_to_state(2.7, 0.15, 8.0, 40.0, 20.0, 30.0, ecliptic=False)
    return np.r_[r, v]


def _tracklets(mjd0, nights, ra, dec, sig, rng):
    from moving import nightlink as N
    trks = []
    cursor = 0
    for nt in nights:
        n = 3
        sl = slice(cursor, cursor + n)
        noise_ra = rng.normal(0, sig / 3600.0, n) / np.cos(np.radians(dec[sl]))
        noise_de = rng.normal(0, sig / 3600.0, n)
        times = mjd0 + nt + np.arange(n) * 0.5 / 24.0
        trks.append(N.Tracklet(times, ra[sl] + noise_ra, dec[sl] + noise_de, sig, code="500", tid=int(nt), night=int(nt), truth=1))
        cursor += n
    return trks


def _recover(span_nights, sig=0.2, seed=7):
    from moving import nbody_refine as NR
    from moving import nightlink as N
    mjd0 = 61000.0
    state = _mba_state()
    times = np.concatenate([mjd0 + nt + np.arange(3) * 0.5 / 24.0 for nt in span_nights])
    ra, dec, _meta = NR.observe(state, mjd0, times, code="500")
    rng = np.random.default_rng(seed)
    trks = _tracklets(mjd0, span_nights, ra, dec, sig, rng)
    gap = float(span_nights[-1] - span_nights[0] + 5)
    res = N.link_nights(trks, max_gap_days=gap)
    raw = None
    if len(res["groups"]) == 1:
        raw = res["groups"][0].get("rms_arcsec")
    out = NR.refine_result(trks, res, allow_assist=False, chi2_max=4.0, resid_max=6.0, timeout=300)
    return trks, out, raw, state, mjd0


def _assert_recovered(nights, seed):
    """Same CODES truth particle.  Bounds were fixed before the first passing run."""
    _need_codes()
    from moving import kepler as K
    from moving import nbody_refine as NR
    _trks, out, _raw, state, mjd0 = _recover(nights, sig=0.2, seed=seed)
    assert out["nbody"]["status"] == "ok", out["nbody"]
    assert out["nbody"]["backend"] == "codes"
    assert len(out["groups"]) == 1, out.get("nbody_rejected")
    g = out["groups"][0]
    assert g["nbody"]["accepted"] is True
    assert g["chi2_red"] < 4.0
    assert g["rms_arcsec"] < 0.5
    assert g["rms_arcsec"] <= g["two_body_raw_rms_arcsec"] + 0.05
    truth = NR.propagate_helio(state, mjd0, g["nbody_mjd"])
    el_t = truth["elements"]
    el = g["elements"]
    assert abs(el["a"] - el_t["a"]) / el_t["a"] < 0.02, (el, el_t, g["nbody"])
    assert abs(el["e"] - el_t["e"]) < 0.01, (el, el_t)
    assert abs(el["inc_deg"] - el_t["inc_deg"]) < 0.2, (el, el_t)
    fit = g["fit"]
    r, v = K.propagate_2body(fit.state[0], fit.state[1], g["nbody_mjd"] - fit.t_ref)
    el2 = K.state_to_elements(r, v)
    print("NBODY_MEASURE nights=%s seed=%s status=%s iters=%s two_body_raw_rms=%.4f nbody_rms=%.4f "
          "chi2_red=%.3f da_nbody=%.6f da_twobody=%.6f de_nbody=%.6f de_twobody=%.6f "
          "di_nbody=%.6f di_twobody=%.6f"
          % (list(nights), seed, g["nbody"]["status"], g["nbody"].get("iterations"),
             g["two_body_raw_rms_arcsec"], g["rms_arcsec"], g["chi2_red"],
             el["a"] - el_t["a"], el2["a"] - el_t["a"],
             el["e"] - el_t["e"], el2["e"] - el_t["e"],
             el["inc_deg"] - el_t["inc_deg"], el2["i"] - el_t["inc_deg"]))


def test_codes_refit_recovers_a_21_day_main_belt_arc():
    """Three nights over 21 days, 0.2 arcsec noise, seed 7."""
    _assert_recovered([0, 7, 21], 7)


def test_codes_refit_recovers_a_45_day_main_belt_arc():
    """Nights 0, 20 and 45, seed 9."""
    _assert_recovered([0, 20, 45], 9)


def test_codes_refit_recovers_a_120_day_main_belt_arc():
    """Nights 0, 40 and 120, seed 11.  In-arc residual stays at the noise; this checks the elements."""
    _assert_recovered([0, 40, 120], 11)
