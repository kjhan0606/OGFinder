"""observe, propagate_helio, and propagate_many on the in-tree integrator.

The position and sky limits below compare the worker with a direct call of the
same ogfmeas propagator.  They are a wiring check.  They are not a Horizons
or DE440 accuracy gate.
"""
import os
import subprocess
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"

FORCE_MODEL = "ogfmeas(builtin planets and Moon, solar Schwarzschild, Sun J2, Earth J2/J3/J4)"
_MJD = 61000.0
_POS_AU = 1e-9
_VEL = 1e-9
_SKY_DEG = 1e-6


def _state():
    from moving import kepler as K
    r, v = K.elements_to_state(2.7, 0.15, 8.0, 40.0, 20.0, 30.0, ecliptic=False)
    return np.r_[r, v]


def _independent_helio(state, mjd_ref, mjd_out):
    from moving.orbit import Propagator
    from moving.util import utc_mjd_to_tdb_jd
    prop = Propagator()
    jd0 = float(utc_mjd_to_tdb_jd([mjd_ref])[0])
    jd1 = float(utc_mjd_to_tdb_jd([mjd_out])[0])
    sun0 = prop.sun(jd0)
    bary = np.r_[np.asarray(state[:3], float) + sun0[:3], np.asarray(state[3:], float) + sun0[3:]]
    if abs(jd1 - jd0) > 1e-12:
        bary = prop.propagate(bary[None, :], jd0, np.array([jd1]))[0, 0]
    sun1 = prop.sun(jd1)
    return np.r_[bary[:3] - sun1[:3], bary[3:] - sun1[3:]]


def _independent_radec(state, mjd_ref, mjd):
    from moving.orbit import Obs, Propagator, predict_radec
    from moving.util import utc_mjd_to_tdb_jd
    prop = Propagator()
    mjd = np.asarray(mjd, float)
    jd0 = float(utc_mjd_to_tdb_jd([mjd_ref])[0])
    sun0 = prop.sun(jd0)
    bary = np.r_[np.asarray(state[:3], float) + sun0[:3], np.asarray(state[3:], float) + sun0[3:]]
    obs = Obs(mjd, np.zeros(len(mjd)), np.zeros(len(mjd)), np.ones(len(mjd)), np.ones(len(mjd)),
              ["500"] * len(mjd))
    ra, dec, _tau, _rho = predict_radec(prop, bary[None, :], jd0, obs)
    return ra[0], dec[0]


def test_import_does_not_load_external_orbit_packages():
    script = (
        "import sys\n"
        "before = set(sys.modules)\n"
        "import moving.nbody_refine, moving.arc_ranging, moving.nbody_worker\n"
        "added = set(sys.modules) - before\n"
        "banned = {'assist', 'rebound', 'spiceypy', 'neo_orbit_calculator'}\n"
        "hit = sorted(n for n in added if n.split('.')[0] in banned)\n"
        "print('\\n'.join(hit))\n"
    )
    proc = subprocess.run([sys.executable, "-c", script], cwd=_ROOT, capture_output=True, text=True, timeout=60)
    assert proc.returncode == 0, proc.stderr[-500:]
    assert proc.stdout.strip() == ""


def test_propagate_helio_matches_the_in_tree_integrator():
    from moving.nbody_worker import dispatch
    state = _state()
    out = dispatch({
        "cmd": "propagate_helio",
        "state_helio": state.tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD + 10.0,
    })
    assert out["status"] == "ok", out
    assert out["backend"] == "ogfmeas"
    assert out["force_model"] == FORCE_MODEL
    got = np.asarray(out["state_helio"], float)
    ref = _independent_helio(state, _MJD, _MJD + 10.0)
    dpos = float(np.max(np.abs(got[:3] - ref[:3])))
    dvel = float(np.max(np.abs(got[3:] - ref[3:])))
    print("OGFMEAS_MEASURE propagate_helio dpos_au=%.3e dvel=%.3e a=%.6f" % (dpos, dvel, out["elements"]["a"]))
    assert dpos < _POS_AU
    assert dvel < _VEL
    assert out["elements"]["a"] > 0
    assert "neo_orbit_calculator" not in sys.modules


def test_observe_matches_predict_radec():
    from moving.nbody_worker import dispatch
    state = _state()
    mjd = np.array([_MJD - 1.0, _MJD + 0.5, _MJD + 7.0])
    out = dispatch({
        "cmd": "observe",
        "state_helio": state.tolist(),
        "mjd_ref": _MJD,
        "mjd": mjd.tolist(),
        "code": ["500", "500", "500"],
    })
    assert out["status"] == "ok", out
    assert out["backend"] == "ogfmeas"
    assert out["force_model"] == FORCE_MODEL
    ra = np.asarray(out["ra"], float)
    dec = np.asarray(out["dec"], float)
    assert np.all(np.isfinite(ra)) and np.all(np.isfinite(dec))
    ra_ref, dec_ref = _independent_radec(state, _MJD, mjd)
    dra = (ra - ra_ref + 180.0) % 360.0 - 180.0
    ddec = dec - dec_ref
    print("OGFMEAS_MEASURE observe dra_deg=%.3e ddec_deg=%.3e" % (float(np.max(np.abs(dra))), float(np.max(np.abs(ddec)))))
    assert np.max(np.abs(dra)) < _SKY_DEG
    assert np.max(np.abs(ddec)) < _SKY_DEG


def test_propagate_many_matches_one_row_and_steps_backward():
    from moving.nbody_worker import dispatch
    state = _state()
    states = np.vstack([
        state,
        state + np.array([0.01, -0.02, 0.005, 0.0, 0.0, 0.0]),
        state + np.array([-0.01, 0.01, -0.004, 1e-4, -1e-4, 5e-5]),
    ])
    many = dispatch({
        "cmd": "propagate_many",
        "states_helio": states.tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD + 10.0,
    })
    one = dispatch({
        "cmd": "propagate_helio",
        "state_helio": states[1].tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD + 10.0,
    })
    assert many["status"] == "ok", many
    assert many["backend"] == "ogfmeas"
    assert many["force_model"] == FORCE_MODEL
    assert many["n"] == 3
    assert many["n_failed"] == 0
    assert len(many["elements"]) == 3
    assert all(row["a"] > 0 for row in many["elements"])
    for key in ("a", "e", "inc_deg", "q", "r_au"):
        assert abs(many["elements"][1][key] - one["elements"][key]) < 1e-8
    back = dispatch({
        "cmd": "propagate_helio",
        "state_helio": state.tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD - 4.0,
    })
    assert back["status"] == "ok", back
    assert back["backend"] == "ogfmeas"
    got = np.asarray(back["state_helio"], float)
    ref = _independent_helio(state, _MJD, _MJD - 4.0)
    dpos = float(np.max(np.abs(got[:3] - ref[:3])))
    print("OGFMEAS_MEASURE backward dpos_au=%.3e" % dpos)
    assert np.all(np.isfinite(got))
    assert dpos < _POS_AU
    many_back = dispatch({
        "cmd": "propagate_many",
        "states_helio": states.tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD - 4.0,
    })
    assert many_back["status"] == "ok", many_back
    assert many_back["backend"] == "ogfmeas"
    assert many_back["n"] == 3
    assert many_back["n_failed"] == 0


def test_bad_states_helio_errors():
    from moving.nbody_worker import dispatch
    out = dispatch({
        "cmd": "propagate_many",
        "states_helio": [[1.0, 2.0, 3.0, 4.0, 5.0]],
        "mjd_ref": _MJD,
        "mjd_out": _MJD + 1.0,
    })
    assert out["status"] == "error"
    assert "states_helio" in out["message"]


def test_parent_default_reports_ogfmeas():
    from moving import nbody_refine as NR
    state = _state()
    out = NR.propagate_helio(state, _MJD, _MJD + 3.0)
    assert out["backend"] == "ogfmeas"
    assert out["force_model"] == FORCE_MODEL
    ref = _independent_helio(state, _MJD, _MJD + 3.0)
    dpos = float(np.max(np.abs(np.asarray(out["state_helio"], float)[:3] - ref[:3])))
    print("OGFMEAS_MEASURE parent dpos_au=%.3e" % dpos)
    assert dpos < _POS_AU
    ra, dec, meta = NR.observe(state, _MJD, [_MJD + 1.0, _MJD + 2.0], code="500")
    assert meta["backend"] == "ogfmeas"
    assert meta["force_model"] == FORCE_MODEL
    assert np.all(np.isfinite(ra)) and np.all(np.isfinite(dec))


def test_codes_path_stays_off_when_assist_is_refused_and_no_tree_is_named():
    from moving import nbody_refine as NR
    out = NR.call_worker({
        "cmd": "propagate_helio",
        "allow_assist": False,
        "use_default": False,
        "codes_root": "",
        "state_helio": _state().tolist(),
        "mjd_ref": _MJD,
        "mjd_out": _MJD + 1.0,
    }, codes_root="", timeout=60)
    assert out["status"] == "unavailable"
    assert out.get("backend") != "ogfmeas"
