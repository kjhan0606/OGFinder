"""Orbit determination with the in-tree solar-system integrator.

The Horizons position gate stays 5e-6 AU after 100 days. That check uses the
built-in planet ephemeris, not a DE440 kernel and not the 16 asteroids.
"""
import numpy as np
from conftest import needs_ephem, needs_net


@needs_ephem
@needs_net
def test_propagator_matches_horizons():
    from moving import orbit as O, horizons as H, kepler as K, util
    p = O.Propagator()
    js = H.sbdb("433"); el, _ = H.sbdb_elements(js)
    jd0 = el["epoch_jd"]
    r, v = K.elements_to_state(el["a"], el["e"], el["i"], el["om"], el["w"], el["ma"])
    sun = p.sun(jd0)
    st = np.r_[K.R_ECL2EQ @ r + sun[:3], K.R_ECL2EQ @ v + sun[3:]]
    jd1 = jd0 + 100.0
    S = p.propagate(st[None, :], jd0, [jd1])[0, 0]
    hz = H.vectors("433;", [jd1], center="500@0")[0]
    assert np.linalg.norm(S[:3] - hz[:3]) < 5e-6


@needs_ephem
def test_differential_correction_synthetic():
    """Generate observations from a known orbit (ASSIST), perturb the initial state, recover it."""
    from moving import orbit as O, kepler as K
    p = O.Propagator()
    jd0 = 2459000.5
    r, v = K.elements_to_state(2.7, 0.15, 8.0, 100.0, 50.0, 30.0)
    sun = p.sun(jd0)
    st = np.r_[K.R_ECL2EQ @ r + sun[:3], K.R_ECL2EQ @ v + sun[3:]]
    mjd = 59000.0 + np.linspace(0, 40, 14) - 0.0
    obs = O.Obs(mjd, np.zeros(14), np.zeros(14), 0.2, 0.2, ["G96"] * 14)
    ra, dec, _, _ = O.predict_radec(p, st[None, :], jd0, obs)
    rng = np.random.default_rng(1)
    obs.ra = ra[0] + rng.normal(0, 0.2 / 3600, 14) / np.cos(np.radians(dec[0])); obs.dec = dec[0] + rng.normal(0, 0.2 / 3600, 14)
    start = st + np.r_[rng.normal(0, 2e-5, 3), rng.normal(0, 2e-7, 3)]
    res = O.differential_correction(p, start, jd0, obs, reject=False)
    assert res["chi2_red"] < 2.5
    sig = np.sqrt(np.diag(res["cov"]))
    pull = (res["state"] - st) / sig
    assert np.all(np.abs(pull) < 4.0), pull
