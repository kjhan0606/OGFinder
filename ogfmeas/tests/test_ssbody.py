"""In-tree solar-system integrator. Gates were fixed before the run."""
import sys

import numpy as np

import ogfmeas.ssbody as S


def test_imports_do_not_load_rebound_or_assist():
    assert "rebound" not in sys.modules
    assert "assist" not in sys.modules


def test_exponential_converges():
    def fun(_t, y):
        return y

    y = S.integrate_ivp(fun, np.array([[1.0]]), 0.0, np.array([1.0]), rtol=1e-12, atol=1e-14, max_step=0.2)
    assert abs(float(y[0, 0, 0]) - np.e) < 1e-10


def test_sun_only_matches_two_body():
    from moving import kepler as K
    mu = S.gm_of("sun")
    r, v = K.elements_to_state(2.7, 0.15, 8.0, 40.0, 20.0, 30.0, ecliptic=False, mu=mu)
    state = np.r_[r, v]
    jd0 = 2459000.5
    got = S.integrate_ivp(
        lambda _t, y: np.concatenate([y[:, 3:6], S.point_mass_acceleration(y[:, 0:3], np.zeros((1, 3)), np.array([mu]))], axis=1),
        state[None, :], jd0, np.array([jd0 + 100.0]), rtol=1e-12, atol=1e-14, max_step=2.0,
    )[0, 0]
    r2, v2 = K.propagate_2body(r, v, 100.0, mu=mu)
    assert np.linalg.norm(got[:3] - r2) < 1e-8
    assert np.linalg.norm(got[3:] - v2) < 1e-10


def test_round_trip_returns_to_the_start():
    mu = S.gm_of("sun")
    state = np.array([1.5, 0.2, 0.0, 0.0, np.sqrt(mu / 1.5), 0.0])

    def fun(_t, y):
        acc = S.point_mass_acceleration(y[:, 0:3], np.zeros((1, 3)), np.array([mu]))
        return np.concatenate([y[:, 3:6], acc], axis=1)

    mid = S.integrate_ivp(fun, state[None, :], 0.0, np.array([40.0]), rtol=1e-12, atol=1e-14)[0, 0]
    back = S.integrate_ivp(fun, mid[None, :], 40.0, np.array([0.0]), rtol=1e-12, atol=1e-14)[0, 0]
    assert np.linalg.norm(back[:3] - state[:3]) < 1e-9


def test_equatorial_j2_is_more_attractive():
    gm = S.gm_of("earth")
    radius = 6378.1366 / 149597870.7
    j2 = 0.00108262539
    r = 10.0 * radius
    acc = S.zonal_acceleration(np.array([[r, 0.0, 0.0]]), gm, radius, j2, 0.0, 0.0)[0]
    expected = -1.5 * j2 * gm * radius ** 2 / r ** 4
    assert abs(acc[0] - expected) / abs(expected) < 1e-4
    assert abs(acc[1]) < 1e-6 * abs(expected)
    assert abs(acc[2]) < 1e-6 * abs(expected)


def test_builtin_earth_matches_astropy_on_the_grid():
    planets = S.BuiltinEphemeris()
    jd = 2459000.5
    planets.ensure(jd, jd)
    got = planets.state("earth", jd)
    from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
    from astropy.time import Time
    with solar_system_ephemeris.set("builtin"):
        p, v = get_body_barycentric_posvel("earth", Time(jd, format="jd", scale="tdb"))
    direct = np.array([c.to_value("AU") for c in p.xyz] + [c.to_value("AU/d") for c in v.xyz])
    assert np.linalg.norm(got[:3] - direct[:3]) < 1e-10
    assert np.linalg.norm(got[3:] - direct[3:]) < 1e-12


def test_propagator_does_not_import_the_gpl_packages():
    for name in ("rebound", "assist"):
        sys.modules.pop(name, None)
    from moving import orbit as O
    prop = O.Propagator()
    assert prop.force_model.startswith("ogfmeas")
    assert O.assist_available() is True
    assert "rebound" not in sys.modules
    assert "assist" not in sys.modules
    sun = prop.sun(2459000.5)
    assert sun.shape == (6,)
    assert np.isfinite(sun).all()
