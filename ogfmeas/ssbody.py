"""Test-particle motion in the solar system.

Positions of the Sun, Moon and planets come from astropy's built-in
ephemeris. The force on a massless particle is Newtonian gravity of those
bodies, the solar Schwarzschild term, the Sun's J2, and the Earth's J2, J3
and J4. Masses and zonal coefficients are the JPL DE440 header constants
(Park et al. 2021, AJ 161, 105, GROUP 1041). The integrator is the
Dormand-Prince 5(4) pair (Dormand & Prince 1980).

This does not read a DE440 or asteroid kernel, does not include the 16
massive asteroids, and does not implement the full Einstein-Infeld-Hoffmann
equations. It does not import rebound or assist.
"""
import numpy as np
from scipy.interpolate import CubicSpline

# JPL DE440 header, GROUP 1041. GM in AU^3 day^-2, radii in km, AU in km.
_AU_KM = 149597870.7
_CLIGHT_KM_S = 299792.458
_C_AUD = _CLIGHT_KM_S * 86400.0 / _AU_KM
_EMRAT = 81.30056822149722
_GMB = 8.997011392947347e-10
_GM_EARTH = _GMB * _EMRAT / (1.0 + _EMRAT)
_GM_MOON = _GMB / (1.0 + _EMRAT)
_GM = {
    "sun": 0.00029591220828411956,
    "mercury": 4.912500194889318e-11,
    "venus": 7.243452332644119e-10,
    "earth": _GM_EARTH,
    "moon": _GM_MOON,
    "mars": 9.549548829725812e-11,
    "jupiter": 2.825345825225792e-07,
    "saturn": 8.45970599337629e-08,
    "uranus": 1.29202656496824e-08,
    "neptune": 1.524357347885194e-08,
}
_NAMES = tuple(_GM)
# IAU solar north pole (Archinal et al. 2018, Celest. Mech. Dyn. Astr. 130, 22).
_SUN_POLE_RA = np.radians(286.13)
_SUN_POLE_DEC = np.radians(63.87)
_SUN_POLE = np.array([
    np.cos(_SUN_POLE_DEC) * np.cos(_SUN_POLE_RA),
    np.cos(_SUN_POLE_DEC) * np.sin(_SUN_POLE_RA),
    np.sin(_SUN_POLE_DEC),
])
_ZONAL = {
    "sun": (696000.0 / _AU_KM, 2.1961391516529825e-07, 0.0, 0.0, _SUN_POLE),
    "earth": (6378.1366 / _AU_KM, 0.00108262539, -2.53241e-06, -1.619898e-06, np.array([0.0, 0.0, 1.0])),
}

FORCE_MODEL = "ogfmeas(builtin planets and Moon, solar Schwarzschild, Sun J2, Earth J2/J3/J4)"

# Dormand & Prince 1980, J. Comput. Appl. Math. 6, 19.
_C = np.array([0.0, 0.2, 0.3, 0.8, 8.0 / 9.0, 1.0, 1.0])
_A = [
    [],
    [0.2],
    [3.0 / 40.0, 9.0 / 40.0],
    [44.0 / 45.0, -56.0 / 15.0, 32.0 / 9.0],
    [19372.0 / 6561.0, -25360.0 / 2187.0, 64448.0 / 6561.0, -212.0 / 729.0],
    [9017.0 / 3168.0, -355.0 / 33.0, 46732.0 / 5247.0, 49.0 / 176.0, -5103.0 / 18656.0],
    [35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0],
]
_B5 = np.array([35.0 / 384.0, 0.0, 500.0 / 1113.0, 125.0 / 192.0, -2187.0 / 6784.0, 11.0 / 84.0, 0.0])
_B4 = np.array([5179.0 / 57600.0, 0.0, 7571.0 / 16695.0, 393.0 / 640.0, -92097.0 / 339200.0, 187.0 / 2100.0, 1.0 / 40.0])
_E = _B5 - _B4


def gm_of(name):
    return _GM[name]


def schwarzschild(rel_sun, vel, mu, c_aud=_C_AUD):
    """Solar Schwarzschild acceleration. ``rel_sun`` and ``vel`` are (n, 3)."""
    r2 = np.sum(rel_sun * rel_sun, axis=1)
    r = np.sqrt(r2)
    v2 = np.sum(vel * vel, axis=1)
    rv = np.sum(rel_sun * vel, axis=1)
    coef = mu / (c_aud * c_aud * r2 * r)
    return coef[:, None] * ((4.0 * mu / r - v2)[:, None] * rel_sun + (4.0 * rv)[:, None] * vel)


def _legendre_potential(rel, gm, radius, j2, j3, j4):
    r = np.linalg.norm(rel, axis=-1)
    s = rel[..., 2] / r
    q = radius / r
    p2 = 0.5 * (3.0 * s * s - 1.0)
    p3 = 0.5 * (5.0 * s * s * s - 3.0 * s)
    p4 = (35.0 * s ** 4 - 30.0 * s * s + 3.0) / 8.0
    return -(gm / r) * (j2 * q ** 2 * p2 + j3 * q ** 3 * p3 + j4 * q ** 4 * p4)


def _pole_frame(pole):
    pole = np.asarray(pole, dtype=np.float64)
    pole = pole / np.linalg.norm(pole)
    helper = np.array([1.0, 0.0, 0.0]) if abs(pole[2]) > 0.9 else np.array([0.0, 0.0, 1.0])
    x_axis = np.cross(helper, pole)
    x_axis /= np.linalg.norm(x_axis)
    y_axis = np.cross(pole, x_axis)
    return np.vstack([x_axis, y_axis, pole])


def zonal_acceleration(rel, gm, radius, j2, j3, j4):
    """Gradient of the J2/J3/J4 force function. ``rel`` is (n, 3) with z along the pole."""
    rel = np.asarray(rel, dtype=np.float64)
    h = np.maximum(np.linalg.norm(rel, axis=1), 1e-8) * 1e-6
    acc = np.empty_like(rel)
    for k in range(3):
        step = np.zeros_like(rel)
        step[:, k] = h
        acc[:, k] = (_legendre_potential(rel + step, gm, radius, j2, j3, j4)
                     - _legendre_potential(rel - step, gm, radius, j2, j3, j4)) / (2.0 * h)
    return acc


def point_mass_acceleration(pos, body_pos, gm):
    """Newtonian acceleration at ``pos`` (n, 3) from bodies at ``body_pos`` (nb, 3)."""
    delta = body_pos[None, :, :] - pos[:, None, :]
    dist2 = np.sum(delta * delta, axis=2)
    dist = np.sqrt(dist2)
    factor = gm[None, :] / (dist2 * dist + 1e-30)
    return np.sum(factor[:, :, None] * delta, axis=1)


def field(pos, vel, body_pos, gm, sun_index=0, zonals=None, gr=True):
    """Acceleration of massless particles. ``zonals`` maps a body index to (radius, j2, j3, j4, pole)."""
    acc = point_mass_acceleration(pos, body_pos, gm)
    if gr:
        rel = pos - body_pos[sun_index]
        acc += schwarzschild(rel, vel, gm[sun_index])
    if zonals:
        for index, (radius, j2, j3, j4, pole) in zonals.items():
            rel = pos - body_pos[index]
            frame = _pole_frame(pole)
            local = rel @ frame.T
            acc += zonal_acceleration(local, gm[index], radius, j2, j3, j4) @ frame
    return acc


def integrate_ivp(fun, y0, t0, times, rtol=1e-11, atol=1e-12, max_step=2.0):
    """Integrate dy/dt = fun(t, y) with Dormand-Prince 5(4). ``y0`` is (n, d).

    ``times`` is one direction relative to ``t0``. Returns y at each time, shape (m, n, d).
    """
    y = np.array(y0, dtype=np.float64, copy=True)
    times = np.asarray(times, dtype=np.float64)
    out = np.empty((times.size,) + y.shape, dtype=np.float64)
    if times.size == 0:
        return out
    direction = 1.0 if times[0] >= t0 else -1.0
    if np.any((times - t0) * direction < -1e-15):
        raise ValueError("times must lie on one side of t0")
    t = float(t0)
    h = direction * min(0.05, max_step)
    atol_v = np.array(atol, dtype=np.float64)
    cursor = 0
    steps = 0
    while cursor < times.size:
        target = float(times[cursor])
        if abs(target - t) <= 1e-14:
            out[cursor] = y
            cursor += 1
            continue
        remain = target - t
        h_try = remain if abs(remain) < abs(h) else h
        accepted = False
        for _reject in range(40):
            stages = [fun(t, y)]
            for i in range(1, 7):
                increment = np.zeros_like(y)
                for j, aij in enumerate(_A[i]):
                    increment += aij * stages[j]
                stages.append(fun(t + _C[i] * h_try, y + h_try * increment))
            err_vec = h_try * sum(_E[j] * stages[j] for j in range(7))
            y5 = y + h_try * sum(_B5[j] * stages[j] for j in range(7))
            scale = atol_v + rtol * np.maximum(np.abs(y), np.abs(y5))
            err = float(np.max(np.abs(err_vec) / scale))
            if not np.isfinite(err):
                raise RuntimeError("test-particle integrator diverged")
            factor = 5.0 if err == 0.0 else min(5.0, max(0.2, 0.9 * err ** -0.2))
            h_next = np.sign(h_try) * min(max_step, abs(h_try) * factor)
            if err <= 1.0:
                t = t + float(h_try)
                y = y5
                h = h_next if abs(h_next) > 1e-12 else direction * 1e-12
                accepted = True
                steps += 1
                break
            h_try = h_next if abs(h_next) > 1e-14 else h_try * 0.2
        if not accepted:
            raise RuntimeError("test-particle integrator rejected every trial step")
        if steps > 200000:
            raise RuntimeError("test-particle integrator step limit")
        if abs(t - target) <= 1e-12:
            out[cursor] = y
            cursor += 1
    return out


class BuiltinEphemeris:
    """Barycentric states from astropy's built-in ephemeris, cached on a cubic grid."""

    def __init__(self, names=_NAMES, step_day=0.5):
        self.names = tuple(names)
        self.index = {name: i for i, name in enumerate(self.names)}
        self.gm = np.array([_GM[name] for name in self.names], dtype=np.float64)
        self.step = float(step_day)
        self.jd = None
        self._pos = None
        self._vel = None
        self.zonals = {}
        for name, spec in _ZONAL.items():
            if name in self.index:
                self.zonals[self.index[name]] = spec

    def ensure(self, jd_lo, jd_hi):
        lo = float(min(jd_lo, jd_hi)) - 1.0
        hi = float(max(jd_lo, jd_hi)) + 1.0
        if self.jd is not None and self.jd[0] <= lo and self.jd[-1] >= hi:
            return
        if self.jd is not None:
            lo = min(lo, float(self.jd[0]))
            hi = max(hi, float(self.jd[-1]))
        grid = np.arange(lo, hi + self.step * 0.5, self.step)
        from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
        from astropy.time import Time
        when = Time(grid, format="jd", scale="tdb")
        pos = np.empty((grid.size, len(self.names), 3), dtype=np.float64)
        vel = np.empty_like(pos)
        with solar_system_ephemeris.set("builtin"):
            for i, name in enumerate(self.names):
                p, v = get_body_barycentric_posvel(name, when)
                pos[:, i, :] = np.stack([c.to_value("AU") for c in p.xyz], axis=1)
                vel[:, i, :] = np.stack([c.to_value("AU/d") for c in v.xyz], axis=1)
        self.jd = grid
        self._pos = [CubicSpline(grid, pos[:, i, :], extrapolate=False) for i in range(len(self.names))]
        self._vel = [CubicSpline(grid, vel[:, i, :], extrapolate=False) for i in range(len(self.names))]

    def state(self, name, jd):
        self.ensure(jd, jd)
        i = self.index[name]
        t = float(jd)
        return np.array([self._pos[i](t)[k] for k in range(3)] + [self._vel[i](t)[k] for k in range(3)], dtype=np.float64)

    def positions(self, jd):
        self.ensure(jd, jd)
        t = float(jd)
        return np.vstack([spline(t) for spline in self._pos])

    def propagate(self, states, jd0, jds, gr=True, zonals=True):
        states = np.atleast_2d(np.asarray(states, dtype=np.float64))
        jds = np.atleast_1d(np.asarray(jds, dtype=np.float64))
        jd0 = float(jd0)
        if jds.size:
            self.ensure(min(jd0, float(np.min(jds))), max(jd0, float(np.max(jds))))
        gm = self.gm
        extra = self.zonals if zonals else {}
        sun = self.index["sun"]

        def fun(_t, y):
            pos = y[:, 0:3]
            vel = y[:, 3:6]
            acc = field(pos, vel, self.positions(_t), gm, sun_index=sun, zonals=extra, gr=gr)
            return np.concatenate([vel, acc], axis=1)

        order = np.argsort(jds)
        out = np.zeros((len(states), len(jds), 6), dtype=np.float64)
        forward = [k for k in order if jds[k] >= jd0]
        backward = [k for k in order[::-1] if jds[k] < jd0]
        for seq in (forward, backward):
            if not seq:
                continue
            targets = jds[np.array(seq)]
            got = integrate_ivp(fun, states, jd0, targets)
            for row, k in enumerate(seq):
                out[:, k, :] = got[row]
        return out
