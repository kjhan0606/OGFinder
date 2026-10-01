"""Orbital mechanics helpers: element <-> state conversions, universal-variable two-body
propagation (vectorised), ecliptic/equatorial rotation.  Units: AU, day, degrees for angles."""
import numpy as np

AU_M = 149597870700.0
GM_SUN_M3S2 = 1.32712440041279419e20          # DE440
MU = GM_SUN_M3S2 * 86400.0 ** 2 / AU_M ** 3     # AU^3/d^2
EPS_J2000 = np.radians(84381.448 / 3600.0)      # IAU 1976 obliquity (as used by JPL/MPC ecliptic frame)
_c, _s = np.cos(EPS_J2000), np.sin(EPS_J2000)
R_ECL2EQ = np.array([[1, 0, 0], [0, _c, -_s], [0, _s, _c]])
R_EQ2ECL = R_ECL2EQ.T


def stumpff(z):
    z = np.asarray(z, dtype=float)
    C = np.empty_like(z); S = np.empty_like(z)
    small = np.abs(z) < 1e-6
    pos = z > 1e-6
    neg = z < -1e-6
    zs = z[small]
    C[small] = 0.5 - zs / 24 + zs ** 2 / 720
    S[small] = 1 / 6 - zs / 120 + zs ** 2 / 5040
    zp = z[pos]; sq = np.sqrt(zp)
    C[pos] = (1 - np.cos(sq)) / zp
    S[pos] = (sq - np.sin(sq)) / sq ** 3
    zn = z[neg]; sq = np.sqrt(-zn)
    C[neg] = (np.cosh(sq) - 1) / (-zn)
    S[neg] = (np.sinh(sq) - sq) / sq ** 3
    return C, S


def propagate_2body(r0, v0, dt, mu=MU):
    """Vectorised universal-variable Kepler propagation.
    r0, v0: (...,3); dt: broadcastable scalar/array (...,).  Returns r, v (...,3)."""
    r0 = np.asarray(r0, float); v0 = np.asarray(v0, float)
    dt = np.asarray(dt, float) * np.ones(r0.shape[:-1])
    rm = np.linalg.norm(r0, axis=-1)
    vr = np.sum(r0 * v0, axis=-1) / rm
    alpha = 2.0 / rm - np.sum(v0 * v0, axis=-1) / mu
    sm = np.sqrt(mu)
    chi = sm * np.abs(alpha) * dt
    # initial guess (Vallado): elliptic/hyperbolic
    chi = np.where(np.abs(alpha) > 1e-12, chi, sm * dt / rm)
    for _ in range(60):
        z = alpha * chi ** 2
        C, S = stumpff(z)
        F = rm * vr / sm * chi ** 2 * C + (1 - alpha * rm) * chi ** 3 * S + rm * chi - sm * dt
        dF = rm * vr / sm * chi * (1 - alpha * chi ** 2 * S) + (1 - alpha * rm) * chi ** 2 * C + rm
        dchi = F / dF
        chi = chi - dchi
        if np.nanmax(np.abs(dchi)) < 1e-12:
            break
    z = alpha * chi ** 2
    C, S = stumpff(z)
    f = 1 - chi ** 2 / rm * C
    g = dt - chi ** 3 / sm * S
    r = f[..., None] * r0 + g[..., None] * v0
    rn = np.linalg.norm(r, axis=-1)
    fd = sm / (rn * rm) * chi * (alpha * chi ** 2 * S - 1)
    gd = 1 - chi ** 2 / rn * C
    v = fd[..., None] * r0 + gd[..., None] * v0
    return r, v


def elements_to_state(a, e, i, om, w, M, ecliptic=True, mu=MU):
    """Heliocentric Keplerian (a [AU], e, angles [deg], M mean anomaly [deg]) -> state.
    Returns (r, v) heliocentric in ecliptic J2000 if ecliptic=True (no rotation) else rotated to equatorial."""
    i, om, w, M = [np.radians(x) for x in (i, om, w, M)]
    # solve Kepler's equation
    if e < 1:
        E = M + e * np.sin(M)
        for _ in range(100):
            dE = (E - e * np.sin(E) - M) / (1 - e * np.cos(E))
            E -= dE
            if abs(dE) < 1e-14:
                break
        x = a * (np.cos(E) - e); y = a * np.sqrt(1 - e * e) * np.sin(E)
        n = np.sqrt(mu / a ** 3)
        vx = -a * n * np.sin(E) / (1 - e * np.cos(E)); vy = a * n * np.sqrt(1 - e * e) * np.cos(E) / (1 - e * np.cos(E))
    else:
        H = np.arcsinh(M / e) if e > 1 else M
        for _ in range(100):
            dH = (e * np.sinh(H) - H - M) / (e * np.cosh(H) - 1)
            H -= dH
            if abs(dH) < 1e-14:
                break
        x = a * (np.cosh(H) - e); y = -a * np.sqrt(e * e - 1) * np.sinh(H)
        n = np.sqrt(mu / (-a) ** 3)
        vx = -a * n * np.sinh(H) / (e * np.cosh(H) - 1)
        vy = -a * n * np.sqrt(e * e - 1) * np.cosh(H) / (e * np.cosh(H) - 1)
    cO, sO, cw, sw, ci, si = np.cos(om), np.sin(om), np.cos(w), np.sin(w), np.cos(i), np.sin(i)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci, sO * si],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci, -cO * si],
                  [sw * si, cw * si, ci]])
    r = R @ np.array([x, y, 0.0]); v = R @ np.array([vx, vy, 0.0])
    if not ecliptic:
        r = R_ECL2EQ @ r; v = R_ECL2EQ @ v
    return r, v


def state_to_elements(r, v, mu=MU, ecliptic_input=False):
    """Heliocentric state -> dict(a,e,i,om,w,M,q,Q,P,n,tp_offset) in AU/deg/day.
    If ecliptic_input False the vectors are taken as equatorial and rotated to ecliptic."""
    r = np.asarray(r, float); v = np.asarray(v, float)
    if not ecliptic_input:
        r = R_EQ2ECL @ r; v = R_EQ2ECL @ v
    rn = np.linalg.norm(r); v2 = v @ v
    h = np.cross(r, v); hn = np.linalg.norm(h)
    ev = np.cross(v, h) / mu - r / rn
    e = np.linalg.norm(ev)
    en = v2 / 2 - mu / rn
    a = -mu / (2 * en)
    i = np.arccos(h[2] / hn)
    nvec = np.array([-h[1], h[0], 0.0]); nn = np.linalg.norm(nvec)
    if nn > 1e-14:
        om = np.arctan2(nvec[1], nvec[0]) % (2 * np.pi)
        w = np.arccos(np.clip(nvec @ ev / (nn * e), -1, 1))
        if ev[2] < 0:
            w = 2 * np.pi - w
    else:
        om = 0.0
        w = np.arctan2(ev[1], ev[0]) % (2 * np.pi)
    nu = np.arccos(np.clip(ev @ r / (e * rn), -1, 1))
    if r @ v < 0:
        nu = 2 * np.pi - nu
    if e < 1:
        E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2), np.sqrt(1 + e) * np.cos(nu / 2))
        M = (E - e * np.sin(E)) % (2 * np.pi)
        n = np.sqrt(mu / a ** 3)
    else:
        F = 2 * np.arctanh(np.sqrt((e - 1) / (e + 1)) * np.tan(nu / 2))
        M = e * np.sinh(F) - F
        n = np.sqrt(mu / (-a) ** 3)
    q = a * (1 - e)
    out = dict(a=a, e=e, i=np.degrees(i), om=np.degrees(om), w=np.degrees(w),
               ma=np.degrees(M) if e < 1 else M, q=q, n=np.degrees(n))
    out["Q"] = a * (1 + e) if e < 1 else np.inf
    out["P_yr"] = 2 * np.pi / n / 365.25 if e < 1 else np.inf
    return out


def orbit_class(a, e, i=None):
    """Coarse dynamical class from heliocentric a, e (used for IOD posterior class probabilities)."""
    if e >= 1 or a <= 0:
        return "HYP"
    q = a * (1 - e)
    if q < 1.3:
        return "NEO"
    if a < 2.0 and q >= 1.3:
        return "MarsCrosser/Hungaria-ish"
    if 2.0 <= a < 3.3 and q >= 1.3:
        return "MBA"
    if 3.3 <= a < 5.5 and q >= 1.3:
        return "Hilda/Cybele/Jupiter-region"
    if 5.5 <= a < 30.1:
        return "Centaur/JFC-region"
    if a >= 30.1:
        return "TNO"
    return "other"


def radec_to_unit(ra_deg, dec_deg):
    ra = np.radians(ra_deg); de = np.radians(dec_deg)
    return np.stack([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)], axis=-1)


def unit_to_radec(u):
    u = np.asarray(u, float)
    ra = np.degrees(np.arctan2(u[..., 1], u[..., 0])) % 360.0
    de = np.degrees(np.arcsin(np.clip(u[..., 2] / np.linalg.norm(u, axis=-1), -1, 1)))
    return ra, de
