"""Orbit-based vetting of tracklets (precision stage of the moving-object linker).

A tracklet is a short arc: positions in >= 3 exposures give the sky-plane position and angular velocity but neither distance nor radial velocity.  What an
orbit fit *can* say about it is whether its motion is compatible with Sun-bound bodies at all, and how typical it is of the asteroid population seen
in that direction on that date.  This module computes the second from a population synthesis:

  population_rates()   draws bound two-body orbits from a small-body element table (JPL SBDB columns a, e, i, node, peri, H), puts each at a random mean
                       anomaly, computes the geocentric apparent position / angular velocity / V magnitude at the observation epoch with the real Earth state
                       (tracklet.earth_helio_state) and keeps the objects that appear within `radius_deg` of the field and brighter than `vmax`.
  RatePrior            Gaussian-kernel density of the apparent rate vector (east, north; arcsec/h) of that sample, log density, and the density level that
                       contains a given fraction of the population (highest-density region).
  vet_tracklets()      adds `log_prior`, `in_hdr` and `bound_feasible` to tracklet dicts and optionally drops the ones outside the region.
Pure numpy.  The prior is a model: it assumes the synthesis (osculating elements treated as constant, light time and perturbations neglected,
random mean anomaly) describes the detectable population; NEOs, comets and TNOs outside the element table range are not represented.
"""
import math

import numpy as np

GM_SUN = 2.9591220828559e-4          # AU^3 / day^2
OBL = math.radians(23.4392911)


def _kepler(M, e):
    E = M + e * np.sin(M)
    for _ in range(12):
        E = E - (E - e * np.sin(E) - M) / (1.0 - e * np.cos(E))
    return E


def _state(a, e, inc, om, w, M):
    """Heliocentric ecliptic position (AU) and velocity (AU/d) of a bound orbit; angles in radians."""
    E = _kepler(M, e)
    n = np.sqrt(GM_SUN / a ** 3)
    cE, sE = np.cos(E), np.sin(E)
    x = a * (cE - e); y = a * np.sqrt(1 - e * e) * sE
    f = n * a / (1 - e * cE)
    vx = -f * sE; vy = f * np.sqrt(1 - e * e) * cE
    cO, sO, ci, si, cw, sw = np.cos(om), np.sin(om), np.cos(inc), np.sin(inc), np.cos(w), np.sin(w)
    R = np.array([[cO * cw - sO * sw * ci, -cO * sw - sO * cw * ci],
                  [sO * cw + cO * sw * ci, -sO * sw + cO * cw * ci],
                  [sw * si, cw * si]])
    pos = np.stack([R[0, 0] * x + R[0, 1] * y, R[1, 0] * x + R[1, 1] * y, R[2, 0] * x + R[2, 1] * y], 1)
    vel = np.stack([R[0, 0] * vx + R[0, 1] * vy, R[1, 0] * vx + R[1, 1] * vy, R[2, 0] * vx + R[2, 1] * vy], 1)
    return pos, vel


def _ecl2eq(v):
    c, s = math.cos(OBL), math.sin(OBL)
    return np.stack([v[:, 0], c * v[:, 1] - s * v[:, 2], s * v[:, 1] + c * v[:, 2]], 1)


def hg_magnitude(H, r, delta, rE, G=0.15):
    """Apparent V of an asteroid (HG system) from heliocentric r, geocentric delta and Earth-Sun distance rE (AU)."""
    cosb = np.clip((r * r + delta * delta - rE * rE) / (2 * r * delta), -1, 1)
    alpha = np.arccos(cosb)
    t = np.tan(alpha / 2)
    phi = (1 - G) * np.exp(-3.33 * t ** 0.63) + G * np.exp(-1.87 * t ** 1.22)
    return H + 5 * np.log10(r * delta) - 2.5 * np.log10(np.maximum(phi, 1e-6))


def population_rates(elements, mjd_utc, ra, dec, radius_deg=3.0, n_draw=20, vmax=22.0, seed=1, earth_state=None, max_orbits=400000):
    """elements: dict of arrays a (AU), e, i, om, w (deg), H.  -> dict(rate_e, rate_n [arcsec/h], V, delta, r, n) of the synthetic objects seen within
    radius_deg of (ra, dec) at the epoch.  earth_state = (pos AU, vel AU/d) equatorial; default tracklet.earth_helio_state."""
    if earth_state is None:
        from .tracklet import earth_helio_state
        earth_state = earth_helio_state(mjd_utc)
    pe, ve = earth_state
    rng = np.random.default_rng(seed)
    n_all = len(elements['a'])
    sel = rng.choice(n_all, min(n_all, max_orbits), replace=False) if n_all > max_orbits else np.arange(n_all)
    a = np.asarray(elements['a'], float)[sel]; e = np.asarray(elements['e'], float)[sel]; inc = np.radians(np.asarray(elements['i'], float)[sel])
    om = np.radians(np.asarray(elements['om'], float)[sel]); w = np.radians(np.asarray(elements['w'], float)[sel]); H = np.asarray(elements['H'], float)[sel]
    ok = (e < 0.95) & (a > 0.5) & np.isfinite(H)
    a, e, inc, om, w, H = a[ok], e[ok], inc[ok], om[ok], w[ok], H[ok]
    al, de = math.radians(ra), math.radians(dec)
    los = np.array([math.cos(de) * math.cos(al), math.cos(de) * math.sin(al), math.sin(de)])
    ex = np.array([-math.sin(al), math.cos(al), 0.0]); ey = np.array([-math.sin(de) * math.cos(al), -math.sin(de) * math.sin(al), math.cos(de)])
    cosr = math.cos(math.radians(radius_deg))
    rE = float(np.linalg.norm(pe))
    out = dict(rate_e=[], rate_n=[], V=[], delta=[], r=[], H=[])
    for _ in range(n_draw):
        M = rng.uniform(0, 2 * np.pi, len(a))
        pos, vel = _state(a, e, inc, om, w, M)
        pos = _ecl2eq(pos); vel = _ecl2eq(vel)
        g = pos - pe[None, :]; gv = vel - ve[None, :]
        d = np.linalg.norm(g, axis=1)
        u = g / d[:, None]
        m = (u @ los) > cosr
        if not m.any():
            continue
        g, gv, d, u = g[m], gv[m], d[m], u[m]
        r = np.linalg.norm(pos[m], axis=1)
        vperp = gv - (gv * u).sum(1)[:, None] * u
        arc = 3600.0 * 180.0 / math.pi / 24.0                       # rad/day -> arcsec/h
        # sky-plane components at the object's own direction, east/north unit vectors there
        ra_o = np.arctan2(u[:, 1], u[:, 0]); de_o = np.arcsin(np.clip(u[:, 2], -1, 1))
        exo = np.stack([-np.sin(ra_o), np.cos(ra_o), np.zeros_like(ra_o)], 1)
        eyo = np.stack([-np.sin(de_o) * np.cos(ra_o), -np.sin(de_o) * np.sin(ra_o), np.cos(de_o)], 1)
        out['rate_e'].append((vperp * exo).sum(1) / d * arc); out['rate_n'].append((vperp * eyo).sum(1) / d * arc)
        out['delta'].append(d); out['r'].append(r); out['H'].append(H[m])
        out['V'].append(hg_magnitude(H[m], r, d, rE))
    for k in out:
        out[k] = np.concatenate(out[k]) if out[k] else np.zeros(0)
    keep = out['V'] < vmax
    for k in out:
        out[k] = out[k][keep]
    out['n'] = int(keep.sum())
    return out


class RatePrior:
    """Gaussian-kernel density of the apparent rate vector (east, north; arcsec/h) of a synthetic population."""

    def __init__(self, rate_e, rate_n, bandwidth=3.0, max_points=20000, seed=1):
        pts = np.c_[rate_e, rate_n]
        if len(pts) > max_points:
            pts = pts[np.random.default_rng(seed).choice(len(pts), max_points, replace=False)]
        self.pts = pts
        self.h = float(bandwidth)
        self.n_pop = len(rate_e)
        self._self_logp = np.sort(self.logp(pts[: min(len(pts), 4000), 0], pts[: min(len(pts), 4000), 1]))

    def logp(self, re, rn):
        """log of the density (per (arcsec/h)^2) of the population at the rate vectors."""
        re = np.atleast_1d(np.asarray(re, float)); rn = np.atleast_1d(np.asarray(rn, float))
        out = np.empty(len(re))
        for s in range(0, len(re), 500):
            d2 = (re[s:s + 500, None] - self.pts[None, :, 0]) ** 2 + (rn[s:s + 500, None] - self.pts[None, :, 1]) ** 2
            dens = np.exp(-0.5 * d2 / self.h ** 2).sum(1) / (len(self.pts) * 2 * math.pi * self.h ** 2)
            out[s:s + 500] = np.log(np.maximum(dens, 1e-300))
        return out

    def hdr_level(self, fraction=0.99):
        """log-density level above which `fraction` of the population lies (highest-density region)."""
        return float(np.quantile(self._self_logp, 1.0 - fraction))


def vet_tracklets(trs, prior, fraction=0.99, drop=False):
    """Adds log_prior / in_hdr (fraction of the population) to tracklet dicts (rate vector from rate_ash and pa_deg: east = rate sin(pa), north = rate cos(pa))."""
    lvl = prior.hdr_level(fraction)
    re = np.array([t['rate_ash'] * math.sin(math.radians(t['pa_deg'])) for t in trs])
    rn = np.array([t['rate_ash'] * math.cos(math.radians(t['pa_deg'])) for t in trs])
    lp = prior.logp(re, rn) if len(trs) else np.zeros(0)
    for t, v in zip(trs, lp):
        t['log_prior'] = float(v); t['in_hdr'] = bool(v >= lvl)
    return [t for t in trs if t['in_hdr']] if drop else trs


# ---------------------------------------------------------------------------------------------------------------- element tables / glue
def fetch_elements(path, hmax=19.5, timeout=300):
    """Download the numbered + unnumbered asteroid elements (JPL SBDB query API: a, e, i, om, w, H with a in 1.5-5.5 AU, H <= hmax; ~145 MB JSON for
    H <= 19.5) once and store them as .npz.  Returns the dict of arrays.  Only the distribution of orbits is used, not the epochs."""
    import json
    import os
    import requests
    if os.path.isfile(path):
        return dict(np.load(path))
    r = requests.get("https://ssd-api.jpl.nasa.gov/sbdb_query.api", timeout=timeout,
                     params={"fields": "a,e,i,om,w,H", "sb-cdata": json.dumps({"AND": ["H|RG|0|%g" % hmax, "a|RG|1.5|5.5"]}), "sb-kind": "a"})
    r.raise_for_status()
    rows = r.json()["data"]
    el = {k: np.array([float(x[j]) if x[j] not in (None, "") else np.nan for x in rows]) for j, k in enumerate(["a", "e", "i", "om", "w", "H"])}
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    np.savez(path, **el)
    return el


def prior_for_tracklets(trs, elements, mjd_utc=None, radius_deg=3.0, vmax=22.0, n_draw=40, bandwidth=3.0, seed=1):
    """RatePrior for the field of a tracklet list: field centre = mean position of the members, epoch = mean member time unless given."""
    ra = float(np.mean([np.mean(t['ra']) for t in trs])); dec = float(np.mean([np.mean(t['dec']) for t in trs]))
    mjd = float(np.mean([np.mean(t['t']) for t in trs])) if mjd_utc is None else mjd_utc
    pop = population_rates(elements, mjd, ra, dec, radius_deg=radius_deg, n_draw=n_draw, vmax=vmax, seed=seed)
    if pop['n'] < 200:
        return None
    return RatePrior(pop['rate_e'], pop['rate_n'], bandwidth=bandwidth)


DEFAULT_ELEMENTS = "~/.cache/ogfinder_regression/sbdb_elements_H19.5.npz"


def resolve_prior(spec, trs):
    """spec: RatePrior | elements dict | path of an .npz element table | 'auto' (cached SBDB table, downloaded once) -> RatePrior or None"""
    import os
    if isinstance(spec, RatePrior):
        return spec
    if isinstance(spec, str):
        path = os.path.expanduser(DEFAULT_ELEMENTS if spec == "auto" else spec)
        spec = fetch_elements(path) if spec == "auto" else dict(np.load(path))
    return prior_for_tracklets(trs, spec) if trs else None
