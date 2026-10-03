"""Cross-night tracklet linking by orbit fitting (2-body, heliocentric).

Input: tracklets (a few detections within one night/visit) from several nights.  Output: groups of tracklets that are
consistent with ONE Sun-bound orbit, together with the orbit and the fit statistics.

Method (a "linear + gravity" linker in the spirit of Gauss/Herget IOD, HelioLinC and Holman et al. 2018):
  1. each tracklet -> attributable (ra, dec, mu_a, mu_d at its mean time) with covariance from the linear tangent-plane fit;
  2. pair gate: for tracklet A a (rho, rho_dot) grid of admissible (energy < 0) heliocentric states is built; each is propagated
     with universal-variable 2-body motion (+ light time) to the epoch of a later-night tracklet B and compared with B's
     position and rate, with the uncertainties of both attributables propagated over the gap (the unknown range makes the
     comparison a one-parameter-family search, not a point test);
  3. the best grid hypotheses of each surviving pair seed a full 6-parameter least-squares fit (Levenberg-Marquardt, state
     at A's epoch) of all observations of A and B; the pair is accepted when chi2/dof and the largest residual pass;
  4. accepted pairs are grown greedily into chains: another tracklet is added when the current orbit predicts it within
     tolerance and the refit of all observations stays acceptable; at most one tracklet per night and station per chain.
Model errors (planetary perturbations, topocentric parallax approximations, light-time) are covered by an additive error
floor sigma_eff^2 = sigma^2 + floor0^2 + (floor_rate * |t - t_ref|)^2 (arcsec, days).

Limits: pure 2-body Sun-only dynamics (arcs of up to a few weeks are fine for main-belt objects with a floor of ~1 arcsec;
NEO close approaches and arcs of months need an N-body propagator -- see moving.orbit with ASSIST); observer = geocentre
(+ station offset from the MPC parallax constants when `code` is given); no covariance-based Mahalanobis test of the final fit
(chi2/dof and max residual only); greedy (not globally optimal) chain assignment.
"""
import numpy as np
from scipy.optimize import least_squares
from . import kepler as K
from .util import C_AUD, AU_KM

ARC = np.pi / 180.0 / 3600.0            # radian per arcsec
MU = K.MU


# ------------------------------------------------------------------------------------------------ ephemeris helpers
def earth_helio(mjd_utc):
    """Heliocentric Earth-Moon barycentre (pos AU, vel AU/d), ICRF equatorial, astropy built-in ephemeris; arrays (n,3)."""
    from astropy.time import Time
    from astropy.coordinates import get_body_barycentric_posvel, solar_system_ephemeris
    t = Time(np.atleast_1d(np.asarray(mjd_utc, float)), format="mjd", scale="utc")
    with solar_system_ephemeris.set("builtin"):
        pe, ve = get_body_barycentric_posvel("earth", t); ps, vs = get_body_barycentric_posvel("sun", t)
    return (pe - ps).xyz.to("au").value.T, (ve - vs).xyz.to("au/day").value.T


def observer_helio(mjd_utc, code=None):
    """Heliocentric observer position (n,3) AU; `code` None/'500' = geocentre, else MPC station (parallax constants)."""
    mjd = np.atleast_1d(np.asarray(mjd_utc, float))
    pos, vel = earth_helio(mjd)
    if code not in (None, "500", ""):
        from . import obs as O
        pos = pos + O.geocentric_gcrs_km(code, mjd) / AU_KM
    return pos, vel


def _los(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.stack([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)], -1)


def _basis(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    ea = np.array([-np.sin(a), np.cos(a), 0.0])
    ed = np.array([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)])
    return ea, ed


def sky_diff(ra1, de1, ra2, de2):
    """(dx east, dy north) in arcsec of point 1 relative to point 2 (small-angle, cos(dec) applied)."""
    dx = ((np.asarray(ra1) - ra2 + 180.0) % 360.0 - 180.0) * np.cos(np.radians(de2)) * 3600.0
    dy = (np.asarray(de1) - de2) * 3600.0
    return dx, dy


def predict(r0, v0, t0, t, obs_pos, lt_iter=3):
    """2-body prediction of apparent (ra, dec) [deg] and range for state(s) r0,v0 (...,3) at epoch t0 [MJD] observed at times t
    [broadcastable] from heliocentric observer positions obs_pos (...,3); light-time corrected."""
    r0 = np.asarray(r0, float); v0 = np.asarray(v0, float)
    tau = np.zeros(np.broadcast(np.asarray(t, float), r0[..., 0]).shape)
    for _ in range(lt_iter):
        r, _ = K.propagate_2body(r0, v0, np.asarray(t, float) - tau - t0)
        d = r - obs_pos
        rng = np.linalg.norm(d, axis=-1)
        tau = rng / C_AUD
    ra = np.degrees(np.arctan2(d[..., 1], d[..., 0])) % 360.0
    de = np.degrees(np.arcsin(np.clip(d[..., 2] / rng, -1, 1)))
    return ra, de, rng


# ------------------------------------------------------------------------------------------------ tracklets
class Tracklet:
    """One night's detections of one candidate moving object.  mjd UTC; ra/dec deg; sig arcsec (scalar or array)."""

    def __init__(self, mjd, ra, dec, sig=0.3, code=None, tid=None, night=None, truth=None):
        o = np.argsort(mjd)
        self.mjd = np.asarray(mjd, float)[o]; self.ra = np.asarray(ra, float)[o]; self.dec = np.asarray(dec, float)[o]
        self.sig = np.broadcast_to(np.asarray(sig, float), self.mjd.shape).copy()[o] if np.ndim(sig) else np.full(len(o), float(sig))
        self.code = code; self.id = tid; self.truth = truth
        if night is None:
            lon = 0.0
            if code not in (None, "500", ""):
                try:
                    from . import obs as O
                    lon = O.obscodes().get(code, (0.0,))[0] or 0.0
                except Exception:
                    lon = 0.0
            night = int(np.floor(np.mean(self.mjd) + lon / 360.0 + 0.5))      # local noon-to-noon "night" index of the station
        self.night = night
        self._fit()
        self._obs = None

    def __len__(self):
        return len(self.mjd)

    def _fit(self):
        t = self.mjd; self.t0 = float(np.average(t, weights=1.0 / self.sig ** 2))
        self.ra0 = float(self.ra[len(t) // 2]); self.dec0 = float(self.dec[len(t) // 2])
        x, y = sky_diff(self.ra, self.dec, self.ra0, self.dec0)
        dt = t - self.t0; w = 1.0 / self.sig ** 2
        sdt2 = float(np.sum(w * dt ** 2)); sw = float(np.sum(w))
        if len(t) >= 2 and sdt2 > 0:
            mx = np.sum(w * dt * x) / sdt2; my = np.sum(w * dt * y) / sdt2
            cx = np.sum(w * x) / sw; cy = np.sum(w * y) / sw
            self.mu_x, self.mu_y = float(mx), float(my)           # arcsec/day
            self.sig_mu = float(1.0 / np.sqrt(sdt2))              # arcsec/day
        else:
            cx = cy = 0.0; self.mu_x = self.mu_y = 0.0; self.sig_mu = 1e4
        self.sig_pos = float(1.0 / np.sqrt(sw))
        cd = np.cos(np.radians(self.dec0))
        self.ra_c = (self.ra0 + cx / 3600.0 / cd) % 360.0; self.dec_c = self.dec0 + cy / 3600.0
        self.arc_days = float(t.max() - t.min())

    def obs_t0(self):
        """(position, velocity) of the observer at the tracklet's reference time (cached)."""
        if getattr(self, "_obs0", None) is None:
            p, v = observer_helio(self.t0, self.code)
            self._obs0 = (p[0], v[0])
        return self._obs0

    def obs(self):
        if self._obs is None:
            self._obs = observer_helio(self.mjd, self.code)
        return self._obs


def tracklets_from_json(trs, night_of=None, sig_default=0.3):
    """Convert moving.pipeline / ds9_moving `tracklets.json` entries (keys t/mjd, ra, dec, sig) into Tracklet objects."""
    out = []
    for k, t in enumerate(trs):
        mjd = t.get("mjd", t.get("t"))
        if mjd is None:
            continue
        out.append(Tracklet(mjd, t["ra"], t["dec"], t.get("sig", sig_default), code=t.get("code"), tid=t.get("id", k),
                            night=None if night_of is None else night_of(t)))
    return out


# ------------------------------------------------------------------------------------------------ hypotheses
def hypotheses(trk, n_rho=20, n_rd=13, rho_range=(0.02, 8.0), amax=100.0):
    """Admissible heliocentric states (energy < 0, a < amax) for tracklet attributable on a (log rho, rho_dot/v_esc) grid.
    Returns r (H,3), v (H,3), rho (H,), rhod (H,)."""
    op, ov = trk.obs_t0()
    los = _los(trk.ra_c, trk.dec_c); ea, ed = _basis(trk.ra_c, trk.dec_c)
    cd = np.cos(np.radians(trk.dec_c))
    mu_a = trk.mu_x * ARC; mu_d = trk.mu_y * ARC                  # rad/day (east incl. cos dec, north)
    rho = np.geomspace(rho_range[0], rho_range[1], n_rho)
    R, V, RHO, RD = [], [], [], []
    for rh in rho:
        r = op + rh * los
        rn = np.linalg.norm(r)
        vesc = np.sqrt(2 * MU / rn)
        for x in np.linspace(-1.0, 1.0, n_rd):
            rd = x * vesc
            v = ov + rd * los + rh * (mu_a * ea + mu_d * ed)
            E = 0.5 * v @ v - MU / rn
            if E >= 0 or rn < 0.15:
                continue
            a = -MU / (2 * E)
            if a > amax:
                continue
            R.append(r); V.append(v); RHO.append(rh); RD.append(rd)
    return np.array(R), np.array(V), np.array(RHO), np.array(RD)


# ------------------------------------------------------------------------------------------------ pair gate
def pair_gate(A, B, hyp=None, model_pos=3.0, model_rate_floor=1.0, top=4, **hkw):
    """Gate for a later tracklet B given earlier A.  Returns (chi2_min, [indices of the best hypotheses], hyp)."""
    if hyp is None:
        hyp = hypotheses(A, **hkw)
    R, V = hyp[0], hyp[1]
    if len(R) == 0:
        return np.inf, [], hyp
    ob, ovb = B.obs_t0()
    h = 0.02                                                     # position and rate from one stacked propagation (t_B and t_B + h)
    tt = np.array([B.t0, B.t0 + h])[:, None]
    oo = np.stack([ob, ob + ovb * h])[:, None, :]
    ra_, de_, _ = predict(R[None], V[None], A.t0, tt, oo)
    ra, de, ra2, de2 = ra_[0], de_[0], ra_[1], de_[1]
    px, py = sky_diff(ra, de, B.ra_c, B.dec_c)                     # predicted minus measured
    qx, qy = sky_diff(ra2, de2, ra, de); qx /= h; qy /= h         # predicted rate, arcsec/day
    dtg = abs(B.t0 - A.t0)
    sp2 = A.sig_pos ** 2 + B.sig_pos ** 2 + (dtg * A.sig_mu) ** 2 + (dtg * B.sig_mu) ** 2 + model_pos ** 2
    sm2 = A.sig_mu ** 2 + B.sig_mu ** 2 + model_rate_floor ** 2 + (0.02 * np.hypot(qx, qy)) ** 2
    chi = (px ** 2 + py ** 2) / sp2 + ((qx - B.mu_x) ** 2 + (qy - B.mu_y) ** 2) / sm2
    idx = np.argsort(chi)[:top]
    return float(chi[idx[0]]), idx.tolist(), hyp


# ------------------------------------------------------------------------------------------------ fitting
def _eff_sigma(trk, t_ref, floor0, floor_rate):
    return np.sqrt(trk.sig ** 2 + floor0 ** 2 + (floor_rate * np.abs(trk.mjd - t_ref)) ** 2)


class Fit:
    pass


def _stack(trks):
    mjd = np.concatenate([t.mjd for t in trks]); ra = np.concatenate([t.ra for t in trks]); de = np.concatenate([t.dec for t in trks])
    sig = np.concatenate([t.sig for t in trks]); op = np.concatenate([t.obs()[0] for t in trks])
    return mjd, ra, de, sig, op


def fit_orbit(trks, state0, t_ref, floor0=0.3, floor_rate=0.05, max_nfev=120):
    """Full 2-body least squares on all observations of `trks` starting from heliocentric state0=(r,v) at MJD t_ref.
    Parameters are the state vector (AU, AU/d); residuals are (dx, dy)/sigma_eff.  Returns Fit or None."""
    mjd, ra, de, sig, op = _stack(trks)
    seff = np.sqrt(sig ** 2 + floor0 ** 2 + (floor_rate * np.abs(mjd - t_ref)) ** 2)
    n = len(mjd)

    def res(p):
        try:
            pra, pde, _ = predict(p[:3], p[3:], t_ref, mjd, op)
        except Exception:
            return np.full(2 * n + 1, 1e3)
        if not np.all(np.isfinite(pra)):
            return np.full(2 * n + 1, 1e3)
        dx, dy = sky_diff(ra, de, pra, pde)
        # soft prior: Sun-bound orbit (v/v_esc < 0.98); keeps short-arc fits from wandering onto hyperbolic solutions
        ratio = np.linalg.norm(p[3:]) / np.sqrt(2.0 * MU / np.linalg.norm(p[:3]))
        return np.r_[dx / seff, dy / seff, max(0.0, ratio - 0.98) / 0.02]

    p0 = np.r_[np.asarray(state0[0], float), np.asarray(state0[1], float)]
    xs = np.r_[np.full(3, 0.01), np.full(3, 1e-4)]
    try:
        sol = least_squares(res, p0, x_scale=xs, method="lm", max_nfev=max_nfev * 7)
    except Exception:
        return None
    f = Fit(); f.state = (sol.x[:3].copy(), sol.x[3:].copy()); f.t_ref = t_ref
    r = sol.fun[:-1]; f.chi2 = float(np.sum(r ** 2)); f.dof = max(2 * n - 6, 1)
    f.chi2_red = f.chi2 / f.dof; f.max_resid = float(np.max(np.abs(r))); f.n_obs = n
    dx, dy = r[:n] * seff, r[n:] * seff
    f.rms_arcsec = float(np.sqrt(np.mean(dx ** 2 + dy ** 2) / 2.0))
    f.seff = seff
    return f


def orbit_summary(f):
    r, v = f.state
    try:
        el = K.state_to_elements(r, v, ecliptic_input=False)
    except Exception:
        return {}
    return dict(a=float(el["a"]), e=float(el["e"]), inc_deg=float(el["i"]), q=float(el["q"]), r_au=float(np.linalg.norm(r)))


def predict_tracklet_resid(f, C, floor0=0.3, floor_rate=0.05):
    """rms normalised residual of tracklet C against fit f's orbit (used before refitting with C)."""
    pra, pde, _ = predict(f.state[0], f.state[1], f.t_ref, C.mjd, C.obs()[0])
    dx, dy = sky_diff(C.ra, C.dec, pra, pde)
    se = np.sqrt(C.sig ** 2 + floor0 ** 2 + (floor_rate * np.abs(C.mjd - f.t_ref)) ** 2)
    return float(np.sqrt(np.mean((dx / se) ** 2 + (dy / se) ** 2)))


def acceptable(f, chi2_max, resid_max):
    if f is None:
        return False
    r, v = f.state
    rn = np.linalg.norm(r)
    if not np.isfinite(f.chi2) or rn < 0.1 or 0.5 * v @ v - MU / rn >= 0:
        return False
    return f.chi2_red <= chi2_max and f.max_resid <= resid_max


def fit_pair(A, B, hyp=None, idx=None, floor0=0.3, floor_rate=0.05, chi2_max=4.0, resid_max=6.0, top=4):
    """Full fit of A+B from the best gate hypotheses.  Returns the best Fit (not necessarily acceptable) or None."""
    if hyp is None or idx is None:
        _, idx, hyp = pair_gate(A, B, top=top)
    best = None
    for i in idx:
        f = fit_orbit([A, B], (hyp[0][i], hyp[1][i]), A.t0, floor0, floor_rate)
        if f is not None and (best is None or f.chi2 < best.chi2):
            best = f
        if best is not None and acceptable(best, chi2_max * 0.5, resid_max):
            break
    return best


# ------------------------------------------------------------------------------------------------ linker
def link_nights(*a, **k):
    with np.errstate(all="ignore"):
        return _link_nights(*a, **k)


def _link_nights(trks, max_gap_days=30.0, min_gap_days=0.3, gate_chi2=60.0, chi2_max=4.0, resid_max=6.0, floor0=0.3, floor_rate=0.05,
                max_motion_deg_day=2.5, n_rho=20, n_rd=13, top=4, one_per_night=True, grow_tol=6.0, min_tracklets=3, progress=None):
    """Link tracklets from different nights into orbit-consistent groups.
    Groups with fewer than `min_tracklets` (default 3) tracklets are returned separately as `pairs`: two short-arc tracklets give only 2
    degrees of freedom of orbit constraint (8 attributable numbers vs 6 orbital parameters), so in a dense field chance pairings pass the fit
    (measured: see docs/moving_objects.md); >= 3 tracklets leave 6 constraint degrees of freedom.
    Returns dict(groups=[dict(ids, nights, fit, elements, chi2_red, rms_arcsec, max_resid, n_obs, pair_scores)], unlinked=[ids],
    stats=dict(n_pairs_tested, n_gate, n_fit, n_accept))."""
    trks = sorted(trks, key=lambda t: t.t0)
    n = len(trks)
    cands = []; partners = {}
    stats = dict(n_tracklets=n, n_pairs_tested=0, n_gate=0, n_fit=0, n_accept=0)
    hyp_cache = {}
    for i in range(n):
        A = trks[i]
        for j in range(i + 1, n):
            B = trks[j]
            dt = B.t0 - A.t0
            if dt > max_gap_days:
                break
            if dt < min_gap_days or (one_per_night and A.night == B.night and A.code == B.code):
                continue
            sep = np.hypot(*sky_diff(B.ra_c, B.dec_c, A.ra_c, A.dec_c))
            if sep / 3600.0 > max_motion_deg_day * dt + 0.2:
                continue
            stats["n_pairs_tested"] += 1
            if i not in hyp_cache:
                hyp_cache[i] = hypotheses(A, n_rho=n_rho, n_rd=n_rd)
            chi, idx, _ = pair_gate(A, B, hyp=hyp_cache[i], top=top)
            if chi > gate_chi2:
                continue
            stats["n_gate"] += 1
            f = fit_pair(A, B, hyp_cache[i], idx, floor0, floor_rate, chi2_max, resid_max, top)
            stats["n_fit"] += 1
            if acceptable(f, chi2_max, resid_max):
                stats["n_accept"] += 1
                cands.append((f.chi2_red, i, j, f))
                partners.setdefault(i, {})[j] = f.chi2_red; partners.setdefault(j, {})[i] = f.chi2_red
        if progress:
            progress(i, n)
    cands.sort(key=lambda c: c[0])
    used = set(); groups = []
    for _, i, j, f in cands:
        if i in used or j in used:
            continue
        members = [i, j]; cur = f; scores = [f.chi2_red]; rejected = set()
        while True:
            nights = {(trks[k].night, trks[k].code) for k in members}
            best = None
            for k in range(n):
                if k in used or k in members or k in rejected:
                    continue
                if one_per_night and (trks[k].night, trks[k].code) in nights:
                    continue
                if abs(trks[k].t0 - cur.t_ref) > max_gap_days:
                    continue
                rr = predict_tracklet_resid(cur, trks[k], floor0, floor_rate)
                pc = [partners[k][m] for m in members if m in partners.get(k, {})]
                if pc:                                  # accepted pair-link with a member: always a candidate, ranked first
                    rr = min(rr, grow_tol) * 1e-3 + float(np.mean(pc)) * 1e-4
                if rr <= grow_tol and (best is None or rr < best[0]):
                    best = (rr, k)
            if best is None:
                break
            k = best[1]
            f2 = fit_orbit([trks[m] for m in members + [k]], cur.state, cur.t_ref, floor0, floor_rate)
            if acceptable(f2, chi2_max, resid_max):
                members.append(k); cur = f2; scores.append(f2.chi2_red)
            else:
                rejected.add(k)
        used.update(members)
        groups.append(dict(ids=[trks[m].id for m in members], nights=[trks[m].night for m in members],
                           elements=orbit_summary(cur), chi2_red=cur.chi2_red, rms_arcsec=cur.rms_arcsec, max_resid=cur.max_resid,
                           n_obs=cur.n_obs, fit=cur, truth=[trks[m].truth for m in members]))
    pairs = [g for g in groups if len(g["ids"]) < min_tracklets]
    groups = [g for g in groups if len(g["ids"]) >= min_tracklets]
    linked = {i for g in groups for i in g["ids"]}
    return dict(groups=groups, pairs=pairs, unlinked=[t.id for t in trks if t.id not in linked], stats=stats)


# ------------------------------------------------------------------------------------------------ vetting
def vet_with_links(trks, result):
    """Per-tracklet status: 'linked' (member of an orbit-consistent group of >= 3 tracklets), 'pair' (orbit-consistent with exactly one other tracklet: a candidate, false-pair rate
    is high in dense fields) else 'single'.  Single-night tracklets
    should still go through moving.orbitlink.vet_tracklets (population rate prior); linked ones have been verified with an
    orbit fit over >= 2 nights, which is a much stronger test than the rate prior."""
    st = {}
    for gi, g in enumerate(result.get("pairs", [])):
        for tid in g["ids"]:
            st[tid] = dict(status="pair", group=gi, chi2_red=g["chi2_red"], a=g["elements"].get("a"), e=g["elements"].get("e"))
    for gi, g in enumerate(result["groups"]):
        for tid in g["ids"]:
            st[tid] = dict(status="linked", group=gi, chi2_red=g["chi2_red"], a=g["elements"].get("a"), e=g["elements"].get("e"))
    return {t.id: st.get(t.id, dict(status="single")) for t in trks}
