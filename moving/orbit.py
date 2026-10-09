"""Orbit propagation and batch least-squares differential correction.

Force model: the in-tree test-particle integrator (ogfmeas.ssbody). Planet and
Moon positions come from astropy's built-in ephemeris. Forces are Newtonian
point masses, the solar Schwarzschild term, solar J2, and Earth J2/J3/J4.
The 16 massive asteroids and a DE440 kernel are not used. The result records
the model name in ``force_model``.

Observation model: barycentric ICRF positions; light-time iteration (first-order Taylor
correction per perturbed particle); astrometric (catalog-reduced) directions, so NO stellar
aberration correction is applied (it is already absorbed by the star-catalog reduction); solar
gravitational light deflection is neglected (<~ 10 mas for elongations > 90 deg).
"""
import os
import pathlib
import sys
import numpy as np
from .util import EPHEM_DIR, C_AUD, ARCSEC, log, utc_mjd_to_tdb_jd
from . import kepler as K
from . import obs as OBS

PLANETS_FILE = "linux_p1550p2650.440"
ASTEROID_FILE = "sb441-n16.bsp"


def _ssbody():
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas.ssbody as ssbody
    return ssbody


def ephem_paths():
    return (os.path.join(EPHEM_DIR, PLANETS_FILE), os.path.join(EPHEM_DIR, ASTEROID_FILE))


def assist_available():
    """True when the in-tree solar-system integrator can be constructed.

    The name is kept for the desktop setup report. This does not import
    rebound or assist, and it does not require an ephemeris kernel.
    """
    try:
        _ssbody()
        return True
    except Exception:
        return False


class Propagator:
    """Integrate test particles. Times are JD TDB."""
    def __init__(self, ephem=None):
        if not assist_available():
            self.force_model = "UNAVAILABLE"
            raise RuntimeError("in-tree solar-system integrator is unavailable")
        self._planets = ephem or _ssbody().BuiltinEphemeris()
        self.force_model = _ssbody().FORCE_MODEL

    def body(self, name, jd):
        """Barycentric position and velocity of a planet or the Moon at JD TDB."""
        return self._planets.state(name, jd)

    def sun(self, jd):
        return self.body("sun", jd)

    def propagate(self, states, jd0, jds):
        """states: (n,6) barycentric at jd0; jds: (m,) target times (any order).
        Returns (n, m, 6) states."""
        return self._planets.propagate(states, jd0, jds)


def propagate_one(prop, state, jd0, jds):
    """Simple helper: states at each jd (one particle) using a fresh simulation per sweep."""
    return prop.propagate(np.array([state]), jd0, jds)[0]


def propagate_assist(states, jd0, jds):
    """Integrate with the user-installed rebound and assist packages.

    Those packages run in a child process, and only when this function is
    called. A missing package or a missing kernel file raises RuntimeError
    and does not return in-tree states. Propagator stays on ogfmeas.ssbody.
    """
    import json
    import subprocess
    root = str(pathlib.Path(__file__).resolve().parents[1])
    payload = json.dumps({
        "states": np.atleast_2d(np.asarray(states, float)).tolist(),
        "jd0": float(jd0),
        "jds": np.asarray(jds, float).tolist(),
    })
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    env.setdefault("MKL_NUM_THREADS", "2")
    env.setdefault("NUMEXPR_NUM_THREADS", "2")
    proc = subprocess.run(
        [sys.executable, "-m", "moving.assist_worker"],
        input=payload,
        capture_output=True,
        text=True,
        cwd=root,
        env=env,
    )
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        if proc.returncode == 2:
            raise RuntimeError(
                "assist was not found. The in-tree orbit model is ogfmeas. "
                + detail)
        raise RuntimeError("assist worker failed: " + detail)
    out = json.loads(proc.stdout)
    return np.asarray(out["states"], float), out["force_model"]


# ------------------------------------------------------------------ observations container
class Obs:
    """Arrays describing astrometric observations. All angles in degrees, sigmas in arcsec."""
    def __init__(self, mjd_utc, ra, dec, sig_ra, sig_dec, stn, sat_xyz_km=None, mag=None, band=None,
                 cat=None, ids=None, obs_pos_bary_au=None):
        self.mjd_utc = np.asarray(mjd_utc, float)
        self.ra = np.asarray(ra, float); self.dec = np.asarray(dec, float)
        n = len(self.ra)
        self.sig_ra = np.broadcast_to(np.asarray(sig_ra, float), (n,)).copy()   # on-sky (RA*cosDec)
        self.sig_dec = np.broadcast_to(np.asarray(sig_dec, float), (n,)).copy()
        self.stn = list(stn) if not isinstance(stn, str) else [stn] * n
        self.sat = sat_xyz_km if sat_xyz_km is not None else [None] * n
        self.mag = np.full(n, np.nan) if mag is None else np.asarray(mag, float)
        self.band = band if band is not None else [""] * n
        self.cat = cat if cat is not None else [""] * n
        self.ids = ids if ids is not None else list(range(n))
        self.jd_tdb = utc_mjd_to_tdb_jd(self.mjd_utc)
        self.n = n
        self._obspos = obs_pos_bary_au
        self.active = np.ones(n, bool)

    def subset(self, mask):
        mask = np.asarray(mask)
        idx = np.where(mask)[0] if mask.dtype == bool else mask
        o = Obs(self.mjd_utc[idx], self.ra[idx], self.dec[idx], self.sig_ra[idx], self.sig_dec[idx],
                [self.stn[i] for i in idx], [self.sat[i] for i in idx], self.mag[idx],
                [self.band[i] for i in idx], [self.cat[i] for i in idx], [self.ids[i] for i in idx])
        if self._obspos is not None:
            o._obspos = self._obspos[idx]
        return o

    def observer_bary(self, prop):
        """Barycentric observer positions & velocities (n,6), AU and AU/d (velocity: Earth only)."""
        if self._obspos is not None and self._obspos.shape[1] == 6:
            return self._obspos
        out = np.zeros((self.n, 6))
        stns = sorted(set(self.stn))
        for stn in stns:
            idx = [i for i in range(self.n) if self.stn[i] == stn]
            sat = self.sat[idx[0]] if stn not in OBS.obscodes() or OBS.obscodes()[stn][0] is None else None
            if sat is not None:
                gc = np.array([self.sat[i] for i in idx])
            else:
                gc = OBS.geocentric_gcrs_km(stn, self.mjd_utc[idx])
            for k, i in enumerate(idx):
                e = prop.body("earth", self.jd_tdb[i])
                out[i, :3] = e[:3] + gc[k] / OBS.AU_KM
                out[i, 3:] = e[3:]
        self._obspos = out
        return out


# ------------------------------------------------------------------ prediction
def _radec_from_vec(d):
    ra = np.degrees(np.arctan2(d[..., 1], d[..., 0])) % 360.0
    dec = np.degrees(np.arcsin(d[..., 2] / np.linalg.norm(d, axis=-1)))
    return ra, dec


def predict_radec(prop, states, jd0, obs, obspos=None, tau0=None):
    """Predict (ra, dec) [deg] for each particle (n,6) at all observations.
    Returns ra (n,m), dec (n,m), tau_nom (m,), rho (n,m)."""
    states = np.atleast_2d(states)
    n, m = len(states), obs.n
    op = obs.observer_bary(prop) if obspos is None else obspos
    jdobs = obs.jd_tdb
    if tau0 is None:
        # 2-body estimate for the first light-time guess
        r, v = K.propagate_2body(states[0, :3][None, :], states[0, 3:][None, :], jdobs - jd0)
        tau0 = np.linalg.norm(r - op[:, :3], axis=-1) / C_AUD
    # two passes on the nominal particle for tau (cheap relative to full set)
    tau = tau0
    for _ in range(2):
        S = prop.propagate(states[:1], jd0, jdobs - tau)[0]
        tau = np.linalg.norm(S[:, :3] - op[:, :3], axis=-1) / C_AUD
    te_nom = jdobs - tau
    SS = prop.propagate(states, jd0, te_nom)          # (n,m,6)
    ra = np.zeros((n, m)); dec = np.zeros((n, m)); rho = np.zeros((n, m))
    for i in range(n):
        x = SS[i, :, :3].copy(); v = SS[i, :, 3:]
        tp = tau.copy()
        for _ in range(3):
            xx = x + v * (-(tp - tau))[:, None]           # move emission time by (tau - tp)
            d = xx - op[:, :3]
            tp = np.linalg.norm(d, axis=-1) / C_AUD
        # first-order Taylor in time: x(te_nom + (tau - tp))
        xx = x + v * (tau - tp)[:, None]
        d = xx - op[:, :3]
        rho[i] = np.linalg.norm(d, axis=-1)
        ra[i], dec[i] = _radec_from_vec(d)
    return ra, dec, tau, rho


def residuals_arcsec(ra_obs, dec_obs, ra_mod, dec_mod):
    d_ra = ((ra_obs - ra_mod + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec_obs)) * 3600.0
    d_de = (dec_obs - dec_mod) * 3600.0
    return d_ra, d_de


def partials(prop, state, jd0, obs, steps=None):
    """Finite-difference partials of (ra*cosdec, dec) [arcsec] w.r.t. the 6 state elements.
    Returns model (ra,dec,tau,rho), J (2m,6)."""
    if steps is None:
        steps = np.array([1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8])
    P = np.vstack([state] + [state + s * np.eye(6)[k] for k in range(6) for s in (steps[k], -steps[k])])
    ra, dec, tau, rho = predict_radec(prop, P, jd0, obs)
    m = obs.n
    J = np.zeros((2 * m, 6))
    cosd = np.cos(np.radians(obs.dec))
    for k in range(6):
        pr, mr = ra[1 + 2 * k], ra[2 + 2 * k]
        dra = ((pr - mr + 180) % 360 - 180) * cosd * 3600.0 / (2 * steps[k])
        dde = (dec[1 + 2 * k] - dec[2 + 2 * k]) * 3600.0 / (2 * steps[k])
        J[0::2, k] = dra; J[1::2, k] = dde
    return (ra[0], dec[0], tau, rho[0]), J


# ------------------------------------------------------------------ differential correction
class FitResult(dict):
    pass


def carpino_chi2(res, J, W, C):
    """Per-observation Carpino et al. (2003) statistic chi_i^2 = r^T (Gamma_i - H C H^T)^-1 r (2x2 blocks)."""
    m = len(res) // 2
    chi2 = np.zeros(m)
    for i in range(m):
        H = J[2 * i:2 * i + 2]
        G = np.linalg.inv(W[2 * i:2 * i + 2, 2 * i:2 * i + 2])
        R = G - H @ C @ H.T
        r = res[2 * i:2 * i + 2]
        try:
            chi2[i] = r @ np.linalg.solve(R, r)
        except np.linalg.LinAlgError:
            chi2[i] = np.inf
        if chi2[i] < 0:
            chi2[i] = np.inf
    return chi2


def differential_correction(prop, state0, jd0, obs, max_iter=12, chi_rej2=9.0, chi_rec2=7.0,
                            reject=True, verbose=False, lm_lambda0=0.0, fit_mask=None):
    """Batch least squares (Gauss-Newton with Levenberg-Marquardt damping).

    Weights: 1/sigma^2 per component (on-sky RA*cosDec and Dec).  Outlier handling follows the
    Carpino et al. (2003) chi_i^2 statistic with rejection / recovery thresholds chi_rej2 / chi_rec2
    (defaults are our choice, configurable).  Returns FitResult with state, covariance (6x6, state
    elements), residuals, chi2, rms, rejected flags.
    """
    state = np.array(state0, float)
    m = obs.n
    active = obs.active.copy() if fit_mask is None else np.asarray(fit_mask, bool).copy()
    sig = np.empty(2 * m); sig[0::2] = obs.sig_ra; sig[1::2] = obs.sig_dec
    lam = lm_lambda0
    prev_chi2 = None
    for it in range(max_iter):
        (ra, dec, tau, rho), J = partials(prop, state, jd0, obs)
        dra, dde = residuals_arcsec(obs.ra, obs.dec, ra, dec)
        res = np.empty(2 * m); res[0::2] = dra; res[1::2] = dde
        act2 = np.repeat(active, 2)
        w = 1.0 / sig ** 2
        A = J[act2]; r = res[act2]; ww = w[act2]
        N = A.T @ (A * ww[:, None]); b = A.T @ (ww * r)
        # column scaling for conditioning
        s = np.sqrt(np.diag(N)); s[s == 0] = 1
        Ns = N / np.outer(s, s)
        chi2 = float(np.sum(ww * r ** 2))
        try:
            dx = np.linalg.solve(Ns + lam * np.diag(np.diag(Ns)), b / s) / s
        except np.linalg.LinAlgError:
            dx = np.linalg.lstsq(Ns, b / s, rcond=None)[0] / s
        # trial step with step halving
        trial = state + dx
        (ra2, dec2, _, _) = predict_radec(prop, trial[None, :], jd0, obs, tau0=tau)[:4]
        d2a, d2d = residuals_arcsec(obs.ra, obs.dec, ra2[0], dec2[0])
        r2 = np.empty(2 * m); r2[0::2] = d2a; r2[1::2] = d2d
        chi2_new = float(np.sum(w[act2] * r2[act2] ** 2))
        if chi2_new <= chi2 * (1 + 1e-9):
            state = trial
            lam = max(lam / 5, 0.0)
            ok = True
        else:
            lam = max(lam * 10, 1e-3)
            ok = False
        if verbose:
            log("DC it=%d chi2 %.3f -> %.3f active=%d lam=%.2g |dx/sig| max=%.2g" % (
                it, chi2, chi2_new, active.sum(), lam, np.max(np.abs(dx) * s)))
        # outlier management using the up-to-date linearisation
        if reject and ok:
            (ra, dec, tau, rho), J = partials(prop, state, jd0, obs)
            dra, dde = residuals_arcsec(obs.ra, obs.dec, ra, dec)
            res = np.empty(2 * m); res[0::2] = dra; res[1::2] = dde
            act2 = np.repeat(active, 2)
            Wm = np.diag(w)
            Nn = J[act2].T @ (J[act2] * w[act2][:, None])
            try:
                C = np.linalg.inv(Nn)
                c2 = carpino_chi2(res, J, np.diag(w), C)
                new_active = active.copy()
                for i in range(m):
                    if active[i] and c2[i] > chi_rej2:
                        new_active[i] = False
                    elif (not active[i]) and c2[i] < chi_rec2:
                        new_active[i] = True
                if new_active.sum() >= 4 and (new_active != active).any():
                    active = new_active
                    prev_chi2 = None
                    continue
            except np.linalg.LinAlgError:
                pass
        if ok and (np.max(np.abs(dx) * s) < 1e-3 or abs(chi2 - chi2_new) < 1e-6 * max(chi2, 1.0)) and it > 0:
            break
    (ra, dec, tau, rho), J = partials(prop, state, jd0, obs)
    dra, dde = residuals_arcsec(obs.ra, obs.dec, ra, dec)
    res = np.empty(2 * m); res[0::2] = dra; res[1::2] = dde
    act2 = np.repeat(active, 2)
    w = 1.0 / sig ** 2
    N = J[act2].T @ (J[act2] * w[act2][:, None])
    C = np.linalg.inv(N)
    chi2 = float(np.sum(w[act2] * res[act2] ** 2))
    ndof = int(act2.sum() - 6)
    return FitResult(state=state, jd0=jd0, cov=C, residual_ra=dra, residual_dec=dde, active=active,
                     chi2=chi2, ndof=ndof, chi2_red=chi2 / max(ndof, 1),
                     rms_arcsec=float(np.sqrt(np.mean(res[act2] ** 2))), n_used=int(active.sum()),
                     n_total=m, force_model=prop.force_model, J=J, tau=tau, rho=rho,
                     model_ra=ra, model_dec=dec)


# ------------------------------------------------------------------ elements & covariance
def helio_elements(prop, state_bary, jd):
    """Osculating heliocentric ecliptic-J2000 elements for a barycentric state at JD TDB."""
    sun = prop.sun(jd)
    rel = np.asarray(state_bary) - sun
    el = K.state_to_elements(rel[:3], rel[3:])
    return el


ELEM_KEYS = ("a", "e", "i", "om", "w", "ma", "q")


def elements_and_sigmas(prop, state, cov, jd, jd_target=None):
    """Elements at jd_target (default jd) with 1-sigma via first-order propagation of the state covariance
    (numerical Jacobian of state->elements, and of the flow if jd_target != jd)."""
    jt = jd if jd_target is None else jd_target
    steps = np.array([1e-6, 1e-6, 1e-6, 1e-8, 1e-8, 1e-8])
    P = np.vstack([state] + [state + s * np.eye(6)[k] for k in range(6) for s in (steps[k], -steps[k])])
    if jt != jd:
        S = prop.propagate(P, jd, [jt])[:, 0, :]
    else:
        S = P
    els = []
    for s in S:
        e = helio_elements(prop, s, jt)
        els.append([e[k] for k in ELEM_KEYS])
    els = np.array(els)
    J = np.zeros((len(ELEM_KEYS), 6))
    for k in range(6):
        J[:, k] = (els[1 + 2 * k] - els[2 + 2 * k]) / (2 * steps[k])
    # angle wrap
    for j in (3, 4, 5):
        d = els[1:, j] - els[0, j]
        els[1:, j] = els[0, j] + (d + 180) % 360 - 180
    for k in range(6):
        J[:, k] = (els[1 + 2 * k] - els[2 + 2 * k]) / (2 * steps[k])
    Ce = J @ cov @ J.T
    e0 = helio_elements(prop, S[0], jt)
    sig = {k: float(np.sqrt(max(Ce[j, j], 0))) for j, k in enumerate(ELEM_KEYS)}
    return e0, sig, Ce
