"""Initial orbit determination.

(A classical Gauss IOD is NOT implemented; ranging() + differential correction in orbit.py is used instead.)
ranging(): admissible-region + statistical-ranging posterior for very short arcs (Milani et al. 2004 admissible
  region; Virtanen et al. 2001 / Granvik et al. 2009 ranging idea), SIMPLIFIED:
    * (rho, rho_dot) sampled from the admissible region defined by: rho>0, heliocentric two-body energy < E_max
      (bound or mildly hyperbolic cut-off set by `amax`), heliocentric distance r>r_min, optional H prior;
    * for each (rho, rho_dot) the four attributable parameters (ra, dec, ra_dot, dec_dot) are optimised against the
      observed astrometry (two-body propagation with the actual observer position of each exposure, so observer
      parallax is handled exactly);
    * weight  w = prior(a,e,i) * exp(-chi2_min/2) / sqrt(det(H))  (Laplace marginalisation over the attributable);
    * the prior is a flat prior in (rho, rho_dot) (the Jacobian to (a,e,i) is NOT applied) multiplied by an optional
      population prior in (a, e, i) that is switched off by default.
  Output: weighted samples of heliocentric (a, e, i) etc. and class probabilities.  Not validated against the MPC
  Orbit Fitting software; validated here only on JPL Horizons-generated synthetic arcs and on real HST data.
"""
import numpy as np
from scipy.optimize import least_squares
from . import kepler as K
from .util import C_AUD


def _los(ra, dec):
    a, d = np.radians(ra), np.radians(dec)
    return np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])


def attributable_to_state(ra, dec, dra, ddec, rho, rhod, op, ov):
    """ra/dec [deg], dra = d(ra*cosdec)/dt, ddec [rad/day] at the reference time; op/ov observer bary pos/vel (AU, AU/d)."""
    a, d = np.radians(ra), np.radians(dec)
    los = np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])
    e_a = np.array([-np.sin(a), np.cos(a), 0.0])
    e_d = np.array([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)])
    r = op + rho * los
    v = ov + rhod * los + rho * (dra * e_a + ddec * e_d)
    return r, v


def propagate_obs(state_helio_ref, jd_ref, jds, obs_helio):
    """Predict ra/dec for 2-body motion; state_helio_ref=(r,v) heliocentric equatorial at jd_ref; obs_helio (n,3) heliocentric
    observer positions at jds.  First-order light-time."""
    r0, v0 = state_helio_ref
    out = np.zeros((len(jds), 2))
    rho = np.zeros(len(jds))
    for i, (t, o) in enumerate(zip(jds, obs_helio)):
        tau = 0.0
        for _ in range(3):
            r, v = K.propagate_2body(r0, v0, t - tau - jd_ref)
            d = r - o
            tau = np.linalg.norm(d) / C_AUD
        a = np.degrees(np.arctan2(d[1], d[0])) % 360.0
        dd = np.degrees(np.arcsin(d[2] / np.linalg.norm(d)))
        out[i] = (a, dd); rho[i] = np.linalg.norm(d)
    return out, rho


def ranging(ra, dec, jd, sig_arcsec, obs_bary, sun_bary, n_samples=3000, rho_range=(0.05, 8.0), amax=60.0,
            rng=None, h_prior=None, jd_ref=None, verbose=False):
    """Statistical ranging posterior.  Inputs: ra, dec arrays [deg] at jd (TDB, days); sig_arcsec per obs; obs_bary (n,6)
    barycentric observer state; sun_bary (n,6) barycentric Sun state at each obs time.  The two-body heliocentric force
    centre is the Sun (barycentre offset neglected -> fine for arcs <~ weeks, flagged in docs).
    Returns dict with arrays of samples: a, e, i, q, rho, rhod, w (normalised weights), chi2, state (n,6 heliocentric at jd_ref),
    jd_ref, class_probs (dict), neff."""
    rng = np.random.default_rng(1 if rng is None else rng)
    ra = np.asarray(ra); dec = np.asarray(dec); jd = np.asarray(jd, float)
    n = len(jd)
    jd_ref = float(np.mean(jd)) if jd_ref is None else jd_ref
    sig = np.broadcast_to(np.asarray(sig_arcsec, float), (n,)) / 3600.0
    obs_h = obs_bary[:, :3] - sun_bary[:, :3]
    k_ref = int(np.argmin(np.abs(jd - jd_ref)))
    op = obs_bary[k_ref, :3] - sun_bary[k_ref, :3]; ov = obs_bary[k_ref, 3:] - sun_bary[k_ref, 3:]
    # initial attributable by linear fit in tangent plane (arcsec/day -> rad/day)
    from .tracklet import fit_tracklet
    f = fit_tracklet(jd, ra, dec, np.broadcast_to(sig * 3600.0, (n,)), t_ref=jd_ref)
    a0, d0 = f["ra0"], f["dec0"]
    dra0 = f["vx_asd"] * np.pi / 180 / 3600.0 / np.cos(np.radians(d0)) * np.cos(np.radians(d0))
    dde0 = f["vy_asd"] * np.pi / 180 / 3600.0
    # admissible-region sampling: draw rho (log-uniform), rhod uniform in the bound-orbit interval
    mu = K.MU
    cand = []
    tries = 0
    while len(cand) < n_samples and tries < 200 * n_samples:
        tries += 1
        rho = np.exp(rng.uniform(np.log(rho_range[0]), np.log(rho_range[1])))
        # energy bound: |v|^2/2 - mu/r < -mu/(2 amax)  -> restrict rhod
        # solve quadratic for rhod
        los = _los(a0, d0)
        r = op + rho * los
        rn = np.linalg.norm(r)
        if rn < 0.3:
            continue
        # v = ov + rhod*los + rho*(mu_vec): |v|^2 = rhod^2 + 2 rhod (los.(ov+rho mu_vec)) + |ov + rho mu_vec|^2
        a_ = np.radians(a0); d_ = np.radians(d0)
        e_a = np.array([-np.sin(a_), np.cos(a_), 0.0]); e_d = np.array([-np.sin(d_) * np.cos(a_), -np.sin(d_) * np.sin(a_), np.cos(d_)])
        w0 = ov + rho * (dra0 * e_a + dde0 * e_d)
        B = los @ w0; C = w0 @ w0 - 2 * mu / rn + mu / amax
        disc = B * B - C
        if disc <= 0:
            continue
        lo, hi = -B - np.sqrt(disc), -B + np.sqrt(disc)
        rhod = rng.uniform(lo, hi)
        cand.append((rho, rhod, hi - lo))
    cand = np.array(cand)
    res = []
    for (rho, rhod, width) in cand:
        def resid(p):
            r_, v_ = attributable_to_state(p[0], p[1], p[2], p[3], rho, rhod, op, ov)
            pred, _ = propagate_obs((r_, v_), jd_ref, jd, obs_h)
            dr = ((ra - pred[:, 0] + 180) % 360 - 180) * np.cos(np.radians(dec)) / sig
            dd = (dec - pred[:, 1]) / sig
            return np.r_[dr, dd]
        p0 = np.array([a0, d0, dra0, dde0])
        sc = np.array([1e-5, 1e-5, 1e-8, 1e-8])
        try:
            sol = least_squares(resid, p0, x_scale=sc, method="lm", max_nfev=60)
        except Exception:
            continue
        chi2 = float(np.sum(sol.fun ** 2))
        Jm = sol.jac
        H = Jm.T @ Jm
        sd = np.linalg.slogdet(H)[1] if np.all(np.isfinite(H)) else np.inf
        r_, v_ = attributable_to_state(sol.x[0], sol.x[1], sol.x[2], sol.x[3], rho, rhod, op, ov)
        el = K.state_to_elements(r_, v_)
        res.append((rho, rhod, chi2, sd, el["a"], el["e"], el["i"], el["q"], np.r_[r_, v_], width))
    if not res:
        raise RuntimeError("ranging produced no admissible samples")
    chi2 = np.array([x[2] for x in res]); sd = np.array([x[3] for x in res])
    # flat prior in (rho, rhod) -> importance weights; the log-uniform rho draw means density ~ 1/rho
    rho_s = np.array([x[0] for x in res])
    width_s = np.array([x[9] for x in res])
    # proposal density is  p(rho) ~ 1/rho (log-uniform) x 1/width(rho) (uniform rhod): importance weight for a flat (rho, rhod) prior
    logw = -0.5 * (chi2 - chi2.min()) - 0.5 * sd + np.log(rho_s) + np.log(width_s)
    w = np.exp(logw - logw.max()); w /= w.sum()
    neff = 1.0 / np.sum(w ** 2)
    a = np.array([x[4] for x in res]); e = np.array([x[5] for x in res]); i = np.array([x[6] for x in res]); q = np.array([x[7] for x in res])
    cls = {}
    for k in range(len(a)):
        c = K.orbit_class(a[k], e[k], i[k])
        cls[c] = cls.get(c, 0.0) + w[k]
    best = int(np.argmin(chi2))
    return dict(a=a, e=e, i=i, q=q, rho=rho_s, rhod=np.array([x[1] for x in res]), w=w, chi2=chi2, neff=float(neff),
                state=np.array([x[8] for x in res]), jd_ref=jd_ref, class_probs=cls, best=best,
                n_samples=len(res), attributable=dict(ra=a0, dec=d0, dra=dra0, ddec=dde0))


def refine_two_body(R, ra, dec, jd, sig_arcsec, obs_bary, sun_bary, top=6):
    """Polish the best ranging samples with a full 6-parameter (attributable + rho, rho_dot) least-squares fit of the 2-body model.
    Returns list of dicts(state_helio_eq(6), chi2) sorted by chi2.  Used to seed the N-body differential correction: a 4-day
    arc with 0.2" astrometry constrains (rho, rho_dot) far better than a few hundred random admissible-region draws can sample."""
    n = len(jd); sig = np.broadcast_to(np.asarray(sig_arcsec, float), (n,)) / 3600.0
    obs_h = obs_bary[:, :3] - sun_bary[:, :3]
    k_ref = int(np.argmin(np.abs(jd - R["jd_ref"])))
    op = obs_bary[k_ref, :3] - sun_bary[k_ref, :3]; ov = obs_bary[k_ref, 3:] - sun_bary[k_ref, 3:]
    a0 = R["attributable"]; out = []
    order = np.argsort(R["chi2"])[:top]
    for i in order:
        def resid(p):
            r_, v_ = attributable_to_state(p[0], p[1], p[2], p[3], p[4], p[5], op, ov)
            if np.linalg.norm(r_) < 0.05:
                return np.full(2 * n, 1e6)
            try:
                pred, _ = propagate_obs((r_, v_), R["jd_ref"], jd, obs_h)
            except Exception:
                return np.full(2 * n, 1e6)
            dr = ((ra - pred[:, 0] + 180) % 360 - 180) * np.cos(np.radians(dec)) / sig
            dd = (dec - pred[:, 1]) / sig
            return np.r_[dr, dd]
        p0 = np.array([a0["ra"], a0["dec"], a0["dra"], a0["ddec"], R["rho"][i], R["rhod"][i]])
        sc = np.array([1e-5, 1e-5, 1e-8, 1e-8, 0.01, 1e-4])
        try:
            sol = least_squares(resid, p0, x_scale=sc, method="lm", max_nfev=400)
        except Exception:
            continue
        r_, v_ = attributable_to_state(*sol.x, op, ov)
        out.append(dict(state=np.r_[r_, v_], chi2=float(np.sum(sol.fun ** 2)), rho=float(sol.x[4]), rhod=float(sol.x[5])))
    out.sort(key=lambda d: d["chi2"])
    return out


def weighted_quantiles(x, w, qs=(0.16, 0.5, 0.84)):
    idx = np.argsort(x); xs = x[idx]; cw = np.cumsum(w[idx]); cw /= cw[-1]
    return [float(np.interp(q, cw, xs)) for q in qs]
