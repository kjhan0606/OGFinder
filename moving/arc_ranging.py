"""Orbital-element distribution for tracklets the night linker already kept.

``moving.iod.ranging`` draws topocentric range and range-rate inside the
Sun-bound admissible region, optimises the attributable, and weights each
draw by its astrometric chi-squared (Laplace marginalisation, flat prior in
range and range-rate).  The Jacobian from that prior to (a, e, i) is not
applied.  This is the simplified sampler already in this tree.  It is not
the OpenOrb program of Granvik et al. 2009, and it does not write ``.sor``
or ``.orb`` files.  On an arc of at least a day, when a 6-parameter fit
determines the orbit, the reported cloud is the local Laplace sample of that
fit.  A uniform draw over the whole admissible region misses that valley.

The two-body point fit stored on the group is left as it is.  A subset of
the accepted samples can be propagated by the child process
(``python -m moving.nbody_worker``, command ``propagate_many``).  The default
integrator is the in-tree model.  ``allow_assist=False`` keeps the external
CODES path when that tree is present.  A missing backend leaves the two-body
distribution in place.  This module does not import rebound, assist,
spiceypy, or neo_orbit_calculator.
"""
import numpy as np
from scipy.optimize import least_squares

from . import iod, kepler as K, nightlink as N

LIKE_MIN = 1e-3
# Steps used to scale the 6-parameter Hessian.  Raw eigenvalues mix degrees with AU.
_SCALE6 = np.array([1e-5, 1e-5, 1e-8, 1e-8, 0.01, 1e-4])
NOTE = (
    "simplified statistical ranging (moving.iod.ranging): flat prior in range and "
    "range-rate, Laplace weight on the attributable; the Jacobian to (a, e, i) is "
    "not applied; this is not OpenOrb"
)
LAPLACE_NOTE = (
    "local Laplace sample of the 6-parameter fit to the same observations. "
    "A uniform draw over the whole admissible region does not land in this valley "
    "once the arc determines the range. Not OpenOrb"
)


def stacked_observations(members):
    """Heliocentric observer states aligned with every detection.  Sun at the origin, as in the nightlink two-body model."""
    mjd, ra, dec, sig, pos, vel = [], [], [], [], [], []
    for t in members:
        p, v = t.obs()
        mjd.append(np.asarray(t.mjd, float))
        ra.append(np.asarray(t.ra, float))
        dec.append(np.asarray(t.dec, float))
        sig.append(np.asarray(t.sig, float))
        pos.append(np.asarray(p, float))
        vel.append(np.asarray(v, float))
    mjd = np.concatenate(mjd)
    ra = np.concatenate(ra)
    dec = np.concatenate(dec)
    sig = np.concatenate(sig)
    pos = np.concatenate(pos, axis=0)
    vel = np.concatenate(vel, axis=0)
    order = np.argsort(mjd)
    obs = np.concatenate([pos[order], vel[order]], axis=1)
    sun = np.zeros((len(order), 6))
    return mjd[order], ra[order], dec[order], sig[order], obs, sun


def _quantiles(x, w):
    return [float(v) for v in iod.weighted_quantiles(np.asarray(x, float), np.asarray(w, float))]


def _best_rms(R, mjd, ra, dec, obs):
    b = int(R["best"])
    state = np.asarray(R["state"][b], float)
    op = obs[:, :3]
    pra, pde, _ = N.predict(state[:3], state[3:], R["jd_ref"], mjd, op)
    dx, dy = N.sky_diff(ra, dec, pra, pde)
    return float(np.sqrt(np.mean(dx ** 2 + dy ** 2) / 2.0))


def _params_from_state(r, v, op, ov):
    los = np.asarray(r, float) - np.asarray(op, float)
    rho = float(np.linalg.norm(los))
    los = los / rho
    rel = np.asarray(v, float) - np.asarray(ov, float)
    rhod = float(los @ rel)
    tang = rel - rhod * los
    ra = float(np.arctan2(los[1], los[0]))
    dec = float(np.arcsin(np.clip(los[2], -1.0, 1.0)))
    ea = np.array([-np.sin(ra), np.cos(ra), 0.0])
    ed = np.array([-np.sin(dec) * np.cos(ra), -np.sin(dec) * np.sin(ra), np.cos(dec)])
    return np.array([
        np.degrees(ra) % 360.0, np.degrees(dec),
        float(tang @ ea) / rho, float(tang @ ed) / rho, rho, rhod,
    ])


def _resid6(p, tref, mjd, ra, dec, sig_deg, op, ov, obs_h):
    n = len(mjd)
    if not np.isfinite(p).all() or p[4] <= 0.02:
        return np.full(2 * n, 1e3)
    try:
        r_, v_ = iod.attributable_to_state(p[0], p[1], p[2], p[3], p[4], p[5], op, ov)
        pred, _ = iod.propagate_obs((r_, v_), tref, mjd, obs_h)
    except Exception:
        return np.full(2 * n, 1e3)
    if not np.all(np.isfinite(pred)):
        return np.full(2 * n, 1e3)
    dr = ((ra - pred[:, 0] + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec)) / sig_deg
    dd = (dec - pred[:, 1]) / sig_deg
    return np.r_[dr, dd]


def _sky_rms(r, v, tref, mjd, ra, dec, op):
    pra, pde, _ = N.predict(r, v, tref, mjd, op)
    dx, dy = N.sky_diff(ra, dec, pra, pde)
    return float(np.sqrt(np.mean(dx ** 2 + dy ** 2) / 2.0))


def _try_local(mjd, ra, dec, sig, obs, n_samples, seed, code="500"):
    """Laplace sample around a 6-parameter fit.  None when that fit does not determine the orbit."""
    span = float(mjd[-1] - mjd[0])
    if span < 1.0 or len(mjd) < 4:
        return None
    tref = float(np.mean(mjd))
    k = int(np.argmin(np.abs(mjd - tref)))
    op, ov = obs[k, :3], obs[k, 3:]
    obs_h = obs[:, :3]
    sig_deg = np.asarray(sig, float) / 3600.0
    trk = N.Tracklet(mjd, ra, dec, sig, code=code or "500", tid=0, night=0)
    R, V = N.hypotheses(trk, n_rho=12, n_rd=9)[:2]
    if len(R) == 0:
        return None
    scores = []
    for i in range(len(R)):
        pra, pde, _ = N.predict(R[i], V[i], trk.t0, mjd, obs_h)
        dx, dy = N.sky_diff(ra, dec, pra, pde)
        scores.append(float(np.mean(dx ** 2 + dy ** 2)))
    order = np.argsort(scores)[:4]
    scale = _SCALE6
    best = None
    for i in order:
        rr, vv = K.propagate_2body(R[i], V[i], tref - trk.t0)
        p0 = _params_from_state(rr, vv, op, ov)
        try:
            sol = least_squares(lambda p: _resid6(p, tref, mjd, ra, dec, sig_deg, op, ov, obs_h),
                                p0, x_scale=scale, method="lm", max_nfev=80)
        except Exception:
            continue
        chi2 = float(np.sum(sol.fun ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, sol)
    if best is None:
        return None
    chi2, sol = best
    r_mode, v_mode = iod.attributable_to_state(*sol.x, op, ov)
    rms = _sky_rms(r_mode, v_mode, tref, mjd, ra, dec, obs_h)
    if rms > 1.0:
        return None
    # Covariance is taken in scaled parameters.  The raw Hessian is ill-conditioned
    # because right ascension in degrees sits next to range-rate in AU/day.
    H = sol.jac.T @ sol.jac
    Hs = (_SCALE6[:, None] * H) * _SCALE6[None, :]
    evals, evec = np.linalg.eigh(0.5 * (Hs + Hs.T))
    if not np.all(np.isfinite(evals)) or np.min(evals) <= 0:
        return None
    draw_scale = 1.0 / np.sqrt(evals)
    rng = np.random.default_rng(int(seed))
    raw = sol.x + (rng.normal(size=(int(n_samples), 6)) @ (evec * draw_scale).T) * _SCALE6
    states, chi, els = [], [], []
    for p in raw:
        if p[4] <= 0.02:
            continue
        try:
            rr, vv = iod.attributable_to_state(p[0], p[1], p[2], p[3], p[4], p[5], op, ov)
            el = K.state_to_elements(rr, vv)
            fun = _resid6(p, tref, mjd, ra, dec, sig_deg, op, ov, obs_h)
        except Exception:
            continue
        if not np.isfinite(el["a"]) or el["a"] <= 0:
            continue
        states.append(np.r_[rr, vv])
        chi.append(float(np.sum(fun ** 2)))
        els.append(el)
    if len(states) < max(10, int(n_samples) // 5):
        return None
    chi = np.asarray(chi, float)
    # The draws already come from the Gaussian approximation, so the reported
    # quantiles are unweighted.  chi2 stays on the pack for the propagation cut.
    w = np.ones(len(states)) / len(states)
    a = np.array([e["a"] for e in els])
    e = np.array([e["e"] for e in els])
    inc = np.array([e["i"] for e in els])
    q = np.array([e["q"] for e in els])
    qa = _quantiles(a, w)
    if qa[2] - qa[0] > 1.0:
        return None
    cls = {}
    for ai, ei, ii, wi in zip(a, e, inc, w):
        c = K.orbit_class(ai, ei, ii)
        cls[c] = cls.get(c, 0.0) + float(wi)
    el_mode = K.state_to_elements(r_mode, v_mode)
    summary = {
        "status": "ok",
        "sampler": "laplace",
        "note": LAPLACE_NOTE,
        "n_obs": int(len(mjd)),
        "n_samples": int(len(states)),
        "neff": float(len(states)),
        "mjd_epoch": tref,
        "chi2_min": float(chi2),
        "best_rms_arcsec": rms,
        "best": {
            "a": float(el_mode["a"]), "e": float(el_mode["e"]),
            "inc_deg": float(el_mode["i"]), "q": float(el_mode["q"]),
        },
        "quantiles_16_50_84": {
            "a": qa, "e": _quantiles(e, w), "inc_deg": _quantiles(inc, w), "q": _quantiles(q, w),
        },
        "class_probs": {str(k): float(v) for k, v in cls.items() if v > 0},
    }
    pack = {
        "state": np.asarray(states, float),
        "chi2": chi,
        "w": w,
        "jd_ref": tref,
        "mjd_last": float(np.max(mjd)),
    }
    return summary, pack


def _ranging_distribution(mjd, ra, dec, sig, obs, sun, n_samples, seed):
    R = iod.ranging(ra, dec, mjd, sig, obs, sun, n_samples=int(n_samples), rng=int(seed))
    w = np.asarray(R["w"], float)
    b = int(R["best"])
    summary = {
        "status": "ok",
        "sampler": "iod.ranging",
        "note": NOTE,
        "n_obs": int(len(mjd)),
        "n_samples": int(R["n_samples"]),
        "neff": float(R["neff"]),
        "mjd_epoch": float(R["jd_ref"]),
        "chi2_min": float(R["chi2"][b]),
        "best_rms_arcsec": _best_rms(R, mjd, ra, dec, obs),
        "best": {
            "a": float(R["a"][b]),
            "e": float(R["e"][b]),
            "inc_deg": float(R["i"][b]),
            "q": float(R["q"][b]),
        },
        "quantiles_16_50_84": {
            "a": _quantiles(R["a"], w),
            "e": _quantiles(R["e"], w),
            "inc_deg": _quantiles(R["i"], w),
            "q": _quantiles(R["q"], w),
        },
        "class_probs": {str(k): float(v) for k, v in R["class_probs"].items() if float(v) > 0},
    }
    if summary["neff"] < 2.0 and summary["best_rms_arcsec"] > 2.0:
        summary["status"] = "missed_valley"
    pack = {
        "state": np.asarray(R["state"], float),
        "chi2": np.asarray(R["chi2"], float),
        "w": w,
        "jd_ref": float(R["jd_ref"]),
        "mjd_last": float(np.max(mjd)),
    }
    return summary, pack


def distribution(members, n_samples=400, seed=1):
    """Weighted (a, e, i, q) distribution at the mean observation time.

    A short arc uses the admissible-region Monte Carlo.  An arc of at least one
    day uses a local Laplace sample when that 6-parameter fit actually determines
    the orbit.  Returns ``(summary, pack)``.  ``pack`` is not written to JSON.
    """
    mjd, ra, dec, sig, obs, sun = stacked_observations(members)
    if len(mjd) < 3:
        raise ValueError("ranging needs at least 3 observations")
    if float(mjd[-1] - mjd[0]) >= 1.0:
        codes = {getattr(t, "code", None) for t in members}
        code = codes.pop() if len(codes) == 1 else "500"
        local = _try_local(mjd, ra, dec, sig, obs, n_samples, seed, code=code or "500")
        if local is not None:
            return local
    return _ranging_distribution(mjd, ra, dec, sig, obs, sun, n_samples, seed)


def accepted_indices(pack, like_min=LIKE_MIN):
    """Samples whose astrometric likelihood is at least ``like_min`` times the best sample."""
    chi2 = np.asarray(pack["chi2"], float)
    like = np.exp(-0.5 * (chi2 - float(np.min(chi2))))
    keep = np.flatnonzero(np.isfinite(like) & (like >= like_min))
    if len(keep) == 0:
        keep = np.array([int(np.argmin(chi2))])
    return keep


def draw_states(pack, n, seed, like_min=LIKE_MIN):
    """Weighted draw, without replacement, of accepted heliocentric states."""
    keep = accepted_indices(pack, like_min)
    w = np.asarray(pack["w"], float)[keep]
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    if float(w.sum()) <= 0:
        w = np.ones(len(keep))
    w = w / w.sum()
    n_take = min(int(n), len(keep))
    pick = np.random.default_rng(int(seed)).choice(keep, size=n_take, replace=False, p=w)
    return np.asarray(pack["state"], float)[pick], keep


def propagate_accepted(pack, n_propagate, seed, codes_root=None, use_default=True, timeout=300,
                       like_min=LIKE_MIN, allow_assist=True):
    """Propagate a weighted subset of accepted samples to the last observation time."""
    from . import nbody_refine as NR

    states, keep = draw_states(pack, n_propagate, seed, like_min)
    mjd_out = float(pack["mjd_last"])
    if mjd_out <= float(pack["jd_ref"]) + 1e-4:
        mjd_out = float(pack["jd_ref"]) + 1.0
    req = {
        "cmd": "propagate_many",
        "allow_assist": bool(allow_assist),
        "use_default": bool(use_default),
        "states_helio": states.tolist(),
        "mjd_ref": float(pack["jd_ref"]),
        "mjd_out": mjd_out,
    }
    if codes_root is not None:
        req["codes_root"] = codes_root
    out = NR.call_worker(req, codes_root=codes_root, timeout=timeout)
    base = {"n_accepted": int(len(keep)), "n_requested": int(min(int(n_propagate), len(keep))), "like_min": float(like_min)}
    if out.get("status") != "ok":
        base.update(status=out.get("status", "error"), message=out.get("message") or out.get("reason") or "")
        return base
    rows = out.get("elements") or []
    if not rows:
        base.update(status="error", message="the integrator returned no elements")
        return base
    ones = np.ones(len(rows))
    base.update(
        status="ok",
        backend=out.get("backend"),
        force_model=out.get("force_model"),
        n_propagated=int(len(rows)),
        n_failed=int(out.get("n_failed") or 0),
        mjd_epoch=float(out.get("mjd_epoch")),
        quantiles_16_50_84={
            "a": _quantiles([r["a"] for r in rows], ones),
            "e": _quantiles([r["e"] for r in rows], ones),
            "inc_deg": _quantiles([r["inc_deg"] for r in rows], ones),
            "q": _quantiles([r["q"] for r in rows], ones),
        },
    )
    return base


def attach(trks, result, n_samples=400, seed=1, n_propagate=32, codes_root=None, use_default=True,
           timeout=None, allow_assist=True):
    """Add a ``ranging`` summary to each linked group and each pair.  Point-fit fields stay unchanged."""
    by_id = {t.id: t for t in trks}
    n_ok = 0
    timeout = 300 if timeout is None else timeout
    if n_propagate:
        timeout = max(int(timeout), 8 * int(n_propagate))
    for kind, bucket in (("group", result.get("groups") or []), ("pair", result.get("pairs") or [])):
        for gi, g in enumerate(bucket):
            try:
                members = [by_id[i] for i in g["ids"]]
            except KeyError as exc:
                g["ranging"] = {"status": "error", "message": "missing tracklet %s" % exc}
                continue
            try:
                summary, pack = distribution(members, n_samples=n_samples, seed=int(seed) + (0 if kind == "group" else 1000) + gi)
            except Exception as exc:
                g["ranging"] = {"status": "error", "message": "%s: %s" % (type(exc).__name__, exc)}
                continue
            if int(n_propagate) > 0:
                summary["codes"] = propagate_accepted(
                    pack, n_propagate, int(seed) + 17 + gi, codes_root=codes_root,
                    use_default=use_default, timeout=timeout, allow_assist=allow_assist)
            else:
                summary["codes"] = {"status": "skipped", "n_accepted": 0, "n_propagated": 0}
            g["ranging"] = summary
            n_ok += 1
    result["ranging"] = {
        "status": "ok" if n_ok else "empty",
        "n": int(n_ok),
        "n_samples": int(n_samples),
        "n_propagate": int(n_propagate),
        "note": NOTE,
    }
    return result
