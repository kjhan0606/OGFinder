"""Child process: N-body differential correction for one or more nightlink groups.

The parent linker (``moving.nightlink``) does not import this module.  Run it as
``python -m moving.nbody_worker`` with one JSON object on stdin and one JSON
object on stdout.  rebound and assist are imported only when they are already
installed.  Otherwise the worker uses the external CODES integrator
(``neo_orbit_calculator`` under ``OGF_CODES_ROOT`` or ``~/BACKUP/3.5ST``):
Fortran Dormand-Prince, JPL DE440s, SB441-N16, full 1PN and J2/J4/J6.
That tree is not shipped with OGFinder and is never imported into the web process.
"""
import json
import os
import sys
import traceback
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")


def _root_of_repo():
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def resolve_codes_root(req):
    if req.get("use_default") is False and not req.get("codes_root"):
        return None
    given = req.get("codes_root") or os.environ.get("OGF_CODES_ROOT")
    if given:
        return os.path.abspath(os.path.expanduser(given))
    default = os.path.expanduser("~/BACKUP/3.5ST")
    kern = os.path.join(default, "neo_orbit_calculator", "kernels", "de440s.bsp")
    if os.path.isfile(kern):
        return default
    return None


def assist_ready():
    try:
        from moving import orbit as O
        return bool(O.assist_available())
    except Exception:
        return False


def _codes(root):
    if root not in sys.path:
        sys.path.insert(0, root)
    from neo_orbit_calculator.core import (  # noqa: E402
        AU_KM, C_KM_S, DAY_S, DE440Environment, ForceModel, _rhs_factory,
    )
    from neo_orbit_calculator.fortran_backend import FortranIntegrator
    kernels = os.path.join(root, "neo_orbit_calculator", "kernels")
    env = DE440Environment(Path(kernels))
    integ = FortranIntegrator(env)
    model = ForceModel()
    return env, integ, model, _rhs_factory, AU_KM, C_KM_S, DAY_S


def _radec(d):
    nrm = np.linalg.norm(d, axis=-1)
    ra = np.degrees(np.arctan2(d[..., 1], d[..., 0])) % 360.0
    dec = np.degrees(np.arcsin(np.clip(d[..., 2] / nrm, -1.0, 1.0)))
    return ra, dec


def _sky(ra_obs, dec_obs, ra_mod, dec_mod):
    dx = ((ra_obs - ra_mod + 180.0) % 360.0 - 180.0) * np.cos(np.radians(dec_obs)) * 3600.0
    dy = (dec_obs - dec_mod) * 3600.0
    return dx, dy


class CodesEngine:
    def __init__(self, root):
        self.env, self.integ, self.model, self.rhs_factory, self.AU, self.C, self.DAY = _codes(root)
        self.force_model = (
            "CODES Fortran real%d + DE440s + SB441-N16 + full 1PN + J2/J4/J6"
            % self.integ.precision_digits
        )

    def _propagate(self, state, et0, ets):
        ets = np.asarray(ets, float)
        rel = ets - et0
        if np.any(rel < -1e-3):
            raise RuntimeError("CODES Fortran integrator only steps forward of the epoch")
        order = np.argsort(rel)
        rel_s = rel[order]
        prepend = rel_s[0] > 1e-6
        times = np.r_[0.0, rel_s] if prepend else rel_s.copy()
        if not prepend:
            times[0] = 0.0
        for i in range(1, len(times)):
            if times[i] < times[i - 1]:
                times[i] = times[i - 1]
        out, _ = self.integ.propagate(
            np.asarray(state, float), et0, np.ascontiguousarray(times, float), self.model,
            rtol=1e-12, atol_position_km=1e-3, atol_velocity_kms=1e-9, max_step_days=2.0)
        if prepend:
            out = out[1:]
        inv = np.empty(len(order), int)
        inv[order] = np.arange(len(order))
        return out[inv]

    def _scipy_to(self, state, et_from, et_to):
        from scipy.integrate import solve_ivp
        rhs = self.rhs_factory(self.env, et_from, self.model)
        dt = float(et_to - et_from)
        sol = solve_ivp(rhs, (0.0, dt), np.asarray(state, float), t_eval=[dt], method="DOP853",
                        rtol=1e-11, atol=1e-4, max_step=0.5 * self.DAY)
        if not sol.success:
            raise RuntimeError(sol.message)
        return sol.y[:, -1]

    def bary_at(self, state_helio, mjd_from, mjd_to):
        from moving.util import utc_mjd_to_tdb_jd
        jd_from = float(utc_mjd_to_tdb_jd([mjd_from])[0])
        jd_to = float(utc_mjd_to_tdb_jd([mjd_to])[0])
        et_from = self.env.jd_to_et(jd_from)
        et_to = self.env.jd_to_et(jd_to)
        sun0 = self.env.state("SUN", et_from)
        bary = np.empty(6)
        bary[:3] = np.asarray(state_helio[:3], float) * self.AU + sun0[:3]
        bary[3:] = np.asarray(state_helio[3:], float) * self.AU / self.DAY + sun0[3:]
        if et_to > et_from + 1e-3:
            bary = self._propagate(bary, et_from, np.array([et_to]))[0]
        elif et_to < et_from - 1e-3:
            bary = self._scipy_to(bary, et_from, et_to)
        return bary, et_to

    def helio_of(self, bary, et):
        sun = self.env.state("SUN", et)
        h = np.empty(6)
        h[:3] = (bary[:3] - sun[:3]) / self.AU
        h[3:] = (bary[3:] - sun[3:]) * self.DAY / self.AU
        return h

    def observers(self, mjd, codes):
        from moving.util import utc_mjd_to_tdb_jd
        mjd = np.asarray(mjd, float)
        jd = utc_mjd_to_tdb_jd(mjd)
        ets = np.array([self.env.jd_to_et(float(j)) for j in jd])
        pos = np.zeros((len(mjd), 3))
        for i, et in enumerate(ets):
            pos[i] = self.env.state("EARTH", float(et))[:3]
        for code in set(codes):
            if code in (None, "", "500"):
                continue
            from moving.obs import geocentric_gcrs_km
            idx = [i for i, c in enumerate(codes) if c == code]
            pos[idx] += geocentric_gcrs_km(code, mjd[idx])
        return ets, pos

    def predict(self, states, et0, ets, obs_pos, tau_seconds=None):
        """states (n,6) km, km/s at et0.  Returns ra (n,m), dec (n,m), tau seconds (m,)."""
        states = np.atleast_2d(np.asarray(states, float))
        ets = np.asarray(ets, float)
        tau = np.full(len(ets), 0.005 * self.DAY) if tau_seconds is None else np.asarray(tau_seconds, float).copy()
        for _ in range(2):
            S = self._propagate(states[0], et0, ets - tau)
            tau = np.linalg.norm(S[:, :3] - obs_pos, axis=1) / self.C
        te = ets - tau
        ra = np.zeros((len(states), len(ets)))
        dec = np.zeros_like(ra)
        for i, st in enumerate(states):
            S = self._propagate(st, et0, te)
            tau_i = np.linalg.norm(S[:, :3] - obs_pos, axis=1) / self.C
            xx = S[:, :3] + S[:, 3:] * (tau - tau_i)[:, None]
            ra[i], dec[i] = _radec(xx - obs_pos)
        return ra, dec, tau


def _elements(helio):
    from moving import kepler as K
    el = K.state_to_elements(helio[:3], helio[3:])
    return dict(a=float(el["a"]), e=float(el["e"]), inc_deg=float(el["i"]), q=float(el["q"]),
                r_au=float(np.linalg.norm(helio[:3])))


def _fit_codes(engine, group, floor_arcsec, max_iter):
    mjd = np.asarray(group["mjd"], float)
    ra = np.asarray(group["ra"], float)
    dec = np.asarray(group["dec"], float)
    sig = np.hypot(np.asarray(group["sig"], float), float(floor_arcsec))
    codes = list(group["code"])
    mjd_epoch = float(np.min(mjd) - 0.1)
    bary, et0 = engine.bary_at(np.asarray(group["state_helio"], float), float(group["mjd_ref"]), mjd_epoch)
    ets, obs_pos = engine.observers(mjd, codes)
    steps = np.array([10.0, 10.0, 10.0, 1e-5, 1e-5, 1e-5])
    eye = np.eye(6)
    state = np.asarray(bary, float)
    tau = None
    lam = 1e-3
    accepted = 0
    converged = False
    last = None
    for it in range(int(max_iter)):
        P = [state]
        for k in range(6):
            P.append(state + steps[k] * eye[k])
            P.append(state - steps[k] * eye[k])
        ra_m, dec_m, tau = engine.predict(np.vstack(P), et0, ets, obs_pos, tau)
        dx, dy = _sky(ra, dec, ra_m[0], dec_m[0])
        res = np.empty(2 * len(mjd))
        res[0::2] = dx
        res[1::2] = dy
        J = np.zeros((2 * len(mjd), 6))
        cosd = np.cos(np.radians(dec))
        for k in range(6):
            pr, mr = ra_m[1 + 2 * k], ra_m[2 + 2 * k]
            J[0::2, k] = ((pr - mr + 180.0) % 360.0 - 180.0) * cosd * 3600.0 / (2.0 * steps[k])
            J[1::2, k] = (dec_m[1 + 2 * k] - dec_m[2 + 2 * k]) * 3600.0 / (2.0 * steps[k])
        w = np.repeat(1.0 / sig ** 2, 2)
        Js = J * steps
        N = Js.T @ (Js * w[:, None])
        b = Js.T @ (w * res)
        chi2 = float(np.sum(w * res ** 2))
        try:
            dy_s = np.linalg.solve(N + lam * np.eye(6), b)
        except np.linalg.LinAlgError:
            dy_s = np.linalg.lstsq(N + lam * np.eye(6), b, rcond=None)[0]
        trial = state + dy_s * steps
        ra_t, dec_t, tau_t = engine.predict(trial[None, :], et0, ets, obs_pos, tau)
        dx2, dy2 = _sky(ra, dec, ra_t[0], dec_t[0])
        r2 = np.empty_like(res)
        r2[0::2] = dx2
        r2[1::2] = dy2
        chi2_new = float(np.sum(w * r2 ** 2))
        if np.isfinite(chi2_new) and chi2_new <= chi2 * (1.0 + 1e-9):
            state = trial
            tau = tau_t
            lam = max(lam / 3.0, 1e-6)
            accepted += 1
            last = (dx2, dy2, chi2_new)
            # dy is in units of the finite-difference step (10 km, 1e-5 km/s).  0.05 of that step is ~0.5 km.
            if it > 0 and (np.max(np.abs(dy_s)) < 0.05 or abs(chi2 - chi2_new) < 1e-3 * max(chi2, 1.0)):
                converged = True
                break
        else:
            lam = min(lam * 10.0, 1e6)
    if last is None:
        return {"id": group.get("id"), "status": "not_converged", "message": "no accepted step"}
    dx, dy, chi2 = last
    n = len(mjd)
    ndof = 2 * n - 6
    helio = engine.helio_of(state, et0)
    mx = float(max(np.max(np.abs(dx / sig)), np.max(np.abs(dy / sig))))
    return {
        "id": group.get("id"),
        "status": "converged" if converged else "not_converged",
        "iterations": int(it + 1),
        "accepted_steps": int(accepted),
        "chi2": float(chi2),
        "ndof": int(ndof),
        "chi2_red": float(chi2 / max(ndof, 1)),
        "rms_arcsec": float(np.sqrt(np.mean(dx ** 2 + dy ** 2) / 2.0)),
        "max_resid": mx,
        "n_used": int(n),
        "n_obs": int(n),
        "mjd_epoch": mjd_epoch,
        "state_helio": helio.tolist(),
        "elements": _elements(helio),
    }


def _fit_assist(group, floor_arcsec, max_iter):
    from moving import kepler as K
    from moving import orbit as O
    from moving.util import utc_mjd_to_tdb_jd
    prop = O.Propagator()
    mjd = np.asarray(group["mjd"], float)
    ra = np.asarray(group["ra"], float)
    dec = np.asarray(group["dec"], float)
    sig = np.hypot(np.asarray(group["sig"], float), float(floor_arcsec))
    codes = [("500" if c in (None, "") else c) for c in group["code"]]
    jd = utc_mjd_to_tdb_jd(mjd)
    jd0 = float(np.mean(jd))
    mjd0 = float(np.mean(mjd))
    r, v = K.propagate_2body(group["state_helio"][:3], group["state_helio"][3:], mjd0 - float(group["mjd_ref"]))
    sun = prop.sun(jd0)
    state = np.r_[np.asarray(r, float) + sun[:3], np.asarray(v, float) + sun[3:]]
    obs = O.Obs(mjd, ra, dec, sig, sig, codes)
    res = O.differential_correction(prop, state, jd0, obs, reject=True, max_iter=int(max_iter))
    active = np.asarray(res["active"], bool)
    if int(active.sum()) < 4:
        return {"id": group.get("id"), "status": "not_converged", "message": "fewer than 4 observations left"}
    sun = prop.sun(jd0)
    helio = np.r_[res["state"][:3] - sun[:3], res["state"][3:] - sun[3:]]
    el = O.helio_elements(prop, res["state"], jd0)
    dra = np.asarray(res["residual_ra"], float)[active]
    dde = np.asarray(res["residual_dec"], float)[active]
    sig_a = sig[active]
    return {
        "id": group.get("id"),
        "status": "converged",
        "iterations": None,
        "chi2": float(res["chi2"]),
        "ndof": int(res["ndof"]),
        "chi2_red": float(res["chi2_red"]),
        "rms_arcsec": float(np.sqrt(np.mean(dra ** 2 + dde ** 2) / 2.0)),
        "max_resid": float(max(np.max(np.abs(dra / sig_a)), np.max(np.abs(dde / sig_a)))),
        "n_used": int(res["n_used"]),
        "n_obs": int(res["n_total"]),
        "mjd_epoch": mjd0,
        "state_helio": helio.tolist(),
        "elements": dict(a=float(el["a"]), e=float(el["e"]), inc_deg=float(el["i"]), q=float(el["q"]),
                         r_au=float(np.linalg.norm(helio[:3]))),
    }


def _observe(engine, req):
    mjd = np.asarray(req["mjd"], float)
    codes = list(req.get("code") or ["500"] * len(mjd))
    if len(codes) == 1 and len(mjd) > 1:
        codes = codes * len(mjd)
    # Epoch sits before every emission time.  The Fortran propagator does not step backward.
    mjd_epoch = float(np.min(mjd) - 0.1)
    bary, et0 = engine.bary_at(np.asarray(req["state_helio"], float), float(req["mjd_ref"]), mjd_epoch)
    ets, obs_pos = engine.observers(mjd, codes)
    ra, dec, _ = engine.predict(bary[None, :], et0, ets, obs_pos)
    return {"status": "ok", "ra": ra[0].tolist(), "dec": dec[0].tolist(), "mjd_epoch": mjd_epoch,
            "state_helio": engine.helio_of(bary, et0).tolist()}


def dispatch(req):
    cmd = req.get("cmd", "fit")
    allow_assist = bool(req.get("allow_assist", True))
    root = resolve_codes_root(req)
    use_assist = allow_assist and assist_ready()
    if cmd == "probe":
        if use_assist:
            return {"status": "ok", "backend": "assist", "force_model": "ASSIST(DE440+16 asteroids+GR+J2/J3/J4)"}
        if root and os.path.isfile(os.path.join(root, "neo_orbit_calculator", "kernels", "de440s.bsp")):
            engine = CodesEngine(root)
            return {"status": "ok", "backend": "codes", "force_model": engine.force_model, "codes_root": root}
        return {"status": "unavailable", "reason": "neither ASSIST+DE440 nor a CODES tree with de440s.bsp is available"}
    if not use_assist and not (root and os.path.isfile(os.path.join(root, "neo_orbit_calculator", "kernels", "de440s.bsp"))):
        return {"status": "unavailable", "reason": "neither ASSIST+DE440 nor a CODES tree with de440s.bsp is available"}
    engine = None if use_assist else CodesEngine(root)
    if cmd == "observe":
        if engine is None:
            return {"status": "error", "message": "observe is implemented for the CODES backend only"}
        return _observe(engine, req)
    if cmd == "propagate_helio":
        if engine is None:
            return {"status": "error", "message": "propagate_helio is implemented for the CODES backend only"}
        bary, et = engine.bary_at(np.asarray(req["state_helio"], float), float(req["mjd_ref"]), float(req["mjd_out"]))
        helio = engine.helio_of(bary, et)
        return {"status": "ok", "mjd_epoch": float(req["mjd_out"]), "state_helio": helio.tolist(), "elements": _elements(helio)}
    if cmd != "fit":
        return {"status": "error", "message": "unknown cmd %s" % cmd}
    floor = float(req.get("floor_arcsec", 0.05))
    max_iter = int(req.get("max_iter", 12))
    rows = []
    for group in req.get("groups", []):
        try:
            if use_assist:
                rows.append(_fit_assist(group, floor, max_iter))
            else:
                rows.append(_fit_codes(engine, group, floor, max_iter))
        except Exception as exc:
            rows.append({"id": group.get("id"), "status": "error", "message": "%s: %s" % (type(exc).__name__, exc)})
    return {
        "status": "ok",
        "backend": "assist" if use_assist else "codes",
        "force_model": "ASSIST(DE440+16 asteroids+GR+J2/J3/J4)" if use_assist else engine.force_model,
        "codes_root": None if use_assist else root,
        "groups": rows,
    }


def main():
    try:
        req = json.loads(sys.stdin.read() or "{}")
        out = dispatch(req)
    except Exception:
        out = {"status": "error", "message": traceback.format_exc(limit=8)}
    json.dump(out, sys.stdout)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
