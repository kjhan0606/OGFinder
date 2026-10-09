"""N-body refinement of nightlink groups, in a separate process.

``link_nights`` stays a Sun-only two-body fit so this module's import does not
load rebound, assist, or spiceypy.  ``refine_result`` sends the linked groups
to ``python -m moving.nbody_worker``.  The worker calls
``moving.orbit.differential_correction`` on the in-tree solar-system model,
and otherwise the external CODES integrator (see the worker module).
Groups the N-body fit rejects leave ``groups`` and are listed in
``nbody_rejected``.  A missing backend leaves the two-body result unchanged.
"""
import json
import os
import subprocess
import sys

import numpy as np

from . import nightlink as N

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _child_env(codes_root):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "2"
    env["OPENBLAS_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    env["NUMEXPR_NUM_THREADS"] = "2"
    parts = [p for p in (codes_root, env.get("PYTHONPATH")) if p]
    if parts:
        env["PYTHONPATH"] = os.pathsep.join(parts)
    return env


def call_worker(request, codes_root=None, timeout=300):
    """Run one worker request.  Returns the JSON object, or an unavailable/error dict."""
    try:
        proc = subprocess.run(
            [sys.executable, "-m", "moving.nbody_worker"],
            input=json.dumps(request),
            text=True,
            capture_output=True,
            timeout=timeout,
            cwd=_ROOT,
            env=_child_env(codes_root or request.get("codes_root")),
        )
    except subprocess.TimeoutExpired:
        return {"status": "error", "message": "n-body worker timed out after %s s" % timeout}
    text = (proc.stdout or "").strip()
    try:
        out = json.loads(text)
    except json.JSONDecodeError:
        err = (proc.stderr or "")[-800:]
        return {"status": "error", "message": "worker exit %s: %s" % (proc.returncode, err or text[-800:])}
    return out


def probe(codes_root=None, allow_assist=True, use_default=True):
    req = {"cmd": "probe", "allow_assist": allow_assist, "use_default": use_default}
    if codes_root is not None:
        req["codes_root"] = codes_root
    return call_worker(req, codes_root=codes_root, timeout=120)


def observe(state_helio, mjd_ref, mjd, code="500", codes_root=None, timeout=180,
            allow_assist=True, use_default=True):
    """Apparent RA/Dec (deg) of a heliocentric equatorial state (AU, AU/day).

    Light-time corrected, geocentre or an MPC station.  The default backend is
    ogfmeas.  ``allow_assist=False`` selects CODES when a tree with de440s.bsp
    is present.
    """
    mjd = np.asarray(mjd, float)
    codes = [code] * len(mjd) if isinstance(code, str) or code is None else list(code)
    req = {
        "cmd": "observe",
        "allow_assist": bool(allow_assist),
        "use_default": bool(use_default),
        "state_helio": np.asarray(state_helio, float).tolist(),
        "mjd_ref": float(mjd_ref),
        "mjd": mjd.tolist(),
        "code": codes,
    }
    if codes_root is not None:
        req["codes_root"] = codes_root
    out = call_worker(req, codes_root=codes_root, timeout=timeout)
    if out.get("status") != "ok":
        raise RuntimeError(out.get("message") or out.get("reason") or out)
    return np.asarray(out["ra"], float), np.asarray(out["dec"], float), out


def propagate_helio(state_helio, mjd_ref, mjd_out, codes_root=None, timeout=180,
                    allow_assist=True, use_default=True):
    """Heliocentric state at ``mjd_out``.  The default backend is ogfmeas."""
    req = {
        "cmd": "propagate_helio",
        "allow_assist": bool(allow_assist),
        "use_default": bool(use_default),
        "state_helio": np.asarray(state_helio, float).tolist(),
        "mjd_ref": float(mjd_ref),
        "mjd_out": float(mjd_out),
    }
    if codes_root is not None:
        req["codes_root"] = codes_root
    out = call_worker(req, codes_root=codes_root, timeout=timeout)
    if out.get("status") != "ok":
        raise RuntimeError(out.get("message") or out.get("reason") or out)
    return out


def _raw_rms(fit, members):
    mjd, ra, de, _sig, op = N._stack(members)
    pra, pde, _ = N.predict(fit.state[0], fit.state[1], fit.t_ref, mjd, op)
    dx, dy = N.sky_diff(ra, de, pra, pde)
    return float(np.sqrt(np.mean(dx ** 2 + dy ** 2) / 2.0))


def _job(group, members):
    fit = group["fit"]
    mjd = np.concatenate([t.mjd for t in members])
    ra = np.concatenate([t.ra for t in members])
    dec = np.concatenate([t.dec for t in members])
    sig = np.concatenate([t.sig for t in members])
    code = []
    for t in members:
        code.extend([t.code or "500"] * len(t))
    state = np.r_[np.asarray(fit.state[0], float), np.asarray(fit.state[1], float)]
    return {
        "id": list(group["ids"]),
        "mjd": mjd.tolist(),
        "ra": ra.tolist(),
        "dec": dec.tolist(),
        "sig": sig.tolist(),
        "code": code,
        "mjd_ref": float(fit.t_ref),
        "state_helio": state.tolist(),
    }


def _residual_ok(rep, chi2_max, resid_max):
    el = rep.get("elements") or {}
    return (rep.get("status") in ("converged", "not_converged")
            and rep.get("chi2_red", np.inf) <= chi2_max
            and rep.get("max_resid", np.inf) <= resid_max
            and rep.get("n_used", 0) >= 4
            and el.get("a", -1) > 0
            and el.get("e", 2) < 1)


def _apply(group, rep, chi2_max, resid_max):
    """Return 'accept', 'reject', or 'keep'.

    A fit that meets the chi2 and residual cuts replaces the two-body elements
    even if the last iteration missed the formal step tolerance.  Only a
    converged fit that fails those cuts removes the group.  An error, or a
    fit that never settled and still fails the cuts, leaves the two-body group.
    """
    group["two_body_rms_arcsec"] = group.get("rms_arcsec")
    group["two_body_elements"] = dict(group.get("elements") or {})
    group["nbody"] = {k: rep.get(k) for k in (
        "status", "chi2_red", "rms_arcsec", "max_resid", "n_used", "n_obs", "iterations", "message")}
    if _residual_ok(rep, chi2_max, resid_max):
        group["elements"] = rep["elements"]
        group["chi2_red"] = rep["chi2_red"]
        group["rms_arcsec"] = rep["rms_arcsec"]
        group["max_resid"] = rep["max_resid"]
        group["n_obs"] = rep["n_used"]
        group["nbody"]["accepted"] = True
        group["nbody_mjd"] = rep["mjd_epoch"]
        group["nbody_state_helio"] = rep["state_helio"]
        return "accept"
    group["nbody"]["accepted"] = False
    if rep.get("status") == "converged":
        return "reject"
    return "keep"


def refine_result(trks, result, codes_root=None, allow_assist=True, use_default=True,
                  chi2_max=4.0, resid_max=6.0, floor_arcsec=0.05, timeout=300):
    """Re-fit each confirmed group.  Mutates and returns ``result``."""
    by_id = {t.id: t for t in trks}
    jobs = []
    usable = []
    for group in result.get("groups", []):
        if group.get("fit") is None:
            continue
        members = [by_id[i] for i in group["ids"]]
        group["two_body_raw_rms_arcsec"] = _raw_rms(group["fit"], members)
        jobs.append(_job(group, members))
        usable.append(group)
    if not jobs:
        result["nbody"] = {"status": "empty"}
        return result
    req = {
        "cmd": "fit",
        "allow_assist": allow_assist,
        "use_default": use_default,
        "floor_arcsec": floor_arcsec,
        "groups": jobs,
    }
    if codes_root is not None:
        req["codes_root"] = codes_root
    reply = call_worker(req, codes_root=codes_root, timeout=timeout)
    if reply.get("status") != "ok":
        result["nbody"] = reply
        return result
    by_key = {tuple(int(x) for x in row.get("id", [])): row for row in reply.get("groups", [])}
    kept = []
    rejected = list(result.get("nbody_rejected", []))
    for group in result.get("groups", []):
        row = by_key.get(tuple(int(x) for x in group["ids"]))
        if row is None or group.get("fit") is None:
            kept.append(group)
            continue
        outcome = _apply(group, row, chi2_max, resid_max)
        if outcome == "reject":
            rejected.append(group)
        else:
            kept.append(group)
    result["groups"] = kept
    result["nbody_rejected"] = rejected
    # Same rule as link_nights: pair members stay out of the confirmed set, so they remain unlinked.
    linked = {i for g in kept for i in g["ids"]}
    result["unlinked"] = [t.id for t in trks if t.id not in linked]
    result["nbody"] = {
        "status": "ok",
        "backend": reply.get("backend"),
        "force_model": reply.get("force_model"),
        "codes_root": reply.get("codes_root"),
        "n_refined": sum(1 for g in kept if g.get("nbody", {}).get("accepted")),
        "n_rejected": len(rejected),
        "floor_arcsec": floor_arcsec,
    }
    return result
