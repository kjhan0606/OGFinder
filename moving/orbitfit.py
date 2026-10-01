"""High-level orbit determination from a set of astrometric observations (Obs container).

fit_orbit(obs):
  1. if no initial state: for arcs <= `short_arc_days` use statistical ranging (iod.ranging) to get a sample; the
     sample members with the best chi2 (top `n_starts`) seed ASSIST differential correction; the DC solution
     with the lowest chi2 is taken.  For longer arcs the same procedure is used with a larger admissible-region search.
  2. differential correction (orbit.differential_correction) with Carpino et al. (2003) outlier rejection.
  3. elements +- sigma at the epoch of the middle observation (or `epoch_jd`), orbit class, covariance.
  4. (optional) emcee posterior sampling of the 6 state parameters for short arcs where the covariance is non-Gaussian.
Photometric H is estimated separately (photometry.hg_absolute).
"""
import numpy as np
from . import orbit as O, iod, kepler as K, util
from .util import log


def _obs_helio_inputs(prop, obs):
    ob = obs.observer_bary(prop)
    sun = np.array([prop.sun(j) for j in obs.jd_tdb])
    return ob, sun


def seed_by_ranging(prop, obs, n_samples=600, **kw):
    ob, sun = _obs_helio_inputs(prop, obs)
    sig = 0.5 * (obs.sig_ra + obs.sig_dec)
    R = iod.ranging(obs.ra, obs.dec, obs.jd_tdb, sig, ob, sun, n_samples=n_samples, **kw)
    return R


def state_bary_from_helio(prop, state_helio_eq, jd):
    s = prop.sun(jd)
    return np.concatenate([state_helio_eq[:3] + s[:3], state_helio_eq[3:] + s[3:]])


def fit_orbit(obs, prop=None, state0=None, jd0=None, n_starts=6, ranging_samples=600, verbose=False, epoch_jd=None, **dc_kw):
    """jd0 = fit epoch (default: mean observation time, so IOD/DC never extrapolate); epoch_jd = epoch of the reported elements
    (default jd0; propagated with ASSIST, covariance mapped to first order)."""
    prop = prop or O.Propagator()
    jd0 = float(np.mean(obs.jd_tdb)) if jd0 is None else jd0
    info = {}
    if state0 is None:
        R = seed_by_ranging(prop, obs, n_samples=ranging_samples, jd_ref=jd0)
        info["ranging"] = dict(neff=R["neff"], class_probs=R["class_probs"], n_samples=R["n_samples"])
        # polish the best samples with a 6-parameter 2-body least-squares fit (needed for long/precise arcs)
        ob, sun = _obs_helio_inputs(prop, obs)
        sig = 0.5 * (obs.sig_ra + obs.sig_dec)
        idx = np.argsort(obs.jd_tdb)
        if len(idx) > 10:
            idx = idx[np.unique(np.round(np.linspace(0, len(idx) - 1, 10)).astype(int))]
        pol = iod.refine_two_body(R, obs.ra[idx], obs.dec[idx], obs.jd_tdb[idx], sig[idx], ob[idx], sun[idx], top=n_starts)
        info["two_body_chi2"] = [p_["chi2"] for p_ in pol]
        starts = [state_bary_from_helio(prop, p_["state"], R["jd_ref"]) for p_ in pol]
        info["ranging_result"] = R
    else:
        starts = [np.asarray(state0, float)]
    best = None
    for s in starts:
        try:
            res = O.differential_correction(prop, s, jd0, obs, verbose=verbose, **dict(dict(lm_lambda0=0.05, max_iter=25), **dc_kw))
        except Exception as e:
            log("DC failed:", e); continue
        if best is None or (res["chi2"] + 2 * (res["n_total"] - res["n_used"]) * 9) < (best["chi2"] + 2 * (best["n_total"] - best["n_used"]) * 9):
            best = res
    if best is None:
        raise RuntimeError("differential correction failed for all starting states")
    el, sg, Ce = O.elements_and_sigmas(prop, best["state"], best["cov"], jd0, jd_target=epoch_jd)
    best["epoch_jd"] = float(jd0 if epoch_jd is None else epoch_jd)
    best["elements"] = el; best["sigmas"] = sg; best["elem_cov"] = Ce
    best["orbit_class"] = K.orbit_class(el["a"], el["e"], el["i"])
    best["info"] = info
    # a unique orbit is only claimed when the fit is well conditioned; for arcs of a fraction of a day (single HST orbit)
    # the 6 parameters are degenerate and the DC "solution" is not meaningful (hyperbolic / enormous sigmas).
    arc = float(np.ptp(obs.mjd_utc))
    rel_a = sg["a"] / abs(el["a"]) if el["a"] != 0 else np.inf
    ok = (el["e"] < 1.0) and (el["a"] > 0) and (rel_a < 0.5) and (arc > 1.0 or best["ndof"] >= 6) and (sg["e"] < 0.5)
    best["determined"] = bool(ok)
    if not ok and "ranging_result" in info:
        R = info["ranging_result"]
        q = {k: iod.weighted_quantiles(R[k], R["w"]) for k in ("a", "e", "i", "q")}
        info["ranging_quantiles_16_50_84"] = q
        best["status_note"] = ("arc %.2f d too short for a unique orbit: elements below are NOT constrained; use the class probabilities and "
                               "the ranging 16/50/84%% quantiles (flat prior in range and range-rate over the admissible region, see docs)" % arc)
    return best


def mcmc_posterior(prop, obs, res, nwalkers=24, nsteps=400, burn=150, seed=1):
    """emcee posterior of the barycentric state (6 params) for a (possibly poorly constrained) short arc, started from the DC
    solution scattered by its covariance.  Likelihood: independent Gaussian residuals with the observation sigmas (rejected
    observations excluded)."""
    import emcee
    jd0 = res["jd0"]; active = res["active"]
    sub = obs.subset(active)
    sub._obspos = obs.observer_bary(prop)[active]
    sig = np.empty(2 * sub.n); sig[0::2] = sub.sig_ra; sig[1::2] = sub.sig_dec
    x0 = res["state"]; Cv = res["cov"]
    scale = np.sqrt(np.diag(Cv))

    def lnp(theta):
        st = x0 + theta * scale
        try:
            ra, dec, _, _ = O.predict_radec(prop, st[None, :], jd0, sub)[:4]
        except Exception:
            return -np.inf
        d1, d2 = O.residuals_arcsec(sub.ra, sub.dec, ra[0], dec[0])
        r = np.empty(2 * sub.n); r[0::2] = d1; r[1::2] = d2
        return -0.5 * np.sum((r / sig) ** 2)
    rng = np.random.default_rng(seed)
    p0 = rng.normal(0, 0.3, (nwalkers, 6))
    S = emcee.EnsembleSampler(nwalkers, 6, lnp)
    S.run_mcmc(p0, nsteps, progress=False)
    ch = S.get_chain(discard=burn, flat=True)
    states = x0 + ch * scale
    return dict(states=states, scale=scale, acceptance=float(np.mean(S.acceptance_fraction)),
                autocorr=(S.get_autocorr_time(tol=0) if hasattr(S, "get_autocorr_time") else None))
