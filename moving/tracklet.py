"""Tracklets: linking of detections between exposures, constant-velocity fits with uncertainties,
and a simplified HelioLinC-style heliocentric linking for longer arcs.

* link_exposures(): for exposures with times t_i, detections (ra, dec, sigma, flux).  All pairs (i<j) of
  detections in different exposures define a velocity; a third+ exposure confirms if a detection lies within
  `tol` of the prediction.  Greedy best-chi2 assignment, each detection used once.  Observer motion
  (HST orbital parallax) is NOT removed from the linear fit by default: the fit is done on the apparent
  topocentric sky motion; rate and PA are reported topocentric.
* fit_tracklet(): weighted least squares of (xi, eta) tangent-plane offsets linear in time -> rate, PA,
  covariance (position at reference time and rates), chi2.
* helio_linc(): HelioLinC (Holman et al. 2018; Heinze et al. 2022), simplified: tracklets -> (r, rdot)
  hypotheses grid -> propagate each tracklet attributable (ra, dec, rates + observer vector) to a heliocentric
  position at a common epoch with 2-body (universal variable) motion -> cluster in 6D with scaled DBSCAN
  (instead of the original KD-tree with velocity tolerance).  Validated only on synthetic data.
"""
import numpy as np
from .util import tangent_offsets, C_AUD, angsep_arcsec
from . import kepler as K


def fit_tracklet(t_days, ra, dec, sig_arcsec, t_ref=None, obs_off_au=None, fit_parallax=False, inv_delta=0.0):
    """Constant angular velocity fit of detections.  t in days (MJD), ra/dec deg, sig arcsec (1-sigma per axis).

    obs_off_au: (n,3) observer displacement from the geocentre (AU, equatorial) -- needed for spacecraft such as
    HST whose orbital parallax (up to ~10 arcsec for a main-belt object over one orbit) dominates the apparent motion.
    The model is  xi(t) = x0 + vx (t-t_ref) - 206265 * d_perp_x(t) * k,  k = 1/Delta [1/AU]; k is fixed to `inv_delta`
    or fitted if fit_parallax=True (needs >= 4 exposures to leave 1 dof, >=3 with a prior on k not implemented).
    Returns dict(ra0, dec0, t_ref, vx_asd, vy_asd [arcsec/day, geocentric motion when k is applied], rate_ash, pa_deg,
    sig_rate_ash, sig_pa_deg, inv_delta, chi2, ndof, rms_resid_arcsec, ...)."""
    t = np.asarray(t_days, float); n = len(t)
    t_ref = float(np.mean(t)) if t_ref is None else t_ref
    ra0 = float(np.mean(ra)); de0 = float(np.mean(dec))
    xi, eta = tangent_offsets(np.asarray(ra), np.asarray(dec), ra0, de0)
    sig = np.broadcast_to(np.asarray(sig_arcsec, float), (n,))
    cols_x = [np.ones(n), t - t_ref]
    # tangent-plane components of the observer offset
    if obs_off_au is not None:
        a0, d0 = np.radians(ra0), np.radians(de0)
        ex = np.array([-np.sin(a0), np.cos(a0), 0.0])
        ey = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
        ox = np.asarray(obs_off_au) @ ex * 206264.806; oy = np.asarray(obs_off_au) @ ey * 206264.806      # arcsec * AU
    else:
        ox = oy = np.zeros(n)
    W = 1.0 / sig ** 2
    if fit_parallax and obs_off_au is not None and n >= 4:
        # unknown: x0, y0, vx, vy, k ; observations stacked (xi then eta)
        A = np.zeros((2 * n, 5)); y = np.r_[xi, eta]; w = np.r_[W, W]
        A[:n, 0] = 1; A[:n, 2] = t - t_ref; A[:n, 4] = -ox
        A[n:, 1] = 1; A[n:, 3] = t - t_ref; A[n:, 4] = -oy
        N = A.T @ (A * w[:, None]); Ci5 = np.linalg.inv(N)
        p = Ci5 @ (A.T @ (w * y)); r = y - A @ p
        cx = np.array([p[0], p[2]]); cy = np.array([p[1], p[3]]); k = float(p[4])
        npar = 5
        Ci = np.array([[Ci5[0, 0], Ci5[0, 2]], [Ci5[2, 0], Ci5[2, 2]]])
        rx, ry = r[:n], r[n:]
        sig_k = float(np.sqrt(Ci5[4, 4]))
    else:
        k = float(inv_delta)
        xi_c = xi + k * ox; eta_c = eta + k * oy               # remove parallax shift for the hypothesised k
        A = np.stack([np.ones(n), t - t_ref], 1)
        N = A.T @ (A * W[:, None]); Ci = np.linalg.inv(N)
        cx = Ci @ (A.T @ (W * xi_c)); cy = Ci @ (A.T @ (W * eta_c))
        rx = xi_c - A @ cx; ry = eta_c - A @ cy
        npar = 4; sig_k = np.nan
    chi2 = float(np.sum(W * (rx ** 2 + ry ** 2)))
    ndof = 2 * n - npar
    vx, vy = cx[1], cy[1]
    rate = np.hypot(vx, vy)
    pa = np.degrees(np.arctan2(vx, vy)) % 360.0
    sv = np.sqrt(Ci[1, 1])
    from .align import tangent_inverse
    ra_r, de_r = tangent_inverse(cx[0], cy[0], ra0, de0)
    return dict(ra0=float(ra_r), dec0=float(de_r), t_ref=t_ref, vx_asd=float(vx), vy_asd=float(vy),
                rate_ash=float(rate / 24.0), sig_rate_ash=float(sv / 24.0), pa_deg=float(pa),
                sig_pa_deg=float(np.degrees(sv / rate)) if rate > 0 else np.nan, inv_delta=k, sig_inv_delta=sig_k,
                chi2=chi2, ndof=int(ndof), chi2_red=chi2 / max(ndof, 1),
                rms_resid_arcsec=float(np.sqrt(np.mean(rx ** 2 + ry ** 2))), n=n,
                sig_pos0_arcsec=float(np.sqrt(Ci[0, 0])), cov_v=float(Ci[1, 1]))


def _trail_consistency(f, m):
    """Compare the intra-exposure trail (orientation, length) of each member with the tracklet motion: a real mover
    leaves a trail of length ~ rate * t_exp along its direction of motion.  Adds f['trail_pa_dev_deg'] (rms of the
    orientation difference, mod 180), f['trail_len_ratio'] (median observed/predicted length) and f['score']
    (lower = better; combines positional rms in units of the expected noise and the trail mismatch)."""
    devs = []; ratios = []
    for x in m:
        if x.get("trail_pa") is None or not x.get("texp"):
            continue
        dpa = abs(((x["trail_pa"] - f["pa_deg"] + 90) % 180) - 90)
        devs.append(dpa)
        pred = f["rate_ash"] * x["texp"] / 3600.0
        if pred > 0 and x.get("trail_len"):
            ratios.append(x["trail_len"] / pred)
    f["trail_pa_dev_deg"] = float(np.sqrt(np.mean(np.square(devs)))) if devs else np.nan
    f["trail_len_ratio"] = float(np.median(ratios)) if ratios else np.nan
    pen = 0.0
    if devs:
        pen += (f["trail_pa_dev_deg"] / 20.0) ** 2
    if ratios:
        # measured extents (4 sigma of the moment ellipse) overestimate the geometric trail length by ~1.5-2.5x (PSF + seeing),
        # so only ratios outside [0.5, 3.5] are penalised
        r = f["trail_len_ratio"]
        if r > 3.5:
            pen += (np.log(r / 3.5) / 0.4) ** 2
        elif r < 0.5:
            pen += (np.log(0.5 / max(r, 1e-3)) / 0.4) ** 2
    # physical consistency: a mover with rate v trails by v * t_exp during an exposure.  A member detected as a compact
    # POINT although the tracklet predicts a trail much longer than the PSF contradicts the hypothesis.
    npt = 0
    for x in m:
        if x.get("texp") and x.get("channel", "point") != "trail":
            pred = f["rate_ash"] * x["texp"] / 3600.0
            if pred > x.get("psf_fwhm_arcsec", 0.09) * 4.0:
                npt += 1
    f["n_point_contradict"] = npt
    pen += 4.0 * npt
    # brightness consistency: compared within the same detection channel only (trail flux is a total, point flux a peak value)
    sc = []
    for ch in {x.get("channel", "") for x in m}:
        fl = np.array([x.get("flux", np.nan) for x in m if x.get("channel", "") == ch], float)
        fl = fl[np.isfinite(fl) & (fl > 0)]
        if len(fl) >= 3:
            sc.append(np.std(np.log10(fl)))
    f["flux_scatter_dex"] = float(max(sc)) if sc else np.nan
    if sc:                                        # a mover keeps its brightness
        pen += (f["flux_scatter_dex"] / 0.12) ** 2
    k = f.get("inv_delta", 0.0)
    if np.isfinite(k) and (k < 0.0 or k > 3.5):          # 1/Delta [1/AU] outside any plausible range for a Solar-System body
        pen += 25.0
    f["score"] = float(f["rms_resid_arcsec"] / 0.1 + pen)


def _spread_triplets(nex):
    """Exposure triplets (indices into the time-sorted exposure list) used to seed tracklets: all C(n,3) for n <= 6, else
    (first, j, last), (first, j, last-1) and (1, j, last) -- every real track of >= 3 detections contains one of them."""
    import itertools
    if nex <= 6:
        return list(itertools.combinations(range(nex), 3))
    L = nex - 1; out = set()
    for j in range(1, L):
        out.add((0, j, L))
        if j != L - 1:
            out.add((0, j, L - 1))
        if j != 1:
            out.add((1, j, L))
    return sorted(out)


def link_exposures(dets, tol_arcsec=0.6, min_exposures=3, rate_min_ash=0.5, rate_max_ash=600.0,
                   sig_floor=0.01, max_tracklets=500, obs_off_au=None, inv_delta_grid=None,
                   max_cand_per_exposure=3000, min_disp_arcsec=1.5, flux_tol_dex=None,
                   k_max=3.5, rate_per_k_ash=300.0, rate_slack_ash=4.0, ext_tol_factor=2.0, rms_max_arcsec=None,
                   sigma_eff_arcsec=0.35, cr_penalty=3.0, missing_penalty=8.0, k_prior_weight=1.0, k_prior_center=0.3,
                   rate_prior_weight=1.0, rate_prior_start_ash=30.0, max_seeds=4_000_000, dedupe=True, stats=None):
    """Link detections of >= 3 exposures into tracklets (constant geocentric velocity + orbital parallax of the observer).

    dets: list of dict(ex=int exposure index, t=MJD, ra, dec [deg], sig [arcsec], id, flux; optional cls, channel ('point'|'trail'),
    trail_pa, trail_len [arcsec], texp [s]).  obs_off_au: dict ex -> geocentric observer offset (AU, equatorial; HST); without
    it the observer is taken as fixed (k = 0).

    Algorithm (pairwise-then-extend, velocity-space consistency; replaces the first/last-pair linker kept as
    `link_exposures_legacy`):
      1. seeds: for every exposure triplet (a, b, c) all detection pairs (a, c) with a displacement >= min_disp; the middle
         exposure b is searched by a sorted 1-D lookup *along the parallax direction*: the observer offset makes the apparent
         position of b move linearly with k = 1/Delta [1/AU], so the b detection fixes k analytically (no k grid), and only
         the perpendicular distance has to be < tol.  Seeds are gated physically: 0 <= k <= k_max and the geocentric rate
         <= rate_per_k * k + slack (an object bound to the Sun moves at most ~ 60 km/s / Delta = 300 arcsec/h * k).
      2. extension: the seed's constant-velocity+parallax prediction is looked up in every other exposure (KD-tree, radius
         ext_tol_factor * tol); identical member sets are merged.
      3. every candidate gets a weighted least-squares fit of (x0, y0, vx, vy, k) (k fixed to the seed value if it has only 3
         members: 3 points are always fit by a free k, so they cannot be validated), a gate on the rms residual, and a
         score (lower = better, chi2-like): positional chi2 + class penalty (members the single-exposure classifier flagged as
         cosmic ray / edge / faint) + trail consistency (orientation, length) + brightness consistency within a detection
         channel + `missing_penalty` per exposure without a member.
      4. candidates are ranked by (score) with the missing-exposure penalty included, de-duplicated (a candidate sharing >= 2
         detections with a better one is the same object seen through the other channel / with a wrong member and is
         dropped; subsets are dropped) and the best `max_tracklets` are returned with the final fit_tracklet() solution.
    `stats` (dict) receives counters (seeds, candidates, gated, kept).  Returns list of tracklet dicts as before."""
    if not dets:
        return []
    min_exposures = max(int(min_exposures), 3)
    byex = {}
    for d in dets:
        byex.setdefault(d["ex"], []).append(d)
    exs = sorted(byex, key=lambda e: byex[e][0]["t"])
    nex = len(exs)
    if nex < min_exposures:
        return []
    for e in exs:                                                   # bound the work per exposure: brightest first
        if len(byex[e]) > max_cand_per_exposure:
            byex[e] = sorted(byex[e], key=lambda d: -(d.get("snr") or 0))[:max_cand_per_exposure]
    ra_c = float(np.mean([d["ra"] for d in dets])); de_c = float(np.mean([d["dec"] for d in dets]))
    a0, d0 = np.radians(ra_c), np.radians(de_c)
    ex_v = np.array([-np.sin(a0), np.cos(a0), 0.0]); ey_v = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
    X = []; Tt = []; Ob = []
    for e in exs:
        x, y = tangent_offsets(np.array([d["ra"] for d in byex[e]]), np.array([d["dec"] for d in byex[e]]), ra_c, de_c)
        X.append(np.c_[x, y]); Tt.append(float(byex[e][0]["t"]))
        o = np.asarray(obs_off_au[e], float) * 206264.806 if obs_off_au is not None else np.zeros(3)
        Ob.append(np.array([o @ ex_v, o @ ey_v]))
    Tt = np.array(Tt); Ob = np.array(Ob)
    from scipy.spatial import cKDTree
    trees = [cKDTree(x) for x in X]

    rows = []; K0 = []; nseed = 0
    for (ia_, ib_, ic_) in _spread_triplets(nex):
        Xa, Xb, Xc = X[ia_], X[ib_], X[ic_]
        if len(Xa) == 0 or len(Xb) == 0 or len(Xc) == 0:
            continue
        ta, tb, tc = Tt[ia_], Tt[ib_], Tt[ic_]
        if tc <= ta:
            continue
        g = (tb - ta) / (tc - ta)
        D = (Ob[ia_] - Ob[ib_]) + g * (Ob[ic_] - Ob[ia_])
        nD = float(np.hypot(*D))
        raw = np.hypot(Xc[None, :, 0] - Xa[:, None, 0], Xc[None, :, 1] - Xa[:, None, 1])
        ja, jc = np.nonzero(raw >= min_disp_arcsec)
        if len(ja) == 0:
            continue
        S0 = Xa[ja] + g * (Xc[jc] - Xa[ja])
        if nD > 0.3:
            Dh = D / nD; nh = np.array([-Dh[1], Dh[0]])
            pb = Xb @ nh; order = np.argsort(pb); pbs = pb[order]
            tp = S0 @ nh
            lo = np.searchsorted(pbs, tp - tol_arcsec); hi = np.searchsorted(pbs, tp + tol_arcsec)
            cnt = hi - lo; m = cnt > 0
            ja, jc, lo, cnt, S0 = ja[m], jc[m], lo[m], cnt[m], S0[m]
            if cnt.sum() == 0:
                continue
            rep = np.repeat(np.arange(len(ja)), cnt)
            off = np.arange(cnt.sum()) - np.repeat(np.cumsum(cnt) - cnt, cnt)
            jb = order[lo[rep] + off]
            R = Xb[jb] - S0[rep]
            k = (R @ Dh) / nD
            kslack = tol_arcsec / nD            # the middle detection fixes k only to +-tol/|D| (small HST baseline: |D| ~ 0.1-0.3")
            ok = (k >= -0.05 - kslack) & (k <= k_max + kslack) & (np.abs(R @ nh) < tol_arcsec)
            ja, jb, jc, k = ja[rep][ok], jb[ok], jc[rep][ok], np.clip(k[ok], 0.0, k_max)
        else:
            # (nearly) linear observer motion over the triplet: k is not constrained by the middle exposure alone (parallax is
            # degenerate with the velocity), so scan a short k grid (only taken for a fixed observer or a straight-line orbit arc)
            kgrid = np.array([0.0]) if obs_off_au is None else np.arange(0.0, k_max + 1e-9, 0.1)
            JA = []; JB = []; JC = []; KK = []
            for kg in kgrid:
                pa_k = Xa[ja] + kg * Ob[ia_]; pc_k = Xc[jc] + kg * Ob[ic_]
                P = pa_k + g * (pc_k - pa_k) - kg * Ob[ib_]
                r = trees[ib_].query_ball_point(P, tol_arcsec)
                cnt = np.array([len(x) for x in r])
                if cnt.sum() == 0:
                    continue
                rep = np.repeat(np.arange(len(ja)), cnt)
                JB.append(np.concatenate([np.array(x, int) for x in r if len(x)])); JA.append(ja[rep]); JC.append(jc[rep]); KK.append(np.full(len(rep), kg))
            if not JA:
                continue
            ja, jb, jc, k = np.concatenate(JA), np.concatenate(JB), np.concatenate(JC), np.concatenate(KK)
        if len(ja) == 0:
            continue
        # physical gate on the geocentric rate implied by (a, c) and k
        va = (Xc[jc] + k[:, None] * Ob[ic_]) - (Xa[ja] + k[:, None] * Ob[ia_])
        rate = np.hypot(va[:, 0], va[:, 1]) / ((tc - ta) * 24.0)
        okr = (rate >= rate_min_ash) & (rate <= min(rate_max_ash, 1e9)) & (rate <= rate_per_k_ash * k + rate_slack_ash)
        ja, jb, jc, k = ja[okr], jb[okr], jc[okr], k[okr]
        nseed += len(ja)
        if nseed > max_seeds:
            from .util import log
            log("link_exposures: seed limit %d reached (raise max_seeds or tighten snr / tol)" % max_seeds)
            break
        if len(ja) == 0:
            continue
        mem = -np.ones((len(ja), nex), int)
        mem[:, ia_] = ja; mem[:, ib_] = jb; mem[:, ic_] = jc
        # extension to the other exposures with the constant-velocity + parallax prediction
        pa_ = Xa[ja] + k[:, None] * Ob[ia_]; pc_ = Xc[jc] + k[:, None] * Ob[ic_]
        vel = (pc_ - pa_) / (tc - ta)
        for q in range(nex):
            if q in (ia_, ib_, ic_) or len(X[q]) == 0:
                continue
            P = pa_ + vel * (Tt[q] - ta) - k[:, None] * Ob[q]
            dd, jj = trees[q].query(P, distance_upper_bound=tol_arcsec * ext_tol_factor)
            okq = np.isfinite(dd)
            mem[okq, q] = jj[okq]
        rows.append(mem); K0.append(k)
    if not rows:
        return []
    M = np.vstack(rows); K0 = np.concatenate(K0)
    M, u = np.unique(M, axis=0, return_index=True); K0 = K0[u]
    nm = (M >= 0).sum(1)
    keep = nm >= min_exposures
    M, K0, nm = M[keep], K0[keep], nm[keep]
    n = len(M)
    if stats is not None:
        stats.update(seeds=nseed, candidates=n)
    if n == 0:
        return []
    has = M >= 0

    def gather(fn, default=np.nan):
        out = np.full((n, nex), default, float)
        for i, e in enumerate(exs):
            a = np.array([fn(d) for d in byex[e]], float)
            h = has[:, i]; out[h, i] = a[M[h, i]]
        return out
    px = np.zeros((n, nex)); py = np.zeros((n, nex))
    for i in range(nex):
        h = has[:, i]; px[h, i] = X[i][M[h, i], 0]; py[h, i] = X[i][M[h, i], 1]
    tref = float(Tt.mean()); dt = (Tt - tref)[None, :] * np.ones((n, 1))
    sg = gather(lambda d: max(d.get("sig", 0.05) or 0.05, sig_floor), 1.0)
    sg = np.sqrt(sg ** 2 + sigma_eff_arcsec ** 2)
    W = has / sg ** 2
    A = np.zeros((n, nex, 2, 5))
    A[:, :, 0, 0] = 1; A[:, :, 1, 1] = 1; A[:, :, 0, 2] = dt; A[:, :, 1, 3] = dt
    A[:, :, 0, 4] = -Ob[None, :, 0]; A[:, :, 1, 4] = -Ob[None, :, 1]
    y = np.stack([px, py], 2); w = W[:, :, None] * np.ones((1, 1, 2))
    N = np.einsum("nepi,nep,nepj->nij", A, w, A)
    rhs = np.einsum("nepi,nep,nep->ni", A, w, y)
    free_k = (nm >= 4) & (obs_off_au is not None)
    prior = np.where(free_k, 1e-6, 1e9)
    N[:, 4, 4] += prior; rhs[:, 4] += prior * K0 * (~free_k)
    if obs_off_au is None:                                          # no parallax column at all
        N[:, 4, 4] = 1.0; rhs[:, 4] = 0.0
    sol = np.linalg.solve(N, rhs[:, :, None])[:, :, 0]
    res = (y - np.einsum("nepi,ni->nep", A, sol)) * has[:, :, None]
    chi2 = (res ** 2 * w).sum((1, 2))
    rms = np.sqrt((res ** 2).sum((1, 2)) / (2.0 * nm))
    kfit = sol[:, 4]; vx, vy = sol[:, 2], sol[:, 3]
    # formal uncertainty of the free k: with the small HST baseline (~0.1-0.4 arcsec of parallax per unit k) k is often only
    # known to +-1, so the physical gate 0 <= k <= k_max is applied to k within 2 sigma, and the score uses k clipped into range
    with np.errstate(all="ignore"):
        sk = np.sqrt(np.maximum(np.linalg.inv(N)[:, 4, 4], 0.0))
    sk = np.where(free_k, np.nan_to_num(sk, nan=0.0, posinf=0.0), 0.0)
    kfit_raw = kfit.copy(); kfit = np.clip(kfit, 0.0, k_max)
    rate = np.hypot(vx, vy) / 24.0
    rms_max = rms_max_arcsec if rms_max_arcsec is not None else max(0.5 * tol_arcsec, 0.3)
    gate = (rms <= rms_max) & (kfit_raw >= -0.05 - 2 * sk) & (kfit_raw <= k_max + 2 * sk) & (rate >= rate_min_ash * 0.5) & (rate <= rate_per_k_ash * np.clip(kfit_raw + 2 * sk, 0.0, k_max) + 1.5 * rate_slack_ash + 1e-9) \
        & (rate <= rate_max_ash)
    if stats is not None:
        stats["gated"] = int(gate.sum())
    if not gate.any():
        return []
    # ---- score (lower = better)
    TR = gather(lambda d: 1.0 if d.get("channel", "point") == "trail" else 0.0, 0.0)
    FL = gather(lambda d: d.get("flux", np.nan) or np.nan)
    PA = gather(lambda d: d["trail_pa"] if d.get("trail_pa") is not None else np.nan)
    LN = gather(lambda d: d["trail_len"] if d.get("trail_len") is not None else np.nan)
    TX = gather(lambda d: d.get("texp") or 0.0, 0.0)
    BAD = gather(lambda d: 1.0 if d.get("cls") in ("artefact_cr", "artefact_edge", "artefact_dipole", "faint") else 0.0, 0.0)
    score = chi2.copy()
    score += cr_penalty * (BAD * has).sum(1)
    score += missing_penalty * (nex - nm)
    # weak population priors (tuned on the HST injection tests; set the weights to 0 for NEO/close-approach searches):
    # log-normal prior on k = 1/Delta around k_prior_center (main belt, Delta ~ 3 AU) and a penalty on geocentric rates
    # above rate_prior_start_ash (most small bodies in a deep field are slow)
    score += k_prior_weight * np.log(np.maximum(kfit, 0.05) / k_prior_center) ** 2
    score += rate_prior_weight * np.maximum(rate - rate_prior_start_ash, 0.0) / 10.0
    # apparent (topocentric) velocity at every exposure: v - k dO/dt
    if nex >= 2:
        dO = np.gradient(Ob, Tt, axis=0) if nex > 2 else np.tile((Ob[1] - Ob[0]) / (Tt[1] - Tt[0]), (2, 1))
    else:
        dO = np.zeros_like(Ob)
    tvx = vx[:, None] - kfit[:, None] * dO[None, :, 0]; tvy = vy[:, None] - kfit[:, None] * dO[None, :, 1]
    pa_e = np.degrees(np.arctan2(tvx, tvy)) % 360.0
    dpa = np.abs(((PA - pa_e + 90.0) % 180.0) - 90.0)
    pred = np.hypot(tvx, tvy) / 24.0 * TX / 3600.0
    with np.errstate(invalid="ignore", divide="ignore"):
        lr = np.log(np.clip(LN / np.maximum(pred, 1e-3), 1e-3, 1e3))
    ist = (TR > 0) & has
    score += np.where(ist & np.isfinite(dpa), np.minimum((dpa / 25.0) ** 2, 9.0), 0.0).sum(1)
    score += np.where(ist & np.isfinite(lr), np.minimum((lr / 0.7) ** 2, 9.0), 0.0).sum(1)
    for chv in (0.0, 1.0):
        m = (TR == chv) & has & (FL > 0)
        L = np.where(m, np.log10(np.where(FL > 0, FL, 1.0)), np.nan); c = m.sum(1)
        mu = np.nansum(L, 1) / np.maximum(c, 1)
        sd = np.sqrt(np.nansum((L - mu[:, None]) ** 2, 1) / np.maximum(c - 1, 1))
        score += np.where(c >= 2, (sd / 0.25) ** 2 * (c - 1), 0.0)
    score = np.where(gate, score, np.inf)

    order = np.argsort(score, kind="stable")
    order = order[np.isfinite(score[order])]
    # ---- de-duplicate: pairs of member detections already claimed by a better tracklet
    kept = []; claimed = set()
    for q in order:
        mem = [(i, int(M[q, i])) for i in range(nex) if M[q, i] >= 0]
        pairs = [(mem[a_], mem[b_]) for a_ in range(len(mem)) for b_ in range(a_ + 1, len(mem))]
        if dedupe and any(p in claimed for p in pairs):
            continue
        claimed.update(pairs)
        kept.append(q)
        if len(kept) >= max_tracklets:
            break
    if stats is not None:
        stats["kept"] = len(kept)
    out = []
    for q in kept:
        mem = [(i, int(M[q, i])) for i in range(nex) if M[q, i] >= 0]
        m = [byex[exs[i]][j] for i, j in mem]
        t = np.array([x["t"] for x in m]); ra = np.array([x["ra"] for x in m]); de = np.array([x["dec"] for x in m])
        sgm = np.array([max(x.get("sig", 0.05), sig_floor) for x in m])
        off = np.array([obs_off_au[x["ex"]] for x in m]) if obs_off_au is not None else None
        f = fit_tracklet(t, ra, de, sgm, obs_off_au=off, fit_parallax=(off is not None and len(m) >= 4),
                         inv_delta=float(kfit[q]) if len(m) < 4 else float(K0[q]))
        f["members"] = [x["id"] for x in m]; f["ex"] = [x["ex"] for x in m]; f["t"] = t.tolist()
        f["ra"] = ra.tolist(); f["dec"] = de.tolist(); f["sig"] = sgm.tolist()
        f["flux"] = [x.get("flux", np.nan) for x in m]
        _trail_consistency(f, m)
        f["link_chi2"] = float(chi2[q]); f["link_rms_arcsec"] = float(rms[q]); f["n_missing"] = int(nex - nm[q])
        f["score_legacy"] = f["score"]
        f["score"] = float(score[q])
        out.append(f)
    return out


def link_exposures_legacy(dets, tol_arcsec=0.6, min_exposures=3, rate_min_ash=0.5, rate_max_ash=600.0,
                   sig_floor=0.01, max_tracklets=500, obs_off_au=None, inv_delta_grid=tuple(np.r_[0.0, np.arange(0.08, 1.2, 0.04), np.arange(1.3, 3.01, 0.2)]),
                   max_cand_per_exposure=3000, min_disp_arcsec=1.5, flux_tol_dex=0.35):
    """[Legacy linker (pre-2026-10): per-parallax-hypothesis first/last pairing, kept for the before/after benchmark in
    moving/validation/link_bench.py and for regression comparison.]
    dets: list of dict(ex=int exposure index, t=MJD, ra, dec, sig [arcsec], id, flux).
    obs_off_au: dict ex -> geocentric observer offset (3-vector, AU); with it, for each hypothesised 1/Delta in
    `inv_delta_grid` the detections are first corrected for orbital parallax, then linked with a constant-velocity
    model (this is what makes HST arcs linkable).  Returns list of tracklets (dicts with members, fit)."""
    if not dets:
        return []
    byex = {}
    for d in dets:
        byex.setdefault(d["ex"], []).append(d)
    exs = sorted(byex, key=lambda e: byex[e][0]["t"])
    if len(exs) < min_exposures:
        return []
    ra_c = np.mean([d["ra"] for d in dets]); de_c = np.mean([d["dec"] for d in dets])
    a0, d0 = np.radians(ra_c), np.radians(de_c)
    ex_v = np.array([-np.sin(a0), np.cos(a0), 0.0]); ey_v = np.array([-np.sin(d0) * np.cos(a0), -np.sin(d0) * np.sin(a0), np.cos(d0)])
    from scipy.spatial import cKDTree
    grids = list(inv_delta_grid) if obs_off_au is not None else [0.0]
    cand = {}
    first, last = exs[0], exs[-1]
    t_first = byex[first][0]["t"]; t_last = byex[last][0]["t"]
    dt_h = (t_last - t_first) * 24.0
    if dt_h <= 0:
        return []
    base = {e: tangent_offsets(np.array([d["ra"] for d in byex[e]]), np.array([d["dec"] for d in byex[e]]), ra_c, de_c) for e in exs}
    for k in grids:
        P = {}
        for e in exs:
            xs, ys = base[e]
            if obs_off_au is not None and k:
                o = np.asarray(obs_off_au[e]) * 206264.806
                xs = xs + k * (o @ ex_v); ys = ys + k * (o @ ey_v)
            P[e] = (xs, ys)
        # raw (uncorrected) displacement test -> not a static source
        DX = P[last][0][None, :] - P[first][0][:, None]; DY = P[last][1][None, :] - P[first][1][:, None]
        RAWX = base[last][0][None, :] - base[first][0][:, None]; RAWY = base[last][1][None, :] - base[first][1][:, None]
        rate = np.hypot(DX, DY) / dt_h
        okp = (rate >= rate_min_ash) & (rate <= rate_max_ash) & (np.hypot(RAWX, RAWY) >= min_disp_arcsec)
        if flux_tol_dex is not None:
            # brightness gate between first and last detection (same detection channel only: trail flux is a total, point flux a peak)
            F0 = np.array([max(d.get("flux", 0) or 0, 1e-6) for d in byex[first]]); F1 = np.array([max(d.get("flux", 0) or 0, 1e-6) for d in byex[last]])
            C0 = np.array([d.get("channel", "") for d in byex[first]]); C1 = np.array([d.get("channel", "") for d in byex[last]])
            same = C0[:, None] == C1[None, :]
            dl = np.abs(np.log10(F0)[:, None] - np.log10(F1)[None, :])
            okp &= (~same) | (dl < flux_tol_dex)
        ia, ib = np.nonzero(okp)
        if len(ia) == 0:
            continue
        counts = np.full(len(ia), 2, int)
        nn = {}                                           # exposure -> index of the matched intermediate detection (-1 none)
        for e in exs[1:-1]:
            tt = byex[e][0]["t"]; f = (tt - t_first) / (t_last - t_first)
            px = P[first][0][ia] + f * DX[ia, ib]; py = P[first][1][ia] + f * DY[ia, ib]
            tree = cKDTree(np.c_[P[e][0], P[e][1]])
            dd, jj = tree.query(np.c_[px, py], distance_upper_bound=tol_arcsec)
            ok_ = np.isfinite(dd)
            counts += ok_
            nn[e] = np.where(ok_, jj, -1)
        for q in np.nonzero(counts >= min_exposures)[0]:
            pk = [(first, int(ia[q])), (last, int(ib[q]))] + [(e, int(nn[e][q])) for e in nn if nn[e][q] >= 0]
            key = tuple(sorted(pk))
            if key not in cand:
                cand[key] = int(counts[q]); cand[("k", key)] = k
    # fit every candidate, then take them greedily by (number of exposures, chi2_red)
    fits = []
    for key in (kk for kk in cand if kk[0] != "k"):
        mem = list(key)
        m = [byex[e][i] for e, i in mem]
        t = np.array([x["t"] for x in m]); ra = np.array([x["ra"] for x in m]); de = np.array([x["dec"] for x in m])
        sg = np.array([max(x.get("sig", 0.05), sig_floor) for x in m])
        off = np.array([obs_off_au[x["ex"]] for x in m]) if obs_off_au is not None else None
        f = fit_tracklet(t, ra, de, sg, obs_off_au=off, fit_parallax=(off is not None and len(m) >= 4),
                         inv_delta=cand.get(("k", key), 0.0))
        f["members"] = [x["id"] for x in m]; f["ex"] = [x["ex"] for x in m]; f["t"] = t.tolist()
        f["ra"] = ra.tolist(); f["dec"] = de.tolist(); f["sig"] = sg.tolist()
        f["flux"] = [x.get("flux", np.nan) for x in m]
        f["_key"] = mem
        _trail_consistency(f, m)
        fits.append(f)
    fits.sort(key=lambda f: (-f["n"], f["score"]))
    # Non-exclusive selection: a detection may belong to several competing tracklets (greedy exclusive assignment lets
    # chance alignments with a smaller residual steal members from the true track in crowded, CR-rich data).  Only
    # tracklets whose members are a subset of an already kept, better tracklet are dropped.  Ranking is by `score`.
    fits.sort(key=lambda f: (-f["n"], f["score"]))        # more exposures first, then consistency score
    tracklets = []
    for f in fits:
        key = set(f["_key"])
        if any(key <= set(t["_members_key"]) for t in tracklets):
            continue
        f["_members_key"] = f.pop("_key")
        tracklets.append(f)
        if len(tracklets) >= max_tracklets:
            break
    for t in tracklets:
        t.pop("_members_key", None)
    return tracklets


def link_trail_pairs(dets, **kw):
    return link_exposures(dets, **kw)


# --------------------------------------------------------------------------- HelioLinC (simplified)
def attributable_state(ra, dec, mu_ra_asd, mu_dec_asd, obs_pos, obs_vel, rho, rho_dot):
    """Heliocentric-equatorial... returns barycentric-equatorial state (AU, AU/d) of the object whose topocentric
    line of sight is (ra,dec) with angular rates mu (arcsec/day) at range rho [AU] and range-rate rho_dot [AU/d]."""
    a = np.radians(ra); d = np.radians(dec)
    los = np.array([np.cos(d) * np.cos(a), np.cos(d) * np.sin(a), np.sin(d)])
    e_ra = np.array([-np.sin(a), np.cos(a), 0.0])
    e_de = np.array([-np.sin(d) * np.cos(a), -np.sin(d) * np.sin(a), np.cos(d)])
    mu_a = mu_ra_asd * (np.pi / 180 / 3600)      # (cosdec * dRA/dt) rad/day
    mu_d = mu_dec_asd * (np.pi / 180 / 3600)
    r = np.asarray(obs_pos) + rho * los
    v = np.asarray(obs_vel) + rho_dot * los + rho * (mu_a * e_ra + mu_d * e_de)
    return r, v


def helio_linc(tracklets, obs, sun_fn, rho_grid=None, rhodot_grid=None, epoch=None, eps_pos=0.02, eps_vel=0.002, min_samples=2):
    """Link tracklets observed at different epochs (nights).  tracklets: dicts with ra0, dec0, vx_asd, vy_asd,
    t_ref (MJD UTC->TDB approx), and obs = list of (obs_pos(3), obs_vel(3)) barycentric in AU, AU/d, at t_ref
    (same order); sun_fn(jd)-> barycentric sun state.  Returns list of clusters (index lists) with the best hypothesis.

    For each (rho, rhodot) hypothesis (the same for all tracklets, heliocentric distance hypotheses in the original
    method; here topocentric range which is equivalent for the grid) every tracklet is converted to a heliocentric
    state at its t_ref, propagated (2-body) to `epoch` and the 6D points are clustered (DBSCAN, scaled by eps_pos [AU],
    eps_vel [AU/d])."""
    from sklearn.cluster import DBSCAN
    n = len(tracklets)
    if n < 2:
        return []
    jd = np.array([t["jd_tdb"] for t in tracklets])
    epoch = float(np.mean(jd)) if epoch is None else epoch
    rho_grid = np.linspace(0.5, 6.0, 12) if rho_grid is None else rho_grid
    rhodot_grid = np.array([-0.01, -0.005, 0, 0.005, 0.01]) if rhodot_grid is None else rhodot_grid
    best = []
    for rho in rho_grid:
        for rd in rhodot_grid:
            X = np.zeros((n, 6))
            for i, t in enumerate(tracklets):
                op, ov = obs[i]
                r, v = attributable_state(t["ra0"], t["dec0"], t["vx_asd"], t["vy_asd"], op, ov, rho, rd)
                s = sun_fn(t["jd_tdb"])
                rh, vh = r - s[:3], v - s[3:]
                rr, vv = K.propagate_2body(rh, vh, epoch - t["jd_tdb"])
                X[i] = np.r_[rr, vv]
            Xs = X / np.array([eps_pos] * 3 + [eps_vel] * 3)
            lab = DBSCAN(eps=1.0, min_samples=min_samples).fit_predict(Xs)
            for L in set(lab) - {-1}:
                idx = np.where(lab == L)[0]
                if len(set(tracklets[i]["night"] for i in idx)) >= 2:
                    best.append(dict(idx=idx.tolist(), rho=float(rho), rhodot=float(rd), size=len(idx)))
    # keep unique largest clusters
    best.sort(key=lambda c: -c["size"])
    seen = set(); out = []
    for c in best:
        key = tuple(c["idx"])
        if any(i in seen for i in key):
            continue
        out.append(c); seen.update(key)
    return out
