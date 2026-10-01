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


def link_exposures(dets, tol_arcsec=0.6, min_exposures=3, rate_min_ash=0.5, rate_max_ash=600.0,
                   sig_floor=0.01, max_tracklets=500, obs_off_au=None, inv_delta_grid=tuple(np.r_[0.0, np.arange(0.08, 1.2, 0.04), np.arange(1.3, 3.01, 0.2)]),
                   max_cand_per_exposure=3000, min_disp_arcsec=1.5, flux_tol_dex=0.35):
    """dets: list of dict(ex=int exposure index, t=MJD, ra, dec, sig [arcsec], id, flux).
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
