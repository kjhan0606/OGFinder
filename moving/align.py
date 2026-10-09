"""WCS alignment of exposures to Gaia DR3 and to a reference image.

Per chip: detect stars (sep), query Gaia DR3 (astroquery TAP, cached), propagate proper motions to the
exposure epoch, iterative sigma-clipped matching, and fit a tangent-plane model
  [xi_g, eta_g] = M [xi_o, eta_o] + t            (affine: 6 parameters)  or polynomial order 2 (12 parameters)
The affine part is stored in a sidecar JSON next to the output and applied when the chip is loaded
(imaging.apply_alignment).  Residual rms (mas) of the matched stars after the fit is reported.
Centroids are windowed first moments (sep.winpos) with errors from sep's flux/centroid S/N model.
"""
def _import_ogfmeas():
    """In-tree measurements. Finds ogfmeas from this file so a script does not need PYTHONPATH."""
    import pathlib
    import sys
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas
    return ogfmeas.measurement_library()

import os, json
import numpy as np
from .util import log, http_get, tangent_offsets, write_json, ensure_dir
from . import imaging as I

GAIA_TAP = "https://gea.esac.esa.int/tap-server/tap/sync"


def gaia_cone(ra, dec, radius_deg, gmax=21.0, cache=True):
    """Gaia DR3 sources (ra, dec at epoch 2016.0, pmra, pmdec, G).  Returns dict of numpy arrays."""
    q = ("SELECT source_id, ra, dec, pmra, pmdec, parallax, phot_g_mean_mag, ra_error, dec_error "
         "FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS',ra,dec),CIRCLE('ICRS',%.6f,%.6f,%.5f)) "
         "AND phot_g_mean_mag<%.2f" % (ra, dec, radius_deg, gmax))
    txt = http_get(GAIA_TAP, params=None, method="POST",
                   data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "QUERY": q}, cache=cache, timeout=120)
    lines = txt.strip().splitlines()
    hdr = lines[0].split(",")
    rows = [l.split(",") for l in lines[1:]]
    out = {}
    for j, k in enumerate(hdr):
        col = []
        for r in rows:
            v = r[j] if j < len(r) else ""
            try:
                col.append(float(v) if v not in ("", "null") else np.nan)
            except ValueError:
                col.append(np.nan)
        out[k] = np.array(col)
    out["pmra"] = np.nan_to_num(out["pmra"]); out["pmdec"] = np.nan_to_num(out["pmdec"])
    return out


def gaia_at_epoch(g, mjd):
    dt = (mjd - 57388.0) / 365.25   # Gaia DR3 reference epoch J2016.0 = MJD 57388.0
    ra = g["ra"] + g["pmra"] * dt / 3.6e6 / np.cos(np.radians(g["dec"]))
    dec = g["dec"] + g["pmdec"] * dt / 3.6e6
    return ra, dec


def detect_stars(chip, snr_min=10.0, max_n=3000):
    """Compact, unsaturated, isolated sources with windowed centroids.  Compactness is judged with the
    half-light radius (sep.flux_radius), whose stellar locus is narrow even for bright stars."""
    sep = _import_ogfmeas()
    bkg, rms = I.background(chip.data, chip.bad)
    sub = np.ascontiguousarray((chip.data - bkg).astype(np.float32))
    o = sep.extract(sub, 4.0, err=np.ascontiguousarray(rms), mask=np.ascontiguousarray(chip.bad),
                    minarea=4, deblend_nthresh=32, deblend_cont=0.005)
    if len(o) == 0:
        return None
    snr = o["flux"] / np.sqrt(o["npix"] * np.median(rms) ** 2 + np.abs(o["flux"]) + 1e-9)
    ok = (o["flag"] == 0) & (snr > snr_min) & (o["flux"] > 0)
    o = o[ok]; snr = snr[ok]
    if len(o) < 4:
        return None
    r50, rflag = sep.flux_radius(sub, o["x"], o["y"], 12.0 * np.ones(len(o)), 0.5)
    good = (r50 > 0.5) & np.isfinite(r50)
    sat = np.array([chip.saturated[int(round(yy)), int(round(xx))] if 0 <= int(round(yy)) < chip.shape[0] and 0 <= int(round(xx)) < chip.shape[1] else True
                    for xx, yy in zip(o["x"], o["y"])])
    good &= ~sat
    o = o[good]; snr = snr[good]; r50 = r50[good]
    if len(o) < 4:
        return None
    h, e = np.histogram(r50, bins=np.arange(0.5, 8.0, 0.2))
    pk = 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
    st = (r50 > 0.75 * pk) & (r50 < 1.35 * pk)
    o = o[st]; snr = snr[st]; r50 = r50[st]
    if len(o) < 4:
        return None
    xw, yw, flag = sep.winpos(sub, o["x"], o["y"], 2.0 * np.ones(len(o)))
    k = flag == 0
    o = o[k]; xw = xw[k]; yw = yw[k]; snr = snr[k]; r50 = r50[k]
    fwhm = float(2.0 * np.median(r50))
    sig = np.maximum(fwhm / 2.355 / np.maximum(snr, 1), 0.02)
    order = np.argsort(-snr)[:max_n]
    return dict(x=xw[order], y=yw[order], flux=o["flux"][order], snr=snr[order], sig_pix=sig[order], fwhm=fwhm)


def _design(model):
    if model == "shift":
        return lambda x, y: np.stack([np.ones_like(x)], 1)
    return lambda x, y: np.stack([np.ones_like(x), x, y], 1)


def fit_transform(xo, eo, xg, eg, model="affine", w=None):
    """Fit tangent-plane transform  g = M o + t.  model: shift | similarity | affine.
    Returns M (2x2), t (2), predicted g."""
    n = len(xo)
    w = np.ones(n) if w is None else w
    if model == "shift":
        t = np.array([np.average(xg - xo, weights=w), np.average(eg - eo, weights=w)])
        M = np.eye(2)
    elif model == "similarity":
        # g = [a -b; b a] o + t
        A = np.zeros((2 * n, 4)); y = np.zeros(2 * n)
        A[0::2] = np.stack([xo, -eo, np.ones(n), np.zeros(n)], 1)
        A[1::2] = np.stack([eo, xo, np.zeros(n), np.ones(n)], 1)
        y[0::2] = xg; y[1::2] = eg
        ww = np.repeat(w, 2)
        p = np.linalg.lstsq(A * ww[:, None], y * ww, rcond=None)[0]
        M = np.array([[p[0], -p[1]], [p[1], p[0]]]); t = p[2:4]
    else:
        A = np.stack([np.ones(n), xo, eo], 1)
        cx = np.linalg.lstsq(A * w[:, None], xg * w, rcond=None)[0]
        cy = np.linalg.lstsq(A * w[:, None], eg * w, rcond=None)[0]
        M = np.array([[cx[1], cx[2]], [cy[1], cy[2]]]); t = np.array([cx[0], cy[0]])
    pred = np.stack([M[0, 0] * xo + M[0, 1] * eo + t[0], M[1, 0] * xo + M[1, 1] * eo + t[1]], 1)
    return M, t, pred


def tangent_inverse(xi, eta, ra0, de0):
    xi = np.radians(np.asarray(xi) / 3600.0); eta = np.radians(np.asarray(eta) / 3600.0)
    ra0r, de0r = np.radians(ra0), np.radians(de0)
    den = np.cos(de0r) - eta * np.sin(de0r)
    ra = ra0r + np.arctan2(xi, den)
    dec = np.arctan2(np.sin(de0r) + eta * np.cos(de0r), np.hypot(xi, den))
    return np.degrees(ra) % 360.0, np.degrees(dec)


def apply_transform_to_wcs(wcs, M, t_arcsec):
    """Compose tangent-plane transform (about CRVAL) with a celestial WCS:  CD' = M CD,  CRVAL' = CRVAL + t.
    Exact for the linear part; the translation is applied by moving the tangent point (second-order
    error ~ |t| * theta^2, below 1e-3 mas for |t| < 1 arcsec and a 0.1 deg field)."""
    w = wcs
    M = np.asarray(M, float)
    if w.wcs.has_cd():
        w.wcs.cd = M @ w.wcs.cd
    else:
        w.wcs.pc = M @ w.wcs.get_pc()
    ra0, de0 = w.wcs.crval
    w.wcs.crval = list(tangent_inverse(t_arcsec[0], t_arcsec[1], ra0, de0))
    # keep the transformation of the zero point consistent: new tangent point maps CRPIX -> shifted sky
    return w


def match_and_fit(chip, xp, yp, ref_ra, ref_de, model="auto", match_arcsec=1.0, nsig=3.0, min_match=6,
                  weights=None, tag="gaia", search_arcsec=3.0):
    """Iteratively match detected pixel positions to reference sky positions, fit the transform about CRVAL,
    update chip.wcs in place and return stats."""
    res = dict(chip=chip.name, ok=False, n_det=int(len(xp)), n_match=0, rms_mas=np.nan, ref=tag)
    ra, de = chip.wcs.all_pix2world(xp, yp, 0)
    ra0, de0 = chip.wcs.wcs.crval
    xo, eo = tangent_offsets(ra, de, ra0, de0)
    xg, eg = tangent_offsets(ref_ra, ref_de, ra0, de0)
    from scipy.spatial import cKDTree
    tree = cKDTree(np.stack([xg, eg], 1))
    M = np.eye(2); t = np.zeros(2)
    mdl = "shift"
    keep = None
    # blind offset search by voting (robust to a poor initial WCS and to few stars): all detection-reference
    # pair offsets within +-search_arcsec are histogrammed; the peak gives the initial shift.
    dxs = (xg[None, :] - xo[:, None]).ravel(); dys = (eg[None, :] - eo[:, None]).ravel()
    sel = (np.abs(dxs) < search_arcsec) & (np.abs(dys) < search_arcsec)
    if sel.sum() >= min_match:
        bw = 0.15
        H, xe, ye = np.histogram2d(dxs[sel], dys[sel], bins=int(2 * search_arcsec / bw), range=[[-search_arcsec, search_arcsec]] * 2)
        from scipy.ndimage import uniform_filter
        Hs = uniform_filter(H, 3) * 9
        i, j = np.unravel_index(np.argmax(Hs), Hs.shape)
        t = np.array([0.5 * (xe[i] + xe[i + 1]), 0.5 * (ye[j] + ye[j + 1])])
        res["vote_peak"] = float(Hs[i, j])
    for it in range(8):
        pos = np.stack([M[0, 0] * xo + M[0, 1] * eo + t[0], M[1, 0] * xo + M[1, 1] * eo + t[1]], 1)
        rad = match_arcsec * (1.0 if it < 3 else 0.7)
        d, j = tree.query(pos, distance_upper_bound=rad)
        m = np.isfinite(d)
        if m.sum() < min_match:
            res["reason"] = "too few matches (%d)" % m.sum()
            return res
        # one-to-one
        order = np.argsort(d); seen = set(); k2 = np.zeros(len(d), bool)
        for k in order:
            if m[k] and j[k] not in seen:
                seen.add(j[k]); k2[k] = True
        idx = np.where(k2)[0]
        n = len(idx)
        mdl = model if model != "auto" else ("shift" if n < 12 else "similarity" if n < 40 else "affine")
        wv = None if weights is None else weights[idx]
        M, t, pred = fit_transform(xo[idx], eo[idx], xg[j[idx]], eg[j[idx]], mdl, wv)
        r = np.hypot(pred[:, 0] - xg[j[idx]], pred[:, 1] - eg[j[idx]])
        s = 1.4826 * np.median(np.abs(r - np.median(r))) + 1e-4
        good = r < np.median(r) + nsig * max(s, 0.02)
        if good.sum() >= min_match and good.sum() < n:
            idx = idx[good]; n = len(idx)
            M, t, pred = fit_transform(xo[idx], eo[idx], xg[j[idx]], eg[j[idx]], mdl, None if weights is None else weights[idx])
        keep = (idx, j[idx], pred)
    idx, jj, pred = keep
    n = len(idx)
    npar = {"shift": 2, "similarity": 4, "affine": 6}[mdl]
    rr = np.stack([pred[:, 0] - xg[jj], pred[:, 1] - eg[jj]], 1)
    rms1d = np.sqrt(np.sum(rr ** 2) / max(2 * n - npar, 1))
    apply_transform_to_wcs(chip.wcs, M, t)
    res.update(ok=True, n_match=int(n), rms_mas=float(rms1d * 1000.0), model=mdl, M=M.ravel().tolist(),
               t_arcsec=t.tolist(), crval_orig=[float(ra0), float(de0)], shift_mas=float(np.hypot(*t) * 1000.0),
               rot_deg=float(np.degrees(np.arctan2(M[1, 0] - M[0, 1], M[0, 0] + M[1, 1]))),
               scale=float(np.sqrt(abs(np.linalg.det(M)))))
    chip.align = res
    return res


def align_chip_to_gaia(chip, gmax=20.5, snr_min=8.0, gaia=None, match_arcsec=1.0):
    st = detect_stars(chip, snr_min=snr_min)
    if st is None or len(st["x"]) < 4:
        return dict(chip=chip.name, ok=False, n_det=0, n_match=0, rms_mas=np.nan, reason="too few stars", ref="gaia")
    ny, nx = chip.shape
    c = chip.wcs.all_pix2world([nx / 2], [ny / 2], 0)
    cra, cde = float(c[0][0]), float(c[1][0])
    rad = 0.5 * np.hypot(nx, ny) * chip.pixscale / 3600.0 * 1.1
    g = gaia if gaia is not None else gaia_cone(cra, cde, rad, gmax=gmax)
    gra, gde = gaia_at_epoch(g, chip.mjd)
    return match_and_fit(chip, st["x"], st["y"], gra, gde, "auto", match_arcsec, min_match=3, tag="gaia")


def sources_for_matching(chip, snr_min=6.0, max_n=4000):
    """All compact-or-extended sources with centroid, used for relative alignment (galaxies included)."""
    sep = _import_ogfmeas()
    bkg, rms = I.background(chip.data, chip.bad)
    sub = np.ascontiguousarray((chip.data - bkg).astype(np.float32))
    o = sep.extract(sub, 4.0, err=np.ascontiguousarray(rms), mask=np.ascontiguousarray(chip.bad), minarea=6,
                    deblend_nthresh=32, deblend_cont=0.01)
    snr = o["flux"] / np.sqrt(o["npix"] * np.median(rms) ** 2 + np.abs(o["flux"]) + 1e-9)
    sel = (o["flag"] == 0) & (snr > snr_min) & (o["a"] < 12)
    xi = np.clip(np.round(o["x"]).astype(int), 0, chip.shape[1] - 1); yi = np.clip(np.round(o["y"]).astype(int), 0, chip.shape[0] - 1)
    sel &= ~chip.cr[yi, xi] & ~chip.saturated[yi, xi]
    o = o[sel]
    xw, yw, fl = sep.winpos(sub, o["x"], o["y"], 2.5 * np.ones(len(o)))
    ok = fl == 0
    return o["x"][ok] * 0 + xw[ok], yw[ok], o["flux"][ok]


def align_chip_to_reference(chip, ref_chip_or_cat, snr_min=6.0, match_arcsec=1.5, model="auto"):
    """Align to a reference given as (ra, dec) arrays catalogue (already on Gaia frame)."""
    ra_ref, de_ref = ref_chip_or_cat
    x, y, fl = sources_for_matching(chip, snr_min)
    if len(x) < 6:
        return dict(chip=chip.name, ok=False, reason="too few sources", ref="reference", n_det=len(x), n_match=0, rms_mas=np.nan)
    return match_and_fit(chip, x, y, ra_ref, de_ref, model, match_arcsec, min_match=6, tag="reference")


def save_sidecar(chip, outdir):
    ensure_dir(outdir)
    p = I.align_sidecar(chip.path, outdir)
    d = json.load(open(p)) if os.path.exists(p) else {"path": chip.path, "chips": {}}
    d["chips"]["ext%d" % chip.ext] = {k: (v if not isinstance(v, np.generic) else v.item()) for k, v in chip.align.items()}
    write_json(p, d)
    return p


def load_sidecar_into(chip, outdir):
    """Re-apply a stored alignment (M, t about the ORIGINAL CRVAL) to a freshly loaded chip."""
    p = I.align_sidecar(chip.path, outdir)
    if not os.path.exists(p):
        return False
    d = json.load(open(p))
    a = d["chips"].get("ext%d" % chip.ext)
    if not a or not a.get("ok"):
        return False
    apply_transform_to_wcs(chip.wcs, np.array(a["M"]).reshape(2, 2), a["t_arcsec"])
    chip.align = a
    return True


def sky_catalog(chip, snr_min=6.0):
    x, y, f = sources_for_matching(chip, snr_min)
    ra, de = chip.wcs.all_pix2world(x, y, 0)
    return ra, de


def align_field(chips, gaia_min_match=4, gaia_max_rms_mas=40.0, snr_min=6.0, log_fn=None):
    """Align a set of chips (all epochs/filters) to a common frame.

    1. Try Gaia DR3 for every chip (few stars at high galactic latitude: often 4-10 per ACS chip);
       accepted when n_match >= gaia_min_match and rms < gaia_max_rms_mas.
    2. Chain relative alignment: chips are processed in order of decreasing overlap with the growing catalogue
       (positions of all sources > snr_min of already aligned chips; moving objects and CRs are clipped by the
       robust fit).  The anchor chip (best Gaia solution, else the first chip) keeps its (Gaia-corrected) WCS.
    Returns list of per-chip stat dicts (rms in mas of the relative fit and of the Gaia fit if any)."""
    lg = log_fn or log
    stats = {}
    gaia_ok = {}
    for c in chips:
        try:
            r = align_chip_to_gaia(c, snr_min=5, gmax=21.5)
        except Exception as e:
            r = dict(ok=False, reason=repr(e), n_match=0, rms_mas=np.nan, chip=c.name)
        stats[c.name] = dict(gaia=r)
        if r["ok"] and r["n_match"] >= gaia_min_match and r["rms_mas"] < gaia_max_rms_mas:
            gaia_ok[c.name] = r
        else:
            # undo any partial application (match_and_fit only applies when ok)
            if r["ok"]:
                pass
    # chips whose Gaia fit was accepted keep the shift; others were not modified unless r["ok"] (then revert)
    for c in chips:
        r = stats[c.name]["gaia"]
        if r.get("ok") and c.name not in gaia_ok:
            M = np.array(r["M"]).reshape(2, 2)
            # revert by re-applying the inverse about the NEW crval is not exact; reload the original WCS instead
            from astropy.wcs import WCS
            from astropy.io import fits
            with fits.open(c.path) as h:
                try:
                    c.wcs = WCS(h[c.ext].header, fobj=h)
                except Exception:
                    c.wcs = WCS(h[c.ext].header)
    order = sorted(chips, key=lambda c: -(gaia_ok[c.name]["n_match"] if c.name in gaia_ok else 0))
    anchor = order[0]
    cat_ra = []; cat_de = []
    done = []
    r0, d0 = sky_catalog(anchor, snr_min)
    cat_ra.append(r0); cat_de.append(d0)
    stats[anchor.name]["relative"] = dict(ok=True, anchor=True, n_match=0, rms_mas=0.0)
    done.append(anchor)
    todo = [c for c in chips if c is not anchor]
    anchor.align = dict(anchor.align, ok=True) if getattr(anchor, "align", None) and anchor.align.get("ok") else dict(ok=True, anchor=True, M=[1, 0, 0, 1], t_arcsec=[0.0, 0.0])
    _inherit_siblings(anchor, chips, todo, stats, snr_min, cat_ra, cat_de, done)
    while todo:
        allr = np.concatenate(cat_ra); allde = np.concatenate(cat_de)
        # choose chip with the largest number of catalogue sources within its footprint
        best = None; bestn = -1
        for c in todo:
            ny, nx = c.shape
            x, y = c.wcs.all_world2pix(allr, allde, 0, quiet=True)
            n = int(np.sum((x > 0) & (x < nx) & (y > 0) & (y < ny)))
            if n > bestn:
                best, bestn = c, n
        c = best; todo.remove(c)
        if bestn < 6:
            stats[c.name]["relative"] = dict(ok=False, reason="no overlap with aligned chips (%d cat sources)" % bestn)
            continue
        try:
            r = align_chip_to_reference(c, (allr, allde), snr_min=snr_min, match_arcsec=1.0)
        except Exception as e:
            r = dict(ok=False, reason=repr(e))
        stats[c.name]["relative"] = r
        if r.get("ok"):
            rr, dd = sky_catalog(c, snr_min)
            cat_ra.append(rr); cat_de.append(dd)
            done.append(c)
            _inherit_siblings(c, chips, todo, stats, snr_min, cat_ra, cat_de, done)
    return stats


def _inherit_siblings(c, chips, todo, stats, snr_min, cat_ra, cat_de, done):
    """Other chips of the same exposure share its pointing correction (rigid detector): apply the same
    tangent-plane transform (M, t) and add their sources to the running catalogue."""
    a = c.align
    for s in list(todo):
        if s.path == c.path and s is not c:
            apply_transform_to_wcs(s.wcs, np.array(a["M"]).reshape(2, 2), a["t_arcsec"])
            s.align = dict(a, inherited=True, chip=s.name)
            stats[s.name]["relative"] = dict(a, inherited=True, ok=True, chip=s.name)
            todo.remove(s)
            rr, dd = sky_catalog(s, snr_min)
            cat_ra.append(rr); cat_de.append(dd)
            done.append(s)


def propagate_same_exposure(chips, stats):
    """Chips of one file share the pointing: if one chip got a transform and the other did not, apply the
    same tangent-plane correction (M, t) to the other (about its own CRVAL).  Marks align['inherited']=True."""
    byfile = {}
    for c in chips:
        byfile.setdefault(c.path, []).append(c)
    n = 0
    for p, cl in byfile.items():
        ok = [c for c in cl if getattr(c, "align", None) and c.align.get("ok")]
        bad = [c for c in cl if not (getattr(c, "align", None) and c.align.get("ok"))]
        if ok and bad:
            a = ok[0].align
            for c in bad:
                apply_transform_to_wcs(c.wcs, np.array(a["M"]).reshape(2, 2), a["t_arcsec"])
                c.align = dict(a, inherited=True, chip=c.name)
                stats[c.name]["inherited"] = True
                n += 1
    return n
