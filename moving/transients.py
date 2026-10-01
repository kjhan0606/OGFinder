"""Static-transient candidate search, host association and light-curve heuristics.

Input: single-exposure detections (pipeline.detect_in_region) AFTER removal of detections that belong to linked mover
tracklets or to known solar-system objects.  A static transient candidate is a positive point-like detection that
  * is not classified as artefact (CR / edge / dipole),
  * is found at the same sky position (within `tol` arcsec) in >= `min_epochs` exposures taken on different exposure
    files (cosmic rays do not repeat) -- or in a single exposure if `allow_single`,
  * is not present in the template (by construction of the difference image).
Host association: nearest catalogue galaxy from the OGFinder catalogue TSV (any ds9_sextract-style table with ALPHA/DELTA
(or RA/DEC) columns and optionally a size column and a redshift column); reports separation (arcsec), separation in units of
the host's effective radius (if a radius column exists), and the host redshift/photo-z if a column of that name exists.
Light-curve: forced aperture photometry on the difference image at the candidate position in every exposure.
SN/AGN heuristic (documented; NOT a classifier): nuclear offset < 0.3 arcsec or host radius ratio < 0.3 -> 'nuclear/AGN-like',
rise/decline slope > 0.05 mag/day over the baseline -> 'fast evolving', else 'consistent with constant flux'.
"""
import numpy as np
from . import util
from .util import angsep_arcsec


def group_static(dets, tol_arcsec=0.4, min_epochs=2, classes=("point", "faint"), snr_min=6.0, exclude=None):
    """dets: list of detection dicts (with ra, dec, ex, cls, snr, file).  exclude: set of detection indices to skip."""
    exclude = exclude or set()
    sel = [(i, d) for i, d in enumerate(dets) if i not in exclude and d["sign"] > 0 and d["cls"] in classes and d["snr"] >= snr_min]
    if not sel:
        return []
    from scipy.spatial import cKDTree
    ra0 = np.mean([d["ra"] for _, d in sel]); de0 = np.mean([d["dec"] for _, d in sel])
    xy = np.array([util.tangent_offsets(np.array([d["ra"]]), np.array([d["dec"]]), ra0, de0) for _, d in sel]).reshape(len(sel), 2)
    tree = cKDTree(xy)
    seen = set(); groups = []
    for a in range(len(sel)):
        if a in seen:
            continue
        mem = tree.query_ball_point(xy[a], tol_arcsec)
        mem = [m for m in mem if m not in seen]
        files = {sel[m][1]["file"] for m in mem if "file" in sel[m][1]}
        if len(files) >= min_epochs or (min_epochs <= 1):
            groups.append([sel[m][0] for m in mem])
            seen.update(mem)
    return groups


def read_catalog(path):
    rows = util.read_tsv(path)
    cols = rows[0].keys() if rows else []
    def pick(*names):
        for n in names:
            for c in cols:
                if c.lower() == n.lower():
                    return c
        return None
    return rows, dict(ra=pick("ALPHA_J2000", "X_WORLD", "RA", "ra_deg"), dec=pick("DELTA_J2000", "Y_WORLD", "DEC", "dec_deg"),
                      rad=pick("FLUX_RADIUS", "R_EFF", "A_IMAGE", "KRON_RADIUS"), z=pick("z_phot", "photo_z", "ZPHOT", "z_best", "z"),
                      id=pick("NUMBER", "ID", "id"), mag=pick("MAG_AUTO", "mag"), pix=pick("PIXSCALE"))


def associate_host(ra, dec, cat_rows, cmap, pixscale_arcsec=None, max_arcsec=5.0):
    if not cat_rows or cmap["ra"] is None:
        return None
    best = None
    for r in cat_rows:
        try:
            s = angsep_arcsec(ra, dec, float(r[cmap["ra"]]), float(r[cmap["dec"]]))
        except (TypeError, ValueError):
            continue
        if s < max_arcsec and (best is None or s < best[0]):
            best = (s, r)
    if best is None:
        return None
    s, r = best
    out = dict(sep_arcsec=float(s), id=r.get(cmap["id"]) if cmap["id"] else None)
    if cmap["rad"] and pixscale_arcsec:
        try:
            rr = float(r[cmap["rad"]]) * pixscale_arcsec        # half-light radius in pixels -> arcsec (approximate)
            out["host_radius_arcsec"] = rr; out["offset_in_re"] = s / rr if rr > 0 else np.nan
        except (TypeError, ValueError):
            pass
    if cmap["z"]:
        try:
            out["z"] = float(r[cmap["z"]])
        except (TypeError, ValueError):
            pass
    if cmap["mag"]:
        try:
            out["host_mag"] = float(r[cmap["mag"]])
        except (TypeError, ValueError):
            pass
    return out


def heuristic_class(host, lc):
    why = []
    cls = "SN-like (offset transient)"
    if host is None:
        cls = "hostless / unassociated"; why.append("no catalogue galaxy within 5 arcsec")
    elif host["sep_arcsec"] < 0.3 or host.get("offset_in_re", 9) < 0.3:
        cls = "nuclear (AGN-like or TDE)"; why.append("offset %.2f arcsec from host centre" % host["sep_arcsec"])
    else:
        why.append("offset %.2f arcsec%s" % (host["sep_arcsec"], (" = %.1f R_e" % host["offset_in_re"]) if "offset_in_re" in host else ""))
    if lc is not None and len(lc["t"]) >= 3:
        t = np.array(lc["t"]); m = np.array(lc["mag"], float); ok = np.isfinite(m)
        if ok.sum() >= 3 and np.ptp(t[ok]) > 0.01:
            sl = np.polyfit(t[ok], m[ok], 1)[0]
            why.append("slope %.3f mag/day over %.2f d" % (sl, np.ptp(t[ok])))
            if abs(sl) > 0.05:
                why.append("fast evolving")
    return cls, why


def lightcurve(cands_positions, diff_images, r_pix=3.0):
    """diff_images: list of (t_mjd, wcs, alpha_img, zp_ab, texp).  Returns dict(t, flux, mag, ...) at each position."""
    from .photometry import aperture_flux, ab_mag
    out = []
    for (ra, dec) in cands_positions:
        t = []; fl = []; er = []; mg = []
        for tm, wcs, img, zp, _ in diff_images:
            x, y = wcs.all_world2pix([ra], [dec], 0)
            ny, nx = img.shape
            if not (5 < x[0] < nx - 5 and 5 < y[0] < ny - 5):
                continue
            f, e = aperture_flux(img, x[0], y[0], r_pix)
            t.append(tm); fl.append(f); er.append(e)
            mg.append(float(ab_mag(f, zp)) if (f > 3 * e and zp is not None) else np.nan)
        out.append(dict(t=t, flux=fl, err=er, mag=mg))
    return out
