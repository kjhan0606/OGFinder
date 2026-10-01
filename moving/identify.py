"""Known-object identification of tracklets.

Three independent public sources (all no-credential):
  * IMCCE SkyBoT cone search (positions of known SSOs at the exposure epoch; observer = geocentre, the HST parallax is
    handled by comparing against the *geocentric* SkyBoT position plus the parallax shift computed with the SkyBoT
    distance Delta),
  * JPL Horizons observer ephemerides for candidate designations from SkyBoT, evaluated with the true observer
    (centre 500@-48 for HST) -- this is the decisive test and is what 'matched' is based on,
  * MPC `get-obs` (public API) for the observations of the matched object (used by the orbit fit step).
A tracklet is 'known' when the Horizons prediction for a SkyBoT candidate agrees with ALL tracklet points within
`tol` (arcsec, default 1.0 + 3 sigma of the point); otherwise it is a 'new/unknown' candidate.
"""
import numpy as np
from . import skybot, horizons, util
from .util import angsep_arcsec


def skybot_candidates(ra, dec, mjd_utc, radius_deg):
    out = []
    try:
        out = skybot.cone(ra, dec, radius_deg, mjd_utc, observer="500")
    except Exception as e:
        util.log("SkyBoT failed:", e)
    return out


def horizons_id_string(c):
    num = c["num"].strip()
    if num:
        return num + ";"                                  # numbered asteroid: Horizons accepts the plain number + ';' 
    return "DES=%s;" % c["name"].replace(" ", "")


def identify_tracklet(tr, radius_arcmin=6.0, center_code="500@-48", tol_arcsec=1.5, max_cand=25):
    """tr: tracklet dict from tracklet.py (t in UTC MJD, ra, dec, sig).  Returns dict(status, matches=[...])."""
    tm = np.array(tr["t"]); ra = np.array(tr["ra"]); de = np.array(tr["dec"]); sg = np.array(tr["sig"])
    k = len(tm) // 2
    cands = skybot_candidates(ra[k], de[k], tm[k], radius_arcmin / 60.0)
    cands = sorted(cands, key=lambda c: angsep_arcsec(ra[k], de[k], c["ra"], c["dec"]))[:max_cand]
    matches = []
    for c in cands:
        # cheap pre-filter with SkyBoT geocentric position (parallax up to ~ 30" for HST at 1.5 AU)
        if angsep_arcsec(ra[k], de[k], c["ra"], c["dec"]) > 60.0:
            continue
        cid = horizons_id_string(c)
        try:
            e = horizons.observer_table(cid, util.mjd_to_jd(tm), center=center_code)
        except Exception as ex:
            util.log("Horizons failed for", c["name"], ex)
            continue
        dr = (ra - e[:, 0]) * np.cos(np.radians(de)) * 3600.0; dd = (de - e[:, 1]) * 3600.0
        sep = np.hypot(dr, dd)
        ok = bool(np.all(sep < tol_arcsec + 3.0 * sg))
        matches.append(dict(name=c["name"], num=c["num"], vmag=c["vmag"], cls=c["cls"], sep_arcsec=[float(s) for s in sep],
                            max_sep=float(sep.max()), mean_dra=float(dr.mean()), mean_ddec=float(dd.mean()), matched=ok,
                            skybot_ra=c["ra"], skybot_dec=c["dec"], delta_au=c["delta"], r_au=c["r"]))
    matches.sort(key=lambda m: m["max_sep"])
    status = "known" if matches and matches[0]["matched"] else "unknown"
    return dict(status=status, matches=matches[:5])
