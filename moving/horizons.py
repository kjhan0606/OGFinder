"""Minimal JPL Horizons REST client (vectors, observer ephemerides) with cache.

Used for: HST (-48) / JWST (-170) geocentric position, ephemeris truth for validation,
and orbit fetches (SBDB).  Everything in ICRF equatorial, AU and AU/day.
"""
import re
import numpy as np
from .util import http_get

HORIZONS = "https://ssd.jpl.nasa.gov/api/horizons.api"
SBDB = "https://ssd-api.jpl.nasa.gov/sbdb.api"


def _hz(**kw):
    p = dict(format="json")
    p.update(kw)
    j = http_get(HORIZONS, params=p, as_json=True)
    if "result" not in j:
        raise RuntimeError("Horizons error: %s" % str(j)[:300])
    return j["result"]


def _block(res):
    i, j = res.find("$$SOE"), res.find("$$EOE")
    if i < 0:
        raise RuntimeError("Horizons returned no ephemeris:\n" + res[:800])
    return res[i + 5:j].strip().splitlines()


def vectors(command, jd_tdb_list, center="500@0", chunk=40):
    """State vectors (x,y,z,vx,vy,vz) [AU, AU/d], ICRF equatorial, at TDB JDs.
    Uses the TLIST option so arbitrary times are allowed."""
    jd_tdb_list = np.atleast_1d(np.asarray(jd_tdb_list, dtype=float))
    out = np.zeros((len(jd_tdb_list), 6))
    for s in range(0, len(jd_tdb_list), chunk):
        sub = jd_tdb_list[s:s + chunk]
        res = _hz(COMMAND="'%s'" % command, OBJ_DATA="NO", MAKE_EPHEM="YES", EPHEM_TYPE="VECTORS",
                  CENTER="'%s'" % center, TLIST="'%s'" % " ".join("%.9f" % t for t in sub),
                  REF_PLANE="FRAME", REF_SYSTEM="ICRF", OUT_UNITS="AU-D", VEC_TABLE="2",
                  CSV_FORMAT="YES")
        rows = _block(res)
        k = 0
        for line in rows:
            f = [x.strip() for x in line.split(",") if x.strip() != ""]
            if len(f) < 8:
                continue
            out[s + k] = [float(x) for x in f[2:8]]
            k += 1
        if k != len(sub):
            raise RuntimeError("Horizons vectors: expected %d rows got %d" % (len(sub), k))
    return out


def observer_table(command, jd_utc_list, center="500", quantities="1", chunk=40):
    """Astrometric RA/Dec (deg, light-time corrected, ICRF) for `command` seen from
    `center` (MPC code like '500' geocentric, '250' HST?, or '@-48'), times in UTC JD.
    Returns array (n,2) [RA, DEC] and extra columns (delta, r, V mag if available)."""
    jd = np.atleast_1d(np.asarray(jd_utc_list, dtype=float))
    out = np.full((len(jd), 5), np.nan)
    for s in range(0, len(jd), chunk):
        sub = jd[s:s + chunk]
        res = _hz(COMMAND="'%s'" % command, OBJ_DATA="NO", MAKE_EPHEM="YES", EPHEM_TYPE="OBSERVER",
                  CENTER="'%s'" % center, TLIST="'%s'" % " ".join("%.9f" % t for t in sub),
                  TLIST_TYPE="JD", TIME_TYPE="UT", QUANTITIES="'1,9,19,20'", ANG_FORMAT="DEG",
                  CSV_FORMAT="YES", EXTRA_PREC="YES", APPARENT="AIRLESS", CAL_FORMAT="JD")
        i = res.find("$$SOE")
        hdr_line = None
        for line in res[:i].splitlines()[::-1]:
            if "R.A." in line or "RA" in line.split(",")[0:4].__str__():
                hdr_line = line
                break
        rows = _block(res)
        cols = [c.strip() for c in hdr_line.split(",")] if hdr_line else []
        k = 0
        for line in rows:
            f = [x.strip() for x in line.split(",")]
            try:
                ira = [j for j, c in enumerate(cols) if c.startswith("R.A.") or c == "R.A._(ICRF)"][0]
            except IndexError:
                ira = 3
            ra = float(f[ira]); de = float(f[ira + 1])
            out[s + k, 0] = ra; out[s + k, 1] = de
            for j, c in enumerate(cols):
                if c.startswith("V") and "mag" in c.lower() or c == "V":
                    try: out[s + k, 2] = float(f[j])
                    except Exception: pass
                if c.startswith("delta"):
                    try: out[s + k, 3] = float(f[j])
                    except Exception: pass
                if c.startswith("r") and "rdot" not in c and c.strip() == "r":
                    try: out[s + k, 4] = float(f[j])
                    except Exception: pass
            k += 1
    return out


def sbdb(sstr, cov=True, phys=True):
    """SBDB query. Returns dict (json)."""
    p = {"sstr": sstr, "full-prec": 1}
    if cov:
        p["cov"] = "mat"
    if phys:
        p["phys-par"] = 1
    return http_get(SBDB, params=p, as_json=True)


def sbdb_elements(j):
    """Return dict of elements from sbdb json: epoch_jd (TDB), a,e,i,om,w,ma,tp,q..."""
    el = {e["name"]: float(e["value"]) for e in j["orbit"]["elements"]}
    el["epoch_jd"] = float(j["orbit"]["epoch"])
    sig = {e["name"]: float(e["sigma"]) for e in j["orbit"]["elements"] if e.get("sigma") not in (None, "")}
    return el, sig


def cartesian_from_elements(el, jd_target=None, mu=None):
    """Heliocentric ecliptic-J2000 Keplerian -> barycentric ICRF equatorial state needs the Sun;
    use heliocentric elements to get a heliocentric equatorial state (2-body).  `el` keys:
    a,e,i,om,w,ma (deg), epoch_jd.  Returns heliocentric equatorial state at epoch."""
    from .kepler import elements_to_state
    return elements_to_state(el["a"], el["e"], el["i"], el["om"], el["w"], el["ma"], ecliptic=True)
