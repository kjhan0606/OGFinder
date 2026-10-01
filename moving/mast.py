"""MAST queries and downloads (astroquery.mast, with plain-HTTP fallback).

* query_candidates(): high-resolution imaging covering a position, HST (ACS/WFC3/WFPC2),
  JWST (NIRCam), PS1 (stacks via the PS1 image service), with the date, filter,
  exposure time, pixel scale and an approximate depth.
* download(): calibrated products (flc/flt/i2d/drc/drz/cal/stk) into ~/.ds9/mast_cache
  (override with OGF_MAST_CACHE), skipping files already present.
"""
import os, re, io, json, math
import numpy as np
from .util import MAST_CACHE, log, ensure_dir, http_get

PIXSCALE = {  # arcsec/pixel of the native (calibrated) grid
    "ACS/WFC": 0.05, "ACS/HRC": 0.027, "WFC3/UVIS": 0.04, "WFC3/IR": 0.128,
    "WFPC2/WFC": 0.0996, "WFPC2/PC": 0.0455, "NIRCAM": 0.031, "NIRCAM/LONG": 0.063,
    "NIRISS": 0.0656, "GPC1": 0.25, "STIS/CCD": 0.05, "WFPC2": 0.0996,
}
# approximate 5-sigma point-source AB depth for a 1000 s exposure (very rough, flagged "approx")
DEPTH1000 = {"ACS/WFC": 27.2, "WFC3/UVIS": 26.8, "WFC3/IR": 26.3, "WFPC2/WFC": 25.5,
             "NIRCAM": 28.3, "GPC1": 23.0}
CAL_PATTERNS = {
    "HST": ["_flc.fits", "_flt.fits", "_drc.fits", "_drz.fits"],
    "JWST": ["_cal.fits", "_i2d.fits", "_rate.fits"],
    "PS1": [".fits"],
}


def approx_depth(inst, texp):
    base = DEPTH1000.get(inst) or DEPTH1000.get(str(inst).split("/")[0])
    if base is None or not texp or texp <= 0:
        return float("nan")
    return base + 1.25 * math.log10(texp / 1000.0)


def pixscale_of(inst):
    inst = str(inst)
    for k, v in PIXSCALE.items():
        if inst.upper().startswith(k.upper()):
            return v
    return float("nan")


def _parse_region(s):
    out = []
    for p in re.split(r"POLYGON", str(s)):
        if not p.strip():
            continue
        p = p.replace("ICRS", " ").replace("J2000", " ").replace("UNKNOWNFrame", " ")
        n = re.findall(r"-?\d+\.\d+(?:[eE][-+]?\d+)?|-?\d+", p)
        v = np.array(n, float)
        if len(v) >= 6:
            out.append(v[:len(v) // 2 * 2].reshape(-1, 2))
    return out


def covers(s_region, ra, dec):
    """True if the footprint polygon(s) in `s_region` contain (ra, dec) [deg]."""
    from matplotlib.path import Path
    for poly in _parse_region(s_region):
        # unwrap RA around the target
        p = poly.copy()
        p[:, 0] = ((p[:, 0] - ra + 180) % 360) - 180 + ra
        if Path(p).contains_point((ra, dec)):
            return True
    return False


def query_candidates(ra, dec, radius_arcmin=1.0, collections=("HST", "JWST", "PS1"),
                     instruments=None, filters=None, calib_level=(2, 3),
                     exposures_only=False, footprint_check=True):
    """Return list of dict candidates sorted by date.  ra, dec, radius in deg/deg/arcmin."""
    from astroquery.mast import Observations
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    c = SkyCoord(ra, dec, unit="deg")
    kw = dict(coordinates=c, radius=radius_arcmin * u.arcmin, obs_collection=list(collections),
              dataproduct_type="image", calib_level=list(calib_level))
    if instruments:
        kw["instrument_name"] = list(instruments)
    if filters:
        kw["filters"] = list(filters)
    tab = Observations.query_criteria(**kw)
    out = []
    for r in tab:
        coll = str(r["obs_collection"]); oid = str(r["obs_id"])
        inst = str(r["instrument_name"])
        if coll == "HST":
            # keep individual-exposure association rows (9-char IPPPSSOOT) and HAP products
            is_exp = bool(re.match(r"^[a-z0-9]{9}$", oid))
            is_hap = oid.startswith("hst_")
            if not (is_exp or is_hap):
                continue
            if exposures_only and not is_exp:
                continue
            if is_hap and (oid.endswith("total") or "_total_" in oid):
                continue  # detection image duplicates the filter products
        if str(r["filters"]).lower() in ("detection", "--"):
            continue
        texp = float(r["t_exptime"]) if r["t_exptime"] is not np.ma.masked else float("nan")
        if footprint_check and r["s_region"] is not np.ma.masked and not covers(r["s_region"], ra, dec):
            continue
        tm = 0.5 * (float(r["t_min"]) + float(r["t_max"]))
        out.append(dict(obs_id=oid, obsid=str(r["obsid"]), collection=coll, instrument=inst,
                        filter=str(r["filters"]), texp=texp, mjd=tm, mjd_start=float(r["t_min"]),
                        pixscale=pixscale_of(inst), depth=approx_depth(inst, texp),
                        proposal=str(r["proposal_id"]), calib_level=int(r["calib_level"]),
                        target=str(r["target_name"]),
                        ra=float(r["s_ra"]), dec=float(r["s_dec"]), s_region=str(r["s_region"]),
                        kind=("exposure" if re.match(r"^[a-z0-9]{9}$", oid) else "coadd")))
    out.sort(key=lambda d: d["mjd"])
    return out


def group_epochs(cands, gap_days=1.0):
    """Assign an epoch index to candidates (clusters in time separated by > gap_days)."""
    ep = -1; last = None
    for c in sorted(cands, key=lambda d: d["mjd"]):
        if last is None or c["mjd"] - last > gap_days:
            ep += 1
        c["epoch"] = ep
        last = c["mjd"]
    return cands


def list_products(obsid, suffixes=None, collection="HST"):
    from astroquery.mast import Observations
    from astropy.table import Table
    tab = Observations.query_criteria(obsid=obsid)
    prods = Observations.get_product_list(tab)
    suff = suffixes or CAL_PATTERNS.get(collection, [".fits"])
    sel = [p for p in prods if any(str(p["productFilename"]).endswith(s) for s in suff)
           and "hlet" not in str(p["productFilename"])]
    plain = [p for p in sel if not str(p["productFilename"]).startswith("hst_")]
    return plain if plain else sel      # drop the HAP-named duplicates of the same exposure


def download(obsid, suffixes=None, collection="HST", dest=None, max_files=None):
    """Download calibrated FITS products of one MAST obsid; returns list of local paths."""
    from astroquery.mast import Observations
    dest = ensure_dir(dest or MAST_CACHE)
    prods = list_products(obsid, suffixes, collection)
    if max_files:
        prods = prods[:max_files]
    paths = []
    for p in prods:
        fn = str(p["productFilename"])
        local = os.path.join(dest, fn)
        if os.path.exists(local) and os.path.getsize(local) > 0:
            paths.append(local); continue
        uri = str(p["dataURI"])
        log("MAST download %s (%.1f MB)" % (fn, float(p["size"]) / 1e6 if p["size"] else -1))
        res = Observations.download_file(uri, local_path=local, cache=False)
        if str(res[0]).upper() not in ("COMPLETE",):
            log("download problem: %s" % (res,))
            continue
        paths.append(local)
    return paths


def ps1_cutout(ra, dec, size_arcsec, filt="i", dest=None):
    """PS1 stack cutout through the ps1images service (FITS).  Returns local path."""
    import requests
    dest = ensure_dir(dest or MAST_CACHE)
    size = int(round(size_arcsec / 0.25))
    url = "https://ps1images.stsci.edu/cgi-bin/ps1cutouts"
    fn = os.path.join(dest, "ps1_%s_%.5f_%.5f_%d.fits" % (filt, ra, dec, size))
    if os.path.exists(fn):
        return fn
    u = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?red=%s&format=fits&size=%d&ra=%.6f&dec=%.6f&"
         "output_size=0&autoscale=0" % (
             _ps1_file(ra, dec, filt), size, ra, dec))
    r = requests.get(u, timeout=120)
    if r.status_code != 200 or len(r.content) < 2880:
        raise RuntimeError("PS1 cutout failed: %s" % r.status_code)
    with open(fn, "wb") as f:
        f.write(r.content)
    return fn


def _ps1_file(ra, dec, filt):
    t = http_get("https://ps1images.stsci.edu/cgi-bin/ps1filenames.py",
                 params=dict(ra=ra, dec=dec, filters=filt, type="stack"))
    lines = t.strip().splitlines()
    return lines[1].split()[7]
