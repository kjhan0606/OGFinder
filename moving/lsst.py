"""Rubin / LSST reference provider (public, no-credential access only).

Findings (probed 2026-10-01 from this box, see docs/moving_objects.md):
  * Rubin DP1 images and catalogs (SIA v2, TAP, SODA cutouts, HiPS tiles) require a Rubin Science
    Platform login (HTTP 401 without a token).  Only data-rights holders may have an RSP account
    (dp1.lsst.io, "Data access").  This module therefore never sends or stores credentials.
  * Public without login: the TAP availability document, the HiPS *list* (metadata only), and the
    CDS "Rubin First Look" HiPS (colour PNG art built from press-release TIFFs; NOT science data).
  * A user who holds data rights can download DP1 calexp/coadd FITS themselves; point
    --lsst-local-dir at that folder and they are used as a reference/epoch source (no network).

So the provider is: probe -> coverage check against the seven DP1 fields -> local files if any ->
otherwise marked unavailable with the reason (callers fall back to MAST / Pan-STARRS).
"""
import os, glob
import numpy as np
from .util import http_get, log, angsep_arcsec

TAP_AVAIL = "https://data.lsst.cloud/api/tap/availability"
HIPS_LIST = "https://data.lsst.cloud/api/hips/v2/dp1/list"
SIA_DP1 = "https://data.lsst.cloud/api/sia/dp1/query"

# DP1 fields (dp1.lsst.io/overview/observations.html): name, RA, Dec [deg]; ~1 deg^2 each.
DP1_FIELDS = [
    ("47 Tuc", 6.02, -72.08), ("Low Ecliptic Latitude", 37.86, 6.98),
    ("Fornax dSph", 40.00, -34.45), ("ECDFS", 53.13, -28.10),
    ("EDFS", 59.10, -48.73), ("Low Galactic Latitude", 95.00, -25.00),
    ("Seagull Nebula", 106.23, -10.51),
]
DP1_FIELD_RADIUS_DEG = 0.6   # half-diagonal of a ~1 deg^2 field; coverage test is only approximate


def dp1_field_for(ra, dec):
    """Name of the DP1 field whose nominal area (0.6 deg radius) contains the position, else None."""
    for name, r0, d0 in DP1_FIELDS:
        if angsep_arcsec(ra, dec, r0, d0) / 3600.0 <= DP1_FIELD_RADIUS_DEG:
            return name
    return None


def probe(timeout=20):
    """Probe the LSST endpoints without credentials.  Returns dict of {name: (http_code, bytes)}."""
    import requests
    res = {}
    for name, url, kw in [
        ("tap_availability", TAP_AVAIL, {}),
        ("hips_list_dp1", HIPS_LIST, {}),
        ("sia_dp1", SIA_DP1, {"params": {"POS": "CIRCLE 53.1 -28.1 0.05", "MAXREC": 1}}),
        ("hips_tile_band_r", "https://data.lsst.cloud/api/hips/v2/dp1/deep_coadd/band_r/properties", {}),
    ]:
        try:
            r = requests.get(url, timeout=timeout, **kw)
            res[name] = (r.status_code, len(r.content))
        except Exception as e:
            res[name] = (0, 0)
    return res


def find_local(ra, dec, local_dir, radius_arcmin=2.0):
    """FITS files under local_dir (user-downloaded DP1 data) whose WCS footprint contains the position."""
    from astropy.io import fits
    from astropy.wcs import WCS
    out = []
    for fn in sorted(glob.glob(os.path.join(os.path.expanduser(local_dir), "**", "*.fits*"), recursive=True)):
        try:
            with fits.open(fn) as h:
                for k, hd in enumerate(h):
                    if hd.data is None or hd.data.ndim != 2:
                        continue
                    w = WCS(hd.header)
                    if not w.has_celestial:
                        continue
                    ny, nx = hd.data.shape
                    x, y = w.all_world2pix(ra, dec, 0)
                    if -radius_arcmin * 60 / 0.2 < x < nx + radius_arcmin * 60 / 0.2 and \
                       -radius_arcmin * 60 / 0.2 < y < ny + radius_arcmin * 60 / 0.2:
                        mjd = hd.header.get("MJD-OBS") or hd.header.get("MJD") or h[0].header.get("MJD-OBS")
                        filt = hd.header.get("FILTER") or h[0].header.get("FILTER") or h[0].header.get("BAND") or "?"
                        out.append(dict(obs_id=os.path.basename(fn), obsid=fn, collection="LSST",
                                        instrument="LSSTComCam", filter=str(filt),
                                        texp=float(h[0].header.get("EXPTIME", float("nan")) or float("nan")),
                                        mjd=float(mjd) if mjd else float("nan"), pixscale=0.2,
                                        depth=float("nan"), proposal="DP1", calib_level=2,
                                        target="local", ra=ra, dec=dec, kind="exposure", path=fn))
                        break
        except Exception:
            continue
    return out


def query_candidates(ra, dec, radius_arcmin=1.0, local_dir=None):
    """Return (candidates, status) for the LSST provider.  status is a dict with the reason
    when the provider cannot serve data."""
    status = {"provider": "LSST", "available": False, "reason": "", "dp1_field": dp1_field_for(ra, dec)}
    cands = []
    if local_dir:
        cands = find_local(ra, dec, local_dir, radius_arcmin)
        if cands:
            status.update(available=True, reason="local user-supplied Rubin FITS (%d)" % len(cands))
            return cands, status
    pr = probe()
    status["probe"] = pr
    if status["dp1_field"] is None:
        status["reason"] = "position is outside the seven Rubin DP1 fields (no public or private DP1 data there)"
    elif pr.get("sia_dp1", (0, 0))[0] == 401:
        status["reason"] = ("inside DP1 field '%s' but Rubin Science Platform login is required (HTTP 401); "
                            "no credentials are used. Download DP1 FITS yourself and pass --lsst-local-dir. "
                            "Falling back to MAST/Pan-STARRS." % status["dp1_field"])
    else:
        status["reason"] = "SIA returned %s; not implemented without authentication" % (pr.get("sia_dp1"),)
    return cands, status
