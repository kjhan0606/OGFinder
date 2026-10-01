"""IMCCE SkyBoT cone search (VO conesearch, text output) with caching, also used for the
expected-position tests and the known-object cross-match."""
import numpy as np
from .util import http_get, utc_mjd_to_tdb_jd

URL = "http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php"


def cone(ra, dec, radius_deg, mjd_utc, observer="500", mime="text"):
    """Known solar-system bodies in a cone at UTC MJD.  Returns list of dicts
    (num, name, ra, dec, cls, vmag, err_arcsec, dra_arcsec_h, ddec_arcsec_h, delta_au, r_au)."""
    jd = float(mjd_utc) + 2400000.5
    txt = http_get(URL, params={"-ra": "%.6f" % ra, "-dec": "%.6f" % dec, "-rd": "%.5f" % radius_deg,
                                "-ep": "%.6f" % jd, "-mime": mime, "-output": "basic", "-loc": observer,
                                "-observer": observer}, timeout=90)
    out = []
    for line in txt.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        c = [x.strip() for x in line.split("|")]
        if len(c) < 12:
            continue
        try:
            h, m, s = c[2].split(); ra_d = 15.0 * (float(h) + float(m) / 60 + float(s) / 3600)
            sg = -1 if c[3].strip().startswith("-") else 1
            dd, mm, ss = c[3].replace("+", "").replace("-", "").split()
            de_d = sg * (float(dd) + float(mm) / 60 + float(ss) / 3600)
            out.append(dict(num=c[0] if c[0] != "-" else "", name=c[1], ra=ra_d, dec=de_d, cls=c[4],
                            vmag=float(c[5]) if c[5] not in ("", "-") else np.nan, err=float(c[6]),
                            dra=float(c[8]), ddec=float(c[9]), delta=float(c[10]), r=float(c[11])))
        except Exception:
            continue
    return out
