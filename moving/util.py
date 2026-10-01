"""Shared helpers: caching HTTP, JSON/TSV I/O, hashing, time conversion."""
import os, sys, json, hashlib, time, io
import numpy as np

CACHE_DIR = os.path.expanduser(os.environ.get("OGF_MOVING_CACHE", "~/.ds9/moving_cache"))
MAST_CACHE = os.path.expanduser(os.environ.get("OGF_MAST_CACHE", "~/.ds9/mast_cache"))
EPHEM_DIR = os.path.expanduser(os.environ.get("OGF_EPHEM_DIR", "~/.ds9/ephem"))

C_KMS = 299792.458
AU_KM = 149597870.7
C_AUD = C_KMS * 86400.0 / AU_KM          # speed of light, AU/day
ARCSEC = np.pi / 180.0 / 3600.0


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def sha256_file(path, block=1 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(block)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def ensure_dir(p):
    os.makedirs(p, exist_ok=True)
    return p


def http_get(url, params=None, cache=True, timeout=90, retries=3, headers=None,
             method="GET", data=None, as_json=False, cache_tag=""):
    """GET/POST with on-disk cache (keyed on url+params+data).  Returns text or json."""
    import requests
    key = hashlib.sha1(json.dumps([url, params, data, method, cache_tag], sort_keys=True,
                                  default=str).encode()).hexdigest()
    ensure_dir(CACHE_DIR)
    cp = os.path.join(CACHE_DIR, key + ".txt")
    if cache and os.path.exists(cp) and os.path.getsize(cp) > 0:
        txt = open(cp, encoding="utf-8").read()
        return json.loads(txt) if as_json else txt
    err = None
    for k in range(retries):
        try:
            if method == "GET":
                r = requests.get(url, params=params, timeout=timeout, headers=headers)
            else:
                r = requests.post(url, params=params, data=data, timeout=timeout, headers=headers)
            if r.status_code == 200:
                txt = r.text
                if cache:
                    with open(cp, "w", encoding="utf-8") as f:
                        f.write(txt)
                return json.loads(txt) if as_json else txt
            err = "HTTP %d: %s" % (r.status_code, r.text[:200])
        except Exception as e:  # network
            err = repr(e)
        time.sleep(1.5 * (k + 1))
    raise RuntimeError("request failed: %s %s: %s" % (url, params, err))


def write_json(path, obj):
    def conv(o):
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, default=conv)


def read_json(path):
    with open(path) as f:
        return json.load(f)


def write_tsv(path, rows, cols=None):
    """rows: list of dicts.  Writes tab-separated with header; returns cols."""
    if cols is None:
        cols = []
        for r in rows:
            for k in r:
                if k not in cols:
                    cols.append(k)
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            out = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, (float, np.floating)):
                    v = "%.10g" % v
                elif isinstance(v, (list, tuple, np.ndarray)):
                    v = ",".join(str(x) for x in v)
                out.append(str(v))
            f.write("\t".join(out) + "\n")
    return cols


def read_tsv(path):
    rows = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            v = line.rstrip("\n").split("\t")
            d = {}
            for k, x in zip(hdr, v):
                d[k] = _auto(x)
            rows.append(d)
    return rows


def _auto(x):
    if x == "":
        return ""
    try:
        return int(x)
    except ValueError:
        pass
    try:
        return float(x)
    except ValueError:
        return x


def versions():
    out = {"python": sys.version.split()[0]}
    for m in ("numpy", "scipy", "astropy", "sep", "reproject", "rebound", "assist",
              "astroquery", "emcee", "sklearn", "jplephem", "requests"):
        try:
            mod = __import__(m)
            out[m] = getattr(mod, "__version__", "?")
        except Exception:
            out[m] = None
    return out


def mjd_to_jd(m):
    return np.asarray(m, dtype=float) + 2400000.5


def utc_mjd_to_tdb_jd(mjd_utc):
    """UTC MJD -> TDB JD (geocentric TDB-TT approximation from erfa via astropy)."""
    from astropy.time import Time
    t = Time(np.atleast_1d(np.asarray(mjd_utc, dtype=float)), format="mjd", scale="utc")
    return t.tdb.jd


def tdb_jd_to_utc_mjd(jd_tdb):
    from astropy.time import Time
    t = Time(np.atleast_1d(np.asarray(jd_tdb, dtype=float)), format="jd", scale="tdb")
    return t.utc.mjd


def angsep_arcsec(ra1, de1, ra2, de2):
    """Great-circle separation (deg inputs) in arcsec."""
    ra1, de1, ra2, de2 = [np.radians(np.asarray(x, dtype=float)) for x in (ra1, de1, ra2, de2)]
    s = np.sin((de2 - de1) / 2) ** 2 + np.cos(de1) * np.cos(de2) * np.sin((ra2 - ra1) / 2) ** 2
    return np.degrees(2 * np.arcsin(np.sqrt(np.clip(s, 0, 1)))) * 3600.0


def tangent_offsets(ra, de, ra0, de0):
    """Gnomonic offsets (xi, eta) in arcsec of (ra,de) relative to (ra0,de0), degrees in."""
    ra, de, ra0, de0 = [np.radians(np.asarray(x, dtype=float)) for x in (ra, de, ra0, de0)]
    cosc = np.sin(de0) * np.sin(de) + np.cos(de0) * np.cos(de) * np.cos(ra - ra0)
    xi = np.cos(de) * np.sin(ra - ra0) / cosc
    eta = (np.cos(de0) * np.sin(de) - np.sin(de0) * np.cos(de) * np.cos(ra - ra0)) / cosc
    return np.degrees(xi) * 3600.0, np.degrees(eta) * 3600.0
