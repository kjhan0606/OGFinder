"""Astrometric observations: MPC 80-column parsing/writing, observer positions
(ground stations via MPC parallax constants + Earth rotation; spacecraft via Horizons or
MPC 's' vector lines), Eggl et al. (2020) catalog debiasing (JPL debias files), weights."""
import os, re
import numpy as np
from .util import EPHEM_DIR, log, utc_mjd_to_tdb_jd, ensure_dir

RE_KM = 6378.137
AU_KM = 149597870.7

# ---------------------------------------------------------------- MPC 80 column
def _parse_ra(s):
    h, m, sec = s.split()
    return 15.0 * (float(h) + float(m) / 60 + float(sec) / 3600)


def _parse_dec(s):
    s = s.strip()
    sign = -1.0 if s[0] == "-" else 1.0
    d, m, sec = s[1:].split()
    return sign * (float(d) + float(m) / 60 + float(sec) / 3600)


def _mpc_date_to_mjd(s):
    yr = int(s[0:4]); mo = int(s[5:7]); dd = float(s[8:])
    from astropy.time import Time
    t = Time({"year": yr, "month": mo, "day": 1}, scale="utc")
    return t.mjd + dd - 1.0


def _unpack_designation(s):
    return s.strip()


def parse_obs80(text):
    """Parse MPC 80-column optical observations (and 's' satellite vector lines).
    Returns list of dicts: desig, mjd_utc, ra, dec [deg], mag, band, cat, stn, sat_xyz_km or None, line."""
    obs = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        L = lines[i]
        if len(L) < 80:
            L = L.ljust(80)
        note2 = L[14]
        if note2 in ("s", "r", "v", "R", "V"):   # continuation lines handled with their parent
            i += 1; continue
        if not re.match(r"^\s*\d{4} \d{2} \d{2}", L[15:32]) and not re.match(r"^\d{4} \d{2} [\d.]+", L[15:32].strip()):
            i += 1; continue
        try:
            mjd = _mpc_date_to_mjd(L[15:32].strip())
            ra = _parse_ra(L[32:44]); dec = _parse_dec(L[44:56])
        except Exception:
            i += 1; continue
        mag = None; band = L[70].strip()
        try:
            mag = float(L[65:70])
        except ValueError:
            pass
        d = dict(desig=L[0:12].strip(), mjd_utc=mjd, ra=ra, dec=dec, mag=mag, band=band,
                 cat=L[71], stn=L[77:80].strip(), note1=L[13], note2=note2, sat_xyz_km=None, line=L.rstrip())
        if note2 == "S" and i + 1 < len(lines) and len(lines[i + 1]) >= 80 and lines[i + 1][14] == "s":
            S = lines[i + 1]
            try:
                unit = S[32]
                x = float(S[34:46]); y = float(S[46:58]); z = float(S[58:70])
                f = AU_KM if unit == "2" else 1.0
                d["sat_xyz_km"] = np.array([x, y, z]) * f
                d["sat_line"] = S.rstrip()
            except ValueError:
                pass
        obs.append(d)
        i += 1
    return obs


def format_ra(ra_deg):
    """HH MM SS.dd (11 chars, MPC 80-col precision 0.01 s)"""
    h = (ra_deg % 360) / 15.0
    tot = round(h * 3600.0, 2)
    hh = int(tot // 3600); mm = int((tot - hh * 3600) // 60); ss = tot - hh * 3600 - mm * 60
    return "%02d %02d %05.2f" % (hh % 24, mm, ss)


def format_dec(dec_deg):
    """+DD MM SS.d (10 chars, 0.1 arcsec)"""
    sg = "-" if dec_deg < 0 else "+"
    tot = round(abs(dec_deg) * 3600.0, 1)
    d = int(tot // 3600); mm = int((tot - d * 3600) // 60); ss = tot - d * 3600 - mm * 60
    return "%s%02d %02d %04.1f" % (sg, d, mm, ss)


def format_date(mjd_utc):
    from astropy.time import Time
    t = Time(mjd_utc, format="mjd", scale="utc")
    y, mo, d = t.ymdhms.year, t.ymdhms.month, t.ymdhms.day
    frac = (mjd_utc - np.floor(mjd_utc))
    return "%04d %02d %08.5f" % (y, mo, d + frac)


def pack_designation(des):
    """Packed 12-column designation, only for numbered asteroids (5 digits, right-aligned in cols 1-5);
    others are written as 7-character provisional if the user supplies them already packed."""
    des = str(des).strip()
    if des.isdigit() and len(des) <= 5:
        return des.rjust(5) + " " * 7
    return des.ljust(12)[:12]


def write_obs80(observations, desig="", stn="250", cat="", band="", note="C"):
    """MPC1992 80-column optical lines. observations: dicts mjd_utc, ra, dec [deg J2000], optional
    mag, band, cat, stn, desig, sat_xyz_km (-> an 's' line, geocentric km, unit code 2? no: km, code 1)."""
    out = []
    for o in observations:
        mag = o.get("mag")
        line = pack_designation(o.get("desig", desig))      # cols 1-12
        line = line.ljust(14) + note[:1]                    # col 15 = note 2 (C = CCD)
        line = line.ljust(15) + format_date(o["mjd_utc"]).ljust(17)   # cols 16-32
        line = line.ljust(32) + format_ra(o["ra"]) + " " + format_dec(o["dec"])   # 33-44, 45-55
        line = line.ljust(65) + (("%4.1f " % mag) if (mag is not None and np.isfinite(mag)) else "     ")
        line = line.ljust(70) + ((o.get("band", band) or " ")[:1])
        line = line.ljust(71) + ((o.get("cat", cat) or " ")[:1])
        line = line.ljust(77) + o.get("stn", stn).ljust(3)[:3]
        out.append(line)
        if o.get("sat_xyz_km") is not None:
            x, y, z = o["sat_xyz_km"]
            s = pack_designation(o.get("desig", desig)).ljust(14) + "s"
            s = s.ljust(15) + format_date(o["mjd_utc"]).ljust(17)
            s = s.ljust(32) + "1 " + ("%+11.4f %+11.4f %+11.4f" % (x, y, z))   # km, cols 33-34 flag
            out.append(s.ljust(77) + o.get("stn", stn).ljust(3)[:3])
    return "\n".join(out) + "\n"


# ------------------------------------------------------------ station table
_OBSCODES = {}


def obscodes():
    """MPC observatory codes: code -> (lon_deg, rho_cos, rho_sin, name).  Spacecraft have None."""
    if _OBSCODES:
        return _OBSCODES
    fn = os.path.join(EPHEM_DIR, "ObsCodes.html")
    if not os.path.exists(fn):
        import requests
        ensure_dir(EPHEM_DIR)
        r = requests.get("https://www.minorplanetcenter.net/iau/lists/ObsCodes.html", timeout=60)
        r.raise_for_status()
        open(fn, "w").write(r.text)
    for line in open(fn, errors="replace"):
        if len(line) < 30 or line.startswith("<") or line.startswith("Code"):
            continue
        code = line[0:3]
        try:
            lon = float(line[4:13]); c = float(line[13:21]); s = float(line[21:30])
            _OBSCODES[code] = (lon, c, s, line[30:].strip())
        except ValueError:
            _OBSCODES[code] = (None, None, None, line[30:].strip())
    return _OBSCODES


def _c2t_matrix(mjd_utc_arr):
    """GCRS->ITRS rotation matrices (n,3,3) from erfa c2t06a (IAU 2006/2000A)."""
    import erfa
    from astropy.time import Time
    from astropy.utils import iers
    iers.conf.auto_download = False
    iers.conf.iers_degraded_accuracy = "ignore"
    t = Time(np.atleast_1d(mjd_utc_arr), format="mjd", scale="utc")
    try:
        dut1 = t.delta_ut1_utc
        xp, yp = iers.earth_orientation_table.get().pm_xy(t)
        xp = xp.to_value("rad"); yp = yp.to_value("rad")
    except Exception:
        t.delta_ut1_utc = np.zeros(len(t)); dut1 = t.delta_ut1_utc; xp = yp = np.zeros(len(t))
    tt = t.tt; ut1 = t.ut1
    M = erfa.c2t06a(tt.jd1, tt.jd2, ut1.jd1, ut1.jd2, xp, yp)
    return M


def geocentric_gcrs_km(code, mjd_utc, sat_xyz_km=None, hz_fallback=True):
    """Geocentric GCRS (~ICRF) position of the station in km for UTC MJD array.
    Ground: MPC parallax constants rotated with erfa.  Spacecraft: sat_xyz_km, else Horizons."""
    mjd_utc = np.atleast_1d(np.asarray(mjd_utc, float))
    oc = obscodes().get(code)
    if code == "500":
        return np.zeros((len(mjd_utc), 3))
    if oc is not None and oc[0] is not None:
        lon, c, s, _ = oc
        itrs = RE_KM * np.array([c * np.cos(np.radians(lon)), c * np.sin(np.radians(lon)), s])
        M = _c2t_matrix(mjd_utc)             # GCRS -> ITRS
        return np.einsum("nji,j->ni", M, itrs)   # M^T @ itrs
    if sat_xyz_km is not None:
        return np.tile(np.asarray(sat_xyz_km, float), (len(mjd_utc), 1))
    sc = {"250": "-48", "274": "-170"}.get(code)
    if sc and hz_fallback:
        from . import horizons
        jd = utc_mjd_to_tdb_jd(mjd_utc)
        return horizons.vectors(sc, jd, center="500@399")[:, :3] * AU_KM
    raise ValueError("no position source for observatory %s" % code)


def spacecraft_geocentric_km(body, mjd_utc):
    from . import horizons
    jd = utc_mjd_to_tdb_jd(mjd_utc)
    return horizons.vectors(body, jd, center="500@399")[:, :3] * AU_KM


# ------------------------------------------------------------ debiasing (Eggl et al. 2020)
_DEBIAS = {}


def _load_debias():
    if _DEBIAS:
        return _DEBIAS
    d = os.path.join(EPHEM_DIR, "debias")
    fb = os.path.join(d, "bias.dat"); ft = os.path.join(d, "tiles.dat")
    if not (os.path.exists(fb) and os.path.exists(ft)):
        return None
    codes = None
    rows = []
    for line in open(fb):
        if line.startswith("!"):
            if codes is None and line.strip().startswith("! a b c"):
                codes = line.strip()[1:].split()
            continue
        if line.strip():
            rows.append([float(x) for x in line.split()])
    tiles = np.loadtxt(ft, comments="!")
    _DEBIAS["codes"] = codes
    _DEBIAS["bias"] = np.array(rows)
    _DEBIAS["tiles"] = tiles
    ra, de = np.radians(tiles[:, -2]), np.radians(tiles[:, -1])
    from scipy.spatial import cKDTree
    _DEBIAS["xyz"] = np.stack([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)], 1)
    _DEBIAS["tree"] = cKDTree(_DEBIAS["xyz"])
    return _DEBIAS


def debias(ra_deg, dec_deg, mjd_utc, catcode):
    """Return (dra_cosdec_arcsec, ddec_arcsec) to SUBTRACT from the observation (Eggl et al. 2020,
    JPL debias files, Gaia-DR2-referenced).  Zero if the catalog is not tabulated."""
    D = _load_debias()
    if D is None or not catcode or catcode.strip() == "" or catcode not in D["codes"]:
        return 0.0, 0.0
    k = D["codes"].index(catcode)
    ra, de = np.radians(ra_deg), np.radians(dec_deg)
    u = np.array([np.cos(de) * np.cos(ra), np.cos(de) * np.sin(ra), np.sin(de)])
    _, ti = D["tree"].query(u)
    dra, dde, pmra, pmde = D["bias"][ti, 4 * k:4 * k + 4]
    yr = 2000.0 + (mjd_utc - 51544.5) / 365.25
    return dra + (yr - 2000.0) * pmra / 1000.0, dde + (yr - 2000.0) * pmde / 1000.0


# ------------------------------------------------------------ weights
# NOT the Vereš et al. (2017) table (not reproduced here).  Defaults are coarse, user-adjustable
# one-sigma values in arcsec; pass --weights-file CSV (stn,sigma_arcsec) to override.
DEFAULT_SIGMA = {"F51": 0.20, "F52": 0.20, "T08": 0.30, "T05": 0.30, "G96": 0.30, "703": 0.50, "E12": 0.50,
                 "704": 0.60, "644": 0.60, "691": 0.50, "699": 0.50, "C57": 0.50, "W68": 0.50, "250": 0.03,
                 "274": 0.03, "568": 0.20, "M22": 0.30, "V00": 0.30, "G45": 0.30, "D29": 0.50}


def station_sigma(stn, year, table=None):
    t = table or DEFAULT_SIGMA
    if stn in t:
        return t[stn]
    return 1.0 if year > 1995 else 2.0 if year > 1950 else 3.0


def load_weights_file(path):
    t = {}
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        a = re.split(r"[,\s]+", line.strip())
        t[a[0]] = float(a[1])
    return t
