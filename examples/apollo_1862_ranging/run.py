"""(1862) Apollo: twelve Catalina observations to an element distribution.

Reads observations.obs80, builds one tracklet per night with the station
sigma, links the nights, and writes result.json.  The nightlink command-line
tool replaces every tracklet's station with one --obs-code and, when a
tracklet has no sigma, uses 0.3 arcsec.  This script keeps station 703 and
obs.station_sigma (0.50 arcsec).
"""
import json
import os
import sys
import time
import urllib.request

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "2"

import numpy as np

_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _ROOT)

from moving import arc_ranging as AR
from moving import nightlink as N
from moving import obs as O

HERE = os.path.dirname(os.path.abspath(__file__))
OBS_PATH = os.path.join(HERE, "observations.obs80")
OUT_PATH = os.path.join(HERE, "result.json")
STATION = "703"
SEED = 1862
N_SAMPLES = 200
N_PROPAGATE = 8
SBDB_URL = "https://ssd-api.jpl.nasa.gov/sbdb.api?sstr=1862&full-prec=true"
SBDB_NAMES = ("e", "a", "q", "i", "om", "w", "ma")


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, float):
        return float(value) if np.isfinite(value) else None
    return value


def _jd_tdb_to_utc(jd):
    from astropy.time import Time
    return Time(float(jd), format="jd", scale="tdb").utc.isot


def _mjd_to_utc(mjd):
    from astropy.time import Time
    return Time(float(mjd), format="mjd", scale="utc").utc.isot


def _tracklets(parsed):
    by = {}
    for row in parsed:
        by.setdefault(row["stn"], []).append(row)
    trks = []
    made = []
    for stn, rows in by.items():
        rows = sorted(rows, key=lambda o: o["mjd_utc"])
        groups = []
        cur = [rows[0]]
        for row in rows[1:]:
            if row["mjd_utc"] - cur[-1]["mjd_utc"] <= 3.0 / 24.0:
                cur.append(row)
            else:
                groups.append(cur)
                cur = [row]
        groups.append(cur)
        for group in groups:
            if len(group) < 2:
                continue
            year = 2000 + int((group[0]["mjd_utc"] - 51544.5) / 365.25)
            sig = float(O.station_sigma(stn, year))
            mjd = [o["mjd_utc"] for o in group]
            trk = N.Tracklet(
                mjd, [o["ra"] for o in group], [o["dec"] for o in group], sig,
                code=stn, tid=len(trks), night=int(np.floor(group[0]["mjd_utc"] + 0.3)),
            )
            trk.obs()
            trks.append(trk)
            made.append({
                "id": int(trk.id),
                "night": int(trk.night),
                "code": stn,
                "sigma_arcsec": sig,
                "year_for_sigma": int(year),
                "n_obs": len(group),
                "mjd_first": float(mjd[0]),
                "mjd_last": float(mjd[-1]),
                "arc_hours": float((mjd[-1] - mjd[0]) * 24.0),
                "t0": float(trk.t0),
            })
    return trks, made


def _sbdb():
    req = urllib.request.Request(SBDB_URL, headers={"User-Agent": "astrafex-example"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        payload = json.loads(resp.read().decode("utf-8"))
    orbit = payload["orbit"]
    wanted = []
    for element in orbit["elements"]:
        if element.get("name") in SBDB_NAMES:
            wanted.append({
                "name": element.get("name"),
                "value": element.get("value"),
                "sigma": element.get("sigma"),
                "units": element.get("units"),
                "title": element.get("title"),
            })
    epoch = orbit.get("epoch")
    cov = orbit.get("cov_epoch")
    obj = payload.get("object") or {}
    return {
        "url": SBDB_URL,
        "fullname": obj.get("fullname"),
        "orbit_class": (obj.get("orbit_class") or {}).get("name"),
        "neo": obj.get("neo"),
        "pha": obj.get("pha"),
        "epoch_jd_tdb": epoch,
        "epoch_utc": _jd_tdb_to_utc(epoch) if epoch else None,
        "cov_epoch_jd_tdb": cov,
        "cov_epoch_utc": _jd_tdb_to_utc(cov) if cov else None,
        "equinox": orbit.get("equinox"),
        "first_obs": orbit.get("first_obs"),
        "last_obs": orbit.get("last_obs"),
        "soln_date": orbit.get("soln_date"),
        "n_obs_used": orbit.get("n_obs_used"),
        "rms": orbit.get("rms"),
        "condition_code": orbit.get("condition_code"),
        "comment": orbit.get("comment"),
        "elements": wanted,
        "note": (
            "JPL SBDB solution at its own epoch. It is not the epoch of these "
            "twelve observations, and it uses a different observation set and force model. "
            "The elements are not subtracted from the two-body fit or from the Laplace sample."
        ),
    }


def main():
    text = open(OBS_PATH).read()
    parsed = O.parse_obs80(text)
    parsed = [o for o in parsed if o["stn"] == STATION and o["note2"] not in ("S", "s")]
    if len(parsed) != 12:
        raise SystemExit("expected 12 station-703 observations, got %d" % len(parsed))
    code = O.obscodes()[STATION]
    trks, made = _tracklets(parsed)
    t0 = time.perf_counter()
    linked = N.link_nights(trks, max_gap_days=12, min_tracklets=3)
    link_seconds = time.perf_counter() - t0
    if len(linked["groups"]) != 1 or linked["pairs"] or linked["unlinked"]:
        raise SystemExit("expected one group and nothing else: %s" % (
            {k: linked[k] if k != "groups" else len(linked["groups"]) for k in ("groups", "pairs", "unlinked")}))
    group = linked["groups"][0]
    point_epoch = float(group["fit"].t_ref)
    point = {
        "ids": [int(i) for i in group["ids"]],
        "nights": [int(n) for n in group["nights"]],
        "n_obs": int(group["n_obs"]),
        "chi2_red": float(group["chi2_red"]),
        "rms_arcsec": float(group["rms_arcsec"]),
        "max_resid": float(group["max_resid"]),
        "elements": _plain(group["elements"]),
        "mjd_epoch": point_epoch,
        "epoch_utc": _mjd_to_utc(point_epoch),
    }
    before = (point["elements"]["a"], point["chi2_red"], point["rms_arcsec"])
    for item in linked["groups"] + linked["pairs"]:
        item.pop("fit", None)
    t1 = time.perf_counter()
    ranged = AR.attach(
        trks, linked, n_samples=N_SAMPLES, seed=SEED, n_propagate=N_PROPAGATE, timeout=300)
    range_seconds = time.perf_counter() - t1
    summary = ranged["groups"][0]["ranging"]
    after = (
        float(ranged["groups"][0]["elements"]["a"]),
        float(ranged["groups"][0]["chi2_red"]),
        float(ranged["groups"][0]["rms_arcsec"]),
    )
    if after != before:
        raise SystemExit("ranging changed the two-body point fit")
    try:
        sbdb = _sbdb()
        sbdb_error = None
    except Exception as exc:
        sbdb = None
        sbdb_error = "%s: %s" % (type(exc).__name__, exc)
    mjd = [float(o["mjd_utc"]) for o in parsed]
    doc = {
        "designation": "1862",
        "fullname_short": "1862 Apollo",
        "observations_file": "observations.obs80",
        "mpc_get_obs": "https://data.minorplanetcenter.net/api/get-obs",
        "selection": (
            "Public MPC OBS80 for designation 1862, station 703, MJD 59250 to 59260, "
            "note2 not S or s. Twelve lines. The rest of the archive is not in this directory."
        ),
        "station": {
            "code": STATION,
            "lon_deg": float(code[0]),
            "rho_cos": float(code[1]),
            "rho_sin": float(code[2]),
            "name": code[3],
            "sigma_arcsec": float(O.station_sigma(STATION, 2021)),
            "sigma_source": "moving/obs.py DEFAULT_SIGMA['703']",
        },
        "n_obs": len(parsed),
        "mjd_first": min(mjd),
        "mjd_last": max(mjd),
        "span_days": max(mjd) - min(mjd),
        "first_utc": _mjd_to_utc(min(mjd)),
        "last_utc": _mjd_to_utc(max(mjd)),
        "observations": [
            {
                "mjd_utc": float(o["mjd_utc"]),
                "ra_deg": float(o["ra"]),
                "dec_deg": float(o["dec"]),
                "mag": None if o["mag"] is None else float(o["mag"]),
                "band": o["band"] or None,
                "stn": o["stn"],
            }
            for o in parsed
        ],
        "tracklets": made,
        "link": {
            "max_gap_days": 12,
            "min_tracklets": 3,
            "n_groups": len(ranged["groups"]),
            "n_pairs": len(ranged["pairs"]),
            "unlinked": list(ranged["unlinked"]),
            "seconds": link_seconds,
            "point_fit": point,
            "note": (
                "Two-body Sun-only fit from moving.nightlink. The elements are osculating "
                "at mjd_epoch, which is the epoch of that fit, not the mean observation time "
                "and not the JPL epoch."
            ),
        },
        "ranging": _plain(summary),
        "ranging_run": {
            "n_samples": N_SAMPLES,
            "seed": SEED,
            "distribution_seed": SEED,
            "propagation_seed": SEED + 17,
            "n_propagate": N_PROPAGATE,
            "seconds": range_seconds,
            "threads": 2,
        },
        "jpl_sbdb": sbdb,
        "jpl_sbdb_error": sbdb_error,
    }
    with open(OUT_PATH, "w") as handle:
        json.dump(doc, handle, indent=2)
        handle.write("\n")
    qa = summary["quantiles_16_50_84"]["a"]
    print(
        "APOLLO sampler=%s status=%s neff=%.2f rms=%.4f a=%.5f/%.5f/%.5f best_a=%.5f "
        "point_a=%.5f class=%s codes=%s seconds=%.1f"
        % (
            summary.get("sampler"), summary.get("status"), float(summary.get("neff")),
            float(summary.get("best_rms_arcsec")), qa[0], qa[1], qa[2],
            float(summary["best"]["a"]), point["elements"]["a"],
            summary.get("class_probs"), (summary.get("codes") or {}).get("status"),
            link_seconds + range_seconds,
        )
    )
    width = float(qa[2]) - float(qa[0])
    if summary.get("sampler") != "laplace" or summary.get("status") != "ok" or float(summary.get("neff")) < 10 or width > 1.0:
        raise SystemExit("distribution is not the local Laplace sample of this arc")


if __name__ == "__main__":
    main()
