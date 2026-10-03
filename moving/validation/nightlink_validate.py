#!/usr/bin/env python3
"""Validation of moving.nightlink (cross-night orbit-fit linking).

  --exp synth     2-body truth, synthetic noise, decoys (SYNTHETIC; dynamics model exact -> optimistic)
  --exp horizons  REAL asteroid orbits: n-body ephemerides from JPL Horizons for objects selected in one sky window (from SBDB
                  elements), synthetic Gaussian noise + decoys made from shifted copies of real tracklets (SEMI-SYNTHETIC)
  --exp mpc       REAL astrometry from the MPC (get-obs) of the same objects: per-station per-night tracklets, topocentric parallax,
                  real astrometric errors; truth = MPC designation (REAL data; objects pre-selected as known, bright, in one window)
Outputs JSON + markdown in --out.  Network needed for horizons/mpc (cached in ~/.cache of moving.util).
"""
import argparse, json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from moving import nightlink as N, nightlink_sim as S, kepler as K, horizons as HZ, mpc as MPC, obs as O
from moving.util import http_get


def synth(args):
    rows = []
    cases = []
    for sig in (0.1, 0.3, 1.0):
        for nights in ([0, 1, 3], [0, 2, 6], [0, 1, 7, 14], [0, 3, 10], [0, 1], [0, 7]):
            cases.append((sig, nights))
    for sig, nights in cases:
        agg = dict(n_obj=0, rec=0, pure=0, wrong=0, groups=0, complete=0, t=0.0, tracklets=0, pair_pure=0, pair_wrong=0)
        for seed in range(args.seeds):
            rng = np.random.default_rng(1000 + seed)
            trks = S.make_field(args.nobj, nights, rng, sig=sig, n_decoy=args.ndecoy, p_detect=0.9)
            t0 = time.time(); res = N.link_nights(trks); dt = time.time() - t0
            sc = S.score(trks, res)
            for k, v in (("n_obj", sc["n_objects"]), ("rec", sc["recovered"]), ("pure", sc["pure"]), ("wrong", sc["wrong"]),
                         ("groups", sc["n_groups"]), ("complete", sc["complete"]), ("pair_pure", sc["pair_pure"]), ("pair_wrong", sc["pair_wrong"]), ("t", dt), ("tracklets", len(trks))):
                agg[k] += v
        rows.append(dict(sigma=sig, nights=nights, **agg))
        print("synth sigma=%.1f nights=%s recall(>=3 tracklets)=%.3f wrong=%d / %d groups; pairs pure %d wrong %d  (%.0f s/field, %d tracklets)" % (
            sig, nights, agg["rec"] / max(agg["n_obj"], 1), agg["wrong"], agg["groups"], agg["pair_pure"], agg["pair_wrong"], agg["t"] / args.seeds, agg["tracklets"] // args.seeds), flush=True)
    return rows


# ------------------------------------------------------------------------------------------------ real objects
def sbdb_table(path, hmax=15.5):
    if os.path.isfile(path):
        return json.load(open(path))
    import requests
    r = requests.get("https://ssd-api.jpl.nasa.gov/sbdb_query.api", timeout=300, params={
        "fields": "pdes,a,e,i,om,w,ma,epoch,H", "sb-kind": "a", "sb-ns": "n", "sb-cdata": json.dumps({"AND": ["H|RG|0|%g" % hmax, "a|RG|1.8|4.5"]})})
    r.raise_for_status(); rows = r.json()["data"]
    json.dump(rows, open(path, "w")); return rows


def pick_window(rows, mjd0, half_deg=4.0, nmax=70):
    """Densest window of ``2*half_deg`` square at mjd0 among 2-body-propagated SBDB objects (position error of a few arcmin)."""
    a = np.array([float(x[1]) for x in rows]); e = np.array([float(x[2]) for x in rows]); inc = np.array([float(x[3]) for x in rows])
    om = np.array([float(x[4]) for x in rows]); w = np.array([float(x[5]) for x in rows]); ma = np.array([float(x[6]) for x in rows])
    ep = np.array([float(x[7]) for x in rows]) - 2400000.5
    n_mean = np.degrees(np.sqrt(K.MU / a ** 3))
    M = (ma + n_mean * (mjd0 - ep)) % 360.0
    st = [K.elements_to_state(a[k], e[k], inc[k], om[k], w[k], M[k], ecliptic=False) for k in range(len(a))]
    r = np.array([s[0] for s in st]); v = np.array([s[1] for s in st])
    op, _ = N.observer_helio(mjd0)
    ra, de, _ = N.predict(r, v, mjd0, mjd0, np.broadcast_to(op[0], r.shape))
    best = (0, None)
    for rc in np.arange(0, 360, 2.0):
        for dc in np.arange(-20, 21, 2.0):
            dx = ((ra - rc + 180) % 360 - 180) * np.cos(np.radians(dc))
            m = (np.abs(dx) < half_deg) & (np.abs(de - dc) < half_deg)
            if m.sum() > best[0]:
                best = (int(m.sum()), (rc, dc, np.where(m)[0]))
    n, (rc, dc, idx) = best
    idx = idx[np.argsort([float(rows[k][8]) for k in idx])][:nmax]
    return rc, dc, [rows[k][0] for k in idx], n


def horizons_tracklets(des_list, nights, t0, rng, sig, n_per_night=3, spacing_h=0.5):
    trks = []; rate_pool = []
    jd0 = t0 + 2400000.5
    for oi, des in enumerate(des_list):
        tl = []
        for nt in nights:
            tl += list(jd0 + nt + np.arange(n_per_night) * spacing_h / 24.0)
        try:
            ot = HZ.observer_table(des + ";", tl, center="500")
        except Exception as ex:
            print("horizons failed", des, ex); continue
        if not np.all(np.isfinite(ot[:, 0])):
            continue
        for k, nt in enumerate(nights):
            sl = slice(k * n_per_night, (k + 1) * n_per_night)
            times = np.array(tl[sl]) - 2400000.5
            ra, de = ot[sl, 0].copy(), ot[sl, 1].copy()
            cd = np.cos(np.radians(de))
            ra = ra + rng.normal(0, sig, len(ra)) / 3600 / cd; de = de + rng.normal(0, sig, len(de)) / 3600
            trks.append(N.Tracklet(times, ra, de, sig, tid=len(trks), night=int(nt), truth=oi))
    return trks


def add_shifted_decoys(trks, rng, nights, per_night, half_deg, rc, dc, mjd0):
    """Decoys = real tracklets copied from other epochs with random sky shifts inside the window (keeps realistic rates)."""
    real = [t for t in trks if t.truth is not None and t.truth >= 0]
    out = list(trks)
    for nt in nights:
        for _ in range(per_night):
            src = real[rng.integers(len(real))]
            sh_ra = rng.uniform(-half_deg, half_deg) / np.cos(np.radians(dc)); sh_de = rng.uniform(-half_deg, half_deg)
            ra = src.ra - src.ra0 + rc + sh_ra
            de = src.dec - src.dec0 + dc + sh_de
            times = mjd0 + nt + (src.mjd - src.mjd.min())
            out.append(N.Tracklet(times, ra % 360, de, src.sig, tid=len(out), night=int(nt), truth=-1))
    return out


def horizons_exp(args):
    os.makedirs(args.out, exist_ok=True)
    mjd0 = args.mjd0
    rows = sbdb_table(os.path.join(args.out, "sbdb_H155.json"))
    rc, dc, des, ncount = pick_window(rows, mjd0)
    print("window centre RA %.1f Dec %.1f: %d candidates, using %d" % (rc, dc, ncount, len(des)), flush=True)
    out = []
    for sig in args.hz_sigmas:
        for nights in args.hz_nights:
            rng = np.random.default_rng(7)
            trks = horizons_tracklets(des, nights, mjd0, rng, sig)
            ndec = len(des)
            trks = add_shifted_decoys(trks, rng, nights, ndec, 4.0, rc, dc, mjd0)
            t0 = time.time(); res = N.link_nights(trks, floor0=args.floor0, floor_rate=args.floor_rate); dt = time.time() - t0
            sc = S.score(trks, res)
            sc.update(sigma=sig, nights=nights, tracklets=len(trks), seconds=dt, n_obj_total=len(des))
            sc["wrong_groups"] = [dict(truth=g["truth"], nights=g["nights"], chi2=round(g["chi2_red"], 2), n_trk=len(g["ids"])) for g in res["groups"] + res["pairs"]
                                  if not (len(set(g["truth"])) == 1 and g["truth"][0] >= 0)]
            out.append(sc)
            print("horizons sigma=%.1f nights=%s objects=%d recall=%.3f pure=%d wrong=%d complete=%d; pairs pure %d wrong %d (%d tracklets, %.0f s)" % (
                sig, nights, sc["n_objects"], sc["recall"], sc["pure"], sc["wrong"], sc["complete"], sc["pair_pure"], sc["pair_wrong"], len(trks), dt), flush=True)
    return out, dict(rc=rc, dc=dc, designations=des)


def mpc_tracklets(des_list, mjd_lo, mjd_hi, max_gap_h=3.0, min_obs=2):
    trks = []; skipped = 0
    for oi, des in enumerate(des_list):
        try:
            j = MPC.get_obs(des)
        except Exception as ex:
            print("mpc failed", des, ex); continue
        ob = O.parse_obs80(j.get("OBS80", ""))
        ob = [o for o in ob if mjd_lo <= o["mjd_utc"] <= mjd_hi and o["note2"] not in ("S", "s") and o["stn"] not in ("247", "270", "")]
        by = {}
        for o in ob:
            by.setdefault(o["stn"], []).append(o)
        for stn, L in by.items():
            L.sort(key=lambda o: o["mjd_utc"])
            cur = [L[0]]
            groups = []
            for o in L[1:]:
                if o["mjd_utc"] - cur[-1]["mjd_utc"] <= max_gap_h / 24.0:
                    cur.append(o)
                else:
                    groups.append(cur); cur = [o]
            groups.append(cur)
            for g in groups:
                if len(g) < min_obs:
                    skipped += 1; continue
                yr = 2000 + int((g[0]["mjd_utc"] - 51544.5) / 365.25)
                sig = O.station_sigma(stn, yr)
                try:
                    trks.append(N.Tracklet([o["mjd_utc"] for o in g], [o["ra"] for o in g], [o["dec"] for o in g], sig, code=stn,
                                           tid=len(trks), night=int(np.floor(g[0]["mjd_utc"] + 0.3)), truth=oi))
                    trks[-1].obs()                      # station known?
                except Exception:
                    trks.pop(); skipped += 1
    return trks, skipped


def mpc_objects(cands, mjd_lo, mjd_hi, want, seed=3):
    """Random window objects (faint end preferred is NOT applied) that really have >= 2 station-nights with >= 2 MPC observations."""
    rng = np.random.default_rng(seed)
    keep = []
    for k in rng.permutation(len(cands)):
        des = cands[k]
        try:
            trks, _ = mpc_tracklets([des], mjd_lo, mjd_hi)
        except Exception:
            continue
        if len({(t.night, t.code) for t in trks}) >= 2:
            keep.append(des)
        if len(keep) >= want:
            break
    return keep


def mpc_exp(args):
    os.makedirs(args.out, exist_ok=True)
    mjd0 = args.mjd0
    rows = sbdb_table(os.path.join(args.out, "sbdb_H155.json"))
    rc, dc, cands, ncount = pick_window(rows, mjd0, nmax=100000)
    des = mpc_objects(cands, mjd0 - 1.0, mjd0 + max(args.spans), args.nmpc)
    print("mpc: window RA %.1f Dec %.1f: %d candidates, %d with >= 2 station-nights of MPC astrometry" % (rc, dc, len(cands), len(des)), flush=True)
    out = []
    for span in args.spans:
        trks, skipped = mpc_tracklets(des, mjd0 - 1.0, mjd0 + span)
        nobj = len({t.truth for t in trks})
        nights = sorted({t.night for t in trks})
        cnt = {}
        for t in trks:
            cnt.setdefault(t.truth, set()).add(t.night)
        n_multi = sum(1 for v in cnt.values() if len(v) >= 2)
        t0 = time.time(); res = N.link_nights(trks, floor0=args.floor0, floor_rate=args.floor_rate, max_gap_days=span + 2); dt = time.time() - t0
        sc = S.score(trks, res)
        # nights-aware score: truth links possible only for objects with tracklets on >= 2 nights
        sc.update(span_days=span, tracklets=len(trks), objects_seen=nobj, objects_multi_night=n_multi, nights=len(nights), seconds=dt, skipped_short=skipped)
        out.append(sc)
        print("mpc span=%d d: %d tracklets of %d objects (%d multi-night, %d nights); groups(>=3): recall=%.3f pure=%d wrong=%d complete=%d; pairs pure %d wrong %d (%.0f s)" % (
            span, len(trks), nobj, n_multi, len(nights), sc["recall"], sc["pure"], sc["wrong"], sc["complete"], sc["pair_pure"], sc["pair_wrong"], dt), flush=True)
        sc["wrong_groups"] = [dict(truth=g["truth"], chi2=g["chi2_red"]) for g in res["groups"] + res["pairs"] if len(set(g["truth"])) > 1]
        sc["unlinked_objs"] = sorted({t.truth for t in trks if t.id in res["unlinked"]} - {tt for g in res["groups"] for tt in g["truth"]})[:20]
    return out, dict(rc=rc, dc=dc, designations=des)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exp", default="synth")
    ap.add_argument("--out", default="/workspace/work/nl/val")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--nobj", type=int, default=40)
    ap.add_argument("--ndecoy", type=int, default=40)
    ap.add_argument("--mjd0", type=float, default=60950.0)
    ap.add_argument("--floor0", type=float, default=0.3)
    ap.add_argument("--floor-rate", dest="floor_rate", type=float, default=0.05)
    ap.add_argument("--hz-sigmas", type=float, nargs="+", default=[0.1, 0.3, 1.0])
    ap.add_argument("--hz-nights", type=lambda x: [int(v) for v in x.split(",")], nargs="+", default=[[0, 1, 3], [0, 2, 6], [0, 1, 7, 14], [0, 7]])
    ap.add_argument("--nmpc", type=int, default=60)
    ap.add_argument("--spans", type=int, nargs="+", default=[4, 8, 14])
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    res = {}
    for e in a.exp.split(","):
        if e == "synth":
            res["synth"] = synth(a)
        elif e == "horizons":
            res["horizons"], res["horizons_meta"] = horizons_exp(a)
        elif e == "mpc":
            res["mpc"], res["mpc_meta"] = mpc_exp(a)
    json.dump(res, open(os.path.join(a.out, "nightlink_validation_%s.json" % a.exp.replace(",", "_")), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
