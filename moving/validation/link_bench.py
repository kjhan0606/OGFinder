"""Injection benchmark for the tracklet linker.

usage: link_bench.py injN.json.pkl [...] [--legacy] [--no-rescore] [--no-veto] [--no-bound] [--max-tracklets 400]
(--no-rescore = the chi2-like first-stage ranking without the logistic score / bound-orbit cut / clustering)

The .pkl files come from validation/inject.py (HST ACS BB89 field with synthetic movers).  An injected object is "recovered"
when one returned tracklet has >= 3 members within 1 arcsec of the true positions.  `ceiling` is the number of injected
objects that have >= 3 detections within 1" in the linker's input pool at all (detection-side upper bound)."""
import sys, os, time, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from moving import pipeline as P, tracklet as T


def _near(pos, e, ra, dec, tol=1.0):
    return np.hypot((ra - pos[e][2]) * np.cos(np.radians(pos[e][3])), dec - pos[e][3]) * 3600.0 < tol


def score(trs, truth, dets=None, tol=1.0):
    rec = 0; ranks = []
    for tk in truth:
        for r, t in enumerate(trs):
            if sum(_near(tk["pos"], e, a, d, tol) for e, a, d in zip(t["ex"], t["ra"], t["dec"])) >= 3:
                rec += 1; ranks.append(r); break
    return rec, ranks


def _n_true_tracklets(trs, truth, tol=1.0):
    """Returned tracklets that contain >= 3 members within `tol` of ONE injected object (a true object can have several)."""
    n = 0
    for t in trs:
        for tk in truth:
            if sum(_near(tk["pos"], e, a, d) for e, a, d in zip(t["ex"], t["ra"], t["dec"])) >= 3:
                n += 1; break
    return n


def ceiling(dets, truth, offs, tol=1.0):
    """injected objects with >= 3 detections within `tol` arcsec in the *actual linker input pool* (after the S/N, class and
    sharpness selection and the per-exposure cap): no linker can recover more than this."""
    cap = {}
    saved = T.link_exposures
    T.link_exposures = lambda sel, **kw: cap.setdefault("sel", sel) and []
    try:
        P.link_detections(dets, offs, snr_min=8.0, tol_arcsec=1.0, min_exposures=3, max_per_exposure=600)
    finally:
        T.link_exposures = saved
    sel = cap.get("sel", [])
    n = 0
    for tk in truth:
        c = sum(any(d["ex"] == e and _near(tk["pos"], e, d["ra"], d["dec"], tol) for d in sel) for e in tk["pos"])
        n += c >= 3
    return n


def run(pk, legacy=False, max_tracklets=400, **kw):
    dets, truth, offs = pickle.load(open(pk, "rb"))
    if legacy:
        saved = T.link_exposures; T.link_exposures = T.link_exposures_legacy
    t0 = time.time(); vst = {}
    try:
        trs = P.link_detections(dets, offs, veto_stats=vst, snr_min=8.0, tol_arcsec=1.0, min_exposures=3, max_per_exposure=600,
                                max_tracklets=max_tracklets, **kw)
    finally:
        if legacy:
            T.link_exposures = saved
    dt = time.time() - t0
    rec, ranks = score(trs, truth)
    top = {N: sum(1 for r in ranks if r < N) for N in (10, 25, 50, 100)}
    ntrue = _n_true_tracklets(trs, truth)
    hi = [t for t in trs if t.get("prob", 0) >= 0.5]; rec_hi, _ = score(hi, truth)
    return dict(file=os.path.basename(pk), n_inj=len(truth), recovered=rec, ceiling=ceiling(dets, truth, offs), n_tracklets=len(trs),
                false=len(trs) - ntrue, true_tracklets=ntrue, precision=round(ntrue / max(len(trs), 1), 3), in_top10=top[10], in_top25=top[25],
                in_top50=top[50], in_top100=top[100], n_p50=len(hi), rec_p50=rec_hi, prec_p50=round(_n_true_tracklets(hi, truth) / max(len(hi), 1), 2), veto=dict(vst), precision_top50=round(sum(1 for r in ranks if r < 50) / 50.0, 3), seconds=round(dt, 1))


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    legacy = "--legacy" in sys.argv
    kw = dict(rescore=False, bound_orbit=False, cluster_dedupe=False) if "--no-rescore" in sys.argv else {}
    if "--no-bound" in sys.argv:
        kw["bound_orbit"] = False
    if "--no-veto" in sys.argv:
        kw["veto"] = False
    mt = 400
    if "--max-tracklets" in sys.argv:
        mt = int(sys.argv[sys.argv.index("--max-tracklets") + 1]); a = [x for x in a if x != str(mt)]
    tot = 0; totn = 0
    for pk in a:
        r = run(pk, legacy, mt, **kw); tot += r["recovered"]; totn += r["n_inj"]
        print(("legacy " if legacy else "stage1 " if not kw.get("rescore", True) else "new    ") + str(r), flush=True)
    print("TOTAL recovered %d / %d" % (tot, totn))
