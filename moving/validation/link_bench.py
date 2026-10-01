"""Injection benchmark for the tracklet linker.

usage: link_bench.py injN.json.pkl [...] [--legacy] [--max-tracklets 400]

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
    t0 = time.time()
    try:
        trs = P.link_detections(dets, offs, snr_min=8.0, tol_arcsec=1.0, min_exposures=3, max_per_exposure=600,
                                max_tracklets=max_tracklets, **kw)
    finally:
        if legacy:
            T.link_exposures = saved
    dt = time.time() - t0
    rec, ranks = score(trs, truth)
    top50 = sum(1 for r in ranks if r < 50)
    return dict(file=os.path.basename(pk), n_inj=len(truth), recovered=rec, ceiling=ceiling(dets, truth, offs), n_tracklets=len(trs),
                false=len(trs) - rec, in_top50=top50, seconds=round(dt, 1))


if __name__ == "__main__":
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    legacy = "--legacy" in sys.argv
    mt = 400
    if "--max-tracklets" in sys.argv:
        mt = int(sys.argv[sys.argv.index("--max-tracklets") + 1]); a = [x for x in a if x != str(mt)]
    tot = 0; totn = 0
    for pk in a:
        r = run(pk, legacy, mt); tot += r["recovered"]; totn += r["n_inj"]
        print(("legacy " if legacy else "new    ") + str(r), flush=True)
    print("TOTAL recovered %d / %d" % (tot, totn))
