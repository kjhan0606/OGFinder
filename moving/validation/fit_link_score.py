"""Fit the tracklet-level logistic score used by tracklet.link_exposures (LLR_MODEL, or LLR_MODEL_RB with --rb).

usage: fit_link_score.py inj1.json.pkl ... [--holdout injA.json.pkl ...] [--cands 40000] [--rb]

For each injection set the linker is run with rescore on, dedupe off and a very large `max_tracklets`, so that the feature vector
(`tracklet['features']`, the same code path as in production) is available for the top candidates; candidates are labelled true when
they contain >= 3 detections within 1 arcsec of one injected object.  A logistic regression (L2, C=1) is fitted on all training sets
and printed as a LLR_MODEL literal, with leave-one-set-out recovery (distinct injected objects among the first 400 kept).

--rb  fit the 9-feature model (`tracklet.LLR_FEATURES_RB`: the 7 base features + mean / min real-bogus logit of the members) used when the
      detections carry `rb`.  Use only injection sets of fields other than the evaluation field (`validation/train_realbogus.py --make-sets`)."""
import sys, os, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from moving import pipeline as P, tracklet as T
import link_bench as LB


def candidates(pk, n=40000, feats=None):
    feats = feats or T.LLR_FEATURES
    dets, truth, offs = pickle.load(open(pk, "rb"))
    trs = P.link_detections(dets, offs, snr_min=8.0, tol_arcsec=1.0, min_exposures=3, max_per_exposure=600, max_tracklets=n, dedupe=False,
                            bound_orbit=False)
    lab = []
    for t in trs:
        l = -1
        for k, tk in enumerate(truth):
            if sum(LB._near(tk["pos"], e, a, d) for e, a, d in zip(t["ex"], t["ra"], t["dec"])) >= 3:
                l = k; break
        lab.append(l)
    X = np.array([[t["features"][k] for k in feats] for t in trs])
    return X, np.array(lab)


def recovered_top(X, lab, model, top=400):
    s = T.llr_logit(X, model); o = np.argsort(-s)[:top]
    return len({lab[i] for i in o if lab[i] >= 0})


if __name__ == "__main__":
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    args = sys.argv[1:]; hold = []
    n = 40000
    rb = "--rb" in args
    if rb:
        args.remove("--rb")
    feats = T.LLR_FEATURES_RB if rb else T.LLR_FEATURES
    if "--cands" in args:
        i = args.index("--cands"); n = int(args[i + 1]); del args[i:i + 2]
    if "--holdout" in args:
        i = args.index("--holdout"); hold = args[i + 1:]; args = args[:i]
    data = {pk: candidates(pk, n, feats) for pk in args}
    hdata = {pk: candidates(pk, n, feats) for pk in hold}

    def fit(sets):
        X = np.vstack([data[s][0] for s in sets]); y = np.concatenate([data[s][1] >= 0 for s in sets]).astype(int)
        sc = StandardScaler().fit(X); m = LogisticRegression(C=1.0, max_iter=5000).fit(sc.transform(X), y)
        coef = m.coef_[0] / sc.scale_
        return dict(coef=tuple(round(float(c), 3) for c in coef), intercept=round(float(m.intercept_[0] - (coef * sc.mean_).sum()), 2))
    for te in args:
        mod = fit([s for s in args if s != te])
        print("leave-out %-16s recovered in top400 (candidates, no dedupe): %d / %d" % (os.path.basename(te), recovered_top(*data[te], mod), len(set(data[te][1]) - {-1})))
    mod = fit(args)
    print("LLR_FEATURES%s =" % ("_RB" if rb else ""), feats)
    print("LLR_MODEL%s =" % ("_RB" if rb else ""), mod)
    shipped = T.LLR_MODEL_RB if rb and T.LLR_MODEL_RB else T.LLR_MODEL
    for pk in hold:
        print("holdout %-16s  fitted model: %d, shipped model: %d, chi2 score only: %d (of %d reachable by any candidate)" % (
            os.path.basename(pk), recovered_top(*hdata[pk], mod), recovered_top(*hdata[pk], shipped),
            len({hdata[pk][1][i] for i in np.argsort(hdata[pk][0][:, 6])[:400] if hdata[pk][1][i] >= 0}), len(set(hdata[pk][1]) - {-1})))
