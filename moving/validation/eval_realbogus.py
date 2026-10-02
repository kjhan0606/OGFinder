"""Held-out evaluation of real/bogus models on the BB89 injection sets (never used for training).

usage: eval_realbogus.py MODEL.json [MODEL2.json ...] --sets inj1.json.pkl inj2.json.pkl ...  [--snr 6]
For each model and each group (plain inj1-8, CR-heavy inj9-12, all) prints AUC and the recall of real detections at 0.2 % / 2 % bogus kept
(threshold taken on the same held-out detections: operating-point comparison, not a calibrated threshold)."""
import sys, os, pickle
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from moving import realbogus as RB


def load_sets(paths, smin):
    D = {}
    for p in paths:
        dets, truth, offs = pickle.load(open(p, "rb"))
        D[os.path.basename(p).split(".")[0]] = [d for d in dets if d["sign"] > 0 and d["snr"] >= smin and "real" in d]
    return D


def score(model_path, dets):
    import json
    m = RB.load_from_dict(json.load(open(model_path)))
    X = np.array([RB.feature_vector(d) for d in dets], float)
    return RB.predict(X, m)[1]


def auc(y, p):
    from sklearn.metrics import roc_auc_score
    return roc_auc_score(y, p)


def recall_at(y, p, frac):
    b = np.sort(p[y == 0])[::-1]
    if len(b) == 0: return float("nan")
    thr = b[min(len(b) - 1, int(frac * len(b)))]
    return float((p[y == 1] > thr).mean())


def main():
    a = sys.argv[1:]
    smin = 6.0
    if "--snr" in a:
        i = a.index("--snr"); smin = float(a[i + 1]); del a[i:i + 2]
    i = a.index("--sets"); models = a[:i]; sets = a[i + 1:]
    D = load_sets(sets, smin)
    groups = {"inj1-8 plain": [k for k in D if k[3:].isdigit() and int(k[3:]) <= 8], "inj9-12 CR-heavy": [k for k in D if k[3:].isdigit() and int(k[3:]) >= 9]}
    groups["all"] = list(D)
    for mp in models:
        print("model", mp)
        for gname, keys in groups.items():
            dets = [d for k in keys for d in D[k]]
            if not dets: continue
            y = np.array([int(d["real"]) for d in dets])
            p = score(mp, dets)
            print("  %-18s n=%6d real=%4d  AUC %.3f  recall@0.2%%bogus %.2f  recall@2%%bogus %.2f" % (gname, len(y), y.sum(), auc(y, p), recall_at(y, p, 0.002), recall_at(y, p, 0.02)))


if __name__ == "__main__":
    main()
