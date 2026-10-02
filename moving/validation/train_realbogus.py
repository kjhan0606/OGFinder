"""Train the real/bogus classifier (`moving/realbogus.py`) from injection sets on fields OTHER than the evaluation field (BB89).

usage:
  train_realbogus.py --make-sets DIR            # inject + detect on the four other COSMOS ACS visits (3 sets each: 2 plain, 1 CR-heavy/no-DQ)
  train_realbogus.py DIR/tr*.json.pkl [--out moving/data/realbogus_model.json] [--trees 60 --depth 3 --seed 0]     # fit + write the model
  train_realbogus.py --oof DIR/tr*.json.pkl                  # write out-of-fold `rb` into the pickles (input of fit_link_score.py --rb)

Training data = every positive detection with S/N >= 6 of the pickles written by `inject.py` (which labels each detection exactly from the
injected-only image: real = the injected flux makes up a significant part of the detected flux).  Reproducible: fixed seeds (injection seeds
1001..1024, model seed 0), the cached ACS exposures named below, scikit-learn GradientBoostingClassifier.  The evaluation sets inj1..inj12
(BB89) are never read here.  Leave-one-FIELD-out AUC / recall at fixed false-positive rate is printed before the final fit on all fields."""
import sys, os, json, pickle, subprocess, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
import numpy as np
from moving import realbogus as RB

FIELDS = {   # field -> exposures (MAST cache); none overlaps the BB89 visit j8pu38*
    "f1": "j8pu1stiq_flc,j8pu1stmq_flc,j8pu1stpq_flc,j8pu1sttq_flc",
    "f2": "j8pu09rbq_flc,j8pu09rsq_flc,j8pu09rvq_flc,j8pu09s0q_flc",
    "f3": "j8pu33blq_flc,j8pu33btq_flc,j8pu33bwq_flc,j8pu33c0q_flc",
    "f4": "j8pu3oqmq_flc,j8pu3oqqq_flc,j8pu3oqtq_flc,j8pu3oqxq_flc",
    # added in round 2 (item 2): more cached COSMOS ACS visits (different pointings/epochs than BB89; jbhm43 = 4 x 120 s, jboa* = 3 x 340-460 s)
    "f5": "jbhm43uzq_flc,jbhm43viq_flc,jbhm43waq_flc,jbhm43wsq_flc",
    "f6": "jboa40aeq_flc,jboa40aiq_flc,jboa40alq_flc",
    "f7": "jboa41xpq_flc,jboa41xtq_flc,jboa41xwq_flc",
    "f8": "jboa84quq_flc,jboa84qyq_flc,jboa84r4q_flc"}
SETS = [("a", 1000, ""), ("b", 1010, ""), ("c", 1020, "--nodq --cr-extra 3000")]


def make_sets(outdir):
    os.makedirs(outdir, exist_ok=True)
    here = os.path.dirname(os.path.abspath(__file__))
    jobs = []
    for fi, (f, stems) in enumerate(FIELDS.items(), 1):
        for tag, base, extra in SETS:
            out = os.path.join(outdir, "tr%d%s.json" % (fi, tag))
            if os.path.exists(out + ".pkl"):
                continue
            jobs.append(subprocess.Popen([sys.executable, os.path.join(here, "inject.py"), str(base + fi), "30", out, "--field", stems] + extra.split(),
                                         stdout=open(out + ".log", "w"), stderr=subprocess.STDOUT))
            if len(jobs) >= 4:
                [j.wait() for j in jobs]; jobs = []
    [j.wait() for j in jobs]


def load(paths):
    X = []; y = []; g = []; snr = []
    for p in paths:
        dets, truth, offs = pickle.load(open(p, "rb"))
        field = os.path.basename(p)[:3]
        for d in dets:
            if d["sign"] < 0 or d["snr"] < 6.0 or "real" not in d:
                continue
            X.append(RB.feature_vector(d)); y.append(int(d["real"])); g.append(field); snr.append(d["snr"])
    return np.array(X, float), np.array(y), np.array(g), np.array(snr)


def fit(X, y, trees, depth, seed, w=None):
    from sklearn.ensemble import GradientBoostingClassifier
    fill = np.nanmedian(np.where(np.isfinite(X), X, np.nan), axis=0)
    fill = np.where(np.isfinite(fill), fill, 0.0)
    Xf = np.where(np.isfinite(X), X, fill)
    gb = GradientBoostingClassifier(n_estimators=trees, max_depth=depth, learning_rate=0.1, subsample=0.8, min_samples_leaf=20, random_state=seed)
    gb.fit(Xf, y, sample_weight=w)
    return gb, fill


def main():
    a = sys.argv[1:]
    if "--oof" in a:
        a.remove("--oof"); score_out_of_fold([x for x in a if x.endswith(".pkl")]); return
    if "--make-sets" in a:
        make_sets(a[a.index("--make-sets") + 1]); return
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "realbogus_model.json")
    trees, depth, seed = 60, 3, 0
    if "--out" in a:
        i = a.index("--out"); out = a[i + 1]; del a[i:i + 2]
    if "--trees" in a:
        i = a.index("--trees"); trees = int(a[i + 1]); del a[i:i + 2]
    if "--depth" in a:
        i = a.index("--depth"); depth = int(a[i + 1]); del a[i:i + 2]
    if "--seed" in a:
        i = a.index("--seed"); seed = int(a[i + 1]); del a[i:i + 2]
    paths = [x for x in a if x.endswith(".pkl")]
    X, y, g, snr = load(paths)
    print("training detections %d (real %d, bogus %d) from %d sets, fields %s" % (len(y), y.sum(), (1 - y).sum(), len(paths), sorted(set(g))))
    from sklearn.metrics import roc_auc_score
    for hold in sorted(set(g)):
        tr = g != hold
        gb, fill = fit(X[tr], y[tr], trees, depth, seed)
        m = {"init": 0, "lr": 0, "fill": fill}
        Xh = np.where(np.isfinite(X[~tr]), X[~tr], fill)
        p = gb.predict_proba(Xh)[:, 1]
        yh = y[~tr]
        thr = np.sort(p[yh == 0])[::-1][int(0.002 * (yh == 0).sum())]
        print("leave-field-out %s: n=%d real=%d AUC %.3f | at 0.2%% bogus kept: recall %.2f | at 2%% bogus kept: recall %.2f" % (
            hold, len(yh), yh.sum(), roc_auc_score(yh, p), (p[yh == 1] > thr).mean(),
            (p[yh == 1] > np.sort(p[yh == 0])[::-1][int(0.02 * (yh == 0).sum())]).mean()))
    gb, fill = fit(X, y, trees, depth, seed)
    model = RB.export_sklearn(gb, fill, meta=dict(trained_on=[os.path.basename(p) for p in paths], n_real=int(y.sum()), n_bogus=int((1 - y).sum()),
                                                  trees=trees, depth=depth, seed=seed, note="trained on injections in non-BB89 fields; see docs/moving_objects.md"))
    json.dump(model, open(out, "w"), separators=(",", ":"))
    print("wrote", out, os.path.getsize(out), "bytes")
    imp = sorted(zip(gb.feature_importances_, RB.FEATURES), reverse=True)[:10]
    print("top features:", ", ".join("%s %.3f" % (n, v) for v, n in imp))



def score_out_of_fold(paths, trees=60, depth=3, seed=0, write=True):
    """Add an out-of-fold `rb` (model trained on the OTHER fields only) to every detection of the training pickles and rewrite them.  Needed to fit
    the tracklet-level model (`fit_link_score.py --rb`) without feeding it in-sample (over-confident) probabilities."""
    X, y, g, snr = load(paths)
    models = {}
    for hold in sorted(set(g)):
        tr = g != hold
        gb, fill = fit(X[tr], y[tr], trees, depth, seed)
        models[hold] = RB.load_from_dict(RB.export_sklearn(gb, fill))
    for p in paths:
        dets, truth, offs = pickle.load(open(p, "rb"))
        RB.score_dets(dets, model=models[os.path.basename(p)[:3]])
        if write:
            pickle.dump((dets, truth, offs), open(p, "wb"))
        print("scored", p, sum(1 for d in dets if "rb" in d), flush=True)


if __name__ == "__main__":
    main()
