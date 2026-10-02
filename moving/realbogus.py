"""Fast real/bogus classifier for difference-image detections.

Features are cheap shape / cosmic-ray / context quantities computed in `detect._shape_features` (flux profile in rings, PSF-fit chi2, fine
structure and Laplacian at the peak, L.A.Cosmic and archive-flag coverage, ...) plus a few existing ones (S/N, elongation, template S/N, negative
fraction).  The model is a small gradient-boosted tree ensemble (scikit-learn `GradientBoostingClassifier`, 60 trees of depth 3) trained on
*injected sources vs. false candidates* of the injection benchmark on fields OTHER than the evaluation field (`validation/train_realbogus.py`
makes the training sets from the other COSMOS ACS visits in the cache and writes `moving/data/realbogus_model.json`).  Inference needs numpy only:
the trees are stored as arrays in the JSON file (about 60 kB) and evaluated by a vectorised walk, ~1 us per detection.

`score_dets(dets)` adds `rb` (probability) to every detection of S/N >= 6 and `rb_raw` (logit); it returns a small summary.  Without a model file
`rb` is not set (`available()` is False) and the linker falls back to its previous S/N ordering."""
import os, json
import numpy as np

FEATURES = ("snr", "log_snr", "trail", "elong", "a_pix", "b_pix", "npix", "tpl_snr", "neg_frac", "pos_frac", "sharp_rel", "on_cr",
            "f_r1", "f_r2", "f_r3", "pk_nb", "psf_chi2", "psf_amp", "fine", "lap", "n_hi", "asym", "lac3", "lac_n7", "arch3", "n_near",
            "lac_frac", "arch_frac", "fill", "trail_peak_frac", "near_bad", "log_trail_len")
MODEL_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "realbogus_model.json")
_MODEL = None


def feature_vector(d):
    """One detection (dict) -> feature row (NaN where not defined; the model handles NaN through the stored fill value)."""
    g = lambda k, dflt=np.nan: (lambda v: float(v) if v is not None and np.isfinite(float(v)) else dflt)(d.get(k, dflt))
    snr = max(g("snr", 0.0), 1e-3)
    trail = 1.0 if d.get("channel") == "trail" else 0.0
    sh, shp = g("sharp"), g("sharp_psf")
    L = g("trail_len_pix", 0.0)
    return [snr, np.log10(snr), trail, g("elong"), g("a_pix"), g("b_pix"), g("npix"), g("tpl_snr"), g("neg_frac"), g("pos_frac"),
            sh / shp if np.isfinite(sh) and np.isfinite(shp) and shp > 0 else np.nan, 1.0 if d.get("on_cr") in (True, 1, "True", "1") else 0.0,
            g("f_r1"), g("f_r2"), g("f_r3"), g("pk_nb"), g("psf_chi2"), g("psf_amp"), g("fine"), g("lap"), g("n_hi"), g("asym"),
            g("lac3", 0.0), g("lac_n7", 0.0), g("arch3", 0.0), g("n_near", 0.0), g("lac_frac"), g("arch_frac"), g("fill"), g("trail_peak_frac"),
            1.0 if d.get("near_bad") in (True, 1, "True", "1") else 0.0, np.log10(1.0 + L)]


def load_from_dict(m):
    m = dict(m)
    m["trees"] = [dict(f=np.array(t["f"], int), thr=np.array(t["thr"], float), l=np.array(t["l"], int), r=np.array(t["r"], int), v=np.array(t["v"], float))
                  for t in m["trees"]]
    m["fill"] = np.array(m["fill"], float)
    return m


def load(path=None, force=False):
    global _MODEL
    if _MODEL is not None and not force and path is None:
        return _MODEL
    p = path or MODEL_PATH
    if not os.path.isfile(p):
        _MODEL = None if path is None else _MODEL
        return None
    m = load_from_dict(json.load(open(p)))
    if path is None:
        _MODEL = m
    return m


def available():
    return load() is not None


def _tree_predict(t, X):
    node = np.zeros(len(X), int)
    for _ in range(64):
        f = t["f"][node]
        leaf = f < 0
        if leaf.all():
            break
        go_left = X[np.arange(len(X)), np.where(leaf, 0, f)] <= t["thr"][node]
        nxt = np.where(go_left, t["l"][node], t["r"][node])
        node = np.where(leaf, node, nxt)
    return t["v"][node]


def predict(X, model=None):
    """Logit and probability for a feature matrix (n x len(FEATURES))."""
    m = model or load()
    X = np.array(X, float)
    bad = ~np.isfinite(X)
    if bad.any():
        X[bad] = np.broadcast_to(m["fill"], X.shape)[bad]
    z = np.full(len(X), m["init"], float)
    for t in m["trees"]:
        z += m["lr"] * _tree_predict(t, X)
    return z, 1.0 / (1.0 + np.exp(-z))


def score_dets(dets, min_snr=6.0, model=None):
    m = model or load()
    idx = [i for i, d in enumerate(dets) if d["sign"] > 0 and d["snr"] >= min_snr]
    if m is None or not idx:
        return dict(available=m is not None, scored=0)
    X = np.array([feature_vector(dets[i]) for i in idx])
    z, p = predict(X, m)
    for i, a, b in zip(idx, z, p):
        dets[i]["rb_raw"] = float(a); dets[i]["rb"] = float(b)
    return dict(available=True, scored=len(idx), n_ge_0p5=int((p >= 0.5).sum()), n_ge_0p1=int((p >= 0.1).sum()))


def export_sklearn(gb, fill, feature_names=FEATURES, meta=None):
    """sklearn GradientBoostingClassifier (binary, init = prior) -> JSON-able dict used by `load`."""
    trees = []
    for est in gb.estimators_[:, 0]:
        t = est.tree_
        trees.append(dict(f=[int(f) if l >= 0 else -1 for f, l in zip(t.feature, t.children_left)],
                          thr=[round(float(x), 6) for x in t.threshold], l=[int(x) for x in t.children_left], r=[int(x) for x in t.children_right],
                          v=[round(float(x[0][0]), 6) for x in t.value]))
    prior = float(gb.init_.class_prior_[1]) if hasattr(gb.init_, "class_prior_") else 0.5
    init = float(np.log(prior / (1 - prior)))
    return dict(schema=1, features=list(feature_names), init=init, lr=float(gb.learning_rate), fill=[float(x) for x in fill], trees=trees, meta=meta or {})
