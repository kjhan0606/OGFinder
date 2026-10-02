"""Light-curve features, Bazin fit, a pure-numpy random-forest predictor and the SN score for transient light curves.

Light curve = arrays t (MJD or days), flux, err (flux units, any zero point; magnitudes are converted with ZP = 25 by `from_mags`).
`features(lc)` -> dict of FEATURES (+ fit results); `predict(model, feats)` -> class probabilities; `classify(lc, model, host_model, ...)`.
Trained models are JSON (see plugins/lightcurves/train.py); this module needs only numpy/scipy (astropy for the Lomb-Scargle feature).
"""
import json
import math
import numpy as np

ZP = 25.0
FEATURES = ['log_baseline', 'n_obs', 'log_snr_peak', 'log_chi2_const', 'eta', 'log_tau_r', 'log_tau_f', 't0_frac', 'log_chi2_bazin', 'log_dchi2_ratio',
            'half_dur_frac', 'decl_rate', 'rise_rate', 'skew', 'kurt', 'beyond1', 'mid_ratio', 'trend', 'ls_power', 'ls_logP', 'max_slope_frac']
HOST_FEATURES = FEATURES + ['log_host_offset']
SN_CLASSES = ('SNIa', 'SNIbc', 'SNII')


def from_mags(t, mag, mag_err):
    mag = np.asarray(mag, float); me = np.asarray(mag_err, float)
    f = 10.0 ** (-0.4 * (mag - ZP))
    return np.asarray(t, float), f, f * me * math.log(10.0) / 2.5


def to_mag(f):
    with np.errstate(divide='ignore', invalid='ignore'):
        return ZP - 2.5 * np.log10(np.asarray(f, float))


def bazin(t, A, t0, tr, tf, B):
    x = np.clip((t - t0) / tr, -60, 60)
    return A * np.exp(-(t - t0) / tf) / (1.0 + np.exp(-x)) + B


def fit_bazin(t, f, e):
    """Multi-start least-squares Bazin fit.  -> dict(A, t0, tau_r, tau_f, B, chi2, ok)."""
    from scipy.optimize import least_squares
    t = np.asarray(t, float); f = np.asarray(f, float); e = np.maximum(np.asarray(e, float), 1e-12)
    ip = int(np.argmax(f)); fp = f[ip]; base = max(float(np.ptp(t)), 1.0)
    best = None
    lo = [-np.inf, t.min() - base, 0.3, 0.5, -np.inf]; hi = [np.inf, t.max() + base, 60.0, 400.0, np.inf]

    def res(p):
        return (bazin(t, *p) - f) / e
    for t0 in (t[ip],):
        for tr, tf in ((2.0, 15.0), (6.0, 60.0), (0.8, 3.0)):
            A0 = max(fp - np.median(f), e.mean()) * (1.0 + math.exp(-1.0) * 0)  # amplitude near the peak
            p0 = [A0 * 1.5, t0, tr, tf, float(np.min(f))]
            p0[4] = min(max(p0[4], -abs(fp)), abs(fp))
            try:
                r = least_squares(res, p0, bounds=(lo, hi), max_nfev=120)
            except Exception:
                continue
            c = float(np.sum(r.fun ** 2))
            if best is None or c < best[0]:
                best = (c, r.x)
    if best is None:
        return dict(ok=False)
    c, x = best
    return dict(ok=True, A=float(x[0]), t0=float(x[1]), tau_r=float(x[2]), tau_f=float(x[3]), B=float(x[4]), chi2=c)


def _ls(t, f, e):
    try:
        from astropy.timeseries import LombScargle
        base = float(np.ptp(t))
        if len(t) < 6 or base <= 0:
            return np.nan, np.nan
        fmin = 1.0 / base; fmax = 1.0 / 0.15
        freq = np.exp(np.linspace(math.log(fmin), math.log(fmax), 500))
        p = LombScargle(t, f, e).power(freq)
        i = int(np.argmax(p))
        return float(p[i]), float(math.log10(1.0 / freq[i]))
    except Exception:
        return np.nan, np.nan


def features(t, f, e):
    """-> (dict of FEATURES, dict of auxiliary quantities).  NaN where undefined."""
    t = np.asarray(t, float); f = np.asarray(f, float); e = np.asarray(e, float)
    o = np.argsort(t); t, f, e = t[o], f[o], np.maximum(e[o], 1e-12)
    n = len(t)
    base = max(float(np.ptp(t)), 1e-3)
    nan = float('nan')
    F = {k: nan for k in FEATURES}
    aux = {}
    F['log_baseline'] = math.log10(base); F['n_obs'] = float(n)
    ip = int(np.argmax(f)); fmax = float(f[ip]); fmin = float(np.min(f))
    F['log_snr_peak'] = math.log10(max(fmax / e[ip], 1e-3))
    w = 1.0 / e ** 2
    mean = float(np.sum(w * f) / np.sum(w))
    chi2c = float(np.sum(((f - mean) / e) ** 2)) / max(n - 1, 1)
    F['log_chi2_const'] = math.log10(max(chi2c, 1e-3))
    d = np.diff(f)
    var = float(np.var(f, ddof=1)) if n > 2 else nan
    F['eta'] = float(np.sum(d ** 2) / max(n - 1, 1) / var) if var and var > 0 else nan
    amp = max(fmax - fmin, 1e-12)
    # distribution
    sd = float(np.std(f)) or 1e-12
    F['skew'] = float(np.mean(((f - f.mean()) / sd) ** 3)); F['kurt'] = float(np.mean(((f - f.mean()) / sd) ** 4) - 3.0)
    F['beyond1'] = float(np.mean(np.abs(f - mean) > np.median(e)))
    p5, p50, p95 = np.percentile(f, [5, 50, 95])
    F['mid_ratio'] = float((p50 - p5) / max(p95 - p5, 1e-12))
    from scipy.stats import spearmanr
    try:
        F['trend'] = float(spearmanr(t, f)[0])
    except Exception:
        pass
    # half-maximum duration (fraction of the baseline above half of the amplitude)
    half = fmin + 0.5 * amp
    above = f >= half
    if above.sum() >= 1:
        ta = t[above]
        F['half_dur_frac'] = float((ta.max() - ta.min()) / base) if above.sum() > 1 else 0.0
    # slopes in magnitudes: decline after the peak (within 25 d), rise before it (within 15 d), maximum |slope|
    mg = to_mag(f); good = (f > 3 * e)
    tp = t[ip]
    post = good & (t > tp) & (t <= tp + 25.0)
    if post.sum() >= 2 and np.ptp(t[post]) > 0.5:
        F['decl_rate'] = float(np.polyfit(t[post], mg[post], 1)[0])
    pre = good & (t < tp) & (t >= tp - 15.0)
    if pre.sum() >= 2 and np.ptp(t[pre]) > 0.5:
        F['rise_rate'] = float(np.polyfit(t[pre], mg[pre], 1)[0])
    dt = np.diff(t); ok = dt > 1e-3
    if ok.any():
        sl = np.abs(d[ok] / dt[ok]) / max(fmax, 1e-12)
        F['max_slope_frac'] = float(np.max(sl))
    pw, lp = _ls(t, f, e)
    F['ls_power'] = pw; F['ls_logP'] = lp
    bz = fit_bazin(t, f, e) if n >= 5 else dict(ok=False)
    if bz['ok']:
        F['log_tau_r'] = math.log10(bz['tau_r']); F['log_tau_f'] = math.log10(bz['tau_f'])
        F['t0_frac'] = float((bz['t0'] - t[0]) / base)
        c2 = bz['chi2'] / max(n - 5, 1)
        F['log_chi2_bazin'] = math.log10(max(c2, 1e-3))
        F['log_dchi2_ratio'] = math.log10(max(chi2c, 1e-3) / max(c2, 1e-3))
    aux.update(bazin=bz, t_peak=float(tp), f_peak=fmax, m_peak=float(to_mag(fmax)), snr_peak=float(fmax / e[ip]), n=int(n), baseline=base, chi2r_const=chi2c)
    return F, aux


# ------------------------------------------------------------------ random-forest predictor (JSON)
def export_forest(rf, features_, fill, classes, meta=None):
    trees = []
    for est in rf.estimators_:
        t = est.tree_
        v = t.value[:, 0, :]
        v = v / np.maximum(v.sum(1, keepdims=True), 1e-12)
        trees.append(dict(f=[int(a) if l >= 0 else -1 for a, l in zip(t.feature, t.children_left)], thr=[round(float(x), 6) for x in t.threshold],
                          l=[int(x) for x in t.children_left], r=[int(x) for x in t.children_right], v=[[round(float(y), 4) for y in row] for row in v]))
    return dict(schema=1, kind='random_forest', classes=list(classes), features=list(features_), fill=[float(x) for x in fill], trees=trees, meta=meta or {})


def load_model(path):
    m = json.load(open(path))
    m['_trees'] = [dict(f=np.array(t['f'], int), thr=np.array(t['thr'], float), l=np.array(t['l'], int), r=np.array(t['r'], int), v=np.array(t['v'], float)) for t in m['trees']]
    m['_fill'] = np.array(m['fill'], float)
    return m


def predict(model, X):
    """X: (n, n_features) -> probabilities (n, n_classes); NaN replaced by the stored fill values."""
    X = np.array(X, float)
    bad = ~np.isfinite(X)
    if bad.any():
        X[bad] = np.broadcast_to(model['_fill'], X.shape)[bad]
    P = 0.0
    for t in model['_trees']:
        node = np.zeros(len(X), int)
        for _ in range(80):
            fe = t['f'][node]
            leaf = fe < 0
            if leaf.all():
                break
            left = X[np.arange(len(X)), np.where(leaf, 0, fe)] <= t['thr'][node]
            node = np.where(leaf, node, np.where(left, t['l'][node], t['r'][node]))
        P = P + t['v'][node]
    return P / len(model['_trees'])


def vector(F, names, host_offset_re=None):
    v = []
    for k in names:
        if k == 'log_host_offset':
            v.append(math.log10(max(host_offset_re, 0.01)) if host_offset_re is not None and np.isfinite(host_offset_re) else float('nan'))
        else:
            v.append(F[k])
    return v


def classify(t, f, e, model, host_model=None, host_offset_re=None, min_prob=0.5, min_points=5, min_snr=3.0):
    """Classify one light curve.  -> dict(class, p_class, probs, p_sn, features, aux, note).
    class 'unclassified' when fewer than min_points points, peak S/N < min_snr, or the top probability < min_prob."""
    t = np.asarray(t, float); f = np.asarray(f, float); e = np.asarray(e, float)
    ok = np.isfinite(t) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    t, f, e = t[ok], f[ok], e[ok]
    res = dict(cls='unclassified', best_guess='', p_class=float('nan'), probs={}, p_sn=float('nan'), n=int(len(t)), used_host=False, note='')
    if len(t) < min_points:
        res['note'] = 'fewer than %d points' % min_points
        return res
    F, aux = features(t, f, e)
    res['features'] = F; res['aux'] = aux
    if aux['snr_peak'] < min_snr:
        res['note'] = 'peak S/N %.1f below %.1f' % (aux['snr_peak'], min_snr)
        return res
    use_host = host_model is not None and host_offset_re is not None and np.isfinite(host_offset_re)
    mdl = host_model if use_host else model
    P = predict(mdl, [vector(F, mdl['features'], host_offset_re)])[0]
    classes = mdl['classes']
    probs = {c: float(p) for c, p in zip(classes, P)}
    top = max(probs, key=probs.get)
    res.update(probs=probs, p_sn=float(sum(probs.get(c, 0.0) for c in SN_CLASSES)), used_host=bool(use_host), p_class=probs[top])
    res['cls'] = top if probs[top] >= min_prob else 'unclassified'
    if res['cls'] == 'unclassified':
        res['note'] = 'top class %s at p=%.2f below %.2f' % (top, probs[top], min_prob)
    res['best_guess'] = top
    return res
