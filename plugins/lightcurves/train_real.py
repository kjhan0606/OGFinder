#!/usr/bin/env python3
"""Validate / retrain the light-curve classifier on REAL ZTF light curves (made by fetch_ztf.py).

    python3 plugins/lightcurves/train_real.py ZTF_LCS.json [--folds 5] [--synthetic 400] [--out plugins/lightcurves/real_validation_report.json]
                                              [--write-model plugins/lightcurves/model_lc_ztf.json] [--jobs 6]

Reports, on real objects only (stratified k-fold, the same objects never in train and test):
  A  shipped synthetic-trained model (model_lc.json), no retraining                     -> accuracy, balanced accuracy, SN-vs-rest AUC
  B  random forest trained on the real training folds only
  C  synthetic library (--synthetic per class) + real training folds (real weighted by --real-weight)
Classes with no real examples ('fast', 'static') are never a true class here; predictions into them count as errors in A.
Also accepts a generic JSON list of dicts with keys label, t, mag, err[, sign] (or t, flux, err) so other surveys (ALeRCE exports, PLAsTiCC-style
tables converted to this layout) can be plugged in."""
import argparse, json, os, sys
from multiprocessing import Pool
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from ogfkit import lcsynth, lcclass  # noqa: E402


def to_flux(o):
    if 'flux' in o:
        return np.array(o['t'], float), np.array(o['flux'], float), np.array(o['err'], float)
    t, f, e = lcclass.from_mags(o['t'], o['mag'], o['err'])
    return t, f * np.array(o.get('sign', [1] * len(f)), float), e


def _feat(o):
    t, f, e = to_flux(o)
    F, aux = lcclass.features(t, f, e)
    return [F[k] for k in lcclass.FEATURES], aux['snr_peak']


def _feat_syn(lc):
    F, aux = lcclass.features(lc['t'], lc['flux'], lc['err'])
    return [F[k] for k in lcclass.FEATURES]


def fill_of(X):
    f = np.nanmedian(np.where(np.isfinite(X), X, np.nan), axis=0)
    return np.where(np.isfinite(f), f, 0.0)


def metrics(y, P, classes, sn_idx):
    pred = P.argmax(1)
    present = sorted(set(y))
    acc = float((pred == y).mean())
    bal = float(np.mean([(pred[y == c] == c).mean() for c in present]))
    from sklearn.metrics import roc_auc_score
    snp = P[:, sn_idx].sum(1); sny = np.isin(y, sn_idx).astype(int)
    auc = float(roc_auc_score(sny, snp)) if 0 < sny.sum() < len(sny) else float('nan')
    sn_pred = snp >= 0.5
    comp = float(sn_pred[sny == 1].mean()); pur = float(sny[sn_pred].mean()) if sn_pred.any() else float('nan')
    conf = np.zeros((len(classes),) * 2, int)
    for t, p in zip(y, pred): conf[t, p] += 1
    return dict(accuracy=acc, balanced_accuracy=bal, sn_auc=auc, sn_completeness=comp, sn_purity=pur, confusion=conf.tolist())


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('lcs'); ap.add_argument('--folds', type=int, default=5); ap.add_argument('--synthetic', type=int, default=400)
    ap.add_argument('--real-weight', type=float, default=5.0)
    ap.add_argument('--out', default=os.path.join(HERE, 'real_validation_report.json')); ap.add_argument('--write-model', default='')
    ap.add_argument('--drop', default='', help='comma list of features removed from B and C (confound ablation, e.g. n_obs,log_baseline)')
    ap.add_argument('--min-snr', type=float, default=3.0); ap.add_argument('--jobs', type=int, default=6); ap.add_argument('--seed', type=int, default=7)
    a = ap.parse_args(argv)
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.model_selection import StratifiedKFold
    objs = json.load(open(a.lcs))
    classes = lcsynth.CLASSES
    with Pool(a.jobs) as p:
        R = p.map(_feat, objs, chunksize=8)
    X = np.array([r[0] for r in R], float); snr = np.array([r[1] for r in R])
    y = np.array([classes.index(o['label']) for o in objs])
    keep = snr >= a.min_snr
    X, y, snr = X[keep], y[keep], snr[keep]
    print('real objects %d (peak S/N >= %.1f), per class %s' % (len(y), a.min_snr, {classes[c]: int((y == c).sum()) for c in sorted(set(y))}))
    sn_idx = [classes.index(c) for c in lcsynth.SN_CLASSES]
    # synthetic library for A (shipped model is trained on it) and C
    syn = lcsynth.make_set(a.synthetic, 20260401)
    with Pool(a.jobs) as p:
        Xs = np.array(p.map(_feat_syn, syn, chunksize=20), float)
    ys = np.array([classes.index(l['cls']) for l in syn])
    shipped = lcclass.load_model(os.path.join(HERE, 'model_lc.json'))
    PA = lcclass.predict(shipped, X)
    report = dict(n_real=int(len(y)), per_class={classes[c]: int((y == c).sum()) for c in sorted(set(y))}, classes=classes, min_snr=a.min_snr,
                  A_shipped_synthetic_model=metrics(y, PA, classes, sn_idx))
    # snr-binned accuracy of A
    pa = PA.argmax(1)
    report['A_by_peak_snr'] = {('%g-%g' % (lo, hi)): dict(n=int(((snr >= lo) & (snr < hi)).sum()), accuracy=float((pa == y)[(snr >= lo) & (snr < hi)].mean()))
                               for lo, hi in ((3, 10), (10, 30), (30, 1e9)) if ((snr >= lo) & (snr < hi)).sum() > 5}
    cols = [i for i, k in enumerate(lcclass.FEATURES) if k not in a.drop.split(',')]
    fnames = [lcclass.FEATURES[i] for i in cols]
    XA, XsA = X, Xs
    X, Xs = X[:, cols], Xs[:, cols]
    skf = StratifiedKFold(a.folds, shuffle=True, random_state=a.seed)
    PB = np.zeros((len(y), len(classes))); PC = np.zeros_like(PB)
    for tr, te in skf.split(X, y):
        f = fill_of(X[tr])
        Xt = np.where(np.isfinite(X[tr]), X[tr], f)
        cw = {c: 1.0 for c in set(y)}
        rf = RandomForestClassifier(n_estimators=200, min_samples_leaf=3, max_depth=12, class_weight='balanced', random_state=a.seed, n_jobs=a.jobs).fit(Xt, y[tr])
        PB[np.ix_(te, rf.classes_)] = rf.predict_proba(np.where(np.isfinite(X[te]), X[te], f))
        XC = np.r_[Xs, X[tr]]; yC = np.r_[ys, y[tr]]; w = np.r_[np.ones(len(ys)), np.full(len(tr), a.real_weight)]
        fC = fill_of(XC)
        rfc = RandomForestClassifier(n_estimators=200, min_samples_leaf=3, max_depth=12, random_state=a.seed, n_jobs=a.jobs).fit(np.where(np.isfinite(XC), XC, fC), yC, sample_weight=w)
        PC[np.ix_(te, rfc.classes_)] = rfc.predict_proba(np.where(np.isfinite(X[te]), X[te], fC))
    report['B_real_only_cv'] = metrics(y, PB, classes, sn_idx)
    report['C_synthetic_plus_real_cv'] = metrics(y, PC, classes, sn_idx)
    report['settings'] = dict(drop=a.drop, folds=a.folds, synthetic_per_class=a.synthetic, real_weight=a.real_weight, seed=a.seed)
    json.dump(report, open(a.out, 'w'), indent=1)
    for k in ('A_shipped_synthetic_model', 'B_real_only_cv', 'C_synthetic_plus_real_cv'):
        m = report[k]
        print('%-28s acc %.3f  balanced %.3f  SN-vs-rest AUC %.3f  SN completeness %.3f purity %.3f' % (k, m['accuracy'], m['balanced_accuracy'], m['sn_auc'], m['sn_completeness'], m['sn_purity']))
    if a.write_model:
        XC = np.r_[Xs, X]; yC = np.r_[ys, y]; w = np.r_[np.ones(len(ys)), np.full(len(y), a.real_weight)]
        fC = fill_of(XC)
        rfc = RandomForestClassifier(n_estimators=50, min_samples_leaf=6, max_depth=12, random_state=a.seed, n_jobs=a.jobs).fit(np.where(np.isfinite(XC), XC, fC), yC, sample_weight=w)
        m = lcclass.export_forest(rfc, fnames, fC, classes, meta=dict(kind='lc', library='ogfkit/lcsynth.py + %d real ZTF objects (%s)' % (len(y), os.path.basename(a.lcs)),
                                                                             real_weight=a.real_weight))
        json.dump(m, open(a.write_model, 'w'), separators=(',', ':'))
        print('wrote', a.write_model, os.path.getsize(a.write_model), 'bytes')


if __name__ == '__main__':
    main()
