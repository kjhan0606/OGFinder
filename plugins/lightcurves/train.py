#!/usr/bin/env python3
"""Reproducible training of the light-curve classifiers (random forests on the synthetic library of ogfkit/lcsynth.py).

    python3 plugins/lightcurves/train.py [--per-class 600] [--seed 20260101] [--out plugins/lightcurves] [--jobs 6]

Writes model_lc.json (light-curve features), model_host.json (+ host offset in R_e), training_report.json (held-out accuracy, confusion matrix).
Same seed + same package versions -> same models.  Needs numpy, scipy, astropy, scikit-learn (training only; prediction is pure numpy).
"""
import argparse
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..')))
from ogfkit import lcsynth, lcclass  # noqa: E402


def _feat(lc):
    F, aux = lcclass.features(lc['t'], lc['flux'], lc['err'])
    return [F[k] for k in lcclass.FEATURES], lc['host_offset_re'], aux['snr_peak']


def build(per_class, seed, jobs):
    lcs = lcsynth.make_set(per_class, seed)
    with Pool(jobs) as p:
        res = p.map(_feat, lcs, chunksize=20)
    X = np.array([r[0] for r in res], float)
    ho = np.array([r[1] if r[1] is not None else np.nan for r in res], float)
    y = np.array([lcsynth.CLASSES.index(l['cls']) for l in lcs])
    return X, ho, y, np.array([r[2] for r in res])


def hostcol(ho):
    return np.log10(np.maximum(ho, 0.01))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--per-class', type=int, default=600)
    ap.add_argument('--seed', type=int, default=20260101)
    ap.add_argument('--out', default=HERE)
    ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args(argv)
    from sklearn.ensemble import RandomForestClassifier
    X, ho, y, snr = build(a.per_class, a.seed, a.jobs)
    Xv, hov, yv, snrv = build(max(a.per_class // 4, 50), a.seed + 1, a.jobs)
    report = dict(seed=a.seed, per_class=a.per_class, classes=lcsynth.CLASSES)
    for name, cols, use_host in (('lc', lcclass.FEATURES, False), ('host', lcclass.HOST_FEATURES, True)):
        A = np.c_[X, hostcol(ho)] if use_host else X
        Av = np.c_[Xv, hostcol(hov)] if use_host else Xv
        fill = np.nanmedian(A, axis=0)
        fill = np.where(np.isfinite(fill), fill, 0.0)
        Af = np.where(np.isfinite(A), A, fill); Avf = np.where(np.isfinite(Av), Av, fill)
        rf = RandomForestClassifier(n_estimators=50, min_samples_leaf=6, max_depth=12, random_state=a.seed, n_jobs=a.jobs).fit(Af, y)
        pv = rf.predict_proba(Avf)
        pred = pv.argmax(1)
        conf = np.zeros((len(lcsynth.CLASSES),) * 2, int)
        for t, p in zip(yv, pred):
            conf[t, p] += 1
        report[name] = dict(accuracy=float((pred == yv).mean()), confusion=conf.tolist(), n_val=int(len(yv)))
        m = lcclass.export_forest(rf, cols, fill, lcsynth.CLASSES, meta=dict(seed=a.seed, per_class=a.per_class, kind=name, library='ogfkit/lcsynth.py'))
        with open(os.path.join(a.out, 'model_%s.json' % name), 'w') as f:
            json.dump(m, f, separators=(',', ':'))
        print('%s: held-out accuracy %.3f (%d trees, %d bytes)' % (name, report[name]['accuracy'], len(m['trees']), os.path.getsize(os.path.join(a.out, 'model_%s.json' % name))))
    with open(os.path.join(a.out, 'training_report.json'), 'w') as f:
        json.dump(report, f, indent=1)


if __name__ == '__main__':
    main()
