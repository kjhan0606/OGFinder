"""Photo-z closure validation on the SDSS spec-z sample (photo_z/data/sdss_specphoto.h5, same split as photoz_validate.py real: held = idx[7000:], test = last 1500).

1. PIT recalibration: map learnt on the 1500 calibration objects (held minus test) -> applied to the 1500 test objects (coverage, KS, CRPS, log score before / after);
   5-fold cross-fit on all 3000; recovery of an injected miscalibration (MDN widths x0.5 and x2, cross-fit).
2. Consistency: EAZY (eazy-py 0.8.7, EAZY templates, native engine through plugins/sedcodes, z 0-1.0) vs the repository MDN photo-z on the same SDSS ugriz photometry:
   width c of the normalised difference learnt on the calibration half, flags on the test half, outlier fractions vs spec-z.
Usage: closure_validate.py OUT.json    (needs torch + the MDN checkpoint; EAZY part skipped when $OGF_EAZY_PYTHON / eazy data are missing)
"""
import json, os, subprocess, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('OGF_ROOT', os.path.abspath(os.path.join(HERE, '..', '..', '..')))
sys.path.append(ROOT)
from ogfkit import photoz_stats as ps, pz_closure as pc  # noqa: E402

H5 = os.path.join(ROOT, 'photo_z', 'data', 'sdss_specphoto.h5')
CKPT = os.path.join(ROOT, 'photo_z', 'data', 'checkpoints', 'photo_z_mdn_best.pt')
EZ_PY = os.environ.get('OGF_EAZY_PYTHON', '/workspace/eazy_venv/bin/python')
EZ_DATA = os.environ.get('EAZYCODE', '/workspace/eazy_data/eazy-photoz')


def r3(d):
    return json.loads(json.dumps(d, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else o.tolist() if hasattr(o, 'tolist') else str(o)))


def main(out):
    import h5py
    from photo_z.predict import predict_photoz
    with h5py.File(H5) as f:
        X = f['features'][:]; z = f['redshift'][:]
    n = len(z)
    idx = np.random.RandomState(42).permutation(n)
    held = idx[7000:]; te = idx[n - 1500:]; cal = idx[7000:n - 1500]
    r = predict_photoz(X[held], CKPT, return_mixture=True)
    pred = ps.Predictive('mixture', pi=r['mix_pi'], mu=r['mix_mu'], sigma=r['mix_sigma'])
    zs = z[held]
    is_te = np.isin(held, te)
    rep = dict(n_cal=int((~is_te).sum()), n_test=int(is_te.sum()))
    # 1. cal -> test
    sub = lambda m: pc._subset(pred, np.where(m)[0])
    pc_cal = sub(~is_te).cdf(zs[~is_te])
    rc = pc.PITRecal.fit(pc_cal)
    ptest = sub(is_te)
    before, _ = pc.summary(ptest, zs[is_te])
    after, _ = pc.summary(pc.recalibrated_predictive(ptest, rc), zs[is_te])
    rep['cal_to_test'] = dict(before=before, after=after, map=rc.to_dict())
    # cross-fit on all 3000
    cf = pc.crossfit(pred, zs, k=5)
    rep['crossfit_all'] = dict(before=cf['before'], after_oof=cf['after_oof'])
    for s in (0.5, 2.0):
        cf2 = pc.crossfit(pred.scaled(s), zs, k=5)
        rep['injected_scale_x%g' % s] = dict(before=cf2['before'], after_oof=cf2['after_oof'])
    # 2. EAZY vs MDN consistency
    if os.path.exists(EZ_PY) and os.path.isdir(EZ_DATA):
        names = ['U', 'G', 'R', 'I', 'Z']
        wd = os.path.join(os.path.dirname(os.path.abspath(out)), 'closure_work'); os.makedirs(wd, exist_ok=True)
        cat = os.path.join(wd, 'cat_held.tsv')
        with open(cat, 'w') as fh:
            fh.write('NUMBER\t' + '\t'.join('MAG_SDSS_' + b for b in names) + '\t' + '\t'.join('MAGERR_SDSS_' + b for b in names) + '\tZ_SPEC\n')
            for i in held:
                fh.write('%d\t%s\t%s\t%.5f\n' % (i + 1, '\t'.join('%.4f' % v for v in X[i, :5]), '\t'.join('%.4f' % v for v in X[i, 5:10]), z[i]))
        p = subprocess.run([sys.executable, os.path.join(ROOT, 'plugins', 'sedcodes', 'sedcodes.py'), '--task', 'photoz', '--code', 'eazy', '--engine', 'native', '--python', EZ_PY, '--eazy-data', EZ_DATA,
                            '--catalog', cat, '--work', os.path.join(wd, 'ez'), '--z-max', '1.0', '--mag-columns', ','.join('MAG_SDSS_' + b for b in names)], capture_output=True, text=True)
        lines = [l for l in p.stdout.strip().split('\n') if not l.startswith('eazy photoz')]
        hdr = lines[0].split('\t'); rows = {l.split('\t')[0]: dict(zip(hdr, l.split('\t'))) for l in lines[1:]}
        zE = np.array([float(rows[str(i + 1)]['EZ_Z']) if str(i + 1) in rows else np.nan for i in held])
        sE = np.array([float(rows[str(i + 1)]['EZ_ZERR']) if str(i + 1) in rows else np.nan for i in held])
        zM = r['z_point']; sM = r['z_err']
        ok = np.isfinite(zE) & np.isfinite(sE) & (sE > 0)
        cal_m = ~is_te & ok; te_m = is_te & ok
        # EAZY errors: ZERR is a PDF width that can be huge for multimodal objects; cap the width used for the consistency at 0.3
        sEc = np.clip(sE, 0.005, 0.3)
        cal_c = pc.calibrate(zE[cal_m], sEc[cal_m], zM[cal_m], sM[cal_m])
        cons = pc.consistency(zE[te_m], sEc[te_m], zM[te_m], sM[te_m], c=cal_c['c'], nsig=3.0, zspec=zs[te_m])
        cons0 = pc.consistency(zE[te_m], sEc[te_m], zM[te_m], sM[te_m], c=1.0, nsig=3.0, zspec=zs[te_m])
        rep['eazy_vs_mdn'] = dict(n_cal=int(cal_m.sum()), n_test=int(te_m.sum()), calibration=cal_c,
                                  flagged_calibrated={k: v for k, v in cons.items() if k not in ('flag', 'd')}, flagged_uncalibrated={k: v for k, v in cons0.items() if k not in ('flag', 'd')},
                                  eazy_point=ps.point_metrics(zE[te_m], zs[te_m]), mdn_point=ps.point_metrics(zM[te_m], zs[te_m]))
    json.dump(r3(rep), open(out, 'w'), indent=1)
    return rep


if __name__ == '__main__':
    rep = main(sys.argv[1] if len(sys.argv) > 1 else 'closure_report.json')
    for k in ('cal_to_test', 'crossfit_all', 'injected_scale_x0.5', 'injected_scale_x2'):
        d = rep[k]; a = d.get('before'); b = d.get('after') or d.get('after_oof')
        print(k, 'KS p %.2g -> %.2g' % (a['ks_p'], b['ks_p']), 'cov68 %.3f -> %.3f' % (a['coverage']['0.68'], b['coverage']['0.68']), 'CRPS %.4f -> %.4f' % (a['crps'], b['crps']), 'logscore %.3f -> %.3f' % (a['logscore'], b['logscore']))
    if 'eazy_vs_mdn' in rep:
        e = rep['eazy_vs_mdn']
        print('eazy_vs_mdn', json.dumps(e['calibration']), {k: (round(v, 4) if isinstance(v, float) else v) for k, v in e['flagged_calibrated'].items()})
        print('eazy', {k: round(v, 4) for k, v in e['eazy_point'].items() if isinstance(v, float)})
        print('mdn', {k: round(v, 4) for k, v in e['mdn_point'].items() if isinstance(v, float)})
