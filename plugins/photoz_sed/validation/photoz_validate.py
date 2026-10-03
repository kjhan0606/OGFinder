#!/usr/bin/env python3
"""Validation of the photo-z distribution tools (writes JSON reports next to this file or into --out).

    photoz_validate.py synth   power / false-positive rate of the PIT tests, accuracy metrics vs known truth (synthetic PDFs of known calibration)
    photoz_validate.py real    SDSS held-out set (photo_z/data/sdss_specphoto.h5) with the repository MDN checkpoint: calibration, percentile accuracy,
                               end-to-end through ds9_photo_z.py --mixture-out and photoz_quality.py; recalibration; representativeness with a
                               bright spec-z selection (truth = the whole held-out sample)
"""
import json
import os
import subprocess
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import photoz_stats as ps  # noqa: E402

H5 = os.path.join(ROOT, 'photo_z', 'data', 'sdss_specphoto.h5')
CKPT = os.path.join(ROOT, 'photo_z', 'data', 'checkpoints', 'photo_z_mdn_best.pt')


def synth(out):
    rep = {}
    n_trials, n = 200, 1000
    cases = {'calibrated': dict(), 'overconfident_x0.8': dict(scale=0.8), 'underconfident_x1.25': dict(scale=1.25), 'bias_+0.01(1+z)': dict(bias=0.01),
             'outliers_3%': dict(frac_out=0.03)}
    tab = {}
    for name, kw in cases.items():
        rej_ks = rej_cvm = rej_chi = 0
        crps, cov68 = [], []
        for t in range(n_trials):
            rng = np.random.default_rng(1000 + t)
            zt = rng.uniform(0.05, 1.0, n)
            sig = 0.03 * (1 + zt) * rng.uniform(0.5, 1.5, n)
            zp = zt + sig * rng.standard_normal(n) + kw.get('bias', 0) * (1 + zt)
            if kw.get('frac_out'):
                o = rng.random(n) < kw['frac_out']
                zp[o] = rng.uniform(0, 1.2, o.sum())
            pred = ps.Predictive('gauss', zphot=zp, sigma=sig * kw.get('scale', 1.0))
            s = ps.pit_report(pred.cdf(zt))
            rej_ks += s['ks_p'] < 0.05
            rej_cvm += s['cvm_p'] < 0.05
            rej_chi += s['chi2_p'] < 0.05
            cov68.append(s['coverage'][1]['observed'])
            crps.append(np.mean(ps.crps(pred, zt)))
        tab[name] = dict(reject_ks=rej_ks / n_trials, reject_cvm=rej_cvm / n_trials, reject_chi2=rej_chi / n_trials, cover68_mean=float(np.mean(cov68)), crps_mean=float(np.mean(crps)))
    rep['pit_tests_n1000_200trials'] = tab
    # truth for point metrics: Gaussian dz with sigma_dz, contaminated
    rng = np.random.default_rng(5)
    N = 200000
    zs = rng.uniform(0.1, 1.2, N)
    dz = 0.025 * rng.standard_normal(N) + 0.004
    cat = rng.random(N) < 0.04
    dz[cat] = rng.uniform(-0.6, 0.6, cat.sum())
    zp = zs + dz * (1 + zs)
    m = ps.point_metrics(zp, zs)
    truth_out = 0.04 * (1 - 0.3 / 0.6) + 0.96 * 2 * stats.norm.sf((0.15 - 0.004) / 0.025) * 0.5 + 0.96 * stats.norm.sf((0.15 + 0.004) / 0.025) * 0.0
    p_out = (0.04 * (1 - (0.15 + 0.15) / 1.2)
             + 0.96 * (stats.norm.sf((0.15 - 0.004) / 0.025) + stats.norm.cdf((-0.15 - 0.004) / 0.025)))
    rep['point_metrics_truth'] = dict(measured=m, true_bias_core=0.004, true_sigma_core=0.025, expected_outlier_frac=float(p_out))
    # optimal width scale
    sc = {}
    for s_true in (0.6, 0.8, 1.0, 1.5):
        rng = np.random.default_rng(11)
        zt = rng.uniform(0.05, 1, 20000)
        sig = 0.03 * (1 + zt)
        zp = zt + sig * rng.standard_normal(len(zt))
        sc[s_true] = ps.optimal_scale(ps.Predictive('gauss', zphot=zp, sigma=sig * s_true), zt)['scale'] * s_true
    rep['optimal_scale_times_applied_scale (expect 1)'] = sc
    json.dump(rep, open(os.path.join(out, 'photoz_synth_report.json'), 'w'), indent=1)
    print(json.dumps(tab, indent=1))


def real(out):
    import h5py
    from photo_z.predict import predict_photoz
    with h5py.File(H5) as f:
        X = f['features'][:]
        z = f['redshift'][:]
        mr = f['extra/mag_r'][:]
    n = len(z)
    idx = np.random.RandomState(42).permutation(n)
    te = idx[n - 1500:]
    held = idx[7000:]
    rep = {}
    r = predict_photoz(X[te], CKPT, return_mixture=True)
    pr = ps.Predictive('mixture', pi=r['mix_pi'], mu=r['mix_mu'], sigma=r['mix_sigma'])
    s, pit, cr, ls, zsc = ps.calibration_summary(pr, z[te])
    rep['test_1500'] = dict(point=ps.point_metrics(r['z_point'], z[te]), pit={k: s['pit'][k] for k in ('ks_p', 'cvm_p', 'mean', 'var', 'shape', 'hist', 'frac_extreme')},
                            coverage=s['pit']['coverage'], crps=s['crps_mean'], logscore=s['logscore_mean'], zscore=s['zscore'])
    rep['test_1500']['gauss_from_zerr'] = ps.calibration_summary(ps.Predictive('gauss', zphot=r['z_point'], sigma=r['z_err']), z[te])[0]['pit']['ks_p']
    sc = ps.optimal_scale(pr, z[te])
    rep['test_1500']['optimal_scale'] = sc
    f = np.arange(len(pit)) % 2
    pc = np.where(f == 0, ps.recalibrate_pit(pit[f == 1], pit), ps.recalibrate_pit(pit[f == 0], pit))
    rep['test_1500']['pit_recal_crossfit_ks_p'] = float(stats.kstest(pc, 'uniform').pvalue)
    # known miscalibration -> detected (MDN widths x0.5 / x2)
    for sfac in (0.5, 2.0):
        q = ps.calibration_summary(pr.scaled(sfac), z[te])[0]
        rep['test_1500']['scaled_x%g' % sfac] = dict(ks_p=q['pit']['ks_p'], shape=q['pit']['shape'], cover68=q['pit']['coverage'][1]['observed'], crps=q['crps_mean'])
    # binned
    rep['binned_zspec'] = ps.binned_metrics(z[te], r['z_point'], z[te], [0, 0.1, 0.2, 0.35, 0.5, 0.7, 1.01], pr)
    rep['binned_mag_r'] = ps.binned_metrics(mr[te], r['z_point'], z[te], [10, 16.5, 17.2, 17.7, 18.5, 25], pr)
    # end-to-end CLI chain: catalog from the raw SDSS photometry -> ds9_photo_z.py --mixture-out -> photoz_quality.py
    tmp = os.path.join(out, 'sdss_e2e')
    os.makedirs(tmp, exist_ok=True)
    cat = os.path.join(tmp, 'cat.tsv')
    with open(cat, 'w') as fh:
        fh.write('NUMBER\tMAG_AUTO\t' + '\t'.join('MAG_' + b for b in 'ugriz') + '\t' + '\t'.join('MAGERR_' + b for b in 'ugriz') + '\tZ_SPEC\n')
        for i in te:
            fh.write('%d\t%.4f\t%s\t%s\t%.5f\n' % (i + 1, mr[i], '\t'.join('%.4f' % v for v in X[i, :5]), '\t'.join('%.4f' % v for v in X[i, 5:10]), z[i]))
    pz = subprocess.run([sys.executable, os.path.join(ROOT, 'ds9', 'library', 'ds9_photo_z.py'), '/dev/null', '--catalog', cat, '--bands', 'u,g,r,i,z', '--checkpoint', CKPT,
                         '--mixture-out', os.path.join(tmp, 'mix.npz')], capture_output=True, text=True)
    cols = pz.stdout.strip().split('\n')
    hdr = [c for c in cols if not c.startswith('#')]
    names = hdr[0].split('\t')
    rows = [dict(zip(names, l.split('\t'))) for l in hdr[1:]]
    zp_pipe = np.array([float(x['PHOTO_Z']) for x in rows])
    rep['e2e'] = dict(mdn_used='Using MDN model' in pz.stderr and 'falling back' not in pz.stderr,
                      max_abs_diff_vs_direct_predict=float(np.max(np.abs(zp_pipe - r['z_point']))), median_abs_diff=float(np.median(np.abs(zp_pipe - r['z_point']))))
    with open(os.path.join(tmp, 'zp.tsv'), 'w') as fh:
        fh.write('NUMBER\tPHOTO_Z\tPHOTO_Z_ERR\n')
        for x in rows:
            fh.write('%s\t%s\t%s\n' % (x['NUMBER'], x['PHOTO_Z'], x['PHOTO_Z_ERR']))
    # merge catalog + photo-z
    with open(cat) as fh:
        L = fh.read().strip().split('\n')
    cc = [l.split('\t') for l in L]
    mm = {x['NUMBER']: x for x in rows}
    with open(os.path.join(tmp, 'full.tsv'), 'w') as fh:
        fh.write('\t'.join(cc[0] + ['PHOTO_Z', 'PHOTO_Z_ERR']) + '\n')
        for l in cc[1:]:
            fh.write('\t'.join(l + [mm[l[0]]['PHOTO_Z'], mm[l[0]]['PHOTO_Z_ERR']]) + '\n')
    q = subprocess.run([sys.executable, os.path.join(HERE, '..', 'photoz_quality.py'), '--catalog', os.path.join(tmp, 'full.tsv'), '--work', os.path.join(tmp, 'q'), '--mode', 'quality',
                        '--mixture-file', os.path.join(tmp, 'mix.npz')], capture_output=True, text=True)
    qr = json.load(open(os.path.join(tmp, 'q', 'pzq_report.json')))
    rep['e2e']['cli_quality'] = dict(n=qr['n_spec'], nmad=qr['point']['sigma_nmad'], outlier=qr['point']['outlier_frac'], pit_ks_p=qr['pit']['ks_p'], cover68=qr['pit']['coverage'][1]['observed'],
                                     crps=qr['crps_mean'])
    # representativeness: bright spec-z selection of the 3000 held-out objects; truth = NMAD / outliers of all 3000
    rh = predict_photoz(X[held], CKPT, return_mixture=False)
    zh = z[held]
    m_r = mr[held]
    F = np.column_stack([m_r, X[held, 0] - X[held, 1], X[held, 1] - X[held, 2], X[held, 2] - X[held, 3], X[held, 3] - X[held, 4]])
    truth = ps.point_metrics(rh['z_point'], zh)
    rep['repr'] = {'truth_all_%d' % len(zh): {k: truth[k] for k in ('sigma_nmad', 'outlier_frac', 'bias_mean', 'absdz_p95')}}
    rngs = np.random.default_rng(21)
    selections = [('mag_r<17.2 (hard cut)', m_r < 17.2),
                  ('soft: p=clip(10^(-0.3 (m-17.2)), 0.05, 1)', rngs.random(len(m_r)) < np.clip(10 ** (-0.3 * (m_r - 17.2)), 0.05, 1)),
                  ('soft: p=clip(10^(-0.3 (m-16.5)), 0.03, 1)', rngs.random(len(m_r)) < np.clip(10 ** (-0.3 * (m_r - 16.5)), 0.03, 1))]
    for label, sel in selections:
        Zt = ps.robust_scale(F)
        Zs = Zt[sel]
        t = ps.knn_classifier_test(Zs, Zt, k=20, nperm=200)
        _, fo, ex = ps.coverage_distance(Zs, Zt)
        d = {}
        for k in (5, 10, 20):
            w = ps.knn_weights(Zs, Zt, k=k)
            wm = ps.point_metrics(rh['z_point'][sel], zh[sel], weights=w)
            d['k%d' % k] = dict(ess=ps.effective_fraction(w), nmad=wm['sigma_nmad'], outlier=wm['outlier_frac'], bias=wm['bias_mean'], p95=wm['absdz_p95'])
        un = ps.point_metrics(rh['z_point'][sel], zh[sel])
        rep['repr'][label] = dict(n_spec=int(sel.sum()), knn_auc=t['auc'], knn_p=t['p'], frac_outside_support=fo, unweighted=dict(nmad=un['sigma_nmad'], outlier=un['outlier_frac'], bias=un['bias_mean'], p95=un['absdz_p95']),
                                  weighted=d, feature_ks=[round(x['ks'], 3) for x in ps.feature_table(F[sel], F, list('abcde'))])
    # no-selection control (random half): must be consistent
    rng = np.random.default_rng(3)
    sel = rng.random(len(zh)) < 0.5
    Zt = ps.robust_scale(F)
    t = ps.knn_classifier_test(Zt[sel], Zt, k=20, nperm=200)
    rep['repr']['random_half_control'] = dict(knn_auc=t['auc'], knn_p=t['p'], frac_outside=ps.coverage_distance(Zt[sel], Zt)[1])
    json.dump(rep, open(os.path.join(out, 'photoz_real_report.json'), 'w'), indent=1, default=float)
    print(json.dumps(rep['test_1500']['pit'], indent=0)[:600])
    print(json.dumps(rep['e2e'], indent=0)[:600])
    print(json.dumps(rep['repr'], indent=0)[:2500])


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'synth'
    out = sys.argv[2] if len(sys.argv) > 2 else HERE
    os.makedirs(out, exist_ok=True)
    {'synth': synth, 'real': real}[mode](out)
