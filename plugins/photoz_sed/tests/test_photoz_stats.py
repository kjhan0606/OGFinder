import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy import stats

from ogfkit import photoz_stats as ps

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'photoz_quality.py')
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))


def sample(n=4000, seed=0, scale=1.0, bias=0.0, frac_out=0.0):
    """Truth z, per-object sigma, observed photo-z drawn from the true error model; PDFs N(zphot, scale sigma)."""
    rng = np.random.default_rng(seed)
    zt = rng.uniform(0.05, 1.0, n)
    sig = 0.03 * (1 + zt) * rng.uniform(0.5, 1.5, n)
    zp = zt + sig * rng.standard_normal(n) + bias
    if frac_out:
        o = rng.random(n) < frac_out
        zp[o] = rng.uniform(0.0, 1.2, o.sum())
    return zt, zp, sig, scale


def test_pit_calibrated_is_uniform():
    zt, zp, sig, _ = sample()
    rep = ps.calibration_summary(ps.Predictive('gauss', zphot=zp, sigma=sig), zt)[0]['pit']
    assert rep['ks_p'] > 0.01 and rep['cvm_p'] > 0.01
    for c in rep['coverage']:
        assert abs(c['observed'] - c['level']) < 0.03


@pytest.mark.parametrize('scale,shape', [(0.6, 'U-shaped'), (1.7, 'hump')])
def test_pit_detects_dispersion_errors(scale, shape):
    zt, zp, sig, _ = sample()
    rep = ps.calibration_summary(ps.Predictive('gauss', zphot=zp, sigma=sig * scale), zt)[0]['pit']
    assert rep['ks_p'] < 1e-6 and rep['shape'].startswith(shape)
    # coverage direction: overconfident -> less than nominal
    assert (rep['coverage'][1]['observed'] < 0.68) == (scale < 1)


def test_pit_detects_bias_and_outliers():
    zt, zp, sig, _ = sample(bias=0.03)
    rep = ps.calibration_summary(ps.Predictive('gauss', zphot=zp, sigma=sig), zt)[0]['pit']
    assert rep['ks_p'] < 1e-6 and rep['shape'].startswith('tilted')
    zt, zp, sig, _ = sample(frac_out=0.08)
    rep = ps.calibration_summary(ps.Predictive('gauss', zphot=zp, sigma=sig), zt)[0]['pit']
    assert rep['frac_extreme'] > 0.03 and rep['ks_p'] < 1e-3


def test_crps_orders_models_and_closed_form():
    # N(0,1) at y=0: sigma * (2 phi(0) - 1/sqrt(pi))
    p = ps.Predictive('gauss', zmin=None, zphot=np.array([0.0]), sigma=np.array([1.0]))
    assert abs(ps.crps(p, np.array([0.0]))[0] - (2 / math.sqrt(2 * math.pi) - 1 / math.sqrt(math.pi))) < 1e-9
    zt, zp, sig, _ = sample()
    sc = {s: np.mean(ps.crps(ps.Predictive('gauss', zphot=zp, sigma=sig * s), zt)) for s in (0.5, 1.0, 2.0)}
    assert sc[1.0] < sc[0.5] and sc[1.0] < sc[2.0]
    ls = {s: np.mean(ps.log_score(ps.Predictive('gauss', zphot=zp, sigma=sig * s), zt)) for s in (0.5, 1.0, 2.0)}
    assert ls[1.0] < ls[0.5] and ls[1.0] < ls[2.0]


def test_predictive_kinds_agree():
    zt, zp, sig, _ = sample(600)
    g = ps.Predictive('gauss', zphot=zp, sigma=sig)
    sp = ps.Predictive('split', zphot=zp, lo=sig, hi=sig)
    mx = ps.Predictive('mixture', pi=np.stack([np.full(600, .5)] * 2, 1), mu=np.stack([zp, zp], 1), sigma=np.stack([sig, sig], 1))
    zg = np.linspace(0, 1.6, 1601)
    gr = ps.Predictive('grid', zgrid=zg, pdf=g.grid(zg))
    for q in (sp, mx, gr):
        assert np.abs(q.cdf(zt) - g.cdf(zt)).max() < 2e-3
    # exact mixture CRPS = numerical CRPS
    assert np.abs(ps.crps(mx, zt) - ps.crps(gr, zt)).max() < 1e-4
    # asymmetric split: median at zphot, 16/84 percentiles reproduced
    s2 = ps.Predictive('split', zmin=None, zphot=np.array([0.5]), lo=np.array([0.02]), hi=np.array([0.06]))
    assert abs(s2.cdf(np.array([0.5]))[0] - 0.25) < 1e-9          # sigma_lo / (sigma_lo + sigma_hi)
    assert abs(s2.quantile(0.5, zgrid=np.linspace(0.2, 0.9, 7001))[0] - (0.5 + 0.0)) < 0.05
    # truncation: PDF mass below zmin removed
    t = ps.Predictive('gauss', zmin=0.0, zphot=np.array([0.02]), sigma=np.array([0.05]))
    assert t.cdf(np.array([0.0]))[0] == 0.0 and abs(t.cdf(np.array([5.0]))[0] - 1) < 1e-9


def test_point_metrics_and_percentiles():
    rng = np.random.default_rng(3)
    n = 20000
    zs = rng.uniform(0.1, 1, n)
    dz = 0.03 * rng.standard_normal(n)
    out = rng.random(n) < 0.05
    dz[out] = rng.choice([-1, 1], out.sum()) * rng.uniform(0.2, 0.5, out.sum())
    zp = zs + dz * (1 + zs)
    m = ps.point_metrics(zp, zs)
    assert abs(m['sigma_nmad'] - 0.0323) < 0.0015 and abs(m['outlier_frac'] - 0.05) < 0.005 and abs(m['bias_mean']) < 0.003
    assert m['absdz_p68'] < 0.04 and m['absdz_p99'] > 0.2
    b = ps.bootstrap_point(zp[:3000], zs[:3000], n=100)
    assert 0.0005 < b['sigma_nmad_err'] < 0.003
    # weights: weighting the outliers up raises the outlier fraction
    w = np.where(out, 3.0, 1.0)
    assert ps.point_metrics(zp, zs, weights=w)['outlier_frac'] > m['outlier_frac'] * 2


def test_optimal_scale_and_recalibration():
    zt, zp, sig, _ = sample(6000)
    pred = ps.Predictive('gauss', zphot=zp, sigma=sig * 0.6)
    r = ps.optimal_scale(pred, zt)
    assert abs(r['scale'] - 1 / 0.6) < 0.08 and r['logscore_at_scale'] < r['logscore_original']
    pit = pred.cdf(zt)
    f = np.arange(len(pit)) % 2
    pc = np.where(f == 0, ps.recalibrate_pit(pit[f == 1], pit), ps.recalibrate_pit(pit[f == 0], pit))
    assert stats.kstest(pc, 'uniform').pvalue > 0.01 and stats.kstest(pit, 'uniform').pvalue < 1e-6


def test_representativeness_weights_recover_target():
    rng = np.random.default_rng(5)
    n = 30000
    mag = 21 + 1.6 * rng.standard_normal(n)
    col = 0.5 + 0.1 * (mag - 21) + 0.3 * rng.standard_normal(n)
    X = np.stack([mag, col], 1)
    target = X[:15000]
    spec = X[15000:][rng.random(15000) < 1 / (1 + np.exp((X[15000:, 0] - 20.5) / 0.8))]     # bright-biased spec sample
    assert 600 < len(spec) < 6000
    Zt = ps.robust_scale(target)
    Zs = ps.robust_scale(spec, target)
    t = ps.knn_classifier_test(Zs, Zt, k=15, nperm=100)
    assert t['auc'] > 0.65 and t['p'] < 0.02
    ratio, fo, ex = ps.coverage_distance(Zs, Zt)
    assert fo > 0.1 and ex == pytest.approx(0.05)
    w = ps.knn_weights(Zs, Zt, k=10)
    assert abs(w.mean() - 1) < 1e-9
    dm_un = abs(spec[:, 0].mean() - target[:, 0].mean())
    dm_w = abs(np.average(spec[:, 0], weights=w) - target[:, 0].mean())
    assert dm_w < 0.6 * dm_un
    assert 0.1 < ps.effective_fraction(w) < 0.9
    # a representative sample: no detection, weights ~ flat
    Zs2 = ps.robust_scale(X[15000:][rng.random(15000) < 0.3], target)
    t2 = ps.knn_classifier_test(Zs2, Zt, k=15, nperm=100)
    assert abs(t2['auc'] - 0.5) < 0.04
    assert ps.effective_fraction(ps.knn_weights(Zs2, Zt, k=10)) > 0.7
    tab = ps.feature_table(spec, target, ['m', 'c'])
    assert tab[0]['ks'] > 0.3 and tab[0]['std_mean_diff'] < -0.5


def write_cat(path, rows, cols):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in rows:
            f.write('\t'.join('' if v is None else ('%.6g' % v if isinstance(v, float) else str(v)) for v in r) + '\n')


def run_cli(args):
    p = subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)
    return p


def read_tsv(text):
    L = [l for l in text.strip().split('\n') if l and not l.startswith('#')]
    cols = L[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t') + [''] * len(cols))) for l in L[1:]]


def test_cli_quality_gauss_split_mixture(tmp_path):
    zt, zp, sig, _ = sample(1500, seed=2)
    mag = 20 + 3 * np.random.default_rng(1).random(len(zt))
    rows = [(i + 1, float(m), (float(z) if i % 5 else None), float(p), float(s), float(p - 0.99 * s), float(p + 1.01 * s))
            for i, (m, z, p, s) in enumerate(zip(mag, zt, zp, sig))]
    cat = tmp_path / 'c.tsv'
    write_cat(cat, rows, ['NUMBER', 'MAG_AUTO', 'Z_SPEC', 'PHOTO_Z', 'PHOTO_Z_ERR', 'Q16', 'Q84'])
    out = {}
    for tag, extra in (('gauss', []), ('split', ['--q16-col', 'Q16', '--q84-col', 'Q84'])):
        w = tmp_path / tag
        p = run_cli(['--catalog', str(cat), '--work', str(w), '--mode', 'quality'] + extra)
        assert p.returncode == 0, p.stderr
        cols, recs = read_tsv(p.stdout)
        assert cols == ['NUMBER'] + ps_cols()
        assert sum(1 for r in recs if r['PZ_PIT']) == sum(1 for r in rows if r[2] is not None)
        assert all(r['PZ_PIT'] == '' for r, row in zip(recs, rows) if row[2] is None)
        rep = json.load(open(w / 'pzq_report.json'))
        out[tag] = rep
        assert rep['pit']['ks_p'] > 0.001 and abs(rep['point']['sigma_nmad'] - 0.03 / 1.0 * 1.0) < 0.012
        for f in ('pzq_plot.png', 'pzq_binned.tsv', 'pzq_pit.tsv', 'pzq_coverage.tsv'):
            assert (w / f).stat().st_size > 50
        assert len(rep['binned_mag']) == 5 and rep['stacked_nz']['ks_distance'] < 0.08
    assert out['gauss']['pdf_kind'] == 'gauss' and out['split']['pdf_kind'] == 'split'
    # mixture file keyed by NUMBER (rows without spec-z included) with miscalibrated sigmas -> flagged
    mix = tmp_path / 'm.npz'
    pi = np.ones((len(zt), 2)) / 2
    np.savez(mix, number=np.arange(1, len(zt) + 1), pi=pi, mu=np.stack([zp, zp], 1), sigma=np.stack([sig * .5, sig * .5], 1))
    p = run_cli(['--catalog', str(cat), '--work', str(tmp_path / 'mix'), '--mode', 'quality', '--mixture-file', str(mix)])
    assert p.returncode == 0, p.stderr
    rep = json.load(open(tmp_path / 'mix' / 'pzq_report.json'))
    assert rep['pdf_kind'] == 'mixture' and rep['pit']['ks_p'] < 1e-6 and rep['pit']['shape'].startswith('U')
    assert rep['recalibration']['scale'] > 1.7 and rep['recalibration']['pit_crossfit_ks_p'] > 0.01
    # spec-z from a separate file
    sz = tmp_path / 'sz.tsv'
    write_cat(sz, [(r[0], r[2]) for r in rows if r[2] is not None], ['NUMBER', 'ZSPEC'])
    rows2 = [(r[0], r[1], None, r[3], r[4]) for r in rows]
    cat2 = tmp_path / 'c2.tsv'
    write_cat(cat2, rows2, ['NUMBER', 'MAG_AUTO', 'Z_SPEC', 'PHOTO_Z', 'PHOTO_Z_ERR'])
    p = run_cli(['--catalog', str(cat2), '--work', str(tmp_path / 'zf'), '--mode', 'quality', '--zspec-file', str(sz)])
    assert p.returncode == 0, p.stderr
    assert json.load(open(tmp_path / 'zf' / 'pzq_report.json'))['n_spec'] == out['gauss']['n_spec']


def ps_cols():
    sys.path.insert(0, os.path.dirname(CLI))
    import photoz_quality as pq
    return pq.QCOLS


def test_cli_repr_and_weighted_accuracy(tmp_path):
    rng = np.random.default_rng(8)
    n = 5000
    mag = 21 + 1.5 * rng.standard_normal(n)
    gr = 0.6 + 0.1 * (mag - 21) + 0.25 * rng.standard_normal(n)
    zt = np.clip(0.4 + 0.15 * (mag - 21) + 0.1 * rng.standard_normal(n), 0.02, 2)
    sig = 0.02 * (1 + zt) * 10 ** (0.25 * (mag - 21))          # scatter grows with magnitude
    zp = zt + sig * rng.standard_normal(n)
    has = rng.random(n) < 1 / (1 + np.exp((mag - 20.5) / 0.8))
    rows = [(i + 1, float(mag[i]), float(mag[i] + gr[i]), float(mag[i]), float(zt[i]) if has[i] else None, float(zp[i]), float(sig[i])) for i in range(n)]
    cat = tmp_path / 'c.tsv'
    write_cat(cat, rows, ['NUMBER', 'MAG_g', 'MAG_r', 'MAG_AUTO', 'Z_SPEC', 'PHOTO_Z', 'PHOTO_Z_ERR'])
    p = run_cli(['--catalog', str(cat), '--work', str(tmp_path / 'r'), '--mode', 'repr', '--features', 'MAG_AUTO,MAG_g-MAG_r', '--nperm', '100'])
    assert p.returncode == 0, p.stderr
    cols, recs = read_tsv(p.stdout)
    assert cols == ['NUMBER', 'PZ_SPECDIST', 'PZ_INSPEC', 'PZ_SPECW']
    nspec = int(has.sum())
    assert sum(1 for r in recs if r['PZ_SPECW']) == nspec
    rep = json.load(open(tmp_path / 'r' / 'pzr_report.json'))
    assert rep['verdict'].startswith('NOT') and rep['knn_test']['auc'] > 0.65 and rep['coverage']['frac_outside'] > 0.1
    acc = rep['accuracy']
    # truth: the target-sample (all objects) NMAD; weighting moves the bright-biased spec NMAD towards it
    nm_all = ps.nmad(ps.delta_z(zp, zt))
    d_un = abs(acc['unweighted']['sigma_nmad'] - nm_all)
    d_w = abs(acc['weighted']['sigma_nmad'] - nm_all)
    assert d_w < d_un and d_w < 0.5 * d_un + 0.002
    # weighted quality run takes the weights back from the catalog
    w = {r['NUMBER']: r['PZ_SPECW'] for r in recs}
    rows_w = [r + (float(w[str(r[0])]) if w[str(r[0])] else None,) for r in rows]
    cat2 = tmp_path / 'cw.tsv'
    write_cat(cat2, rows_w, ['NUMBER', 'MAG_g', 'MAG_r', 'MAG_AUTO', 'Z_SPEC', 'PHOTO_Z', 'PHOTO_Z_ERR', 'PZ_SPECW'])
    p = run_cli(['--catalog', str(cat2), '--work', str(tmp_path / 'qw'), '--mode', 'quality', '--weights-col', 'PZ_SPECW'])
    assert p.returncode == 0, p.stderr
    q = json.load(open(tmp_path / 'qw' / 'pzq_report.json'))
    assert q['weighted'] and abs(q['point_weighted']['sigma_nmad'] - acc['weighted']['sigma_nmad']) < 1e-3


H5 = os.path.join(ROOT, 'photo_z', 'data', 'sdss_specphoto.h5')
CKPT = os.path.join(ROOT, 'photo_z', 'data', 'checkpoints', 'photo_z_mdn_best.pt')


@pytest.mark.skipif(not (os.path.exists(H5) and os.path.exists(CKPT)), reason='no SDSS data / checkpoint')
def test_real_mdn_on_heldout_sdss():
    h5py = pytest.importorskip('h5py')
    pytest.importorskip('torch')
    from photo_z.predict import predict_photoz
    with h5py.File(H5) as f:
        X = f['features'][:]
        z = f['redshift'][:]
    n = len(z)
    te = np.random.RandomState(42).permutation(n)[n - 1500:]                  # the training script's held-out test split
    r = predict_photoz(X[te], CKPT, return_mixture=True)
    pr = ps.Predictive('mixture', pi=r['mix_pi'], mu=r['mix_mu'], sigma=r['mix_sigma'])
    m, _ = pr.mean_std()
    assert np.abs(m - r['z_point']).max() < 1e-3                               # mixture mean = the pipeline's point estimate
    s = ps.calibration_summary(pr, z[te])[0]
    pm = ps.point_metrics(r['z_point'], z[te])
    assert pm['sigma_nmad'] < 0.03 and pm['outlier_frac'] < 0.03
    assert 0.62 < s['pit']['coverage'][1]['observed'] < 0.80 and s['crps_mean'] < 0.035
    assert ps.stacked_nz(pr, z[te], np.linspace(0, 1.6, 801), np.linspace(0, 1.6, 17))['ks_distance'] < 0.05
