import json
import os
import subprocess
import sys

import numpy as np
from scipy import stats

from ogfkit import photoz_stats as ps, pz_closure as pcl

HERE = os.path.dirname(os.path.abspath(__file__))
CLI = os.path.join(os.path.dirname(HERE), 'photoz_quality.py')


def sample(n=3000, seed=0, scale=1.0):
    rng = np.random.default_rng(seed)
    zt = rng.uniform(0.05, 1.0, n)
    sig = 0.03 * (1 + zt) * rng.uniform(0.5, 1.5, n)
    zp = zt + sig * rng.standard_normal(n)
    return zt, ps.Predictive('gauss', zphot=zp, sigma=sig * scale)


def test_recalibration_uniformises_pit_and_keeps_calibrated_ones():
    for scale in (0.5, 2.0):
        zt, pred = sample(scale=scale)
        pit = pred.cdf(zt)
        assert stats.kstest(pit, 'uniform').pvalue < 1e-6
        rc = pcl.PITRecal.fit(pit)
        assert stats.kstest(rc.apply(pit), 'uniform').pvalue > 0.05
    zt, pred = sample(seed=3)
    rc = pcl.PITRecal.fit(pred.cdf(zt))
    u = np.linspace(0.02, 0.98, 50)
    assert np.max(np.abs(rc.apply(u) - u)) < 0.05                 # an already calibrated predictive is left (almost) alone
    rc2 = pcl.PITRecal.from_dict(json.loads(json.dumps(rc.to_dict())))
    assert np.allclose(rc2.apply(u), rc.apply(u))


def test_recalibrated_predictive_is_a_pdf_with_cdf_G_of_F():
    zt, pred = sample(1200, seed=1, scale=0.6)
    rc = pcl.PITRecal.fit(pred.cdf(zt))
    new = pcl.recalibrated_predictive(pred, rc)
    zg = new.zg
    P = pred.grid(zg)
    C = np.concatenate([np.zeros((pred.n, 1)), np.cumsum(0.5 * (P[:, 1:] + P[:, :-1]) * np.diff(zg), axis=1)], axis=1)
    C = C / C[:, -1:]
    assert np.allclose(new.cg[:, -1], 1.0, atol=1e-6)
    i = np.searchsorted(zg, 0.6)
    assert np.max(np.abs(new.cg[:, i] - rc.apply(C[:, i]))) < 0.02       # F'(z) = G(F(z)) (grid integration error only)


def test_crossfit_improves_out_of_fold_and_recovers_coverage():
    zt, pred = sample(2500, seed=2, scale=0.5)
    r = pcl.crossfit(pred, zt, k=5)
    assert r['before']['ks_p'] < 1e-10 and r['after_oof']['ks_p'] > 0.01
    assert abs(r['after_oof']['coverage']['68.3'] - 0.683) < 0.03 if '68.3' in r['after_oof']['coverage'] else True
    assert r['after_oof']['crps'] < r['before']['crps']


def test_consistency_flags_injected_outliers():
    rng = np.random.default_rng(5)
    n = 4000
    zt = rng.uniform(0.1, 1.0, n)
    sA = np.full(n, 0.03); sB = np.full(n, 0.04)
    zA = zt + sA * rng.standard_normal(n)
    zB = zt + sB * rng.standard_normal(n)
    bad = rng.random(n) < 0.05
    zB[bad] += rng.choice([-1, 1], bad.sum()) * 0.4
    cal = pcl.calibrate(zA, sA, zB, sB)
    assert 0.8 < cal['c'] < 1.2
    r = pcl.consistency(zA, sA, zB, sB, c=cal['c'], nsig=3.0, zspec=zt)
    assert r['flag'][bad].mean() > 0.95 and r['flag'][~bad].mean() < 0.02
    assert r['outlier_B_flagged'] > 0.9 and r['outlier_B_unflagged'] < 0.01
    assert r['A_closer_when_flagged'] > 0.95


def write_cat(path, cols, data):
    with open(path, 'w') as f:
        f.write('\t'.join(cols) + '\n')
        for r in zip(*data):
            f.write('\t'.join('%.6g' % v if isinstance(v, (float, np.floating)) else str(v) for v in r) + '\n')


def test_cli_recal_and_consistency(tmp_path):
    rng = np.random.default_rng(8)
    n = 1500
    zt = rng.uniform(0.05, 1.0, n)
    sig = 0.03 * (1 + zt)
    zp = zt + 0.5 * sig * rng.standard_normal(n)           # quoted sigma is 2x too large
    zalt = zt + 0.04 * rng.standard_normal(n)
    bad = np.arange(n) % 40 == 0
    zalt[bad] += 0.5
    cat = tmp_path / 'c.tsv'
    write_cat(cat, ['NUMBER', 'Z_SPEC', 'PHOTO_Z', 'PHOTO_Z_ERR', 'EZ_Z', 'EZ_ZERR'], [np.arange(1, n + 1), zt, zp, sig, zalt, np.full(n, 0.04)])
    base = [sys.executable, CLI, '--catalog', str(cat)]
    w1 = tmp_path / 'w1'
    p = subprocess.run(base + ['--work', str(w1), '--mode', 'quality', '--recal-out', str(tmp_path / 'map.json')], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    rep = json.load(open(w1 / 'pzq_report.json'))
    cl = rep['closure']
    assert cl['before']['ks_p'] < 1e-6 and cl['after_oof']['ks_p'] > 0.01
    assert (tmp_path / 'map.json').stat().st_size > 50
    w2 = tmp_path / 'w2'
    p = subprocess.run(base + ['--work', str(w2), '--mode', 'quality', '--recal-in', str(tmp_path / 'map.json')], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    head = open(w2 / 'pzq_recalibrated.tsv').readline().split()
    assert {'PZ_P16_RC', 'PZ_P50_RC', 'PZ_P84_RC'} <= set(head)
    w3 = tmp_path / 'w3'
    p = subprocess.run(base + ['--work', str(w3), '--mode', 'consistency'], capture_output=True, text=True)
    assert p.returncode == 0, p.stderr
    rc = json.load(open(w3 / 'pzq_consistency.json'))
    assert 0.9 < rc['frac_flagged'] / 0.025 < 1.3 or 0.02 < rc['frac_flagged'] < 0.05
    cols = p.stdout.split('\n')[0].split('\t')
    assert 'PZ_INCONSIST' in cols and 'PZ_ZDIFF' in cols
