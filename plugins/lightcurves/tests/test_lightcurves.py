"""Light-curve classification tests: held-out synthetic validation, calibration, injection through the real difference-image forced photometry, readers, CLI, reproducible training."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PLUG = os.path.abspath(os.path.join(HERE, '..'))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
sys.path.insert(0, ROOT)
from ogfkit import lcsynth, lcclass  # noqa: E402

M = lcclass.load_model(os.path.join(PLUG, 'model_lc.json'))
H = lcclass.load_model(os.path.join(PLUG, 'model_host.json'))
CLI = os.path.join(PLUG, 'lightcurves.py')
SN = list(lcsynth.SN_CLASSES)


def run(args):
    return subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)


@pytest.fixture(scope='module')
def heldout():
    lcs = lcsynth.make_set(40, 777)            # seed differs from the training seeds 20260101 / 20260102
    out = []
    for lc in lcs:
        r1 = lcclass.classify(lc['t'], lc['flux'], lc['err'], M, None, None, min_prob=0.0)
        r2 = lcclass.classify(lc['t'], lc['flux'], lc['err'], M, H, lc['host_offset_re'], min_prob=0.0)
        out.append((lc['cls'], r1, r2))
    return out


def test_heldout_accuracy_and_sn_selection(heldout):
    cls = np.array([h[0] for h in heldout]); sn = np.isin(cls, SN)
    for k, name in ((1, 'lc'), (2, 'host')):
        g = np.array([h[k]['best_guess'] for h in heldout]); ps = np.array([h[k]['p_sn'] for h in heldout])
        acc = float((g == cls).mean()); sel = ps >= 0.5
        compl = float((sel & sn).sum() / sn.sum()); pur = float((sel & sn).sum() / sel.sum())
        print('%s: 7-class accuracy %.3f; SN selection (P_SN>=0.5) completeness %.3f purity %.3f; per class %s' % (name, acc, compl, pur, {c: round(float((g[cls == c] == c).mean()), 2) for c in lcsynth.CLASSES}))
        assert acc > 0.70
        assert compl > 0.85 and pur > 0.88
    # the three SN types are only partly separable from one light curve
    g = np.array([h[1]['best_guess'] for h in heldout])
    sn_sn = np.mean([g[i] in SN for i in range(len(g)) if sn[i]])
    assert sn_sn > 0.85


def test_calibration_of_sn_score(heldout):
    sn = np.isin([h[0] for h in heldout], SN); ps = np.array([h[1]['p_sn'] for h in heldout])
    ece = 0.0
    for a, b in zip(np.linspace(0, 1, 6)[:-1], np.linspace(0, 1, 6)[1:]):
        m = (ps >= a) & ((ps <= b) if b == 1 else (ps < b))
        if m.any():
            ece += m.sum() / len(ps) * abs(ps[m].mean() - sn[m].mean())
    print('expected calibration error of P_SN: %.3f' % ece)
    assert ece < 0.08


def test_accuracy_grows_with_snr_and_host_helps_agn(heldout):
    cls = np.array([h[0] for h in heldout]); snr = np.array([h[1].get('aux', {}).get('snr_peak', np.nan) for h in heldout])
    g = np.array([h[1]['best_guess'] for h in heldout])
    lo = float((g[snr < 15] == cls[snr < 15]).mean()); hi = float((g[snr >= 30] == cls[snr >= 30]).mean())
    print('accuracy S/N<15: %.2f ; S/N>=30: %.2f' % (lo, hi))
    assert hi > lo + 0.05
    g2 = np.array([h[2]['best_guess'] for h in heldout])
    a1 = float((g[cls == 'AGN'] == 'AGN').mean()); a2 = float((g2[cls == 'AGN'] == 'AGN').mean())
    print('AGN recall: light curve only %.2f, with host offset %.2f' % (a1, a2))
    assert a2 > a1


def test_unclassified_rules():
    t = np.array([0., 1, 2]); f = np.array([1., 2, 1]); e = np.ones(3)
    r = lcclass.classify(t, f, e, M)
    assert r['cls'] == 'unclassified' and 'fewer' in r['note']
    t = np.arange(10.); f = np.random.RandomState(1).randn(10); e = np.ones(10) * 5
    r = lcclass.classify(t, f, e, M)
    assert r['cls'] == 'unclassified' and 'S/N' in r['note']


# ---- injection into difference images, forced photometry with the existing moving/ code
def inject(n_per_class=10, seed=5, n_epoch=16, baseline=90.0):
    from astropy.wcs import WCS
    from moving.transients import lightcurve
    rng = np.random.RandomState(seed)
    classes = list(lcsynth.CLASSES)
    n = n_per_class * len(classes)
    side = 40 * int(np.ceil(np.sqrt(n))) + 40
    w = WCS(naxis=2); w.wcs.ctype = ['RA---TAN', 'DEC--TAN']; w.wcs.crval = [150.0, 2.0]; w.wcs.crpix = [side / 2, side / 2]; w.wcs.cdelt = [-0.2 / 3600, 0.2 / 3600]
    k = int(np.ceil(np.sqrt(n)))
    pos = [(40 + 40 * (i % k), 40 + 40 * (i // k)) for i in range(n)]
    tt = np.sort(rng.rand(n_epoch) * baseline)
    mlim = 25.3                                    # 5 sigma point-source depth in the aperture
    sig_ap = lcsynth.mag2flux(mlim) / 5.0
    r_ap = 3.0; psf_sig = 1.5
    frac = 1 - np.exp(-r_ap ** 2 / (2 * psf_sig ** 2))
    sig_pix = sig_ap / np.sqrt(np.pi * r_ap ** 2)
    truth = []
    fl = np.zeros((n, n_epoch))
    for i in range(n):
        c = classes[i // n_per_class]
        peak = rng.uniform(mlim - 3.8, mlim - 0.8)
        z = rng.uniform(0.03, 0.8) if c in SN else 0.0
        t0 = rng.uniform(-15, baseline * 0.9) if c in SN + ['fast'] else 0.0
        f, p = lcsynth.model_flux(c, tt, rng, peak, z, t0)
        fl[i] = f; truth.append(c)
    yy, xx = np.mgrid[0:side, 0:side]
    diffs = []
    for e in range(n_epoch):
        img = rng.randn(side, side) * sig_pix
        for i, (x, y) in enumerate(pos):
            sl = (slice(y - 12, y + 13), slice(x - 12, x + 13))
            g = np.exp(-((xx[sl] - x) ** 2 + (yy[sl] - y) ** 2) / (2 * psf_sig ** 2))
            img[sl] += fl[i, e] * g / (2 * np.pi * psf_sig ** 2)
        diffs.append((59000.0 + tt[e], w, img, 25.0, 1.0))
    cand = [tuple(float(v) for v in w.all_pix2world([x], [y], 0)) for x, y in pos]
    cand = [(c[0][0] if hasattr(c[0], '__len__') else c[0], c[1][0] if hasattr(c[1], '__len__') else c[1]) for c in cand]
    lcs = lightcurve(cand, diffs, r_pix=r_ap)
    return truth, lcs, fl, frac, cand, tt


def test_injection_through_difference_image_photometry(tmp_path):
    truth, lcs, fl, frac, cand, tt = inject()
    ratio = []
    for l, f in zip(lcs, fl):
        ff = np.array(l['flux']); ok = f > 5 * np.array(l['err'])
        ratio += list(ff[ok] / f[ok])
    print('forced-photometry flux / injected flux: median %.3f (aperture fraction %.3f)' % (np.median(ratio), frac))
    assert abs(np.median(ratio) - frac) < 0.03
    # write the Moving-objects files and classify through the CLI
    wd = tmp_path / 'mov'; wd.mkdir()
    json.dump([dict(id=i, t=l['t'], flux=l['flux'], err=l['err'], mag=l['mag']) for i, l in enumerate(lcs)], open(wd / 'lightcurves.json', 'w'))
    with open(wd / 'transients.tsv', 'w') as f:
        f.write('id\tra\tdec\thost_id\toffset_re\n')
        for i, c in enumerate(cand):
            f.write('%d\t%.7f\t%.7f\t\t\n' % (i, c[0], c[1]))
    cat = tmp_path / 'cat.tsv'; cat.write_text('NUMBER\tALPHA_J2000\tDELTA_J2000\n1\t150.0\t2.0\n')
    r = run(['--catalog', str(cat), '--work', str(tmp_path / 'w'), '--source', 'moving', '--moving-dir', str(wd), '--min-prob', '0.0'])
    assert r.returncode == 0, r.stderr
    res = {}
    for l in open(tmp_path / 'w' / 'lc_results.tsv').read().strip().split('\n')[1:]:
        t = l.split('\t'); res[t[0]] = t
    cols = open(tmp_path / 'w' / 'lc_results.tsv').readline().strip().split('\t')
    gi = cols.index('BEST_GUESS'); pi = cols.index('PSN')
    pred = [res['T%d' % i][gi] for i in range(len(truth))]
    ps = np.array([float(res['T%d' % i][pi]) if res['T%d' % i][pi] else np.nan for i in range(len(truth))])
    tr = np.array(truth); sn = np.isin(tr, SN); ok = np.isfinite(ps)
    acc = float(np.mean([p == t for p, t in zip(pred, truth)]))
    sel = ps >= 0.5
    compl = float((sel & sn).sum() / sn.sum()); pur = float((sel & sn).sum() / max(sel.sum(), 1))
    print('injection (%d light curves): 7-class accuracy %.2f, SN selection completeness %.2f purity %.2f' % (len(truth), acc, compl, pur))
    assert acc > 0.60 and compl > 0.80 and pur > 0.80


def test_readers_flux_mag_and_zp(tmp_path):
    sys.path.insert(0, PLUG)
    import lightcurves as L
    t = np.arange(8) * 3.0 + 59000
    f = np.array([1, 3, 8, 12, 9, 6, 4, 2.5], float) * 10; e = np.ones(8) * 0.5
    p = tmp_path / 'a.csv'
    p.write_text('name,mjd,flux,flux_err\n' + ''.join('x,%g,%g,%g\n' % a for a in zip(t, f, e)))
    r = L.read_long(str(p), 'auto', zp=25.0)
    np.testing.assert_allclose(r[0]['flux'], f)
    r2 = L.read_long(str(p), 'auto', zp=27.5)                       # ZP 27.5 -> ZP 25 scale: x 10^(-1)
    np.testing.assert_allclose(r2[0]['flux'], f * 0.1)
    m = 25 - 2.5 * np.log10(f)
    p2 = tmp_path / 'b.tsv'
    p2.write_text('id\tmjd\tmag\tmag_err\n' + ''.join('y\t%g\t%g\t0.05\n' % a for a in zip(t, m)))
    r3 = L.read_long(str(p2))
    np.testing.assert_allclose(r3[0]['flux'], f, rtol=2e-4)
    np.testing.assert_allclose(r3[0]['err'], f * 0.05 * np.log(10) / 2.5, rtol=2e-4)
    with pytest.raises(SystemExit):
        (tmp_path / 'c.csv').write_text('a,b\n1,2\n'); L.read_long(str(tmp_path / 'c.csv'))


def test_cli_links_and_columns(tmp_path):
    cat = tmp_path / 'cat.tsv'
    cat.write_text('NUMBER\tALPHA_J2000\tDELTA_J2000\n' + ''.join('%d\t%.6f\t%.6f\n' % (i, 150 + i * 0.001, 2.0) for i in range(1, 31)))
    r = subprocess.run([sys.executable, os.path.join(PLUG, 'make_demo.py'), str(tmp_path / 'd'), ','.join(str(i) for i in range(1, 29))], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    r = run(['--catalog', str(cat), '--work', str(tmp_path / 'w'), '--lc-file', str(tmp_path / 'd' / 'lc_demo.tsv'), '--meta-out', str(tmp_path / 'm.json')])
    assert r.returncode == 0, r.stderr
    lines = r.stdout.strip().split('\n')
    cols = lines[0].split('\t')
    assert cols == ['NUMBER', 'LC_ID', 'LC_N', 'LC_CLASS', 'LC_PCLASS', 'LC_PSN', 'LC_PEAK_MJD', 'LC_PEAK_MAG', 'LC_RISE_D', 'LC_DECL', 'LC_CHI2R']
    rows = [dict(zip(cols, l.split('\t') + [''] * 12)) for l in lines[1:]]
    assert len(rows) == 30 and sum(1 for x in rows if x['LC_ID']) == 28 and rows[29]['LC_ID'] == ''
    truth = {l.split('\t')[0]: l.split('\t')[1] for l in open(tmp_path / 'd' / 'truth.tsv').read().strip().split('\n')[1:]}
    ok = sum(1 for x in rows if x['LC_ID'] and x['LC_CLASS'] == truth[x['NUMBER']])
    print('CLI demo: %d/28 classes match truth (unclassified count %d)' % (ok, sum(1 for x in rows if x['LC_CLASS'] == 'unclassified')))
    assert ok >= 17
    assert 'lightcurves' in json.load(open(tmp_path / 'm.json'))
    # position linking: ra/dec columns, no number column
    src = open(tmp_path / 'd' / 'lc_demo.tsv').read().split('\n')
    out = ['id\tmjd\tflux\tflux_err\tra\tdec']
    for l in src[1:]:
        if not l:
            continue
        a = l.split('\t'); n = int(a[1])
        out.append('\t'.join([a[0], a[2], a[3], a[4], '%.6f' % (150 + n * 0.001 + 0.5 / 3600), '2.000100']))
    (tmp_path / 'pos.tsv').write_text('\n'.join(out) + '\n')
    r = run(['--catalog', str(cat), '--work', str(tmp_path / 'w2'), '--lc-file', str(tmp_path / 'pos.tsv'), '--link', 'position', '--link-radius-arcsec', '2'])
    assert r.returncode == 0, r.stderr
    rows2 = [dict(zip(cols, l.split('\t') + [''] * 12)) for l in r.stdout.strip().split('\n')[1:]]
    assert sum(1 for x in rows2 if x['LC_ID']) == 28
    assert all(x['LC_ID'] == 'LC%s' % x['NUMBER'] for x in rows2 if x['LC_ID'])
    r = run(['--catalog', str(cat), '--work', str(tmp_path / 'w3'), '--lc-file', str(tmp_path / 'nope.tsv')])
    assert r.returncode != 0 and 'not found' in r.stderr


def test_plot_png(tmp_path):
    cat = tmp_path / 'cat.tsv'; cat.write_text('NUMBER\n1\n2\n')
    subprocess.run([sys.executable, os.path.join(PLUG, 'make_demo.py'), str(tmp_path / 'd'), '1,2'], check=True)
    out = tmp_path / 'p.png'
    r = run(['--task', 'plot', '--catalog', str(cat), '--work', str(tmp_path / 'w'), '--lc-file', str(tmp_path / 'd' / 'lc_demo.tsv'), '--id', 'LC1', '--out', str(out)])
    assert r.returncode == 0, r.stderr
    from PIL import Image
    assert Image.open(out).size == (640, 340)


def test_training_is_reproducible(tmp_path):
    pytest.importorskip('sklearn')
    outs = []
    for k in range(2):
        d = tmp_path / ('t%d' % k); d.mkdir()
        r = subprocess.run([sys.executable, os.path.join(PLUG, 'train.py'), '--per-class', '24', '--seed', '11', '--out', str(d), '--jobs', '2'], capture_output=True, text=True)
        assert r.returncode == 0, r.stderr
        outs.append((d / 'model_lc.json').read_bytes())
    assert outs[0] == outs[1] and len(outs[0]) > 10000
    rep = json.load(open(PLUG + '/training_report.json'))
    print('shipped models: held-out accuracy lc %.3f host %.3f (seed %d, %d per class)' % (rep['lc']['accuracy'], rep['host']['accuracy'], rep['seed'], rep['per_class']))
    assert rep['lc']['accuracy'] > 0.7
