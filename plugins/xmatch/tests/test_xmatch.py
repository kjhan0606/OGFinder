"""Cross-match tests: synthetic truth (completeness / purity / shift / chance estimate), readers, CLI (file, pixel, fake TAP server), live Gaia (skipped offline)."""
import csv
import json
import os
import subprocess
import sys
import threading
import http.server
import urllib.parse

import numpy as np
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
CLI = os.path.join(HERE, '..', 'xmatch.py')
sys.path.insert(0, ROOT)
from ogfkit import xmatch as xm  # noqa: E402


def field(n=400, n_ref_extra=300, sigma=0.15, shift=(0.0, 0.0), seed=1, size_arcmin=6.0):
    """Catalog 1 (n objects, truth positions + sigma noise) and catalog 2 (same sources with independent noise + shift, plus unrelated extras)."""
    rng = np.random.RandomState(seed)
    ra0, dec0 = 150.0, 2.0
    x = (rng.rand(n) - .5) * size_arcmin / 60.0; y = (rng.rand(n) - .5) * size_arcmin / 60.0
    ra_t = ra0 + x / np.cos(np.radians(dec0 + y)); dec_t = dec0 + y
    ra1 = ra_t + rng.randn(n) * sigma / 3600 / np.cos(np.radians(dec_t)); dec1 = dec_t + rng.randn(n) * sigma / 3600
    ra2 = ra_t + (rng.randn(n) * sigma + shift[0]) / 3600 / np.cos(np.radians(dec_t)); dec2 = dec_t + (rng.randn(n) * sigma + shift[1]) / 3600
    xe = (rng.rand(n_ref_extra) - .5) * size_arcmin / 60.0; ye = (rng.rand(n_ref_extra) - .5) * size_arcmin / 60.0
    ra2 = np.r_[ra2, ra0 + xe / np.cos(np.radians(dec0 + ye))]; dec2 = np.r_[dec2, dec0 + ye]
    return ra1, dec1, ra2, dec2, n


def test_completeness_purity_vs_radius():
    ra1, dec1, ra2, dec2, n = field()
    out = []
    for rad in (0.3, 0.6, 1.0, 2.0, 4.0):
        r = xm.crossmatch(ra1, dec1, ra2, dec2, rad)['res']
        b = r['best']
        good = (b == np.arange(n))                       # true counterpart is index i in catalog 2
        comp = float(good.mean()); purity = float(good.sum() / max(np.sum(b >= 0), 1))
        out.append((rad, comp, purity))
    # sigma 0.15" per axis on each catalog -> pair separation Rayleigh with 0.21": 0.6" is >99 % complete
    d = {o[0]: o for o in out}
    assert abs(d[0.6][1] - (1 - np.exp(-0.6 ** 2 / (2 * 0.15 ** 2 * 2)))) < 0.02 and d[0.6][2] > 0.97     # Rayleigh expectation 0.983
    assert d[0.3][1] < d[0.6][1]                          # tighter radius loses real pairs
    assert d[4.0][2] >= 0.9                               # bigger radius: impurity grows only mildly at this density
    print('completeness/purity', out)


def test_chance_estimate_matches_truth():
    # catalog 2 has NO true counterparts: every match is by chance
    rng = np.random.RandomState(5)
    n = 500
    ra1 = 150 + (rng.rand(n) - .5) * .1; dec1 = 2 + (rng.rand(n) - .5) * .1
    est, true = [], []
    for s in range(12):
        r2 = np.random.RandomState(100 + s)
        ra2 = 150 + (r2.rand(800) - .5) * .1; dec2 = 2 + (r2.rand(800) - .5) * .1
        r = xm.crossmatch(ra1, dec1, ra2, dec2, 3.0)
        true.append(r['n_matched']); est.append(r['chance']['mean'])
    t, e = float(np.mean(true)), float(np.mean(est))
    print('chance matches: true %.1f estimated %.1f' % (t, e))
    assert abs(e - t) / t < 0.1


def test_systematic_shift_recovered_and_applied():
    ra1, dec1, ra2, dec2, n = field(shift=(0.45, -0.30), sigma=0.12, seed=3)
    r = xm.crossmatch(ra1, dec1, ra2, dec2, 2.0)
    s = r['shift']
    print('shift', s)
    assert abs(s['dra'] - 0.45) < 0.03 and abs(s['ddec'] + 0.30) < 0.03
    r2 = xm.crossmatch(ra1, dec1, ra2, dec2, 0.6, apply_shift=True)
    r3 = xm.crossmatch(ra1, dec1, ra2, dec2, 0.6)
    ok2 = float((r2['res']['best'] == np.arange(n)).mean()); ok3 = float((r3['res']['best'] == np.arange(n)).mean())
    print('completeness at 0.6": raw %.3f shifted %.3f' % (ok3, ok2))
    assert ok2 > 0.97 > ok3 + 0.2


def test_ambiguity_and_flags():
    ra1 = np.array([10.0, 10.0 + 20 / 3600, 11.0]); dec1 = np.array([0.0, 0.0, 0.0])
    ra2 = np.array([10.0 + 0.4 / 3600, 10.0 - 0.8 / 3600, 10.0 + 20.3 / 3600]); dec2 = np.zeros(3)
    r = xm.crossmatch(ra1, dec1, ra2, dec2, 1.5)['res']
    assert list(r['best']) == [0, 2, -1]
    assert r['n'][0] == 2 and r['flag'][0] & 2 and r['flag'][0] & 4          # ambiguous but mutual nearest
    assert r['flag'][1] == 1 | 4 and r['flag'][2] == 0
    assert abs(r['second_sep'][0] - 0.8) < 0.01


def test_pixel_match():
    rng = np.random.RandomState(2)
    x = rng.rand(100) * 500; y = rng.rand(100) * 500
    x2 = x + rng.randn(100) * 0.3 + 1.0; y2 = y + rng.randn(100) * 0.3 - 0.5
    r = xm.crossmatch(None, None, None, None, 3.0, pixel=True, x1=x, y1=y, x2=x2, y2=y2)
    assert r['n_matched'] == 100
    assert abs(r['shift']['dra'] - 1.0) < 0.15 and abs(r['shift']['ddec'] + 0.5) < 0.15


def test_readers_and_sexagesimal(tmp_path):
    from astropy.table import Table
    ra = np.array([150.1, 150.2, 150.3]); dec = np.array([2.1, 2.2, 2.3]); mag = np.array([18.0, 19.5, 20.1])
    t = Table({'RAJ2000': ra, 'DEJ2000': dec, 'gmag': mag, 'name': ['a', 'b', 'c']})
    p = {}
    for ext in ('fits', 'csv', 'ecsv', 'vot'):
        p[ext] = str(tmp_path / ('t.' + ext))
        t.write(p[ext], format={'fits': 'fits', 'csv': 'ascii.csv', 'ecsv': 'ascii.ecsv', 'vot': 'votable'}[ext])
    with open(tmp_path / 't.tsv', 'w') as f:
        f.write('# comment\nRA\tDEC\tgmag\n150.1\t2.1\t18.0\n150.2\t2.2\t19.5\n150.3\t2.3\t20.1\n')
    p['tsv'] = str(tmp_path / 't.tsv')
    for ext, path in p.items():
        cols, d = xm.read_table(path)
        rc, dc = xm.find_columns(cols)
        assert rc and dc, (ext, cols)
        a, b = xm.parse_sky(d[rc], d[dc])
        np.testing.assert_allclose(a, ra, atol=1e-6)
        np.testing.assert_allclose(b, dec, atol=1e-6)
    a, b = xm.parse_sky(np.array(['10:00:24.0', '10:00:48.0'], dtype=object), np.array(['+02:06:00', '-02:12:00'], dtype=object))
    np.testing.assert_allclose(a, [150.1, 150.2], atol=1e-6)         # sexagesimal RA is hours
    np.testing.assert_allclose(b, [2.1, -2.2], atol=1e-6)


def write_cat(path, ra, dec, x=None, y=None):
    with open(path, 'w') as f:
        f.write('NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\n')
        for i in range(len(ra)):
            f.write('%d\t%g\t%g\t%.8f\t%.8f\n' % (i + 1, (x[i] if x is not None else 0), (y[i] if y is not None else 0), ra[i], dec[i]))


def run(args, **kw):
    return subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True, **kw)


def read_out(txt):
    lines = txt.strip().split('\n')
    cols = lines[0].split('\t')
    return cols, [dict(zip(cols, l.split('\t') + [''] * 20)) for l in lines[1:]]


def test_cli_file_mode(tmp_path):
    ra1, dec1, ra2, dec2, n = field(n=200, n_ref_extra=100, shift=(0.3, 0.2), seed=7)
    cat = str(tmp_path / 'cat.tsv'); write_cat(cat, ra1, dec1)
    ref = str(tmp_path / 'ref.csv')
    with open(ref, 'w', newline='') as f:
        w = csv.writer(f); w.writerow(['objid', 'ra', 'dec', 'gmag', 'plx'])
        for i in range(len(ra2)):
            w.writerow(['S%d' % i, '%.8f' % ra2[i], '%.8f' % dec2[i], 15 + i * 0.01, 0.5])
    r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--source', 'file', '--ref-file', ref, '--radius-arcsec', '1.5', '--ref-id-col', 'objid',
             '--copy-columns', 'gmag,plx', '--apply-shift', '--meta-out', str(tmp_path / 'm.json')])
    assert r.returncode == 0, r.stderr
    cols, rows = read_out(r.stdout)
    assert cols == ['NUMBER', 'XM_SEP', 'XM_N', 'XM_FLAG', 'XM_DRA', 'XM_DDEC', 'XM_ID', 'XM_V1', 'XM_V2', 'XM_V3', 'XM_V4']
    assert len(rows) == 200
    good = sum(1 for i, rw in enumerate(rows) if rw['XM_ID'] == 'S%d' % i)
    print('CLI file mode: %d/200 correct ids' % good)
    assert good >= 195
    assert abs(float(rows[0]['XM_V1']) - 15.0) < 1e-6 and float(rows[0]['XM_V2']) == 0.5
    s = json.load(open(tmp_path / 'w' / 'xmatch_summary.json'))
    assert abs(s['shift']['dra'] - 0.3) < 0.05 and s['n_matched'] >= 195
    assert os.path.getsize(tmp_path / 'w' / 'xmatch_pairs.tsv') > 100
    assert 'xmatch' in json.load(open(tmp_path / 'm.json'))


def test_cli_pixel_and_errors(tmp_path):
    rng = np.random.RandomState(3)
    x = rng.rand(50) * 400; y = rng.rand(50) * 400
    cat = str(tmp_path / 'cat.tsv'); write_cat(cat, np.zeros(50), np.zeros(50), x, y)
    ref = str(tmp_path / 'ref.csv')
    with open(ref, 'w') as f:
        f.write('px,py,m\n' + ''.join('%g,%g,%g\n' % (x[i] + .2, y[i] - .1, i) for i in range(50)))
    r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--mode', 'pixel', '--ref-file', ref, '--ref-x-col', 'px', '--ref-y-col', 'py', '--radius-arcsec', '1.0', '--copy-columns', 'm'])
    assert r.returncode == 0, r.stderr
    _, rows = read_out(r.stdout)
    assert all(float(rw['XM_V1']) == i for i, rw in enumerate(rows))
    assert 0.2 < float(rows[0]['XM_SEP']) < 0.3
    r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--source', 'file', '--ref-file', str(tmp_path / 'nope.csv')])
    assert r.returncode != 0 and 'not found' in r.stderr
    r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--source', 'tap', '--tap-url', 'http://127.0.0.1:9/tap'])
    assert r.returncode != 0 and 'allow' in r.stderr.lower()


class TapHandler(http.server.BaseHTTPRequestHandler):
    table = ''
    seen = []

    def do_POST(self):
        n = int(self.headers.get('Content-Length', 0))
        form = urllib.parse.parse_qs(self.rfile.read(n).decode())
        TapHandler.seen.append((self.path, form))
        body = TapHandler.table.encode()
        self.send_response(200); self.send_header('Content-Type', 'text/csv'); self.send_header('Content-Length', str(len(body))); self.end_headers(); self.wfile.write(body)

    def log_message(self, *a):
        pass


def test_cli_tap_fake_server(tmp_path):
    ra1, dec1, ra2, dec2, n = field(n=100, n_ref_extra=40, seed=11)
    cat = str(tmp_path / 'cat.tsv'); write_cat(cat, ra1, dec1)
    TapHandler.table = 'source_id,ra,dec,phot_g_mean_mag\n' + ''.join('%d,%.8f,%.8f,%.2f\n' % (1000 + i, ra2[i], dec2[i], 16 + i * .01) for i in range(len(ra2)))
    TapHandler.seen = []
    srv = http.server.HTTPServer(('127.0.0.1', 0), TapHandler)
    th = threading.Thread(target=srv.serve_forever, daemon=True); th.start()
    try:
        r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--source', 'tap', '--allow-network', '--tap-url', 'http://127.0.0.1:%d/tap' % srv.server_address[1], '--radius-arcsec', '1.0',
                 '--ref-id-col', 'source_id', '--copy-columns', 'phot_g_mean_mag'], timeout=60)
    finally:
        srv.shutdown()
    assert r.returncode == 0, r.stderr
    path, form = TapHandler.seen[0]
    assert path == '/tap/sync' and form['LANG'] == ['ADQL'] and 'CIRCLE(' in form['QUERY'][0] and 'gaiadr3' in form['QUERY'][0]
    _, rows = read_out(r.stdout)
    good = sum(1 for i, rw in enumerate(rows) if rw['XM_ID'] == str(1000 + i))
    print('fake TAP: %d/100 correct' % good)
    assert good >= 97
    assert os.path.exists(tmp_path / 'w' / 'xmatch_tap.csv')


def _live_ok():
    try:
        import urllib.request
        urllib.request.urlopen('https://gea.esac.esa.int/tap-server/tap/availability', timeout=10)
        return True
    except Exception:
        return False


@pytest.mark.skipif(not os.environ.get('OGF_LIVE_TESTS') and not _live_ok(), reason='no network')
def test_live_gaia_dr3_recovers_injected_offset(tmp_path):
    """Real Gaia DR3 sources of a field (M51) -> perturbed copy (0.1" noise, +0.35"/-0.25" offset) as the 'catalog' -> match through the TAP path."""
    import urllib.request, urllib.parse, io
    q = ("SELECT TOP 600 source_id, ra, dec, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', 202.4696, 47.1952, 0.15)) AND phot_g_mean_mag < 20.5 ORDER BY phot_g_mean_mag")
    data = urllib.request.urlopen('https://gea.esac.esa.int/tap-server/tap/sync', urllib.parse.urlencode({'REQUEST': 'doQuery', 'LANG': 'ADQL', 'FORMAT': 'csv', 'QUERY': q}).encode(), timeout=60).read().decode()
    rows = list(csv.DictReader(io.StringIO(data)))
    assert len(rows) > 100, len(rows)
    rng = np.random.RandomState(4)
    ra = np.array([float(r['ra']) for r in rows]); dec = np.array([float(r['dec']) for r in rows])
    ra_c = ra - (0.35 + rng.randn(len(ra)) * 0.1) / 3600 / np.cos(np.radians(dec)); dec_c = dec - (-0.25 + rng.randn(len(ra)) * 0.1) / 3600
    cat = str(tmp_path / 'cat.tsv'); write_cat(cat, ra_c, dec_c)
    r = run(['--catalog', cat, '--work', str(tmp_path / 'w'), '--source', 'tap', '--allow-network', '--radius-arcsec', '2.0', '--ref-id-col', 'source_id', '--copy-columns', 'phot_g_mean_mag'], timeout=180)
    assert r.returncode == 0, r.stderr
    _, out = read_out(r.stdout)
    good = sum(1 for rw, rr in zip(out, rows) if rw['XM_ID'] == rr['source_id'])
    s = json.load(open(tmp_path / 'w' / 'xmatch_summary.json'))
    print('live Gaia: %d/%d correct, shift %.3f %.3f, chance %.1f, reference rows %d' % (good, len(rows), s['shift']['dra'], s['shift']['ddec'], s['chance']['mean'], s['n_reference']))
    assert good / len(rows) > 0.95
    assert abs(s['shift']['dra'] - 0.35) < 0.05 and abs(s['shift']['ddec'] + 0.25) < 0.05
