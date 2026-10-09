"""plugins/report/report.py on real data (m51, HUDF F160W crop): contents of the generated HTML are checked from the Python side.

The catalogs come from ogfmeas/sextract.py (skipped when the script or the test images are absent); the review
columns are written the way the GUI writes them (REVIEW / REVIEW_NOTE / REVIEW_TIME).
"""
import base64
import html as htmlmod
import io
import json
import os
import re
import subprocess
import sys

import numpy as np
import pytest

import report as R

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')
SEX = os.path.join(ROOT, 'ogfmeas', 'sextract.py')
M51 = os.path.join(FITS, 'm51.fits')
HUDF = os.path.join(FITS, 'hudf_f160w.fits')
SCRIPT = os.path.join(ROOT, 'plugins', 'report', 'report.py')

needs_m51 = pytest.mark.skipif(not (os.path.exists(SEX) and os.path.exists(M51)), reason='ogfmeas/sextract.py or m51.fits missing')
needs_hudf = pytest.mark.skipif(not (os.path.exists(SEX) and os.path.exists(HUDF)), reason='ogfmeas/sextract.py or hudf_f160w.fits missing')


def extract(image, tmp, *args):
    out = subprocess.run([sys.executable, SEX, image, *args], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, timeout=300).stdout.decode()
    p = os.path.join(str(tmp), os.path.basename(image) + '.tsv')
    with open(p, 'w') as f:
        f.write(out)
    return p


def review_catalog(src, dst, decisions):
    """decisions: {NUMBER: (status, note)}; adds the three REVIEW columns like the GUI does"""
    cols, rows = R.read_tsv(src)
    ni = cols.index('NUMBER')
    lines = ['\t'.join(cols + ['REVIEW', 'REVIEW_NOTE', 'REVIEW_TIME'])]
    for r in rows:
        st, note = decisions.get(r[ni], ('', ''))
        lines.append('\t'.join(r + [st, note, '2026-10-01 12:00:00' if st else '']))
    with open(dst, 'w') as f:
        f.write('\n'.join(lines))
    return dst


def thumbs_of(text):
    return [base64.b64decode(m) for m in re.findall(r'<img class="thumb"[^>]*src="data:image/png;base64,([^"]+)"', text)]


def png_array(b):
    from PIL import Image
    return np.asarray(Image.open(io.BytesIO(b)))


@pytest.fixture(scope='module')
def m51_cat(tmp_path_factory):
    tmp = tmp_path_factory.mktemp('m51')
    return extract(M51, tmp, '--detect-thresh', '5'), tmp


@needs_m51
def test_m51_report_contents(m51_cat, tmp_path):
    cat, _ = m51_cat
    cols, rows = R.read_tsv(cat)
    nums = [r[cols.index('NUMBER')] for r in rows]
    dec = {nums[0]: ('accept', 'bright core & "ring" <b>x</b>'), nums[1]: ('reject', ''), nums[2]: ('uncertain', 'check PSF'),
           nums[3]: ('accept', '')}
    rc = review_catalog(cat, str(tmp_path / 'rev.tsv'), dec)
    prov = {'git_head': 'abc123', 'git_dirty_sources': 0, 'ds9_version': '8.7', 'python': 'python3',
            'steps': [{'seq': 1, 'step': 'extract', 'title': 'Extract sources', 'class': 'auto', 'ms': '240', 'failed': '0',
                       'argv': ['ds9_sextract', M51, '--detect-thresh', '5'], 'payload': {}, 'wall': '2026-10-01T12:00:00+0900'},
                      {'seq': 2, 'step': 'review.set', 'title': 'Review: accept 2', 'class': 'note', 'ms': '', 'failed': '0',
                       'argv': [], 'payload': {'numbers': [nums[0], nums[3]], 'status': 'accept', 'note': ''}, 'wall': '2026-10-01T12:01:00+0900'}],
            'plugin_params': {'extract': {'detect-thresh': '5'}}, 'extraction_params': {'param,detect-thresh': '5'}}
    pj = tmp_path / 'prov.json'
    pj.write_text(json.dumps(prov))
    out = str(tmp_path / 'r.html')
    res = R.generate(rc, M51, out, provenance=str(pj), include='all', title='M51 <test>', author='J. Kim', max_rows=50)
    t = open(out, encoding='utf-8').read()
    assert t.startswith('<!DOCTYPE html>') and '</html>' in t
    assert htmlmod.escape('M51 <test>') in t and 'M51 <test>' not in t
    assert res['counts'] == {'accept': 2, 'reject': 1, 'uncertain': 1, 'none': len(rows) - 4}
    assert res['listed'] == 50 and res['selected'] == len(rows)
    assert 'accepted: <b>2</b>' in t and 'rejected: <b>1</b>' in t and 'uncertain: <b>1</b>' in t
    assert t.count('class="cand st-accept" data-review="accept"') >= 2
    assert 'bright core &amp; &quot;ring&quot; &lt;b&gt;x&lt;/b&gt;' in t and '<b>x</b>' not in t
    th = thumbs_of(t)
    assert len(th) == 50 == res['thumbnails']
    for b in th[:5]:
        assert b[:8] == b'\x89PNG\r\n\x1a\n'
        assert png_array(b).shape[:2] == (128, 128)
    assert 'abc123' in t and 'review.set' in t and '--detect-thresh 5' in t
    assert len(re.findall(r'<code>[0-9a-f]{64}</code>', t)) == 2          # sha256 of the image and of the catalog
    assert 'Review history' in t and '2026-10-01T12:01:00+0900' in t
    assert not re.findall(r'(?:src|href)="(?!data:|#)[^"]*"', t)           # self-contained
    assert 'http://' not in t and 'https://' not in t


@needs_m51
def test_cutout_is_centred_on_the_candidate(m51_cat):
    """real image + real catalog: the 5x5 centre of the cut-out is brighter than the cut-out as a whole for compact sources, and the
    synthetic case pins the 1-based pixel convention exactly"""
    cat, _ = m51_cat
    data, _h = R.load_image(M51)
    cols, rows = R.read_tsv(cat)
    ci = {c: i for i, c in enumerate(cols)}
    compact = [r for r in rows if float(r[ci['FWHM_IMAGE']]) < 6 and float(r[ci['FLAGS']]) == 0]
    compact = sorted(compact, key=lambda r: float(r[ci['MAG_AUTO']]))[:8]
    assert len(compact) >= 3
    good = 0
    for r in compact:
        c = R.cutout(data, float(r[ci['X_IMAGE']]), float(r[ci['Y_IMAGE']]), 33)
        if np.nanmean(c[14:19, 14:19]) > np.nanmean(c) + 0.5 * np.nanstd(c):
            good += 1
    assert good >= len(compact) - 1, (good, len(compact))
    # synthetic: a spike at 1-based (x=40, y=25) is at the centre pixel (16,16) of a 33 px cut-out
    d = np.zeros((60, 80), dtype='float32')
    d[25 - 1, 40 - 1] = 9.0
    c = R.cutout(d, 40, 25, 33)
    assert np.unravel_index(np.nanargmax(c), c.shape) == (16, 16)


def test_cutout_edges_are_nan_padded():
    d = np.arange(100, dtype='float32').reshape(10, 10)
    c = R.cutout(d, 1, 1, 5)                  # corner: 1-based (1,1) is d[0,0]
    assert c.shape == (5, 5) and np.isnan(c[0, 0]) and c[2, 2] == d[0, 0] and c[4, 4] == d[2, 2]
    assert np.isnan(R.cutout(d, 500, 500, 5)).all()
    b = R.render_png(R.to_uint8(R.cutout(d, 500, 500, 5), 0, 1, 'zscale'), 2, True)
    assert b[:4] == b'\x89PNG'


@needs_m51
def test_include_modes(m51_cat, tmp_path):
    cat, _ = m51_cat
    cols, rows = R.read_tsv(cat)
    nums = [r[cols.index('NUMBER')] for r in rows]
    dec = {nums[0]: ('accept', ''), nums[1]: ('reject', ''), nums[2]: ('uncertain', ''), nums[3]: ('accept', '')}
    rc = review_catalog(cat, str(tmp_path / 'rev.tsv'), dec)
    n = len(rows)
    want = {'all': n, 'accepted': 2, 'accepted_uncertain': 3, 'not_rejected': n - 1, 'uncertain': 1, 'rejected': 1,
            'unreviewed': n - 4, 'table': n}
    for inc, k in want.items():
        res = R.generate(rc, M51, str(tmp_path / ('r_%s.html' % inc)), include=inc, thumbs=False, max_rows=100000)
        assert res['selected'] == k, inc
    res = R.generate(rc, M51, str(tmp_path / 'o.html'), include='table', thumbs=False, only_numbers={nums[0], nums[5]})
    assert res['selected'] == 2


@needs_m51
def test_sort_max_rows_columns_and_errors(m51_cat, tmp_path):
    cat, _ = m51_cat
    out = str(tmp_path / 's.html')
    cols, rows = R.read_tsv(cat)
    res = R.generate(cat, M51, out, sort_by='MAG_AUTO', max_rows=5, thumbs=False, columns=['NUMBER', 'MAG_AUTO'])
    t = open(out, encoding='utf-8').read()
    assert res['listed'] == 5 and 'further candidates were not listed' in t
    first = re.findall(r"<tr class=\"cand st-none\" data-review=\"none\">\s*<td class=\"num\">[^<]*</td>\s*<td class=\"num\">(-?\d+\.\d+)</td>", t)[0]
    assert float(first) == min(float(r[cols.index('MAG_AUTO')]) for r in rows)
    with pytest.raises(ValueError):
        R.generate(cat, M51, out, sort_by='NOPE')
    with pytest.raises(ValueError):
        R.generate(cat, M51, out, columns=['NOPE'])
    with pytest.raises(ValueError):
        R.generate(cat, M51, out, include='bogus')


@needs_hudf
def test_hudf_report_with_real_catalog(tmp_path):
    from astropy.io import fits
    with fits.open(HUDF, memmap=True) as h:
        sub = np.array(h[0].data[1500:2100, 1500:2100], dtype='float32')
    crop = str(tmp_path / 'hudf_crop.fits')
    fits.writeto(crop, sub)
    cat = extract(crop, tmp_path, '--detect-thresh', '3')
    cols, rows = R.read_tsv(cat)
    assert len(rows) > 20
    nums = [r[cols.index('NUMBER')] for r in rows]
    rc = review_catalog(cat, str(tmp_path / 'rev.tsv'), {nums[0]: ('accept', 'star-like'), nums[1]: ('reject', 'artifact')})
    out = str(tmp_path / 'h.html')
    res = R.generate(rc, crop, out, include='accepted_uncertain', cutout_size=48, thumb_zoom=3, stretch='asinh',
                     png_dir=str(tmp_path / 'png'))
    t = open(out, encoding='utf-8').read()
    assert res['selected'] == 1 and res['thumbnails'] == 1
    th = thumbs_of(t)
    arr = png_array(th[0])
    assert arr.shape[:2] == (144, 144)
    assert int(arr[..., 0].max()) - int(arr[..., 0].min()) > 50            # the cut-out of a source is not flat
    pngs = os.listdir(str(tmp_path / 'png'))
    assert len(pngs) == 1 and pngs[0].startswith('cutout_')
    assert 'star-like' in t and 'artifact' not in t                        # the rejected one is not listed


def test_html_escaping_of_catalog_cells(tmp_path):
    p = tmp_path / 'c.tsv'
    p.write_text('NUMBER\tX_IMAGE\tY_IMAGE\tFOO\tREVIEW\tREVIEW_NOTE\n'
                 '1\t10\t10\t<script>alert(1)</script>\taccept\t"><img src=x onerror=alert(2)>\n')
    out = str(tmp_path / 'e.html')
    R.generate(str(p), '', out, thumbs=False, columns=['NUMBER', 'FOO'])
    t = open(out, encoding='utf-8').read()
    assert '<script>alert(1)' not in t and '&lt;script&gt;alert(1)' in t
    assert '<img src=x' not in t
    assert t.count('<script>') == 1                                        # only the report's own filter script


def test_cli_stdout_contract(tmp_path):
    p = tmp_path / 'c.tsv'
    p.write_text('NUMBER\tX_IMAGE\tY_IMAGE\tMAG_AUTO\n1\t10\t10\t20.1\n2\t12\t12\t21.0\n')
    out = str(tmp_path / 'x.html')
    r = subprocess.run([sys.executable, SCRIPT, '--catalog', str(p), '--output', out, '--no-thumbs'],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert r.returncode == 0, r.stderr
    f = r.stdout.decode().strip().split('\t')
    assert f[0] == 'REPORT' and f[1] == out and f[2] == '2' and f[3] == '2' and f[4] == '0' and f[5] == '-'
    r = subprocess.run([sys.executable, SCRIPT, '--catalog', str(tmp_path / 'none.tsv'), '--output', out],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert r.returncode == 1 and b'ERROR' in r.stderr


@pytest.mark.skipif(R.find_chrome() is None, reason='no headless Chrome/Chromium for the optional PDF')
def test_optional_pdf(tmp_path):
    p = tmp_path / 'c.tsv'
    p.write_text('NUMBER\tX_IMAGE\tY_IMAGE\tMAG_AUTO\n1\t10\t10\t20.1\n')
    res = R.generate(str(p), '', str(tmp_path / 'x.html'), thumbs=False, pdf=True)
    assert res['pdf'].endswith('.pdf') and open(res['pdf'], 'rb').read(5) == b'%PDF-'


def test_pdf_missing_chrome_is_a_warning(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(R, 'find_chrome', lambda: None)
    p = tmp_path / 'c.tsv'
    p.write_text('NUMBER\tX_IMAGE\tY_IMAGE\n1\t1\t1\n')
    res = R.generate(str(p), '', str(tmp_path / 'x.html'), thumbs=False, pdf=True)
    assert res['pdf'] == '-' and os.path.exists(str(tmp_path / 'x.html'))
    assert 'no Chrome' in capsys.readouterr().err
