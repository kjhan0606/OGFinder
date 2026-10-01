#!/usr/bin/env python3
"""OGFinder candidate report: one self-contained HTML file (optionally a PDF) from a catalog TSV and the image.

    python3 report.py --catalog cat.tsv --image image.fits --output report.html
                      [--provenance prov.json] [--include table|all|accepted|accepted_uncertain|not_rejected|uncertain|rejected|unreviewed]
                      [--title T] [--author A] [--sort-by COL|-COL] [--max-rows N] [--cutout-size PX] [--thumb-zoom K]
                      [--stretch zscale|asinh|linear] [--no-thumbs] [--no-marker] [--columns A,B,...] [--no-all-columns]
                      [--png-dir DIR] [--pdf] [--table-filter TEXT] [--only-numbers FILE]

The catalog is the table the GUI shows (TSV, header line first).  Review state is read from the columns REVIEW
(accept | reject | uncertain | empty), REVIEW_NOTE and REVIEW_TIME if present.  Thumbnails are PNG cut-outs of the image around
X_IMAGE/Y_IMAGE (FITS 1-based pixel coordinates), embedded as data: URIs, so the HTML has no external files.  The provenance JSON
is written by the GUI (plugins/report/report.tcl): session steps, versions, plugin parameters, input checksums.

stdout (one line, tab separated):  REPORT <html> <n_candidates_listed> <n_total_after_filter> <n_thumbnails> <pdf-or-'-'>
Only numpy, astropy and PIL are needed.  Every value from the catalog or the provenance is HTML-escaped.
"""
import argparse
import base64
import hashlib
import html
import io
import json
import os
import platform
import shutil
import subprocess
import sys
import time
import zlib

REVIEW_VALUES = ('accept', 'reject', 'uncertain')
INCLUDE_CHOICES = ('table', 'all', 'accepted', 'accepted_uncertain', 'not_rejected', 'uncertain', 'rejected', 'unreviewed')
PREFERRED = ['NUMBER', 'kind', 'X_IMAGE', 'Y_IMAGE', 'ALPHA_J2000', 'DELTA_J2000', 'MAG_AUTO', 'MAGERR_AUTO', 'FLUX_RADIUS',
             'FWHM_IMAGE', 'ELLIPTICITY', 'CLASS_STAR', 'n', 'score', 'rate', 'host_id']
# the standard ds9_sextract columns: anything else in the catalog is an analysis result and is shown as an extra column
SEXTRACTOR_COLUMNS = set('''NUMBER X_IMAGE Y_IMAGE ALPHA_J2000 DELTA_J2000 MAG_AUTO MAG_ISOCOR MAG_APER FLUX_AUTO FLUXERR_AUTO
FLUXERR_APER MAGERR_AUTO MAGERR_APER FLUX_APER_2 FLUX_APER_3 FLUX_APER_5 FLUXERR_APER_2 FLUXERR_APER_3 FLUXERR_APER_5 MAG_APER_2
MAG_APER_3 MAG_APER_5 FLUX_RADIUS A_IMAGE B_IMAGE THETA_IMAGE ELLIPTICITY KRON_RADIUS FWHM_IMAGE ISO_RADIUS NPIX_ISO MU_MAX
MU_THRESHOLD CLASS_STAR FLAGS'''.split())
MAX_EXTRA_COLUMNS = 10


# ------------------------------------------------------------------------------------------------ catalog
def read_tsv(path):
    """-> (columns, rows); rows are lists of str padded to len(columns)"""
    with open(path, 'r', encoding='utf-8', errors='replace') as f:
        lines = [ln.rstrip('\r\n') for ln in f]
    lines = [ln for ln in lines if ln.strip()]
    if not lines:
        raise ValueError('empty catalog: %s' % path)
    cols = [c.strip() for c in lines[0].split('\t')]
    rows = []
    for ln in lines[1:]:
        r = [v.strip() for v in ln.split('\t')]
        r += [''] * (len(cols) - len(r))
        rows.append(r[:len(cols)])
    return cols, rows


def review_of(row, idx):
    v = row[idx['REVIEW']].strip().lower() if 'REVIEW' in idx else ''
    return v if v in REVIEW_VALUES else ''


def select_rows(cols, rows, include, only_numbers=None):
    idx = {c: i for i, c in enumerate(cols)}
    if only_numbers is not None:
        if 'NUMBER' not in idx:
            raise ValueError('--only-numbers needs a NUMBER column in the catalog')
        ni = idx['NUMBER']
        rows = [r for r in rows if r[ni] in only_numbers]
    keep = {
        'table': lambda s: True, 'all': lambda s: True,
        'accepted': lambda s: s == 'accept',
        'accepted_uncertain': lambda s: s in ('accept', 'uncertain'),
        'not_rejected': lambda s: s != 'reject',
        'uncertain': lambda s: s == 'uncertain',
        'rejected': lambda s: s == 'reject',
        'unreviewed': lambda s: s == '',
    }[include]
    return [r for r in rows if keep(review_of(r, idx))]


def sort_rows(cols, rows, sort_by):
    if not sort_by:
        return rows
    desc = sort_by.startswith('-')
    name = sort_by.lstrip('-+')
    if name not in cols:
        raise ValueError('sort column %r not in the catalog (columns: %s)' % (name, ', '.join(cols[:12])))
    i = cols.index(name)

    def num(v):
        try:
            return (0, float(v))
        except ValueError:
            return (1, v)
    return sorted(rows, key=lambda r: num(r[i]), reverse=desc)


def choose_columns(cols, requested):
    if requested:
        miss = [c for c in requested if c not in cols]
        if miss:
            raise ValueError('requested column(s) not in the catalog: %s' % ', '.join(miss))
        shown = list(requested)
    else:
        shown = [c for c in PREFERRED if c in cols]
        extras = [c for c in cols if c not in SEXTRACTOR_COLUMNS and c not in PREFERRED and not c.startswith('REVIEW')
                  and c not in ('kind', 'id', 'ra', 'dec', 'mag', '_pos')]
        shown += extras[:MAX_EXTRA_COLUMNS]
    return shown


# ------------------------------------------------------------------------------------------------ image / thumbnails
def load_image(path):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        for h in hdul:
            d = h.data
            if d is not None and getattr(d, 'ndim', 0) >= 2:
                import numpy as np
                while d.ndim > 2:
                    d = d[0]
                return np.asarray(d, dtype='float32'), h.header
    raise ValueError('no 2-D image in %s' % path)


def cutout(data, x, y, size):
    """size x size array centred on the 1-based FITS pixel (x, y); NaN outside the image"""
    import numpy as np
    ny, nx = data.shape
    half = size // 2
    cx = int(round(x)) - 1
    cy = int(round(y)) - 1
    out = np.full((size, size), np.nan, dtype='float32')
    x0, y0 = cx - half, cy - half
    xs0, ys0 = max(x0, 0), max(y0, 0)
    xs1, ys1 = min(x0 + size, nx), min(y0 + size, ny)
    if xs1 <= xs0 or ys1 <= ys0:
        return out
    out[ys0 - y0:ys1 - y0, xs0 - x0:xs1 - x0] = data[ys0:ys1, xs0:xs1]
    return out


def global_limits(data, mode):
    """display limits from the whole image (so all thumbnails share one stretch)"""
    import numpy as np
    sub = data[::max(1, data.shape[0] // 400), ::max(1, data.shape[1] // 400)]
    sub = sub[np.isfinite(sub)]
    if sub.size == 0:
        return 0.0, 1.0
    if mode == 'linear':
        return float(np.percentile(sub, 1)), float(np.percentile(sub, 99.5))
    try:
        from astropy.visualization import ZScaleInterval
        lo, hi = ZScaleInterval().get_limits(sub)
    except Exception:
        lo, hi = np.percentile(sub, 5), np.percentile(sub, 95)
    if mode == 'asinh':
        hi = float(np.percentile(sub, 99.9)) if hi <= lo else max(float(hi), float(np.percentile(sub, 99.9)))
    if not hi > lo:
        hi = lo + 1.0
    return float(lo), float(hi)


def to_uint8(arr, lo, hi, mode):
    import numpy as np
    a = np.nan_to_num((arr - lo) / (hi - lo), nan=0.0, posinf=1.0, neginf=0.0)
    a = np.clip(a, 0.0, 1.0)
    if mode == 'asinh':
        a = np.arcsinh(10.0 * a) / np.arcsinh(10.0)
    return (a * 255.0 + 0.5).astype('uint8')


def render_png(a8, zoom=2, marker=True):
    from PIL import Image, ImageDraw
    im = Image.fromarray(a8[::-1], mode='L')              # FITS y points up
    if zoom > 1:
        im = im.resize((im.width * zoom, im.height * zoom), Image.NEAREST)
    im = im.convert('RGB')
    if marker:
        d = ImageDraw.Draw(im)
        c = im.width // 2
        r = max(6, im.width // 8)
        d.ellipse([c - r, c - r, c + r, c + r], outline=(255, 210, 0))
    buf = io.BytesIO()
    im.save(buf, format='PNG', optimize=True)
    return buf.getvalue()


# ------------------------------------------------------------------------------------------------ provenance
def load_provenance(path):
    if not path:
        return {}
    with open(path, 'r', encoding='utf-8') as f:
        return json.load(f)


def file_facts(path):
    if not path or not os.path.exists(path):
        return {'path': path or '', 'exists': False}
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 22), b''):
            h.update(chunk)
    st = os.stat(path)
    return {'path': os.path.abspath(path), 'exists': True, 'bytes': st.st_size, 'sha256': h.hexdigest(),
            'mtime': time.strftime('%Y-%m-%dT%H:%M:%S%z', time.localtime(st.st_mtime))}


def python_versions():
    out = {'python': platform.python_version(), 'platform': platform.platform()}
    for mod in ('numpy', 'astropy', 'PIL', 'sep', 'scipy'):
        try:
            m = __import__(mod)
            out[mod] = getattr(m, '__version__', getattr(m, 'PILLOW_VERSION', '?'))
        except Exception:
            out[mod] = 'not installed'
    return out


# ------------------------------------------------------------------------------------------------ HTML
CSS = '''
body{font-family:Helvetica,Arial,sans-serif;margin:24px;color:#1b1b1b;background:#fff}
h1{margin:0 0 4px 0;font-size:24px} h2{margin-top:28px;border-bottom:1px solid #bbb;padding-bottom:3px;font-size:18px}
.meta{color:#555;font-size:13px} table{border-collapse:collapse;font-size:12px} th,td{border:1px solid #ccc;padding:3px 6px;vertical-align:top}
th{background:#eee;text-align:left} td.num{text-align:right;font-variant-numeric:tabular-nums} img.thumb{display:block;image-rendering:pixelated}
tr.st-accept td.rv{background:#cfeccf} tr.st-reject td.rv{background:#f3c9c9} tr.st-uncertain td.rv{background:#fbe7b0}
.counts span{display:inline-block;margin-right:14px;padding:2px 8px;border:1px solid #ccc;border-radius:3px;font-size:13px}
details{margin:2px 0} summary{cursor:pointer;color:#345} pre{white-space:pre-wrap;word-break:break-all;margin:2px 0;font-size:11px}
.note{white-space:pre-wrap;max-width:260px} .small{font-size:11px;color:#666} .warn{color:#8a4b00}
@media print{.noprint{display:none} body{margin:8px}}
'''

FILTER_JS = '''
function ogfFilter(v){var rows=document.querySelectorAll('tr.cand');for(var i=0;i<rows.length;i++){
var s=rows[i].getAttribute('data-review')||'none';rows[i].style.display=(v==='all'||v===s)?'':'none';}}
'''


def esc(v):
    return html.escape(str(v), quote=True)


def fmt_cell(v):
    return esc(v)


def build_html(ctx):
    """ctx: dict produced by generate(); returns the HTML text"""
    o = []
    A = o.append
    A('<!DOCTYPE html>\n<html lang="en"><head><meta charset="utf-8">')
    A('<title>%s</title><style>%s</style></head><body>' % (esc(ctx['title']), CSS))
    A('<h1>%s</h1>' % esc(ctx['title']))
    A('<div class="meta">generated %s%s &middot; OGFinder report v%s &middot; image <code>%s</code></div>' % (
        esc(ctx['generated']), (' &middot; prepared by %s' % esc(ctx['author'])) if ctx['author'] else '',
        esc(ctx['version']), esc(os.path.basename(ctx['image'])) if ctx['image'] else '(none)'))
    # -------------------------------------------------------------------- summary
    c = ctx['counts']
    A('<h2>Summary</h2><div class="counts">')
    A('<span>catalog rows: <b>%d</b></span><span>in this report: <b>%d</b>%s</span>' % (
        ctx['n_catalog'], ctx['n_selected'], ' (first %d listed)' % ctx['n_listed'] if ctx['n_listed'] < ctx['n_selected'] else ''))
    A('<span>accepted: <b>%d</b></span><span>rejected: <b>%d</b></span><span>uncertain: <b>%d</b></span>'
      '<span>unreviewed: <b>%d</b></span>' % (c['accept'], c['reject'], c['uncertain'], c['none']))
    A('</div><div class="meta">row selection: <b>%s</b>%s</div>' % (esc(ctx['include']), (
        ' &middot; sorted by %s' % esc(ctx['sort_by'])) if ctx['sort_by'] else ''))
    if ctx['table_filter']:
        A('<div class="meta">table filter when exported: %s</div>' % esc(ctx['table_filter']))
    if ctx['warnings']:
        A('<div class="warn">%s</div>' % '<br>'.join(esc(w) for w in ctx['warnings']))
    # -------------------------------------------------------------------- candidates
    A('<h2>Candidates</h2><div class="noprint meta">show: ')
    for v, lab in (('all', 'all'), ('accept', 'accepted'), ('uncertain', 'uncertain'), ('reject', 'rejected'), ('none', 'unreviewed')):
        A('<a href="#" onclick="ogfFilter(\'%s\');return false">%s</a> ' % (v, lab))
    A('</div>')
    shown = ctx['shown_columns']
    A('<table><tr>%s%s<th>review</th><th>note</th></tr>' % ('<th>cut-out</th>' if ctx['thumbs'] else '',
                                                              ''.join('<th>%s</th>' % esc(c) for c in shown)))
    idx = ctx['index']
    for n, row in enumerate(ctx['rows']):
        st = review_of(row, idx)
        A('<tr class="cand st-%s" data-review="%s">' % (st or 'none', st or 'none'))
        if ctx['thumbs']:
            png = ctx['png'].get(n)
            if png is None:
                A('<td class="small">no cut-out</td>')
            else:
                z = ctx['thumb_px']
                A('<td><img class="thumb" width="%d" height="%d" alt="cut-out %s" src="data:image/png;base64,%s"></td>' % (
                    z, z, esc(row[idx['NUMBER']]) if 'NUMBER' in idx else n, base64.b64encode(png).decode('ascii')))
        for cname in shown:
            v = row[idx[cname]]
            A('<td%s>%s</td>' % (' class="num"' if _isnum(v) else '', fmt_cell(v)))
        A('<td class="rv">%s%s</td>' % (esc(st or '-'), (('<div class="small">%s</div>' % esc(row[idx['REVIEW_TIME']]))
                                                       if 'REVIEW_TIME' in idx and row[idx['REVIEW_TIME']] else '')))
        note = row[idx['REVIEW_NOTE']] if 'REVIEW_NOTE' in idx else ''
        A('<td class="note">%s</td>' % esc(note))
        A('</tr>')
        if ctx['all_columns']:
            rest = [(c, row[i]) for c, i in idx.items() if c not in shown and not c.startswith('REVIEW')]
            if rest:
                A('<tr class="cand st-%s" data-review="%s"><td colspan="%d"><details><summary>all measurements of %s</summary>'
                  '<pre>%s</pre></details></td></tr>' % (
                      st or 'none', st or 'none', len(shown) + 3 + (0 if ctx['thumbs'] else -1),
                      esc(row[idx['NUMBER']]) if 'NUMBER' in idx else n,
                      esc('  '.join('%s=%s' % (c, v) for c, v in rest))))
    A('</table>')
    if ctx['n_listed'] < ctx['n_selected']:
        A('<div class="small">%d further candidates were not listed (--max-rows %d); the full table is in the catalog file.</div>' % (
            ctx['n_selected'] - ctx['n_listed'], ctx['max_rows']))
    # -------------------------------------------------------------------- review history
    hist = ctx['review_history']
    A('<h2>Review history (from the session log)</h2>')
    if hist:
        A('<table><tr><th>#</th><th>time</th><th>candidates</th><th>decision</th><th>note</th></tr>')
        for h in hist:
            A('<tr><td class="num">%s</td><td>%s</td><td>%s</td><td>%s</td><td class="note">%s</td></tr>' % (
                esc(h['seq']), esc(h['time']), esc(h['numbers']), esc(h['status']), esc(h['note'])))
        A('</table>')
    else:
        A('<div class="small">no review steps in the session log (the decisions above come from the catalog columns REVIEW / '
          'REVIEW_NOTE / REVIEW_TIME).</div>')
    # -------------------------------------------------------------------- provenance
    A('<h2>Processing provenance</h2>')
    prov = ctx['prov']
    A('<h3>Inputs</h3><table><tr><th>file</th><th>bytes</th><th>sha256</th><th>modified</th></tr>')
    for f in ctx['files']:
        if f['exists']:
            A('<tr><td><code>%s</code></td><td class="num">%d</td><td><code>%s</code></td><td>%s</td></tr>' % (
                esc(f['path']), f['bytes'], esc(f['sha256']), esc(f['mtime'])))
        else:
            A('<tr><td><code>%s</code></td><td colspan="3">not found</td></tr>' % esc(f['path']))
    A('</table>')
    A('<h3>Software</h3><table>')
    soft = [('OGFinder (git HEAD)', prov.get('git_head', '')), ('sources modified since HEAD', _yn(prov.get('git_dirty_sources', ''))),
            ('SAOImageDS9 version', prov.get('ds9_version', '')), ('Python used by the GUI', prov.get('python', '')),
            ('ds9 executable', prov.get('ds9_exe', '')), ('session started', prov.get('session_started', ''))]
    soft += list(ctx['versions'].items())
    for k, v in soft:
        A('<tr><th>%s</th><td>%s</td></tr>' % (esc(k), esc(v)))
    A('</table>')
    steps = prov.get('steps', [])
    A('<h3>Session log summary</h3>')
    if steps:
        cnt = {}
        for s in steps:
            cls = 'failed' if str(s.get('failed', '0')) not in ('0', '', 'False') else s.get('class', '?')
            cnt[cls] = cnt.get(cls, 0) + 1
        A('<div class="meta">%d recorded steps: %s</div>' % (len(steps), ', '.join('%s %d' % (esc(k), v) for k, v in sorted(cnt.items()))))
        A('<table><tr><th>#</th><th>class</th><th>step</th><th>title</th><th>seconds</th><th>command</th></tr>')
        for s in steps:
            ms = s.get('ms', '')
            sec = '%.1f' % (float(ms) / 1000.0) if _isnum(ms) else ''
            argv = s.get('argv', [])
            cmd = ' '.join(str(a) for a in argv[1:]) if argv else ''
            if len(cmd) > 400:
                cmd = cmd[:400] + ' ...'
            A('<tr><td class="num">%s</td><td>%s</td><td>%s</td><td>%s</td><td class="num">%s</td><td><pre>%s</pre></td></tr>' % (
                esc(s.get('seq', '')), esc(s.get('class', '')), esc(s.get('step', '')), esc(s.get('title', '')), esc(sec), esc(cmd)))
        A('</table>')
    else:
        A('<div class="small">no provenance file was given: this report lists no session steps.</div>')
    A('<h3>Parameters</h3>')
    pp = prov.get('plugin_params', {})
    if pp:
        A('<table><tr><th>plugin</th><th>parameters</th></tr>')
        for pid, d in pp.items():
            A('<tr><td>%s</td><td><pre>%s</pre></td></tr>' % (esc(pid), esc('  '.join('%s=%s' % (k, v) for k, v in d.items()))))
        A('</table>')
    ex = prov.get('extraction_params', {})
    if ex:
        A('<div class="meta">extraction / analysis parameters (catalog panel):</div><pre>%s</pre>' % esc(
            '  '.join('%s=%s' % (k, v) for k, v in ex.items())))
    if not pp and not ex:
        A('<div class="small">no parameters recorded.</div>')
    A('<h3>Report options</h3><pre>%s</pre>' % esc(json.dumps(ctx['options'], indent=1, sort_keys=True)))
    A('<div class="small">This file is self-contained (no external resources).  The step list and parameters are what the GUI '
      'recorded; export the session as a Python script (Workflow &gt; Session) to re-run the AUTO steps on other data.</div>')
    A('<script>%s</script></body></html>' % FILTER_JS)
    return '\n'.join(o)


def _yn(v):
    return {'1': 'yes', '0': 'no', 1: 'yes', 0: 'no', True: 'yes', False: 'no'}.get(v, v)


def _isnum(v):
    try:
        float(v)
        return True
    except (TypeError, ValueError):
        return False


def review_history(prov):
    out = []
    for s in prov.get('steps', []):
        if not str(s.get('step', '')).startswith('review.'):
            continue
        pl = s.get('payload', {}) or {}
        nums = pl.get('numbers', [])
        out.append({'seq': s.get('seq', ''), 'time': s.get('wall', ''), 'numbers': ', '.join(str(n) for n in nums) if isinstance(nums, list) else nums,
                    'status': pl.get('status', ''), 'note': pl.get('note', ''), 'step': s.get('step', '')})
    return out


# ------------------------------------------------------------------------------------------------ driver
def find_chrome():
    for n in ('google-chrome', 'chromium', 'chromium-browser', 'chrome', 'microsoft-edge'):
        p = shutil.which(n)
        if p:
            return p
    return None


def make_pdf(html_path, pdf_path):
    """returns (ok, message).  Needs a headless Chrome / Chromium / Edge; nothing else is tried."""
    chrome = find_chrome()
    if not chrome:
        return False, 'no Chrome/Chromium found on PATH: open the HTML in a browser and print it to PDF'
    cmd = [chrome, '--headless', '--no-sandbox', '--disable-gpu', '--no-pdf-header-footer', '--print-to-pdf=%s' % pdf_path,
           'file://' + os.path.abspath(html_path)]
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=120)
    except Exception as e:
        return False, 'PDF conversion failed: %s' % e
    if os.path.exists(pdf_path) and os.path.getsize(pdf_path) > 0:
        return True, ''
    return False, 'PDF conversion failed (rc=%s): %s' % (p.returncode, p.stderr.decode('utf-8', 'replace')[-200:])


def generate(catalog, image, output, provenance=None, include='table', title='OGFinder candidate report', author='',
             sort_by='', max_rows=500, cutout_size=64, thumb_zoom=2, stretch='zscale', thumbs=True, marker=True,
             columns=None, all_columns=True, png_dir=None, pdf=False, table_filter='', only_numbers=None):
    if include not in INCLUDE_CHOICES:
        raise ValueError('include must be one of %s' % ', '.join(INCLUDE_CHOICES))
    cols, rows = read_tsv(catalog)
    idx_all = {c: i for i, c in enumerate(cols)}
    prov = load_provenance(provenance)
    n_catalog = len(rows)
    counts = {'accept': 0, 'reject': 0, 'uncertain': 0, 'none': 0}
    for r in rows:
        counts[review_of(r, idx_all) or 'none'] += 1
    sel = select_rows(cols, rows, include, only_numbers)
    sel = sort_rows(cols, sel, sort_by)
    listed = sel[:max_rows] if max_rows and max_rows > 0 else sel
    warnings = []
    shown = choose_columns(cols, columns)
    png = {}
    data = None
    if thumbs:
        if not image or not os.path.exists(image):
            warnings.append('image %r not found: no cut-outs.' % image)
        elif 'X_IMAGE' not in idx_all or 'Y_IMAGE' not in idx_all:
            warnings.append('the table has no X_IMAGE/Y_IMAGE columns: no cut-outs.')
        else:
            data, _hdr = load_image(image)
            lo, hi = global_limits(data, stretch)
            for n, r in enumerate(listed):
                try:
                    x = float(r[idx_all['X_IMAGE']])
                    y = float(r[idx_all['Y_IMAGE']])
                except ValueError:
                    continue
                png[n] = render_png(to_uint8(cutout(data, x, y, cutout_size), lo, hi, stretch), thumb_zoom, marker)
    if png_dir and png:
        os.makedirs(png_dir, exist_ok=True)
        for n, b in png.items():
            nm = listed[n][idx_all['NUMBER']] if 'NUMBER' in idx_all else str(n)
            safe = ''.join(ch if ch.isalnum() or ch in '._-' else '_' for ch in nm)
            with open(os.path.join(png_dir, 'cutout_%s.png' % safe), 'wb') as f:
                f.write(b)
    options = {'only_numbers': (len(only_numbers) if only_numbers is not None else 'all rows'), 'include': include, 'sort_by': sort_by, 'max_rows': max_rows, 'cutout_size': cutout_size, 'thumb_zoom': thumb_zoom,
               'stretch': stretch, 'thumbnails': bool(thumbs), 'marker': bool(marker), 'columns': columns or 'auto',
               'all_columns': bool(all_columns), 'catalog_file': os.path.abspath(catalog)}
    ctx = dict(title=title, author=author, generated=time.strftime('%Y-%m-%d %H:%M:%S %z'), version='1', image=image or '',
               counts=counts, n_catalog=n_catalog, n_selected=len(sel), n_listed=len(listed), include=include, sort_by=sort_by,
               table_filter=table_filter or prov.get('table_filter', ''), warnings=warnings, shown_columns=shown, index=idx_all,
               rows=listed, png=png, thumbs=bool(thumbs and png), thumb_px=cutout_size * thumb_zoom, all_columns=all_columns,
               review_history=review_history(prov), prov=prov, versions=python_versions(), max_rows=max_rows, options=options,
               files=[file_facts(image), file_facts(catalog)])
    text = build_html(ctx)
    with open(output, 'w', encoding='utf-8') as f:
        f.write(text)
    pdf_name = '-'
    if pdf:
        pdf_path = os.path.splitext(output)[0] + '.pdf'
        ok, msg = make_pdf(output, pdf_path)
        if ok:
            pdf_name = pdf_path
        else:
            print('WARNING: ' + msg, file=sys.stderr)
    return dict(html=output, listed=len(listed), selected=len(sel), thumbnails=len(png), pdf=pdf_name, counts=counts,
                crc=zlib.crc32(text.encode('utf-8')) & 0xffffffff)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--image', default='')
    ap.add_argument('--output', required=True)
    ap.add_argument('--provenance', default='')
    ap.add_argument('--include', default='table', choices=INCLUDE_CHOICES)
    ap.add_argument('--title', default='OGFinder candidate report')
    ap.add_argument('--author', default='')
    ap.add_argument('--sort-by', default='')
    ap.add_argument('--max-rows', type=int, default=500)
    ap.add_argument('--cutout-size', type=int, default=64)
    ap.add_argument('--thumb-zoom', type=int, default=2)
    ap.add_argument('--stretch', default='zscale', choices=('zscale', 'asinh', 'linear'))
    ap.add_argument('--no-thumbs', action='store_true')
    ap.add_argument('--no-marker', action='store_true')
    ap.add_argument('--columns', default='')
    ap.add_argument('--no-all-columns', action='store_true')
    ap.add_argument('--png-dir', default='')
    ap.add_argument('--pdf', action='store_true')
    ap.add_argument('--table-filter', default='', help='text describing the table filter that was active (shown in the report)')
    ap.add_argument('--only-numbers', default='', help='file with one NUMBER per line: list only these catalog rows (the rows the table shows)')
    a = ap.parse_args(argv)
    only = None
    if a.only_numbers:
        with open(a.only_numbers, 'r', encoding='utf-8') as f:
            only = set(ln.strip() for ln in f if ln.strip())
    try:
        r = generate(a.catalog, a.image, a.output, provenance=a.provenance or None, include=a.include, title=a.title,
                     author=a.author, sort_by=a.sort_by, max_rows=a.max_rows, cutout_size=a.cutout_size, thumb_zoom=a.thumb_zoom,
                     stretch=a.stretch, thumbs=not a.no_thumbs, marker=not a.no_marker,
                     columns=[c for c in a.columns.split(',') if c] or None, all_columns=not a.no_all_columns,
                     png_dir=a.png_dir or None, pdf=a.pdf, table_filter=a.table_filter, only_numbers=only)
    except (ValueError, OSError) as e:
        print('ERROR: %s' % e, file=sys.stderr)
        return 1
    print('REPORT\t%s\t%d\t%d\t%d\t%s' % (r['html'], r['listed'], r['selected'], r['thumbnails'], r['pdf']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
