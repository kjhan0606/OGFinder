#!/usr/bin/env python3
"""Spectra linked to catalog objects: 1D spectra, 2D spectra, IFU cubes; line fits and redshifts.

    spectra.py --task link --catalog TSV --work DIR [--spec-dir D --spec-pattern spec_{NUMBER}.fits | --link-file TSV]
    spectra.py --task fit  --catalog TSV --work DIR [same link options] [--snr-min 4 --z-column Z_SPEC --fit-line Ha ...]
    spectra.py --task plot --number N --catalog TSV --work DIR --out PNG [same link options]

link/fit follow the add_columns contract.  The link table (DIR/spectra_links.tsv: NUMBER FILE KIND X Y ROW) is written by every task; a user-supplied --link-file
(columns NUMBER FILE, optional KIND X Y ROW) overrides the file-name pattern.  fit also writes DIR/spectra_results.json (per object: redshift, lines).
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, spectra as sp, meta as ometa  # noqa: E402

COLUMNS = {
    'link': ['SP_FILE', 'SP_KIND', 'SP_OK'],
    'fit': ['SP_Z', 'SP_ZERR', 'SP_ZQ', 'SP_NLINES', 'SP_SNR', 'SP_LINE_FLUX', 'SP_LINE_FLUXERR', 'SP_LINE_FWHM_KMS', 'SP_LINE_EW'],
}
C_KMS = 299792.458


def detect_kind(path):
    if not path.lower().endswith(('.fits', '.fit', '.fits.gz', '.fz')):
        return '1d'
    from astropy.io import fits
    try:
        with fits.open(path) as hl:
            for h in hl:
                if h.data is not None and getattr(h.data, 'ndim', 0) == 3:
                    return 'cube'
            for h in hl:
                if hasattr(h, 'columns') and h.data is not None and len(h.columns) >= 2:
                    return '1d'
            for h in hl:
                if h.data is not None and getattr(h.data, 'ndim', 0) == 2:
                    return '2d'
    except Exception:
        return '1d'
    return '1d'


def build_links(a, rows):
    """-> list of dicts per catalog row: NUMBER, FILE (absolute or ''), KIND, X, Y, ROW, OK."""
    user = {}
    if a.link_file and os.path.isfile(a.link_file):
        cols, lr = tsvio.read_catalog(a.link_file)
        base = os.path.dirname(os.path.abspath(a.link_file))
        for r in lr:
            f = r.get('FILE', '')
            if f and not os.path.isabs(f):
                f = os.path.join(a.spec_dir or base, f) if a.spec_dir else os.path.join(base, f)
            user[str(r.get('NUMBER', '')).strip()] = dict(FILE=f, KIND=r.get('KIND', ''), X=tsvio.fnum(r.get('X')), Y=tsvio.fnum(r.get('Y')), ROW=tsvio.fnum(r.get('ROW')))
    links = []
    for r in rows:
        num = str(r['NUMBER']).strip()
        if user:
            u = user.get(num)
            d = dict(FILE=u['FILE'], KIND=u['KIND'], X=u['X'], Y=u['Y'], ROW=u['ROW']) if u else dict(FILE='', KIND='', X=float('nan'), Y=float('nan'), ROW=float('nan'))
        else:
            f = os.path.join(a.spec_dir, a.spec_pattern.replace('{NUMBER}', num)) if a.spec_dir else ''
            d = dict(FILE=f, KIND='', X=float('nan'), Y=float('nan'), ROW=float('nan'))
        d['NUMBER'] = num
        d['OK'] = int(bool(d['FILE']) and os.path.isfile(d['FILE']))
        if d['OK'] and (not d['KIND'] or a.kind != 'auto'):
            d['KIND'] = detect_kind(d['FILE']) if a.kind == 'auto' else a.kind
        # catalog positions are the default spaxel (cube) / row (2D) coordinates
        if not np.isfinite(d['X']):
            d['X'] = tsvio.fnum(r.get('X_IMAGE')) - 1 if a.cube_xy else float('nan')
        if not np.isfinite(d['Y']):
            d['Y'] = tsvio.fnum(r.get('Y_IMAGE')) - 1 if a.cube_xy else float('nan')
        if not np.isfinite(d['ROW']) and a.row_column and r.get(a.row_column, '') != '':
            d['ROW'] = tsvio.fnum(r.get(a.row_column))
        links.append(d)
    return links


def load_1d(a, lk, rows_by_num=None):
    """Spectrum of one linked object as dict(wave, flux, err, info) whatever the file kind."""
    wu = None if a.wave_unit == 'auto' else a.wave_unit
    k = lk['KIND']
    if k == '1d':
        s = sp.read_spectrum(lk['FILE'], wu)
        s['info'] = dict(kind='1d')
        return s
    if k == '2d':
        t = sp.read_2d(lk['FILE'], wu)
        ny = t['data'].shape[0]
        row = lk['ROW'] if np.isfinite(lk['ROW']) else (ny - 1) / 2.0
        sky = (a.sky_inner, a.sky_outer) if a.sky_outer > a.sky_inner else None
        f, e, inf = sp.extract_2d(t['data'], row, a.extract_halfwidth, t['err'], sky, a.extract_mode)
        return dict(wave=t['wave'], flux=f, err=e, info=dict(kind='2d', row=float(row), data=t['data'], sky=inf.get('sky'), profile_sigma=inf.get('profile_sigma')))
    if k == 'cube':
        c = sp.read_cube(lk['FILE'])
        nz, ny, nx = c['data'].shape
        x = lk['X'] if np.isfinite(lk['X']) else (nx - 1) / 2.0
        y = lk['Y'] if np.isfinite(lk['Y']) else (ny - 1) / 2.0
        f, e = sp.extract_cube(c['data'], x, y, a.aperture_radius, c['err'])
        return dict(wave=c['wave'], flux=f, err=e, info=dict(kind='cube', x=float(x), y=float(y), cube=c['data'], wave_cube=c['wave']))
    raise ValueError('unknown kind %r' % k)


def line_pick(res, name):
    for ft in res['lines']:
        if ft.get('ok') and ft['name'] == name:
            return ft
    return None


def task_fit(a, cols, rows, W):
    links = build_links(a, rows)
    out = []
    results = {}
    nz = 0
    for r, lk in zip(rows, links):
        d = {}
        if lk['OK']:
            try:
                s = load_1d(a, lk)
                err = s['err']
                if err is not None:
                    err = np.where(np.isfinite(err) & (err > 0), err, np.nan)
                zk = tsvio.fnum(r.get(a.z_column)) if a.z_column and a.z_column in r else float('nan')
                an = sp.analyse_spectrum(s['wave'], s['flux'], err, a.snr_min, a.kernel_sigma_px, a.cont_width, z_known=zk if np.isfinite(zk) and zk > 0 else None, zmax=a.zmax)
                z = an['redshift']
                d.update(SP_Z=z['z'] if z['quality'] > 0 or z['quality'] == 4 else None, SP_ZERR=z.get('z_err'), SP_ZQ=z['quality'], SP_NLINES=z.get('n_lines_fit', 0), SP_SNR=an['cont_snr_per_pix'])
                if z['quality'] == 0:
                    d['SP_Z'] = None; d['SP_ZERR'] = None
                ft = line_pick(an, a.fit_line)
                if ft is not None:
                    d.update(SP_LINE_FLUX=ft['flux'], SP_LINE_FLUXERR=ft['flux_err'], SP_LINE_FWHM_KMS=ft['fwhm'] / ft['mu'] * C_KMS, SP_LINE_EW=ft['ew'] / (1 + z['z']) if z['quality'] else None)
                nz += z['quality'] > 0
                results[lk['NUMBER']] = dict(kind=lk['KIND'], file=lk['FILE'], redshift={k: v for k, v in z.items()}, lines=an['lines'], noise=an['noise'], cont_snr=an['cont_snr_per_pix'], detections=an['detections'])
            except Exception as ex:                                  # a broken file must not stop the batch
                results[lk['NUMBER']] = dict(error=str(ex))
                sys.stderr.write('spectra: object %s: %s\n' % (lk['NUMBER'], ex))
        out.append((r['NUMBER'], d))
    return out, dict(n_linked=sum(l['OK'] for l in links), n_redshift=int(nz)), links, results


def task_link(a, cols, rows, W):
    links = build_links(a, rows)
    out = [(r['NUMBER'], dict(SP_FILE=os.path.basename(l['FILE']) if l['OK'] else None, SP_KIND=l['KIND'] if l['OK'] else None, SP_OK=l['OK'])) for r, l in zip(rows, links)]
    return out, dict(n_linked=sum(l['OK'] for l in links), n_missing=sum(1 - l['OK'] for l in links)), links, None


def plot_object(a, rows, links, number, out_png, results=None):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    lk = next((l for l in links if l['NUMBER'] == str(number)), None)
    if lk is None or not lk['OK']:
        raise SystemExit('spectra: object %s has no linked spectrum' % number)
    s = load_1d(a, lk)
    err = s['err']
    if err is not None:
        err = np.where(np.isfinite(err) & (err > 0), err, np.nan)
    row = next((r for r in rows if str(r['NUMBER']) == str(number)), {})
    zk = tsvio.fnum(row.get(a.z_column)) if a.z_column and a.z_column in row else float('nan')
    an = sp.analyse_spectrum(s['wave'], s['flux'], err, a.snr_min, a.kernel_sigma_px, a.cont_width, z_known=zk if np.isfinite(zk) and zk > 0 else None, zmax=a.zmax)
    z = an['redshift']
    info = s['info']
    extra = info['kind'] in ('2d', 'cube')
    fig = plt.figure(figsize=(8.0, 5.2 if extra else 3.6))
    gs = fig.add_gridspec(2 if extra else 1, 1, height_ratios=[1, 1.6] if extra else [1])
    ax = fig.add_subplot(gs[-1])
    w = s['wave']
    cont, _ = sp.continuum(w, s['flux'], a.cont_width)
    ax.plot(w, s['flux'], color='0.3', lw=0.6, drawstyle='steps-mid')
    ax.plot(w, cont, 'b-', lw=0.8)
    if err is not None:
        ax.plot(w, err, color='tab:orange', lw=0.5)
    for ft in an['lines']:
        if ft.get('ok'):
            ax.axvline(ft['mu'], color='tab:red', lw=0.6, alpha=0.7)
            ax.text(ft['mu'], ax.get_ylim()[1], ft['name'], rotation=90, va='top', ha='right', fontsize=6, color='tab:red')
    lo, hi = np.nanpercentile(s['flux'], [0.5, 99.7])
    ax.set_ylim(lo - 0.1 * (hi - lo), hi + 0.15 * (hi - lo))
    ax.set_xlabel('wavelength (A)'); ax.set_ylabel('flux')
    ax.set_title('object %s (%s)  z = %s  quality %d, %d lines' % (number, info['kind'], ('%.5f' % z['z']) if z['quality'] else 'n/a', z['quality'], z.get('n_lines_fit', 0)), fontsize=8)
    if info['kind'] == '2d':
        a2 = fig.add_subplot(gs[0])
        d = info['data']
        v1, v2 = np.nanpercentile(d, [5, 99.5])
        a2.imshow(d, aspect='auto', origin='lower', cmap='gray', vmin=v1, vmax=v2, extent=[w[0], w[-1], -0.5, d.shape[0] - 0.5])
        a2.axhline(info['row'] - a.extract_halfwidth, color='y', lw=0.5); a2.axhline(info['row'] + a.extract_halfwidth, color='y', lw=0.5)
        a2.set_ylabel('row'); a2.set_title('2D spectrum (%s extraction, half width %.1f)' % (a.extract_mode, a.extract_halfwidth), fontsize=7)
    elif info['kind'] == 'cube':
        a2 = fig.add_subplot(gs[0])
        wl = sp.white_light(info['cube'])
        a2.imshow(wl, origin='lower', cmap='gray')
        th = np.linspace(0, 2 * math.pi, 60)
        a2.plot(info['x'] + a.aperture_radius * np.cos(th), info['y'] + a.aperture_radius * np.sin(th), 'y-', lw=0.8)
        a2.set_title('white-light image, aperture r=%.1f spaxels' % a.aperture_radius, fontsize=7)
    fig.tight_layout()
    fig.savefig(out_png, dpi=80)
    plt.close(fig)
    return an


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--task', required=True, choices=['link', 'fit', 'plot'])
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--work', default='.')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--spec-dir', default='')
    ap.add_argument('--spec-pattern', default='spec_{NUMBER}.fits')
    ap.add_argument('--link-file', default='')
    ap.add_argument('--kind', default='auto', choices=['auto', '1d', '2d', 'cube'])
    ap.add_argument('--wave-unit', default='auto', choices=['auto', 'angstrom', 'nm', 'um'])
    ap.add_argument('--row-column', default='')
    ap.add_argument('--cube-xy', type=int, default=0)
    ap.add_argument('--extract-mode', default='boxcar', choices=['boxcar', 'optimal'])
    ap.add_argument('--extract-halfwidth', type=float, default=3.0)
    ap.add_argument('--sky-inner', type=float, default=0.0)
    ap.add_argument('--sky-outer', type=float, default=0.0)
    ap.add_argument('--aperture-radius', type=float, default=2.5)
    ap.add_argument('--snr-min', type=float, default=4.0)
    ap.add_argument('--kernel-sigma-px', type=float, default=2.0)
    ap.add_argument('--cont-width', type=int, default=101)
    ap.add_argument('--z-column', default='')
    ap.add_argument('--zmax', type=float, default=7.0)
    ap.add_argument('--fit-line', default='Ha')
    ap.add_argument('--number', default='')
    ap.add_argument('--out', default='')
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    cols, rows = tsvio.read_catalog(a.catalog)
    if a.task == 'plot':
        links = build_links(a, rows)
        plot_object(a, rows, links, a.number, a.out or os.path.join(a.work, 'spectra_view.png'))
        return 0
    if a.task == 'fit':
        out, summ, links, results = task_fit(a, cols, rows, a.work)
        with open(os.path.join(a.work, 'spectra_results.json'), 'w') as fh:
            json.dump(results, fh, default=lambda o: None)
    else:
        out, summ, links, _ = task_link(a, cols, rows, a.work)
    tsvio.write_table(os.path.join(a.work, 'spectra_links.tsv'), ['NUMBER', 'FILE', 'KIND', 'X', 'Y', 'ROW', 'OK'], links)
    if a.meta_out:
        ometa.update(a.meta_out, 'spectra', dict(summ, task=a.task), nrows=len(rows))
    tsvio.write_columns(COLUMNS[a.task], out)
    sys.stderr.write('spectra %s: %s\n' % (a.task, json.dumps(summ)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
