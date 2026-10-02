#!/usr/bin/env python3
"""Cross-match the current catalog with a local table (FITS / CSV / TSV / VOTable) or with a TAP service (through the AI-bridge tap_query adapter).

    xmatch.py --catalog TSV --work DIR --source file --ref-file TABLE [--radius-arcsec 1.0] [--ref-ra-col ra --ref-dec-col dec] [--copy-columns a,b,c,d]
    xmatch.py --catalog TSV --work DIR --source tap --allow-network --tap-url URL --tap-adql "SELECT ... {ra} {dec} {radius_deg}" [...]

stdout follows the add_columns contract: NUMBER XM_SEP XM_N XM_FLAG XM_DRA XM_DDEC XM_ID XM_V1..XM_V4 (separation in arcsec or pixels; flags in ogfkit/xmatch.py).
Files in DIR: xmatch_summary.json, xmatch_pairs.tsv (every pair within the radius), xmatch_tap.csv (the raw TAP answer).
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, xmatch as xm, meta as ometa  # noqa: E402

COLUMNS = ['XM_SEP', 'XM_N', 'XM_FLAG', 'XM_DRA', 'XM_DDEC', 'XM_ID', 'XM_V1', 'XM_V2', 'XM_V3', 'XM_V4']
DEFAULT_ADQL = ("SELECT TOP 50000 source_id, ra, dec, phot_g_mean_mag, parallax FROM gaiadr3.gaia_source "
                "WHERE 1=CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra}, {dec}, {radius_deg}))")


def fetch_tap(a, ra, dec, radius_deg, work):
    """One cone query over the whole field through ai_bridge's tap_query adapter -> (cols, data, info)."""
    if not a.allow_network:
        raise SystemExit('xmatch: TAP queries need the network: enable "Allow network access" (--allow-network)')
    from ai_bridge import adapters
    prof = {'name': 'xmatch_tap', 'task': 'generic', 'transport': 'tap_query', 'base_url': a.tap_url, 'auth': {'scheme': 'none'},
            'params': {'radius_arcsec': radius_deg * 3600.0}, 'request': {'adql': a.tap_adql or DEFAULT_ADQL, 'format': 'csv'}, 'response': {'fields': {}},
            'timeout_s': a.tap_timeout}
    ad = adapters.make_adapter(prof, {'params': dict(prof['params']), 'bands': []})
    rec = {'id': 'field', 'ra': ra, 'dec': dec, 'x': 0, 'y': 0, 'mags': {}, 'mag_errs': {}, 'catalog_row': {}}
    prepared = ad.prepare([rec], None)
    obj, meta = ad.send(prepared, a.tap_timeout)
    rows = obj['rows']
    cols = obj['columns']
    raw = {}
    with open(os.path.join(work, 'xmatch_tap.csv'), 'w') as fh:
        fh.write(','.join(cols) + '\n')
        for r in rows:
            fh.write(','.join(str(r.get(c, '')) for c in cols) + '\n')
    data = {}
    for c in cols:
        vals = [r.get(c, '') for r in rows]
        raw[c] = vals
        try:
            data[c] = np.array([float(v) if v not in ('', None) else np.nan for v in vals])
        except ValueError:
            data[c] = np.array(vals, dtype=object)
    return cols, data, raw, dict(n_rows=len(rows), http=meta.get('http_status'), url=a.tap_url, truncated=bool(a.tap_limit and len(rows) >= a.tap_limit))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--work', default='.')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--source', default='file', choices=['file', 'tap'])
    ap.add_argument('--ref-file', default='')
    ap.add_argument('--ref-format', default='auto', choices=['auto', 'fits', 'csv', 'ascii', 'votable'])
    ap.add_argument('--ref-ra-col', default='')
    ap.add_argument('--ref-dec-col', default='')
    ap.add_argument('--ref-x-col', default='')
    ap.add_argument('--ref-y-col', default='')
    ap.add_argument('--ref-id-col', default='')
    ap.add_argument('--copy-columns', default='')
    ap.add_argument('--cat-ra-col', default='ALPHA_J2000')
    ap.add_argument('--cat-dec-col', default='DELTA_J2000')
    ap.add_argument('--mode', default='sky', choices=['sky', 'pixel'])
    ap.add_argument('--radius-arcsec', type=float, default=1.0)
    ap.add_argument('--apply-shift', action='store_true')
    ap.add_argument('--allow-network', action='store_true')
    ap.add_argument('--tap-url', default='https://gea.esac.esa.int/tap-server/tap')
    ap.add_argument('--tap-adql', default='')
    ap.add_argument('--tap-radius-arcmin', type=float, default=0.0)
    ap.add_argument('--tap-timeout', type=float, default=120.0)
    ap.add_argument('--tap-limit', type=int, default=50000)
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    cols, rows = tsvio.read_catalog(a.catalog)
    n1 = len(rows)
    pixel = a.mode == 'pixel'
    x1 = np.array([tsvio.fnum(r.get('X_IMAGE')) for r in rows]); y1 = np.array([tsvio.fnum(r.get('Y_IMAGE')) for r in rows])
    ra1 = np.array([tsvio.fnum(r.get(a.cat_ra_col)) for r in rows]); dec1 = np.array([tsvio.fnum(r.get(a.cat_dec_col)) for r in rows])
    if not pixel and not np.isfinite(ra1).any():
        raise SystemExit('xmatch: the catalog has no sky coordinates (%s, %s); use pixel mode or add them' % (a.cat_ra_col, a.cat_dec_col))
    info = {}
    rraw = {}
    if a.source == 'tap':
        ok = np.isfinite(ra1) & np.isfinite(dec1)
        if not ok.any():
            raise SystemExit('xmatch: TAP needs sky coordinates in the catalog')
        # field centre and radius from the catalog (unit vectors, so RA wrap-around is safe)
        v = np.c_[np.cos(np.radians(dec1[ok])) * np.cos(np.radians(ra1[ok])), np.cos(np.radians(dec1[ok])) * np.sin(np.radians(ra1[ok])), np.sin(np.radians(dec1[ok]))].mean(0)
        rac = float(np.degrees(np.arctan2(v[1], v[0])) % 360.0); decc = float(np.degrees(np.arctan2(v[2], np.hypot(v[0], v[1]))))
        pr = xm.sky_match(np.array([rac]), np.array([decc]), ra1[ok], dec1[ok], 3600.0 * 90)
        rad = (float(pr['sep'].max()) / 3600.0 if len(pr['sep']) else 0.1) * 1.05 + a.radius_arcsec / 3600.0
        if a.tap_radius_arcmin > 0:
            rad = a.tap_radius_arcmin / 60.0
        rcols, rdata, rraw, tinfo = fetch_tap(a, rac, decc, rad, a.work)
        info.update(tinfo, centre=[rac, decc], radius_deg=rad)
        a.tap_limit = a.tap_limit
    else:
        if not a.ref_file or not os.path.isfile(a.ref_file):
            raise SystemExit('xmatch: reference table not found: %r' % a.ref_file)
        rcols, rdata = xm.read_table(a.ref_file, a.ref_format)
    if pixel:
        xc, yc = a.ref_x_col, a.ref_y_col
        if xc not in rdata or yc not in rdata:
            raise SystemExit('xmatch: pixel mode needs --ref-x-col and --ref-y-col present in the reference (columns: %s)' % ', '.join(rcols[:20]))
        x2 = xm._to_float_array(rdata[xc]); y2 = xm._to_float_array(rdata[yc])
        r = xm.crossmatch(None, None, None, None, a.radius_arcsec, pixel=True, x1=x1, y1=y1, x2=x2, y2=y2)
        n2 = len(x2)
    else:
        rc, dc = xm.find_columns(rcols, a.ref_ra_col, a.ref_dec_col)
        if rc is None or dc is None:
            raise SystemExit('xmatch: cannot find RA/Dec columns in the reference (columns: %s); set --ref-ra-col / --ref-dec-col' % ', '.join(rcols[:20]))
        ra2, dec2 = xm.parse_sky(rdata[rc], rdata[dc])
        r = xm.crossmatch(ra1, dec1, ra2, dec2, a.radius_arcsec, apply_shift=a.apply_shift)
        n2 = len(ra2)
    res = r['res']
    idc = a.ref_id_col if a.ref_id_col in rdata else None
    copy = [c.strip() for c in a.copy_columns.split(',') if c.strip()][:4]
    for c in copy:
        if c not in rdata:
            raise SystemExit('xmatch: copy column %r not in the reference (columns: %s)' % (c, ', '.join(rcols[:30])))
    out = []
    pairs_out = []
    for i, row in enumerate(rows):
        d = {}
        j = res['best'][i]
        d['XM_N'] = int(res['n'][i])
        if j >= 0:
            d.update(XM_SEP=res['sep'][i], XM_FLAG=int(res['flag'][i]), XM_DRA=res['dra'][i], XM_DDEC=res['ddec'][i])
            if idc:
                d['XM_ID'] = str(rraw[idc][j]) if idc in rraw else str(rdata[idc][j])   # exact text (Gaia ids exceed float precision)
            for k, c in enumerate(copy):
                v = rdata[c][j]
                d['XM_V%d' % (k + 1)] = float(v) if isinstance(v, (float, np.floating)) else str(v)
            pairs_out.append(dict(NUMBER=row['NUMBER'], REF_INDEX=int(j), SEP=res['sep'][i], FLAG=int(res['flag'][i])))
        else:
            d['XM_FLAG'] = 0
        out.append((row['NUMBER'], d))
    summary = dict(source=a.source, mode=a.mode, radius=a.radius_arcsec, n_catalog=n1, n_reference=int(n2), n_matched=r['n_matched'], n_unique=int(np.sum((res['flag'] & 3) == 1)),
                   n_ambiguous=int(np.sum((res['flag'] & 2) > 0)), n_mutual=int(np.sum((res['flag'] & 4) > 0)), shift=r['shift'], chance=r['chance'], tap=info or None,
                   median_sep=float(np.nanmedian(res['sep'])) if r['n_matched'] else None, copy_columns=copy, id_column=idc)
    with open(os.path.join(a.work, 'xmatch_summary.json'), 'w') as fh:
        json.dump(summary, fh, indent=1, default=lambda o: None)
    tsvio.write_table(os.path.join(a.work, 'xmatch_pairs.tsv'), ['NUMBER', 'REF_INDEX', 'SEP', 'FLAG'], pairs_out)
    if a.meta_out:
        ometa.update(a.meta_out, 'xmatch', dict(source=a.source, radius=a.radius_arcsec, n_matched=r['n_matched'], n_reference=int(n2), shift_dra=r['shift']['dra'], shift_ddec=r['shift']['ddec'],
                                                chance=(r['chance'] or {}).get('mean')), nrows=n1)
    tsvio.write_columns(COLUMNS, out)
    ch = r['chance']
    sys.stderr.write('xmatch: %d of %d matched within %.3g (%d unique, %d ambiguous, %d mutual); median offset %.3f, %.3f; expected chance matches %s\n' % (
        r['n_matched'], n1, a.radius_arcsec, summary['n_unique'], summary['n_ambiguous'], summary['n_mutual'], r['shift']['dra'], r['shift']['ddec'], ('%.1f' % ch['mean']) if ch else 'n/a'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
