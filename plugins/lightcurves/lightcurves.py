#!/usr/bin/env python3
"""Classify transient light curves and score them as supernovae.

    lightcurves.py --task classify --catalog CAT.tsv --work DIR (--source file --lc-file LC.tsv | --source moving --moving-dir DIR) [...]
    lightcurves.py --task plot --catalog CAT.tsv --work DIR --id ID --out PNG [same input options]

Input: a long-form table (id, mjd, flux | mag, err; optional number, ra, dec, host_offset_re) or the light curves of the Moving-objects pipeline
(lightcurves.json + transients.tsv in its work directory: forced photometry on the difference images, host offsets from the host association).
stdout (classify) = NUMBER + LC_* columns of the catalog objects that have a light curve (highest SN score when several).  Files in DIR:
lc_results.tsv (one row per light curve, including those without a catalog host), lc_results.json.
"""
import argparse
import csv
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

from ogfkit import tsvio, lcclass, lcsynth, meta as ometa  # noqa: E402

COLUMNS = ['LC_ID', 'LC_N', 'LC_CLASS', 'LC_PCLASS', 'LC_PSN', 'LC_PEAK_MJD', 'LC_PEAK_MAG', 'LC_RISE_D', 'LC_DECL', 'LC_CHI2R']
ALIASES = dict(id=('id', 'name', 'object', 'objid', 'obj', 'lcid', 'number'), t=('mjd', 't', 'time', 'jd', 'bjd'), flux=('flux', 'f', 'fluxcal', 'flux_ujy'),
               ferr=('flux_err', 'fluxerr', 'flux_error', 'ferr', 'e_flux', 'err', 'error', 'sigma'), mag=('mag', 'magnitude', 'mag_ab', 'm'),
               merr=('mag_err', 'magerr', 'e_mag', 'emag', 'mag_error', 'err', 'error', 'sigma'), number=('number', 'catalog_number', 'host_number', 'host_id', 'hostid'),
               ra=('ra', 'ra_deg', 'alpha_j2000'), dec=('dec', 'dec_deg', 'delta_j2000'), ho=('host_offset_re', 'offset_re', 'offset_in_re'), band=('band', 'filter', 'passband'))


def _pick(cols, key, exclude=()):
    low = {c.lower(): c for c in cols}
    for a in ALIASES[key]:
        if a in low and low[a] not in exclude:
            return low[a]
    return None


def _num(s):
    try:
        return float(s)
    except (TypeError, ValueError):
        return float('nan')


def read_long(path, flux_unit='auto', zp=25.0):
    """Long-form table -> list of dict(id, t, flux, err, number, ra, dec, ho) (flux on the zero point 25 scale)."""
    with open(path, newline='') as fh:
        txt = fh.read()
    lines = [l for l in txt.splitlines() if l.strip() and not l.startswith('#')]
    delim = '\t' if '\t' in lines[0] else ','
    rd = csv.DictReader(lines, delimiter=delim)
    cols = rd.fieldnames or []
    rows = list(rd)
    ci = _pick(cols, 'id'); ct = _pick(cols, 't')
    if ci is None or ct is None:
        raise SystemExit('lightcurves: need an id column (id/name/object/number) and a time column (mjd/t/time); found %s' % ', '.join(cols))
    cf = _pick(cols, 'flux'); cm = _pick(cols, 'mag')
    use_flux = (flux_unit == 'flux') or (flux_unit == 'auto' and cf is not None)
    if use_flux:
        if cf is None:
            raise SystemExit('lightcurves: no flux column')
        ce = _pick(cols, 'ferr', exclude=(cf,))
    else:
        if cm is None:
            raise SystemExit('lightcurves: no mag column')
        ce = _pick(cols, 'merr', exclude=(cm,))
    if ce is None:
        raise SystemExit('lightcurves: no error column (flux_err / mag_err)')
    cn = _pick(cols, 'number', exclude=(ci,)) if ci.lower() != 'number' else None
    cra, cde, cho = _pick(cols, 'ra'), _pick(cols, 'dec'), _pick(cols, 'ho')
    groups = {}
    for r in rows:
        groups.setdefault(r[ci].strip(), []).append(r)
    out = []
    for k, rr in groups.items():
        t = np.array([_num(r[ct]) for r in rr])
        if use_flux:
            sc = 10.0 ** (-0.4 * (zp - 25.0))
            f = np.array([_num(r[cf]) for r in rr]) * sc; e = np.array([_num(r[ce]) for r in rr]) * sc
        else:
            _, f, e = lcclass.from_mags(t, [_num(r[cm]) for r in rr], [_num(r[ce]) for r in rr])
        first = rr[0]
        out.append(dict(id=k, t=t, flux=f, err=e, number=(first.get(cn, '').strip() if cn else ''), ra=_num(first.get(cra)) if cra else float('nan'),
                        dec=_num(first.get(cde)) if cde else float('nan'), ho=_num(first.get(cho)) if cho else float('nan')))
    return out


def read_moving(wd):
    """lightcurves.json (+ transients.tsv) of the Moving-objects pipeline -> same dicts.  Flux (e/s, any zero point) is put on the ZP=25 scale with the zero point
    implied by the detected points (mag + 2.5 log10 flux) of each epoch."""
    fn = os.path.join(wd, 'lightcurves.json')
    if not os.path.isfile(fn):
        raise SystemExit('lightcurves: %s not found (run Moving objects > Light Curve first)' % fn)
    L = json.load(open(fn))
    tr = {}
    tf = os.path.join(wd, 'transients.tsv')
    if os.path.isfile(tf):
        c, rows = tsvio.read_catalog(tf)
        tr = {r['id']: r for r in rows}
    out = []
    for lc in L:
        t = np.array(lc['t'], float); f = np.array(lc['flux'], float); e = np.array(lc['err'], float)
        m = np.array([np.nan if v is None else v for v in lc['mag']], float)
        with np.errstate(invalid='ignore', divide='ignore'):
            zpi = m + 2.5 * np.log10(f)
        zpm = np.nanmedian(zpi) if np.isfinite(zpi).any() else 25.0
        zpi = np.where(np.isfinite(zpi), zpi, zpm)
        sc = 10.0 ** (-0.4 * (zpi - 25.0))
        r = tr.get(str(lc['id']), {})
        out.append(dict(id='T%s' % lc['id'], t=t, flux=f * sc, err=e * sc, number=str(r.get('host_id', '') or ''), ra=_num(r.get('ra')), dec=_num(r.get('dec')),
                        ho=_num(r.get('offset_re'))))
    return out


def load_models(model_dir):
    m = lcclass.load_model(os.path.join(model_dir, 'model_lc.json'))
    hp = os.path.join(model_dir, 'model_host.json')
    return m, (lcclass.load_model(hp) if os.path.isfile(hp) else None)


def analyse(lcs, model, host_model, use_host=True, min_prob=0.5, min_points=5, min_snr=3.0):
    res = []
    for lc in lcs:
        ho = lc['ho'] if use_host else None
        r = lcclass.classify(lc['t'], lc['flux'], lc['err'], model, host_model if use_host else None, ho if (ho is not None and np.isfinite(ho)) else None, min_prob, min_points, min_snr)
        a = r.get('aux') or {}
        F = r.get('features') or {}
        t = np.asarray(lc['t'], float); f = np.asarray(lc['flux'], float)
        rise = float('nan')
        if a:
            o = np.argsort(t)
            if f[o][0] < 0.5 * a['f_peak']:
                rise = a['t_peak'] - float(t[o][0])
        row = dict(ID=lc['id'], N=r['n'], CLASS=r['cls'], BEST_GUESS=r.get('best_guess', ''), PCLASS=r['p_class'], PSN=r['p_sn'], USED_HOST=int(r['used_host']),
                   PEAK_MJD=a.get('t_peak', float('nan')), PEAK_MAG=a.get('m_peak', float('nan')), RISE_D=rise, DECL=F.get('decl_rate', float('nan')),
                   CHI2R=a.get('chi2r_const', float('nan')), HOST_OFFSET_RE=lc['ho'], NUMBER=lc['number'], NOTE=r['note'])
        for c in lcsynth.CLASSES:
            row['P_' + c] = r['probs'].get(c, float('nan'))
        res.append(row)
    return res


def link(res, lcs, rows, mode, radius):
    """-> dict catalog NUMBER -> index of the best (highest PSN) light curve."""
    nums = {r['NUMBER']: i for i, r in enumerate(rows)}
    ra = np.array([tsvio.fnum(r.get('ALPHA_J2000')) for r in rows]); de = np.array([tsvio.fnum(r.get('DELTA_J2000')) for r in rows])
    best = {}
    for k, (r, lc) in enumerate(zip(res, lcs)):
        tgt = None
        if mode in ('auto', 'number') and lc['number'] in nums:
            tgt = lc['number']
        if tgt is None and mode in ('auto', 'number'):
            d = ''.join(ch for ch in str(lc['id']) if ch.isdigit())
            if mode == 'number' and lc['number'] == '' and d in nums:
                tgt = d
            elif mode == 'auto' and lc['number'] == '' and str(lc['id']) in nums:
                tgt = str(lc['id'])
        if tgt is None and mode in ('auto', 'position') and np.isfinite(lc['ra']) and np.isfinite(lc['dec']) and np.isfinite(ra).any():
            dd = np.hypot((ra - lc['ra']) * math.cos(math.radians(lc['dec'])), de - lc['dec']) * 3600.0
            dd = np.where(np.isfinite(dd), dd, 1e9)
            j = int(np.argmin(dd))
            if dd[j] <= radius:
                tgt = rows[j]['NUMBER']
        if tgt is None:
            continue
        r['NUMBER'] = tgt
        ps = r['PSN'] if np.isfinite(r['PSN']) else -1
        if tgt not in best or ps > (res[best[tgt]]['PSN'] if np.isfinite(res[best[tgt]]['PSN']) else -1):
            best[tgt] = k
    return best


def plot(lc, row, model, host_model, out, mdl_use_host=True):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.4, 3.4), dpi=100)
    t = np.asarray(lc['t']); f = np.asarray(lc['flux']); e = np.asarray(lc['err'])
    ax.errorbar(t, f, e, fmt='o', ms=3, color='k', lw=0.8)
    if len(t) >= 5:
        b = lcclass.fit_bazin(t, f, e)
        if b.get('ok'):
            tt = np.linspace(t.min(), t.max(), 400)
            ax.plot(tt, lcclass.bazin(tt, b['A'], b['t0'], b['tau_r'], b['tau_f'], b['B']), color='tab:red', lw=1, label='Bazin fit')
    ax.set_xlabel('MJD'); ax.set_ylabel('flux (ZP=25)')
    pc = row['PCLASS']
    ax.set_title('%s: %s  p=%s  P(SN)=%s  N=%d' % (row['ID'], row['CLASS'] if row['CLASS'] != 'unclassified' else 'unclassified (%s)' % (row['BEST_GUESS'] or '-'),
                 '%.2f' % pc if np.isfinite(pc) else '-', '%.2f' % row['PSN'] if np.isfinite(row['PSN']) else '-', row['N']), fontsize=9)
    ax.axhline(0, color='0.7', lw=0.5)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--task', default='classify', choices=['classify', 'plot'])
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--work', default='.')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--source', default='file', choices=['file', 'moving'])
    ap.add_argument('--lc-file', default='')
    ap.add_argument('--moving-dir', default='')
    ap.add_argument('--flux-unit', default='auto', choices=['auto', 'flux', 'mag'])
    ap.add_argument('--zp', type=float, default=25.0)
    ap.add_argument('--link', default='auto', choices=['auto', 'number', 'position'])
    ap.add_argument('--link-radius-arcsec', type=float, default=3.0)
    ap.add_argument('--min-prob', type=float, default=0.5)
    ap.add_argument('--min-points', type=int, default=5)
    ap.add_argument('--min-snr', type=float, default=3.0)
    ap.add_argument('--no-host', action='store_true')
    ap.add_argument('--model-dir', default='')
    ap.add_argument('--id', default='')
    ap.add_argument('--out', default='')
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    if a.source == 'moving':
        if not a.moving_dir:
            raise SystemExit('lightcurves: set the Moving-objects work directory')
        lcs = read_moving(a.moving_dir)
    else:
        if not a.lc_file or not os.path.isfile(a.lc_file):
            raise SystemExit('lightcurves: light-curve file not found: %r' % a.lc_file)
        if a.lc_file.endswith('.json'):
            import shutil
            lcs = read_moving(os.path.dirname(os.path.abspath(a.lc_file)))
        else:
            lcs = read_long(a.lc_file, a.flux_unit, a.zp)
    model, hmodel = load_models(a.model_dir or HERE)
    res = analyse(lcs, model, hmodel, not a.no_host, a.min_prob, a.min_points, a.min_snr)
    if a.task == 'plot':
        k = [i for i, r in enumerate(res) if r['ID'] == a.id]
        if not k:
            raise SystemExit('lightcurves: no light curve with id %r' % a.id)
        plot(lcs[k[0]], res[k[0]], model, hmodel, a.out)
        return 0
    cols, rows = tsvio.read_catalog(a.catalog)
    best = link(res, lcs, rows, a.link, a.link_radius_arcsec)
    out = []
    for r in rows:
        k = best.get(r['NUMBER'])
        if k is None:
            out.append((r['NUMBER'], {}))
            continue
        q = res[k]
        out.append((r['NUMBER'], dict(LC_ID=q['ID'], LC_N=q['N'], LC_CLASS=q['CLASS'], LC_PCLASS=q['PCLASS'], LC_PSN=q['PSN'], LC_PEAK_MJD=q['PEAK_MJD'],
                                     LC_PEAK_MAG=q['PEAK_MAG'], LC_RISE_D=q['RISE_D'], LC_DECL=q['DECL'], LC_CHI2R=q['CHI2R'])))
    rcols = ['ID', 'N', 'CLASS', 'BEST_GUESS', 'PCLASS', 'PSN', 'USED_HOST', 'PEAK_MJD', 'PEAK_MAG', 'RISE_D', 'DECL', 'CHI2R', 'HOST_OFFSET_RE', 'NUMBER', 'NOTE'] + ['P_' + c for c in lcsynth.CLASSES]
    tsvio.write_table(os.path.join(a.work, 'lc_results.tsv'), rcols, res)
    counts = {}
    for r in res:
        counts[r['CLASS']] = counts.get(r['CLASS'], 0) + 1
    summ = dict(n_lightcurves=len(res), n_linked=len(best), class_counts=counts, n_sn_like=int(sum(1 for r in res if np.isfinite(r['PSN']) and r['PSN'] >= 0.5)), source=a.source,
                model=dict(classes=model['classes'], meta=model.get('meta')), min_prob=a.min_prob)
    with open(os.path.join(a.work, 'lc_results.json'), 'w') as fh:
        json.dump(dict(summary=summ, results=[{k: (None if isinstance(v, float) and not math.isfinite(v) else v) for k, v in r.items()} for r in res]), fh, indent=1)
    if a.meta_out:
        ometa.update(a.meta_out, 'lightcurves', dict(n_lightcurves=len(res), n_linked=len(best), n_sn_like=summ['n_sn_like']), nrows=len(rows))
    tsvio.write_columns(COLUMNS, out)
    sys.stderr.write('lightcurves: %d light curves (%s), %d linked to catalog objects, %d with P(SN) >= 0.5\n' % (len(res), ', '.join('%s %d' % kv for kv in sorted(counts.items())), len(best), summ['n_sn_like']))
    return 0


if __name__ == '__main__':
    sys.exit(main())
