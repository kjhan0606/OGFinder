#!/usr/bin/env python3
"""Download REAL ZTF light curves with labels for validating / retraining the light-curve classifier (network needed).

    python3 plugins/lightcurves/fetch_ztf.py OUTDIR [--n-agn 150] [--n-var 60] [--window 220] [--jobs 6]

Sources (public, no account):
  * SN types: VizieR J/ApJ/895/32 (Fremling et al. 2020, ZTF Bright Transient Survey, 764 SNe with SPECTROSCOPIC types) via the TAP service
    tapvizier.cds.unistra.fr.   Ia* -> SNIa, II/IIb/IIn/II-87A -> SNII, Ib/Ic/Ib/c/Ic-BL/Ibn/Ic-pec -> SNIbc; SLSN and 'ambiguous' are skipped.
  * AGN / periodic variables: the ALeRCE lc_classifier labels (api.alerce.online) at probability >= 0.9.  THESE LABELS ARE A MACHINE
    CLASSIFIER'S, NOT SPECTROSCOPIC: AGN -> AGN; RRL, E, CEP, DSCT, LPV -> variable.
  * Light curves: ALeRCE `objects/<oid>/lightcurve` (ZTF alert detections: magpsf, sigmapsf, isdiffpos; difference-image PSF photometry).
Output: OUTDIR/ztf_lcs.json = list of dict(oid, label, source, band, t, mag, err, sign).  The band with more detections is used; every object is cut to a
`--window`-day window (SNe: first detection - 20 d .. + window; others: the window with most detections) so that the observed duration is not a
class feature.  Convert to flux with ogfkit.lcclass.from_mags (sign gives the flux sign of the difference image)."""
import argparse, csv, io, json, os, sys, time, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor

TAP = 'https://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
API = 'https://api.alerce.online/ztf/v1'
MAP = {'Ia': 'SNIa', 'Ia-02cx': 'SNIa', 'Ia-91T': 'SNIa', 'Ia-91bg': 'SNIa', 'Ia-csm': 'SNIa', 'Ia-SC': 'SNIa',
       'II': 'SNII', 'IIb': 'SNII', 'IIn': 'SNII', 'II-87A': 'SNII',
       'Ib': 'SNIbc', 'Ic': 'SNIbc', 'Ib/c': 'SNIbc', 'Ic-BL': 'SNIbc', 'Ibn': 'SNIbc', 'Ic-pec': 'SNIbc'}
ALERCE = {'AGN': 'AGN', 'RRL': 'variable', 'E': 'variable', 'CEP': 'variable', 'DSCT': 'variable', 'LPV': 'variable'}


def get(url, data=None, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers={'User-Agent': 'OGFinder-fetch_ztf'})
            return urllib.request.urlopen(req, timeout=60).read().decode()
        except Exception as e:
            err = e; time.sleep(1.5 * (k + 1))
    raise err


def bts_list():
    q = urllib.parse.urlencode(dict(REQUEST='doQuery', LANG='ADQL', FORMAT='csv', QUERY='SELECT ZTF, SNtype, z FROM "J/ApJ/895/32/table1"')).encode()
    rows = list(csv.DictReader(io.StringIO(get(TAP, q))))
    return [(r['ZTF'], MAP[r['SNtype']], 'BTS spectroscopic %s' % r['SNtype']) for r in rows if r['SNtype'] in MAP]


def alerce_list(cls, n):
    out = []
    page = 1
    while len(out) < n and page < 20:
        u = '%s/objects/?classifier=lc_classifier&class=%s&ranking=1&page_size=100&page=%d&probability=0.9&format=json' % (API, cls, page)
        its = json.loads(get(u)).get('items', [])
        if not its: break
        out += [(i['oid'], ALERCE[cls], 'ALeRCE lc_classifier %s p=%.2f' % (cls, i.get('probability') or 0)) for i in its
                if (i.get('probability') or 0) >= 0.9 and int(i.get('ndet') or 0) >= 15]
        page += 1
    return out[:n]


def lightcurve(item, window):
    oid, label, src = item
    try:
        det = json.loads(get('%s/objects/%s/lightcurve' % (API, oid))).get('detections', [])
    except Exception:
        return None
    best = None
    for fid in (1, 2):
        d = sorted((x for x in det if x.get('fid') == fid and x.get('magpsf') is not None and x.get('sigmapsf')), key=lambda x: x['mjd'])
        if best is None or len(d) > len(best): best = d
    if not best or len(best) < 8: return None
    t = [x['mjd'] for x in best]
    if label.startswith('SN'):
        lo = t[0] - 20.0
    else:   # densest window
        lo = max(t, key=lambda a: sum(1 for b in t if a <= b < a + window))
    sel = [x for x in best if lo <= x['mjd'] < lo + window]
    if len(sel) < 8: return None
    return dict(oid=oid, label=label, source=src, band='g' if best[0]['fid'] == 1 else 'r', t=[x['mjd'] for x in sel], mag=[x['magpsf'] for x in sel],
                err=[x['sigmapsf'] for x in sel], sign=[int(x.get('isdiffpos') or 1) for x in sel])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('outdir'); ap.add_argument('--n-agn', type=int, default=150); ap.add_argument('--n-var', type=int, default=60)
    ap.add_argument('--window', type=float, default=220.0); ap.add_argument('--jobs', type=int, default=6)
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    items = bts_list()
    items += alerce_list('AGN', a.n_agn)
    for c in ('RRL', 'E', 'CEP', 'DSCT', 'LPV'):
        items += alerce_list(c, a.n_var)
    print('objects to fetch:', len(items), flush=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        res = list(ex.map(lambda it: lightcurve(it, a.window), items))
    res = [r for r in res if r]
    json.dump(res, open(os.path.join(a.outdir, 'ztf_lcs.json'), 'w'))
    from collections import Counter
    print('saved', len(res), dict(Counter(r['label'] for r in res)))


if __name__ == '__main__':
    main()
