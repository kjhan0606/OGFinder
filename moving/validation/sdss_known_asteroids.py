"""Real-field linking precision / recall: SDSS run 94 (Stripe 82, 1998-09-19, near opposition on the ecliptic), camcol N, fields 136-167.

Each field is imaged in r, i, u, z, g in that order, ~72 s apart (five exposures of the same sky, every asteroid moves ~2.5 arcsec across them), so every
mover has up to 5 detections.  Steps (everything cached under $OGF_DATA_CACHE or ~/.cache/ogfinder_regression/sdss_run94):
  1. download the frames (data.sdss.org), detect sources with sep (4 sigma, >= 4 px) per band, time stamp = frame TAI + 0.5 exposure + row * 26.3 ms (drift scan)
  2. truth: IMCCE SkyBoT cone searches (known asteroids of the MPC database with predicted positions and rates, observer 645) around every field;
     a known object is "detectable" (the ceiling) when detections in >= 3 bands follow its predicted *motion* (offset from the prediction allowed up to 3 arcsec)
     and move by >= 1.2 arcsec
  3. baseline: moving.tracklet.link_exposures (tile by tile, fixed observer), a tracklet counts as true when >= 3 members are detections of the same known object
  4. orbit stage: moving.orbitlink population prior (JPL SBDB orbits, random mean anomaly, real Earth state, V < 22) -> tracklets outside the 99 % highest-density
     region of the apparent rate vector are dropped.  The prior sees no data of the field.
Precision is a LOWER bound: a tracklet that is not a known object may be an unknown real asteroid (fainter than the catalogue completeness).
Usage: sdss_known_asteroids.py OUT.json [--camcol 3] [--fraction 0.99]
"""
import argparse, glob, json, os, pickle, subprocess, sys, time
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from scipy.spatial import cKDTree
from moving import tracklet as T, orbitlink as OL, skybot

CACHE = os.path.expanduser(os.environ.get('OGF_DATA_CACHE', '~/.cache/ogfinder_regression')) + '/sdss_run94'
BANDS = ['r', 'i', 'u', 'z', 'g']
ROWT = 53.9013 / 2048.0
TAI_UTC = 31.0


def download(camcol, fields=range(136, 168)):
    d = '%s/c%d' % (CACHE, camcol)
    os.makedirs(d, exist_ok=True)
    todo = []
    for f in fields:
        for b in BANDS:
            n = 'frame-%s-000094-%d-%04d' % (b, camcol, f)
            if not os.path.isfile('%s/%s.fits' % (d, n)):
                todo.append(n)
    if todo:
        cmd = "cat | xargs -P 4 -I{} sh -c 'curl -sf -m 300 -o %s/{}.fits.bz2 https://data.sdss.org/sas/dr17/eboss/photoObj/frames/301/94/%d/{}.fits.bz2 && bunzip2 -f %s/{}.fits.bz2'" % (d, camcol, d)
        subprocess.run(cmd, shell=True, input='\n'.join(todo), text=True, check=False)
    return d


def detect(path):
    import sep
    from astropy.io import fits
    from astropy.wcs import WCS
    import warnings
    warnings.simplefilter('ignore')
    b = os.path.basename(path).split('-')[1]
    h = fits.getheader(path)
    d = np.ascontiguousarray(fits.getdata(path).astype(np.float64))
    w = WCS(h)
    bkg = sep.Background(d, bw=64, bh=64)
    sub = d - bkg.back()
    o = sep.extract(sub, 4.0, err=bkg.rms(), minarea=4, deblend_nthresh=32, deblend_cont=0.005)
    fl, fe, _ = sep.sum_circle(sub, o['x'], o['y'], 3.0, err=bkg.rms())
    ra, dec = w.all_pix2world(o['x'], o['y'], 0)
    mjd = (h['TAI'] + 0.5 * 53.9013 + o['y'] * ROWT - TAI_UTC) / 86400.0
    return [dict(ex=BANDS.index(b), ra=float(ra[i]), dec=float(dec[i]), t=float(mjd[i]), flux=float(fl[i]), snr=float(fl[i] / max(fe[i], 1e-9)), sig=0.25)
            for i in range(len(o))], dict(ra0=float(ra.mean()), dec0=float(dec.mean()), t0=float(mjd.mean()), band=b)


def load_dets(d, snr_min=5.0):
    from concurrent.futures import ProcessPoolExecutor
    pk = d + '/dets.pkl'
    if os.path.isfile(pk):
        raw, meta = pickle.load(open(pk, 'rb'))
    else:
        raw, meta = [], []
        with ProcessPoolExecutor(4) as ex:
            for o, m in ex.map(detect, sorted(glob.glob(d + '/frame-?-000094-*.fits'))):
                raw += o; meta.append(m)
        pickle.dump((raw, meta), open(pk, 'wb'))
    dets = [x for x in raw if x['snr'] >= snr_min]
    out = []
    for ex in range(5):                                           # the 128-row frame overlap detects sources twice
        dd = sorted([x for x in dets if x['ex'] == ex], key=lambda x: -x['snr'])
        xy = np.array([[x['ra'] * np.cos(np.radians(x['dec'])), x['dec']] for x in dd])
        drop = set(max(i, j) for i, j in cKDTree(xy).query_pairs(0.5 / 3600.0))
        out += [x for k, x in enumerate(dd) if k not in drop]
    for k, x in enumerate(out):
        x['id'] = k
    return out, meta


def known_objects(d, meta):
    pk = d + '/known.pkl'
    if os.path.isfile(pk):
        return pickle.load(open(pk, 'rb'))
    known = {}
    for m in [m for m in meta if m['band'] == 'r']:
        for _ in range(3):
            try:
                res = skybot.cone(m['ra0'], m['dec0'], 0.3, m['t0'], observer='645')
                break
            except Exception:
                res = []
        for o in res:
            known.setdefault((o['num'], o['name']), dict(o, t0=m['t0']))
    pickle.dump(known, open(pk, 'wb'))
    return known


def match_truth(dets, known, tol=0.6, search=3.0, min_move=1.2):
    by_id = {x['id']: x for x in dets}
    by_ex = {ex: [x for x in dets if x['ex'] == ex] for ex in range(5)}
    cosd = np.cos(np.radians(np.mean([x['dec'] for x in dets])))
    trees = {ex: cKDTree(np.array([[x['ra'] * cosd, x['dec']] for x in by_ex[ex]])) for ex in by_ex}

    def pred(o, x):
        dt = (x['t'] - o['t0']) * 24.0
        return o['ra'] + o['dra'] * dt / 3600.0 / np.cos(np.radians(o['dec'])), o['dec'] + o['ddec'] * dt / 3600.0
    obj = {}
    for key, o in known.items():
        best = {}
        for sex in range(5):
            for c in trees[sex].query_ball_point([o['ra'] * cosd, o['dec']], (search + 1.0) / 3600.0 + 0.003):
                sd = by_ex[sex][c]
                pr, pd = pred(o, sd)
                off = ((sd['ra'] - pr) * cosd, sd['dec'] - pd)
                if np.hypot(*off) * 3600 > search:
                    continue
                m = {sex: sd['id']}
                for ex in range(5):
                    if ex == sex:
                        continue
                    bb = None
                    for c2 in trees[ex].query_ball_point([o['ra'] * cosd + off[0], o['dec'] + off[1]], 0.03):
                        x = by_ex[ex][c2]
                        p2r, p2d = pred(o, x)
                        r = np.hypot((x['ra'] - p2r) * cosd - off[0], x['dec'] - p2d - off[1]) * 3600
                        if r < tol and (bb is None or r < bb[0]):
                            bb = (r, x['id'])
                    if bb:
                        m[ex] = bb[1]
                if len(m) < 3:
                    continue
                ids = [by_id[v] for v in m.values()]
                if max(np.hypot((a['ra'] - b['ra']) * cosd, a['dec'] - b['dec']) for a in ids for b in ids) * 3600 < min_move:
                    continue
                if len(m) > len(best):
                    best = m
        obj[key] = best
    return obj


KW = dict(tol_arcsec=0.8, min_exposures=3, rate_min_ash=10.0, rate_max_ash=80.0, min_disp_arcsec=1.0, rescore=False, bound_orbit=False, cluster_dedupe=False,
          max_tracklets=100000, k_max=0.0, rate_per_k_ash=1e6, rate_slack_ash=1e6, max_cand_per_exposure=20000)


def link_tiles(dets, tile=0.5, pad=0.06, **kw):
    allt = []
    ras = np.array([x['ra'] for x in dets])
    for a in np.arange(ras.min(), ras.max() + tile, tile)[:-1]:
        sub = [x for x in dets if a - pad <= x['ra'] < a + tile + pad]
        if len(sub) >= 30:
            allt += T.link_exposures(sub, **kw)
    seen, out = [], []
    for t in sorted(allt, key=lambda t: t.get('score', 0.0)):
        ids = frozenset(t['members'])
        if any(len(ids & s) >= 2 for s in seen):
            continue
        seen.append(ids); out.append(t)
    return out


def classify(trs, obj):
    d2o = {}
    for k, m in obj.items():
        for i in m.values():
            d2o.setdefault(i, set()).add(k)
    res = []
    for t in trs:
        cnt = {}
        for i in t['members']:
            for k in d2o.get(i, ()):
                cnt[k] = cnt.get(k, 0) + 1
        k, n = max(cnt.items(), key=lambda kv: kv[1]) if cnt else (None, 0)
        res.append((n >= 3, k))
    return res


def summarize(trs, cls, ceiling):
    ok = np.array([c[0] for c in cls], bool)
    rec = set(c[1] for c in cls if c[0])
    return dict(n_tracklets=len(trs), n_true=int(ok.sum()), n_other=int((~ok).sum()), precision_lower_bound=float(ok.mean()) if len(ok) else None,
                recovered=len(rec), ceiling=ceiling, recall=len(rec) / max(ceiling, 1))


def run(camcol, fraction, elements_path):
    d = download(camcol)
    dets, meta = load_dets(d)
    known = known_objects(d, meta)
    obj = match_truth(dets, known)
    ceil = [k for k, m in obj.items() if len(m) >= 3]
    t0 = time.time()
    trs = link_tiles(dets, **KW)
    cls = classify(trs, obj)
    out = dict(camcol=camcol, n_detections=len(dets), n_known=len(known), ceiling=len(ceil), ceiling_ge4=sum(len(m) >= 4 for m in obj.values()), ceiling_5=sum(len(m) == 5 for m in obj.values()),
               link_seconds=time.time() - t0, baseline=summarize(trs, cls, len(ceil)))
    el = OL.fetch_elements(elements_path)
    prior = OL.prior_for_tracklets(trs, el)
    for fr in (0.95, fraction, 0.999):
        kept = OL.vet_tracklets([dict(t) for t in trs], prior, fr)
        m = np.array([t['in_hdr'] for t in kept])
        sub = [t for t, k in zip(trs, m) if k]
        out['orbit_stage_%g' % fr] = summarize(sub, [c for c, k in zip(cls, m) if k], len(ceil))
    kept = OL.vet_tracklets([dict(t) for t in trs], prior, fraction)
    byn = {}
    for t, c in zip(kept, cls):
        if t['in_hdr']:
            n = len(t['members']); a = byn.setdefault(str(n), [0, 0]); a[0] += int(c[0]); a[1] += 1
    out['kept_by_members_true_total'] = byn
    out['prior_population_size'] = int(prior.n_pop)
    return out


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('out')
    ap.add_argument('--camcol', type=int, default=3)
    ap.add_argument('--fraction', type=float, default=0.99)
    ap.add_argument('--elements', default=CACHE + '/../sbdb_elements_H19.5.npz')
    a = ap.parse_args()
    r = run(a.camcol, a.fraction, a.elements)
    json.dump(r, open(a.out, 'w'), indent=1)
    print(json.dumps(r)[:1500])
