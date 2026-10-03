#!/usr/bin/env python3
"""Multi-component galaxy fitting (GALFIT-like) for the catalog objects of the current image.

    multifit.py IMAGE --catalog TSV --work DIR [--mask FITS] [--psf-model JSON | --psf FITS | --psf-fwhm PX]
        [--model sersic|exp|dev|psf|bulge+disk|psf+sersic|auto|decomp|bar|ring|spiral|bulge+disk+bar|bulge+disk+ring|structure] [--neighbours fit|mask|ignore] [--objects "1;5;9"] [--max-objects 100] ...
    multifit.py IMAGE --config FILE.json --work DIR          (explicit components of one object, GALFIT-feedme-like)

With --catalog the add_columns contract is followed: stdout = `NUMBER` + GF_* columns, logs on stderr.
Every cutout is fitted with PSF-convolved components + a common sky; neighbours are fitted simultaneously (`fit`), masked (`mask`) or ignored.
The PSF is the model given by --psf-model (spatially varying, evaluated at the object), else the --psf image, else --psf-fwhm (Gaussian), else a model built from the field stars.
Files in DIR/multifit_*: results.tsv (one row per fitted component), model.fits / residual.fits (full frame, sky not included), montage.png, psf.json.
"""
import argparse
import json
import math
import multiprocessing as mp
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio, psfmodel as pm, multifit as mf, meta as ometa, autodecomp as ad, structfit as sf  # noqa: E402

COLUMNS = ['GF_X', 'GF_Y', 'GF_MAG', 'GF_MAGERR', 'GF_RE', 'GF_REERR', 'GF_N', 'GF_NERR', 'GF_Q', 'GF_PA', 'GF_BT', 'GF_MAG2', 'GF_SKY', 'GF_CHI2', 'GF_NCOMP', 'GF_NNEIGH', 'GF_RESFRAC', 'GF_FLAG']
ST_COLUMNS = ['ST_TYPE', 'ST_MAG', 'ST_BT', 'ST_BARFRAC', 'ST_BARLEN', 'ST_BARPA', 'ST_BARQ', 'ST_BARC0', 'ST_RINGR', 'ST_RINGW', 'ST_RINGFRAC', 'ST_SPTHETA', 'ST_SPAMP', 'ST_DBIC', 'ST_CHI2', 'ST_FLAG']
AD_COLUMNS = ['AD_TYPE', 'AD_MAG', 'AD_BT', 'AD_MAGB', 'AD_MAGD', 'AD_MAGN', 'AD_REB', 'AD_NB', 'AD_QB', 'AD_PAB', 'AD_RED', 'AD_QD', 'AD_PAD', 'AD_RE1', 'AD_N1', 'AD_DBIC', 'AD_CHI2', 'AD_BAR', 'AD_FLAG']
G = {}


def num(r, k, default=float('nan')):
    return tsvio.fnum(r.get(k), default)


def start_values(r, zp):
    x, y = num(r, 'X_IMAGE') - 1, num(r, 'Y_IMAGE') - 1
    fl = num(r, 'FLUX_AUTO')
    if not (fl > 0):
        m = num(r, 'MAG_AUTO')
        fl = 10 ** (-0.4 * (m - zp)) if np.isfinite(m) and m < 90 else float('nan')
    A, B = num(r, 'A_IMAGE', 2.0), num(r, 'B_IMAGE', 2.0)
    if not (A > 0):
        A = 2.0
    if not (B > 0):
        B = A
    th = num(r, 'THETA_IMAGE', 0.0)
    fr = num(r, 'FLUX_RADIUS')
    re0 = max(fr * 1.1 if fr > 0 else 1.5 * A, 1.0)
    return dict(x=x, y=y, flux=fl, A=A, B=min(B, A), pa=th if np.isfinite(th) else 0.0, re=re0, q=min(max(B / A, 0.15), 1.0), cs=num(r, 'CLASS_STAR', 0.0))


def ellipse_mask(shape, x, y, a, b, pa_deg):
    yy, xx = np.mgrid[:shape[0], :shape[1]]
    t = math.radians(pa_deg)
    dx, dy = xx - x, yy - y
    u = dx * math.cos(t) + dy * math.sin(t)
    v = -dx * math.sin(t) + dy * math.cos(t)
    return (u / max(a, 0.5)) ** 2 + (v / max(b, 0.5)) ** 2 <= 1.0


def get_psf(x, y):
    return G['psf']


def fit_object(i):
    """Worker: fit catalog row i.  Returns (row dict, per-component results, cutout geometry, model/residual stamps)."""
    a = G['args']
    rows, sv, data, mask, rms = G['rows'], G['sv'], G['data'], G['mask'], G['rms']
    s = sv[i]
    ny, nx = data.shape
    zp = a.mag_zeropoint
    r = rows[i]
    out = dict(NUMBER=r['NUMBER'], i=i)
    if not (np.isfinite(s['x']) and np.isfinite(s['y']) and s['flux'] > 0):
        out['GF_FLAG'] = 1024
        return out
    half = int(min(max(a.size_scale * s['A'], a.min_size / 2.0), a.max_size / 2.0))
    x0, x1 = int(round(s['x'])) - half, int(round(s['x'])) + half + 1
    y0, y1 = int(round(s['y'])) - half, int(round(s['y'])) + half + 1
    cx0, cx1, cy0, cy1 = max(0, x0), min(nx, x1), max(0, y0), min(ny, y1)
    cut = data[cy0:cy1, cx0:cx1].astype(float)
    cm = (mask[cy0:cy1, cx0:cx1].copy() if mask is not None else np.zeros(cut.shape, bool)) | ~np.isfinite(cut)
    rcut = rms[cy0:cy1, cx0:cx1] if np.ndim(rms) == 2 else rms
    lx, ly = s['x'] - cx0, s['y'] - cy0
    flag = 0
    if cx0 != x0 or cx1 != x1 or cy0 != y0 or cy1 != y1:
        flag |= mf.FLAGS['EDGE']
    # neighbours
    extra, nn = [], 0
    if a.neighbours != 'ignore':
        tree = G['tree']
        cand = tree.query_ball_point([s['x'], s['y']], math.hypot(half, half) + 3 * 30)
        near = []
        for j in cand:
            if j == i:
                continue
            t = sv[j]
            if not (np.isfinite(t['x']) and t['flux'] > 0):
                continue
            inside = cx0 - 0.5 <= t['x'] <= cx1 - 0.5 and cy0 - 0.5 <= t['y'] <= cy1 - 0.5
            near.append((t['flux'], j, inside))
        near.sort(reverse=True)
        nfit = 0
        for fl, j, inside in near:
            t = sv[j]
            if a.neighbours == 'fit' and inside and nfit < a.max_neighbours and t['flux'] > 0.01 * s['flux']:
                nfit += 1
                stellar = t['cs'] > 0.8 and t['A'] < 1.6 * G['psf_sigma']
                bx = (t['x'] - cx0 - 3.0, t['x'] - cx0 + 3.0); by = (t['y'] - cy0 - 3.0, t['y'] - cy0 + 3.0)
                if stellar:
                    extra.append(dict(kind='psf', x=t['x'] - cx0, y=t['y'] - cy0, flux=t['flux'], bounds=dict(x=bx, y=by)))
                else:
                    extra.append(dict(kind='sersic', x=t['x'] - cx0, y=t['y'] - cy0, flux=t['flux'], re=t['re'], q=t['q'], pa=t['pa'], n=1.5, bounds=dict(x=bx, y=by)))
            else:
                em = ellipse_mask(cut.shape, t['x'] - cx0, t['y'] - cy0, a.neighbour_radius * t['A'], a.neighbour_radius * t['B'], t['pa'])
                cm |= em
                nn += 1
        nn += nfit
        if nfit:
            flag |= mf.FLAGS['NEIGHBOUR']
    # never mask the target itself
    tgt = ellipse_mask(cut.shape, lx, ly, 1.5, 1.5, 0.0)
    cm &= ~tgt
    psf = G['psf']
    kw = dict(psf=psf, rms=rcut, mask=cm, sky=a.sky, zp=zp, gain=a.gain or None, max_nfev=a.max_nfev)
    name = a.model
    try:
        if name == 'auto':
            name, res, allr = mf.auto_select(cut, lx, ly, s['flux'], s['re'], s['q'], s['pa'], models=('psf', 'sersic', 'bulge+disk'), bic_margin=a.bic_margin, **kw) if not extra else (None, None, None)
            if res is None:
                name = 'sersic'
                res, _ = mf.fit_preset(cut, name, lx, ly, s['flux'], s['re'], s['q'], s['pa'], extra=extra, **kw)
        elif name == 'decomp':
            fw = max(2.355 * G['psf_sigma'], 1.5)
            models = ('sersic', 'bulge+disk') + (() if a.nucleus == 'off' else ('nucleus+bulge+disk',))
            d = ad.decompose(cut, lx, ly, psf=psf, rms=rcut, mask=cm, sky=a.sky, zp=zp, flux0=s['flux'], re0=s['re'], q0=s['q'], pa0=s['pa'], fwhm=fw, extra=extra,
                             models=models, bic_margin=a.bic_margin, nucleus=a.nucleus, restarts=a.decomp_restarts, max_nfev=a.max_nfev, gain=a.gain or None)
            if d['res'] is None:
                raise RuntimeError('; '.join(d['notes']) or 'no decomposition converged')
            res, name = d['res'], d['name']
            out['decomp'] = d
        elif name == 'structure':
            if extra:
                raise RuntimeError('--model structure fits one galaxy: use --neighbours mask')
            sr = sf.fit_structure(cut, lx, ly, s['flux'], s['re'], s['q'], s['pa'], features=tuple(a.st_features.split(',')) if a.st_features != 'all' else ('bar', 'ring', 'spiral'),
                                  bic_margin=a.bic_margin, geometry=a.st_geometry, psf=psf, rms=rcut, mask=cm, sky=a.sky, zp=zp, gain=a.gain or None, max_nfev=a.max_nfev)
            if sr['res'] is None:
                raise RuntimeError('; '.join(sr['notes']) or 'no structure fit converged')
            res, name = sr['res'], 'structure:' + sr['name']
            out['structure'] = sr
        else:
            res, _ = mf.fit_preset(cut, name, lx, ly, s['flux'], s['re'], s['q'], s['pa'], extra=extra, **kw)
    except Exception as e:
        out['GF_FLAG'] = 2048
        out['err'] = str(e)
        return out
    ncomp_t = out['decomp']['ntarget'] if 'decomp' in out else len(res['components']) if 'structure' in out else len(mf.preset(name, 0, 0, 1, 1))
    comps = res['components'][:ncomp_t]
    for c in comps:
        if 'flux' not in c:                                              # surface-brightness normalised components (Ferrers bar): integrated flux from the magnitude
            c['flux'] = 10 ** (-0.4 * (c['mag'] - zp)) if c.get('mag') == c.get('mag') else 0.0
    flag |= res['flags']
    if res['chi2_red'] > a.chi2_max:
        flag |= mf.FLAGS['CHI2']
    tot = sum(max(c['flux'], 0) for c in comps)
    err2 = sum((c['errors'].get('flux', 0.0)) ** 2 for c in comps)
    ext = [c for c in comps if c['kind'] != 'psf']
    dom = max(ext, key=lambda c: c['flux']) if ext else comps[0]
    first = comps[0]
    # model of the target components only (for the full-frame model image), residual of the whole fit
    tmod = np.zeros(cut.shape)
    for c, cc in zip(comps, res['components']):
        tmod += mf.render_component(dict(c, flux=c['flux']), cut.shape, psf)
    good = res['good']
    resfrac = float(np.sum(np.abs(res['residual'][good])) / max(tot, 1e-9))
    nameid = name
    out.update(GF_X=comps[0]['x'] + cx0 + 1, GF_Y=comps[0]['y'] + cy0 + 1, GF_MAG=(zp - 2.5 * math.log10(tot)) if tot > 0 else None,
               GF_MAGERR=(1.0857 * math.sqrt(err2) / tot) if tot > 0 else None, GF_RE=dom.get('re') if dom['kind'] != 'psf' else None,
               GF_REERR=dom['errors'].get('re') if dom['kind'] != 'psf' else None, GF_N=dom.get('n') if dom['kind'] != 'psf' else None,
               GF_NERR=dom['errors'].get('n') if dom['kind'] == 'sersic' else None, GF_Q=dom.get('q') if dom['kind'] != 'psf' else None,
               GF_PA=dom.get('pa') if dom['kind'] != 'psf' else None, GF_BT=(first['flux'] / tot) if (len(comps) > 1 and tot > 0) else None,
               GF_MAG2=comps[1]['mag'] if len(comps) > 1 else None, GF_SKY=res['sky'], GF_CHI2=res['chi2_red'], GF_NCOMP=len(comps), GF_NNEIGH=nn, GF_RESFRAC=resfrac, GF_FLAG=flag,
               model_name=nameid, comps=[{k: v for k, v in c.items()} for c in comps], all_comps=[{k: v for k, v in c.items()} for c in res['components']], bbox=(cx0, cx1, cy0, cy1), tmod=tmod.astype(np.float32), nfev=res['nfev'], bic=res['bic'])
    if 'decomp' in out:
        out.update(decomp_row(out['decomp'], comps, tot, res, zp, flag))
        out['decomp_plot'] = decomp_plot_data(out['decomp'], comps, a.mag_zeropoint)
        out.pop('decomp')
    if 'structure' in out:
        out.update(structure_row(out['structure'], comps, tot, res, flag, zp))
        out.pop('structure')
    if a.keep_stamps:
        out['stamps'] = (cut.astype(np.float32), (res['model']).astype(np.float32), (res['residual']).astype(np.float32))
    return out


def structure_row(r, comps, tot, res, flag, zp=25.0):
    """ST_* catalog columns of one bar / ring / spiral fit (multifit.py --model structure)."""
    name = r['name']
    fl = [max(c['flux'], 0.0) for c in comps]
    row = dict(ST_MAG=(zp - 2.5 * math.log10(tot)) if tot > 0 else None, ST_TYPE={'base': 0, 'bar': 1, 'ring': 2, 'spiral': 3, 'bar+ring': 4}[name], ST_BT=fl[0] / tot if tot > 0 else None, ST_CHI2=res['chi2_red'], ST_FLAG=r['flag'] + 256 * (res['flags'] & 255),
               ST_DBIC=r['bic']['base'] - r['bic'][name] if 'base' in r['bic'] else None)
    if name in ('bar', 'bar+ring'):
        b = comps[2]
        row.update(ST_BARFRAC=fl[2] / tot, ST_BARLEN=b['rout'], ST_BARPA=b['pa'], ST_BARQ=b['q'], ST_BARC0=b.get('c0'))
    if name in ('ring', 'bar+ring'):
        b = comps[-1]
        row.update(ST_RINGR=b['rring'], ST_RINGW=b['sring'], ST_RINGFRAC=fl[-1] / tot)
    elif name == 'spiral':
        row.update(ST_SPTHETA=comps[1].get('rot_theta'), ST_SPAMP=comps[1].get('f2a'))
    row['ST_NOTE'] = '; '.join(r['notes'])
    return row


def decomp_row(d, comps, tot, res, zp, flag):
    """AD_* catalog columns and the rows of multifit_decomp.tsv for one decomposed object."""
    name = d['name']
    g = d.get('guess', {})
    f = d.get('features', {})
    row = dict(AD_TYPE=ad.TYPES[name], AD_MAG=(zp - 2.5 * math.log10(tot)) if tot > 0 else None, AD_CHI2=res['chi2_red'], AD_FLAG=d['flag'] + 256 * (res['flags'] & 255),
               AD_BAR=1 if f.get('bar') else 0, DEC_MODEL=name, DEC_NOTE='; '.join(d['notes']))
    b = d['bic']
    row['AD_DBIC'] = (b['sersic'] - b[name]) if 'sersic' in b else None
    if name == 'sersic':
        c = comps[0]
        row.update(AD_RE1=c.get('re'), AD_N1=c.get('n'))
    else:
        bulge, disc = comps[-2], comps[-1]
        fb, fd = max(bulge['flux'], 0), max(disc['flux'], 0)
        row.update(AD_BT=fb / (fb + fd) if fb + fd > 0 else None, AD_MAGB=bulge['mag'], AD_MAGD=disc['mag'], AD_REB=bulge['re'], AD_NB=bulge.get('n'), AD_QB=bulge['q'], AD_PAB=bulge['pa'],
                   AD_RED=disc['re'], AD_QD=disc['q'], AD_PAD=disc['pa'])
        if name == 'nucleus+bulge+disk':
            row['AD_MAGN'] = comps[0]['mag']
    gcs = g.get('comps', {}).get('bulge+disk') if g else None
    row.update(DEC_GUESS_BT=g.get('bt_guess') if g else None, DEC_GUESS_BT1D=g.get('bt_1d') if g else None, DEC_NUC_RATIO=g.get('nucleus_ratio') if g else None,
               DEC_GUESS_REB=gcs[0]['re'] if gcs else None, DEC_GUESS_NB=gcs[0]['n'] if gcs else None, DEC_GUESS_RED=gcs[1]['re'] if gcs else None,
               DEC_EPS_OUT=f.get('eps_out'), DEC_EPS_IN=f.get('eps_in'), DEC_CONC=f.get('conc'), DEC_R50=g.get('r50') if g else None, DEC_BIC_SERSIC=b.get('sersic'),
               DEC_BIC_BD=b.get('bulge+disk'), DEC_BIC_NBD=b.get('nucleus+bulge+disk'), DEC_SCALE=d.get('scale'), DEC_EDGEON=1 if f.get('edgeon') else 0, DEC_PROFILE=0 if (d['flag'] & ad.FLAGS['NOPROFILE']) else 1)
    return row


def decomp_plot_data(d, comps, zp):
    g = d.get('guess', {})
    gc = g.get('comps', {}).get('bulge+disk') if g else None
    pr = d.get('profile')
    return dict(name=d['name'], sma=None if pr is None else [float(v) for v in pr['sma']], I=None if pr is None else [float(v) for v in pr['I']], comps=[{k: c[k] for k in ('kind', 'flux', 're', 'n', 'q', 'pa') if k in c} for c in comps], guess=gc, bic=d['bic'], notes=d['notes'])


def pick_objects(rows, sv, a):
    idx = list(range(len(rows)))
    if a.objects:
        want = set(s.strip() for s in a.objects.replace(',', ';').split(';') if s.strip())
        idx = [i for i in idx if str(rows[i]['NUMBER']) in want]
    idx = [i for i in idx if sv[i]['flux'] > 0 and np.isfinite(sv[i]['x'])]
    idx.sort(key=lambda i: -sv[i]['flux'])
    if a.max_objects and not a.objects:
        idx = idx[:a.max_objects]
    return idx


def build_psf(data, mask, a):
    if a.psf_model and os.path.isfile(a.psf_model):
        return pm.load_model(a.psf_model), 'model:' + a.psf_model
    if a.psf and os.path.isfile(a.psf):
        from astropy.io import fits
        with fits.open(a.psf) as h:
            arr = next(np.asarray(x.data, float) for x in h if x.data is not None and x.data.ndim == 2)
        return pm.from_image(arr), 'image:' + a.psf
    if a.psf_fwhm > 0:
        return pm.analytic(max(15, int(math.ceil(7 * a.psf_fwhm)) | 1), a.psf_fwhm, 'gaussian'), 'gaussian %.2f px' % a.psf_fwhm
    mdl, info = pm.build_psf_model(data, mask=mask, fwhm_prior=3.0, snr_min=15, degree=1)
    return mdl, 'field stars (%s, %s stars)' % (mdl.meta.get('mode'), mdl.meta.get('n_stars'))


def montage(results, path, n=6):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        return
    res = [r for r in results if 'stamps' in r][:n]
    if not res:
        return
    fig, ax = plt.subplots(len(res), 3, figsize=(7.5, 2.4 * len(res)), squeeze=False)
    for k, r in enumerate(res):
        d, m, rr = r['stamps']
        v = np.nanpercentile(d, [2, 99.5])
        for j, (im, t, lim) in enumerate(((d, 'data', v), (m, 'model', v), (rr, 'residual', (-3 * np.nanstd(rr), 3 * np.nanstd(rr))))):
            ax[k][j].imshow(np.arcsinh(im - lim[0]) if j < 2 else im, origin='lower', cmap='gray', vmin=None if j < 2 else lim[0], vmax=None if j < 2 else lim[1])
            ax[k][j].set_xticks([]); ax[k][j].set_yticks([])
            if k == 0:
                ax[k][j].set_title(t, fontsize=9)
        ax[k][0].set_ylabel('#%s %s\nchi2 %.2f' % (r['NUMBER'], r.get('model_name', ''), r.get('GF_CHI2') or float('nan')), fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def load_config(a, hdr):
    """JSON config or GALFIT feedme (detected by content) -> config dict (1-based image coordinates)."""
    txt = open(a.config).read()
    if txt.lstrip().startswith('{'):
        return json.loads(txt)
    from ogfkit import galfitio
    ex = 1.0
    for k in ('EXPTIME', 'EXPOSURE', 'ITIME'):
        try:
            ex = float(hdr[k]); break
        except Exception:
            pass
    cfg = galfitio.parse_feedme(a.config, exptime=ex)
    for w in cfg.get('warnings', ()):
        sys.stderr.write('multifit: GALFIT feedme: %s\n' % w)
    if cfg.get('psf_sampling', 1) != 1:
        raise SystemExit('multifit: PSF fine sampling E) = %d is not supported (provide a PSF at the data pixel scale)' % cfg['psf_sampling'])
    return cfg


def export_feedme(a, data_path, comps_img, bbox, cfg, sky, grad, res, psf_file='none', tag='', tie=None, zp=25.0):
    """Write a GALFIT feedme (+ constraints when bounds / ties exist) that reproduces this fit's model; comps_img in 1-based image coordinates; bbox = (x0, x1, y0, y1) 0-based half-open."""
    from ogfkit import galfitio
    out = dict(cfg or {}, components=comps_img, bbox=list(bbox), sky_value=sky, sky_grad=list(grad), tie=tie or [], zp=zp)
    out['sky'] = (cfg or {}).get('sky', 'const')
    out.pop('region', None)
    shape = (bbox[3] - bbox[2], bbox[1] - bbox[0])
    cons = galfitio.write_constraints(out)
    has_cons = len(cons.strip().splitlines()) > 1
    base = a.export_feedme
    if not base.endswith('.feedme'):
        os.makedirs(base, exist_ok=True)
        base = os.path.join(base, 'galfit%s.feedme' % tag)
    cpath = os.path.splitext(base)[0] + '.constraints'
    text = galfitio.write_feedme(out, image=os.path.abspath(data_path), output='imgblock%s.fits' % tag, psf=psf_file, mask=(os.path.abspath(a.mask) if a.mask and os.path.isfile(a.mask) else 'none'),
                                 constraints=(cpath if has_cons else 'none'), shape=shape, zp=zp, exptime=(cfg or {}).get('exptime', 1.0), comment='chi2/dof %.3f' % res['chi2_red'] if res else '')
    with open(base, 'w') as fh:
        fh.write(text)
    if has_cons:
        with open(cpath, 'w') as fh:
            fh.write(cons)
    return base


def neighbour_mask(a, shape, cm, comps, x0, y0):
    """Add the ellipses (--neighbour-radius x A_IMAGE, B_IMAGE, THETA_IMAGE) of the --mask-catalog objects that fall in the fitted box to the pixel mask,
    except objects within --mask-exclude pixels of a fitted component (those are the targets).  -> (mask, number of masked objects)."""
    cols, rows = tsvio.read_catalog(a.mask_catalog)
    cm = np.zeros(shape, bool) if cm is None else cm.copy()
    pos = [(c['x'], c['y']) for c in comps if 'x' in c]
    n = 0
    for r in rows:
        x, y = num(r, 'X_IMAGE') - 1 - x0, num(r, 'Y_IMAGE') - 1 - y0
        if not (np.isfinite(x) and np.isfinite(y)):
            continue
        A, B = num(r, 'A_IMAGE', 2.0), num(r, 'B_IMAGE', 2.0)
        rad = max(a.neighbour_radius * max(A, 1.0), 1.0)
        if not (-rad <= x < shape[1] + rad and -rad <= y < shape[0] + rad):
            continue
        if any(math.hypot(x - px, y - py) <= max(a.mask_exclude, 0.0) for px, py in pos):
            continue
        cm |= ellipse_mask(shape, x, y, a.neighbour_radius * max(A, 1.0), a.neighbour_radius * max(B, 1.0), num(r, 'THETA_IMAGE', 0.0))
        n += 1
    return cm, n


def run_config(a, data, mask, cfg, hdr=None):
    comps = cfg['components']
    bb = cfg.get('bbox') or [0, data.shape[1], 0, data.shape[0]]
    x0, x1, y0, y1 = [int(v) for v in bb]
    cut = data[y0:y1, x0:x1].astype(float)
    cm = (mask[y0:y1, x0:x1] if mask is not None else None)
    one = cfg.get('one_based', True)
    comps_img = [json.loads(json.dumps(c)) for c in comps]              # untouched copy (image coordinates) for the export
    for c in comps:
        if 'x' in c:                                                   # truncation components may inherit the position of their parent
            c['x'] = c['x'] - x0 - 1 if one else c['x'] - x0
            c['y'] = c['y'] - y0 - 1 if one else c['y'] - y0
        for k, off in (('x', x0 + (1 if one else 0)), ('y', y0 + (1 if one else 0))):         # bounds on x / y are given in image coordinates like x / y
            if k in c.get('bounds', {}):
                c['bounds'][k] = [v - off for v in c['bounds'][k]]
    psf, desc = build_psf(data, mask, a)
    nmask = 0
    if a.mask_catalog:                                                  # neighbour masking: ellipses of the catalog objects inside the box, except those at a fitted component
        cm, nmask = neighbour_mask(a, cut.shape, cm, comps, x0, y0)
        desc += ', %d neighbour(s) masked' % nmask
    rms = cfg.get('rms') or None
    if cfg.get('sigma') and os.path.isfile(cfg['sigma']):
        rms = imageio.load_image(cfg['sigma'])[0][y0:y1, x0:x1].astype(float)
        desc += ', sigma image ' + os.path.basename(cfg['sigma'])
    if rms is None:
        rms = imageio.robust_sigma(cut, cm)
    g0 = cfg.get('sky_grad') or [0.0, 0.0]
    res = mf.fit_multistart(cut, comps, restarts=a.restarts, psf=psf, rms=rms, mask=cm, sky=cfg.get('sky', 'const'), sky_value=cfg.get('sky_value'), gain=cfg.get('gain') or a.gain or None,
                            zp=cfg.get('zp', a.mag_zeropoint), tie=cfg.get('tie'), max_nfev=a.max_nfev, sky_grad=g0)
    full_m = np.zeros(data.shape, np.float32); full_r = np.zeros(data.shape, np.float32)
    full_m[y0:y1, x0:x1] = res['model']; full_r[y0:y1, x0:x1] = res['residual']
    imageio.save_fits(os.path.join(a.work, 'multifit_model.fits'), full_m)
    imageio.save_fits(os.path.join(a.work, 'multifit_residual.fits'), full_r)
    lines = ['multi-component fit of %s, PSF: %s' % (os.path.basename(a.image), desc) + (' (start %d of %d won: chi2 %s)' % (res['start_used'], len(res['starts_chi2']), ', '.join('%.1f' % v for v in res['starts_chi2'])) if len(res.get('starts_chi2', ())) > 1 else ''),
             'chi2/dof = %.3f (dof %d), BIC %.1f, sky %.4g +- %.2g, flags %d, converged %s, nfev %d' % (res['chi2_red'], res['dof'], res['bic'], res['sky'], res['sky_err'], res['flags'], res['converged'], res['nfev'])]
    recs = []
    for k, c in enumerate(res['components']):
        er = c['errors']
        pos = ('x=%.3f+-%.3f y=%.3f+-%.3f ' % (c['x'] + x0 + 1, er.get('x', 0), c['y'] + y0 + 1, er.get('y', 0))) if 'x' in c else ''
        names = [n_ for n_ in mf.param_names(c) if n_ not in ('x', 'y', 'flux')]
        lines.append('  [%d] %-6s %smag=%.3f+-%.3f ' % (k, c['kind'], pos, c['mag'], er.get('mag', 0)) + ' '.join('%s=%.4g%s' % (n_, c[n_], '+-%.2g' % er[n_] if n_ in er else '') for n_ in names))
        extra = [n_ for n_ in names if n_ not in ('re', 'n', 'q', 'pa')]
        recs.append(dict(COMP=k, KIND=c['kind'], X=(c['x'] + x0 + 1) if 'x' in c else None, Y=(c['y'] + y0 + 1) if 'y' in c else None, MAG=c['mag'], MAGERR=er.get('mag'),
                         RE=c.get('re'), REERR=er.get('re'), N=c.get('n'), Q=c.get('q'), PA=c.get('pa'), EXTRA=';'.join('%s=%.6g' % (n_, c[n_]) for n_ in extra)))
    tsvio.write_table(os.path.join(a.work, 'multifit_results.tsv'), ['COMP', 'KIND', 'X', 'Y', 'MAG', 'MAGERR', 'RE', 'REERR', 'N', 'Q', 'PA', 'EXTRA'], recs)
    with open(os.path.join(a.work, 'multifit_config_result.json'), 'w') as fh:
        json.dump({k: v for k, v in res.items() if k not in ('model', 'residual', 'good')}, fh, default=lambda o: None)
    if a.export_feedme:
        fitted = []
        for c0, c in zip(comps_img, res['components']):
            d = dict(c0)
            for n_ in mf.param_names(c):
                if n_ in ('x', 'y'):
                    d[n_] = c[n_] + (x0 if n_ == 'x' else y0) + (1 if one else 0)
                elif n_ == 'flux':
                    d['mag'] = c['mag']
                else:
                    d[n_] = c[n_]
            if 'flux' in mf.param_names(c):
                d.pop('flux', None)
            else:
                d.pop('mag', None)
                d.pop('mu', None)
            fitted.append(d)
        for d, c0 in zip(fitted, comps_img):                                  # bounds were given in image coordinates: keep them
            d['bounds'] = c0.get('bounds', {})
        fn = export_feedme(a, a.image, fitted, (x0, x1, y0, y1), cfg, res['sky'], res['sky_grad'], res, psf_file=os.path.abspath(a.psf) if a.psf and os.path.isfile(a.psf) else (cfg.get('psf') or 'none'),
                           tie=cfg.get('tie'), zp=cfg.get('zp', a.mag_zeropoint))
        lines.append('GALFIT feedme written: %s' % fn)
    print('\n'.join(lines))
    return 0


DEC_COLS = ['NUMBER', 'DEC_MODEL', 'AD_TYPE', 'AD_MAG', 'AD_BT', 'AD_MAGB', 'AD_MAGD', 'AD_MAGN', 'AD_REB', 'AD_NB', 'AD_QB', 'AD_PAB', 'AD_RED', 'AD_QD', 'AD_PAD', 'AD_RE1', 'AD_N1', 'AD_DBIC', 'AD_CHI2',
            'AD_BAR', 'AD_FLAG', 'DEC_BIC_SERSIC', 'DEC_BIC_BD', 'DEC_BIC_NBD', 'DEC_SCALE', 'DEC_GUESS_BT', 'DEC_GUESS_BT1D', 'DEC_GUESS_REB', 'DEC_GUESS_NB', 'DEC_GUESS_RED', 'DEC_NUC_RATIO', 'DEC_EPS_OUT',
            'DEC_EPS_IN', 'DEC_CONC', 'DEC_R50', 'DEC_EDGEON', 'DEC_PROFILE', 'DEC_NOTE']


def write_decomp(a, W, results):
    rec = [r for r in results if 'DEC_MODEL' in r]
    rec.sort(key=lambda r: r['i'])
    tsvio.write_table(W('decomp.tsv'), DEC_COLS, rec)
    counts = {}
    for r in rec:
        counts[r['DEC_MODEL']] = counts.get(r['DEC_MODEL'], 0) + 1
    sys.stderr.write('autodecomp: %d objects: %s\n' % (len(rec), ', '.join('%s %d' % kv for kv in sorted(counts.items()))))
    try:
        decomp_plot(W('decomp_plot.png'), sorted([r for r in results if r.get('decomp_plot')], key=lambda r: r['i']), a.montage or 6)
    except Exception as e:
        sys.stderr.write('decomp plot failed: %s\n' % e)
    try:
        ometa.update(os.path.join(os.path.dirname(os.path.abspath(a.work)), 'catalog_meta.json'), 'autodecomp',
                     dict(n=len(rec), **{k.replace('+', '_'): v for k, v in counts.items()}), nrows=None)
    except Exception:
        pass


def decomp_plot(path, results, n=6):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    results = [r for r in results if r['decomp_plot'].get('sma')][:max(n, 1)]
    if not results:
        return
    nc = min(3, len(results))
    nr = (len(results) + nc - 1) // nc
    fig, ax = plt.subplots(nr, nc, figsize=(4.4 * nc, 3.4 * nr), squeeze=False)
    for k, r in enumerate(results):
        A = ax[k // nc][k % nc]
        d = r['decomp_plot']
        sma = np.array(d['sma']); I = np.array(d['I'])
        m = I > 0
        A.semilogy(sma[m], I[m], 'k.', ms=4, label='isophotes')
        rr = np.linspace(max(sma.min(), 0.3), sma.max(), 200)
        tot = np.zeros_like(rr)
        for c in d['comps']:
            if c['kind'] == 'sersic':
                ie = ad.sersic_ie_from_flux(c['flux'], c['re'], c['n'], c['q']); y = ad.sersic_profile(rr, ie, c['re'], c['n'])
            elif c['kind'] == 'exp':
                y = ad.exp_profile(rr, c['flux'] / (2 * math.pi * (c['re'] / 1.678) ** 2 * c['q']), c['re'] / 1.678)
            else:
                continue
            tot += y
            A.semilogy(rr, np.maximum(y, 1e-12), '-', lw=1, label='%s re=%.1f' % (c['kind'], c['re']))
        A.semilogy(rr, np.maximum(tot, 1e-12), 'r-', lw=1.5, label='sum (no PSF)')
        A.set_ylim(max(I[m].min() * 0.5, I.max() * 1e-5), I.max() * 2)
        A.set_title('#%s  %s  dBIC(sersic)=%.0f' % (r['NUMBER'], d['name'], (d['bic'].get('sersic', np.nan) - d['bic'][d['name']])), fontsize=8)
        A.set_xlabel('sma [px]', fontsize=7); A.legend(fontsize=6)
    for k in range(len(results), nr * nc):
        ax[k // nc][k % nc].axis('off')
    fig.tight_layout(); fig.savefig(path, dpi=75); plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--work', default='.')
    ap.add_argument('--catalog', default='')
    ap.add_argument('--config', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--psf-model', default='')
    ap.add_argument('--psf', default='')
    ap.add_argument('--psf-fwhm', type=float, default=0.0)
    ap.add_argument('--model', default='sersic', choices=['sersic', 'exp', 'dev', 'psf', 'bulge+disk', 'psf+sersic', 'auto', 'decomp', 'bar', 'ring', 'spiral', 'bulge+disk+bar', 'bulge+disk+ring', 'structure'])
    ap.add_argument('--neighbours', default='fit', choices=['fit', 'mask', 'ignore'])
    ap.add_argument('--objects', default='')
    ap.add_argument('--max-objects', type=int, default=100)
    ap.add_argument('--max-neighbours', type=int, default=6)
    ap.add_argument('--neighbour-radius', type=float, default=2.5, help='masked neighbours: ellipse of this many A_IMAGE x B_IMAGE (sigma-like semi-axes)')
    ap.add_argument('--size-scale', type=float, default=6.0, help='cutout half-size = size-scale x A_IMAGE')
    ap.add_argument('--min-size', type=int, default=31)
    ap.add_argument('--max-size', type=int, default=121)
    ap.add_argument('--sky', default='const', choices=['const', 'plane', 'fixed'])
    ap.add_argument('--gain', type=float, default=0.0)
    ap.add_argument('--rms', type=float, default=0.0, help='background sigma (0 = measured from the image)')
    ap.add_argument('--chi2-max', type=float, default=3.0)
    ap.add_argument('--bic-margin', type=float, default=10.0)
    ap.add_argument('--max-nfev', type=int, default=150)
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--n-workers', type=int, default=0)
    ap.add_argument('--montage', type=int, default=6)
    ap.add_argument('--no-residual', action='store_true')
    ap.add_argument('--restarts', type=int, default=2, help='config / feedme fits: number of extra starts (R_e and n rescaled) tried after the given one; the lowest chi2 is kept (0 = fit only from the given start)')
    ap.add_argument('--mask-catalog', default='', help='config / feedme fits: neighbour masking - TSV catalog (X_IMAGE Y_IMAGE A_IMAGE B_IMAGE THETA_IMAGE); the ellipses of its objects inside the fitted box are masked, except those within --mask-exclude px of a fitted component')
    ap.add_argument('--mask-exclude', type=float, default=4.0, help='--mask-catalog: objects closer than this (pixels) to a fitted component centre are not masked')
    ap.add_argument('--columns', default='gf', choices=['gf', 'ad', 'st'], help='catalog columns of the result: GF_* (components of the fit) or AD_* (structure decomposition: type, B/T, bulge and disc parameters)')
    ap.add_argument('--nucleus', default='auto', choices=['auto', 'always', 'off'], help='--model decomp: try a point-source nucleus (auto = when the profile shows a central excess)')
    ap.add_argument('--st-geometry', default='given', choices=['given', 'outer'], help='--model structure: start the disc axis ratio / PA from the catalog (given) or from the outer isophotes (outer; for galaxies whose bar / ring dominate the inner light)')
    ap.add_argument('--st-features', default='all', help='--model structure: features to try (all or a comma list of bar,ring,spiral)')
    ap.add_argument('--decomp-restarts', type=int, default=1, help='--model decomp: extra starts (R_e / n rescaled) per candidate model')
    ap.add_argument('--export-feedme', default='', help='write the fitted model as a GALFIT feedme (file; with --config) or into this directory (catalog mode: one galfit_<NUMBER>.feedme per fitted object)')
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, 'multifit_' + n)
    cfg = None
    if a.config and a.image in ('-', 'feedme'):                              # image / PSF / mask named in the feedme
        from ogfkit import galfitio
        cfg0 = galfitio.parse_feedme(a.config)
        a.image = cfg0['input']
        if cfg0.get('psf') and not (a.psf or a.psf_model):
            a.psf = cfg0['psf']
        if cfg0.get('mask') and not a.mask:
            a.mask = cfg0['mask']
    data, hdr = imageio.load_image(a.image)
    mask = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    if a.config:
        cfg = load_config(a, hdr)
        if cfg.get('galfit'):
            if cfg.get('psf') and os.path.isfile(cfg['psf']) and not (a.psf or a.psf_model or a.psf_fwhm > 0):
                a.psf = cfg['psf']
            if cfg.get('mask') and os.path.isfile(cfg['mask']) and mask is None:
                mask = imageio.load_mask(cfg['mask'], data.shape)
                a.mask = cfg['mask']
        return run_config(a, data, mask, cfg, hdr)
    cols, rows = tsvio.read_catalog(a.catalog)
    zp = a.mag_zeropoint
    sv = [start_values(r, zp) for r in rows]
    psf, desc = build_psf(data, mask, a)
    try:
        psf.save(W('psf.json'))
    except Exception:
        pass
    rms = a.rms if a.rms > 0 else imageio.robust_sigma(data, mask)
    from scipy.spatial import cKDTree
    pts = np.array([[s['x'] if np.isfinite(s['x']) else -1e9, s['y'] if np.isfinite(s['y']) else -1e9] for s in sv]) if sv else np.zeros((0, 2))
    idx = pick_objects(rows, sv, a)
    sys.stderr.write('multifit: %d objects of %d, model %s, neighbours %s, PSF %s, rms %.4g\n' % (len(idx), len(rows), a.model, a.neighbours, desc, rms))
    G.update(args=a, rows=rows, sv=sv, data=data, mask=mask, rms=rms, psf=psf, tree=cKDTree(pts) if len(pts) else None,
             psf_sigma=max(psf.fwhm_estimate() / 2.355, 0.8) if hasattr(psf, 'fwhm_estimate') else 1.5)
    a.keep_stamps = a.montage > 0
    nw = a.n_workers if a.n_workers > 0 else min(8, os.cpu_count() or 1)
    if nw > 1 and len(idx) > 3:
        with mp.get_context('fork').Pool(nw) as pool:
            results = pool.map(fit_object, idx, chunksize=1)
    else:
        results = [fit_object(i) for i in idx]
    by = {r['i']: r for r in results}
    # outputs
    out_rows = []
    for i, r in enumerate(rows):
        res = by.get(i)
        out_rows.append((r['NUMBER'], {c: (res.get(c) if res else None) for c in (AD_COLUMNS if a.columns == 'ad' else ST_COLUMNS if a.columns == 'st' else COLUMNS)}))
    full_m = np.zeros(data.shape, np.float32)
    recs = []
    for r in results:
        if 'bbox' not in r:
            continue
        x0, x1, y0, y1 = r['bbox']
        full_m[y0:y1, x0:x1] += r['tmod']
        for k, c in enumerate(r['comps']):
            er = c['errors']
            recs.append(dict(NUMBER=r['NUMBER'], COMP=k, KIND=c['kind'], X=c['x'] + x0 + 1, Y=c['y'] + y0 + 1, MAG=c['mag'], MAGERR=er.get('mag'), RE=c.get('re') if c['kind'] != 'psf' else None,
                             REERR=er.get('re'), N=c.get('n') if c['kind'] != 'psf' else None, NERR=er.get('n'), Q=c.get('q') if c['kind'] != 'psf' else None, PA=c.get('pa') if c['kind'] != 'psf' else None,
                             CHI2=r['GF_CHI2'], FLAG=r['GF_FLAG'], MODEL=r['model_name']))
    tsvio.write_table(W('results.tsv'), ['NUMBER', 'COMP', 'KIND', 'MODEL', 'X', 'Y', 'MAG', 'MAGERR', 'RE', 'REERR', 'N', 'NERR', 'Q', 'PA', 'CHI2', 'FLAG'], recs)
    if not a.no_residual:
        imageio.save_fits(W('model.fits'), full_m)
        imageio.save_fits(W('residual.fits'), (data - full_m).astype(np.float32))
    if a.export_feedme:
        os.makedirs(a.export_feedme, exist_ok=True)
        nexp = 0
        for r in results:
            if 'bbox' not in r:
                continue
            x0, x1, y0, y1 = r['bbox']
            cs = []
            for c in r['all_comps']:
                d = dict(kind=c['kind'], x=c['x'] + x0 + 1, y=c['y'] + y0 + 1, mag=c['mag'], fixed=sorted(c.get('fixed', ())))
                if c['kind'] != 'psf':
                    d.update(re=c['re'], q=c['q'], pa=c['pa'])
                    if c['kind'] == 'sersic':
                        d['n'] = c['n']
                cs.append(d)
            tag = '_%s' % r['NUMBER']
            pf = os.path.join(a.export_feedme, 'galfit%s_psf.fits' % tag)
            try:
                st = mf.psf_stamp(psf, cs[0]['x'] - 1, cs[0]['y'] - 1)
                imageio.save_fits(pf, np.asarray(st, np.float32))
            except Exception:
                pf = 'none'
            export_feedme(a, a.image, cs, (x0, x1, y0, y1), None, r['GF_SKY'], (0.0, 0.0), dict(chi2_red=r['GF_CHI2']), psf_file=os.path.abspath(pf) if pf != 'none' else 'none', tag=tag, zp=zp)
            nexp += 1
        sys.stderr.write('multifit: %d GALFIT feedme files written to %s (target + fitted neighbours; masked neighbours are not exported)\n' % (nexp, a.export_feedme))
    montage(sorted([r for r in results if 'stamps' in r], key=lambda r: r['i']), W('montage.png'), a.montage)
    ok = [r for r in results if r.get('GF_FLAG', 0) == 0 or 'bbox' in r]
    nconv = sum(1 for r in results if 'bbox' in r and not (r['GF_FLAG'] & (mf.FLAGS['NOCONV'] | mf.FLAGS['CHI2'])))
    sys.stderr.write('multifit: %d fitted, %d clean (converged, chi2_red <= %.1f); residual/model: %s\n' % (sum(1 for r in results if 'bbox' in r), nconv, a.chi2_max, W('residual.fits')))
    try:
        ometa.update(os.path.join(os.path.dirname(os.path.abspath(a.work)), 'catalog_meta.json'), 'multifit',
                     dict(model=a.model, neighbours=a.neighbours, n_fitted=len(results), n_clean=nconv, psf=desc), nrows=len(rows))
    except Exception:
        pass
    if any('DEC_MODEL' in r for r in results):
        write_decomp(a, W, results)
    tsvio.write_columns(AD_COLUMNS if a.columns == 'ad' else ST_COLUMNS if a.columns == 'st' else COLUMNS, out_rows)
    return 0


if __name__ == '__main__':
    sys.exit(main())
