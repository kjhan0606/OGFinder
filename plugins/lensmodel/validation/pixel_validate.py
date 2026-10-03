"""Validation of the pixelated source inversion / pixel-level lens fit on synthetic lensed Sersic + clump sources."""
import os, sys, json, math, time
import numpy as np
ROOT = os.environ.get('OGF_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
sys.path.insert(0, ROOT)
from ogfkit import lensmodel as L, lenssynth as S, lensextra as X


def true_source(bx, by, src):
    t = S.sersic_source(bx, by, src['x0'], src['y0'], src['reff'], src.get('n', 1.0), src.get('q', 1.0), src.get('phi', 0.0), src.get('amp', 1.0))
    for c in src.get('clumps', []):
        t = t + S.sersic_source(bx, by, c['x0'], c['y0'], c['reff'], 1.0, 1.0, 0.0, c['amp'])
    return t


def make_case(rng, noise=0.02, npix=120, pixscale=0.05, psf=1.5, clumps=True):
    d = S.random_lens(rng, quad=bool(rng.random() < 0.7), theta_E=rng.uniform(1.0, 1.6))
    m = L.build_sie_shear(d['params'])
    bx, by = d['source']
    src = dict(x0=bx, y0=by, reff=0.12, n=1.0, q=rng.uniform(0.6, 0.9), phi=rng.uniform(0, 180), amp=1.0,
               clumps=[dict(x0=bx + 0.12 * math.cos(a), y0=by + 0.12 * math.sin(a), reff=0.03, amp=1.5) for a in rng.uniform(0, 6.28, 2)] if clumps else [])
    org = ((npix - 1) / 2.0, (npix - 1) / 2.0)
    o = 3
    yy, xx = np.mgrid[0:npix * o, 0:npix * o]
    px = (xx + 0.5) / o - 0.5; py = (yy + 0.5) / o - 0.5
    tx = (px - org[0]) * pixscale; ty = (py - org[1]) * pixscale
    sbx, sby = m.beta(tx, ty)
    im = true_source(sbx, sby, src).reshape(npix, o, npix, o).mean(axis=(1, 3))
    from scipy.ndimage import gaussian_filter
    im = gaussian_filter(im, psf)
    img = im + rng.normal(0, noise, im.shape)
    return dict(lens=d['params'], src=src, img=img, clean=im, org=org, pixscale=pixscale, psf=psf, noise=noise, model=m)


def circ_mask(shape, org, r):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    return np.hypot(xx - org[0], yy - org[1]) < r


def one(rng, noise, fit_lens=True, clumps=True, dlne=True):
    c = make_case(rng, noise, clumps=clumps)
    shape = c['img'].shape
    mask = X.arc_mask(c['img'], noise, c['org'], 2.4 * c['lens']['theta_E'] / c['pixscale'], rmin_pix=0.35 * c['lens']['theta_E'] / c['pixscale'])
    k = X.gaussian_kernel(c['psf'])
    out = {}
    SI = X.SourceInversion(c['model'], c['img'], noise, mask, c['pixscale'], c['org'], (0, 0), kernel=k, n=36, reg='gradient')
    r = SI.best()
    src, mod = SI.images(r)
    g = SI.grid; n = g['n']
    xs = g['x0'] + g['pix'] * np.arange(n)
    ys = g['y0'] + g['pix'] * np.arange(n)
    BX, BY = np.meshgrid(xs, ys)
    ts = true_source(BX, BY, c['src'])
    # compare where the data constrain (>= 3 image-pixel rays), at the resolution of the data: both maps smoothed with a 0.05 arcsec Gaussian
    from scipy.ndimage import gaussian_filter
    L_ = X.lens_operator(c['model'], shape, c['pixscale'], c['org'], (0, 0), g)
    cover = np.asarray((L_[np.flatnonzero(mask.ravel())] > 0).sum(axis=0)).reshape(n, n) >= 3

    def metrics(rec, ts_, pix):
        sg = 0.05 / pix
        a_ = gaussian_filter(np.where(cover, np.nan_to_num(rec), 0.0), sg)
        t_ = gaussian_filter(np.where(cover, ts_, 0.0), sg)
        sel_ = cover & (t_ > 0.1 * t_.max())
        cor = float(np.corrcoef(a_[sel_], t_[sel_])[0, 1]) if sel_.sum() > 5 else float('nan')
        yy_, xx_ = np.mgrid[0:n, 0:n] if a_.shape == (n, n) else np.mgrid[0:a_.shape[0], 0:a_.shape[1]]
        fr_ = float(rec[cover].sum() / ts_[cover].sum()) if False else float(np.nan_to_num(rec)[cover].sum() / ts_[cover].sum())
        cen_t = (t_ * xx_).sum() / t_.sum(), (t_ * yy_).sum() / t_.sum()
        cen_a = (np.clip(a_, 0, None) * xx_).sum() / np.clip(a_, 0, None).sum(), (np.clip(a_, 0, None) * yy_).sum() / np.clip(a_, 0, None).sum()
        return cor, fr_, float(math.hypot(cen_a[0] - cen_t[0], cen_a[1] - cen_t[1]) * pix), float(np.sqrt(np.mean((a_[sel_] - t_[sel_]) ** 2)) / t_.max())
    cor, fr, cen, rms = metrics(src, ts, g['pix'])
    neff = SI.n_eff(r)
    out['inv'] = dict(corr=cor, flux_ratio=fr, centroid_err_arcsec=cen, rms_over_peak=rms, chi2=r['chi2'], ndata=SI.ndata, neff=neff, chi2_red=r['chi2'] / (SI.ndata - neff), lam=r['lam'], pix=g['pix'])
    rt = L.ray_trace_image(c['model'], np.where(mask, c['img'], np.nan), c['pixscale'], c['org'], (0, 0), src_pix=g['pix'])
    # back-projection map lives on its own grid: resample the truth on that grid, same cover logic from counts
    xs2 = rt['x0'] + rt['pix'] * np.arange(rt['src'].shape[1]); ys2 = rt['y0'] + rt['pix'] * np.arange(rt['src'].shape[0]); BX2, BY2 = np.meshgrid(xs2, ys2)
    ts2 = true_source(BX2, BY2, c['src']); s2 = np.nan_to_num(rt['src'])
    cov2 = rt['count'] >= 3
    sg2 = 0.05 / rt['pix']
    a2 = gaussian_filter(np.where(cov2, s2, 0.0), sg2); t2 = gaussian_filter(np.where(cov2, ts2, 0.0), sg2)
    sel2 = cov2 & (t2 > 0.1 * t2.max())
    out['backproj_corr'] = float(np.corrcoef(a2[sel2], t2[sel2])[0, 1]) if sel2.sum() > 5 else float('nan')
    # evidence: true vs perturbed lens
    out['quad'] = len(c['model'].solve_images(c['src']['x0'], c['src']['y0'], half=4.0, n=201)) >= 4
    for tag, dp in ((('theta_E+2%', dict(theta_E=1.02)), ('q-0.05', dict(q=-0.05))) if dlne else ()):
        p = dict(c['lens'])
        for kk, v in dp.items():
            p[kk] = p[kk] * v if kk == 'theta_E' else p[kk] + v
        m2 = L.build_sie_shear(p)
        S2 = X.SourceInversion(m2, c['img'], noise, mask, c['pixscale'], c['org'], (0, 0), kernel=k, n=36, reg='gradient')
        out['dlnE_' + tag] = float(S2.best()['evidence'] - r['evidence'])
    if fit_lens:
        st = dict(c['lens'])
        st['theta_E'] *= 1.03; st['q'] = min(0.99, st['q'] + 0.05); st['phi'] += 5; st['gamma'] = max(0.0, st['gamma'] + 0.02)
        t0 = time.time()
        f = X.fit_lens_pixels(c['img'], noise, mask, c['pixscale'], c['org'], (0, 0), st, kernel=k, n=28, maxiter=250)
        p = f['params']
        out['lensfit'] = dict(dthE=(p['theta_E'] - c['lens']['theta_E']) / c['lens']['theta_E'], dq=p['q'] - c['lens']['q'], dphi=(p['phi'] - c['lens']['phi'] + 90) % 180 - 90,
                              dgam=p['gamma'] - c['lens']['gamma'], nfev=f['nfev'], sec=time.time() - t0, start_dthE=0.03)
    return out


if __name__ == '__main__':
    outp = sys.argv[1]; n = int(sys.argv[2]); noise = float(sys.argv[3]) if len(sys.argv) > 3 else 0.02
    rng = np.random.default_rng(int(sys.argv[4]) if len(sys.argv) > 4 else 11)
    rows = []
    for i in range(n):
        rows.append(one(rng, noise)); print(i, json.dumps(rows[-1])[:300], flush=True)
    json.dump(dict(noise=noise, rows=rows), open(outp, 'w'), indent=1)
