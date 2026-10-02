"""Synthetic strong lenses with known parameters (truth for the lensmodel validation).

random_lens(rng) -> dict(params, source, images) for an SIE + external shear lens and a point source inside the caustic (quad) or outside the inner
astroid (double).  render_arcs() paints a lensed Sersic source on a pixel grid (with Gaussian PSF and noise) for the source-plane reconstruction test.
"""
import math

import numpy as np

from . import lensmodel as L


def random_lens(rng, quad=True, theta_E=None, max_tries=200):
    for _ in range(max_tries):
        p = dict(theta_E=float(theta_E if theta_E else rng.uniform(0.9, 2.2)), q=float(rng.uniform(0.55, 0.95)), phi=float(rng.uniform(0, 180)),
                 gamma=float(rng.uniform(0.0, 0.10)), phi_g=float(rng.uniform(0, 180)), x0=0.0, y0=0.0)
        m = L.build_sie_shear(p)
        cc = m.critical_curves(half=3 * p['theta_E'], n=401)
        tang = [c for c in cc if c['kind'] == 'tangential']
        if not tang:
            continue
        cx, cy = tang[0]['cx'], tang[0]['cy']
        # caustic extent (astroid) -> pick a source inside it for a quad, between astroid and cut for a double
        rc = np.hypot(cx - np.mean(cx), cy - np.mean(cy)).max()
        for _t in range(60):
            r = rc * math.sqrt(rng.uniform(0.02, 0.5 if quad else 1.0)) if quad else rc * rng.uniform(1.15, 1.8)
            a = rng.uniform(0, 2 * math.pi)
            bx, by = r * math.cos(a), r * math.sin(a)
            ims = m.solve_images(bx, by, half=3 * p['theta_E'], n=301, min_abs_mu=0.0)
            if (quad and len(ims) == 4) or ((not quad) and len(ims) == 2):
                return dict(params=p, source=(bx, by), images=ims)
    raise RuntimeError('no lens found')


def sersic_source(bx, by, x0, y0, reff, n=1.0, q=1.0, phi=0.0, amp=1.0):
    k = 1.9992 * n - 0.3271
    c, s = math.cos(math.radians(phi)), math.sin(math.radians(phi))
    dx, dy = bx - x0, by - y0
    xp, yp = c * dx + s * dy, -s * dx + c * dy
    r = np.sqrt(xp * xp + (yp / q) ** 2)
    return amp * np.exp(-k * ((r / reff) ** (1.0 / n) - 1.0))


def render_arcs(model, shape, pixscale, origin_xy, centre_xy, src, psf_sigma_pix=1.5, noise=0.0, seed=0, oversample=3):
    """lensed image of an analytic Sersic source ``src`` = dict(x0, y0, reff, n, q, phi, amp); oversampled, PSF-convolved, + Gaussian noise"""
    from scipy.ndimage import gaussian_filter
    ny, nx = shape
    o = oversample
    yy, xx = np.mgrid[0:ny * o, 0:nx * o]
    px = (xx + 0.5) / o - 0.5
    py = (yy + 0.5) / o - 0.5
    tx = (px - origin_xy[0]) * pixscale + centre_xy[0]
    ty = (py - origin_xy[1]) * pixscale + centre_xy[1]
    bx, by = model.beta(tx, ty)
    im = sersic_source(bx, by, src['x0'], src['y0'], src['reff'], src.get('n', 1.0), src.get('q', 1.0), src.get('phi', 0.0), src.get('amp', 1.0))
    im = im.reshape(ny, o, nx, o).mean(axis=(1, 3))
    im = gaussian_filter(im, psf_sigma_pix)
    if noise > 0:
        im = im + np.random.default_rng(seed).normal(0, noise, im.shape)
    return im


def validation_sample(n_lens=40, noise_arcsec=0.003, quad_fraction=0.75, seed=1, fit_kw=None):
    """Monte-Carlo: random SIE+shear lenses, images with Gaussian position noise, fit with the true lens centre fixed.
    -> list of per-lens dicts with truth, fit and differences"""
    rng = np.random.default_rng(seed)
    out = []
    for i in range(n_lens):
        quad = rng.uniform() < quad_fraction
        lens = random_lens(rng, quad=quad)
        ims = lens['images']
        obs = np.array([[im['x'], im['y']] for im in ims])
        obs_n = obs + rng.normal(0, noise_arcsec, obs.shape)
        kw = dict(fit_kw or {})
        if not quad:
            kw.setdefault('free', ['theta_E'])
            kw.setdefault('fixed', dict(q=1.0, gamma=0.0))
        fit = L.fit_sie_shear(obs_n, max(noise_arcsec, 1e-4), (0.0, 0.0), seed=i, **kw)
        out.append(dict(quad=bool(quad), truth=lens['params'], source=lens['source'], n_images=len(ims), fit=fit, obs=obs_n.tolist(), mu=[im['mu'] for im in ims]))
    return out


def angle_diff(a, b, period=180.0):
    d = (a - b + period / 2) % period - period / 2
    return d


def source_recon_test(seed=0, theta_E_error=0.0, noise=0.01, pixscale=0.05, n=241):
    """Render a lensed Sersic source with the true SIE+shear model, reconstruct it with a (possibly biased) model, compare with the analytic source.
    -> dict(corr, flux_ratio, resid_rms_over_noise, centroid_err_arcsec)"""
    rng = np.random.default_rng(seed)
    lens = random_lens(rng, quad=True)
    p = lens['params']
    bx, by = lens['source']
    ref = ((n - 1) / 2.0, (n - 1) / 2.0)
    truth = L.build_sie_shear(p)
    src = dict(x0=bx, y0=by, reff=0.10, n=1.0, q=0.7, phi=30.0, amp=1.0)
    img = render_arcs(truth, (n, n), pixscale, ref, (0.0, 0.0), src, psf_sigma_pix=0.0, noise=noise, seed=seed + 1)
    pm = dict(p)
    pm['theta_E'] = p['theta_E'] * (1 + theta_E_error)
    model = L.build_sie_shear(pm)
    rec = L.ray_trace_image(model, img.astype(float), pixscale, ref, (0.0, 0.0), src_pix=0.5 * pixscale)
    ny, nx = rec['src'].shape
    gx = rec['x0'] + np.arange(nx) * rec['pix']
    gy = rec['y0'] + np.arange(ny) * rec['pix']
    GX, GY = np.meshgrid(gx, gy)
    tr = sersic_source(GX, GY, src['x0'], src['y0'], src['reff'], 1.0, src['q'], src['phi'], src['amp'])
    good = np.isfinite(rec['src']) & (rec['count'] >= 3) & (tr > 0.05)          # where the source is bright and well sampled
    corr = float(np.corrcoef(rec['src'][good], tr[good])[0, 1])
    flux_ratio = float(rec['src'][good].sum() / tr[good].sum())
    fwd = L.lens_source_map(model, rec, (n, n), pixscale, ref, (0.0, 0.0))
    bright = img > 5 * noise
    resid = float(np.sqrt(np.mean((img - fwd)[bright] ** 2)) / noise)
    s = np.nan_to_num(rec['src'], nan=0.0) * (rec['count'] >= 3)
    w = s.clip(0)
    cx = float((w * GX).sum() / w.sum())
    cy = float((w * GY).sum() / w.sum())
    return dict(corr=corr, flux_ratio=flux_ratio, resid_over_noise=resid, centroid_err=float(math.hypot(cx - bx, cy - by)), n_good=int(good.sum()), theta_E=p['theta_E'])
