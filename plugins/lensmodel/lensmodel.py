#!/usr/bin/env python3
"""Strong-lens modelling for OGFinder.

    lensmodel.py --task fit     --catalog TSV --image FITS --work DIR --image-numbers 3,4,7,9 --lens-number 1 [--z-lens 0.5 --z-source 2.0]
    lensmodel.py --task curves  --work DIR [--image FITS]          (critical curves + caustics -> lens_curves.json, lens_overlay.reg)
    lensmodel.py --task magmap  --work DIR --image FITS            (magnification map FITS on the image grid)
    lensmodel.py --task source  --work DIR --image FITS            (source-plane reconstruction + forward model + residual)
    lensmodel.py --task source  --work DIR --image FITS --source-method inversion [--psf-sigma PX --noise SIG --regularisation gradient]
    lensmodel.py --task multiplane --work DIR --pixscale S --lens-x X --lens-y Y --mp-groups "x,y;x,y;x,y;x,y|x,y;..." --mp-z 0.6,2.0 --z-lens 0.22
    lensmodel.py --task predict --work DIR [--source-x PX --source-y PY]   (counter-images, magnifications, time delays)
    lensmodel.py --task simulate --work DIR [--seed N ...]         (synthetic lens: FITS + catalog + truth)

`fit` follows the add_columns contract (stdout = NUMBER + LENS_* columns for the image rows) and writes DIR/lens_model.json, lens_images.tsv,
lens_curves.json, lens_overlay.reg (ds9 regions, image coordinates).  Pixel coordinates in/out are 1-based (X_IMAGE / Y_IMAGE, ds9 image).
The model lives in the PIXEL frame (x right, y up, arcsec = pixel * scale); position angles are measured from +x counter-clockwise, and the sky PA (east of
north) is added when the header carries a WCS.  See docs/lensmodel.md.
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

from ogfkit import tsvio, imageio, lensmodel as lm, lenssynth  # noqa: E402

COLUMNS = ['LENS_ROLE', 'LENS_MU', 'LENS_RES', 'LENS_DT', 'LENS_PARITY']


def pa_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--task', required=True, choices=['fit', 'curves', 'magmap', 'source', 'predict', 'simulate', 'multiplane'])
    p.add_argument('--catalog', default='')
    p.add_argument('--image', default='')
    p.add_argument('--work', required=True)
    p.add_argument('--meta-out', default='')
    p.add_argument('--pixscale', type=float, default=0.0, help='arcsec / pixel (0 = from the FITS WCS)')
    p.add_argument('--image-numbers', default='', help='catalog NUMBERs of the multiple images')
    p.add_argument('--positions', default='', help='"x1,y1;x2,y2;..." 1-based pixels (used when --image-numbers is empty)')
    p.add_argument('--lens-number', type=int, default=0, help='catalog NUMBER of the lens galaxy (its X_IMAGE/Y_IMAGE is the lens centre)')
    p.add_argument('--lens-x', type=float, default=0.0)
    p.add_argument('--lens-y', type=float, default=0.0)
    p.add_argument('--sigma-pos', type=float, default=0.01, help='image position uncertainty (arcsec)')
    p.add_argument('--z-lens', type=float, default=0.0)
    p.add_argument('--z-source', type=float, default=0.0)
    p.add_argument('--h0', type=float, default=70.0)
    p.add_argument('--omega-m', type=float, default=0.3)
    p.add_argument('--free', default='auto', help='comma list of theta_E,q,phi,gamma,phi_g,x0,y0 or auto (by number of images)')
    p.add_argument('--fixed', default='', help='"q=1,gamma=0" values of non-free parameters')
    p.add_argument('--n-starts', type=int, default=8)
    p.add_argument('--flux-column', default='', help='catalog column with image fluxes (adds a magnification-ratio constraint)')
    p.add_argument('--flux-sigma', type=float, default=0.3, help='sigma of ln(flux ratio)')
    p.add_argument('--nfw-kappa-s', type=float, default=0.0, help='fixed NFW halo component (0 = off)')
    p.add_argument('--nfw-theta-s', type=float, default=20.0)
    p.add_argument('--member-column', default='', help='catalog column flagging cluster members (> 0), e.g. CL_MEMBER')
    p.add_argument('--member-mag-column', default='MAG_AUTO')
    p.add_argument('--member-mag-star', type=float, default=20.0)
    p.add_argument('--member-b-star', type=float, default=0.1, help='PIEMD b0 (arcsec) of an L* member (b0 ~ L^0.5)')
    p.add_argument('--member-core', type=float, default=0.05)
    p.add_argument('--member-cut', type=float, default=10.0)
    p.add_argument('--source-x', type=float, default=float('nan'), help='predict: source position (1-based pixels of the source-plane frame, same origin)')
    p.add_argument('--source-y', type=float, default=float('nan'))
    p.add_argument('--grid', type=int, default=401, help='lens-equation grid size')
    p.add_argument('--source-pix', type=float, default=0.0, help='source-plane pixel (arcsec; 0 = half the image pixel)')
    p.add_argument('--cutout', type=float, default=0.0, help='source: image radius (arcsec) around the lens centre (0 = 2.5 theta_E)')
    p.add_argument('--lens-mask', type=float, default=0.0, help='source: mask pixels closer than this (arcsec) to the lens centre (lens-galaxy light must be removed or masked)')
    p.add_argument('--source-method', choices=['backproject', 'inversion'], default='backproject',
                   help='source: back-projection (default) or regularised linear inversion on a pixel grid with PSF and noise (ogfkit.lensextra)')
    p.add_argument('--psf-sigma', type=float, default=0.0, help='source/inversion: Gaussian PSF sigma (image pixels)')
    p.add_argument('--noise', type=float, default=0.0, help='source/inversion: pixel noise sigma (0 = rms of the blank pixels)')
    p.add_argument('--regularisation', choices=['gradient', 'curvature', 'zero'], default='gradient')
    p.add_argument('--source-n', type=int, default=40, help='source/inversion: source grid size (n x n)')
    p.add_argument('--mp-groups', default='', help='multiplane: image sets "x,y;x,y;...|x,y;..." (1-based pixels), first set = reference source plane')
    p.add_argument('--mp-z', default='', help='multiplane: source redshift of each set, comma list; "free" fits the weight of that plane')
    p.add_argument('--max-mu', type=float, default=100.0, help='magmap: |mu| is clipped to this value')
    p.add_argument('--check-lenstronomy', action='store_true', help='compare the deflection with lenstronomy when it is installed')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--sim-theta-e', type=float, default=1.5)
    p.add_argument('--sim-quad', type=int, default=1)
    p.add_argument('--sim-noise', type=float, default=0.005, help='simulate: Gaussian pixel noise of the image')
    p.add_argument('--sim-lens-light', type=int, default=0, help='simulate: add a bright lens galaxy (then use --lens-mask in the source task)')
    return p.parse_args()


def jpath(a, name):
    return os.path.join(a.work, name)


def read_header_scale(a):
    """-> (pixscale arcsec, header or None)"""
    hdr = None
    if a.image and os.path.isfile(a.image):
        _, hdr = imageio.load_image(a.image)
    if a.pixscale > 0:
        return a.pixscale, hdr
    if hdr is None:
        raise SystemExit('lensmodel: need --pixscale or an image with a WCS')
    from astropy.wcs import WCS
    from astropy.wcs.utils import proj_plane_pixel_scales
    w = WCS(hdr)
    if not w.has_celestial:
        raise SystemExit('lensmodel: the image has no celestial WCS; give --pixscale')
    return float(np.mean(proj_plane_pixel_scales(w.celestial)) * 3600.0), hdr


def sky_pa(hdr, ref_pix0, phi_deg):
    """position angle east of north (deg) of the direction phi (from +x, CCW) at pixel ref_pix0 = (x0, y0) 0-based; None without WCS"""
    try:
        from astropy.wcs import WCS
        w = WCS(hdr).celestial
        if not w.has_celestial:
            return None
        c = math.cos(math.radians(phi_deg))
        s = math.sin(math.radians(phi_deg))
        a = w.pixel_to_world(ref_pix0[0], ref_pix0[1])
        b = w.pixel_to_world(ref_pix0[0] + 20 * c, ref_pix0[1] + 20 * s)
        return float(a.position_angle(b).deg % 180.0)
    except Exception:
        return None


def parse_kv(text):
    out = {}
    for t in (text or '').split(','):
        if '=' in t:
            k, v = t.split('=', 1)
            out[k.strip()] = float(v)
    return out


def inputs(a, scale):
    """-> numbers, pix (N x 2, 0-based), lens pixel (0-based), catalog (cols, rows)"""
    cols, rows = ([], [])
    if a.catalog and os.path.isfile(a.catalog):
        cols, rows = tsvio.read_catalog(a.catalog)
    byn = {}
    for r in rows:
        try:
            byn[int(float(r.get('NUMBER', 'nan')))] = r
        except ValueError:
            pass
    nums, pix = [], []
    if a.image_numbers.strip():
        for n in tsvio.parse_numbers(a.image_numbers):
            if n not in byn:
                raise SystemExit('lensmodel: image NUMBER %d is not in the catalog' % n)
            nums.append(n)
            pix.append([tsvio.fnum(byn[n]['X_IMAGE']) - 1.0, tsvio.fnum(byn[n]['Y_IMAGE']) - 1.0])
    elif a.positions.strip():
        for i, t in enumerate(a.positions.split(';')):
            if t.strip():
                x, y = [float(v) for v in t.replace(' ', '').split(',')]
                nums.append(-(i + 1))
                pix.append([x - 1.0, y - 1.0])
    pix = np.array(pix, float).reshape(-1, 2)
    lens = None
    if a.lens_number and a.lens_number in byn:
        lens = [tsvio.fnum(byn[a.lens_number]['X_IMAGE']) - 1.0, tsvio.fnum(byn[a.lens_number]['Y_IMAGE']) - 1.0]
    elif a.lens_x > 0 and a.lens_y > 0:
        lens = [a.lens_x - 1.0, a.lens_y - 1.0]
    return nums, pix, lens, (cols, byn)


def perturbers_from(a, cols, byn, ref, scale, skip_numbers):
    out = []
    desc = []
    if a.nfw_kappa_s > 0:
        out.append(lm.NFW(a.nfw_kappa_s, a.nfw_theta_s, 0.0, 0.0))
        desc.append(dict(type='NFW', kappa_s=a.nfw_kappa_s, theta_s=a.nfw_theta_s))
    if a.member_column and a.member_column in cols:
        for n, r in byn.items():
            if n in skip_numbers or n == a.lens_number or not tsvio.fnum(r.get(a.member_column), 0) > 0:
                continue
            m = tsvio.fnum(r.get(a.member_mag_column))
            if not np.isfinite(m):
                continue
            ratio = 10 ** (-0.4 * (m - a.member_mag_star))
            sc = math.sqrt(ratio)
            x = (tsvio.fnum(r['X_IMAGE']) - 1.0 - ref[0]) * scale
            y = (tsvio.fnum(r['Y_IMAGE']) - 1.0 - ref[1]) * scale
            out.append(lm.PIEMD(a.member_b_star * sc, a.member_core * sc, a.member_cut * sc, x, y))
            desc.append(dict(type='PIEMD', number=n, b0=a.member_b_star * sc, a=a.member_core * sc, s=a.member_cut * sc, x=x, y=y))
    return out, desc


def to_pix(ref, scale, x, y):
    return ref[0] + np.asarray(x) / scale + 1.0, ref[1] + np.asarray(y) / scale + 1.0


def overlay_regions(path, ref, scale, curves, images, source=None):
    L = ['# Region file format: DS9 version 4.1', 'global font="helvetica 9 normal" width=1', 'image']
    for c in curves:
        x, y = to_pix(ref, scale, c['x'], c['y'])
        pts = ','.join('%.2f,%.2f' % (a, b) for a, b in zip(x[::max(1, len(x) // 180)], y[::max(1, len(y) // 180)]))
        L.append('polygon(%s) # color=red width=2 text={critical curve (%s)} tag={lensmodel}' % (pts, c['kind']))
        cx, cy = to_pix(ref, scale, c['cx'], c['cy'])
        pts = ','.join('%.2f,%.2f' % (a, b) for a, b in zip(cx[::max(1, len(cx) // 120)], cy[::max(1, len(cy) // 120)]))
        L.append('polygon(%s) # color=cyan dash=1 text={caustic (source plane, same pixel frame)} tag={lensmodel}' % pts)
    for im in images:
        col = 'green' if im['matched'] >= 0 else 'yellow'
        txt = ('img %d' % im['matched'] if im['matched'] >= 0 else 'counter-image') + ' mu=%+.1f' % im['mu']
        L.append('circle(%.2f,%.2f,%.1f) # color=%s width=2 text={%s} tag={lensmodel}' % (im['px'], im['py'], 8.0, col, txt))
    if source:
        L.append('point(%.2f,%.2f) # point=cross 12 color=magenta text={source} tag={lensmodel}' % source)
    open(path, 'w').write('\n'.join(L) + '\n')


def build_result_model(res):
    p = res['params']
    return lm.build_sie_shear(p)


def curves_payload(model, half):
    cc = model.critical_curves(half=half, n=801)
    return cc


def physical(res_params, theta_E_eff, z_l, z_s, h0, om):
    out = {}
    if z_l > 0 and z_s > z_l:
        Dl, Ds, Dls = lm.distances(z_l, z_s, h0, om)
        out.update(D_l_Mpc=float(Dl), D_s_Mpc=float(Ds), D_ls_Mpc=float(Dls), sigma_crit_Msun_Mpc2=float(lm.sigma_crit(z_l, z_s, h0, om)),
                   mass_Einstein_Msun=float(lm.mass_in_radius(theta_E_eff, z_l, z_s, h0, om)), theta_E_kpc=float(theta_E_eff * lm.ARCSEC * Dl * 1000.0),
                   sigma_SIS_kms=float(lm.sis_sigma_v(theta_E_eff, z_l, z_s, h0, om)))
    return out


def predictions(a, model, bx, by, obs_arc, ref, scale, z_l, z_s):
    half = max(4.0 * max(model.einstein_guess(), 0.5), 2 * float(np.max(np.hypot(obs_arc[:, 0], obs_arc[:, 1]))) if len(obs_arc) else 0.0)
    ims = lm.predict(model, bx, by, (0.0, 0.0), observed=obs_arc, match_tol=max(0.15, 3 * a.sigma_pos), half=half, n=a.grid)
    for im in ims:
        im['px'], im['py'] = [float(v) for v in to_pix(ref, scale, im['x'], im['y'])]
    if z_l > 0 and z_s > z_l and ims:
        dts = lm.time_delays(model, ims, bx, by, z_l, z_s, a.h0, a.omega_m)
        for im, d in zip(ims, dts):
            im['delay_days'] = float(d)
    return ims


def task_fit(a):
    scale, hdr = read_header_scale(a)
    nums, pix, lens, (cols, byn) = inputs(a, scale)
    if len(pix) < 2:
        raise SystemExit('lensmodel: need at least 2 multiple images (--image-numbers or --positions)')
    if lens is None:
        raise SystemExit('lensmodel: the lens centre is required (--lens-number or --lens-x/--lens-y)')
    ref = np.array(lens)
    obs = (pix - ref) * scale
    n = len(pix)
    free = ['theta_E', 'q', 'phi', 'gamma', 'phi_g'] if a.free == 'auto' and n >= 4 else (['theta_E', 'q', 'phi'] if a.free == 'auto' and n == 3 else
            (['theta_E'] if a.free == 'auto' else [t.strip() for t in a.free.split(',') if t.strip()]))
    fixed = parse_kv(a.fixed)
    if a.free == 'auto' and n < 4:
        fixed.setdefault('gamma', 0.0)
        if n == 2:
            fixed.setdefault('q', 1.0)
    perts, pdesc = perturbers_from(a, cols, byn, ref, scale, set(nums))
    fl = None
    if a.flux_column:
        if not byn or len(nums) != n or min(nums) < 0:
            raise SystemExit('lensmodel: --flux-column needs --image-numbers')
        fl = [tsvio.fnum(byn[k].get(a.flux_column)) for k in nums]
        if not all(np.isfinite(fl)) or min(fl) <= 0:
            raise SystemExit('lensmodel: invalid fluxes in %s' % a.flux_column)
    fit = lm.fit_sie_shear(obs, a.sigma_pos, (0.0, 0.0), perturbers=perts, free=free, fixed=fixed or None, fluxes=fl, flux_sigma=a.flux_sigma,
                           n_starts=a.n_starts, seed=a.seed)
    if not fit.get('success'):
        raise SystemExit('lensmodel: fit failed (%s)' % fit.get('message', ''))
    p = fit['params']
    p['phi'] = float(p['phi'] % 180.0)
    p['phi_g'] = float(p['phi_g'] % 180.0)
    model = lm.build_sie_shear(p, perts)
    bx, by = fit['source']
    cc = model.critical_curves(half=max(3 * p['theta_E'], 2 * float(np.max(np.hypot(obs[:, 0], obs[:, 1])))), n=801)
    theta_eff = cc[0]['theta_eff'] if cc else float('nan')
    ims = predictions(a, model, bx, by, obs, ref, scale, a.z_lens, a.z_source)
    counters = [im for im in ims if im['matched'] < 0]
    res = dict(task='fit', pixscale=scale, lens_pixel=[float(ref[0] + 1), float(ref[1] + 1)], image_numbers=nums, observed_pixel=(pix + 1).tolist(), observed_arcsec=obs.tolist(),
               sigma_pos=a.sigma_pos, free=fit['free'], fixed=fixed, params=p, errors=fit['errors'], chi2=fit['chi2'], dof=fit['dof'], n_data=fit['n_data'],
               rms_arcsec=fit['rms_arcsec'], residuals_arcsec=fit['residuals_arcsec'], source_arcsec=[bx, by],
               source_pixel=[float(v) for v in to_pix(ref, scale, bx, by)], theta_E_critical_curve=theta_eff, ellipticity=1 - p['q'],
               pa_sky_deg=sky_pa(hdr, ref, p['phi']) if hdr is not None else None, perturbers=pdesc,
               z_lens=a.z_lens, z_source=a.z_source, images=ims, n_images_predicted=len(ims), n_counter_images=len(counters),
               physical=physical(p, theta_eff, a.z_lens, a.z_source, a.h0, a.omega_m), warnings=[])
    if fit['dof'] <= 0:
        res['warnings'].append('dof = %d: the model has as many (or more) free parameters as data; chi2 = 0 is not evidence for the model' % fit['dof'])
    if len(ims) < n:
        res['warnings'].append('the model predicts %d images but %d were given' % (len(ims), n))
    if a.check_lenstronomy:
        try:
            ax_n, ay_n = lm.SIE(p['theta_E'], p['q'], p['phi'], 0.0, 0.0).alpha(obs[:, 0], obs[:, 1])
            ax_l, ay_l = lm.lenstronomy_alpha('SIE', p, obs[:, 0], obs[:, 1])
            res['lenstronomy_check'] = dict(max_abs_diff_arcsec=float(max(np.abs(ax_n - ax_l).max(), np.abs(ay_n - ay_l).max())))
        except ImportError:
            res['lenstronomy_check'] = dict(skipped='lenstronomy not installed')
    os.makedirs(a.work, exist_ok=True)
    cur = []
    for c in cc:
        cur.append(dict(kind=c['kind'], theta_eff=c['theta_eff'], area=c['area'], x=c['x'].tolist(), y=c['y'].tolist(), cx=c['cx'].tolist(), cy=c['cy'].tolist()))
    json.dump(dict(curves=cur), open(jpath(a, 'lens_curves.json'), 'w'))
    json.dump(res, open(jpath(a, 'lens_model.json'), 'w'), indent=1, default=float)
    write_images_tsv(a, res)
    overlay_regions(jpath(a, 'lens_overlay.reg'), ref, scale, cc, ims, tuple(to_pix(ref, scale, bx, by)))
    # add_columns output for the image rows
    rows = []
    for k, num in enumerate(nums):
        if num < 0:
            continue
        match = [im for im in ims if im['matched'] == k]
        d = dict(LENS_ROLE='image', LENS_RES=fit['residuals_arcsec'][k])
        if match:
            d.update(LENS_MU=match[0]['mu'], LENS_DT=match[0].get('delay_days'), LENS_PARITY=match[0]['parity'])
        rows.append((num, d))
    if a.meta_out:
        try:
            from ogfkit import meta as ometa
            ometa.update(a.meta_out, 'lensmodel', dict(theta_E=p['theta_E'], q=p['q'], phi=p['phi'], gamma=p['gamma'], phi_g=p['phi_g'], chi2=fit['chi2'], dof=fit['dof'],
                                                          rms_arcsec=fit['rms_arcsec'], n_counter_images=len(counters)), nrows=len(rows))
        except Exception:
            pass
    print_summary(res, sys.stderr)
    tsvio.write_columns(COLUMNS, rows)


def write_images_tsv(a, res):
    recs = []
    for k, im in enumerate(res['images']):
        recs.append(dict(IMAGE=k + 1, X_IMAGE=im['px'], Y_IMAGE=im['py'], X_ARCSEC=im['x'], Y_ARCSEC=im['y'], MU=im['mu'], PARITY=im['parity'],
                         MATCHED=(res['image_numbers'][im['matched']] if im['matched'] >= 0 else ''), DELAY_DAYS=im.get('delay_days', '')))
    tsvio.write_table(jpath(a, 'lens_images.tsv'), ['IMAGE', 'X_IMAGE', 'Y_IMAGE', 'X_ARCSEC', 'Y_ARCSEC', 'MU', 'PARITY', 'MATCHED', 'DELAY_DAYS'], recs)


def print_summary(res, fh):
    p = res['params']
    e = res['errors']
    L = []
    L.append('SIE+shear lens model: %d images, chi2 = %.3f for %d dof, image-plane rms = %.4f arcsec (sigma_pos %.4f)' % (len(res['image_numbers']), res['chi2'], res['dof'],
                                                                                                                    res['rms_arcsec'], res['sigma_pos']))
    for k in ('theta_E', 'q', 'phi', 'gamma', 'phi_g'):
        L.append('  %-8s = %9.4f%s' % (k, p[k], (' +- %.4f' % e[k]) if k in e else ' (fixed)'))
    L.append('  source   = (%.4f, %.4f) arcsec, pixel (%.2f, %.2f)' % (res['source_arcsec'][0], res['source_arcsec'][1], res['source_pixel'][0], res['source_pixel'][1]))
    L.append('  theta_E from the tangential critical curve = %.4f arcsec; sky PA (E of N) = %s' % (res['theta_E_critical_curve'], ('%.1f' % res['pa_sky_deg']) if res['pa_sky_deg'] is not None else 'n/a'))
    for k, v in res['physical'].items():
        L.append('  %s = %.4g' % (k, v))
    L.append('  predicted images: %d (counter-images not given: %d)' % (res['n_images_predicted'], res['n_counter_images']))
    for im in res['images']:
        L.append('    %-12s px (%8.2f, %8.2f)  mu = %+8.2f  %-6s%s' % ('image %d' % res['image_numbers'][im['matched']] if im['matched'] >= 0 else 'COUNTER-IMAGE',
                                                                     im['px'], im['py'], im['mu'], im['parity'], ('  delay %.2f d' % im['delay_days']) if 'delay_days' in im else ''))
    for w in res['warnings']:
        L.append('  WARNING: ' + w)
    fh.write('\n'.join(L) + '\n')


def load_result(a):
    f = jpath(a, 'lens_model.json')
    if not os.path.isfile(f):
        raise SystemExit('lensmodel: no lens_model.json in %s - run Fit first' % a.work)
    return json.load(open(f))


def model_from_result(res):
    p = res['params']
    perts = []
    for d in res.get('perturbers', []):
        if d['type'] == 'NFW':
            perts.append(lm.NFW(d['kappa_s'], d['theta_s'], 0.0, 0.0))
        elif d['type'] == 'PIEMD':
            perts.append(lm.PIEMD(d['b0'], d['a'], d['s'], d['x'], d['y']))
    return lm.build_sie_shear(p, perts)


def task_curves(a):
    res = load_result(a)
    model = model_from_result(res)
    scale = res['pixscale']
    ref = np.array(res['lens_pixel']) - 1.0
    half = max(3 * res['params']['theta_E'], 2 * float(np.max(np.hypot(*np.array(res['observed_arcsec']).T))))
    cc = model.critical_curves(half=half, n=801)
    cur = [dict(kind=c['kind'], theta_eff=c['theta_eff'], area=c['area'], x=c['x'].tolist(), y=c['y'].tolist(), cx=c['cx'].tolist(), cy=c['cy'].tolist()) for c in cc]
    json.dump(dict(curves=cur), open(jpath(a, 'lens_curves.json'), 'w'))
    overlay_regions(jpath(a, 'lens_overlay.reg'), ref, scale, cc, res['images'], tuple(res['source_pixel']))
    print('curves: %d critical curve(s); tangential theta_eff = %.4f arcsec' % (len(cc), cc[0]['theta_eff'] if cc else float('nan')))


def task_magmap(a):
    res = load_result(a)
    model = model_from_result(res)
    data, hdr = imageio.load_image(a.image)
    ny, nx = data.shape
    scale = res['pixscale']
    ref = np.array(res['lens_pixel']) - 1.0
    out = np.zeros((ny, nx), np.float32)
    xs = (np.arange(nx) - ref[0]) * scale
    for y0 in range(0, ny, 64):
        y1 = min(ny, y0 + 64)
        X, Y = np.meshgrid(xs, (np.arange(y0, y1) - ref[1]) * scale)
        mu = model.mu(X, Y)
        out[y0:y1] = np.clip(mu, -a.max_mu, a.max_mu)
    f = jpath(a, 'lens_magnification.fits')
    imageio.save_fits(f, out, hdr, extra=dict(LENSMAP='signed magnification 1/detA clipped at +-%g' % a.max_mu, LENSTHE=res['params']['theta_E']))
    print('magnification map: %s  (|mu|>10 on %.2f%% of pixels)' % (os.path.basename(f), 100.0 * float(np.mean(np.abs(out) > 10))))


def task_source(a):
    res = load_result(a)
    model = model_from_result(res)
    data, hdr = imageio.load_image(a.image)
    scale = res['pixscale']
    ref = np.array(res['lens_pixel']) - 1.0
    rad = a.cutout if a.cutout > 0 else 2.5 * res['params']['theta_E'] + 0.5
    R = int(rad / scale)
    x0, x1 = int(max(0, ref[0] - R)), int(min(data.shape[1], ref[0] + R + 1))
    y0, y1 = int(max(0, ref[1] - R)), int(min(data.shape[0], ref[1] + R + 1))
    cut = data[y0:y1, x0:x1].astype(float)
    sky = float(np.nanmedian(cut))
    cut = cut - sky
    origin = (ref[0] - x0, ref[1] - y0)
    if a.lens_mask > 0:
        yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
        cut[np.hypot(xx - origin[0], yy - origin[1]) * scale < a.lens_mask] = np.nan
    if a.source_method == 'inversion':
        return source_inversion(a, res, model, cut, scale, origin, hdr)
    sp = a.source_pix if a.source_pix > 0 else 0.5 * scale
    rec = lm.ray_trace_image(model, cut, scale, origin, (0.0, 0.0), src_pix=sp)
    fwd = lm.lens_source_map(model, rec, cut.shape, scale, origin, (0.0, 0.0))
    ok = np.isfinite(cut) & (fwd != 0)
    resid = np.where(ok, cut - fwd, 0.0)
    from astropy.io import fits
    h = fits.Header()
    h['CDELT1'] = rec['pix']; h['CDELT2'] = rec['pix']; h['CRPIX1'] = 1; h['CRPIX2'] = 1
    h['CRVAL1'] = rec['x0']; h['CRVAL2'] = rec['y0']; h['CUNIT1'] = 'arcsec'; h['CUNIT2'] = 'arcsec'
    h['COMMENT'] = 'source-plane reconstruction; coordinates in arcsec of the lens frame; surface brightness conserved'
    fits.PrimaryHDU(np.nan_to_num(rec['src'], nan=0.0).astype(np.float32), header=h).writeto(jpath(a, 'lens_source.fits'), overwrite=True)
    fits.PrimaryHDU(fwd.astype(np.float32)).writeto(jpath(a, 'lens_model_image.fits'), overwrite=True)
    fits.PrimaryHDU(resid.astype(np.float32)).writeto(jpath(a, 'lens_residual.fits'), overwrite=True)
    sig = float(np.nanstd(cut[~ok])) if (~ok).sum() > 50 else float(np.nanstd(cut))
    s = np.nan_to_num(rec['src'], nan=0.0)
    tot = float(s.sum())
    cx = float((s.sum(0) * (np.arange(s.shape[1]))).sum() / tot * rec['pix'] + rec['x0']) if tot != 0 else float('nan')
    cy = float((s.sum(1) * (np.arange(s.shape[0]))).sum() / tot * rec['pix'] + rec['y0']) if tot != 0 else float('nan')
    rms = float(np.sqrt(np.mean(resid[ok] ** 2))) if ok.any() else float('nan')
    out = dict(sky_subtracted=sky, source_pixel_arcsec=rec['pix'], source_shape=list(s.shape), source_centroid_arcsec=[cx, cy], n_pixels_mapped=int(np.sum(rec['count'])),
               residual_rms=rms, image_rms_in_mapped=float(np.sqrt(np.mean(cut[ok] ** 2))) if ok.any() else float('nan'), blank_pixel_rms=sig)
    res['source_reconstruction'] = out
    json.dump(res, open(jpath(a, 'lens_model.json'), 'w'), indent=1, default=float)
    print('source plane: %s  %dx%d px of %.4f arcsec, centroid (%.4f, %.4f) arcsec, residual rms %.4g vs image rms %.4g' % (
        os.path.basename(jpath(a, 'lens_source.fits')), s.shape[1], s.shape[0], rec['pix'], cx, cy, rms, out['image_rms_in_mapped']))


def source_inversion(a, res, model, cut, scale, origin, hdr):
    from ogfkit import lensextra as lx
    from astropy.io import fits
    ok = np.isfinite(cut)
    sig = a.noise
    if sig <= 0:
        yy, xx = np.mgrid[0:cut.shape[0], 0:cut.shape[1]]
        blank = ok & (np.hypot(xx - origin[0], yy - origin[1]) * scale > 0.9 * 2.5 * res['params']['theta_E'])
        sig = float(np.nanstd(cut[blank])) if blank.sum() > 100 else float(np.nanstd(cut[ok]))
    work = np.nan_to_num(cut, nan=0.0)
    mask = lx.arc_mask(work, sig, origin, 1e9) & ok                       # NaN pixels (--lens-mask) are excluded
    kern = lx.gaussian_kernel(a.psf_sigma) if a.psf_sigma > 0 else None
    S = lx.SourceInversion(model, work, sig, mask, scale, origin, (0.0, 0.0), kernel=kern, n=a.source_n, reg=a.regularisation)
    r = S.best()
    src, mod = S.images(r)
    resid = np.where(mask, work - mod, 0.0)
    g = S.grid
    h = fits.Header()
    h['CDELT1'] = g['pix']; h['CDELT2'] = g['pix']; h['CRPIX1'] = 1; h['CRPIX2'] = 1
    h['CRVAL1'] = g['x0']; h['CRVAL2'] = g['y0']; h['CUNIT1'] = 'arcsec'; h['CUNIT2'] = 'arcsec'
    h['COMMENT'] = 'pixelated source inversion (%s regularisation, lambda %.4g); coordinates in arcsec of the lens frame' % (a.regularisation, r['lam'])
    fits.PrimaryHDU(src.astype(np.float32), header=h).writeto(jpath(a, 'lens_source.fits'), overwrite=True)
    fits.PrimaryHDU(mod.astype(np.float32)).writeto(jpath(a, 'lens_model_image.fits'), overwrite=True)
    fits.PrimaryHDU(resid.astype(np.float32)).writeto(jpath(a, 'lens_residual.fits'), overwrite=True)
    neff = S.n_eff(r)
    out = dict(method='inversion', regularisation=a.regularisation, lam=float(r['lam']), noise=sig, psf_sigma_pix=a.psf_sigma, source_pixel_arcsec=g['pix'], source_shape=[g['n'], g['n']],
               n_data=int(S.ndata), n_eff=neff, chi2=r['chi2'], chi2_red=r['chi2'] / max(S.ndata - neff, 1), log_evidence=r['evidence'],
               source_flux_sum=float(src.sum() * g['pix'] ** 2), residual_rms=float(np.sqrt(np.mean(resid[mask] ** 2))))
    res['source_reconstruction'] = out
    json.dump(res, open(jpath(a, 'lens_model.json'), 'w'), indent=1, default=float)
    print('source inversion: %dx%d px of %.4f arcsec, lambda %.3g, chi2/dof %.3f (N_eff %.0f of %d pixels), residual rms %.4g (noise %.4g)' % (
        g['n'], g['n'], g['pix'], r['lam'], out['chi2_red'], neff, S.ndata, out['residual_rms'], sig))


def task_multiplane(a):
    from ogfkit import lensextra as lx
    scale, hdr = read_header_scale(a)
    if not a.mp_groups or not (a.lens_x or a.lens_y or a.lens_number):
        raise SystemExit('lensmodel: multiplane needs --mp-groups and the lens centre (--lens-x/--lens-y)')
    ref = np.array([a.lens_x - 1.0, a.lens_y - 1.0])
    groups = []
    zs = [t.strip() for t in a.mp_z.split(',')] if a.mp_z else []
    for k, g in enumerate(a.mp_groups.split('|')):
        pts = np.array([[float(v) for v in t.split(',')] for t in g.split(';') if t.strip()]) - 1.0
        z = zs[k] if k < len(zs) else 'free'
        groups.append(dict(xy=(pts - ref) * scale, z=(None if z in ('free', '', '-') else float(z))))
    if groups[0]['z'] is None:
        raise SystemExit('lensmodel: the first set defines the reference plane and needs a redshift')
    r = lx.fit_multiplane(groups, (0.0, 0.0), sigma=a.sigma_pos, z_l=a.z_lens or None, H0=a.h0, Om=a.omega_m, n_starts=a.n_starts)
    out = dict(params={k: float(v) for k, v in r['params'].items()}, weights=[float(w) for w in r['weights']], errors=r['errors'], chi2=r['chi2'], dof=r['dof'],
               rms_arcsec=r['rms_arcsec'], sources_arcsec=[[float(x), float(y)] for x, y in r['sources']], pixscale=scale, success=r['success'])
    json.dump(out, open(jpath(a, 'lens_multiplane.json'), 'w'), indent=1)
    p = out['params']
    print('multiplane: theta_E(ref) %.4f  q %.3f  phi %.1f  gamma %.4f  weights %s  chi2 %.2f / dof %d' % (p['theta_E'], p['q'], p['phi'], p['gamma'], ','.join('%.4f' % w for w in out['weights']), r['chi2'], r['dof']))


def task_predict(a):
    res = load_result(a)
    model = model_from_result(res)
    scale = res['pixscale']
    ref = np.array(res['lens_pixel']) - 1.0
    if np.isfinite(a.source_x) and np.isfinite(a.source_y):
        bx, by = (a.source_x - 1.0 - ref[0]) * scale, (a.source_y - 1.0 - ref[1]) * scale
    else:
        bx, by = res['source_arcsec']
    zl = a.z_lens or res.get('z_lens', 0.0)
    zs = a.z_source or res.get('z_source', 0.0)
    obs = np.array(res['observed_arcsec'])
    ims = predictions(a, model, bx, by, obs, ref, scale, zl, zs)
    res['images'] = ims
    res['source_arcsec'] = [bx, by]
    res['source_pixel'] = [float(v) for v in to_pix(ref, scale, bx, by)]
    res['n_images_predicted'] = len(ims)
    res['n_counter_images'] = len([i for i in ims if i['matched'] < 0])
    res['z_lens'], res['z_source'] = zl, zs
    res['physical'] = physical(res['params'], res['theta_E_critical_curve'], zl, zs, a.h0, a.omega_m)
    json.dump(res, open(jpath(a, 'lens_model.json'), 'w'), indent=1, default=float)
    write_images_tsv(a, res)
    cur = json.load(open(jpath(a, 'lens_curves.json')))['curves'] if os.path.isfile(jpath(a, 'lens_curves.json')) else []
    overlay_regions(jpath(a, 'lens_overlay.reg'), ref, scale, [dict(kind=c['kind'], x=np.array(c['x']), y=np.array(c['y']), cx=np.array(c['cx']), cy=np.array(c['cy'])) for c in cur],
                    ims, tuple(res['source_pixel']))
    print_summary(res, sys.stdout)


def task_simulate(a):
    """a synthetic SIE+shear lens: image FITS (lens galaxy + lensed Sersic arcs + point-source images), catalog and truth json"""
    rng = np.random.default_rng(a.seed)
    lens = lenssynth.random_lens(rng, quad=bool(a.sim_quad), theta_E=a.sim_theta_e)
    p = lens['params']
    scale = a.pixscale if a.pixscale > 0 else 0.05
    n = 301
    ref = ((n - 1) / 2.0, (n - 1) / 2.0)
    model = lm.build_sie_shear(p)
    bx, by = lens['source']
    arcs = lenssynth.render_arcs(model, (n, n), scale, ref, (0.0, 0.0), dict(x0=bx + 0.05, y0=by - 0.03, reff=0.12, n=1.0, q=0.7, phi=30.0, amp=0.5), psf_sigma_pix=1.5, noise=0.0)
    yy, xx = np.mgrid[0:n, 0:n]
    img = arcs.copy()
    r2 = ((xx - ref[0]) * scale) ** 2 + ((yy - ref[1]) * scale) ** 2
    if a.sim_lens_light:
        img += 3.0 * np.exp(-7.669 * ((np.sqrt(r2) / 0.6) ** 0.25 - 1))      # lens galaxy (de Vaucouleurs-like)
    flux0 = 100.0
    from scipy.ndimage import gaussian_filter
    pts = np.zeros_like(img)
    cat = [dict(NUMBER=1, X_IMAGE=ref[0] + 1, Y_IMAGE=ref[1] + 1, FLUX=0.0, MAG_AUTO=19.0, CL_MEMBER=0)]
    for k, im in enumerate(lens['images']):
        px, py = to_pix(np.array(ref), scale, im['x'], im['y'])
        f = flux0 * abs(im['mu'])
        pts[int(round(py - 1)), int(round(px - 1))] += f
        cat.append(dict(NUMBER=k + 2, X_IMAGE=float(px), Y_IMAGE=float(py), FLUX=f, MAG_AUTO=25 - 2.5 * math.log10(f), CL_MEMBER=0))
    img += gaussian_filter(pts, 1.5) * 1.0
    img += rng.normal(0, a.sim_noise, img.shape)
    os.makedirs(a.work, exist_ok=True)
    from astropy.io import fits
    h = fits.Header()
    h['CTYPE1'] = 'RA---TAN'; h['CTYPE2'] = 'DEC--TAN'; h['CRPIX1'] = ref[0] + 1; h['CRPIX2'] = ref[1] + 1; h['CRVAL1'] = 150.0; h['CRVAL2'] = 2.0
    h['CD1_1'] = -scale / 3600.0; h['CD2_2'] = scale / 3600.0; h['CD1_2'] = 0.0; h['CD2_1'] = 0.0
    fits.PrimaryHDU(img.astype(np.float32), header=h).writeto(jpath(a, 'sim_lens.fits'), overwrite=True)
    tsvio.write_table(jpath(a, 'sim_catalog.tsv'), ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'FLUX', 'MAG_AUTO', 'CL_MEMBER'], cat)
    truth = dict(params=p, source_arcsec=[float(bx), float(by)], pixscale=scale, lens_pixel=[ref[0] + 1, ref[1] + 1],
                 images=[dict(x=im['x'], y=im['y'], mu=im['mu'], parity=im['parity']) for im in lens['images']],
                 image_numbers=list(range(2, 2 + len(lens['images']))), source_sersic=dict(reff=0.12, n=1.0))
    json.dump(truth, open(jpath(a, 'sim_truth.json'), 'w'), indent=1, default=float)
    print('simulated lens in %s: theta_E %.3f, q %.3f, phi %.1f, gamma %.3f, %d images' % (a.work, p['theta_E'], p['q'], p['phi'], p['gamma'], len(lens['images'])))


def main():
    a = pa_args()
    os.makedirs(a.work, exist_ok=True)
    {'fit': task_fit, 'curves': task_curves, 'magmap': task_magmap, 'source': task_source, 'predict': task_predict, 'simulate': task_simulate, 'multiplane': task_multiplane}[a.task](a)


if __name__ == '__main__':
    main()
