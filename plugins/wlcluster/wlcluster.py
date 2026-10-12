#!/usr/bin/env python3
"""Cluster weak lensing for the current catalogue.

    wlcluster.py --task shapes --catalog TSV --image FITS --work DIR
    wlcluster.py --task mass   --catalog TSV --work DIR --z-lens Z
    wlcluster.py --task mcal   --work DIR

shapes measures a shear on each stamp (metacalibration) and writes WL_E1, WL_E2,
WL_G1, WL_G2, WL_RESP, WL_FLAG and WL_W. mass converts a g1, g2 catalogue in the
tangent plane into a tangential profile and an NFW mass. Both follow the
add_columns contract. mass writes DIR/wl_profile.tsv and DIR/wl_mass.json.
Positions in the catalogue are 1-based. A sky catalogue uses
x = -(RA-RA_c) cos(Dec_c), y = Dec-Dec_c, in arcsec, so east is to the left.
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

from ogfkit import tsvio, meta as ometa  # noqa: E402
from ogfkit import wlcluster as wl  # noqa: E402

SHAPE_COLUMNS = ['WL_E1', 'WL_E2', 'WL_G1', 'WL_G2', 'WL_RESP', 'WL_FLAG', 'WL_W']
MASS_COLUMNS = ['WL_R_MPC', 'WL_GT', 'WL_GX', 'WL_W', 'WL_IN']
PROFILE_COLUMNS = ['R_MPC', 'R_ARCSEC', 'N', 'GT', 'GT_ERR', 'GX', 'GX_ERR',
                   'DELTASIGMA', 'DELTASIGMA_ERR', 'BOOST', 'GT_MODEL', 'DELTASIGMA_MODEL']


def col(rows, name, default=float('nan')):
    return np.array([tsvio.fnum(r.get(name), default) for r in rows])


def numbers(rows):
    out = []
    for i, r in enumerate(rows):
        n = str(r.get('NUMBER', '')).strip()
        out.append(n if n else str(i + 1))
    return out


def require_column(cols, name):
    if name not in cols:
        raise SystemExit('wlcluster: column was not found: %s' % name)


def plane_arcsec(a, rows, cols):
    """Arcsec relative to the cluster centre. x increases to the right, y up."""
    if a.coord == 'sky':
        require_column(cols, a.ra_column)
        require_column(cols, a.dec_column)
        ra = col(rows, a.ra_column)
        dec = col(rows, a.dec_column)
        x = -(ra - a.ra) * np.cos(np.deg2rad(a.dec)) * 3600.0
        y = (dec - a.dec) * 3600.0
        return x, y
    if a.center_x <= 0 or a.center_y <= 0:
        raise SystemExit('wlcluster: coord xy needs center-x and center-y, 1-based pixels, both positive')
    if a.pixel_scale <= 0:
        raise SystemExit('wlcluster: coord xy needs a positive pixel-scale in arcsec per pixel')
    require_column(cols, a.x_column)
    require_column(cols, a.y_column)
    # Catalogue pixels and the centre are both 1-based, so the difference is the offset.
    x = (col(rows, a.x_column) - a.center_x) * a.pixel_scale
    y = (col(rows, a.y_column) - a.center_y) * a.pixel_scale
    return x, y


def load_nz(path):
    if not path:
        return None
    if not os.path.isfile(path):
        raise SystemExit('wlcluster: n(z) file was not found: %s' % path)
    z, n = [], []
    with open(path, encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#'):
                continue
            parts = s.replace(',', ' ').split()
            if len(parts) < 2:
                continue
            z.append(float(parts[0]))
            n.append(float(parts[1]))
    if len(z) < 2:
        raise SystemExit('wlcluster: n(z) file has fewer than two rows: %s' % path)
    return np.asarray(z, float), np.asarray(n, float)


def source_redshifts(a, cols, rows, cosmo):
    n = len(rows)
    if a.z_source > 0:
        return np.full(n, float(a.z_source))
    if a.z_column and a.z_column in cols:
        z = col(rows, a.z_column)
        if np.isfinite(z).any():
            return z
    nz = load_nz(a.nz)
    if nz is not None:
        zeff = wl.effective_source_z(cosmo, a.z_lens, nz[0], nz[1], a.z_gap)
        return np.full(n, zeff)
    raise SystemExit('wlcluster: a source redshift was not found (z-source, a z column, or an n(z) file)')


def weights(a, cols, rows):
    if a.w_column and a.w_column in cols:
        w = col(rows, a.w_column)
        return np.where(np.isfinite(w), w, 0.0)
    return np.ones(len(rows))


def random_radii(a, Dl):
    if not a.random_catalog:
        return None
    if not os.path.isfile(a.random_catalog):
        raise SystemExit('wlcluster: random catalogue was not found: %s' % a.random_catalog)
    rcols, rrows = tsvio.read_catalog(a.random_catalog)
    if not rrows:
        return np.zeros(0)
    x, y = plane_arcsec(a, rrows, rcols)
    R = np.hypot(x, y) * wl.ARCSEC * Dl
    return R[np.isfinite(R)]


def write_meta(path, value, nrows):
    if not path:
        return
    folder = os.path.dirname(os.path.abspath(path))
    if folder:
        os.makedirs(folder, exist_ok=True)
    prev = ometa.load(path).get('wlcluster')
    merged = dict(prev) if isinstance(prev, dict) else {}
    merged.update(value)
    ometa.update(path, 'wlcluster', merged, nrows=nrows)


def task_shapes(a, cols, rows):
    if not a.image or not os.path.isfile(a.image):
        raise SystemExit('wlcluster: image was not found: %s' % (a.image or ''))
    require_column(cols, a.x_column)
    require_column(cols, a.y_column)
    from ogfkit import imageio
    data, _hdr = imageio.load_image(a.image)
    x = col(rows, a.x_column) - 1.0
    y = col(rows, a.y_column) - 1.0
    measured = wl.measure_stamps(data, x, y, psf_fwhm_pix=a.psf_fwhm, stamp=a.stamp, dg=a.dg)
    selection = wl.selection_response(measured, a.dg, min_resolution=a.min_resolution, emax=a.emax)
    out = []
    n_ok = 0
    for num, rec in zip(numbers(rows), measured):
        r1, r2 = rec['R1'], rec['R2']
        resp = 0.5 * (r1 + r2) if np.isfinite(r1) and np.isfinite(r2) else None
        ok = int(rec['flag']) == 0
        n_ok += int(ok)
        out.append((num, dict(
            WL_E1=rec['e1'], WL_E2=rec['e2'], WL_G1=rec['g1'], WL_G2=rec['g2'],
            WL_RESP=resp, WL_FLAG=int(rec['flag']), WL_W=(1.0 if ok else None))))
    meta = dict(task='shapes', n_shapes=len(rows), n_shapes_ok=n_ok, psf_fwhm_pix=a.psf_fwhm,
                min_resolution=a.min_resolution, emax=a.emax)
    meta.update(selection)
    write_meta(a.meta_out, meta, len(rows))
    return out


def _profile_records(rows):
    recs = []
    for r in rows:
        recs.append(dict(
            R_MPC=r['R'], R_ARCSEC=r.get('R_arcsec'), N=r['n'],
            GT=r['gt'], GT_ERR=r['gt_err'], GX=r['gx'], GX_ERR=r.get('gx_err', r['gt_err']),
            DELTASIGMA=r['deltasigma'], DELTASIGMA_ERR=r['deltasigma_err'], BOOST=r['boost'],
            GT_MODEL=r.get('gt_model'), DELTASIGMA_MODEL=r.get('deltasigma_model')))
    return recs


def task_mass(a, cols, rows):
    if a.z_lens <= 0:
        raise SystemExit('wlcluster: z-lens must be positive')
    require_column(cols, a.g1_column)
    require_column(cols, a.g2_column)
    if a.rmin <= 0 or a.rmax <= a.rmin:
        raise SystemExit('wlcluster: rmin and rmax must satisfy 0 < rmin < rmax')
    cosmo = wl.Cosmology(H0=a.h0, Om=a.omega_m, ns=a.ns, sigma8=a.sigma8)
    x, y = plane_arcsec(a, rows, cols)
    z = source_redshifts(a, cols, rows, cosmo)
    w = weights(a, cols, rows)
    g1 = col(rows, a.g1_column)
    g2 = col(rows, a.g2_column)
    Dl = float(cosmo.angular_diameter(a.z_lens))
    rand = random_radii(a, Dl)
    g1c, g2c = wl.correct_shear(g1, g2, a.m_bias, a.c1, a.c2)
    g1c, g2c, fac1, fac2 = wl.apply_selection_response(
        g1c, g2c, a.mean_response, a.selection_r1, a.selection_r2)
    gt, gx = wl.tangential(x, y, g1c, g2c)
    R = np.hypot(x, y) * wl.ARCSEC * Dl
    # m and c were already removed. select_sources must not remove them again.
    sel = wl.select_sources(
        x, y, g1c, g2c, z, w, a.z_lens, cosmo, z_gap=a.z_gap, z_shift=a.z_shift,
        rmin=a.rmin, rmax=a.rmax, gmax=a.gmax, m=0.0, c1=0.0, c2=0.0)
    edges = wl.bin_edges(a.rmin, a.rmax, a.nbins)
    summary, prof = wl.fit_mass(
        sel, edges, cosmo, a.z_lens, sigma_e=a.sigma_e, concentration=a.concentration,
        c_fixed=a.c_fixed, two_halo=a.two_halo, reduced=a.reduced, miscenter=a.miscenter,
        f_cen=a.f_cen, sigma_off=a.sigma_off, random_R=rand, samples=a.samples, seed=a.seed,
        min_per_bin=a.min_per_bin, power=a.power, lss_cov=a.lss_cov, n_jackknife=a.n_jackknife,
        count_slope=a.count_slope, fit_miscenter=a.fit_miscenter,
        f_cen_sigma=a.f_cen_sigma, sigma_off_sigma=a.sigma_off_sigma)
    summary['selection_factor_g1'] = fac1
    summary['selection_factor_g2'] = fac2
    summary['boost'] = 'not supplied' if rand is None else 'normalised to the outer radial bin'
    summary['coord'] = a.coord
    summary['task'] = 'mass'
    os.makedirs(a.work, exist_ok=True)
    tsvio.write_table(os.path.join(a.work, 'wl_profile.tsv'), PROFILE_COLUMNS, _profile_records(prof))
    with open(os.path.join(a.work, 'wl_mass.json'), 'w', encoding='utf-8') as f:
        json.dump(summary, f, indent=1)
        f.write('\n')
    write_meta(a.meta_out, {k: v for k, v in summary.items() if k != 'limits'}, len(rows))
    used = np.zeros(len(rows), int)
    if len(sel['index']):
        used[sel['index']] = 1
    out = []
    for i, num in enumerate(numbers(rows)):
        out.append((num, dict(
            WL_R_MPC=R[i] if np.isfinite(R[i]) else None,
            WL_GT=gt[i] if np.isfinite(gt[i]) else None,
            WL_GX=gx[i] if np.isfinite(gx[i]) else None,
            WL_W=w[i] if np.isfinite(w[i]) and w[i] > 0 else None,
            WL_IN=int(used[i]))))
    return out


def task_mcal(a):
    """Internal m and c on galaxies drawn in this program. No catalogue is read."""
    if a.stamp < 16:
        raise SystemExit('wlcluster: stamp must be at least 16 pixels')
    sigma = 0.0
    if a.psf_fwhm > 0:
        sigma = a.psf_fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    doc = wl.calibrate_shear(kind=a.galaxy, n_sersic=a.sersic_n, size=a.stamp,
                             psf_sigma=sigma, noise=a.mcal_noise, seed=a.seed, dg=a.dg)
    doc['task'] = 'mcal'
    os.makedirs(a.work, exist_ok=True)
    with open(os.path.join(a.work, 'wl_mcal.json'), 'w', encoding='utf-8') as f:
        json.dump(doc, f, indent=1)
        f.write('\n')
    write_meta(a.meta_out, doc, 0)
    if doc.get('success'):
        print('m %.6f c1 %.6f c2 %.6f' % (doc['m'], doc['c1'], doc['c2']))
    else:
        print('mcal failed')
    return doc


def build_parser():
    ap = argparse.ArgumentParser(description='Cluster weak lensing')
    ap.add_argument('--task', required=True, choices=['shapes', 'mass', 'mcal'])
    ap.add_argument('--catalog', default='')
    ap.add_argument('--work', default='')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--image', default='')
    ap.add_argument('--x-column', default='X_IMAGE')
    ap.add_argument('--y-column', default='Y_IMAGE')
    ap.add_argument('--psf-fwhm', type=float, default=0.0)
    ap.add_argument('--stamp', type=int, default=32)
    ap.add_argument('--dg', type=float, default=0.01)
    ap.add_argument('--min-resolution', type=float, default=0.0)
    ap.add_argument('--emax', type=float, default=0.0)
    ap.add_argument('--coord', default='xy', choices=['xy', 'sky'])
    ap.add_argument('--center-x', type=float, default=0.0)
    ap.add_argument('--center-y', type=float, default=0.0)
    ap.add_argument('--pixel-scale', type=float, default=0.0)
    ap.add_argument('--ra-column', default='RA')
    ap.add_argument('--dec-column', default='DEC')
    ap.add_argument('--ra', type=float, default=0.0)
    ap.add_argument('--dec', type=float, default=0.0)
    ap.add_argument('--g1-column', default='WL_G1')
    ap.add_argument('--g2-column', default='WL_G2')
    ap.add_argument('--w-column', default='WL_W')
    ap.add_argument('--z-column', default='Z')
    ap.add_argument('--z-lens', type=float, default=0.0)
    ap.add_argument('--z-source', type=float, default=0.0)
    ap.add_argument('--z-shift', type=float, default=0.0)
    ap.add_argument('--z-gap', type=float, default=0.1)
    ap.add_argument('--nz', default='')
    ap.add_argument('--random-catalog', default='')
    ap.add_argument('--rmin', type=float, default=0.3)
    ap.add_argument('--rmax', type=float, default=3.0)
    ap.add_argument('--nbins', type=int, default=8)
    ap.add_argument('--gmax', type=float, default=1.5)
    ap.add_argument('--sigma-e', type=float, default=0.26)
    ap.add_argument('--power', default='halofit', choices=['halofit', 'linear'])
    ap.add_argument('--lss-cov', dest='lss_cov', action='store_true')
    ap.add_argument('--no-lss-cov', dest='lss_cov', action='store_false')
    ap.add_argument('--n-jackknife', type=int, default=0)
    ap.add_argument('--count-slope', type=float, default=0.0)
    ap.add_argument('--mean-response', type=float, default=0.0)
    ap.add_argument('--selection-r1', type=float, default=0.0)
    ap.add_argument('--selection-r2', type=float, default=0.0)
    ap.add_argument('--m-bias', type=float, default=0.0)
    ap.add_argument('--c1', type=float, default=0.0)
    ap.add_argument('--c2', type=float, default=0.0)
    ap.add_argument('--h0', type=float, default=70.0)
    ap.add_argument('--omega-m', type=float, default=0.3)
    ap.add_argument('--ns', type=float, default=0.96)
    ap.add_argument('--sigma8', type=float, default=0.8)
    ap.add_argument('--concentration', default='duffy', choices=['duffy', 'free', 'fixed'])
    ap.add_argument('--c-fixed', type=float, default=4.0)
    ap.add_argument('--reduced', dest='reduced', action='store_true')
    ap.add_argument('--no-reduced', dest='reduced', action='store_false')
    ap.add_argument('--two-halo', dest='two_halo', action='store_true')
    ap.add_argument('--no-two-halo', dest='two_halo', action='store_false')
    ap.add_argument('--miscenter', action='store_true')
    ap.add_argument('--fit-miscenter', action='store_true')
    ap.add_argument('--f-cen', type=float, default=0.8)
    ap.add_argument('--sigma-off', type=float, default=0.1)
    ap.add_argument('--f-cen-sigma', type=float, default=0.0)
    ap.add_argument('--sigma-off-sigma', type=float, default=0.0)
    ap.add_argument('--galaxy', default='gaussian', choices=['gaussian', 'sersic'])
    ap.add_argument('--sersic-n', type=float, default=1.0)
    ap.add_argument('--mcal-noise', type=float, default=0.0)
    ap.add_argument('--samples', type=int, default=200)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--min-per-bin', type=int, default=8)
    ap.set_defaults(reduced=True, two_halo=True, lss_cov=True)
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.task == 'mcal':
        if not a.work:
            raise SystemExit('wlcluster: work directory was not found')
        task_mcal(a)
        return
    if not a.catalog or not os.path.isfile(a.catalog):
        raise SystemExit('wlcluster: catalog was not found: %s' % (a.catalog or ''))
    if not a.work:
        raise SystemExit('wlcluster: work directory was not found')
    cols, rows = tsvio.read_catalog(a.catalog)
    if a.task == 'shapes':
        out = task_shapes(a, cols, rows)
        columns = SHAPE_COLUMNS
    else:
        out = task_mass(a, cols, rows)
        columns = MASS_COLUMNS
    tsvio.write_columns(columns, out)


if __name__ == '__main__':
    main()
