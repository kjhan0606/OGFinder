#!/usr/bin/env python3
"""Spectra linked to catalog objects: 1D spectra, 2D spectra, IFU cubes; line fits and redshifts.

    spectra.py --task link --catalog TSV --work DIR [--spec-dir D --spec-pattern spec_{NUMBER}.fits | --link-file TSV]
    spectra.py --task fit  --catalog TSV --work DIR [same link options] [--snr-min 4 --z-column Z_SPEC --fit-line Ha ...]
    spectra.py --task plot --number N --catalog TSV --work DIR --out PNG [same link options]
    spectra.py --task kin  --catalog TSV --work DIR [same link options] [--kin-line Ha --kin-z Z --kin-inc DEG ...]   (2D slit / IFU-cube kinematics, see docs/spectra.md)
    spectra.py --task science --catalog TSV --work DIR [same link options] [--wave-frame vacuum --ebv-mw 0 --inst-fwhm 0 --h0 70 --omega-m 0.3]

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

from ogfkit import tsvio, spectra as sp, meta as ometa, kinematics as kin, specscience as sci  # noqa: E402

COLUMNS = {
    'link': ['SP_FILE', 'SP_KIND', 'SP_OK'],
    'fit': ['SP_Z', 'SP_ZERR', 'SP_ZQ', 'SP_NLINES', 'SP_SNR', 'SP_LINE_FLUX', 'SP_LINE_FLUXERR', 'SP_LINE_FWHM_KMS', 'SP_LINE_EW'],
    'kin': ['SP_KIN_VSYS', 'SP_KIN_VSINI', 'SP_KIN_VC', 'SP_KIN_VC_ERR', 'SP_KIN_RT', 'SP_KIN_PA', 'SP_KIN_INC', 'SP_KIN_SIGMA', 'SP_KIN_CHI2R', 'SP_KIN_N'],
    'science': ['SP_SCI_Z', 'SP_SCI_ZERR', 'SP_SCI_ZQ', 'SP_SCI_ZSRC', 'SP_SIGMA', 'SP_SIGMA_ERR', 'SP_SIGMA_KIND', 'SP_TYPE', 'SP_SFR', 'SP_SFR_EBV', 'SP_HA_FLUX'],
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


def json_clean(o):
    """numpy -> python, NaN/inf -> null (the Tcl JSON parser and strict readers reject NaN)"""
    if isinstance(o, dict):
        return {str(k): json_clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [json_clean(v) for v in o]
    if isinstance(o, np.ndarray):
        return json_clean(o.tolist())
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (float, np.floating)):
        return float(o) if np.isfinite(o) else None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return o


def kin_redshift(a, r, lk):
    """redshift for the kinematics: --kin-z, else the z column, else the blind redshift of the extracted 1D spectrum (quality >= 1)"""
    if a.kin_z > 0:
        return a.kin_z, 'kin-z'
    zk = tsvio.fnum(r.get(a.z_column)) if a.z_column and a.z_column in r else float('nan')
    if np.isfinite(zk) and zk > 0:
        return zk, 'catalog column %s' % a.z_column
    s = load_1d(a, lk)
    err = s['err']
    if err is not None:
        err = np.where(np.isfinite(err) & (err > 0), err, np.nan)
    an = sp.analyse_spectrum(s['wave'], s['flux'], err, a.snr_min, a.kernel_sigma_px, a.cont_width, z_known=None, zmax=a.zmax)
    z = an['redshift']
    if z['quality'] >= 1:
        return float(z['z']), 'blind redshift (quality %d)' % z['quality']
    return None, 'no redshift'


def kin_plot_2d(path, t, lam_sys, prof, rc, number):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.8), gridspec_kw=dict(width_ratios=[1.3, 1]))
    d = t['data']
    w = t['wave']
    m = np.abs(w - lam_sys) < 25 * lam_sys / 6563.0
    v1, v2 = np.nanpercentile(d[:, m], [5, 99.7])
    ax[0].imshow(d[:, m], origin='lower', aspect='auto', cmap='gray', vmin=v1, vmax=v2, extent=[w[m][0], w[m][-1], -0.5, d.shape[0] - 0.5])
    if len(prof['y']):
        ax[0].plot(lam_sys * (1 + prof['v'] / kin.C_KMS), prof['y'], 'r.', ms=4)
    ax[0].set_xlabel('wavelength (A)'); ax[0].set_ylabel('row')
    ax[1].errorbar(prof['y'], prof['v'], prof['v_err'], fmt='ko', ms=3)
    if rc.get('ok'):
        yy = np.linspace(prof['y'].min(), prof['y'].max(), 200)
        ax[1].plot(yy, rc['vsys'] + kin.arctan_rc(yy - rc['y0'], rc['vc_obs'], rc['rt']), 'r-')
    ax[1].set_xlabel('row'); ax[1].set_ylabel('v (km/s)')
    ax[1].set_title('object %s: arctan fit%s' % (number, (' v_obs=%.0f r_t=%.2f px' % (abs(rc['vc_obs']), rc['rt'])) if rc.get('ok') else ' failed'), fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)


def kin_plot_cube(path, maps, vf, number):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2, 3, figsize=(9, 6))
    vmax = np.nanpercentile(np.abs(maps['vel']), 98) if np.isfinite(maps['vel']).any() else 100
    panels = [('flux', maps['flux'], 'gray', None), ('velocity (km/s)', maps['vel'], 'RdBu_r', (-vmax, vmax)), ('dispersion (km/s)', maps['sigma'], 'viridis', None)]
    if vf.get('ok'):
        panels += [('model velocity', np.where(np.isfinite(maps['vel']), vf['model'], np.nan), 'RdBu_r', (-vmax, vmax)), ('residual', maps['vel'] - vf['model'], 'PuOr', (-0.5 * vmax, 0.5 * vmax)),
                   ('S/N', maps['snr'], 'magma', None)]
    for ax, (ti, im, cm, lim) in zip(axs.ravel(), panels):
        h = ax.imshow(im, origin='lower', cmap=cm, vmin=lim[0] if lim else None, vmax=lim[1] if lim else None)
        ax.set_title(ti, fontsize=8)
        plt.colorbar(h, ax=ax, fraction=0.046)
    if vf.get('ok'):
        for ax in axs.ravel()[:5]:
            th = math.radians(vf['pa'])
            ax.plot([vf['x0'] - 8 * math.cos(th), vf['x0'] + 8 * math.cos(th)], [vf['y0'] - 8 * math.sin(th), vf['y0'] + 8 * math.sin(th)], 'k--', lw=0.6)
        fig.suptitle('object %s: PA %.1f  inc %.1f  vc %.0f  rt %.2f px  chi2r %.2f' % (number, vf['pa'], vf['inc'], vf['vc'], vf['rt'], vf['chi2r']), fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)


def save_maps(path, maps, vf, hdr_extra):
    from astropy.io import fits
    hl = fits.HDUList([fits.PrimaryHDU()])
    ex = dict(hdr_extra)
    for name, key in (('FLUX', 'flux'), ('VEL', 'vel'), ('VELERR', 'vel_err'), ('SIGMA', 'sigma'), ('SIGERR', 'sigma_err'), ('SNR', 'snr')):
        h = fits.ImageHDU(np.asarray(maps[key], np.float32), name=name)
        for k, v in ex.items():
            h.header[k] = v
        hl.append(h)
    if vf.get('ok'):
        hl.append(fits.ImageHDU(np.asarray(vf['model'], np.float32), name='MODEL'))
        hl.append(fits.ImageHDU(np.asarray(maps['vel'] - vf['model'], np.float32), name='RESID'))
    hl.append(fits.ImageHDU(np.asarray(maps['bin'], np.int32), name='BIN'))
    hl.writeto(path, overwrite=True)
    base = path[:-len('_maps.fits')] if path.endswith('_maps.fits') else path[:-5]
    for suffix, key in (('vel', 'vel'), ('sigma', 'sigma'), ('flux', 'flux')):       # single-HDU files for the ds9 frames
        h = fits.PrimaryHDU(np.asarray(maps[key], np.float32))
        for k, v in ex.items():
            h.header[k] = v
        h.writeto('%s_%s.fits' % (base, suffix), overwrite=True)


def task_kin(a, cols, rows, W):
    links = build_links(a, rows)
    out, results = [], {}
    nfit = 0
    wu = None if a.wave_unit == 'auto' else a.wave_unit
    if a.kin_line not in sp.LINE_DICT:
        raise SystemExit('spectra: unknown line %r (known: %s)' % (a.kin_line, ', '.join(sp.LINE_DICT)))
    rest = sp.LINE_DICT[a.kin_line]
    for r, lk in zip(rows, links):
        d = {}
        num = lk['NUMBER']
        if lk['OK'] and lk['KIND'] in ('2d', 'cube'):
            try:
                z, how = kin_redshift(a, r, lk)
                if z is None:
                    results[num] = dict(error='no redshift (give kin-z, a redshift column or a spectrum with a blind redshift)')
                    out.append((r['NUMBER'], d))
                    continue
                lam_sys = rest * (1 + z)
                inc = a.kin_inc if a.kin_inc > 0 else None
                if lk['KIND'] == '2d':
                    t = sp.read_2d(lk['FILE'], wu)
                    prof = kin.fit_slit_line(t['wave'], t['data'], t['err'], lam_sys, a.kin_window_kms, a.kin_bin_snr, a.kin_snr_pix, a.kin_inst_fwhm)
                    rc = kin.fit_rotation_curve(prof['y'], prof['v'], prof['v_err'], inc_deg=inc, slit_offset_deg=a.kin_slit_psi, free_centre=True)
                    sdisp = kin.dispersion_summary(prof['sigma'], prof['sigma_err'], prof['flux'])
                    if rc.get('ok'):
                        d.update(SP_KIN_VSYS=rc['vsys'], SP_KIN_VSINI=rc['vflat_obs'], SP_KIN_RT=rc['rt'], SP_KIN_CHI2R=rc['chi2r'], SP_KIN_N=len(prof['y']), SP_KIN_SIGMA=sdisp['mean'] if sdisp['n'] else None)
                        if inc:
                            d.update(SP_KIN_VC=abs(rc['vc']) if rc.get('vc') is not None else None, SP_KIN_VC_ERR=rc.get('vc_err'), SP_KIN_INC=inc)
                        nfit += 1
                    kin_plot_2d(os.path.join(W, 'kin_%s.png' % num), t, lam_sys, prof, rc, num)
                    tsvio.write_table(os.path.join(W, 'kin_%s.tsv' % num), ['ROW', 'V', 'V_ERR', 'SIGMA', 'SIGMA_ERR', 'FLUX', 'FLUX_ERR', 'SNR', 'NROWS'],
                                      [dict(ROW=prof['y'][i], V=prof['v'][i], V_ERR=prof['v_err'][i], SIGMA=prof['sigma'][i], SIGMA_ERR=prof['sigma_err'][i], FLUX=prof['flux'][i],
                                            FLUX_ERR=prof['flux_err'][i], SNR=prof['snr'][i], NROWS=prof['nrows'][i]) for i in range(len(prof['y']))])
                    results[num] = dict(kind='2d', z=z, z_source=how, line=a.kin_line, lam_sys=lam_sys, rotation_curve={k: v for k, v in rc.items()}, dispersion=sdisp, n_bins=int(len(prof['y'])))
                else:
                    c = sp.read_cube(lk['FILE'])
                    maps = kin.fit_cube_line(c['data'], c['wave'], lam_sys, c['err'], a.kin_window_kms, a.kin_snr_pix, a.kin_bin_snr, a.kin_inst_fwhm)
                    vf = kin.fit_velocity_field(maps['vel'], maps['vel_err'], inc_fixed=inc, free_centre=True)
                    sdisp = kin.dispersion_summary(maps['sigma'], maps['sigma_err'], maps['flux'])
                    if vf.get('ok'):
                        d.update(SP_KIN_VSYS=vf['vsys'], SP_KIN_VSINI=vf['vsini'], SP_KIN_VC=vf['vc'], SP_KIN_VC_ERR=vf['errors'].get('vc'), SP_KIN_RT=vf['rt'], SP_KIN_PA=vf['pa'], SP_KIN_INC=vf['inc'],
                                 SP_KIN_SIGMA=sdisp['mean'] if sdisp['n'] else None, SP_KIN_CHI2R=vf['chi2r'], SP_KIN_N=vf['n'])
                        nfit += 1
                    save_maps(os.path.join(W, 'kin_%s_maps.fits' % num), maps, vf, dict(LINE=a.kin_line, LAMSYS=lam_sys, BUNIT='km/s (VEL, SIGMA)'))
                    kin_plot_cube(os.path.join(W, 'kin_%s.png' % num), maps, vf, num)
                    results[num] = dict(kind='cube', z=z, z_source=how, line=a.kin_line, lam_sys=lam_sys, velocity_field={k: v for k, v in vf.items() if k != 'model'}, dispersion=sdisp,
                                        n_fits=int(maps['n_fits']), n_spaxels=int(maps['n_spaxels']))
            except Exception as ex:
                results[num] = dict(error=str(ex))
                sys.stderr.write('spectra kin: object %s: %s\n' % (num, ex))
        out.append((r['NUMBER'], d))
    return out, dict(n_fit=nfit, n_2d_or_cube=sum(1 for l in links if l['OK'] and l['KIND'] in ('2d', 'cube'))), links, results


def task_science(a, cols, rows, W):
    links = build_links(a, rows)
    out = []
    results = {}
    n_ok = 0
    for r, lk in zip(rows, links):
        d = {}
        if lk['OK']:
            try:
                s = load_1d(a, lk)
                err = s['err']
                if err is not None:
                    err = np.where(np.isfinite(err) & (err > 0), err, np.nan)
                zk = tsvio.fnum(r.get(a.z_column)) if a.z_column and a.z_column in r else float('nan')
                an = sci.galaxy_analysis(s['wave'], s['flux'], err, frame=a.wave_frame, ebv_mw=a.ebv_mw, inst_fwhm_a=a.inst_fwhm,
                                         snr_min=a.snr_min, kernel_sigma_px=a.kernel_sigma_px, cont_width=a.cont_width,
                                         z_known=zk if np.isfinite(zk) and zk > 0 else None, zmax=a.zmax, h0=a.h0, om0=a.omega_m)
                z = an.get('redshift') or {}
                disp = an.get('dispersion') or {}
                sfr = an.get('sfr') or {}
                balmer = an.get('balmer') or {}
                d.update(SP_SCI_Z=z.get('z'), SP_SCI_ZERR=z.get('z_err'), SP_SCI_ZQ=z.get('quality'), SP_SCI_ZSRC=z.get('source'),
                         SP_SIGMA=disp.get('sigma'), SP_SIGMA_ERR=disp.get('sigma_err'), SP_SIGMA_KIND=disp.get('kind'),
                         SP_TYPE=an.get('type'), SP_SFR=sfr.get('sfr') if sfr.get('applies') else None,
                         SP_SFR_EBV=balmer.get('ebv') if balmer.get('used') else None, SP_HA_FLUX=an.get('ha_flux'))
                n_ok += 1
                results[lk['NUMBER']] = an
            except Exception as ex:
                results[lk['NUMBER']] = dict(error=str(ex))
                sys.stderr.write('spectra science: object %s: %s\n' % (lk['NUMBER'], ex))
        out.append((r['NUMBER'], d))
    return out, dict(n_science=n_ok, n_linked=sum(l['OK'] for l in links)), links, results


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
    ap.add_argument('--task', required=True, choices=['link', 'fit', 'plot', 'kin', 'science'])
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
    ap.add_argument('--wave-frame', default='vacuum', choices=['vacuum', 'air'])
    ap.add_argument('--ebv-mw', type=float, default=0.0)
    ap.add_argument('--inst-fwhm', type=float, default=0.0)
    ap.add_argument('--h0', type=float, default=70.0)
    ap.add_argument('--omega-m', type=float, default=0.3)
    ap.add_argument('--kin-line', default='Ha')
    ap.add_argument('--kin-z', type=float, default=0.0)
    ap.add_argument('--kin-window-kms', type=float, default=600.0)
    ap.add_argument('--kin-snr-pix', type=float, default=3.0)
    ap.add_argument('--kin-bin-snr', type=float, default=0.0)
    ap.add_argument('--kin-inst-fwhm', type=float, default=0.0)
    ap.add_argument('--kin-inc', type=float, default=0.0)
    ap.add_argument('--kin-slit-psi', type=float, default=0.0)
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
    elif a.task == 'kin':
        out, summ, links, results = task_kin(a, cols, rows, a.work)
        with open(os.path.join(a.work, 'spectra_kin.json'), 'w') as fh:
            json.dump(json_clean(results), fh)
    elif a.task == 'science':
        out, summ, links, results = task_science(a, cols, rows, a.work)
        with open(os.path.join(a.work, 'spectra_science.json'), 'w') as fh:
            json.dump(json_clean(results), fh)
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
