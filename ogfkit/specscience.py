"""Paper measurements on one 1D spectrum.

Preprocess (sort, drop bad pixels, air to vacuum, Milky Way extinction), then the
existing emission-line engine plus an absorption-line redshift.  Gaussian velocity
dispersion, a spectral type, and an H-alpha star-formation rate.

Wavelengths are vacuum Angstrom after prepare_spectrum.  Flux stays in the file's
units.  This is not a stellar-population fit.  The dispersion is a Gaussian line
width, not a template line-of-sight velocity distribution.
"""
import math

import numpy as np

from .spectra import FWHM, LINE_DICT, LINES, analyse_spectrum, continuum, detect_lines, estimate_redshift, fit_line, refine_redshift

C_KMS = 299792.458
SFR_HA = 5.37e-42
HA_HB_CASEB = 2.86
BROAD_FWHM_KMS = 1200.0
BPT_SNR = 3.0
SIGMA_SNR = 5.0
RV_CCM = 3.1
RV_CALZETTI = 4.05

# Laboratory air wavelengths (Angstrom).  Vacuum values are derived with air_to_vacuum.
ABS_AIR = [
    ('CaK', 3933.663, 1.0),
    ('CaH', 3968.469, 0.8),
    ('MgI5167', 5167.321, 0.6),
    ('MgI5173', 5172.684, 0.6),
    ('MgI5184', 5183.604, 0.6),
    ('NaD5890', 5889.951, 0.7),
    ('NaD5896', 5895.924, 0.7),
    ('CaT8498', 8498.02, 0.5),
    ('CaT8542', 8542.09, 0.5),
    ('CaT8662', 8662.14, 0.5),
]


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def air_to_vacuum(wave_air):
    """Air to vacuum, Angstrom.  SDSS / Morton formula as used by VALD.

    s = 1e4/lambda_air (um^-1);
    n = 1 + 0.00008336624212083 + 0.02408926869968/(130.1065924522 - s^2)
        + 0.0001599740894897/(38.92568793293 - s^2);
    lambda_vac = lambda_air * n.
    """
    w = np.asarray(wave_air, float)
    s2 = (1e4 / w) ** 2
    n = (1.0 + 0.00008336624212083
         + 0.02408926869968 / (130.1065924522 - s2)
         + 0.0001599740894897 / (38.92568793293 - s2))
    out = w * n
    if np.ndim(wave_air) == 0:
        return float(out)
    return out


def vacuum_to_air(wave_vac, niter=6):
    """Inverse of air_to_vacuum.  Iterates lambda_air = lambda_vac / n(lambda_air)."""
    w = np.asarray(wave_vac, float)
    air = w.copy()
    for _ in range(niter):
        air = w / (np.asarray(air_to_vacuum(air), float) / air)
    if np.ndim(wave_vac) == 0:
        return float(air)
    return air


ABS_LINES = [(n, float(air_to_vacuum(lam)), wt) for n, lam, wt in ABS_AIR]


def ccm_extinction(wave_angstrom, ebv, rv=RV_CCM):
    """A(lambda) in magnitudes.  Cardelli, Clayton & Mathis 1989.

    Infrared 0.3 < x < 1.1 and optical 1.1 <= x < 3.3, with x = 1/lambda_um.
    A(lambda)/A(V) = a + b/Rv, and A(lambda) = that * Rv * E(B-V).
    Wavelengths outside that range are NaN.
    """
    w = np.atleast_1d(np.asarray(wave_angstrom, float))
    x = 1e4 / w
    a = np.full(w.shape, np.nan)
    b = np.full(w.shape, np.nan)
    ir = (x > 0.3) & (x < 1.1)
    op = (x >= 1.1) & (x < 3.3)
    a[ir] = 0.574 * x[ir] ** 1.61
    b[ir] = -0.527 * x[ir] ** 1.61
    y = x[op] - 1.82
    a[op] = (1 + 0.17699 * y - 0.50447 * y ** 2 - 0.02427 * y ** 3
             + 0.72085 * y ** 4 + 0.01979 * y ** 5 - 0.77530 * y ** 6 + 0.32999 * y ** 7)
    b[op] = (1.41338 * y + 2.28305 * y ** 2 + 1.07233 * y ** 3
             - 5.38434 * y ** 4 - 0.62251 * y ** 5 + 5.30260 * y ** 6 - 2.09002 * y ** 7)
    out = (a + b / rv) * rv * float(ebv)
    if np.ndim(wave_angstrom) == 0:
        return float(out[0])
    return out


def calzetti_k(wave_angstrom):
    """Calzetti (2000) k(lambda), Rv = 4.05.  A(lambda) = E(B-V) * k(lambda).

    The polynomial is the same one as sed_fit.grid.generate.calzetti_extinction:
    that function returns Av * k / Rv, so k = calzetti_extinction(wave, Av=Rv).
    """
    lam_um = np.asarray(wave_angstrom, float) / 1e4
    k = np.where(
        lam_um < 0.63,
        2.659 * (-2.156 + 1.509 / lam_um - 0.198 / lam_um ** 2 + 0.011 / lam_um ** 3) + RV_CALZETTI,
        2.659 * (-1.857 + 1.040 / lam_um) + RV_CALZETTI,
    )
    return np.maximum(k, 0.0)


def calzetti_a(wave_angstrom, ebv):
    """A(lambda) in magnitudes for a nebular E(B-V) on the Calzetti 2000 curve."""
    return float(ebv) * calzetti_k(wave_angstrom)


def prepare_spectrum(wave, flux, err=None, frame='vacuum', ebv_mw=0.0):
    """Sort, drop non-finite samples, convert air to vacuum, optionally deredden the Milky Way.

    ebv_mw = 0 leaves the flux unchanged.  Dereddening multiplies flux and err by 10^(0.4 A).
    Pixels outside the Cardelli range are left unchanged and counted.
    """
    w = np.asarray(wave, float)
    f = np.asarray(flux, float)
    e = None if err is None else np.asarray(err, float)
    m = np.isfinite(w) & np.isfinite(f)
    if e is not None:
        m = m & np.isfinite(e) & (e > 0)
    n_dropped = int((~m).sum())
    w, f = w[m], f[m]
    if e is not None:
        e = e[m]
    order = np.argsort(w, kind='mergesort')
    w, f = w[order], f[order]
    if e is not None:
        e = e[order]
    frame = (frame or 'vacuum').strip().lower()
    if frame == 'air':
        w = np.asarray(air_to_vacuum(w), float)
    elif frame != 'vacuum':
        raise ValueError('frame must be vacuum or air')
    n_no_extinction = 0
    ebv = float(ebv_mw or 0.0)
    if ebv > 0 and len(w):
        A = np.atleast_1d(np.asarray(ccm_extinction(w, ebv), float))
        ok = np.isfinite(A)
        factor = np.ones(len(w))
        factor[ok] = 10.0 ** (0.4 * A[ok])
        f = f * factor
        if e is not None:
            e = e * factor
        n_no_extinction = int((~ok).sum())
    return dict(wave=w, flux=f, err=e, n_dropped=n_dropped, n_no_extinction=n_no_extinction,
                frame='vacuum', input_frame=frame, ebv_mw=ebv)


def balmer_ebv(f_ha, f_hb, rest_ha=6564.61, rest_hb=4862.68):
    """Nebular E(B-V) from the Balmer decrement on the Calzetti 2000 curve.

    Case B intrinsic Ha/Hb = 2.86 at 10^4 K.
    A(Hb) - A(Ha) = 2.5 log10((Fha/Fhb) / 2.86).
    E(B-V) = that / (k_Hb - k_Ha).  The 0.44 stellar/gas factor is not applied.
    A ratio below 2.86 gives E(B-V) = 0 and used = False.
    """
    if f_hb is None or f_ha is None or f_hb <= 0 or f_ha <= 0:
        return dict(ebv=None, used=False, reason='non-positive flux', curve='Calzetti2000', intrinsic=HA_HB_CASEB, rv=RV_CALZETTI)
    ratio = float(f_ha) / float(f_hb)
    if ratio < HA_HB_CASEB:
        return dict(ebv=0.0, used=False, reason='below case B', ratio=ratio, curve='Calzetti2000', intrinsic=HA_HB_CASEB, rv=RV_CALZETTI)
    k_ha = float(np.asarray(calzetti_k(rest_ha)))
    k_hb = float(np.asarray(calzetti_k(rest_hb)))
    diff = 2.5 * math.log10(ratio / HA_HB_CASEB)
    ebv = diff / (k_hb - k_ha)
    return dict(ebv=float(ebv), used=True, ratio=ratio, a_ha=float(ebv * k_ha), a_hb=float(ebv * k_hb),
                curve='Calzetti2000', intrinsic=HA_HB_CASEB, rv=RV_CALZETTI)


def luminosity_distance_cm(z, h0=70.0, om0=0.3):
    """Luminosity distance in cm for FlatLambdaCDM(H0, Om0)."""
    from astropy.cosmology import FlatLambdaCDM
    cos = FlatLambdaCDM(H0=h0, Om0=om0)
    return float(cos.luminosity_distance(z).to('cm').value)


def halpha_sfr(flux_ha, z, h0=70.0, om0=0.3, flux_err=None):
    """H-alpha star-formation rate.

    Murphy et al. 2011 as adopted by Kennicutt & Evans 2012, Kroupa IMF:
    SFR (Msun/yr) = 5.37e-42 * L_Ha (erg/s).
    L = 4 pi D_L^2 F, with no extra (1+z).  z <= 0 or a non-positive flux gives SFR None.
    The uncertainty scales only the H-alpha flux error.  The Balmer-ratio uncertainty is not included.
    flux_unit is 'file': the number is physical only when the file flux is erg/s/cm^2.
    """
    out = dict(sfr=None, sfr_err=None, l_ha=None, calibration='Murphy2011/KennicuttEvans2012',
               constant=SFR_HA, imf='Kroupa', h0=float(h0), om0=float(om0), cosmology='FlatLambdaCDM',
               flux_unit='file', error_includes_balmer=False)
    z = _num(z)
    flux_ha = _num(flux_ha)
    if z is None or z <= 0 or flux_ha is None or flux_ha <= 0:
        return out
    dl = luminosity_distance_cm(z, h0, om0)
    l_ha = 4.0 * math.pi * dl * dl * flux_ha
    sfr = SFR_HA * l_ha
    out.update(sfr=float(sfr), l_ha=float(l_ha), d_l_cm=float(dl))
    fe = _num(flux_err)
    if fe is not None:
        out['sfr_err'] = float(abs(sfr * fe / flux_ha))
    return out


def kewley2001(log_nii_ha):
    """Kewley et al. 2001 maximum-starburst line: log([OIII]/Hb) = 0.61/(log([NII]/Ha) - 0.47) + 1.19."""
    return 0.61 / (log_nii_ha - 0.47) + 1.19


def kauffmann2003(log_nii_ha):
    """Kauffmann et al. 2003 empirical line: log([OIII]/Hb) = 0.61/(log([NII]/Ha) - 0.05) + 1.3."""
    return 0.61 / (log_nii_ha - 0.05) + 1.3


def classify_bpt(log_nii_ha, log_oiii_hb, log_sii_ha=None, log_oi_ha=None):
    """star-forming, composite, AGN, and (when [SII] or [OI] is given) Seyfert or LINER.

    Below Kauffmann: star-forming.  Between Kauffmann and Kewley: composite.  Above Kewley: AGN.
    An AGN with both [SII] lines is split by Kewley et al. 2006 (MNRAS 372, 961):
    Seyfert if log([OIII]/Hb) > 1.89 log([SII]/Ha) + 0.76, else LINER.
    The [OI] line of the same paper, 1.18 log([OI]/Ha) + 1.30, is used only when [OI] is measured
    and [SII] is not.
    """
    x = float(log_nii_ha)
    y = float(log_oiii_hb)
    if x < 0.05 and y < kauffmann2003(x):
        kind = 'star-forming'
    elif x < 0.47 and y < kewley2001(x):
        kind = 'composite'
    else:
        kind = 'AGN'
        if log_sii_ha is not None:
            kind = 'Seyfert' if y > 1.89 * float(log_sii_ha) + 0.76 else 'LINER'
        elif log_oi_ha is not None:
            kind = 'Seyfert' if y > 1.18 * float(log_oi_ha) + 1.30 else 'LINER'
    return dict(class_=kind, log_nii_ha=x, log_oiii_hb=y,
                log_sii_ha=None if log_sii_ha is None else float(log_sii_ha),
                log_oi_ha=None if log_oi_ha is None else float(log_oi_ha),
                kewley2001=None if x >= 0.47 else float(kewley2001(x)),
                kauffmann2003=None if x >= 0.05 else float(kauffmann2003(x)))


def line_sigma_kms(ft, inst_fwhm_a):
    """Gaussian sigma in km/s.  Instrumental FWHM is removed in quadrature when inst_fwhm_a > 0.

    sigma_kms = sigma_A / mu * c.
    sigma_inst = (FWHM_inst_A / 2.35482) / mu * c.
    When the line is narrower than the instrument, sigma_int = 0 and resolved is False.
    err_int = err_obs * sigma_obs / sigma_int.
    inst_fwhm_a == 0 returns the observed sigma and corrected = False.
    """
    if not ft.get('ok'):
        return None
    mu = _num(ft.get('mu'))
    sig = _num(ft.get('sigma'))
    sig_err = _num(ft.get('sigma_err'))
    if mu is None or sig is None or mu <= 0 or sig <= 0:
        return None
    sig_kms = sig / mu * C_KMS
    err_kms = (sig_err / mu * C_KMS) if sig_err is not None else None
    inst = _num(inst_fwhm_a) or 0.0
    if inst <= 0:
        return dict(sigma=sig_kms, sigma_err=err_kms, sigma_obs=sig_kms, sigma_err_obs=err_kms,
                    corrected=False, resolved=True, sigma_inst=None)
    sig_inst = (inst / FWHM) / mu * C_KMS
    if sig_kms > sig_inst:
        sint = math.sqrt(sig_kms ** 2 - sig_inst ** 2)
        err = (err_kms * sig_kms / sint) if (err_kms is not None and sint > 0) else None
        return dict(sigma=sint, sigma_err=err, sigma_obs=sig_kms, sigma_err_obs=err_kms,
                    corrected=True, resolved=True, sigma_inst=sig_inst)
    return dict(sigma=0.0, sigma_err=None, sigma_obs=sig_kms, sigma_err_obs=err_kms,
                corrected=True, resolved=False, sigma_inst=sig_inst)


def _weighted_mean(rows):
    use = [r for r in rows if r and r.get('resolved') and r.get('sigma_err') and r['sigma_err'] > 0 and _num(r.get('sigma')) is not None]
    if not use:
        return None, None, 0
    w = np.array([1.0 / r['sigma_err'] ** 2 for r in use])
    s = np.array([r['sigma'] for r in use])
    return float(np.sum(w * s) / np.sum(w)), float(1.0 / math.sqrt(np.sum(w))), len(use)


def _fit_growing(wave, flux, err, mu, dl, sigma0):
    """fit_line, widening the window when the Gaussian sigma sits on the upper bound.

    Narrow lines are unchanged: the first fit is returned.  A line broader than the
    default kernel (a broad-line AGN) is refit up to three times.
    """
    s0 = float(sigma0)
    ft = fit_line(wave, flux, err, mu, sigma0=s0)
    for _ in range(3):
        if not ft.get('ok'):
            return ft
        win = 6.0 * s0 + 4.0 * dl
        if ft['sigma'] < 0.85 * (win / 2.0):
            return ft
        s0 = max(ft['sigma'] * 2.5, s0 * 2.0)
        nxt = fit_line(wave, flux, err, mu, sigma0=s0)
        if not nxt.get('ok'):
            return ft
        ft = nxt
    return ft


def _snr(ft):
    if not ft or not ft.get('ok'):
        return 0.0
    fe = ft.get('flux_err') or 0.0
    if fe <= 0:
        return 0.0
    return ft['flux'] / fe


def _usable(ft, z, dl, snr_min=2.5):
    if not ft.get('ok') or not ft.get('mu_err') or ft['mu_err'] <= 0:
        return False
    if abs(_snr(ft)) <= snr_min:
        return False
    if z is None or not math.isfinite(z):
        return False
    return abs(ft['mu'] - ft['rest'] * (1.0 + z)) <= 6.0 * dl


def _slim(ft, inst_fwhm_a):
    if not ft.get('ok'):
        return dict(name=ft.get('name'), rest=_num(ft.get('rest')), ok=False)
    sig = line_sigma_kms(ft, inst_fwhm_a)
    fwhm_obs = ft['fwhm'] / ft['mu'] * C_KMS
    return dict(name=ft.get('name'), rest=_num(ft.get('rest')), ok=True, mu=ft['mu'], mu_err=ft['mu_err'],
                sigma_A=ft['sigma'], sigma_err_A=ft['sigma_err'], flux=ft['flux'], flux_err=ft['flux_err'],
                snr=_snr(ft), fwhm_obs_kms=fwhm_obs, ew=ft.get('ew'), chi2r=ft.get('chi2r'),
                dispersion=sig)


def _index(fits):
    out = {}
    for ft in fits:
        name = ft.get('name')
        if not name or not ft.get('ok'):
            continue
        prev = out.get(name)
        if prev is None or abs(_snr(ft)) > abs(_snr(prev)):
            out[name] = ft
    return out


def _positive(ft, cut):
    return bool(ft and ft.get('ok') and ft['flux'] > 0 and _snr(ft) > cut)


def _fit_matched(wave, flux_for_fit, err, matched, dl, sigma0, z):
    fits = []
    for name, lam, mu_obs, _s in matched:
        ft = _fit_growing(wave, flux_for_fit, err, mu_obs, dl, sigma0)
        ft['name'] = name
        ft['rest'] = lam
        if _usable(ft, z, dl):
            fits.append(ft)
    return fits


UV_CORE = ('CIV', 'CIII]', 'MgII')
UV_LINES = [row for row in LINES if row[0] in ('Lya', 'CIV', 'HeII', 'CIII]', 'MgII')]
NEBULAR_LINES = ('Hb', '[OIII]4960', '[OIII]5008', '[NII]6550', 'Ha', '[NII]6585', '[SII]6718', '[SII]6733')
BROAD_LINES = ('Ha', 'Hb', 'CIV', 'CIII]', 'MgII')
AGREE_KMS = 1500.0


def _fit_q(n):
    if n <= 0:
        return 0
    return 1 if n == 1 else (2 if n == 2 else 3)


def _close(z1, z2, tol_kms=AGREE_KMS):
    z1, z2 = _num(z1), _num(z2)
    if z1 is None or z2 is None:
        return False
    return abs(z1 - z2) / (1.0 + z1) * C_KMS < tol_kms


def choose_redshift(cands):
    """Pick one redshift from emission, ultraviolet-line, and absorption solutions.

    Each candidate needs z and a quality equal to the number of lines that actually fitted
    (1, 2, or 3+).  Solutions within 1500 km/s are one group.  The group that contains the
    highest quality is kept, and inside that group an emission solution of quality >= 2
    is kept ahead of an ultraviolet solution, then absorption.  A quality-1 emission
    redshift therefore does not replace a quality-2 or 3 absorption redshift at a
    different redshift.
    """
    cands = [c for c in cands if c and c.get('quality', 0) > 0 and _num(c.get('z')) is not None]
    if not cands:
        return dict(z=None, z_err=None, quality=0, source='none')
    anchor = max(cands, key=lambda c: (c['quality'], c.get('n') or 0))
    group = [c for c in cands if _close(c['z'], anchor['z'])]
    for source in ('emission', 'uv', 'absorption'):
        hits = [c for c in group if c['source'] == source and c['quality'] >= 2]
        if hits:
            win = max(hits, key=lambda c: (c['quality'], c.get('n') or 0))
            return dict(z=float(win['z']), z_err=_num(win.get('z_err')), quality=int(win['quality']), source=win['source'])
    win = anchor
    return dict(z=float(win['z']), z_err=_num(win.get('z_err')), quality=int(win['quality']), source=win['source'])


def _fit_at_redshift(wave, flux_for_fit, err, z, dl, sigma0, names):
    matched = []
    for name in names:
        lam = LINE_DICT[name]
        mu = lam * (1.0 + z)
        if wave[0] < mu < wave[-1]:
            matched.append((name, lam, mu, None))
    return _fit_matched(wave, flux_for_fit, err, matched, dl, sigma0, z)


def _uv_redshift(wave, flux, err, zmax, snr_min):
    """Redshift from CIV, C III] and Mg II after a wide continuum.

    The 101-pixel continuum follows a quasar broad line, so this pass uses a 401-pixel
    continuum and an 8-pixel kernel.  One detection is not enough, and Lyman-alpha does
    not count: two of CIV, C III] and Mg II must fit at one redshift.  The optical
    doublets are not in this list, so one broad bump cannot be counted as H-alpha plus
    both [N II] lines.
    """
    if len(wave) < 20:
        return None
    cont, _mask = continuum(wave, flux, 401)
    dl = float(np.median(np.diff(wave)))
    kernel = 8.0
    det = detect_lines(wave, flux, err, cont, kernel_sigma_px=kernel, snr_min=snr_min)
    tol_px = max(3.0, kernel * 1.5)
    zr = estimate_redshift(det, (float(wave[0]), float(wave[-1])), tol_px=tol_px, dl=dl, zmax=zmax, lines=UV_LINES)
    zgrid = _num(zr.get('z'))
    if zgrid is None:
        return None
    med = float(np.nanmedian(cont))
    ff = flux - cont + med
    tol = tol_px * dl
    core = []
    for name, lam, mu_obs, _s in zr.get('matched') or []:
        if name not in UV_CORE:
            continue
        ft = _fit_growing(wave, ff, err, mu_obs, dl, kernel * dl)
        ft['name'] = name
        ft['rest'] = lam
        if not ft.get('ok') or _snr(ft) <= 3.0 or ft['flux'] <= 0:
            continue
        if abs(ft['mu'] - lam * (1.0 + zgrid)) > max(tol, 6.0 * dl):
            continue
        core.append(ft)
    if len(core) < 2:
        return None
    z, ze, _c2 = refine_redshift([(ft['rest'], ft['mu'], max(ft['mu_err'], 1e-6)) for ft in core])
    keep = [ft for ft in core if _close(ft['mu'] / ft['rest'] - 1.0, z, 2000.0)]
    if len(keep) < 2:
        return None
    return dict(z=z, z_err=ze, quality=_fit_q(len(keep)), n=len(keep), source='uv', lines=keep)


def _absorption(wave, flux, err, cont, dl, kernel_sigma_px, snr_min, z_known, zmax):
    det = detect_lines(wave, flux, err, cont, kernel_sigma_px, snr_min, sign=-1.0)
    if z_known is None:
        zr = estimate_redshift(det, (float(wave[0]), float(wave[-1])), tol_px=max(3.0, kernel_sigma_px * 1.5),
                               dl=dl, zmax=zmax, lines=ABS_LINES)
    else:
        zr = dict(z=float(z_known), quality=4, n_match=0, matched=[], score=0.0, z_err=None)
        zr['matched'] = [(n, lam, lam * (1.0 + z_known), float('nan')) for n, lam, _w in ABS_LINES if wave[0] < lam * (1.0 + z_known) < wave[-1]]
        zr['n_match'] = len(zr['matched'])
    return zr, det


def _empty(reason):
    return dict(ok=False, reason=reason, type='unknown',
                redshift=dict(z=None, z_err=None, quality=0, source='none'),
                dispersion=dict(sigma=None, sigma_err=None, kind=None, corrected=False),
                sfr=dict(sfr=None, applies=False), balmer=None, bpt=None, ha_flux=None, broad=dict(flag=False))


def galaxy_analysis(wave, flux, err=None, *, frame='vacuum', ebv_mw=0.0, inst_fwhm_a=0.0, snr_min=4.0,
                    kernel_sigma_px=2.0, cont_width=101, z_known=None, zmax=7.0, h0=70.0, om0=0.3):
    """One object: preprocess, emission and absorption lines, dispersion, type, redshift, H-alpha SFR.

    Redshift: a catalog value if z_known is given.  Otherwise emission, a CIV/C III]/Mg II
    redshift, and absorption are compared by how many lines actually fitted.  Solutions within
    1500 km/s stay together and emission is kept.  A weaker emission redshift does not replace
    a stronger absorption redshift at a different redshift.
    Type, in order: broad-line AGN (FWHM above 1200 km/s and line S/N > 5 on Ha, Hb, CIV,
    C III] or Mg II; the Balmer cut is Hao et al. 2005); else the BPT class when [OIII], Hb,
    [NII] and Ha each have flux S/N > 3; else emission-line galaxy when the adopted emission
    or ultraviolet quality is at least 2; else star when the absorption quality is at least 2
    and |z| < 0.002; else absorption galaxy; else unknown.
    SP_SIGMA is the gas dispersion when emission lines with flux S/N > 5 exist, else the stellar
    dispersion from absorption lines.  With inst_fwhm_a = 0 the value is the observed Gaussian
    sigma and the corrected value is left empty.
    """
    prep = prepare_spectrum(wave, flux, err, frame=frame, ebv_mw=ebv_mw)
    w, f, e = prep['wave'], prep['flux'], prep['err']
    if len(w) < 20:
        out = _empty('too few pixels')
        out['preprocess'] = {k: prep[k] for k in ('n_dropped', 'n_no_extinction', 'frame', 'input_frame', 'ebv_mw')}
        return out
    zk = _num(z_known)
    if zk is not None and zk <= 0:
        zk = None
    em = analyse_spectrum(w, f, e, snr_min=snr_min, kernel_sigma_px=kernel_sigma_px, cont_width=cont_width,
                          z_known=zk, zmax=zmax, line_list=LINES)
    cont, _mask = continuum(w, f, cont_width)
    dl = float(np.median(np.diff(w)))
    med = float(np.nanmedian(cont))
    ff = f - cont + med
    sigma0 = kernel_sigma_px * dl
    ez = em['redshift']
    em_q = int(ez.get('quality') or 0)
    em_z = _num(ez.get('z'))
    em_fits = []
    if em_q > 0 and em_z is not None:
        em_fits = _fit_matched(w, ff, e, ez.get('matched') or [], dl, sigma0, em_z)
        # A broad Balmer line is wider than the detection kernel, so fit Ha and Hb at the
        # redshift even when the narrow matched filter did not list them.
        have = {ft['name'] for ft in em_fits}
        extra = []
        for name in ('Ha', 'Hb'):
            if name in have:
                continue
            lam = LINE_DICT[name]
            mu = lam * (1.0 + em_z)
            if w[0] < mu < w[-1]:
                extra.append((name, lam, mu, None))
        em_fits.extend(_fit_matched(w, ff, e, extra, dl, sigma0, em_z))
    az, adet = _absorption(w, f, e, cont, dl, kernel_sigma_px, snr_min, zk, zmax)
    ab_q = int(az.get('quality') or 0)
    ab_z = _num(az.get('z'))
    ab_fits = []
    if ab_q > 0 and ab_z is not None and zk is None:
        ab_fits = _fit_matched(w, ff, e, az.get('matched') or [], dl, sigma0, ab_z)
    elif zk is not None:
        ab_fits = _fit_matched(w, ff, e, az.get('matched') or [], dl, sigma0, zk)
        ab_fits = [ft for ft in ab_fits if ft['flux'] < 0]
    if zk is None and ab_fits:
        z_fit, ze_fit, c2_fit = refine_redshift([(ft['rest'], ft['mu'], ft['mu_err']) for ft in ab_fits])
        ab_z = z_fit
        az['z'] = z_fit
        az['z_err'] = ze_fit
        az['z_chi2r'] = c2_fit

    uv = None if zk is not None else _uv_redshift(w, f, e, zmax, snr_min)
    uv_fits = list((uv or {}).get('lines') or [])
    if zk is not None:
        adopted = dict(z=zk, z_err=None, quality=4, source='catalog')
    else:
        cands = []
        if em_fits:
            cands.append(dict(z=em_z, z_err=_num(ez.get('z_err')), quality=_fit_q(len(em_fits)), n=len(em_fits), source='emission'))
        if uv:
            cands.append(dict(z=uv['z'], z_err=uv.get('z_err'), quality=uv['quality'], n=uv['n'], source='uv'))
        if ab_fits and ab_q > 0:
            cands.append(dict(z=ab_z, z_err=_num(az.get('z_err')), quality=_fit_q(len(ab_fits)), n=len(ab_fits), source='absorption'))
        adopted = choose_redshift(cands)
    if adopted['source'] in ('catalog', 'emission'):
        class_fits = em_fits
    elif adopted['source'] == 'uv':
        class_fits = uv_fits + _fit_at_redshift(w, ff, e, adopted['z'], dl, max(sigma0, 8.0 * dl), NEBULAR_LINES)
    elif adopted['source'] == 'absorption' and adopted.get('z') is not None and abs(adopted['z']) >= 0.002:
        class_fits = _fit_at_redshift(w, ff, e, adopted['z'], dl, sigma0, NEBULAR_LINES)
    else:
        class_fits = []

    lines = _index(class_fits)
    ha, hb = lines.get('Ha'), lines.get('Hb')
    o3, n2 = lines.get('[OIII]5008'), lines.get('[NII]6585')
    bpt = None
    if _positive(o3, BPT_SNR) and _positive(hb, BPT_SNR) and _positive(n2, BPT_SNR) and _positive(ha, BPT_SNR):
        log_s2 = None
        s18, s33 = lines.get('[SII]6718'), lines.get('[SII]6733')
        if _positive(s18, BPT_SNR) and _positive(s33, BPT_SNR):
            log_s2 = math.log10((s18['flux'] + s33['flux']) / ha['flux'])
        bpt = classify_bpt(math.log10(n2['flux'] / ha['flux']), math.log10(o3['flux'] / hb['flux']), log_s2)
        bpt['class'] = bpt.pop('class_')

    broad_flag = False
    broad_line = None
    broad_fwhm = None
    for name in BROAD_LINES:
        ft = lines.get(name)
        if not _positive(ft, SIGMA_SNR):
            continue
        info = line_sigma_kms(ft, inst_fwhm_a)
        if info is None:
            continue
        if info['corrected']:
            width = FWHM * info['sigma'] if info['resolved'] else 0.0
        else:
            width = ft['fwhm'] / ft['mu'] * C_KMS
        if width > BROAD_FWHM_KMS:
            broad_flag = True
            broad_line = name
            broad_fwhm = float(width)
            break

    z_ad = adopted['z']
    if broad_flag:
        kind = 'broad-line AGN'
    elif bpt is not None:
        kind = bpt['class']
    elif adopted['source'] in ('emission', 'uv', 'catalog') and adopted['quality'] >= 2:
        kind = 'emission-line galaxy'
    elif adopted['source'] == 'absorption' and z_ad is not None and abs(z_ad) < 0.002 and adopted['quality'] >= 2:
        kind = 'star'
    elif adopted['source'] == 'absorption' and adopted['quality'] >= 2:
        kind = 'absorption galaxy'
    else:
        kind = 'unknown'

    def _sigma_rows(fits, emission):
        rows = []
        for ft in fits:
            if emission and not _positive(ft, SIGMA_SNR):
                continue
            if not emission and not (ft.get('ok') and ft['flux'] < 0 and abs(_snr(ft)) > SIGMA_SNR):
                continue
            info = line_sigma_kms(ft, inst_fwhm_a)
            if info is not None:
                info = dict(info, name=ft.get('name'))
                rows.append(info)
        return rows

    gas_rows = _sigma_rows(class_fits, True)
    star_rows = _sigma_rows(ab_fits, False)
    gas_mean, gas_err, gas_n = _weighted_mean(gas_rows)
    star_mean, star_err, star_n = _weighted_mean(star_rows)
    inst = _num(inst_fwhm_a) or 0.0
    corrected = inst > 0
    if gas_rows:
        disp = dict(sigma=gas_mean, sigma_err=gas_err, kind='gas', corrected=corrected, n=gas_n)
    elif star_rows:
        disp = dict(sigma=star_mean, sigma_err=star_err, kind='stellar', corrected=corrected, n=star_n)
    else:
        disp = dict(sigma=None, sigma_err=None, kind=None, corrected=corrected, n=0)
    if not corrected:
        disp['sigma_corrected'] = None
    else:
        disp['sigma_corrected'] = disp['sigma']
    disp['gas'] = dict(sigma=gas_mean, sigma_err=gas_err, n=gas_n)
    disp['stellar'] = dict(sigma=star_mean, sigma_err=star_err, n=star_n)
    disp['inst_fwhm_a'] = inst
    disp['note'] = 'Gaussian line width. Not a stellar-template LOSVD.'

    balmer = None
    ha_flux = ha['flux'] if _positive(ha, 0) else None
    ha_err = ha['flux_err'] if ha_flux is not None else None
    f_sfr, e_sfr = ha_flux, ha_err
    if _positive(ha, BPT_SNR) and _positive(hb, BPT_SNR):
        balmer = balmer_ebv(ha['flux'], hb['flux'])
        if balmer.get('used'):
            scale = 10.0 ** (0.4 * balmer['a_ha'])
            f_sfr = ha['flux'] * scale
            e_sfr = ha['flux_err'] * scale
    sfr = halpha_sfr(f_sfr, z_ad, h0=h0, om0=om0, flux_err=e_sfr)
    sfr['applies'] = kind in ('star-forming', 'composite', 'emission-line galaxy') and sfr.get('sfr') is not None
    sfr['ha_flux_file'] = ha_flux
    sfr['ha_flux_dustcorr'] = f_sfr if balmer and balmer.get('used') else None
    if kind == 'composite':
        sfr['note'] = 'Composite: the H-alpha luminosity includes AGN light.'
    elif not sfr['applies']:
        sfr['note'] = 'H-alpha calibration is for star-forming gas. The number is withheld from the catalog column.'

    return dict(
        ok=True,
        preprocess=dict(n_pix=int(len(w)), dlambda=dl, n_dropped=prep['n_dropped'], n_no_extinction=prep['n_no_extinction'],
                        frame='vacuum', input_frame=prep['input_frame'], ebv_mw=prep['ebv_mw'],
                        milky_way='Cardelli, Clayton & Mathis 1989, Rv=3.1' if prep['ebv_mw'] > 0 else 'not applied'),
        redshift=adopted,
        emission=dict(z=em_z, z_err=_num(ez.get('z_err')), quality=em_q, n_lines=len(em_fits),
                      lines=[_slim(ft, inst_fwhm_a) for ft in em_fits]),
        absorption=dict(z=ab_z, z_err=_num(az.get('z_err')), quality=ab_q, n_lines=len(ab_fits),
                        n_detections=len(adet), lines=[_slim(ft, inst_fwhm_a) for ft in ab_fits],
                        rest_frame='vacuum, converted from laboratory air'),
        uv=None if not uv else dict(z=uv['z'], z_err=_num(uv.get('z_err')), quality=uv['quality'], n_lines=uv['n'],
                                    lines=[_slim(ft, inst_fwhm_a) for ft in uv_fits]),
        dispersion=disp,
        type=kind,
        bpt=bpt,
        broad=dict(flag=broad_flag, line=broad_line, fwhm_kms=broad_fwhm, cut_kms=BROAD_FWHM_KMS,
                   reference=('Hao et al. 2005, ApJ 629, 61' if broad_line in ('Ha', 'Hb', None)
                              else '1200 km/s width cut; Hao et al. 2005 measured Ha and Hb only')),
        balmer=balmer,
        ha_flux=ha_flux,
        sfr=sfr,
        calibrations=dict(
            air_vacuum='SDSS/Morton (VALD form)',
            milky_way='Cardelli, Clayton & Mathis 1989, Rv=3.1',
            nebular_dust='Calzetti 2000, Rv=4.05; case B Ha/Hb=2.86 at 1e4 K; the 0.44 stellar/gas factor is not applied',
            bpt='Kewley et al. 2001; Kauffmann et al. 2003; Seyfert/LINER Kewley et al. 2006, MNRAS 372, 961',
            broad_line='FWHM > 1200 km/s. Ha and Hb: Hao et al. 2005, ApJ 629, 61. CIV, C III] and Mg II use the same width cut.',
            sfr='Murphy et al. 2011 as adopted by Kennicutt & Evans 2012; Kroupa IMF; 5.37e-42 * L_Ha',
            dispersion='Gaussian sigma; instrumental FWHM removed in quadrature when given; not a template LOSVD',
            cosmology='FlatLambdaCDM', h0=float(h0), om0=float(om0), flux_unit='file',
        ),
    )
