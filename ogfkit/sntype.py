"""Supernova type from one spectrum and a host redshift.

No template library and no learned classifier.  One photospheric velocity is
chosen on a grid.  Each ion is measured at the wavelength rest*(1-v/c), so a
single trough cannot be counted as both H-alpha and Si II 6355.

Classical types follow Filippenko 1997, ARA&A, 35, 309.  The Si II depth cut
0.35, the Ic ratio, and Ic-BL at 0.1c follow Gal-Yam 2017 (arXiv:1611.09353;
the depth notation is Silverman et al. 2012).  Narrow Balmer emission is IIn
(Schlegel 1990; Filippenko 1997).  Narrow He I without hydrogen is Ibn
(Pastorello et al. 2007, Nature, 447, 829).  II-P and II-L are not split.
91bg and 91T are not assigned.
"""
import math

import numpy as np

from .specscience import C_KMS, SIGMA_SNR, air_to_vacuum, prepare_spectrum
from .spectra import LINE_DICT, continuum, fit_line

# The galaxy menu continuum is 101 pixels and follows a broad supernova line.
CONT_WIDTH = 401
V_MIN = 2000.0
V_MAX = 35000.0
V_STEP = 500.0
VOTE_KMS = 1000.0
SI_STRONG = 0.35
NARROW_FWHM_KMS = 2000.0
NARROW_CENTROID_KMS = 1000.0
NARROW_SIGMA_KMS = 600.0
ICBL = 0.1
MIN_SAMPLES = 20
HE_NAMES = ('HeI5876', 'HeI6678', 'HeI7065')

# Vacuum Angstrom.  H and He I 5876 use the engine line list.  The others are
# the conventional air wavelengths, converted once.
FEATURES = (
    ('Ha', float(LINE_DICT['Ha']), 1.0),
    ('Hb', float(LINE_DICT['Hb']), 1.0),
    ('SiII6355', float(air_to_vacuum(6355.0)), 1.0),
    ('HeI5876', float(LINE_DICT['HeI']), 0.8),
    ('HeI6678', float(air_to_vacuum(6678.15)), 0.8),
    ('HeI7065', float(air_to_vacuum(7065.19)), 0.8),
    ('OI7774', float(air_to_vacuum(7774.0)), 0.8),
    ('SII5454', float(air_to_vacuum(5454.0)), 0.3),
    ('SII5640', float(air_to_vacuum(5640.0)), 0.3),
)
FEAT = {name: (lam, wt) for name, lam, wt in FEATURES}
SI5972 = float(air_to_vacuum(5972.0))

CALIBRATION = (
    "Continuum is a 401-pixel running median. Velocity grid 2000 to 35000 km/s "
    "in steps of 500. Depth is the median of the three pixels at rest*(1-v/c), "
    "so the window cannot slide onto a neighbouring ion. Each ion votes only "
    "within 1000 km/s of its own depth peaks, so a wing is not a second line. "
    "H-alpha and Si II 6355 therefore do not share one trough. A line is detected when at least two of "
    "those pixels have depth signal-to-noise >= 5. Hydrogen requires both "
    "H-alpha and H-beta. Helium requires at least two of He I 5876, 6678, and 7065. "
    "Si II is strong when a > 0.35 (Gal-Yam 2017; near-peak normal Type Ia, "
    "Silverman et al. 2012 notation). Ic requires Si not strong and O I 7774, "
    "including the case a(Si)/a(O I) < 1. Ic-BL requires type Ic and an Si II or "
    "O I minimum, or its half-depth blue edge, at >= 0.1c. A strong Si II line "
    "stays Ia. Narrow emission is a separate Gaussian at the host redshift: "
    "flux/err > 5, FWHM < 2000 km/s, centroid within 1000 km/s. Gal-Yam's narrow "
    "component is a few hundred km/s and the intermediate component is about "
    "2000 km/s; the 1000 km/s window rejects a photospheric peak. IIn: Schlegel 1990 "
    "and Filippenko 1997. Ibn: Pastorello et al. 2007, Nature, 447, 829, on at "
    "least two He I lines and no hydrogen. Type II is not split into II-P and II-L. "
    "Si II 5972 is recorded and is not a 91bg or 91T subtype. quality is 0 when "
    "no type is assigned and 1 when these rules assign one. Equal photospheric "
    "scores keep the lower velocity. The host redshift is required. Supernova "
    "line velocities are not used as a redshift. This is not a template "
    "cross-correlation."
)
LIMITATIONS = (
    "Phase is unknown, so the 0.35 Si II cut is the near-peak normal Type Ia cut. "
    "A trough near 6150 A is assigned at one photospheric velocity: deep Si II "
    "without H-beta is Type Ia, although that absorption in SNe Ib can be hydrogen. "
    "Depth is not a pseudo-equivalent width and not a SNID match. No template library is used."
)
NO_Z_NOTE = "A host redshift is required. Supernova line velocities are not used as a redshift."


def _num(x):
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def _empty(type_, z, quality, note):
    return dict(type=type_, z=z, quality=quality, v=None, v_line=None, a_si=None, a_ha=None,
                a_he=None, a_oi=None, hydrogen=False, helium_n=0, features={}, narrow={},
                v_grid=None, note=note, calibration=CALIBRATION)


def _measure(wr, flux, cont, err, lam, v):
    """Depth at the predicted wavelength.  Three pixels, so one spike fails."""
    pred = lam * (1.0 - v / C_KMS)
    if not np.isfinite(pred) or pred <= 0 or len(wr) < 3:
        return None
    if pred < wr[0] or pred > wr[-1]:
        return None
    j = int(np.argmin(np.abs(wr - pred)))
    i0, i1 = max(0, j - 1), min(len(wr), j + 2)
    c = cont[i0:i1]
    f = flux[i0:i1]
    e = err[i0:i1]
    ok = np.isfinite(c) & (c > 0) & np.isfinite(f) & np.isfinite(e) & (e > 0)
    if not np.isfinite(cont[j]) or cont[j] <= 0 or not np.isfinite(err[j]) or err[j] <= 0:
        return None
    depth_pix = 1.0 - f / c
    snr_pix = (c - f) / e
    good = ok & np.isfinite(depth_pix)
    if int(good.sum()) < 2:
        return None
    depth = float(np.median(depth_pix[good]))
    n_sig = int(np.sum(good & (snr_pix >= SIGMA_SNR) & (depth_pix > 0)))
    detected = bool(depth > 0 and n_sig >= 2 and float(snr_pix[j - i0]) >= SIGMA_SNR)
    return dict(a=depth, snr=float(snr_pix[j - i0]), v=float(v), n_sig=n_sig, detected=detected)


def _peaks(grid, scores, rel=0.95):
    """Separate score peaks.  Each peak is the centre of its upper core.

    rel keeps a later peak when its height is at least this fraction of the
    highest.  A feature uses a low rel so a second identification of the same
    trough still votes.  The wings between those peaks do not.
    """
    s = np.array(scores, float)
    found = []
    guard = 0
    while guard < len(grid):
        guard += 1
        if not np.any(s > 0):
            break
        i = int(np.argmax(s))
        best = float(s[i])
        half = 0.5 * best
        lo = hi = i
        while lo > 0 and s[lo - 1] >= half:
            lo -= 1
        while hi + 1 < len(s) and s[hi + 1] >= half:
            hi += 1
        core = s[lo:hi + 1] >= best * 0.98
        if not np.any(core):
            core = np.zeros(hi - lo + 1, dtype=bool)
            core[i - lo] = True
        vv = grid[lo:hi + 1][core]
        ww = s[lo:hi + 1][core]
        found.append((float(np.average(vv, weights=ww)), best))
        s[lo:hi + 1] = 0.0
    if not found:
        return []
    top = found[0][1]
    return [p for p in found if p[1] >= rel * top]


def _refine(recs, grid, v0):
    """Sub-step velocity from a parabola through the defining line's depth peak."""
    a = np.array([r['a'] if r is not None and r.get('detected') else 0.0 for r in recs])
    if not np.any(a > 0):
        return float(v0)
    i = int(np.argmin(np.abs(grid - v0)))
    lo, hi = max(0, i - 2), min(len(grid), i + 3)
    j = lo + int(np.argmax(a[lo:hi]))
    if a[j] <= 0 or j == 0 or j == len(grid) - 1:
        return float(grid[j])
    coef = np.polyfit(grid[j - 1:j + 2].astype(float), a[j - 1:j + 2], 2)
    if coef[0] >= 0:
        return float(grid[j])
    v = float(-coef[1] / (2.0 * coef[0]))
    if not np.isfinite(v) or abs(v - grid[j]) > V_STEP:
        return float(grid[j])
    return v


def _blue_edge(recs, grid, v0):
    """Velocity on the blue side where this ion's depth falls to half its peak."""
    a = np.array([r['a'] if r is not None and r.get('detected') else 0.0 for r in recs])
    i = int(np.argmin(np.abs(grid - v0)))
    lo, hi = max(0, i - 2), min(len(grid), i + 3)
    i = lo + int(np.argmax(a[lo:hi]))
    peak = float(a[i])
    if peak <= 0:
        return None
    target = 0.5 * peak
    for k in range(i, len(a) - 1):
        if a[k + 1] <= target and a[k] > target:
            frac = (a[k] - target) / (a[k] - a[k + 1])
            return float(grid[k] + frac * (grid[k + 1] - grid[k]))
    if a[-1] >= target:
        return float(grid[-1])
    return None


def _narrow(wave, flux, err, lam, z):
    """Emission Gaussian at the host redshift."""
    mu0 = lam * (1.0 + z)
    if len(wave) < 8 or wave[0] > mu0 or wave[-1] < mu0:
        return None
    sigma0 = mu0 * NARROW_SIGMA_KMS / C_KMS
    ft = fit_line(wave, flux, err, mu0, sigma0=sigma0)
    if not ft.get('ok') or not (ft.get('flux_err') and ft['flux_err'] > 0 and ft.get('mu')):
        return None
    snr = float(ft['flux'] / ft['flux_err'])
    fwhm = float(ft['fwhm'] / ft['mu'] * C_KMS)
    dv = float((ft['mu'] / mu0 - 1.0) * C_KMS)
    ok = bool(ft['flux'] > 0 and snr > SIGMA_SNR and fwhm < NARROW_FWHM_KMS and abs(dv) <= NARROW_CENTROID_KMS)
    return dict(ok=ok, snr=snr, fwhm=fwhm, dv=dv, flux=float(ft['flux']))


def _finish(type_, z, quality, v, v_line, feat, helium_n, hydrogen, narrow, v_grid, extra_note):
    def depth(name):
        rec = feat.get(name)
        if rec and rec.get('detected'):
            return rec['a']
        return None
    he = [a for a in (depth(n) for n in HE_NAMES) if a is not None]
    note = LIMITATIONS if not extra_note else LIMITATIONS + ' ' + extra_note
    return dict(type=type_, z=z, quality=quality, v=v, v_line=v_line,
                a_si=depth('SiII6355'), a_ha=depth('Ha'), a_he=(max(he) if he else None), a_oi=depth('OI7774'),
                hydrogen=hydrogen, helium_n=helium_n, features=feat, narrow=narrow, v_grid=v_grid,
                note=note, calibration=CALIBRATION)


def _photospheric_type(feat):
    """Type from absorption at one velocity.  Narrow emission is handled by the caller."""
    def detected(name):
        rec = feat.get(name)
        return bool(rec and rec.get('detected'))

    def aval(name):
        rec = feat.get(name)
        return rec['a'] if rec and rec.get('detected') else None

    hydrogen = detected('Ha') and detected('Hb')
    helium_n = sum(1 for n in HE_NAMES if detected(n))
    si = aval('SiII6355')
    oi = aval('OI7774')
    si_strong = si is not None and si > SI_STRONG
    v = v_line = None
    if hydrogen and helium_n >= 2:
        type_ = 'IIb'
        v, v_line = feat['Ha']['v'], 'Ha'
    elif hydrogen:
        type_ = 'II'
        v, v_line = feat['Ha']['v'], 'Ha'
    elif si_strong:
        type_ = 'Ia'
        v, v_line = feat['SiII6355']['v'], 'SiII6355'
    elif helium_n >= 2:
        type_ = 'Ib'
        best = max((n for n in HE_NAMES if detected(n)), key=lambda n: feat[n]['a'])
        v, v_line = feat[best]['v'], best
    elif si is not None and oi is not None and oi > 0 and si / oi < 1.0:
        type_ = 'Ic'
        v, v_line = feat['OI7774']['v'], 'OI7774'
    elif oi is not None and helium_n < 2 and not si_strong:
        type_ = 'Ic'
        v, v_line = feat['OI7774']['v'], 'OI7774'
    else:
        type_ = 'uncertain'
    fast = False
    if type_ == 'Ic':
        for name in ('SiII6355', 'OI7774'):
            rec = feat.get(name)
            if not rec or not rec.get('detected'):
                continue
            if rec['v'] >= ICBL * C_KMS or (rec.get('v_blue') is not None and rec['v_blue'] >= ICBL * C_KMS):
                fast = True
        if fast:
            type_ = 'Ic-BL'
            if detected('OI7774'):
                v, v_line = feat['OI7774']['v'], 'OI7774'
            else:
                v, v_line = feat['SiII6355']['v'], 'SiII6355'
    return dict(type=type_, v=v, v_line=v_line, hydrogen=hydrogen, helium_n=helium_n, si_strong=si_strong)


def sn_analysis(wave, flux, err=None, z_host=None, frame='vacuum', ebv_mw=0.0):
    """Type one supernova spectrum.

    z_host must be finite and > 0.  Zero, blank, and non-finite values mean the
    host redshift is missing.  The supernova expansion is not used as a redshift.
    """
    z = _num(z_host)
    if z is None or z <= 0:
        return _empty('unknown', None, 0, NO_Z_NOTE)
    prep = prepare_spectrum(wave, flux, err, frame=frame, ebv_mw=ebv_mw)
    w = prep['wave']
    f = prep['flux']
    e = prep['err']
    if len(w) < MIN_SAMPLES:
        return _empty('unknown', z, 0, 'Fewer than 20 pixels remain after preprocessing.')
    cont, _mask = continuum(w, f, CONT_WIDTH)
    if e is None:
        resid = f - cont
        sig = float(1.4826 * np.nanmedian(np.abs(resid - np.nanmedian(resid))))
        e = np.full(len(w), sig if sig > 0 else np.nan)
    else:
        e = np.asarray(e, float)
    wr = w / (1.0 + z)
    grid = np.arange(V_MIN, V_MAX + 0.5 * V_STEP, V_STEP)
    series = {}
    votes = {}
    for name, lam, _wt in FEATURES:
        recs = [_measure(wr, f, cont, e, lam, float(v)) for v in grid]
        series[name] = recs
        depth = np.array([r['a'] if r is not None and r.get('detected') else 0.0 for r in recs])
        # Every separate peak of this ion votes.  The wings between a true
        # line and a misidentification of another line do not.
        votes[name] = _peaks(grid, depth, rel=0.5)
    scores = np.zeros(len(grid))
    for i, v in enumerate(grid):
        for name, _lam, wt in FEATURES:
            rec = series[name][i]
            if rec is None or not rec['detected']:
                continue
            if any(abs(float(v) - vp) <= VOTE_KMS for vp, _a in votes[name]):
                scores[i] += wt * rec['a']
    peaks = _peaks(grid, scores)

    narrow = {}
    for name in ('Ha',) + HE_NAMES:
        got = _narrow(w, f, e, FEAT[name][0], z)
        if got is not None:
            narrow[name] = got
    narrow_h = bool(narrow.get('Ha') and narrow['Ha']['ok'])
    narrow_he = [n for n in HE_NAMES if narrow.get(n) and narrow[n]['ok']]

    candidates = []
    for v_star, score in peaks:
        feat = {}
        for name, lam, _wt in FEATURES:
            rec = _measure(wr, f, cont, e, lam, v_star)
            if rec is None:
                continue
            if rec['detected'] and name in ('SiII6355', 'OI7774', 'Ha') + HE_NAMES:
                rec = dict(rec)
                rec['v'] = _refine(series[name], grid, v_star)
                if name in ('SiII6355', 'OI7774'):
                    rec['v_blue'] = _blue_edge(series[name], grid, v_star)
            feat[name] = rec
        si_rec = _measure(wr, f, cont, e, SI5972, v_star)
        if si_rec is not None:
            feat['SiII5972'] = si_rec
        got = _photospheric_type(feat)
        n_det = sum(1 for n, _lam, _wt in FEATURES if feat.get(n) and feat[n].get('detected'))
        definite = got['type'] != 'uncertain'
        # Bin the score so a 0.01 noise difference does not beat a second line
        # or the lower-velocity reading of the same trough.
        candidates.append((definite, int(round(score / 0.02)), n_det, -v_star, got, feat, v_star))

    if narrow_h:
        feat = candidates[0][5] if candidates else {}
        si_strong = bool(candidates and candidates[0][4]['si_strong'])
        # The rival Si identification can sit on a different peak from the top score.
        if not si_strong:
            for _d, _s, _n, _vv, got, _f, _v in candidates:
                if got['si_strong']:
                    si_strong = True
                    break
        extra = ''
        if si_strong:
            extra = 'Narrow hydrogen is present with a(Si II) > 0.35; typed IIn. Ia-CSM is not a separate class here.'
        return _finish('IIn', z, 1, narrow['Ha']['dv'], 'Ha', feat, 0, True, narrow,
                       candidates[0][6] if candidates else None, extra)

    primary_h = bool(candidates and candidates[0][4]['hydrogen'])
    if len(narrow_he) >= 2 and not primary_h:
        best = max(narrow_he, key=lambda n: narrow[n]['snr'])
        feat = candidates[0][5] if candidates else {}
        return _finish('Ibn', z, 1, narrow[best]['dv'], best, feat, 0, False, narrow,
                       candidates[0][6] if candidates else None, '')

    if not candidates:
        return _finish('uncertain', z, 0, None, None, {}, 0, False, narrow, None, '')
    best = max(candidates, key=lambda c: c[:4])
    got, feat, v_star = best[4], best[5], best[6]
    # Report the refined velocity of the defining line when that line set the type.
    v = got['v']
    quality = 0 if got['type'] == 'uncertain' else 1
    return _finish(got['type'], z, quality, v, got['v_line'], feat, got['helium_n'], got['hydrogen'],
                   narrow, v_star, '')
