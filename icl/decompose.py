"""Double Sérsic BCG+ICL decomposition.

Fits a two-component 1D Sérsic model to the azimuthally-averaged
surface brightness profile, separating the BCG and ICL contributions.

References:
    Kluge+2020, Montes+2018, Brough+2024
"""

import sys
import numpy as np
from scipy.optimize import curve_fit, brentq

from lsbg.sersic import sersic_1d, fit_sersic_1d, sersic_total_flux


def double_sersic_1d(r, Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl):
    """Double 1D Sérsic profile: I(r) = sersic(BCG) + sersic(ICL)."""
    return sersic_1d(r, Ie_bcg, re_bcg, n_bcg) + sersic_1d(r, Ie_icl, re_icl, n_icl)


def fit_double_sersic(profile, pixel_scale=0.06, mag_zeropoint=25.0):
    """Fit double 1D Sérsic to surface brightness profile.

    Parameters
    ----------
    profile : list of dict
        SB profile with keys R_PIX, R_ARCSEC, MU, MU_ERR, FLUX.
    pixel_scale : float
        Arcsec per pixel.
    mag_zeropoint : float
        Photometric zero-point.

    Returns
    -------
    dict with Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl,
         R_transition, f_ICL_sersic, flux_bcg, flux_icl, chi2.
    None on failure.
    """
    # Extract valid profile points
    r_arr = []
    mu_arr = []
    mu_err_arr = []
    for p in profile:
        if p['MU'] < 99.0 and p['FLUX'] > 0:
            r_arr.append(p['R_PIX'])
            mu_arr.append(p['MU'])
            mu_err_arr.append(max(p.get('MU_ERR', 0.1), 0.01))

    if len(r_arr) < 8:
        print("decompose: too few valid profile points", file=sys.stderr)
        return None

    r_arr = np.array(r_arr)
    mu_arr = np.array(mu_arr)
    mu_err_arr = np.array(mu_err_arr)

    # Convert mu (mag/arcsec^2) → linear intensity
    # I = 10^(-0.4 * (mu - zp - 2.5*log10(pixel_scale^2)))
    zp_area = mag_zeropoint + 2.5 * np.log10(pixel_scale ** 2)
    I_arr = 10.0 ** (-0.4 * (mu_arr - zp_area))
    # Error propagation: dI/dmu = -0.4 * ln(10) * I
    I_err = 0.4 * np.log(10.0) * I_arr * mu_err_arr

    # Step 1: single Sérsic fit for initial guess
    single = fit_sersic_1d(r_arr, I_arr, I_err)
    if single is None:
        # Fallback initial guess
        Ie0 = max(I_arr[0], 1e-10)
        re0 = max(r_arr[len(r_arr) // 3], 5.0)
        n0 = 2.0
    else:
        Ie0 = single['Ie']
        re0 = single['re']
        n0 = single['n']

    # Step 2: double Sérsic initial values
    # BCG: brighter, more concentrated
    Ie_bcg0 = Ie0 * 1.5
    re_bcg0 = max(re0 * 0.5, 1.0)
    n_bcg0 = max(n0, 2.0)
    # ICL: fainter, more extended
    Ie_icl0 = Ie0 * 0.3
    re_icl0 = max(re0 * 3.0, 10.0)
    n_icl0 = max(n0 * 0.5, 0.5)

    p0 = [Ie_bcg0, re_bcg0, n_bcg0, Ie_icl0, re_icl0, n_icl0]

    bounds_lower = [0, 0.5, 1.0, 0, 0.5, 0.2]
    bounds_upper = [np.inf, r_arr[-1] * 2, 10.0, np.inf, r_arr[-1] * 5, 6.0]

    try:
        popt, pcov = curve_fit(
            double_sersic_1d, r_arr, I_arr,
            p0=p0,
            bounds=(bounds_lower, bounds_upper),
            sigma=np.maximum(I_err, 1e-20),
            absolute_sigma=True,
            maxfev=2000
        )
    except Exception as e:
        print(f"decompose: curve_fit failed: {e}", file=sys.stderr)
        return None

    Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl = popt

    # Post-fit: if re_icl < re_bcg, swap components
    if re_icl < re_bcg:
        Ie_bcg, Ie_icl = Ie_icl, Ie_bcg
        re_bcg, re_icl = re_icl, re_bcg
        n_bcg, n_icl = n_icl, n_bcg

    # Compute chi2
    I_model = double_sersic_1d(r_arr, Ie_bcg, re_bcg, n_bcg,
                                Ie_icl, re_icl, n_icl)
    residuals = (I_arr - I_model)
    chi2 = float(np.sum((residuals / np.maximum(I_err, 1e-20)) ** 2)
                 / max(1, len(r_arr) - 6))

    # Find R_transition: where BCG(r) == ICL(r)
    R_transition = float('nan')
    try:
        def diff(r):
            return float(sersic_1d(r, Ie_bcg, re_bcg, n_bcg)
                         - sersic_1d(r, Ie_icl, re_icl, n_icl))

        # Search between 1 pixel and max radius
        r_search = np.logspace(0, np.log10(r_arr[-1] * 2), 200)
        d_vals = np.array([diff(r) for r in r_search])
        sign_changes = np.where(np.diff(np.sign(d_vals)))[0]

        if len(sign_changes) > 0:
            idx = sign_changes[0]
            R_transition = float(brentq(diff, r_search[idx], r_search[idx + 1]))
    except Exception:
        pass

    # Compute total fluxes (circular, ellip=0)
    flux_bcg = sersic_total_flux(Ie_bcg, re_bcg, n_bcg, 0.0)
    flux_icl = sersic_total_flux(Ie_icl, re_icl, n_icl, 0.0)

    if flux_bcg is None or not np.isfinite(flux_bcg):
        flux_bcg = float('nan')
    if flux_icl is None or not np.isfinite(flux_icl):
        flux_icl = float('nan')

    f_total = flux_bcg + flux_icl
    f_ICL_sersic = flux_icl / f_total if f_total > 0 and np.isfinite(f_total) else float('nan')

    # Convert re to arcsec
    re_bcg_arcsec = re_bcg * pixel_scale
    re_icl_arcsec = re_icl * pixel_scale
    R_trans_arcsec = R_transition * pixel_scale if np.isfinite(R_transition) else float('nan')

    return {
        'Ie_bcg': float(Ie_bcg),
        're_bcg': float(re_bcg),
        're_bcg_arcsec': float(re_bcg_arcsec),
        'n_bcg': float(n_bcg),
        'Ie_icl': float(Ie_icl),
        're_icl': float(re_icl),
        're_icl_arcsec': float(re_icl_arcsec),
        'n_icl': float(n_icl),
        'R_transition': float(R_transition),
        'R_transition_arcsec': float(R_trans_arcsec),
        'f_ICL_sersic': float(f_ICL_sersic),
        'flux_bcg': float(flux_bcg),
        'flux_icl': float(flux_icl),
        'chi2': chi2,
    }
