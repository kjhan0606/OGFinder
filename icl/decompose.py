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
from sersic_fit.utils import sersic_bn


def double_sersic_1d(r, Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl):
    """Double 1D Sérsic profile: I(r) = sersic(BCG) + sersic(ICL)."""
    return sersic_1d(r, Ie_bcg, re_bcg, n_bcg) + sersic_1d(r, Ie_icl, re_icl, n_icl)


def _central_intensity(Ie, n):
    """Sérsic central intensity I(0) = Ie * exp(bn)."""
    bn = sersic_bn(n)
    return Ie * np.exp(bn)


def _try_double_fit(r_arr, I_arr, I_err, p0, bounds_lower, bounds_upper,
                    maxfev=3000):
    """Attempt a double Sérsic curve_fit. Returns (popt, chi2) or None."""
    try:
        popt, _ = curve_fit(
            double_sersic_1d, r_arr, I_arr,
            p0=p0,
            bounds=(bounds_lower, bounds_upper),
            sigma=np.maximum(I_err, 1e-20),
            absolute_sigma=True,
            maxfev=maxfev
        )
    except Exception:
        return None

    Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl = popt
    # Ensure BCG is the compact component
    if re_icl < re_bcg:
        popt = np.array([Ie_icl, re_icl, n_icl, Ie_bcg, re_bcg, n_bcg])

    I_model = double_sersic_1d(r_arr, *popt)
    chi2 = float(np.sum(((I_arr - I_model) / np.maximum(I_err, 1e-20)) ** 2)
                 / max(1, len(r_arr) - 6))

    # Check that BCG dominates at centre
    Ie_b, re_b, n_b, Ie_i, re_i, n_i = popt
    I0_bcg = _central_intensity(Ie_b, n_b)
    I0_icl = _central_intensity(Ie_i, n_i)
    bcg_dominant = I0_bcg > I0_icl

    return popt, chi2, bcg_dominant


def fit_double_sersic(profile, pixel_scale=0.06, mag_zeropoint=25.0):
    """Fit double 1D Sérsic to surface brightness profile.

    Uses multi-start strategy:
    1. Inner/outer profile split for BCG/ICL initial guess
    2. Single-Sérsic based initial guess
    3. Picks the best fit where BCG dominates at the centre

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

    # Convert mu (mag/arcsec^2) -> linear intensity
    zp_area = mag_zeropoint + 2.5 * np.log10(pixel_scale ** 2)
    I_arr = 10.0 ** (-0.4 * (mu_arr - zp_area))
    I_err = 0.4 * np.log(10.0) * I_arr * mu_err_arr

    # Bounds for [Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl]
    bounds_lower = [0, 0.5, 1.0, 0, 0.5, 0.2]
    bounds_upper = [np.inf, r_arr[-1] * 2, 10.0, np.inf, r_arr[-1] * 5, 6.0]

    # --- Strategy A: inner/outer profile split ---
    # BCG dominates inner region, ICL dominates outer
    r_split = np.percentile(r_arr, 35)
    inner = r_arr <= r_split
    outer = r_arr > r_split

    candidates = []

    if inner.sum() >= 4 and outer.sum() >= 4:
        bcg_fit = fit_sersic_1d(r_arr[inner], I_arr[inner], I_err[inner])
        if bcg_fit is not None:
            # Subtract BCG from outer region to get ICL residual
            I_residual = I_arr[outer] - sersic_1d(
                r_arr[outer], bcg_fit['Ie'], bcg_fit['re'], bcg_fit['n'])
            I_residual = np.maximum(I_residual, I_arr[outer] * 0.01)
            icl_fit = fit_sersic_1d(r_arr[outer], I_residual, I_err[outer])

            if icl_fit is not None:
                p0_a = [bcg_fit['Ie'], max(bcg_fit['re'], 1.0),
                        max(bcg_fit['n'], 1.0),
                        icl_fit['Ie'], max(icl_fit['re'], 10.0),
                        max(min(icl_fit['n'], 5.0), 0.3)]
                res = _try_double_fit(r_arr, I_arr, I_err, p0_a,
                                      bounds_lower, bounds_upper)
                if res is not None:
                    popt, chi2, bcg_dom = res
                    candidates.append((popt, chi2, bcg_dom, 'inner/outer'))
                    print(f"  decompose trial A (inner/outer): chi2={chi2:.4f}, "
                          f"bcg_dominant={bcg_dom}", file=sys.stderr)

    # --- Strategy B: single Sérsic scaled ---
    single = fit_sersic_1d(r_arr, I_arr, I_err)
    if single is not None:
        Ie0, re0, n0 = single['Ie'], single['re'], single['n']
    else:
        Ie0 = max(I_arr[0], 1e-10)
        re0 = max(r_arr[len(r_arr) // 3], 5.0)
        n0 = 2.0

    p0_b = [Ie0 * 2.0, max(re0 * 0.3, 1.0), max(n0, 2.0),
            Ie0 * 0.2, max(re0 * 4.0, 20.0), max(n0 * 0.4, 0.4)]
    res = _try_double_fit(r_arr, I_arr, I_err, p0_b,
                          bounds_lower, bounds_upper)
    if res is not None:
        popt, chi2, bcg_dom = res
        candidates.append((popt, chi2, bcg_dom, 'single-scaled'))
        print(f"  decompose trial B (single-scaled): chi2={chi2:.4f}, "
              f"bcg_dominant={bcg_dom}", file=sys.stderr)

    # --- Strategy C: high-n BCG + low-n ICL ---
    I_peak = I_arr[0]
    p0_c = [I_peak * 0.7, max(re0 * 0.2, 2.0), 4.0,
            I_peak * 0.05, max(re0 * 5.0, 50.0), 0.5]
    res = _try_double_fit(r_arr, I_arr, I_err, p0_c,
                          bounds_lower, bounds_upper)
    if res is not None:
        popt, chi2, bcg_dom = res
        candidates.append((popt, chi2, bcg_dom, 'high-n BCG'))
        print(f"  decompose trial C (high-n BCG): chi2={chi2:.4f}, "
              f"bcg_dominant={bcg_dom}", file=sys.stderr)

    if not candidates:
        print("decompose: all fitting trials failed", file=sys.stderr)
        return None

    # Pick best: prefer BCG-dominant solutions; among those, lowest chi2
    dominant = [(p, c, d, lbl) for p, c, d, lbl in candidates if d]
    if dominant:
        best = min(dominant, key=lambda x: x[1])
    else:
        # No BCG-dominant solution — pick lowest chi2
        best = min(candidates, key=lambda x: x[1])
        print("  decompose: WARNING — no fit where BCG dominates at centre",
              file=sys.stderr)
    popt, chi2, bcg_dom, label = best
    print(f"  decompose: selected '{label}' trial", file=sys.stderr)

    Ie_bcg, re_bcg, n_bcg, Ie_icl, re_icl, n_icl = popt

    # Find R_transition: where BCG(r) == ICL(r)
    R_transition = float('nan')
    try:
        def diff(r):
            return float(sersic_1d(r, Ie_bcg, re_bcg, n_bcg)
                         - sersic_1d(r, Ie_icl, re_icl, n_icl))

        r_search = np.logspace(0, np.log10(r_arr[-1] * 2), 500)
        d_vals = np.array([diff(r) for r in r_search])
        sign_changes = np.where(np.diff(np.sign(d_vals)))[0]

        if len(sign_changes) > 0:
            idx = sign_changes[0]
            R_transition = float(brentq(diff, r_search[idx],
                                        r_search[idx + 1]))
        else:
            # No exact crossing — find closest approach as approximate
            abs_d = np.abs(d_vals)
            idx_min = int(np.argmin(abs_d))
            if abs_d[idx_min] < 0.1 * max(I_arr[0], 1e-20):
                R_transition = float(r_search[idx_min])
                print(f"  decompose: no exact crossing, using closest "
                      f"approach at r={R_transition:.1f}px", file=sys.stderr)
    except Exception:
        pass

    # Compute total fluxes
    flux_bcg = sersic_total_flux(Ie_bcg, re_bcg, n_bcg, 0.0)
    flux_icl = sersic_total_flux(Ie_icl, re_icl, n_icl, 0.0)

    if flux_bcg is None or not np.isfinite(flux_bcg):
        flux_bcg = float('nan')
    if flux_icl is None or not np.isfinite(flux_icl):
        flux_icl = float('nan')

    f_total = flux_bcg + flux_icl
    f_ICL_sersic = (flux_icl / f_total
                    if f_total > 0 and np.isfinite(f_total)
                    else float('nan'))

    re_bcg_arcsec = re_bcg * pixel_scale
    re_icl_arcsec = re_icl * pixel_scale
    R_trans_arcsec = (R_transition * pixel_scale
                      if np.isfinite(R_transition) else float('nan'))

    print(f"  BCG: n={n_bcg:.2f}, re={re_bcg:.1f}px ({re_bcg_arcsec:.2f}\")",
          file=sys.stderr)
    print(f"  ICL: n={n_icl:.2f}, re={re_icl:.1f}px ({re_icl_arcsec:.2f}\")",
          file=sys.stderr)
    print(f"  R_transition={R_transition:.1f}px, f_ICL={f_ICL_sersic:.4f}, "
          f"chi2={chi2:.4f}", file=sys.stderr)

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
