"""Ensemble HSC Y1 distortion calibration for one DeltaSigma bin.

See Mandelbaum et al. (2018), Appendix C, and Murata et al. (2019),
equations 4, 7, 10 and 11. Inputs are already selected and rotated into
the same tangential/cross frame. Inverse critical density uses proper units.
"""
from __future__ import annotations

import numpy as np


def source_frame_angle(lens_ra, lens_dec, source_ra, source_dec):
    """Angle of the pair axis in the source tangent plane, from east to north.

    Coordinates are in degrees. The head and tail give the same spin-2
    rotation. The source-local direction is calculated from source to lens.
    """
    angles = np.broadcast_arrays(*[np.asarray(v, dtype=float) for v in (lens_ra, lens_dec, source_ra, source_dec)])
    if any(not np.all(np.isfinite(v)) for v in angles):
        raise ValueError('nonfinite sky coordinates')
    if np.any(np.abs(angles[1]) > 90) or np.any(np.abs(angles[3]) >= 90):
        raise ValueError('invalid declination or undefined sky basis at the source pole')
    lr, ld, sr, sd = np.deg2rad(angles)
    dra = lr-sr
    # Equivalent to cos(sd)*sin(ld)-sin(sd)*cos(ld)*cos(dra),
    # but avoids differencing large terms for close pairs.
    north = np.sin(ld-sd) + 2*np.sin(sd)*np.cos(ld)*np.sin(dra/2)**2
    east = np.cos(ld)*np.sin(dra)
    if np.any(np.hypot(east, north) <= 8*np.finfo(float).eps):
        raise ValueError('pair axis undefined at coincidence or antipode')
    return np.arctan2(north, east)


def rotate_sky_components(e1, e2, phi):
    """Rotate a spin-2 field in the east/north basis, with tangential positive."""
    cos2, sin2 = np.cos(2*phi), np.sin(2*phi)
    return np.stack((-e1*cos2-e2*sin2, e1*sin2-e2*cos2), axis=-1)


def y1_shape_response(pair_weights, rms_e, m):
    """Compute ensemble responsivity and multiplicative calibration."""
    weights, rms, bias = (np.asarray(v, dtype=float) for v in (pair_weights, rms_e, m))
    if (weights.ndim != 1 or rms.shape != weights.shape or bias.shape != weights.shape
            or not np.all(np.isfinite(weights) & (weights >= 0))):
        raise ValueError('expected matching vectors and nonnegative finite pair weights')
    used = weights > 0
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError('no positive pair weight')
    if not np.all(np.isfinite(rms[used]) & (rms[used] >= 0) & np.isfinite(bias[used])):
        raise ValueError('invalid weighted shape calibration')
    normalized = weights[used]/total
    responsivity = float(1-normalized @ rms[used]**2)
    mean_m = float(normalized @ bias[used])
    if not (np.isfinite(responsivity) and responsivity > 0 and np.isfinite(mean_m) and 1+mean_m > 0):
        raise ValueError('nonpositive or invalid ensemble calibration')
    return dict(responsivity=responsivity, mean_m=mean_m, weight_sum=total, n_pairs=int(used.sum()))


def y1_bin_delta_sigma(distortion, additive, shape_weight, inverse_critical_density, rms_e, m,
                       *, lens_weight, selection_m, photoz_response):
    """Apply catalog calibration, then supplied selection and photo-z responses.

    ``distortion`` and ``additive`` have shape (pairs, 2) with tangential and
    cross columns. Supply (1+b_z) in ``photoz_response``. Each response is
    defined for the exact selection and weight convention of this bin.
    """
    arrays = [np.asarray(v, dtype=float) for v in (shape_weight, inverse_critical_density, lens_weight)]
    weight, beta, lens = arrays
    if weight.ndim != 1 or any(v.shape != weight.shape for v in arrays):
        raise ValueError('expected matching one-dimensional weight and efficiency arrays')
    if any(not np.all(np.isfinite(v) & (v >= 0)) for v in arrays):
        raise ValueError('weights and efficiencies must be finite and nonnegative')
    if not (np.isfinite(selection_m) and 1+selection_m > 0
            and np.isfinite(photoz_response) and photoz_response > 0):
        raise ValueError('selection and photo-z responses must be finite and positive')
    e, c = (np.asarray(v, dtype=float) for v in (distortion, additive))
    if e.shape != (len(beta), 2) or c.shape != e.shape:
        raise ValueError('expected two rotated components per pair')
    pair_weights = lens*weight*beta**2
    result = y1_shape_response(pair_weights, rms_e, m)
    used = pair_weights > 0
    if not np.all(np.isfinite(e[used]) & np.isfinite(c[used])):
        raise ValueError('invalid weighted distortion or additive calibration')

    # Algebraically w_pair / beta, evaluated without dividing by beta.
    numerator_weights = lens[used]*weight[used]*beta[used]/result['weight_sum']
    signal = numerator_weights @ (e[used]/(2*result['responsivity']) - c[used])
    corrected = signal / ((1+result['mean_m'])*(1+selection_m)*photoz_response)
    if not np.all(np.isfinite(corrected)):
        raise ValueError('invalid calibrated DeltaSigma')
    return dict(**result, selection_m=float(selection_m), photoz_response=float(photoz_response),
                tangential=float(corrected[0]), cross=float(corrected[1]),
                unit='Msun/physical Mpc^2')
