"""2-component Sérsic model: Bulge (n~4) + Disk (n~1) + sky."""

import numpy as np
from sersic_fit.utils import sersic_bn, ellipse_radius


def bulge_disk_2d(params, x, y, n_bulge=4.0, n_disk=1.0):
    """Evaluate 2D Bulge + Disk model.

    Parameters
    ----------
    params : array-like, length 11
        [Ie_b, re_b, e_b, theta_b,
         I0_d, rs_d, e_d, theta_d,
         xc, yc, sky]
    x, y : 2D arrays of pixel coordinates
    n_bulge : float
        Sérsic index for bulge component (default 4.0 = de Vaucouleurs)
    n_disk : float
        Sérsic index for disk component (default 1.0 = exponential)

    Returns
    -------
    model : 2D array, same shape as x
    """
    Ie_b, re_b, e_b, theta_b, I0_d, rs_d, e_d, theta_d, xc, yc, sky = params

    # Bulge component: Sérsic with n = n_bulge
    bn_b = sersic_bn(n_bulge)
    r_b = ellipse_radius(x, y, xc, yc, e_b, theta_b)
    exp_b = bn_b * ((r_b / max(re_b, 0.1)) ** (1.0 / max(n_bulge, 0.1)) - 1.0)
    exp_b = np.clip(exp_b, -50, 50)
    bulge = Ie_b * np.exp(-exp_b)

    # Disk component: Sérsic with n = n_disk
    bn_d = sersic_bn(n_disk)
    r_d = ellipse_radius(x, y, xc, yc, e_d, theta_d)
    exp_d = bn_d * ((r_d / max(rs_d, 0.1)) ** (1.0 / max(n_disk, 0.1)) - 1.0)
    exp_d = np.clip(exp_d, -50, 50)
    disk = I0_d * np.exp(-exp_d)

    return bulge + disk + sky


def compute_bt_ratio(params, x, y, n_bulge=4.0, n_disk=1.0):
    """Compute bulge-to-total (B/T) flux ratio.

    Parameters
    ----------
    params : array-like, length 11
        Same as bulge_disk_2d.
    x, y : 2D arrays of pixel coordinates
    n_bulge, n_disk : float
        Sérsic indices for bulge and disk.

    Returns
    -------
    bt : float
        B/T ratio in [0, 1]. Returns 0.0 if total flux <= 0.
    """
    Ie_b, re_b, e_b, theta_b, I0_d, rs_d, e_d, theta_d, xc, yc, sky = params

    # Bulge flux
    bn_b = sersic_bn(n_bulge)
    r_b = ellipse_radius(x, y, xc, yc, e_b, theta_b)
    exp_b = bn_b * ((r_b / max(re_b, 0.1)) ** (1.0 / max(n_bulge, 0.1)) - 1.0)
    exp_b = np.clip(exp_b, -50, 50)
    flux_bulge = np.sum(Ie_b * np.exp(-exp_b))

    # Disk flux
    bn_d = sersic_bn(n_disk)
    r_d = ellipse_radius(x, y, xc, yc, e_d, theta_d)
    exp_d = bn_d * ((r_d / max(rs_d, 0.1)) ** (1.0 / max(n_disk, 0.1)) - 1.0)
    exp_d = np.clip(exp_d, -50, 50)
    flux_disk = np.sum(I0_d * np.exp(-exp_d))

    total = flux_bulge + flux_disk
    if total <= 0:
        return 0.0

    return float(np.clip(flux_bulge / total, 0.0, 1.0))
