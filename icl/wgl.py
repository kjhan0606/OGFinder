"""Scale cluster images by apparent size and estimated mass, then stack surface brightness.

Flat ΛCDM, with H0 and Ωm taken from the caller. The defaults, 70 and 0.3, are the
same settings as the weak-lensing menu. They are not a published cosmology fit.

A pixel offset from the cluster centre has apparent angular size

    θ = (pixel scale) × (pixel offset)

in arcsec. The angular-diameter distance D_A(z) turns that angle into a proper length,

    R = θ D_A(z)

with θ in radians. The estimated mass M_Δ is the mass inside the radius where the
mean density is Δ times the critical density at the cluster redshift,

    M_Δ = (4π/3) Δ ρ_crit(z) R_Δ³ ,

so

    R_Δ = [3 M_Δ / (4π Δ ρ_crit(z))]^(1/3) .

That radius subtends

    θ_Δ = R_Δ / D_A(z) .

The scaled coordinate used for the stack is the apparent size in units of this angle,

    u = θ / θ_Δ = R / R_Δ .

Mass changes the stack through R_Δ. The apparent size changes it through θ. Redshift
enters both D_A and ρ_crit. With Δ = 200 the mass is M_200c, the same definition as
the weak-lensing mass. The output of a warp is surface brightness, counts per
arcsec², so a change of pixel scale does not rescale the light.
"""

import math

import numpy as np
from scipy.ndimage import map_coordinates

from ogfkit.wlcluster import ARCSEC, Cosmology, rho_crit


def _cosmo(cosmo):
    return Cosmology() if cosmo is None else cosmo


def _positive(name, value):
    value = float(value)
    if not np.isfinite(value) or value <= 0:
        raise ValueError('%s must be positive' % name)
    return value


def r_delta_mpc(mass, z, cosmo=None, delta=200.0):
    """Proper radius, in Mpc, of the sphere that holds `mass` at density Δ ρ_crit(z)."""
    mass = _positive('mass', mass)
    delta = _positive('delta', delta)
    z = float(z)
    if not np.isfinite(z) or z < 0:
        raise ValueError('redshift must be non-negative')
    cosmo = _cosmo(cosmo)
    rho = rho_crit(z, cosmo.H0, cosmo.Om)
    if not np.isfinite(rho) or rho <= 0:
        raise ValueError('critical density is not positive')
    return (3.0 * mass / (4.0 * math.pi * delta * rho)) ** (1.0 / 3.0)


def theta_delta_arcsec(mass, z, cosmo=None, delta=200.0):
    """Apparent angular size, in arcsec, of R_Δ at the cluster redshift."""
    cosmo = _cosmo(cosmo)
    z = float(z)
    distance = float(np.asarray(cosmo.angular_diameter(z)))
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError('angular-diameter distance is not positive at this redshift')
    return r_delta_mpc(mass, z, cosmo, delta) / distance / ARCSEC


def scaled_radius(theta_arcsec, mass, z, cosmo=None, delta=200.0):
    """u = θ / θ_Δ. θ is an apparent angular size in arcsec."""
    theta = np.asarray(theta_arcsec, float)
    scale = theta_delta_arcsec(mass, z, cosmo, delta)
    return theta / scale


def _u_axis(u_max, n_side):
    u_max = _positive('u_max', u_max)
    n_side = int(n_side)
    if n_side < 3:
        raise ValueError('n_side must be at least 3')
    return np.linspace(-u_max, u_max, n_side)


def warp_surface_brightness(image, center, pixel_scale, mass, z, cosmo=None,
                            mask=None, delta=200.0, u_max=2.0, n_side=201):
    """Resample one cluster onto u = θ / θ_Δ.

    `image` is counts per pixel. `center` is 0-based (x, y), the same frame as
    icl.profile. `pixel_scale` is arcsec per pixel. A true mask pixel, and any
    non-finite input pixel, is left out. Linear interpolation from a masked
    neighbour is left out as well, so masked light cannot leak into the mean.
    """
    image = np.asarray(image, float)
    if image.ndim != 2:
        raise ValueError('image must be two-dimensional')
    pixel_scale = _positive('pixel_scale', pixel_scale)
    cx, cy = float(center[0]), float(center[1])
    bad = ~np.isfinite(image)
    if mask is not None:
        mask = np.asarray(mask, bool)
        if mask.shape != image.shape:
            raise ValueError('mask shape does not match the image')
        bad |= mask
    good = (~bad).astype(float)
    sb = np.zeros(image.shape, float)
    np.divide(image, pixel_scale * pixel_scale, out=sb, where=good > 0)

    cosmo = _cosmo(cosmo)
    axis = _u_axis(u_max, n_side)
    uy, ux = np.meshgrid(axis, axis, indexing='ij')
    theta = theta_delta_arcsec(mass, z, cosmo, delta)
    xs = cx + (ux * theta) / pixel_scale
    ys = cy + (uy * theta) / pixel_scale
    coords = np.stack([ys.ravel(), xs.ravel()])
    sampled = map_coordinates(sb, coords, order=1, mode='constant', cval=0.0, prefilter=False)
    covered = map_coordinates(good, coords, order=1, mode='constant', cval=0.0, prefilter=False)
    sampled = sampled.reshape(ux.shape)
    covered = covered.reshape(ux.shape)
    # A value below 1 means a masked or off-image input pixel entered the interpolation.
    valid = covered >= 1.0 - 1e-8
    out = np.where(valid, sampled, np.nan)
    return {
        'sb': out,
        'weight': valid.astype(float),
        'u': axis,
        'theta_delta_arcsec': theta,
        'r_delta_mpc': r_delta_mpc(mass, z, cosmo, delta),
        'pixel_scale': pixel_scale,
        'delta': float(delta),
    }


def mean_stack(frames):
    """Inverse-mask mean of warped surface-brightness frames on one u grid."""
    if not frames:
        raise ValueError('no frames to stack')
    shape = frames[0]['sb'].shape
    axis = frames[0]['u']
    acc = np.zeros(shape, float)
    weight = np.zeros(shape, float)
    for frame in frames:
        if frame['sb'].shape != shape or frame['u'].shape != axis.shape or not np.allclose(frame['u'], axis):
            raise ValueError('frames are not on the same scaled grid')
        ww = np.asarray(frame['weight'], float)
        acc += np.where(ww > 0, frame['sb'], 0.0) * ww
        weight += ww
    mean = np.full(shape, np.nan)
    np.divide(acc, weight, out=mean, where=weight > 0)
    return {'sb': mean, 'weight': weight, 'u': axis, 'n': len(frames)}


def stack_clusters(records, cosmo=None, delta=200.0, u_max=2.0, n_side=201):
    """Warp each cluster by its own mass and apparent pixel scale, then average.

    Each record needs image, center, pixel_scale, mass and z. mask is optional.
    """
    if not records:
        raise ValueError('no clusters to stack')
    cosmo = _cosmo(cosmo)
    frames = []
    for rec in records:
        frames.append(warp_surface_brightness(
            rec['image'], rec['center'], rec['pixel_scale'], rec['mass'], rec['z'],
            cosmo=cosmo, mask=rec.get('mask'), delta=delta, u_max=u_max, n_side=n_side))
    stacked = mean_stack(frames)
    stacked['frames'] = frames
    stacked['delta'] = float(delta)
    return stacked
