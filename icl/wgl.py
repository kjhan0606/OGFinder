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

The stacked image still contains the BCG, so the product name is BCG+diffuse light.
A pure diffuse fraction names the centre cut, the surface-brightness limit, or the
model decomposition that defined it. At fixed redshift and cosmology a mass error
moves the alignment radius by one third of that logarithmic mass error. The u stack
uses the point mass; the forward model averages the profile over the caller's mass
draws. A stack on a common kpc grid is returned beside it as the comparison
coordinate. The default mean is an average of valid pixels, not an inverse-variance
mean.
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


def _prepare_image(image, center, pixel_scale, mask):
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
    return sb, good, cx, cy, pixel_scale


def _sample_plane(plane, good, cx, cy, x_pix, y_pix):
    xs = cx + np.asarray(x_pix, float)
    ys = cy + np.asarray(y_pix, float)
    coords = np.stack([ys.ravel(), xs.ravel()])
    sampled = map_coordinates(plane, coords, order=1, mode='constant', cval=0.0, prefilter=False)
    covered = map_coordinates(good, coords, order=1, mode='constant', cval=0.0, prefilter=False)
    sampled = sampled.reshape(xs.shape)
    covered = covered.reshape(xs.shape)
    # A value below 1 means a masked or off-image input pixel entered the interpolation.
    valid = covered >= 1.0 - 1e-8
    return np.where(valid, sampled, np.nan), valid


def _variance_plane(variance, image_shape, pixel_scale, good):
    """Transport a count variance onto the surface-brightness plane.

    The result is the input variance divided by the pixel area squared. Linear
    resampling of that plane does not include the covariance created by the
    interpolation, so an inverse-variance mean that uses it stays a separate product.
    """
    if variance is None:
        return None
    variance = np.asarray(variance, float)
    if variance.shape != image_shape:
        raise ValueError('variance shape does not match the image')
    out = np.zeros(image_shape, float)
    area2 = pixel_scale ** 4
    usable = (good > 0) & np.isfinite(variance) & (variance > 0)
    np.divide(variance, area2, out=out, where=usable)
    return out


def warp_surface_brightness(image, center, pixel_scale, mass, z, cosmo=None,
                            mask=None, delta=200.0, u_max=2.0, n_side=201, variance=None):
    """Resample one cluster onto u = θ / θ_Δ.

    `image` is counts per pixel. `center` is 0-based (x, y), the same frame as
    icl.profile. `pixel_scale` is arcsec per pixel. A true mask pixel, and any
    non-finite input pixel, is left out. Linear interpolation from a masked
    neighbour is left out as well, so masked light cannot leak into the mean.
    `variance`, when given, is the variance of `image` in counts squared.
    """
    sb, good, cx, cy, pixel_scale = _prepare_image(image, center, pixel_scale, mask)
    cosmo = _cosmo(cosmo)
    axis = _u_axis(u_max, n_side)
    uy, ux = np.meshgrid(axis, axis, indexing='ij')
    theta = theta_delta_arcsec(mass, z, cosmo, delta)
    out, valid = _sample_plane(sb, good, cx, cy, (ux * theta) / pixel_scale, (uy * theta) / pixel_scale)
    frame = {
        'sb': out,
        'weight': valid.astype(float),
        'u': axis,
        'theta_delta_arcsec': theta,
        'r_delta_mpc': r_delta_mpc(mass, z, cosmo, delta),
        'pixel_scale': pixel_scale,
        'delta': float(delta),
    }
    var_plane = _variance_plane(variance, sb.shape, pixel_scale, good)
    if var_plane is not None:
        sampled, _valid = _sample_plane(var_plane, good, cx, cy, (ux * theta) / pixel_scale, (uy * theta) / pixel_scale)
        frame['variance'] = sampled
    return frame


def warp_metric(image, center, pixel_scale, z, cosmo=None, mask=None,
                kpc_max=500.0, n_side=201, variance=None):
    """Resample one cluster onto a common proper-distance grid, in kpc.

    The grid does not use the estimated mass. R = θ D_A(z), with θ in radians,
    then converted from Mpc to kpc. This is the comparison coordinate for a stack
    that was also aligned on R_Δ(M).
    """
    sb, good, cx, cy, pixel_scale = _prepare_image(image, center, pixel_scale, mask)
    cosmo = _cosmo(cosmo)
    z = float(z)
    distance = float(np.asarray(cosmo.angular_diameter(z)))
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError('angular-diameter distance is not positive at this redshift')
    kpc_max = _positive('kpc_max', kpc_max)
    n_side = int(n_side)
    if n_side < 3:
        raise ValueError('n_side must be at least 3')
    axis = np.linspace(-kpc_max, kpc_max, n_side)
    ky, kx = np.meshgrid(axis, axis, indexing='ij')
    kpc_per_pixel = pixel_scale * distance * ARCSEC * 1000.0
    out, valid = _sample_plane(sb, good, cx, cy, kx / kpc_per_pixel, ky / kpc_per_pixel)
    frame = {
        'sb': out,
        'weight': valid.astype(float),
        'kpc': axis,
        'kpc_per_pixel': float(kpc_per_pixel),
        'pixel_scale': pixel_scale,
        'z': z,
    }
    var_plane = _variance_plane(variance, sb.shape, pixel_scale, good)
    if var_plane is not None:
        sampled, _valid = _sample_plane(var_plane, good, cx, cy, kx / kpc_per_pixel, ky / kpc_per_pixel)
        frame['variance'] = sampled
    return frame


def _axis_name(frame):
    has_u = 'u' in frame
    has_kpc = 'kpc' in frame
    if has_u and has_kpc:
        raise ValueError('frame has both u and kpc axes')
    if has_u:
        return 'u'
    if has_kpc:
        return 'kpc'
    raise ValueError('frame has no scaled axis')


def radial_sample(frames, n_bins=4):
    """Clusters that contribute at least one valid pixel in each radial bin.

    The bin edge is the absolute value of the scaled axis. A change between the
    innermost and outermost occupied bins means the outer profile averages a
    different set of clusters.
    """
    if not frames:
        raise ValueError('no frames to stack')
    n_bins = int(n_bins)
    if n_bins < 1:
        raise ValueError('n_bins must be positive')
    axis_name = _axis_name(frames[0])
    axis = np.asarray(frames[0][axis_name], float)
    yy, xx = np.meshgrid(axis, axis, indexing='ij')
    radius = np.hypot(xx, yy)
    limit = float(np.max(np.abs(axis)))
    edges = np.linspace(0.0, limit, n_bins + 1)
    ids = [str(frame['id']) if 'id' in frame else str(i) for i, frame in enumerate(frames)]
    rows = []
    occupied = []
    for i in range(n_bins):
        if i == n_bins - 1:
            sel = (radius >= edges[i]) & (radius <= edges[i + 1])
        else:
            sel = (radius >= edges[i]) & (radius < edges[i + 1])
        present = []
        for frame, name in zip(frames, ids):
            ww = np.asarray(frame['weight'], float)
            if np.any((ww > 0) & sel):
                present.append(name)
        if present:
            occupied.append(tuple(present))
        rows.append({
            axis_name + '_min': float(edges[i]),
            axis_name + '_max': float(edges[i + 1]),
            'n_clusters': len(present),
            'ids': present,
        })
    changes = len(occupied) >= 2 and occupied[0] != occupied[-1]
    result = {'bins': rows, 'sample_changes_with_radius': bool(changes)}
    if changes:
        result['note'] = (
            'The outermost occupied bin does not contain the same clusters as the innermost occupied bin.')
    return result


def mean_stack(frames, n_bins=4):
    """Mean of valid warped pixels. Weights are the 0/1 coverage, summed.

    The returned `sb` is that mean. When every frame carries a `variance` map on
    the same grid, `inverse_variance_sb` is computed separately and is not written
    into `sb`.
    """
    if not frames:
        raise ValueError('no frames to stack')
    axis_name = _axis_name(frames[0])
    shape = frames[0]['sb'].shape
    axis = np.asarray(frames[0][axis_name], float)
    has_var = ['variance' in frame for frame in frames]
    if any(has_var) and not all(has_var):
        raise ValueError('variance must be supplied for every frame or for none')
    acc = np.zeros(shape, float)
    weight = np.zeros(shape, float)
    iv_acc = np.zeros(shape, float)
    iv_weight = np.zeros(shape, float)
    for frame in frames:
        if _axis_name(frame) != axis_name:
            raise ValueError('frames are not on the same scaled grid')
        other = np.asarray(frame[axis_name], float)
        if frame['sb'].shape != shape or other.shape != axis.shape or not np.allclose(other, axis):
            raise ValueError('frames are not on the same scaled grid')
        ww = np.asarray(frame['weight'], float)
        acc += np.where(ww > 0, frame['sb'], 0.0) * ww
        weight += ww
        if has_var[0]:
            var = np.asarray(frame['variance'], float)
            if var.shape != shape:
                raise ValueError('variance shape does not match the warped frame')
            use = (ww > 0) & np.isfinite(var) & (var > 0) & np.isfinite(frame['sb'])
            iv = np.zeros(shape, float)
            np.divide(ww, var, out=iv, where=use)
            iv_acc += np.where(use, frame['sb'], 0.0) * iv
            iv_weight += iv
    mean = np.full(shape, np.nan)
    np.divide(acc, weight, out=mean, where=weight > 0)
    result = {'sb': mean, 'weight': weight, axis_name: axis, 'n': len(frames)}
    result['ids'] = [str(frame['id']) if 'id' in frame else str(i) for i, frame in enumerate(frames)]
    result['radial'] = radial_sample(frames, n_bins=n_bins)
    if has_var[0]:
        iv_mean = np.full(shape, np.nan)
        np.divide(iv_acc, iv_weight, out=iv_mean, where=iv_weight > 0)
        result['inverse_variance_sb'] = iv_mean
        result['inverse_variance_weight'] = iv_weight
        result['inverse_variance_note'] = (
            'Separate from sb. sb is the mean of valid pixels. '
            'A variance map that was linearly resampled does not include interpolation covariance.')
    return result


_PRODUCT = 'bcg+diffuse'
_PRODUCT_DESCRIPTION = (
    'BCG+diffuse light. The mean keeps the BCG. A pure diffuse fraction names the '
    'centre cut, surface-brightness limit, or model decomposition that defined it.')
_LIGHT_DEFINITIONS = ('centre_exclude', 'surface_brightness', 'model_decomposition')
_SAME_LENS_NOTE = (
    'The alignment mass comes from the same lensing data as the mass being compared. '
    'Agreement is not an independent check. The mass posterior and the lensing covariance '
    'belong in the same result.')
_RICHNESS_NOTE = (
    'A richness mass shares satellite counts and satellite light with the denominator '
    'of a BCG+satellite light fraction, so an apparent mass trend can repeat that selection.')
_EXTERNAL_NOTE = (
    'An external mass still needs the scatter, redshift trend and selection of its calibration. '
    'A mass defined as M_500c or M_200m needs the conversion to M_200c written beside the radius.')
_EDDINGTON_NOTE = (
    'Splitting the sample on a noisy mass moves clusters across bins and can produce an Eddington bias. '
    'This stack does not remove that shift.')
_CONVERSION_NOTE = 'The supplied mass is not M_200c. Record the conversion before comparing radii.'
_DEFINITION_NOTE = (
    'M_200c and M_200m are different mass definitions. '
    'An input mass of M_500c or M_200m keeps its conversion beside the radius.')
_SURFACE_BRIGHTNESS_NOTE = (
    'Frame values are counts per arcsec^2. '
    'That conversion leaves exposure time and zeropoint as they were. '
    'The survey gate records observed surface brightness or rest-frame luminosity surface density.')
_SURVEY_REQUIRED = (
    'luminosity_unit', 'band', 'extinction', 'sky_residual', 'enlarged_mask', 'psf',
    'radial_sample', 'covariance')


def _filled(value):
    if value is None:
        return False
    if isinstance(value, str):
        return bool(value.strip())
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    if isinstance(value, (int, float, np.integer, np.floating)):
        return bool(np.isfinite(value))
    return True


def survey_checklist(items=None, rest_frame=None):
    """Gates that have to be filled before the stack is marked as a science measurement.

    `rest_frame` is `observed` or `rest`. A rest-frame mean also needs a K correction
    and a surface-brightness dimming correction. Empty gates leave `science_measurement`
    false; the stack itself is still returned.
    """
    items = {} if items is None else dict(items)
    empty = [key for key in _SURVEY_REQUIRED if not _filled(items.get(key))]
    frame = None if rest_frame is None else str(rest_frame)
    if frame not in ('observed', 'rest'):
        empty.append('rest_frame')
    elif frame == 'rest':
        for key in ('k_correction', 'surface_brightness_dimming'):
            if not _filled(items.get(key)):
                empty.append(key)
    recorded = {key: items.get(key) for key in _SURVEY_REQUIRED}
    recorded['rest_frame'] = frame
    recorded['k_correction'] = items.get('k_correction')
    recorded['surface_brightness_dimming'] = items.get('surface_brightness_dimming')
    return {'items': recorded, 'empty': empty, 'science_measurement': len(empty) == 0}


def named_light_fraction(numerator=None, denominator=None, definition=None,
                         centre_exclude_radius=None, sb_limit=None, decomposition=None,
                         numerator_includes_bcg=False):
    """Record a light fraction under the definition the caller actually used.

    A fraction whose numerator includes the BCG is named bcg+icl. Any other fraction
    keeps the name diffuse and repeats the centre cut, surface-brightness limit, or
    model decomposition that defined it.
    """
    if definition not in _LIGHT_DEFINITIONS:
        raise ValueError('definition must be centre_exclude, surface_brightness, or model_decomposition')
    if definition == 'centre_exclude':
        if centre_exclude_radius is None:
            raise ValueError('centre_exclude needs centre_exclude_radius')
        radius = float(centre_exclude_radius)
        if not np.isfinite(radius) or radius < 0:
            raise ValueError('centre_exclude_radius must be non-negative')
        parameters = {'centre_exclude_radius': radius}
    elif definition == 'surface_brightness':
        if sb_limit is None or (isinstance(sb_limit, str) and not sb_limit.strip()):
            raise ValueError('surface_brightness needs sb_limit')
        parameters = {'sb_limit': float(sb_limit) if isinstance(sb_limit, (int, float, np.floating, np.integer)) else str(sb_limit)}
    else:
        if decomposition is None or not str(decomposition).strip():
            raise ValueError('model_decomposition needs decomposition')
        parameters = {'decomposition': str(decomposition)}
    includes = bool(numerator_includes_bcg)
    if includes:
        description = 'BCG+ICL fraction. The numerator includes the BCG.'
    else:
        description = (
            'Diffuse fraction under the caller definition. '
            'The definition parameters name the centre cut, surface-brightness limit, or model decomposition.')
    ratio = None
    if numerator is not None or denominator is not None:
        if numerator is None or denominator is None:
            raise ValueError('numerator and denominator are both required to form a ratio')
        den = float(denominator)
        num = float(numerator)
        if not np.isfinite(num) or not np.isfinite(den) or den == 0:
            raise ValueError('light fraction needs a finite numerator and a non-zero denominator')
        ratio = num / den
    return {
        'product': 'bcg+icl' if includes else 'diffuse',
        'description': description,
        'definition': definition,
        'definition_parameters': parameters,
        'numerator_includes_bcg': includes,
        'ratio': ratio,
    }


def radius_log_shift(delta_ln_mass):
    """δln R_Δ = δln M_Δ / 3 at fixed redshift, cosmology and overdensity."""
    return np.asarray(delta_ln_mass, float) / 3.0


def gaussian_ln_mass_nodes(sigma, n=7):
    """Nodes and weights for a normal distribution of δln M, from Gauss–Hermite quadrature."""
    sigma = _positive('sigma', sigma)
    n = int(n)
    if n < 1:
        raise ValueError('n must be positive')
    nodes, weights = np.polynomial.hermite.hermgauss(n)
    return sigma * np.sqrt(2.0) * nodes, weights / np.sqrt(np.pi)


def forward_scaled_profile(u, sb, delta_ln_mass, weights=None):
    """Average a radial model after each mass draw stretches R_Δ.

    The stack coordinate uses one estimated mass, u_hat = R / R(M_hat). A draw with
    δln M = ln M_draw − ln M_hat has R_draw = R_hat exp(δln M / 3), so the model
    written on the estimated grid is I(u_hat exp(−δln M / 3)). Draws are averaged
    with the caller's weights. Samples that fall outside `u` become non-finite.
    """
    u = np.asarray(u, float).reshape(-1)
    sb = np.asarray(sb, float).reshape(-1)
    if u.shape != sb.shape:
        raise ValueError('u and sb must have the same length')
    if u.size < 2 or np.any(np.diff(u) <= 0):
        raise ValueError('u must be strictly increasing')
    shifts = np.atleast_1d(np.asarray(delta_ln_mass, float))
    if weights is None:
        weights = np.ones(shifts.shape, float)
    weights = np.asarray(weights, float)
    if weights.shape != shifts.shape or not np.isfinite(weights).all() or np.any(weights < 0) or float(weights.sum()) <= 0:
        raise ValueError('mass weights must be non-negative and match the mass draws')
    if not np.isfinite(shifts).all():
        raise ValueError('mass draws must be finite')
    acc = np.zeros(u.shape, float)
    for shift, weight in zip(shifts, weights):
        drawn = np.interp(u * math.exp(-float(shift) / 3.0), u, sb, left=np.nan, right=np.nan)
        acc += float(weight) * drawn
    return acc / float(weights.sum())


def mass_interpretation(mass_origin=None, splits_sample=False, delta_ln_mass=None, weights=None,
                        lensing_covariance=None, mass_definition='200c', conversion=None,
                        u=None, model_sb=None, stack_sb=None):
    """Say what an agreement between this stack and a mass is allowed to mean.

    `independent_evidence` stays false. A same-lens mass is not an independent check.
    A supplied covariance is recorded, not propagated through a joint likelihood.
    The array is the one the caller passed; this function does not calculate the
    covariance between a mass-based alignment and the light profile.
    """
    origin = None if mass_origin is None else str(mass_origin)
    if origin not in (None, 'same_lens', 'richness', 'external'):
        raise ValueError('mass_origin must be same_lens, richness, external, or omitted')
    notes = []
    if origin == 'same_lens':
        notes.append(_SAME_LENS_NOTE)
    elif origin == 'richness':
        notes.append(_RICHNESS_NOTE)
    elif origin == 'external':
        notes.append(_EXTERNAL_NOTE)
    if splits_sample:
        notes.append(_EDDINGTON_NOTE)
    definition = str(mass_definition)
    recorded_conversion = None
    if definition not in ('200c', 'M200c'):
        recorded_conversion = None if conversion is None else str(conversion)
        if not recorded_conversion:
            notes.append(_CONVERSION_NOTE)
    smeared = None
    residual = None
    if delta_ln_mass is not None or u is not None or model_sb is not None:
        if u is None or model_sb is None or delta_ln_mass is None:
            raise ValueError('a forward model needs u, model_sb and delta_ln_mass together')
        smeared = forward_scaled_profile(u, model_sb, delta_ln_mass, weights)
        if stack_sb is not None:
            residual = np.asarray(stack_sb, float) - smeared
    covariance = None if lensing_covariance is None else np.asarray(lensing_covariance, float)
    return {
        'independent_evidence': False,
        'mass_origin': origin,
        'splits_sample': bool(splits_sample),
        'notes': notes,
        'mass_definition': definition,
        'definition_note': _DEFINITION_NOTE,
        'conversion': recorded_conversion,
        'smeared_model': smeared,
        'residual': residual,
        'lensing_covariance': covariance,
        'covariance_recorded': covariance is not None,
        'covariance_propagated': False,
    }


def area_corrected_residual(diffuse, observed, mask):
    """Extra light left after dividing by the unmasked area.

    `diffuse` is the light the measurement wanted. `observed` is that light plus
    wings and undetected satellites. `mask` is true on pixels removed from the mean.
    The residual equals the mean extra light on open pixels. Masked diffuse light is
    absent from the mean, and the open-area division does not put it back.
    """
    diffuse = np.asarray(diffuse, float)
    observed = np.asarray(observed, float)
    mask = np.asarray(mask, bool)
    if diffuse.shape != observed.shape or diffuse.shape != mask.shape:
        raise ValueError('diffuse, observed and mask must share a shape')
    open_pix = ~mask & np.isfinite(observed) & np.isfinite(diffuse)
    if not np.any(open_pix):
        raise ValueError('no open pixels')
    extra = observed - diffuse
    measured = float(np.mean(observed[open_pix]))
    truth = float(np.mean(diffuse[open_pix]))
    masked_sum = float(np.sum(diffuse[mask & np.isfinite(diffuse)])) if np.any(mask) else 0.0
    return {
        'measured': measured,
        'diffuse_on_open_pixels': truth,
        'residual': measured - truth,
        'extra_light_mean': float(np.mean(extra[open_pix])),
        'masked_diffuse_sum': masked_sum,
        'n_open': int(np.count_nonzero(open_pix)),
        'note': (
            'The residual is extra light on open pixels. '
            'Dividing by the open area does not remove wings or undetected satellites, '
            'and it does not put masked diffuse light back.'),
    }


def stack_clusters(records, cosmo=None, delta=200.0, u_max=2.0, n_side=201, kpc_max=None,
                   survey=None, rest_frame=None, mass_origin=None, splits_sample=False,
                   delta_ln_mass=None, mass_weights=None, lensing_covariance=None,
                   mass_definition='200c', conversion=None, model_u=None, model_sb=None,
                   stack_profile=None):
    """Warp each cluster by its own mass and apparent pixel scale, then average.

    Each record needs image, center, pixel_scale, mass and z. mask, id and variance
    are optional. `variance` is the variance of the image in counts squared.
    The returned object is the u stack. `kpc` is the same images on one proper-distance
    grid. `product` is bcg+diffuse. `survey` is the science-measurement gate.
    """
    if not records:
        raise ValueError('no clusters to stack')
    cosmo = _cosmo(cosmo)
    frames = []
    for i, rec in enumerate(records):
        frame = warp_surface_brightness(
            rec['image'], rec['center'], rec['pixel_scale'], rec['mass'], rec['z'],
            cosmo=cosmo, mask=rec.get('mask'), delta=delta, u_max=u_max, n_side=n_side,
            variance=rec.get('variance'))
        frame['id'] = str(rec['id']) if 'id' in rec and rec['id'] is not None else str(i)
        frames.append(frame)
    stacked = mean_stack(frames)
    stacked['frames'] = frames
    stacked['delta'] = float(delta)
    stacked['product'] = _PRODUCT
    stacked['description'] = _PRODUCT_DESCRIPTION
    stacked['surface_brightness_note'] = _SURFACE_BRIGHTNESS_NOTE
    radii = [r_delta_mpc(rec['mass'], rec['z'], cosmo, delta) for rec in records]
    if kpc_max is None:
        kpc_limit = float(u_max) * max(radii) * 1000.0
        kpc_source = 'u_max * max(R_delta)'
    else:
        kpc_limit = _positive('kpc_max', kpc_max)
        kpc_source = 'caller'
    kpc_frames = []
    for rec, frame in zip(records, frames):
        metric = warp_metric(
            rec['image'], rec['center'], rec['pixel_scale'], rec['z'],
            cosmo=cosmo, mask=rec.get('mask'), kpc_max=kpc_limit, n_side=n_side,
            variance=rec.get('variance'))
        metric['id'] = frame['id']
        kpc_frames.append(metric)
    kpc = mean_stack(kpc_frames)
    kpc['frames'] = kpc_frames
    kpc['product'] = _PRODUCT
    kpc['description'] = _PRODUCT_DESCRIPTION
    kpc['kpc_max'] = float(kpc_limit)
    kpc['kpc_max_source'] = kpc_source
    stacked['kpc'] = kpc
    stacked['survey'] = survey_checklist(survey, rest_frame)
    stacked['interpretation'] = mass_interpretation(
        mass_origin=mass_origin, splits_sample=splits_sample, delta_ln_mass=delta_ln_mass,
        weights=mass_weights, lensing_covariance=lensing_covariance, mass_definition=mass_definition,
        conversion=conversion, u=model_u, model_sb=model_sb, stack_sb=stack_profile)
    return stacked
