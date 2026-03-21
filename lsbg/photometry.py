"""Parallel aperture photometry and curve-of-growth for LSBG candidates.

Measures Kron flux, fixed aperture fluxes, and curve-of-growth
half-light radii for LSBG surface brightness estimation.
"""

import sys
import numpy as np
try:
    import sep
except ImportError:
    import sep_pjw as sep


def curve_of_growth(data, x, y, max_radius=100, n_steps=20):
    """Measure cumulative flux as a function of aperture radius.

    Parameters
    ----------
    data : 2D array
        Image data (background-subtracted).
    x, y : float
        Source center (0-based).
    max_radius : float
        Maximum aperture radius in pixels.
    n_steps : int
        Number of radial steps.

    Returns
    -------
    radii : 1D array
        Aperture radii.
    cumflux : 1D array
        Cumulative flux at each radius.
    r_eff : float
        Half-light radius in pixels (NaN if unmeasurable).
    """
    if data.dtype == np.float64 and data.flags['C_CONTIGUOUS']:
        data_c = data
    else:
        data_c = np.ascontiguousarray(data, dtype=np.float64)
    radii = np.linspace(2.0, max_radius, n_steps)
    cumflux = np.zeros(n_steps)

    for i, r in enumerate(radii):
        try:
            flux, _, _ = sep.sum_circle(data_c, [x], [y], r)
            cumflux[i] = flux[0]
        except Exception:
            cumflux[i] = cumflux[i - 1] if i > 0 else 0.0

    # Find half-light radius
    total = cumflux[-1]
    r_eff = np.nan
    if total > 0:
        half = total * 0.5
        idx = np.searchsorted(cumflux, half)
        if 0 < idx < n_steps:
            # Linear interpolation
            f0, f1 = cumflux[idx - 1], cumflux[idx]
            r0, r1 = radii[idx - 1], radii[idx]
            if f1 > f0:
                r_eff = r0 + (r1 - r0) * (half - f0) / (f1 - f0)
            else:
                r_eff = r0
        elif idx == 0:
            r_eff = radii[0]

    return radii, cumflux, r_eff


def _measure_one_source(data, rms, obj, segmap, label_idx, config):
    """Measure photometry for a single source.

    Returns a dict of measured quantities.
    """
    x = float(obj['x'])
    y = float(obj['y'])
    a = float(obj['a'])
    b = float(obj['b'])
    theta = float(obj['theta'])

    # Avoid copy if already contiguous float64 (e.g. from SharedArray)
    if data.dtype == np.float64 and data.flags['C_CONTIGUOUS']:
        data_c = data
    else:
        data_c = np.ascontiguousarray(data, dtype=np.float64)
    if rms.dtype == np.float64 and rms.flags['C_CONTIGUOUS']:
        rms_c = rms
    else:
        rms_c = np.ascontiguousarray(rms, dtype=np.float64)

    result = {
        'x': x, 'y': y, 'a': a, 'b': b, 'theta': theta,
        'ellipticity': 1.0 - b / a if a > 0 else 0.0,
    }

    # Kron radius and flux
    try:
        kronrad, krflag = sep.kron_radius(
            data_c, [x], [y], [a], [b], [theta], 6.0
        )
        kr = float(kronrad[0])
        if kr < 3.5:
            kr = 3.5
        kron_flux, kron_fluxerr, kron_flag = sep.sum_ellipse(
            data_c, [x], [y], [a], [b], [theta], 2.5 * kr,
            err=rms_c
        )
        result['kron_radius'] = kr
        result['flux_auto'] = float(kron_flux[0])
        result['fluxerr_auto'] = float(kron_fluxerr[0])
        result['flags'] = int(kron_flag[0]) | int(krflag[0])
    except Exception:
        result['kron_radius'] = 0.0
        result['flux_auto'] = 0.0
        result['fluxerr_auto'] = 0.0
        result['flags'] = 1

    # Fixed aperture fluxes
    apertures = [float(s) for s in config.phot_apertures.split(',')]
    for ap in apertures:
        key = f'flux_aper_{int(ap)}'
        try:
            flux, _, _ = sep.sum_circle(data_c, [x], [y], ap, err=rms_c)
            result[key] = float(flux[0])
        except Exception:
            result[key] = 0.0

    # Half-light radius via sep.flux_radius with local normalization.
    # Using rmax based on source shape (not Kron aperture) avoids
    # background contamination that inflates r_eff for diffuse sources.
    rmax = max(10.0, min(60.0, 3.0 * max(a, b)))
    try:
        r_half, _ = sep.flux_radius(data_c, [x], [y], [rmax], 0.5)
        r_eff = float(r_half[0])
        if not np.isfinite(r_eff) or r_eff <= 0:
            r_eff = np.nan
    except Exception:
        r_eff = np.nan
    result['r_eff_pix'] = r_eff

    # Surface brightness from flux within r_eff (self-consistent)
    r_eff_arcsec = r_eff * config.pixel_scale if np.isfinite(r_eff) else np.nan
    result['r_eff_arcsec'] = r_eff_arcsec

    if np.isfinite(r_eff) and r_eff > 0:
        try:
            flux_in_reff, _, _ = sep.sum_circle(
                data_c, [x], [y], [r_eff], err=rms_c)
            flux_half = float(flux_in_reff[0])
        except Exception:
            flux_half = 0.0
    else:
        flux_half = 0.0

    if np.isfinite(r_eff_arcsec) and r_eff_arcsec > 0 and flux_half > 0:
        area = np.pi * r_eff_arcsec ** 2
        mu_eff = (config.mag_zeropoint
                  - 2.5 * np.log10(flux_half)
                  + 2.5 * np.log10(area))
        result['mu_eff'] = mu_eff
    else:
        result['mu_eff'] = np.nan

    # Total SNR
    if result['fluxerr_auto'] > 0:
        result['snr_total'] = result['flux_auto'] / result['fluxerr_auto']
    else:
        result['snr_total'] = 0.0

    # Magnitudes
    zp = config.mag_zeropoint
    f = result['flux_auto']
    fe = result['fluxerr_auto']
    if f > 0:
        result['mag_auto'] = zp - 2.5 * np.log10(f)
        if fe > 0:
            result['magerr_auto'] = 2.5 / np.log(10) * fe / f
        else:
            result['magerr_auto'] = 99.0
    else:
        result['mag_auto'] = 99.0
        result['magerr_auto'] = 99.0

    return result


def _worker_photometry(args):
    """Parallel worker for source photometry."""
    from parallel import SharedArraySpec

    data_spec, rms_spec, obj_dict, segmap_spec, label_idx, config = args

    data_arr, data_shm = data_spec.attach()
    rms_arr, rms_shm = rms_spec.attach()
    segmap_arr, segmap_shm = segmap_spec.attach()

    try:
        result = _measure_one_source(data_arr, rms_arr, obj_dict,
                                     segmap_arr, label_idx, config)
    finally:
        data_shm.close()
        rms_shm.close()
        segmap_shm.close()

    return result


def measure_photometry(data, rms, objects, segmap, config):
    """Measure photometry for all detected sources.

    Parallelized using SharedArray + parallel_map.

    Parameters
    ----------
    data : 2D array
        Background-subtracted image.
    rms : 2D array
        RMS map.
    objects : structured array
        SEP detection catalog.
    segmap : 2D int array
        Segmentation map.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    catalog : list of dict
        Photometry results for each source.
    """
    from parallel import SharedArray, parallel_map, resolve_n_workers

    n_sources = len(objects)
    if n_sources == 0:
        return []

    actual_workers = resolve_n_workers(config.n_workers)

    # Sequential mode
    if actual_workers <= 1 or n_sources < 10:
        results = []
        for i in range(n_sources):
            r = _measure_one_source(data, rms, objects[i], segmap, i, config)
            results.append(r)
            if (i + 1) % 100 == 0:
                print(f"  photometry: {i + 1}/{n_sources}", file=sys.stderr)
        return results

    print(f"LSBG photometry: {n_sources} sources, {actual_workers} workers",
          file=sys.stderr)

    data_c = np.ascontiguousarray(data, dtype=np.float64)
    rms_c = np.ascontiguousarray(rms, dtype=np.float64)
    segmap_c = np.ascontiguousarray(segmap, dtype=np.int32)

    with (SharedArray(data_c) as data_spec,
          SharedArray(rms_c) as rms_spec,
          SharedArray(segmap_c) as segmap_spec):
        tasks = [
            (data_spec, rms_spec,
             {name: float(objects[i][name]) if objects.dtype[name] != np.int32
              else int(objects[i][name])
              for name in objects.dtype.names},
             segmap_spec, i, config)
            for i in range(n_sources)
        ]
        results = parallel_map(_worker_photometry, tasks,
                               n_workers=config.n_workers,
                               label="LSBG photometry")

    return results
