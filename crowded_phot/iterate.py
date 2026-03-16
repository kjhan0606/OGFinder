"""Iterative detection-fitting-subtraction loop for crowded field photometry."""

import numpy as np
from .config import CrowdedPhotConfig
from .group import group_sources
from .nstar import fit_group
from .psf_subtract import subtract_sources


def crowded_photometry(data, psf, sources, cfg=None):
    """Perform crowded field photometry with iteration.

    Parameters
    ----------
    data : 2D array, image data
    psf : 2D array, PSF image
    sources : list of dicts with x, y, number, kron_radius, a
    cfg : CrowdedPhotConfig

    Returns
    -------
    list of result dicts
    """
    if cfg is None:
        cfg = CrowdedPhotConfig()

    try:
        import sep
    except ImportError:
        return []

    psf_sum = np.sum(psf)
    all_results = {}
    current_data = data.copy()

    for iteration in range(cfg.max_iterations):
        # Group sources
        groups = group_sources(sources, radius_scale=cfg.group_radius_scale)

        # Fit each group
        iter_results = []
        for group_idx in groups:
            group_sources_list = [sources[i] for i in group_idx]

            if len(group_sources_list) == 0:
                continue

            # Limit group size
            if len(group_sources_list) > cfg.max_group_size:
                group_sources_list = group_sources_list[:cfg.max_group_size]

            results = fit_group(current_data, psf, group_sources_list,
                                 max_shift=cfg.max_shift)
            iter_results.extend(results)

        # Compute magnitudes
        for r in iter_results:
            flux = r.get('FLUX_CROWD', np.nan)
            if np.isfinite(flux) and flux > 0:
                r['MAG_CROWD'] = -2.5 * np.log10(flux) + cfg.mag_zeropoint

        # Store results
        for r in iter_results:
            all_results[r['NUMBER']] = r

        # Subtract and detect new sources
        if iteration < cfg.max_iterations - 1:
            residual = subtract_sources(current_data, psf, iter_results)

            # Detect new sources in residual
            mask = ~np.isfinite(residual)
            residual[mask] = 0.0
            try:
                bkg = sep.Background(residual, bw=64, bh=64)
                res_sub = residual - bkg
                new_objects = sep.extract(res_sub, cfg.detect_thresh,
                                           err=bkg.rms(),
                                           minarea=cfg.detect_minarea)
            except Exception:
                break

            if len(new_objects) == 0:
                break

            # Add new sources
            max_num = max(s['number'] for s in sources)
            n_new = 0
            for obj in new_objects:
                # Skip if too close to existing
                x_new = obj['x']
                y_new = obj['y']
                too_close = False
                for s in sources:
                    if (s['x'] - x_new)**2 + (s['y'] - y_new)**2 < 4.0:
                        too_close = True
                        break
                if too_close:
                    continue

                max_num += 1
                sources.append({
                    'number': max_num,
                    'x': x_new,
                    'y': y_new,
                    'kron_radius': 3.5,
                    'a': max(obj['a'], 1.0),
                })
                n_new += 1

            if n_new == 0:
                break

            current_data = residual

    return list(all_results.values())
