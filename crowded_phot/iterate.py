"""Iterative detection-fitting-subtraction loop for crowded field photometry."""

import numpy as np
from .config import CrowdedPhotConfig
from .group import group_sources
from .nstar import fit_group
from .psf_subtract import subtract_sources


def _worker_fit_group(args):
    """Worker for parallel_map: fit one group from shared data + PSF."""
    data_spec, psf_spec, group_sources_list, max_shift = args
    shm_data = None
    shm_psf = None
    try:
        data, shm_data = data_spec.attach()
        psf, shm_psf = psf_spec.attach()
        return fit_group(data, psf, group_sources_list, max_shift=max_shift)
    except Exception:
        return []
    finally:
        if shm_data is not None:
            shm_data.close()
        if shm_psf is not None:
            shm_psf.close()


def crowded_photometry(data, psf, sources, cfg=None, n_workers=0):
    """Perform crowded field photometry with iteration.

    Parameters
    ----------
    data : 2D array, image data
    psf : 2D array, PSF image
    sources : list of dicts with x, y, number, kron_radius, a
    cfg : CrowdedPhotConfig
    n_workers : int
        0 = auto, 1 = sequential, N = N workers

    Returns
    -------
    list of result dicts
    """
    if cfg is None:
        cfg = CrowdedPhotConfig()

    try:
        import sep
    except ImportError:
        try:
            import sep_pjw as sep
        except ImportError:
            return []

    from parallel import parallel_map, resolve_n_workers, SharedArray

    actual = resolve_n_workers(n_workers)
    psf_sum = np.sum(psf)
    all_results = {}
    current_data = data.copy()

    for iteration in range(cfg.max_iterations):
        # Group sources
        groups = group_sources(sources, radius_scale=cfg.group_radius_scale)

        # Build group task list
        group_tasks = []
        for group_idx in groups:
            group_sources_list = [sources[i] for i in group_idx]
            if len(group_sources_list) == 0:
                continue
            if len(group_sources_list) > cfg.max_group_size:
                group_sources_list = group_sources_list[:cfg.max_group_size]
            group_tasks.append(group_sources_list)

        # Fit groups — parallel within each iteration
        if actual == 1 or len(group_tasks) <= 1:
            iter_results = []
            for group_sources_list in group_tasks:
                results = fit_group(current_data, psf, group_sources_list,
                                     max_shift=cfg.max_shift)
                iter_results.extend(results)
        else:
            with SharedArray(current_data) as data_spec, \
                 SharedArray(psf) as psf_spec:
                tasks = [(data_spec, psf_spec, gs, cfg.max_shift)
                         for gs in group_tasks]
                group_results = parallel_map(
                    _worker_fit_group, tasks, n_workers=n_workers,
                    label=f"Crowded iter {iteration+1}")
            iter_results = []
            for gr in group_results:
                iter_results.extend(gr)

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
