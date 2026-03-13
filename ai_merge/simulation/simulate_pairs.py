"""
Generate simulated merge/separate pairs for training.

Label strategy (Approach B: Ground Truth based):
- Merge (0): primary galaxy + internal knot/clump (correlated params)
- Separate (1): two independent Sersic galaxies

Rejection sampling: both components must be detected by SEP.
Uses multiprocessing for parallel generation across CPU cores.
"""

import os
import sys
import time
import numpy as np
import h5py
from multiprocessing import Pool, cpu_count
from functools import partial

from ..config import SimConfig
from .profiles import make_galaxy_stamp, make_pair_stamp
from ..utils.sep_runner import run_sep_on_stamp
from ..catalog.extract_features import extract_pair_features


def _rand_range(lo, hi):
    return np.random.uniform(lo, hi)


def _rand_choice(choices):
    return choices[np.random.randint(len(choices))]


def _rand_galaxy_params(cfg, stamp_size):
    """Generate random galaxy parameters centered in stamp."""
    return {
        'flux': _rand_range(*cfg.flux_range),
        're': _rand_range(*cfg.re_range),
        'sersic_n': _rand_choice(cfg.sersic_n_choices),
        'ellip': _rand_range(*cfg.ellip_range),
        'theta_rad': _rand_range(0, np.pi),
        'x_offset': 0.0,
        'y_offset': 0.0,
    }


def generate_merge_pair(cfg):
    """
    Generate a merge pair: primary galaxy + internal knot/clump.
    The secondary is correlated with the primary (sub-feature).

    Returns (image, features, params_dict) or None on failure.
    """
    stamp = cfg.stamp_size
    psf_fwhm = _rand_choice(cfg.psf_fwhm_choices)
    noise_std = _rand_range(*cfg.noise_std_range)

    for _ in range(cfg.max_retries):
        # Primary galaxy
        p1 = _rand_galaxy_params(cfg, stamp)

        # Secondary: correlated knot within primary
        flux_frac = _rand_range(*cfg.merge_flux_frac_range)
        size_frac = _rand_range(*cfg.merge_size_frac_range)
        offset_max = cfg.merge_offset_frac * p1['re']

        angle = _rand_range(0, 2 * np.pi)
        offset_r = _rand_range(0.5, offset_max)

        p2 = {
            'flux': p1['flux'] * flux_frac,
            're': max(1.0, p1['re'] * size_frac),
            'sersic_n': _rand_choice([0.5, 1.0, 2.0]),  # clumps tend to be disk-like
            'ellip': _rand_range(0.0, 0.5),
            'theta_rad': _rand_range(0, np.pi),
            'x_offset': offset_r * np.cos(angle),
            'y_offset': offset_r * np.sin(angle),
        }

        # Generate stamp
        image, clean = make_pair_stamp(stamp, p1, p2, psf_fwhm, noise_std)

        # Run SEP
        result = run_sep_on_stamp(image,
                                  detect_thresh=cfg.detect_thresh,
                                  detect_minarea=cfg.detect_minarea,
                                  deblend_nthresh=cfg.deblend_nthresh,
                                  deblend_mincont=cfg.deblend_mincont)

        if result is None or result['n'] < 2:
            continue

        # Find the two sources closest to expected positions
        cx, cy = stamp / 2.0, stamp / 2.0
        obj = result['objects']
        d_primary = (obj['x'] - cx)**2 + (obj['y'] - cy)**2
        d_secondary = (obj['x'] - (cx + p2['x_offset']))**2 + \
                      (obj['y'] - (cy + p2['y_offset']))**2

        idx_i = np.argmin(d_primary)
        idx_j = np.argmin(d_secondary)

        if idx_i == idx_j:
            d_secondary[idx_i] = np.inf
            if result['n'] < 2:
                continue
            idx_j = np.argmin(d_secondary)

        if idx_i == idx_j:
            continue

        # Extract features
        feats, _ = extract_pair_features(idx_i, idx_j, result)
        if feats is None:
            continue

        if not np.all(np.isfinite(feats)):
            continue

        return image, feats

    return None


def generate_separate_pair(cfg):
    """
    Generate a separate pair: two independent Sersic galaxies.

    Returns (image, features) or None on failure.
    """
    stamp = cfg.stamp_size
    psf_fwhm = _rand_choice(cfg.psf_fwhm_choices)
    noise_std = _rand_range(*cfg.noise_std_range)

    for _ in range(cfg.max_retries):
        # Two independent galaxies
        p1 = _rand_galaxy_params(cfg, stamp)
        p2 = _rand_galaxy_params(cfg, stamp)

        # Place them with realistic separation
        r1 = max(p1['re'], 2.0)
        r2 = max(p2['re'], 2.0)
        sep_factor = _rand_range(*cfg.sep_ratio_range)
        dist = sep_factor * (r1 + r2)

        angle = _rand_range(0, 2 * np.pi)
        p1['x_offset'] = -dist / 2 * np.cos(angle)
        p1['y_offset'] = -dist / 2 * np.sin(angle)
        p2['x_offset'] = dist / 2 * np.cos(angle)
        p2['y_offset'] = dist / 2 * np.sin(angle)

        # Check offsets don't go too far off stamp
        max_off = stamp / 2 - 10
        if abs(p1['x_offset']) > max_off or abs(p1['y_offset']) > max_off:
            continue
        if abs(p2['x_offset']) > max_off or abs(p2['y_offset']) > max_off:
            continue

        # Generate stamp
        image, clean = make_pair_stamp(stamp, p1, p2, psf_fwhm, noise_std)

        # Run SEP
        result = run_sep_on_stamp(image,
                                  detect_thresh=cfg.detect_thresh,
                                  detect_minarea=cfg.detect_minarea,
                                  deblend_nthresh=cfg.deblend_nthresh,
                                  deblend_mincont=cfg.deblend_mincont)

        if result is None or result['n'] < 2:
            continue

        # Find sources closest to expected positions
        cx, cy = stamp / 2.0, stamp / 2.0
        obj = result['objects']
        d1 = (obj['x'] - (cx + p1['x_offset']))**2 + \
             (obj['y'] - (cy + p1['y_offset']))**2
        d2 = (obj['x'] - (cx + p2['x_offset']))**2 + \
             (obj['y'] - (cy + p2['y_offset']))**2

        idx_i = np.argmin(d1)
        idx_j = np.argmin(d2)

        if idx_i == idx_j:
            d2[idx_i] = np.inf
            if result['n'] < 2:
                continue
            idx_j = np.argmin(d2)

        if idx_i == idx_j:
            continue

        # Extract features
        feats, _ = extract_pair_features(idx_i, idx_j, result)
        if feats is None:
            continue

        if not np.all(np.isfinite(feats)):
            continue

        return image, feats

    return None


# ========== Worker functions for multiprocessing ==========

def _worker_merge(seed_and_cfg):
    """Worker: generate one merge pair with unique seed."""
    seed, cfg = seed_and_cfg
    np.random.seed(seed)
    result = generate_merge_pair(cfg)
    if result is None:
        return None
    img, feats = result
    return img.astype(np.float32), feats


def _worker_separate(seed_and_cfg):
    """Worker: generate one separate pair with unique seed."""
    seed, cfg = seed_and_cfg
    np.random.seed(seed)
    result = generate_separate_pair(cfg)
    if result is None:
        return None
    img, feats = result
    return img.astype(np.float32), feats


def generate_dataset(cfg=None, verbose=True, n_workers=None):
    """
    Generate the full training dataset using multiprocessing.

    Parameters
    ----------
    cfg : SimConfig, optional
    verbose : bool
    n_workers : int, optional
        Number of parallel workers. Default: cpu_count().

    Returns
    -------
    output_path : str
        Path to the generated HDF5 file.
    """
    if cfg is None:
        cfg = SimConfig()
    if n_workers is None:
        n_workers = min(cpu_count(), 16)

    os.makedirs(cfg.output_dir, exist_ok=True)

    total = cfg.n_merge + cfg.n_separate
    images_list = []
    features_list = []
    labels_list = []

    t0 = time.time()

    if verbose:
        print(f"Using {n_workers} workers for parallel generation")

    # === Generate merge pairs (label=0) ===
    if verbose:
        print(f"Generating {cfg.n_merge} merge pairs...")

    images_m, features_m = _parallel_generate(
        _worker_merge, cfg, cfg.n_merge, n_workers, verbose, "Merge"
    )

    n_merge_actual = len(images_m)
    if verbose:
        print(f"  Merge done: {n_merge_actual} pairs in {time.time()-t0:.0f}s")

    # === Generate separate pairs (label=1) ===
    t1 = time.time()
    if verbose:
        print(f"Generating {cfg.n_separate} separate pairs...")

    images_s, features_s = _parallel_generate(
        _worker_separate, cfg, cfg.n_separate, n_workers, verbose, "Separate"
    )

    n_sep_actual = len(images_s)
    if verbose:
        print(f"  Separate done: {n_sep_actual} pairs in {time.time()-t1:.0f}s")

    n_total = n_merge_actual + n_sep_actual
    if verbose:
        print(f"Total: {n_total} pairs in {time.time()-t0:.0f}s")

    # Combine
    all_images = np.array(images_m + images_s)
    all_features = np.array(features_m + features_s)
    all_labels = np.array([0]*n_merge_actual + [1]*n_sep_actual, dtype=np.int64)

    # Save to HDF5
    output_path = cfg.output_path
    with h5py.File(output_path, 'w') as f:
        f.create_dataset('images', data=all_images, compression='gzip', compression_opts=4)
        f.create_dataset('features', data=all_features, compression='gzip')
        f.create_dataset('labels', data=all_labels)
        f.attrs['n_merge'] = n_merge_actual
        f.attrs['n_separate'] = n_sep_actual
        f.attrs['stamp_size'] = cfg.stamp_size
        f.attrs['n_features'] = 20

    if verbose:
        file_size = os.path.getsize(output_path) / (1024**2)
        print(f"Saved to {output_path} ({file_size:.1f} MB)")
        print(f"  Merge: {n_merge_actual}, Separate: {n_sep_actual}")

    return output_path


def _parallel_generate(worker_fn, cfg, n_target, n_workers, verbose, label):
    """
    Generate samples in parallel using multiprocessing Pool.
    Submits batches and collects results until n_target is reached.
    """
    images = []
    features = []
    rng = np.random.RandomState(42 if label == "Merge" else 123)
    batch_size = n_workers * 8  # oversample to account for rejection
    total_attempts = 0
    max_attempts = n_target * 10

    while len(images) < n_target and total_attempts < max_attempts:
        remaining = n_target - len(images)
        # Submit more than needed to account for failures
        n_submit = min(batch_size, remaining * 3)
        seeds = rng.randint(0, 2**31, size=n_submit)
        args = [(int(s), cfg) for s in seeds]

        with Pool(processes=n_workers) as pool:
            results = pool.map(worker_fn, args)

        for r in results:
            if r is not None and len(images) < n_target:
                img, feats = r
                images.append(img)
                features.append(feats)

        total_attempts += n_submit

        if verbose and len(images) > 0:
            elapsed = time.time() - (time.time() - 1)  # approximate
            rate = len(images) / max(total_attempts, 1)
            print(f"  {label}: {len(images)}/{n_target} "
                  f"(attempts: {total_attempts}, success rate: {rate:.1%})",
                  flush=True)

    return images, features


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--n-merge', type=int, default=None)
    parser.add_argument('--n-separate', type=int, default=None)
    args = parser.parse_args()

    cfg = SimConfig()
    if args.n_merge is not None:
        cfg.n_merge = args.n_merge
    if args.n_separate is not None:
        cfg.n_separate = args.n_separate

    generate_dataset(cfg=cfg, n_workers=args.workers)
