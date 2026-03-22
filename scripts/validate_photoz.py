#!/usr/bin/env python3
"""
Photo-z validation: compare photo-z estimates against spectroscopic redshifts.

Computes standard photo-z quality metrics and generates diagnostic plots.

Metrics:
- bias: median(Δz/(1+z_spec))
- σ_MAD = 1.4826 × MAD(Δz/(1+z_spec))  (normalized MAD)
- σ_NMAD: same as σ_MAD (standard naming)
- η (outlier fraction): |Δz/(1+z_spec)| > 0.15
- σ_68: 68th percentile of |Δz/(1+z_spec)|

Usage:
    # Validate empirical method (baseline)
    python3 scripts/validate_photoz.py --data photo_z/data/sdss_specphoto.h5 --method empirical

    # Validate trained MDN model
    python3 scripts/validate_photoz.py --data photo_z/data/sdss_specphoto.h5 \
        --checkpoint photo_z/data/checkpoints/photo_z_mdn_best.pt

    # Compare both methods
    python3 scripts/validate_photoz.py --data photo_z/data/sdss_specphoto.h5 --compare
"""

import os
import sys
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def load_data(h5_path):
    """Load validation data from HDF5."""
    import h5py
    with h5py.File(h5_path, 'r') as f:
        features = f['features'][:]
        spec_z = f['redshift'][:]
        feat_names = list(f.attrs.get('feature_names', []))
        bands = list(f.attrs.get('bands', ['u', 'g', 'r', 'i', 'z']))

        extra = {}
        if 'extra' in f:
            for key in f['extra']:
                extra[key] = f['extra'][key][:]

    return features, spec_z, feat_names, bands, extra


def compute_metrics(z_phot, z_spec, label=""):
    """Compute standard photo-z quality metrics."""
    valid = np.isfinite(z_phot) & np.isfinite(z_spec) & (z_spec > 0)
    z_phot = z_phot[valid]
    z_spec = z_spec[valid]
    n = len(z_spec)

    if n == 0:
        print(f"  No valid photo-z estimates")
        return {}

    # Δz / (1 + z_spec)
    dz_norm = (z_phot - z_spec) / (1 + z_spec)

    # Metrics
    bias = np.median(dz_norm)
    sigma_mad = 1.4826 * np.median(np.abs(dz_norm - bias))
    sigma_68 = np.percentile(np.abs(dz_norm), 68)
    outlier_frac = np.mean(np.abs(dz_norm) > 0.15)
    mean_bias = np.mean(dz_norm)
    rms = np.sqrt(np.mean(dz_norm**2))

    # Print results
    header = f"=== Photo-z Metrics{': ' + label if label else ''} ==="
    print(f"\n{header}")
    print(f"  N galaxies:       {n}")
    print(f"  Bias (median):    {bias:.5f}")
    print(f"  Bias (mean):      {mean_bias:.5f}")
    print(f"  σ_MAD (NMAD):     {sigma_mad:.5f}")
    print(f"  σ_68:             {sigma_68:.5f}")
    print(f"  RMS:              {rms:.5f}")
    print(f"  η (outlier >0.15): {outlier_frac:.4f} ({int(outlier_frac*n)}/{n})")

    # Metrics by redshift bin
    z_bins = [(0.0, 0.1), (0.1, 0.2), (0.2, 0.4), (0.4, 0.7), (0.7, 1.0)]
    print(f"\n  Per redshift bin:")
    print(f"  {'z range':>12s}  {'N':>6s}  {'bias':>8s}  {'σ_MAD':>8s}  {'η':>8s}")
    for zlo, zhi in z_bins:
        mask = (z_spec >= zlo) & (z_spec < zhi)
        if np.sum(mask) < 10:
            continue
        dz_bin = dz_norm[mask]
        b = np.median(dz_bin)
        s = 1.4826 * np.median(np.abs(dz_bin - b))
        o = np.mean(np.abs(dz_bin) > 0.15)
        print(f"  [{zlo:.1f}, {zhi:.1f})  {np.sum(mask):6d}  {b:+8.5f}  {s:8.5f}  {o:8.4f}")

    # Metrics by magnitude bin
    print(f"\n  Per r-mag bin (using feature index 2):")

    return {
        'n': n, 'bias': bias, 'sigma_mad': sigma_mad, 'sigma_68': sigma_68,
        'outlier_frac': outlier_frac, 'rms': rms, 'z_phot': z_phot,
        'z_spec': z_spec, 'dz_norm': dz_norm,
    }


def run_empirical(features, bands):
    """Run empirical photo-z estimation."""
    from ds9.library.ds9_photo_z import estimate_photoz_empirical

    # Extract magnitudes (first 5 features) and errors (next 5)
    n_bands = len(bands)
    mags = features[:, :n_bands]
    mag_errs = features[:, n_bands:2*n_bands]

    z_point, z_err, z_q68, z_outlier = estimate_photoz_empirical(
        mags, mag_errs, bands)

    return z_point, z_err


def run_mdn(features, bands, checkpoint):
    """Run MDN photo-z estimation."""
    from photo_z.models.mlp_mdn import PhotoZMDN
    from photo_z.config import PhotoZConfig
    import torch

    cfg = PhotoZConfig()
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load checkpoint
    ckpt = torch.load(checkpoint, map_location=device, weights_only=False)
    ckpt_cfg = ckpt.get('config', {})

    n_features = features.shape[1]
    n_components = ckpt_cfg.get('n_components', cfg.n_components)
    quality_head = ckpt_cfg.get('quality_head', False)

    model = PhotoZMDN(
        n_features=n_features,
        n_hidden=ckpt_cfg.get('n_hidden', cfg.n_hidden),
        n_components=n_components,
        dropout=0.0,  # No dropout at inference
        quality_head=quality_head,
    ).to(device)

    model.load_state_dict(ckpt['model_state'])
    model.eval()

    # Normalize features
    feat_mean = np.array(ckpt.get('feat_mean', np.zeros(n_features)), dtype=np.float32)
    feat_std = np.array(ckpt.get('feat_std', np.ones(n_features)), dtype=np.float32)
    feat_std = np.maximum(feat_std, 1e-8)

    features_norm = (features.astype(np.float32) - feat_mean) / feat_std
    features_tensor = torch.from_numpy(features_norm).to(device)

    # Inference
    batch_size = 1024
    z_points = []
    z_errs = []

    with torch.no_grad():
        for i in range(0, len(features_tensor), batch_size):
            batch = features_tensor[i:i+batch_size]
            z_pt, z_er, _, _ = model.predict(batch)
            z_points.append(z_pt.cpu().numpy())
            z_errs.append(z_er.cpu().numpy())

    z_point = np.concatenate(z_points)
    z_err = np.concatenate(z_errs)

    return z_point, z_err


def make_plots(results, output_dir, label=""):
    """Generate diagnostic plots."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError:
        print("WARNING: matplotlib not available, skipping plots")
        return

    z_phot = results['z_phot']
    z_spec = results['z_spec']
    dz_norm = results['dz_norm']
    tag = label.replace(' ', '_').lower() if label else 'photoz'

    os.makedirs(output_dir, exist_ok=True)

    # 1. z_phot vs z_spec scatter
    fig, ax = plt.subplots(1, 1, figsize=(7, 6))
    ax.hexbin(z_spec, z_phot, gridsize=80, cmap='viridis', mincnt=1)
    ax.plot([0, 1.2], [0, 1.2], 'r--', lw=1.5, label='1:1')
    ax.plot([0, 1.2], [0.15, 1.2*1.15+0.15], 'r:', lw=0.8)
    ax.plot([0, 1.2], [-0.15, 1.2*0.85-0.15], 'r:', lw=0.8)
    ax.set_xlabel('Spectroscopic Redshift', fontsize=12)
    ax.set_ylabel('Photometric Redshift', fontsize=12)
    ax.set_title(f'Photo-z vs Spec-z ({label})\n'
                 f'σ_MAD={results["sigma_mad"]:.4f}, '
                 f'η={results["outlier_frac"]:.3f}', fontsize=11)
    ax.set_xlim(0, min(1.2, z_spec.max() * 1.1))
    ax.set_ylim(0, min(1.5, z_phot.max() * 1.1))
    ax.legend()
    plt.colorbar(ax.collections[0], ax=ax, label='Count')
    plt.tight_layout()
    path1 = os.path.join(output_dir, f'{tag}_scatter.png')
    plt.savefig(path1, dpi=150)
    plt.close()
    print(f"  Saved: {path1}")

    # 2. Δz/(1+z) histogram
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    ax.hist(dz_norm, bins=100, range=(-0.3, 0.3), density=True,
            alpha=0.7, color='steelblue', edgecolor='navy', linewidth=0.5)
    ax.axvline(0, color='red', ls='--', lw=1.5)
    ax.axvline(results['bias'], color='orange', ls='-', lw=1.5,
               label=f'bias={results["bias"]:.4f}')
    ax.axvline(0.15, color='red', ls=':', lw=0.8)
    ax.axvline(-0.15, color='red', ls=':', lw=0.8)
    ax.set_xlabel('Δz / (1 + z_spec)', fontsize=12)
    ax.set_ylabel('Density', fontsize=12)
    ax.set_title(f'Photo-z Residuals ({label})', fontsize=11)
    ax.legend()
    plt.tight_layout()
    path2 = os.path.join(output_dir, f'{tag}_residuals.png')
    plt.savefig(path2, dpi=150)
    plt.close()
    print(f"  Saved: {path2}")

    # 3. σ_MAD vs z_spec
    z_edges = np.arange(0, 1.05, 0.05)
    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    sigma_bins = []
    for j in range(len(z_edges) - 1):
        mask = (z_spec >= z_edges[j]) & (z_spec < z_edges[j+1])
        if np.sum(mask) > 20:
            dz_bin = dz_norm[mask]
            b = np.median(dz_bin)
            sigma_bins.append(1.4826 * np.median(np.abs(dz_bin - b)))
        else:
            sigma_bins.append(np.nan)

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    ax.plot(z_centers, sigma_bins, 'o-', color='steelblue', ms=5)
    ax.axhline(results['sigma_mad'], color='red', ls='--', lw=1,
               label=f'overall σ_MAD={results["sigma_mad"]:.4f}')
    ax.set_xlabel('Spectroscopic Redshift', fontsize=12)
    ax.set_ylabel('σ_MAD', fontsize=12)
    ax.set_title(f'σ_MAD vs Redshift ({label})', fontsize=11)
    ax.set_ylim(0, max(0.1, np.nanmax(sigma_bins) * 1.3))
    ax.legend()
    plt.tight_layout()
    path3 = os.path.join(output_dir, f'{tag}_sigma_vs_z.png')
    plt.savefig(path3, dpi=150)
    plt.close()
    print(f"  Saved: {path3}")


def main():
    parser = argparse.ArgumentParser(
        description='Validate photo-z accuracy against spec-z')
    parser.add_argument('--data', required=True,
                        help='Path to HDF5 with features + spec-z')
    parser.add_argument('--method', default='empirical',
                        choices=['empirical', 'mdn'],
                        help='Photo-z method to evaluate')
    parser.add_argument('--checkpoint', default=None,
                        help='MDN checkpoint (required for --method mdn)')
    parser.add_argument('--compare', action='store_true',
                        help='Run both methods and compare')
    parser.add_argument('--plot-dir', default=None,
                        help='Output directory for plots')
    parser.add_argument('--test-fraction', type=float, default=0.2,
                        help='Fraction of data to use for testing (default: 0.2)')
    args = parser.parse_args()

    if args.plot_dir is None:
        args.plot_dir = os.path.join(_project_root, 'results', 'photoz_validation')

    # Load data
    print(f"Loading data: {args.data}")
    features, spec_z, feat_names, bands, extra = load_data(args.data)
    print(f"  Total: {len(spec_z)} galaxies, {len(feat_names)} features")
    print(f"  Bands: {bands}")

    # Train/test split (use last test_fraction as test set)
    n = len(spec_z)
    np.random.seed(42)
    idx = np.random.permutation(n)
    n_test = int(n * args.test_fraction)
    test_idx = idx[:n_test]

    test_features = features[test_idx]
    test_spec_z = spec_z[test_idx]
    print(f"  Test set: {n_test} galaxies")

    if args.compare:
        # Run both methods
        print("\n--- Empirical Method ---")
        z_emp, z_emp_err = run_empirical(test_features, bands)
        res_emp = compute_metrics(z_emp, test_spec_z, label="Empirical")
        if res_emp:
            make_plots(res_emp, args.plot_dir, label="Empirical")

        # Check for MDN checkpoint
        ckpt_path = args.checkpoint
        if ckpt_path is None:
            ckpt_path = os.path.join(_project_root, 'photo_z', 'data',
                                      'checkpoints', 'photo_z_mdn_best.pt')

        if os.path.exists(ckpt_path):
            print("\n--- MDN Model ---")
            z_mdn, z_mdn_err = run_mdn(test_features, bands, ckpt_path)
            res_mdn = compute_metrics(z_mdn, test_spec_z, label="MDN")
            if res_mdn:
                make_plots(res_mdn, args.plot_dir, label="MDN")

            # Comparison summary
            print("\n=== Comparison ===")
            print(f"  {'Metric':<20s}  {'Empirical':>12s}  {'MDN':>12s}  {'Improvement':>12s}")
            for m in ['bias', 'sigma_mad', 'sigma_68', 'outlier_frac', 'rms']:
                e = res_emp[m]
                d = res_mdn[m]
                if m == 'bias':
                    imp = f"{abs(e)-abs(d):+.5f}"
                else:
                    imp = f"{(1-d/e)*100:+.1f}%"
                print(f"  {m:<20s}  {e:12.5f}  {d:12.5f}  {imp:>12s}")
        else:
            print(f"\nMDN checkpoint not found: {ckpt_path}")
            print("Train with: cd OGFinder && python3 -m photo_z.training.train "
                  f"--data {args.data}")

    elif args.method == 'empirical':
        z_phot, z_err = run_empirical(test_features, bands)
        results = compute_metrics(z_phot, test_spec_z, label="Empirical")
        if results:
            make_plots(results, args.plot_dir, label="Empirical")

    elif args.method == 'mdn':
        if args.checkpoint is None:
            args.checkpoint = os.path.join(
                _project_root, 'photo_z', 'data', 'checkpoints',
                'photo_z_mdn_best.pt')
        if not os.path.exists(args.checkpoint):
            print(f"ERROR: Checkpoint not found: {args.checkpoint}")
            print("Train first with:")
            print(f"  python3 -m photo_z.training.train --data {args.data}")
            sys.exit(1)

        z_phot, z_err = run_mdn(test_features, bands, args.checkpoint)
        results = compute_metrics(z_phot, test_spec_z, label="MDN")
        if results:
            make_plots(results, args.plot_dir, label="MDN")

    print(f"\nPlots saved to: {args.plot_dir}")


if __name__ == '__main__':
    main()
