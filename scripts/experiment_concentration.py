#!/usr/bin/env python3
"""
Experiment: test concentration index (C=R90/R50) and size-magnitude
relation as additional features for photo-z.

Features to test:
1. Concentration index: C = petroR90_r / petroR50_r
2. Size-mag relation: log10(petroR50_r) - 0.2 * mag_r
   (at fixed absolute mag, nearby galaxies appear larger)

Usage:
    python3 scripts/experiment_concentration.py
"""

import os
import sys
import numpy as np
import time

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

import torch
from photo_z.models.mlp_mdn import PhotoZMDN
from photo_z.training.data import SDSSPhotoZDataset
from torch.utils.data import DataLoader, Subset


def load_data_with_extras(h5_path):
    """Load HDF5 and return features, redshifts, extra columns."""
    import h5py
    with h5py.File(h5_path, 'r') as f:
        features = f['features'][:]
        spec_z = f['redshift'][:]
        feat_names = list(f.attrs.get('feature_names', []))
        extra = {}
        if 'extra' in f:
            for key in f['extra']:
                extra[key] = f['extra'][key][:]
    return features, spec_z, feat_names, extra


def train_and_evaluate(features, spec_z, label, seed=42, epochs=100,
                       batch_size=256, lr=1e-3, patience=20):
    """Train MDN and return test metrics."""
    n = len(spec_z)
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)

    n_test = int(n * 0.2)
    n_val = int(n * 0.1)
    n_train = n - n_val - n_test

    train_idx = idx[:n_train]
    val_idx = idx[n_train:n_train + n_val]
    test_idx = idx[n_train + n_val:]

    # Normalize from training set
    feat_mean = features[train_idx].mean(axis=0)
    feat_std = features[train_idx].std(axis=0)
    feat_std[feat_std < 1e-8] = 1.0
    features_norm = (features - feat_mean) / feat_std

    full_ds = SDSSPhotoZDataset(features_norm, spec_z)
    train_ds = Subset(full_ds, train_idx.tolist())
    val_ds = Subset(full_ds, val_idx.tolist())
    test_ds = Subset(full_ds, test_idx.tolist())

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size * 2, shuffle=False)
    test_loader = DataLoader(test_ds, batch_size=batch_size * 2, shuffle=False)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_features = features.shape[1]

    model = PhotoZMDN(
        n_features=n_features, n_hidden=128, n_components=5, dropout=0.1
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_val_loss = float('inf')
    best_state = None
    patience_counter = 0

    for epoch in range(1, epochs + 1):
        model.train()
        for feats, zs in train_loader:
            feats, zs = feats.to(device), zs.to(device)
            optimizer.zero_grad()
            pi, mu, sigma = model(feats)
            loss = PhotoZMDN.mdn_loss(pi, mu, sigma, zs)
            loss.backward()
            optimizer.step()
        scheduler.step()

        # Validate
        model.eval()
        val_loss = 0.0
        val_n = 0
        with torch.no_grad():
            for feats, zs in val_loader:
                feats, zs = feats.to(device), zs.to(device)
                pi, mu, sigma = model(feats)
                loss = PhotoZMDN.mdn_loss(pi, mu, sigma, zs)
                val_loss += loss.item() * feats.size(0)
                val_n += feats.size(0)
        val_loss /= val_n

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                break

    # Evaluate on test set
    model.load_state_dict(best_state)
    model.eval()

    all_z_true, all_z_pred = [], []
    with torch.no_grad():
        for feats, zs in test_loader:
            feats = feats.to(device)
            z_pt, _, _ = model.predict(feats)
            all_z_true.append(zs.numpy())
            all_z_pred.append(z_pt.cpu().numpy())

    z_true = np.concatenate(all_z_true)
    z_pred = np.concatenate(all_z_pred)

    dz = (z_pred - z_true) / (1.0 + z_true)
    bias = np.median(dz)
    sigma_mad = 1.4826 * np.median(np.abs(dz - bias))
    outlier_frac = np.mean(np.abs(dz) > 0.15)

    # Outlier analysis
    outlier_mask = np.abs(dz) > 0.15
    n_outliers = int(np.sum(outlier_mask))

    return {
        'label': label,
        'n_features': n_features,
        'sigma_mad': sigma_mad,
        'outlier_frac': outlier_frac,
        'bias': bias,
        'n_outliers': n_outliers,
        'n_test': len(z_true),
        'best_epoch': epoch - patience_counter,
    }


def main():
    h5_path = os.path.join(_project_root, 'photo_z', 'data', 'sdss_specphoto.h5')

    print(f"Loading data: {h5_path}")
    features, spec_z, feat_names, extra = load_data_with_extras(h5_path)
    print(f"  Samples: {len(spec_z)}, Features: {features.shape[1]}")
    print(f"  Extra keys: {list(extra.keys())}")

    # Check required columns
    if 'petroR90_r' not in extra:
        print("ERROR: petroR90_r not in extra. Re-run fetch_photoz_training.py first.")
        sys.exit(1)

    petroR50 = features[:, feat_names.index('petroR50_r')]  # already in features
    petroR90 = extra['petroR90_r']
    mag_r = extra['mag_r']

    # Clean petroR90 (SDSS uses -9999 for missing)
    bad_r90 = (petroR90 < 0) | ~np.isfinite(petroR90)
    petroR90[bad_r90] = petroR50[bad_r90] * 3.0  # default C~3
    print(f"  petroR90_r bad values replaced: {np.sum(bad_r90)}")

    # Compute new features
    # 1. Concentration index: C = R90/R50 (typical range 2-5)
    concentration = petroR90 / np.maximum(petroR50, 0.1)
    concentration = np.clip(concentration, 1.0, 10.0)

    # 2. Size-mag relation: log10(R50) - 0.2*mag_r
    #    This encodes "at given apparent magnitude, how large is the galaxy"
    #    Nearby galaxies are larger for same luminosity
    size_mag = np.log10(np.maximum(petroR50, 0.1)) - 0.2 * mag_r

    print(f"\n  Concentration C = R90/R50:")
    print(f"    range: [{concentration.min():.2f}, {concentration.max():.2f}]")
    print(f"    median: {np.median(concentration):.2f}")
    print(f"  Size-mag relation:")
    print(f"    range: [{size_mag.min():.2f}, {size_mag.max():.2f}]")
    print(f"    median: {np.median(size_mag):.2f}")

    # Experiments
    experiments = [
        ("Baseline (20 feat)", features),
        ("+ Concentration (C=R90/R50)",
         np.hstack([features, concentration.reshape(-1, 1)])),
        ("+ Size-Mag (log R50 - 0.2*mag_r)",
         np.hstack([features, size_mag.reshape(-1, 1)])),
        ("+ Both (C + Size-Mag)",
         np.hstack([features, concentration.reshape(-1, 1),
                    size_mag.reshape(-1, 1)])),
    ]

    results = []
    for label, feat in experiments:
        print(f"\n{'='*60}")
        print(f"Experiment: {label}")
        print(f"  Features: {feat.shape[1]}")
        t0 = time.time()
        res = train_and_evaluate(feat, spec_z, label)
        elapsed = time.time() - t0
        res['elapsed'] = elapsed
        results.append(res)
        print(f"  σ_MAD = {res['sigma_mad']:.4f}")
        print(f"  η = {res['outlier_frac']*100:.2f}% "
              f"({res['n_outliers']}/{res['n_test']})")
        print(f"  bias = {res['bias']:.5f}")
        print(f"  Time: {elapsed:.1f}s")

    # Summary table
    print(f"\n{'='*60}")
    print(f"{'='*60}")
    print(f"\n=== SUMMARY ===\n")
    print(f"{'Feature set':<35s}  {'N_feat':>6s}  {'σ_MAD':>7s}  "
          f"{'η(%)':>7s}  {'n_outlier':>9s}")
    print("-" * 75)

    baseline_eta = results[0]['outlier_frac']
    baseline_smad = results[0]['sigma_mad']

    for r in results:
        eta_pct = r['outlier_frac'] * 100
        delta_eta = (r['outlier_frac'] - baseline_eta) * 100
        delta_smad = (r['sigma_mad'] - baseline_smad)
        mark = ""
        if r != results[0]:
            mark = f"  (Δη={delta_eta:+.2f}%p, Δσ={delta_smad:+.4f})"
        print(f"{r['label']:<35s}  {r['n_features']:>6d}  "
              f"{r['sigma_mad']:.4f}  {eta_pct:>6.2f}%  "
              f"{r['n_outliers']:>4d}/{r['n_test']}{mark}")


if __name__ == '__main__':
    main()
