#!/usr/bin/env python3
"""
Experiment: MDN with Quality Classification Head.

Idea: add a binary head that predicts P(reliable) — whether the
photo-z estimate for this galaxy will be an outlier or not.

Combined loss:
    L = L_MDN + λ * L_quality

where L_quality = BinaryCrossEntropy(predicted_reliable, actual_reliable)
and actual_reliable = (|Δz/(1+z)| < 0.15) computed from spec-z during training.

At inference, use P(reliable) instead of z_err for filtering.
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
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader, Subset


class PhotoZMDN_Quality(nn.Module):
    """MDN with added quality classification head."""

    def __init__(self, n_features, n_hidden=128, n_components=5, dropout=0.1):
        super().__init__()
        self.n_components = n_components

        # Shared backbone
        self.backbone = nn.Sequential(
            nn.Linear(n_features, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(n_hidden, n_hidden),
            nn.BatchNorm1d(n_hidden),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(n_hidden, 64),
            nn.BatchNorm1d(64),
            nn.ReLU(),
        )

        # MDN head
        self.pi_layer = nn.Linear(64, n_components)
        self.mu_layer = nn.Linear(64, n_components)
        self.sigma_layer = nn.Linear(64, n_components)

        # Quality head: P(reliable)
        self.quality_head = nn.Sequential(
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        h = self.backbone(x)

        # MDN outputs
        pi = torch.softmax(self.pi_layer(h), dim=-1)
        mu = self.mu_layer(h)
        sigma = torch.nn.functional.softplus(self.sigma_layer(h)) + 1e-4

        # Quality output
        q_logit = self.quality_head(h).squeeze(-1)

        return pi, mu, sigma, q_logit

    def predict(self, x):
        pi, mu, sigma, q_logit = self.forward(x)

        # Weighted mean as point estimate
        z_point = (pi * mu).sum(dim=-1)

        # Uncertainty from mixture
        z_err = torch.sqrt(
            (pi * (sigma**2 + mu**2)).sum(dim=-1) - z_point**2
        ).clamp(min=1e-4)

        # Quality probability
        p_reliable = torch.sigmoid(q_logit)

        return z_point, z_err, p_reliable

    @staticmethod
    def mdn_loss(pi, mu, sigma, target):
        target = target.unsqueeze(-1)
        log_prob = -0.5 * ((target - mu) / sigma)**2 \
                   - torch.log(sigma) - 0.5 * np.log(2 * np.pi)
        log_prob = log_prob + torch.log(pi + 1e-10)
        return -torch.logsumexp(log_prob, dim=-1).mean()


class PhotoZQualityDataset(Dataset):
    """Dataset that also provides quality labels."""

    def __init__(self, features, redshifts):
        self.features = torch.from_numpy(features.astype(np.float32))
        self.redshifts = torch.from_numpy(redshifts.astype(np.float32))

    def __len__(self):
        return len(self.redshifts)

    def __getitem__(self, idx):
        return self.features[idx], self.redshifts[idx]


def train_quality_model(features, spec_z, seed=42, epochs=100,
                        batch_size=256, lr=1e-3, patience=20,
                        lambda_q=1.0):
    """Train MDN+Quality model."""
    n = len(spec_z)
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)

    n_test = int(n * 0.2)
    n_val = int(n * 0.1)
    n_train = n - n_val - n_test

    train_idx = idx[:n_train]
    val_idx = idx[n_train:n_train + n_val]
    test_idx = idx[n_train + n_val:]

    # Normalize
    feat_mean = features[train_idx].mean(axis=0)
    feat_std = features[train_idx].std(axis=0)
    feat_std[feat_std < 1e-8] = 1.0
    features_norm = (features - feat_mean) / feat_std

    full_ds = PhotoZQualityDataset(features_norm, spec_z)
    train_ds = Subset(full_ds, train_idx.tolist())
    val_ds = Subset(full_ds, val_idx.tolist())
    test_ds = Subset(full_ds, test_idx.tolist())

    train_loader = DataLoader(train_ds, batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(val_ds, batch_size=batch_size * 2, shuffle=False)
    test_loader = DataLoader(test_ds, batch_size=batch_size * 2, shuffle=False)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    n_features = features.shape[1]

    model = PhotoZMDN_Quality(
        n_features=n_features, n_hidden=128, n_components=5, dropout=0.1
    ).to(device)

    n_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f"  Model params: {n_params:,}")

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)
    bce_loss_fn = nn.BCEWithLogitsLoss()

    # Outlier class weight (outliers are rare ~2%)
    # Use pos_weight to upweight outlier detection
    pos_weight = torch.tensor([1.0]).to(device)  # reliable is majority class

    best_val_loss = float('inf')
    best_state = None
    patience_counter = 0

    for epoch in range(1, epochs + 1):
        model.train()
        for feats, zs in train_loader:
            feats, zs = feats.to(device), zs.to(device)
            optimizer.zero_grad()

            pi, mu, sigma, q_logit = model(feats)

            # MDN loss
            loss_mdn = PhotoZMDN_Quality.mdn_loss(pi, mu, sigma, zs)

            # Quality loss: compute z_pred, check if outlier
            with torch.no_grad():
                z_pred = (pi * mu).sum(dim=-1)
                dz_norm = (z_pred - zs) / (1.0 + zs)
                is_reliable = (torch.abs(dz_norm) < 0.15).float()

            loss_quality = bce_loss_fn(q_logit, is_reliable)

            loss = loss_mdn + lambda_q * loss_quality
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
                pi, mu, sigma, q_logit = model(feats)
                l_mdn = PhotoZMDN_Quality.mdn_loss(pi, mu, sigma, zs)
                z_pred = (pi * mu).sum(dim=-1)
                dz = (z_pred - zs) / (1.0 + zs)
                is_rel = (torch.abs(dz) < 0.15).float()
                l_q = bce_loss_fn(q_logit, is_rel)
                val_loss += (l_mdn.item() + lambda_q * l_q.item()) * feats.size(0)
                val_n += feats.size(0)
        val_loss /= val_n

        if val_loss < best_val_loss:
            best_val_loss = val_loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}
            patience_counter = 0
        else:
            patience_counter += 1
            if patience_counter >= patience:
                print(f"  Early stop at epoch {epoch}")
                break

    # Evaluate on test set
    model.load_state_dict(best_state)
    model.eval()

    all_z_true, all_z_pred, all_z_err, all_p_rel = [], [], [], []
    with torch.no_grad():
        for feats, zs in test_loader:
            feats = feats.to(device)
            z_pt, z_er, p_rel = model.predict(feats)
            all_z_true.append(zs.numpy())
            all_z_pred.append(z_pt.cpu().numpy())
            all_z_err.append(z_er.cpu().numpy())
            all_p_rel.append(p_rel.cpu().numpy())

    z_true = np.concatenate(all_z_true)
    z_pred = np.concatenate(all_z_pred)
    z_err = np.concatenate(all_z_err)
    p_reliable = np.concatenate(all_p_rel)

    return z_true, z_pred, z_err, p_reliable


def compute_filter_metrics(z_true, z_pred, mask, label):
    """Compute metrics for a filtered subset."""
    zt = z_true[mask]
    zp = z_pred[mask]
    dz = (zp - zt) / (1.0 + zt)
    bias = np.median(dz)
    sigma_mad = 1.4826 * np.median(np.abs(dz - bias))
    outlier = np.mean(np.abs(dz) > 0.15)
    n_out = int(np.sum(np.abs(dz) > 0.15))
    n_total = len(zt)
    completeness = n_total / len(z_true) * 100
    return {
        'label': label,
        'n': n_total,
        'completeness': completeness,
        'sigma_mad': sigma_mad,
        'outlier_frac': outlier,
        'n_outliers': n_out,
        'bias': bias,
    }


def main():
    import h5py

    h5_path = os.path.join(_project_root, 'photo_z', 'data', 'sdss_specphoto.h5')
    print(f"Loading data: {h5_path}")
    with h5py.File(h5_path, 'r') as f:
        features = f['features'][:]
        spec_z = f['redshift'][:]
    print(f"  Samples: {len(spec_z)}, Features: {features.shape[1]}")

    # Train MDN+Quality model
    print(f"\n{'='*60}")
    print("Training MDN + Quality Head (λ=1.0)")
    print('='*60)
    t0 = time.time()
    z_true, z_pred, z_err, p_reliable = train_quality_model(
        features, spec_z, lambda_q=1.0)
    elapsed = time.time() - t0
    print(f"  Training time: {elapsed:.1f}s")

    dz_norm = (z_pred - z_true) / (1.0 + z_true)
    is_outlier = np.abs(dz_norm) > 0.15

    # Overall metrics
    n_total = len(z_true)
    n_outlier = int(np.sum(is_outlier))
    print(f"\n  Test set: {n_total}, Outliers: {n_outlier} ({n_outlier/n_total*100:.1f}%)")

    # Distribution of p_reliable
    print(f"\n--- P(reliable) Distribution ---")
    print(f"  Outliers:  median={np.median(p_reliable[is_outlier]):.3f}  "
          f"mean={np.mean(p_reliable[is_outlier]):.3f}")
    print(f"  Normal:    median={np.median(p_reliable[~is_outlier]):.3f}  "
          f"mean={np.mean(p_reliable[~is_outlier]):.3f}")

    # Compare filtering: z_err vs p_reliable
    print(f"\n{'='*60}")
    print("=== Comparison: z_err filter vs P(reliable) filter ===")
    print('='*60)

    print(f"\n--- Method 1: z_err threshold ---")
    print(f"{'Filter':<20s}  {'N':>5s}  {'Complete':>8s}  {'σ_MAD':>7s}  "
          f"{'η(%)':>7s}  {'n_out':>5s}")
    print("-" * 65)

    for thresh in [1.0, 0.20, 0.15, 0.10, 0.08, 0.06, 0.05]:
        mask = z_err < thresh if thresh < 1.0 else np.ones(n_total, dtype=bool)
        label = f"z_err < {thresh}" if thresh < 1.0 else "All (no filter)"
        r = compute_filter_metrics(z_true, z_pred, mask, label)
        print(f"{label:<20s}  {r['n']:>5d}  {r['completeness']:>7.1f}%  "
              f"{r['sigma_mad']:.4f}  {r['outlier_frac']*100:>6.2f}%  "
              f"{r['n_outliers']:>5d}")

    print(f"\n--- Method 2: P(reliable) threshold ---")
    print(f"{'Filter':<20s}  {'N':>5s}  {'Complete':>8s}  {'σ_MAD':>7s}  "
          f"{'η(%)':>7s}  {'n_out':>5s}")
    print("-" * 65)

    for thresh in [0.0, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95]:
        mask = p_reliable > thresh if thresh > 0 else np.ones(n_total, dtype=bool)
        label = f"P(rel) > {thresh}" if thresh > 0 else "All (no filter)"
        r = compute_filter_metrics(z_true, z_pred, mask, label)
        print(f"{label:<20s}  {r['n']:>5d}  {r['completeness']:>7.1f}%  "
              f"{r['sigma_mad']:.4f}  {r['outlier_frac']*100:>6.2f}%  "
              f"{r['n_outliers']:>5d}")

    # Head-to-head: at similar completeness levels
    print(f"\n{'='*60}")
    print("=== Head-to-head at similar completeness ===")
    print('='*60)
    print(f"\n{'Completeness':>12s}  {'z_err σ_MAD':>12s}  {'z_err η':>8s}  "
          f"{'P(rel) σ_MAD':>13s}  {'P(rel) η':>9s}")
    print("-" * 65)

    target_comps = [95, 90, 85, 80, 75]
    for target in target_comps:
        # Find z_err threshold for target completeness
        z_err_sorted = np.sort(z_err)
        idx_ze = min(int(target / 100.0 * n_total), n_total - 1)
        ze_thresh = z_err_sorted[idx_ze]
        mask_ze = z_err <= ze_thresh
        r_ze = compute_filter_metrics(z_true, z_pred, mask_ze, "")

        # Find p_reliable threshold for target completeness
        p_sorted = np.sort(p_reliable)[::-1]
        idx_pr = min(int(target / 100.0 * n_total), n_total - 1)
        pr_thresh = p_sorted[idx_pr]
        mask_pr = p_reliable >= pr_thresh
        r_pr = compute_filter_metrics(z_true, z_pred, mask_pr, "")

        print(f"{target:>10d}%  {r_ze['sigma_mad']:>11.4f}  "
              f"{r_ze['outlier_frac']*100:>7.2f}%  "
              f"{r_pr['sigma_mad']:>12.4f}  "
              f"{r_pr['outlier_frac']*100:>8.2f}%")

    # ROC-like analysis: outlier detection
    print(f"\n{'='*60}")
    print("=== Outlier Detection (ROC-like) ===")
    print('='*60)

    # For z_err: higher z_err → more likely outlier
    # For p_reliable: lower p_reliable → more likely outlier
    from sklearn.metrics import roc_auc_score, average_precision_score
    auc_zerr = roc_auc_score(is_outlier, z_err)
    auc_prel = roc_auc_score(is_outlier, 1 - p_reliable)
    ap_zerr = average_precision_score(is_outlier, z_err)
    ap_prel = average_precision_score(is_outlier, 1 - p_reliable)

    print(f"\n  Outlier detection as binary classification:")
    print(f"  {'Metric':<25s}  {'z_err':>10s}  {'P(reliable)':>12s}")
    print(f"  {'ROC AUC':<25s}  {auc_zerr:>10.4f}  {auc_prel:>12.4f}")
    print(f"  {'Average Precision':<25s}  {ap_zerr:>10.4f}  {ap_prel:>12.4f}")


if __name__ == '__main__':
    main()
