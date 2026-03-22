#!/usr/bin/env python3
"""
Comprehensive uncertainty experiment for photo-z outlier detection.

Methods tested:
1. z_err (baseline) — MDN mixture std
2. Quality Head (v1) — binary head on shared backbone
3. Quality Head (v2) — binary head with MDN internals (π,μ,σ) as extra input
4. Entropy — entropy of the mixture P(z) distribution
5. Multimodality score — separation between GMM modes
6. Ensemble disagreement — variance across multiple trained models
7. Combined score — weighted combination of all indicators

Goal: maximize outlier detection (AUROC, AP) while preserving completeness.
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
from sklearn.metrics import roc_auc_score, average_precision_score


# ============================================================
# Models
# ============================================================

class PhotoZMDN_Base(nn.Module):
    """Standard MDN (baseline)."""

    def __init__(self, n_features, n_hidden=128, n_components=5, dropout=0.1):
        super().__init__()
        self.n_components = n_components
        self.backbone = nn.Sequential(
            nn.Linear(n_features, n_hidden),
            nn.BatchNorm1d(n_hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(n_hidden, n_hidden),
            nn.BatchNorm1d(n_hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(n_hidden, 64),
            nn.BatchNorm1d(64), nn.ReLU(),
        )
        self.pi_layer = nn.Linear(64, n_components)
        self.mu_layer = nn.Linear(64, n_components)
        self.sigma_layer = nn.Linear(64, n_components)

    def forward(self, x):
        h = self.backbone(x)
        pi = torch.softmax(self.pi_layer(h), dim=-1)
        mu = self.mu_layer(h)
        sigma = nn.functional.softplus(self.sigma_layer(h)) + 1e-4
        return pi, mu, sigma

    def predict(self, x):
        pi, mu, sigma = self.forward(x)
        z_point = (pi * mu).sum(dim=-1)
        z_err = torch.sqrt(
            (pi * (sigma**2 + mu**2)).sum(dim=-1) - z_point**2
        ).clamp(min=1e-4)
        return z_point, z_err, pi, mu, sigma

    @staticmethod
    def mdn_loss(pi, mu, sigma, target):
        target = target.unsqueeze(-1)
        log_prob = -0.5 * ((target - mu) / sigma)**2 \
                   - torch.log(sigma) - 0.5 * np.log(2 * np.pi)
        log_prob = log_prob + torch.log(pi + 1e-10)
        return -torch.logsumexp(log_prob, dim=-1).mean()


class PhotoZMDN_QualityV2(nn.Module):
    """MDN + Quality head that also sees MDN internals (π,μ,σ)."""

    def __init__(self, n_features, n_hidden=128, n_components=5, dropout=0.1):
        super().__init__()
        self.n_components = n_components
        self.backbone = nn.Sequential(
            nn.Linear(n_features, n_hidden),
            nn.BatchNorm1d(n_hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(n_hidden, n_hidden),
            nn.BatchNorm1d(n_hidden), nn.ReLU(), nn.Dropout(dropout),
            nn.Linear(n_hidden, 64),
            nn.BatchNorm1d(64), nn.ReLU(),
        )
        self.pi_layer = nn.Linear(64, n_components)
        self.mu_layer = nn.Linear(64, n_components)
        self.sigma_layer = nn.Linear(64, n_components)

        # Quality head v2: backbone features (64) + MDN outputs (3*K)
        q_input_dim = 64 + 3 * n_components
        self.quality_head = nn.Sequential(
            nn.Linear(q_input_dim, 64),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(64, 32),
            nn.ReLU(),
            nn.Linear(32, 1),
        )

    def forward(self, x):
        h = self.backbone(x)
        pi = torch.softmax(self.pi_layer(h), dim=-1)
        mu = self.mu_layer(h)
        sigma = nn.functional.softplus(self.sigma_layer(h)) + 1e-4

        # Quality head gets backbone + MDN outputs
        q_input = torch.cat([h, pi, mu, sigma], dim=-1)
        q_logit = self.quality_head(q_input).squeeze(-1)

        return pi, mu, sigma, q_logit

    def predict(self, x):
        pi, mu, sigma, q_logit = self.forward(x)
        z_point = (pi * mu).sum(dim=-1)
        z_err = torch.sqrt(
            (pi * (sigma**2 + mu**2)).sum(dim=-1) - z_point**2
        ).clamp(min=1e-4)
        p_reliable = torch.sigmoid(q_logit)
        return z_point, z_err, p_reliable, pi, mu, sigma

    @staticmethod
    def mdn_loss(pi, mu, sigma, target):
        target = target.unsqueeze(-1)
        log_prob = -0.5 * ((target - mu) / sigma)**2 \
                   - torch.log(sigma) - 0.5 * np.log(2 * np.pi)
        log_prob = log_prob + torch.log(pi + 1e-10)
        return -torch.logsumexp(log_prob, dim=-1).mean()


# ============================================================
# Uncertainty metrics (computed from MDN outputs)
# ============================================================

def compute_entropy(pi, mu, sigma, z_grid=None):
    """Compute entropy of the mixture distribution P(z).
    Higher entropy = more uncertain / spread out."""
    if z_grid is None:
        z_grid = np.linspace(0, 1.5, 300)

    n = pi.shape[0]
    K = pi.shape[1]
    entropy = np.zeros(n)

    for i in range(n):
        # Evaluate mixture PDF on grid
        pdf = np.zeros(len(z_grid))
        for k in range(K):
            pdf += pi[i, k] * np.exp(-0.5 * ((z_grid - mu[i, k]) / sigma[i, k])**2) \
                   / (sigma[i, k] * np.sqrt(2 * np.pi))
        pdf = np.maximum(pdf, 1e-30)
        dz = z_grid[1] - z_grid[0]
        pdf /= (pdf.sum() * dz)  # normalize
        entropy[i] = -np.sum(pdf * np.log(pdf + 1e-30)) * dz

    return entropy


def compute_multimodality(pi, mu, sigma):
    """Multimodality score: how bimodal is the mixture?
    Uses separation between top-2 weighted modes relative to their widths."""
    n = pi.shape[0]
    score = np.zeros(n)

    for i in range(n):
        # Sort components by weight
        order = np.argsort(pi[i])[::-1]
        if pi[i, order[0]] > 0.95:
            # Single dominant mode
            score[i] = 0.0
            continue

        # Separation between top-2 modes
        k1, k2 = order[0], order[1]
        sep = abs(mu[i, k1] - mu[i, k2])
        avg_sigma = 0.5 * (sigma[i, k1] + sigma[i, k2])
        w2 = pi[i, k2]

        # Score: large separation + significant 2nd component = high multimodality
        score[i] = (sep / max(avg_sigma, 0.01)) * min(w2 / pi[i, k1], 1.0)

    return score


def compute_pit(z_true, pi, mu, sigma):
    """PIT (Probability Integral Transform) values.
    Well-calibrated model → PIT ~ Uniform(0,1).
    Extreme PIT values (near 0 or 1) indicate outliers."""
    n = pi.shape[0]
    K = pi.shape[1]
    pit = np.zeros(n)

    for i in range(n):
        # CDF of mixture at z_true
        cdf = 0.0
        for k in range(K):
            from scipy.stats import norm
            cdf += pi[i, k] * norm.cdf(z_true[i], loc=mu[i, k], scale=sigma[i, k])
        pit[i] = cdf

    return pit


# ============================================================
# Training
# ============================================================

class SimpleDataset(Dataset):
    def __init__(self, features, redshifts):
        self.features = torch.from_numpy(features.astype(np.float32))
        self.redshifts = torch.from_numpy(redshifts.astype(np.float32))
    def __len__(self):
        return len(self.redshifts)
    def __getitem__(self, idx):
        return self.features[idx], self.redshifts[idx]


def split_data(features, spec_z, seed=42):
    """Split into train/val/test."""
    n = len(spec_z)
    rng = np.random.RandomState(seed)
    idx = rng.permutation(n)
    n_test = int(n * 0.2)
    n_val = int(n * 0.1)
    n_train = n - n_val - n_test
    return idx[:n_train], idx[n_train:n_train+n_val], idx[n_train+n_val:]


def train_base_mdn(features, spec_z, seed=42, epochs=100, lr=1e-3, patience=20):
    """Train standard MDN, return model + test predictions."""
    train_idx, val_idx, test_idx = split_data(features, spec_z, seed)

    feat_mean = features[train_idx].mean(axis=0)
    feat_std = features[train_idx].std(axis=0)
    feat_std[feat_std < 1e-8] = 1.0
    features_norm = (features - feat_mean) / feat_std

    ds = SimpleDataset(features_norm, spec_z)
    train_loader = DataLoader(Subset(ds, train_idx.tolist()), batch_size=256, shuffle=True)
    val_loader = DataLoader(Subset(ds, val_idx.tolist()), batch_size=512, shuffle=False)
    test_loader = DataLoader(Subset(ds, test_idx.tolist()), batch_size=512, shuffle=False)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = PhotoZMDN_Base(n_features=features.shape[1]).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_val_loss = float('inf')
    best_state = None
    wait = 0

    for ep in range(1, epochs + 1):
        model.train()
        for f, z in train_loader:
            f, z = f.to(device), z.to(device)
            optimizer.zero_grad()
            pi, mu, sigma = model(f)
            PhotoZMDN_Base.mdn_loss(pi, mu, sigma, z).backward()
            optimizer.step()
        scheduler.step()

        model.eval()
        vl = 0; vn = 0
        with torch.no_grad():
            for f, z in val_loader:
                f, z = f.to(device), z.to(device)
                pi, mu, sigma = model(f)
                vl += PhotoZMDN_Base.mdn_loss(pi, mu, sigma, z).item() * f.size(0)
                vn += f.size(0)
        vl /= vn

        if vl < best_val_loss:
            best_val_loss = vl; best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}; wait = 0
        else:
            wait += 1
            if wait >= patience: break

    model.load_state_dict(best_state)
    model.eval()

    all_zt, all_zp, all_ze, all_pi, all_mu, all_sig = [], [], [], [], [], []
    with torch.no_grad():
        for f, z in test_loader:
            f = f.to(device)
            zp, ze, pi, mu, sig = model.predict(f)
            all_zt.append(z.numpy())
            all_zp.append(zp.cpu().numpy())
            all_ze.append(ze.cpu().numpy())
            all_pi.append(pi.cpu().numpy())
            all_mu.append(mu.cpu().numpy())
            all_sig.append(sig.cpu().numpy())

    return (np.concatenate(all_zt), np.concatenate(all_zp), np.concatenate(all_ze),
            np.concatenate(all_pi), np.concatenate(all_mu), np.concatenate(all_sig),
            model, feat_mean, feat_std)


def train_quality_v2(features, spec_z, seed=42, epochs=100, lr=1e-3,
                     patience=20, lambda_q=1.0):
    """Train MDN + Quality Head v2 (with MDN internals)."""
    train_idx, val_idx, test_idx = split_data(features, spec_z, seed)

    feat_mean = features[train_idx].mean(axis=0)
    feat_std = features[train_idx].std(axis=0)
    feat_std[feat_std < 1e-8] = 1.0
    features_norm = (features - feat_mean) / feat_std

    ds = SimpleDataset(features_norm, spec_z)
    train_loader = DataLoader(Subset(ds, train_idx.tolist()), batch_size=256, shuffle=True)
    val_loader = DataLoader(Subset(ds, val_idx.tolist()), batch_size=512, shuffle=False)
    test_loader = DataLoader(Subset(ds, test_idx.tolist()), batch_size=512, shuffle=False)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model = PhotoZMDN_QualityV2(n_features=features.shape[1]).to(device)
    bce = nn.BCEWithLogitsLoss()

    optimizer = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

    best_val_loss = float('inf')
    best_state = None
    wait = 0

    for ep in range(1, epochs + 1):
        model.train()
        for f, z in train_loader:
            f, z = f.to(device), z.to(device)
            optimizer.zero_grad()
            pi, mu, sigma, q_logit = model(f)
            l_mdn = PhotoZMDN_QualityV2.mdn_loss(pi, mu, sigma, z)
            with torch.no_grad():
                zp = (pi * mu).sum(-1)
                dz = (zp - z) / (1.0 + z)
                is_rel = (torch.abs(dz) < 0.15).float()
            l_q = bce(q_logit, is_rel)
            (l_mdn + lambda_q * l_q).backward()
            optimizer.step()
        scheduler.step()

        model.eval()
        vl = 0; vn = 0
        with torch.no_grad():
            for f, z in val_loader:
                f, z = f.to(device), z.to(device)
                pi, mu, sigma, q_logit = model(f)
                l_mdn = PhotoZMDN_QualityV2.mdn_loss(pi, mu, sigma, z)
                zp = (pi * mu).sum(-1)
                dz = (zp - z) / (1.0 + z)
                is_rel = (torch.abs(dz) < 0.15).float()
                l_q = bce(q_logit, is_rel)
                vl += (l_mdn.item() + lambda_q * l_q.item()) * f.size(0)
                vn += f.size(0)
        vl /= vn

        if vl < best_val_loss:
            best_val_loss = vl; best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}; wait = 0
        else:
            wait += 1
            if wait >= patience: break

    model.load_state_dict(best_state)
    model.eval()

    all_zt, all_zp, all_ze, all_prel = [], [], [], []
    all_pi, all_mu, all_sig = [], [], []
    with torch.no_grad():
        for f, z in test_loader:
            f = f.to(device)
            zp, ze, prel, pi, mu, sig = model.predict(f)
            all_zt.append(z.numpy())
            all_zp.append(zp.cpu().numpy())
            all_ze.append(ze.cpu().numpy())
            all_prel.append(prel.cpu().numpy())
            all_pi.append(pi.cpu().numpy())
            all_mu.append(mu.cpu().numpy())
            all_sig.append(sig.cpu().numpy())

    return (np.concatenate(all_zt), np.concatenate(all_zp), np.concatenate(all_ze),
            np.concatenate(all_prel), np.concatenate(all_pi),
            np.concatenate(all_mu), np.concatenate(all_sig))


# ============================================================
# Main
# ============================================================

def main():
    import h5py

    h5_path = os.path.join(_project_root, 'photo_z', 'data', 'sdss_specphoto.h5')
    print(f"Loading: {h5_path}")
    with h5py.File(h5_path, 'r') as f:
        features = f['features'][:]
        spec_z = f['redshift'][:]
    print(f"  N={len(spec_z)}, Features={features.shape[1]}")

    # ========================================
    # 1. Train base MDN (for z_err, entropy, multimodality)
    # ========================================
    print(f"\n{'='*70}")
    print("Step 1: Training base MDN")
    print('='*70)
    t0 = time.time()
    (z_true, z_pred, z_err, pi, mu, sigma,
     base_model, feat_mean, feat_std) = train_base_mdn(features, spec_z)
    print(f"  Time: {time.time()-t0:.1f}s")

    dz_norm = (z_pred - z_true) / (1.0 + z_true)
    is_outlier = np.abs(dz_norm) > 0.15
    n_total = len(z_true)
    n_outlier = int(np.sum(is_outlier))
    print(f"  Test: {n_total}, Outliers: {n_outlier} ({n_outlier/n_total*100:.1f}%)")

    # Compute uncertainty metrics from base MDN
    print("\n  Computing entropy...")
    entropy = compute_entropy(pi, mu, sigma)

    print("  Computing multimodality score...")
    multimod = compute_multimodality(pi, mu, sigma)

    # Additional derived metrics
    # Max component weight (low = uncertain)
    pi_max = np.max(pi, axis=1)
    # Std of component means (high = multimodal)
    mu_std = np.std(mu, axis=1)
    # Weighted mean sigma
    weighted_sigma = np.sum(pi * sigma, axis=1)

    # ========================================
    # 2. Train MDN + Quality Head v2
    # ========================================
    print(f"\n{'='*70}")
    print("Step 2: Training MDN + Quality Head v2 (with MDN internals)")
    print('='*70)
    t0 = time.time()
    (z_true_q, z_pred_q, z_err_q, p_reliable_q,
     pi_q, mu_q, sigma_q) = train_quality_v2(features, spec_z)
    print(f"  Time: {time.time()-t0:.1f}s")

    dz_norm_q = (z_pred_q - z_true_q) / (1.0 + z_true_q)
    is_outlier_q = np.abs(dz_norm_q) > 0.15

    # ========================================
    # 3. Ensemble (train 5 models with different seeds)
    # ========================================
    print(f"\n{'='*70}")
    print("Step 3: Training ensemble (5 models)")
    print('='*70)
    ensemble_preds = []
    t0 = time.time()
    for i, seed in enumerate([42, 123, 7, 999, 2024]):
        print(f"  Model {i+1}/5 (seed={seed})...", end=" ", flush=True)
        zt, zp, ze, _, _, _, _, _, _ = train_base_mdn(features, spec_z, seed=seed)
        ensemble_preds.append(zp)
        print("done")
    print(f"  Total time: {time.time()-t0:.1f}s")

    # Ensemble disagreement = std of predictions across models
    ensemble_preds = np.array(ensemble_preds)  # shape (5, n_test)
    ensemble_std = np.std(ensemble_preds, axis=0)
    ensemble_mean = np.mean(ensemble_preds, axis=0)
    # Use same z_true from first model (same test split with seed=42)

    # ========================================
    # 4. Combined score (logistic regression on all indicators)
    # ========================================
    print(f"\n{'='*70}")
    print("Step 4: Training combined classifier")
    print('='*70)

    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler

    # Features for the meta-classifier (all from base MDN)
    meta_features = np.column_stack([
        z_err,
        entropy,
        multimod,
        pi_max,
        mu_std,
        weighted_sigma,
        ensemble_std,
    ])
    meta_names = ['z_err', 'entropy', 'multimod', 'pi_max', 'mu_std',
                  'weighted_sigma', 'ensemble_std']

    # Train/test split for meta-classifier (use first half to train, second half to eval)
    n_meta = len(z_true)
    mid = n_meta // 2

    scaler = StandardScaler()
    X_train_meta = scaler.fit_transform(meta_features[:mid])
    X_test_meta = scaler.transform(meta_features[mid:])
    y_train_meta = is_outlier[:mid].astype(int)
    y_test_meta = is_outlier[mid:].astype(int)

    clf = LogisticRegression(class_weight='balanced', max_iter=1000)
    clf.fit(X_train_meta, y_train_meta)
    combined_score_test = clf.predict_proba(X_test_meta)[:, 1]  # P(outlier)

    # For fair comparison, also compute individual scores on test half
    z_err_test = z_err[mid:]
    entropy_test = entropy[mid:]
    multimod_test = multimod[mid:]
    ensemble_std_test = ensemble_std[mid:]
    is_outlier_test = is_outlier[mid:]

    print(f"  Meta-classifier trained on {mid} samples")
    print(f"  Feature importances (logistic regression coefficients):")
    for name, coef in zip(meta_names, clf.coef_[0]):
        print(f"    {name:<20s}: {coef:+.4f}")

    # ========================================
    # 5. Comprehensive comparison
    # ========================================
    print(f"\n{'='*70}")
    print("="*70)
    print("=== COMPREHENSIVE COMPARISON: Outlier Detection ===")
    print("="*70)

    # All methods — compute AUROC and AP
    methods = {
        'z_err':             (z_err, is_outlier, False),      # higher = more outlier
        'Entropy':           (entropy, is_outlier, False),
        'Multimodality':     (multimod, is_outlier, False),
        '1 - π_max':         (1 - pi_max, is_outlier, False),
        'μ_std':             (mu_std, is_outlier, False),
        'Weighted σ':        (weighted_sigma, is_outlier, False),
        'Ensemble std':      (ensemble_std, is_outlier, False),
        'Quality v2 (1-P)':  (1 - p_reliable_q, is_outlier_q, False),
    }

    print(f"\n{'Method':<22s}  {'ROC AUC':>8s}  {'Avg Prec':>9s}  "
          f"{'Outlier median':>14s}  {'Normal median':>14s}")
    print("-" * 78)

    scores_for_combine = {}
    for name, (score, outlier_flag, _) in methods.items():
        try:
            auc = roc_auc_score(outlier_flag, score)
            ap = average_precision_score(outlier_flag, score)
        except:
            auc = ap = 0.0

        med_out = np.median(score[outlier_flag]) if np.sum(outlier_flag) > 0 else 0
        med_norm = np.median(score[~outlier_flag])

        print(f"{name:<22s}  {auc:>8.4f}  {ap:>9.4f}  "
              f"{med_out:>14.4f}  {med_norm:>14.4f}")
        scores_for_combine[name] = score

    # Combined (only on test half)
    try:
        auc_comb = roc_auc_score(is_outlier_test, combined_score_test)
        ap_comb = average_precision_score(is_outlier_test, combined_score_test)
        med_out_c = np.median(combined_score_test[is_outlier_test]) if np.sum(is_outlier_test) > 0 else 0
        med_norm_c = np.median(combined_score_test[~is_outlier_test])
        print(f"{'Combined (LR)*':<22s}  {auc_comb:>8.4f}  {ap_comb:>9.4f}  "
              f"{med_out_c:>14.4f}  {med_norm_c:>14.4f}")
    except:
        pass
    print("\n  * Combined evaluated on held-out half of test set")

    # ========================================
    # 6. Completeness vs Purity comparison (best methods)
    # ========================================
    print(f"\n{'='*70}")
    print("=== Completeness vs Outlier Rate (top methods) ===")
    print("="*70)

    best_methods = {
        'z_err':            z_err,
        'Entropy':          entropy,
        'Ensemble std':     ensemble_std,
        'Quality v2 (1-P)': 1 - p_reliable_q,
    }

    target_comps = [98, 95, 90, 85, 80, 75]

    for comp_target in target_comps:
        print(f"\n  Completeness = {comp_target}%:")
        print(f"    {'Method':<22s}  {'σ_MAD':>7s}  {'η(%)':>7s}  {'n_out':>5s}")

        for mname, mscore in best_methods.items():
            # Use corresponding z_pred and z_true
            if mname == 'Quality v2 (1-P)':
                zt, zp = z_true_q, z_pred_q
                out_flag = is_outlier_q
            else:
                zt, zp = z_true, z_pred
                out_flag = is_outlier

            n = len(zt)
            k = int(comp_target / 100.0 * n)
            thresh = np.sort(mscore)[k] if k < n else mscore.max() + 1
            mask = mscore <= thresh

            dz = (zp[mask] - zt[mask]) / (1.0 + zt[mask])
            bias = np.median(dz)
            smad = 1.4826 * np.median(np.abs(dz - bias))
            eta = np.mean(np.abs(dz) > 0.15) * 100
            nout = int(np.sum(np.abs(dz) > 0.15))

            print(f"    {mname:<22s}  {smad:.4f}  {eta:>6.2f}%  {nout:>5d}")

    # ========================================
    # 7. Correlation between methods
    # ========================================
    print(f"\n{'='*70}")
    print("=== Correlation between uncertainty indicators ===")
    print("="*70)

    corr_methods = ['z_err', 'Entropy', 'Multimodality', 'Ensemble std']
    corr_data = [z_err, entropy, multimod, ensemble_std]

    print(f"\n    {'':>14s}", end="")
    for name in corr_methods:
        print(f"  {name:>12s}", end="")
    print()

    for i, (name_i, data_i) in enumerate(zip(corr_methods, corr_data)):
        print(f"    {name_i:>14s}", end="")
        for j, (name_j, data_j) in enumerate(zip(corr_methods, corr_data)):
            r = np.corrcoef(data_i, data_j)[0, 1]
            print(f"  {r:>12.3f}", end="")
        print()

    # ========================================
    # Summary
    # ========================================
    print(f"\n{'='*70}")
    print("=== FINAL SUMMARY ===")
    print("="*70)
    print("""
    Best single indicators (by AUROC):
    - Check ranking above

    Recommendations:
    1. For simple use: z_err threshold (no extra computation)
    2. For better detection: Quality Head v2 (uses MDN internals)
    3. For best detection: Combined score (ensemble + all indicators)
    4. For production: ensemble_std (robust, model-agnostic)
    """)


if __name__ == '__main__':
    main()
