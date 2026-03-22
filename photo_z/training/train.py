"""Training loop for PhotoZ MDN model.

Uses SDSS spectro-photometric data from HDF5 file.
Adam + CosineAnnealingLR + early stopping on validation loss.
"""

import os
import sys
import time
import numpy as np
import torch
from torch.utils.data import DataLoader

from ..config import TrainConfig
from ..models.mlp_mdn import PhotoZMDN
from .data import create_datasets


def evaluate(model, loader, device):
    """
    Evaluate model on a data loader.

    Returns
    -------
    avg_loss : float
        Mean negative log-likelihood (MDN loss only).
    metrics : dict
        Dictionary with bias, sigma_MAD, outlier_fraction.
    """
    model.eval()
    total_loss = 0.0
    total_n = 0

    all_z_true = []
    all_z_pred = []
    all_z_err = []

    with torch.no_grad():
        for features, redshifts in loader:
            features = features.to(device)
            redshifts = redshifts.to(device)

            fwd = model(features)
            if model.has_quality_head:
                pi, mu, sigma, _ = fwd
            else:
                pi, mu, sigma = fwd
            loss = PhotoZMDN.mdn_loss(pi, mu, sigma, redshifts)
            total_loss += loss.item() * features.size(0)
            total_n += features.size(0)

            z_point, z_err, _, _ = model.predict(features)
            all_z_true.append(redshifts.cpu().numpy())
            all_z_pred.append(z_point.cpu().numpy())
            all_z_err.append(z_err.cpu().numpy())

    avg_loss = total_loss / total_n

    z_true = np.concatenate(all_z_true)
    z_pred = np.concatenate(all_z_pred)

    # Normalized residuals: delta_z / (1 + z_spec)
    dz = (z_pred - z_true) / (1.0 + z_true)

    # Bias: mean of normalized residuals
    bias = np.mean(dz)

    # sigma_MAD: 1.4826 * MAD(dz)
    sigma_mad = 1.4826 * np.median(np.abs(dz - np.median(dz)))

    # Outlier fraction: |dz| > 0.15
    outlier_frac = np.mean(np.abs(dz) > 0.15)

    metrics = {
        'bias': bias,
        'sigma_MAD': sigma_mad,
        'outlier_frac': outlier_frac,
    }

    return avg_loss, metrics


def train(config=None, data_path=None):
    """
    Main training loop with early stopping.

    Parameters
    ----------
    config : TrainConfig or None
        Training configuration.
    data_path : str or None
        Path to HDF5 training data.
    """
    if config is None:
        config = TrainConfig()

    if data_path is None:
        print("ERROR: --data path to HDF5 training data is required.",
              file=sys.stderr)
        sys.exit(1)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Set random seed
    torch.manual_seed(config.seed)
    np.random.seed(config.seed)

    # Load data
    print(f"Loading data from {data_path} ...")
    train_ds, val_ds, test_ds, metadata = create_datasets(data_path, config)
    n_features = metadata['n_features']
    print(f"Train: {len(train_ds)}, Val: {len(val_ds)}, Test: {len(test_ds)}")
    print(f"Features: {n_features}")

    train_loader = DataLoader(train_ds, batch_size=config.batch_size,
                              shuffle=True, num_workers=config.num_workers,
                              pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size * 2,
                            shuffle=False, num_workers=config.num_workers,
                            pin_memory=True)

    # Model
    use_quality = config.quality_head
    model = PhotoZMDN(
        n_features=n_features,
        n_hidden=128,
        n_components=5,
        dropout=0.1,
        quality_head=use_quality,
    ).to(device)
    print(f"Model parameters: {model.count_parameters():,}")
    if use_quality:
        print(f"Quality Head v2: ENABLED (λ={config.quality_lambda})")
        bce_loss_fn = torch.nn.BCEWithLogitsLoss()

    # Optimizer + scheduler
    optimizer = torch.optim.Adam(model.parameters(), lr=config.lr,
                                weight_decay=config.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.epochs)

    # Checkpoint dir
    os.makedirs(config.checkpoint_dir, exist_ok=True)

    best_val_loss = float('inf')
    patience_counter = 0

    for epoch in range(1, config.epochs + 1):
        t0 = time.time()

        # Train
        model.train()
        train_loss = 0.0
        train_total = 0

        for features, redshifts in train_loader:
            features = features.to(device)
            redshifts = redshifts.to(device)

            optimizer.zero_grad()
            fwd = model(features)

            if use_quality:
                pi, mu, sigma, q_logit = fwd
                loss_mdn = PhotoZMDN.mdn_loss(pi, mu, sigma, redshifts)
                # Quality label: is this prediction reliable?
                with torch.no_grad():
                    z_pred = (pi * mu).sum(dim=1)
                    dz = (z_pred - redshifts) / (1.0 + redshifts)
                    is_reliable = (torch.abs(dz) < 0.15).float()
                loss_q = bce_loss_fn(q_logit, is_reliable)
                loss = loss_mdn + config.quality_lambda * loss_q
            else:
                pi, mu, sigma = fwd
                loss = PhotoZMDN.mdn_loss(pi, mu, sigma, redshifts)

            loss.backward()
            optimizer.step()

            train_loss += loss.item() * features.size(0)
            train_total += features.size(0)

        scheduler.step()

        train_loss /= train_total

        # Validate
        val_loss, val_metrics = evaluate(model, val_loader, device)

        elapsed = time.time() - t0
        lr = optimizer.param_groups[0]['lr']
        print(f"Epoch {epoch:3d}/{config.epochs} | "
              f"train_loss={train_loss:.4f} | "
              f"val_loss={val_loss:.4f} "
              f"sigma_MAD={val_metrics['sigma_MAD']:.4f} "
              f"outlier={val_metrics['outlier_frac']:.4f} | "
              f"lr={lr:.6f} | {elapsed:.1f}s")

        # Early stopping on validation loss
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            ckpt = {
                'model_state': model.state_dict(),
                'epoch': epoch,
                'val_loss': val_loss,
                'val_metrics': val_metrics,
                'config': {
                    'n_features': n_features,
                    'n_hidden': 128,
                    'n_components': 5,
                    'dropout': 0.1,
                    'quality_head': use_quality,
                },
                'feat_mean': metadata['feat_mean'].tolist(),
                'feat_std': metadata['feat_std'].tolist(),
            }
            ckpt_path = config.checkpoint_path
            torch.save(ckpt, ckpt_path)
            print(f"  -> Saved best model (val_loss={val_loss:.4f})")
        else:
            patience_counter += 1
            if patience_counter >= config.early_stop_patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(best val_loss={best_val_loss:.4f})")
                break

    print(f"\nBest validation loss: {best_val_loss:.4f}")
    print(f"Checkpoint: {config.checkpoint_path}")

    # Evaluate on test set
    print("\n--- Test Set Evaluation ---")
    test_loader = DataLoader(test_ds, batch_size=config.batch_size * 2,
                             shuffle=False, num_workers=config.num_workers)

    ckpt = torch.load(config.checkpoint_path, map_location=device,
                      weights_only=False)
    model.load_state_dict(ckpt['model_state'])

    test_loss, test_metrics = evaluate(model, test_loader, device)
    print(f"Test NLL:       {test_loss:.4f}")
    print(f"Test bias:      {test_metrics['bias']:.6f}")
    print(f"Test sigma_MAD: {test_metrics['sigma_MAD']:.4f}")
    print(f"Test outlier:   {test_metrics['outlier_frac']:.4f}")

    return best_val_loss


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description='Train PhotoZ MDN model')
    parser.add_argument('--data', type=str, required=True,
                        help='Path to HDF5 training data')
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--patience', type=int, default=20)
    parser.add_argument('--quality-head', action='store_true',
                        help='Enable Quality Head v2 for outlier detection')
    parser.add_argument('--quality-lambda', type=float, default=1.0,
                        help='Weight for quality loss (default: 1.0)')
    args = parser.parse_args()

    cfg = TrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        early_stop_patience=args.patience,
        quality_head=args.quality_head,
        quality_lambda=args.quality_lambda,
    )
    train(config=cfg, data_path=args.data)
