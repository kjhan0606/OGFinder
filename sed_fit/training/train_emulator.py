"""Training script for the SPS Emulator MLP (forward model).

Trains the emulator to predict broadband magnitudes from physical
parameters using an HDF5 grid of SPS models.

Usage:
    cd OGFinder && python3 -m sed_fit.training.train_emulator
    cd OGFinder && python3 -m sed_fit.training.train_emulator --data sed_fit/data/sps_grid.h5
"""

import os
import sys
import time
import argparse
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset

from ..config import TrainConfig
from ..models.emulator import SPSEmulatorMLP


def load_grid(h5_path):
    """
    Load SPS grid from HDF5 file.

    Returns
    -------
    params : np.ndarray, shape (N, 6)
    mags : np.ndarray, shape (N, n_bands)
    meta : dict with band_names, wavelengths, etc.
    """
    import h5py
    with h5py.File(h5_path, 'r') as f:
        params = f['params'][:]
        mags = f['mags'][:]
        meta = {
            'band_names': list(f.attrs.get('band_names', [])),
            'n_bands': int(f.attrs.get('n_bands', mags.shape[1])),
            'n_samples': int(f.attrs.get('n_samples', params.shape[0])),
        }
        if 'wavelengths' in f.attrs:
            meta['wavelengths'] = f.attrs['wavelengths']
    return params, mags, meta


def train(data_path=None, cfg=None, verbose=True):
    """
    Train the SPS emulator (forward model).

    Parameters
    ----------
    data_path : str
        Path to HDF5 SPS grid.
    cfg : TrainConfig, optional
    verbose : bool

    Returns
    -------
    model : SPSEmulatorMLP
    history : dict
    """
    if cfg is None:
        cfg = TrainConfig()

    if data_path is None:
        data_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data", "sps_grid.h5"
        )

    if not os.path.exists(data_path):
        print(f"ERROR: SPS grid not found at {data_path}", file=sys.stderr)
        print("Generate it with:")
        print(f"  python3 -m sed_fit.grid.generate --output {data_path}")
        sys.exit(1)

    # Seed
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if verbose:
        print(f"Device: {device}")

    # Load data
    params, mags, meta = load_grid(data_path)
    n_bands = meta['n_bands']
    N = params.shape[0]
    if verbose:
        print(f"Loaded {N} SPS models with {n_bands} bands")

    # Train/val split
    n_val = int(N * cfg.val_fraction)
    n_train = N - n_val
    idx = np.random.permutation(N)
    train_idx = idx[:n_train]
    val_idx = idx[n_train:]

    # Compute normalization statistics on training set
    params_mean = params[train_idx].mean(axis=0)
    params_std = params[train_idx].std(axis=0)
    params_std = np.maximum(params_std, 1e-8)

    mags_mean = mags[train_idx].mean(axis=0)
    mags_std = mags[train_idx].std(axis=0)
    mags_std = np.maximum(mags_std, 1e-8)

    # Normalize
    params_norm = (params - params_mean) / params_std
    mags_norm = (mags - mags_mean) / mags_std

    # Datasets
    X_train = torch.from_numpy(params_norm[train_idx].astype(np.float32))
    Y_train = torch.from_numpy(mags_norm[train_idx].astype(np.float32))
    X_val = torch.from_numpy(params_norm[val_idx].astype(np.float32))
    Y_val = torch.from_numpy(mags_norm[val_idx].astype(np.float32))

    train_ds = TensorDataset(X_train, Y_train)
    val_ds = TensorDataset(X_val, Y_val)

    use_cuda = device.type == 'cuda'
    n_workers = cfg.num_workers if use_cuda else 0
    train_loader = DataLoader(train_ds, batch_size=cfg.batch_size,
                              shuffle=True, num_workers=n_workers,
                              pin_memory=use_cuda)
    val_loader = DataLoader(val_ds, batch_size=cfg.batch_size * 2,
                            shuffle=False, num_workers=n_workers,
                            pin_memory=use_cuda)

    if verbose:
        print(f"Train: {n_train}, Val: {n_val}")

    # Model
    model = SPSEmulatorMLP(
        n_phys_params=6,
        n_bands=n_bands,
        n_hidden=256,
    ).to(device)

    if verbose:
        print(f"Model parameters: {model.count_parameters():,}")

    # Loss, optimizer, scheduler
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(model.parameters(), lr=cfg.lr,
                                 weight_decay=cfg.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=cfg.epochs)

    # Training loop
    best_val_loss = float('inf')
    patience_counter = 0
    history = {'train_loss': [], 'val_loss': []}

    os.makedirs(cfg.checkpoint_dir, exist_ok=True)

    t0 = time.time()
    for epoch in range(1, cfg.epochs + 1):
        # Train
        model.train()
        train_loss = 0.0
        n_batches = 0
        for X_batch, Y_batch in train_loader:
            X_batch = X_batch.to(device, non_blocking=True)
            Y_batch = Y_batch.to(device, non_blocking=True)

            optimizer.zero_grad(set_to_none=True)
            pred = model(X_batch)
            loss = criterion(pred, Y_batch)
            loss.backward()
            optimizer.step()

            train_loss += loss.item()
            n_batches += 1

        train_loss /= max(n_batches, 1)
        scheduler.step()

        # Validate
        model.eval()
        val_loss = 0.0
        n_val_batches = 0
        with torch.no_grad():
            for X_batch, Y_batch in val_loader:
                X_batch = X_batch.to(device, non_blocking=True)
                Y_batch = Y_batch.to(device, non_blocking=True)
                pred = model(X_batch)
                loss = criterion(pred, Y_batch)
                val_loss += loss.item()
                n_val_batches += 1

        val_loss /= max(n_val_batches, 1)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)

        if verbose and (epoch % 10 == 0 or epoch == 1):
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch:3d}/{cfg.epochs}: "
                  f"train_loss={train_loss:.6f} val_loss={val_loss:.6f} "
                  f"lr={lr:.1e}")

        # Early stopping / checkpointing
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            checkpoint = {
                'model_state': model.state_dict(),
                'epoch': epoch,
                'val_loss': float(best_val_loss),
                'params_mean': params_mean.tolist(),
                'params_std': params_std.tolist(),
                'mags_mean': mags_mean.tolist(),
                'mags_std': mags_std.tolist(),
                'config': {
                    'n_phys_params': 6,
                    'n_bands': n_bands,
                    'n_hidden': 256,
                },
                'band_names': meta.get('band_names', []),
            }
            torch.save(checkpoint, cfg.checkpoint_emulator_path)
        else:
            patience_counter += 1
            if patience_counter >= cfg.early_stop_patience:
                if verbose:
                    print(f"Early stopping at epoch {epoch}")
                break

    elapsed = time.time() - t0
    if verbose:
        print(f"\nTraining complete in {elapsed:.1f}s")
        print(f"Best val loss: {best_val_loss:.6f}")
        print(f"Checkpoint: {cfg.checkpoint_emulator_path}")

    # Reload best model
    ckpt = torch.load(cfg.checkpoint_emulator_path, map_location=device,
                      weights_only=True)
    model.load_state_dict(ckpt['model_state'])

    # Final evaluation: compute residuals in mag space
    model.eval()
    with torch.no_grad():
        pred_val = model(X_val.to(device))
        pred_val = pred_val.cpu().numpy()

    # Un-normalize
    pred_mags = pred_val * mags_std + mags_mean
    true_mags = mags[val_idx]
    residuals = pred_mags - true_mags

    if verbose:
        print(f"\n=== Validation Residuals (mag) ===")
        print(f"  Mean:   {residuals.mean():.4f}")
        print(f"  Std:    {residuals.std():.4f}")
        print(f"  Median: {np.median(residuals):.4f}")
        print(f"  MAD:    {np.median(np.abs(residuals)):.4f}")
        if meta.get('band_names'):
            print(f"\n  Per-band RMS:")
            for i, band in enumerate(meta['band_names']):
                rms = np.sqrt(np.mean(residuals[:, i] ** 2))
                print(f"    {band:8s}: {rms:.4f} mag")

    return model, history


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Train SPS Emulator MLP (forward model)')
    parser.add_argument('--data', type=str, default=None,
                        help='Path to HDF5 SPS grid')
    parser.add_argument('--epochs', type=int, default=None)
    parser.add_argument('--batch-size', type=int, default=None)
    parser.add_argument('--lr', type=float, default=None)
    parser.add_argument('--patience', type=int, default=None)
    args = parser.parse_args()

    cfg = TrainConfig()
    if args.epochs is not None:
        cfg.epochs = args.epochs
    if args.batch_size is not None:
        cfg.batch_size = args.batch_size
    if args.lr is not None:
        cfg.lr = args.lr
    if args.patience is not None:
        cfg.early_stop_patience = args.patience

    train(data_path=args.data, cfg=cfg)


if __name__ == '__main__':
    main()
