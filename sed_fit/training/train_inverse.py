"""Training script for the Inverse SED MLP+MDN model.

Trains the inverse model to predict physical parameters from observed
photometry + photo-z, using a Mixture Density Network for uncertainty
estimation.

Usage:
    cd OGFinder && python3 -m sed_fit.training.train_inverse
    cd OGFinder && python3 -m sed_fit.training.train_inverse --data sed_fit/data/sps_grid.h5
"""

import os
import sys
import time
import argparse
import numpy as np
import torch
from torch.utils.data import DataLoader, TensorDataset

from ..config import TrainConfig, SEDConfig, OUTPUT_PARAM_NAMES
from ..models.inverse import InverseSEDMLP
from ..features import normalize_magnitudes


def load_grid(h5_path):
    """
    Load SPS grid from HDF5 file.

    Returns
    -------
    params : np.ndarray, shape (N, 6)
        Physical parameters: z, log_mass, log_age, log_Z, Av, log_tau.
    mags : np.ndarray, shape (N, n_bands)
        Synthetic magnitudes.
    meta : dict
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
    return params, mags, meta


def prepare_training_data(params, mags, val_fraction=0.15, seed=42):
    """
    Prepare input/target pairs for the inverse model.

    Input:  normalized magnitudes + redshift  (n_bands + 1)
    Target: log_mass, log_age, log_Z, Av, log_tau  (5 params)

    Parameters
    ----------
    params : np.ndarray, shape (N, 6)
    mags : np.ndarray, shape (N, n_bands)
    val_fraction : float
    seed : int

    Returns
    -------
    train_X, train_Y, val_X, val_Y : np.ndarray
    norm_stats : dict with input normalization stats
    target_stats : dict with target normalization stats
    """
    rng = np.random.RandomState(seed)
    N = params.shape[0]

    # Input: normalized mags + photo_z
    mags_norm = normalize_magnitudes(mags)
    photo_z = params[:, 0]  # redshift is first param

    # Add small noise to photo_z to simulate photo-z uncertainty
    photo_z_noisy = photo_z + rng.normal(0, 0.02 * (1 + photo_z)).astype(np.float32)
    photo_z_noisy = np.maximum(photo_z_noisy, 0.0)

    X = np.column_stack([mags_norm, photo_z_noisy]).astype(np.float32)

    # Target: 5 physical params (excluding z, which is input)
    Y = params[:, 1:6].astype(np.float32)  # log_mass, log_age, log_Z, Av, log_tau

    # Train/val split
    n_val = int(N * val_fraction)
    idx = rng.permutation(N)
    train_idx = idx[n_val:]
    val_idx = idx[:n_val]

    # Compute normalization on training set
    X_mean = X[train_idx].mean(axis=0)
    X_std = X[train_idx].std(axis=0)
    X_std = np.maximum(X_std, 1e-8)

    Y_mean = Y[train_idx].mean(axis=0)
    Y_std = Y[train_idx].std(axis=0)
    Y_std = np.maximum(Y_std, 1e-8)

    # Normalize inputs
    X_norm = (X - X_mean) / X_std

    # Normalize targets (for better MDN training)
    Y_norm = (Y - Y_mean) / Y_std

    norm_stats = {'mean': X_mean, 'std': X_std}
    target_stats = {'mean': Y_mean, 'std': Y_std}

    return (X_norm[train_idx], Y_norm[train_idx],
            X_norm[val_idx], Y_norm[val_idx],
            norm_stats, target_stats)


def train(data_path=None, cfg=None, verbose=True):
    """
    Train the inverse SED model (MDN).

    Parameters
    ----------
    data_path : str
        Path to HDF5 SPS grid.
    cfg : TrainConfig, optional
    verbose : bool

    Returns
    -------
    model : InverseSEDMLP
    history : dict
    """
    if cfg is None:
        cfg = TrainConfig()

    sed_cfg = SEDConfig()

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

    # Prepare training data
    (train_X, train_Y, val_X, val_Y,
     norm_stats, target_stats) = prepare_training_data(
        params, mags, val_fraction=cfg.val_fraction, seed=cfg.seed)

    if verbose:
        print(f"Train: {len(train_X)}, Val: {len(val_X)}")
        print(f"Input dim: {train_X.shape[1]} (n_bands={n_bands} + photo_z)")
        print(f"Target dim: {train_Y.shape[1]} ({', '.join(OUTPUT_PARAM_NAMES)})")

    # Datasets
    train_ds = TensorDataset(
        torch.from_numpy(train_X),
        torch.from_numpy(train_Y))
    val_ds = TensorDataset(
        torch.from_numpy(val_X),
        torch.from_numpy(val_Y))

    use_cuda = device.type == 'cuda'
    n_workers = cfg.num_workers if use_cuda else 0
    train_loader = DataLoader(train_ds, batch_size=cfg.batch_size,
                              shuffle=True, num_workers=n_workers,
                              pin_memory=use_cuda)
    val_loader = DataLoader(val_ds, batch_size=cfg.batch_size * 2,
                            shuffle=False, num_workers=n_workers,
                            pin_memory=use_cuda)

    # Model
    n_input = n_bands + 1
    n_output = sed_cfg.n_output_params
    n_hidden = sed_cfg.n_hidden_inverse
    n_components = sed_cfg.n_components
    dropout = sed_cfg.dropout

    model = InverseSEDMLP(
        n_input=n_input,
        n_output=n_output,
        n_hidden=n_hidden,
        n_components=n_components,
        dropout=dropout,
    ).to(device)

    if verbose:
        print(f"Model parameters: {model.count_parameters():,}")
        print(f"MDN: {n_components} components, {n_output} output dims")

    # Optimizer + scheduler
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
            pi, mu, sigma = model(X_batch)
            loss = InverseSEDMLP.mdn_loss(pi, mu, sigma, Y_batch)
            loss.backward()

            # Gradient clipping for MDN stability
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)

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
                pi, mu, sigma = model(X_batch)
                loss = InverseSEDMLP.mdn_loss(pi, mu, sigma, Y_batch)
                val_loss += loss.item()
                n_val_batches += 1

        val_loss /= max(n_val_batches, 1)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)

        if verbose and (epoch % 10 == 0 or epoch == 1):
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch:3d}/{cfg.epochs}: "
                  f"train_nll={train_loss:.4f} val_nll={val_loss:.4f} "
                  f"lr={lr:.1e}")

        # Early stopping / checkpointing
        if val_loss < best_val_loss:
            best_val_loss = val_loss
            patience_counter = 0
            checkpoint = {
                'model_state': model.state_dict(),
                'epoch': epoch,
                'val_loss': float(best_val_loss),
                'norm_mean': norm_stats['mean'].tolist(),
                'norm_std': norm_stats['std'].tolist(),
                'target_mean': target_stats['mean'].tolist(),
                'target_std': target_stats['std'].tolist(),
                'config': {
                    'n_input': n_input,
                    'n_output': n_output,
                    'n_hidden': n_hidden,
                    'n_components': n_components,
                    'dropout': dropout,
                },
                'band_names': meta.get('band_names', []),
                'output_param_names': OUTPUT_PARAM_NAMES,
            }
            torch.save(checkpoint, cfg.checkpoint_inverse_path)
        else:
            patience_counter += 1
            if patience_counter >= cfg.early_stop_patience:
                if verbose:
                    print(f"Early stopping at epoch {epoch}")
                break

    elapsed = time.time() - t0
    if verbose:
        print(f"\nTraining complete in {elapsed:.1f}s")
        print(f"Best val NLL: {best_val_loss:.4f}")
        print(f"Checkpoint: {cfg.checkpoint_inverse_path}")

    # Reload best model and evaluate
    ckpt = torch.load(cfg.checkpoint_inverse_path, map_location=device,
                      weights_only=True)
    model.load_state_dict(ckpt['model_state'])
    model.eval()

    # Evaluate on validation set
    val_tensor = torch.from_numpy(val_X).to(device)
    pred = model.predict(val_tensor)

    # Un-normalize targets for comparison
    Y_mean = target_stats['mean']
    Y_std = target_stats['std']
    val_Y_orig = val_Y * Y_std + Y_mean

    if verbose:
        print(f"\n=== Validation Residuals ===")
        for i, name in enumerate(OUTPUT_PARAM_NAMES):
            pred_val = pred[name]['value'] * Y_std[i] + Y_mean[i]
            true_val = val_Y_orig[:, i]
            residuals = pred_val - true_val
            rms = np.sqrt(np.mean(residuals ** 2))
            bias = np.mean(residuals)
            print(f"  {name:10s}: bias={bias:+.4f}  RMS={rms:.4f}  "
                  f"median_err={np.median(pred[name]['error'] * Y_std[i]):.4f}")

    return model, history


def main():
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description='Train Inverse SED MLP+MDN model')
    parser.add_argument('--data', type=str, default=None,
                        help='Path to HDF5 SPS grid')
    parser.add_argument('--epochs', type=int, default=None)
    parser.add_argument('--batch-size', type=int, default=None)
    parser.add_argument('--lr', type=float, default=None)
    parser.add_argument('--patience', type=int, default=None)
    parser.add_argument('--n-components', type=int, default=None,
                        help='Number of MDN mixture components')
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
