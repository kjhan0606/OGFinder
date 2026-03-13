"""Training script for MLP merge classifier."""

import os
import sys
import time
import argparse
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader
from sklearn.metrics import f1_score, classification_report, confusion_matrix

from ..config import TrainConfig
from ..models.mlp import MergeMLP
from .dataset import MergePairDataset, make_train_val_split


def train(data_path=None, cfg=None, verbose=True):
    """
    Train the MLP merge classifier.

    Parameters
    ----------
    data_path : str
        Path to HDF5 training data.
    cfg : TrainConfig, optional
    verbose : bool

    Returns
    -------
    model : MergeMLP
    stats : dict
        Training statistics.
    """
    if cfg is None:
        cfg = TrainConfig()
    if data_path is None:
        data_path = "SAOImageDS9/ai_merge/data/simulated/pairs_v1.h5"

    # Set seed
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)

    # Device
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    if verbose:
        print(f"Device: {device}")

    # Data split
    if verbose:
        print(f"Loading data from {data_path}...")
    train_idx, val_idx, norm_stats = make_train_val_split(
        data_path, val_fraction=cfg.val_fraction, seed=cfg.seed
    )

    train_ds = MergePairDataset(data_path, indices=train_idx, normalize_stats=norm_stats)
    val_ds = MergePairDataset(data_path, indices=val_idx, normalize_stats=norm_stats)

    use_cuda = device.type == 'cuda'
    n_data_workers = 4 if use_cuda else 0
    train_loader = DataLoader(train_ds, batch_size=cfg.batch_size, shuffle=True,
                              num_workers=n_data_workers, pin_memory=use_cuda,
                              persistent_workers=(n_data_workers > 0))
    val_loader = DataLoader(val_ds, batch_size=cfg.batch_size, shuffle=False,
                            num_workers=n_data_workers, pin_memory=use_cuda,
                            persistent_workers=(n_data_workers > 0))

    if verbose:
        print(f"Train: {len(train_ds)}, Val: {len(val_ds)}")

    # Class weights
    train_labels = train_ds.labels
    class_counts = np.bincount(train_labels, minlength=cfg.n_classes).astype(float)
    class_weights = len(train_labels) / (cfg.n_classes * np.maximum(class_counts, 1))
    class_weights = torch.FloatTensor(class_weights).to(device)
    if verbose:
        print(f"Class counts: {class_counts.astype(int)}, Weights: {class_weights.cpu().numpy()}")

    # Model
    model = MergeMLP(
        n_features=cfg.n_features,
        n_classes=cfg.n_classes,
        hidden_dims=cfg.hidden_dims,
        dropout_rates=cfg.dropout_rates,
    ).to(device)

    if verbose:
        print(f"Model parameters: {model.count_parameters():,}")

    # Try torch.compile for faster execution (PyTorch 2.0+)
    compiled_model = model
    try:
        if hasattr(torch, 'compile') and use_cuda:
            compiled_model = torch.compile(model)
            if verbose:
                print("Using torch.compile for CUDA acceleration")
    except Exception:
        compiled_model = model

    # Loss, optimizer, scheduler
    criterion = nn.CrossEntropyLoss(weight=class_weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=cfg.lr,
                                   weight_decay=cfg.weight_decay)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode='max', factor=0.5, patience=5, verbose=verbose
    )

    # Mixed precision for CUDA
    use_amp = use_cuda
    scaler = torch.amp.GradScaler('cuda') if use_amp else None

    # Training loop
    best_val_f1 = 0.0
    patience_counter = 0
    history = {'train_loss': [], 'val_loss': [], 'val_acc': [], 'val_f1': []}

    os.makedirs(cfg.checkpoint_dir, exist_ok=True)

    t0 = time.time()
    for epoch in range(cfg.epochs):
        # Train
        compiled_model.train()
        train_loss = 0.0
        n_batches = 0
        for feats, labels in train_loader:
            feats = feats.to(device, non_blocking=True)
            labels = labels.to(device, non_blocking=True)
            optimizer.zero_grad(set_to_none=True)

            if use_amp:
                with torch.amp.autocast('cuda'):
                    logits = compiled_model(feats)
                    loss = criterion(logits, labels)
                scaler.scale(loss).backward()
                scaler.step(optimizer)
                scaler.update()
            else:
                logits = compiled_model(feats)
                loss = criterion(logits, labels)
                loss.backward()
                optimizer.step()

            train_loss += loss.item()
            n_batches += 1

        train_loss /= max(n_batches, 1)

        # Validate
        compiled_model.eval()
        val_loss = 0.0
        all_preds = []
        all_labels = []
        n_val_batches = 0
        with torch.no_grad():
            for feats, labels in val_loader:
                feats = feats.to(device, non_blocking=True)
                labels = labels.to(device, non_blocking=True)
                if use_amp:
                    with torch.amp.autocast('cuda'):
                        logits = compiled_model(feats)
                        loss = criterion(logits, labels)
                else:
                    logits = compiled_model(feats)
                    loss = criterion(logits, labels)
                val_loss += loss.item()
                preds = logits.argmax(dim=1)
                all_preds.extend(preds.cpu().numpy())
                all_labels.extend(labels.cpu().numpy())
                n_val_batches += 1

        val_loss /= max(n_val_batches, 1)
        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        val_acc = (all_preds == all_labels).mean()
        val_f1 = f1_score(all_labels, all_preds, average='macro')

        scheduler.step(val_f1)

        history['train_loss'].append(train_loss)
        history['val_loss'].append(val_loss)
        history['val_acc'].append(val_acc)
        history['val_f1'].append(val_f1)

        if verbose and (epoch + 1) % 5 == 0:
            lr = optimizer.param_groups[0]['lr']
            print(f"Epoch {epoch+1:3d}/{cfg.epochs}: "
                  f"loss={train_loss:.4f} val_loss={val_loss:.4f} "
                  f"acc={val_acc:.4f} f1={val_f1:.4f} lr={lr:.1e}")

        # Early stopping / checkpointing
        if val_f1 > best_val_f1:
            best_val_f1 = val_f1
            patience_counter = 0
            checkpoint = {
                'model_state': model.state_dict(),
                'norm_mean': torch.from_numpy(norm_stats['mean'].copy()),
                'norm_std': torch.from_numpy(norm_stats['std'].copy()),
                'best_val_f1': float(best_val_f1),
                'best_val_acc': float(val_acc),
                'epoch': int(epoch + 1),
                'config': {
                    'n_features': cfg.n_features,
                    'n_classes': cfg.n_classes,
                    'hidden_dims': cfg.hidden_dims,
                    'dropout_rates': cfg.dropout_rates,
                },
            }
            torch.save(checkpoint, cfg.checkpoint_path)
        else:
            patience_counter += 1
            if patience_counter >= cfg.patience:
                if verbose:
                    print(f"Early stopping at epoch {epoch+1}")
                break

    elapsed = time.time() - t0
    if verbose:
        print(f"\nTraining complete in {elapsed:.1f}s")
        print(f"Best val F1: {best_val_f1:.4f}")

        # Final evaluation on val set with best model
        ckpt = torch.load(cfg.checkpoint_path, map_location=device, weights_only=True)
        model.load_state_dict(ckpt['model_state'])
        model.eval()

        all_preds = []
        all_labels = []
        with torch.no_grad():
            for feats, labels in val_loader:
                feats = feats.to(device)
                logits = model(feats)
                all_preds.extend(logits.argmax(1).cpu().numpy())
                all_labels.extend(labels.numpy())

        all_preds = np.array(all_preds)
        all_labels = np.array(all_labels)
        print("\n=== Best Model Evaluation ===")
        print(classification_report(all_labels, all_preds,
                                     target_names=['Merge', 'Separate']))
        print("Confusion Matrix:")
        print(confusion_matrix(all_labels, all_preds))

    return model, history


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Train MLP merge classifier')
    parser.add_argument('--data', type=str,
                        default='SAOImageDS9/ai_merge/data/simulated/pairs_v1.h5',
                        help='Path to HDF5 training data')
    parser.add_argument('--epochs', type=int, default=None)
    parser.add_argument('--batch-size', type=int, default=None)
    parser.add_argument('--lr', type=float, default=None)
    args = parser.parse_args()

    cfg = TrainConfig()
    if args.epochs:
        cfg.epochs = args.epochs
    if args.batch_size:
        cfg.batch_size = args.batch_size
    if args.lr:
        cfg.lr = args.lr

    train(data_path=args.data, cfg=cfg)
