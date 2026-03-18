"""Training loop for Galaxy Morphology CNN."""

import os
import sys
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from ..config import TrainConfig, GALAXY10_CLASSES
from ..models.cnn import GalaxyMorphCNN
from .dataset import make_train_val_test_split


def compute_class_weights(dataset):
    """Compute inverse-frequency class weights for imbalanced data."""
    labels = dataset.labels
    counts = np.bincount(labels, minlength=10).astype(np.float64)
    counts = np.maximum(counts, 1.0)  # avoid div by zero
    weights = 1.0 / counts
    weights = weights / weights.sum() * len(counts)  # normalize
    return torch.FloatTensor(weights)


def train(h5_path=None, config=None):
    """Main training loop with early stopping.

    Args:
        h5_path: path to Galaxy10_DECals.h5
        config: TrainConfig instance
    """
    if config is None:
        config = TrainConfig()

    if h5_path is None:
        h5_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data", "galaxy10", "Galaxy10_DECals.h5"
        )

    if not os.path.exists(h5_path):
        print(f"ERROR: Galaxy10 dataset not found at {h5_path}")
        print("Download it with:")
        print(f"  wget -P {os.path.dirname(h5_path)} "
              "https://astro.utoronto.ca/~hleung/shared/Galaxy10/Galaxy10_DECals.h5")
        sys.exit(1)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Load data
    train_ds, val_ds, test_ds, norm_stats = make_train_val_test_split(
        h5_path, size=config.cutout_size,
        val_frac=config.val_fraction, test_frac=config.test_fraction,
        seed=config.seed
    )

    train_loader = DataLoader(train_ds, batch_size=config.batch_size,
                               shuffle=True, num_workers=config.num_workers,
                               pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size * 2,
                             shuffle=False, num_workers=config.num_workers,
                             pin_memory=True)

    # Model (3 channels for RGB training)
    model = GalaxyMorphCNN(
        n_classes=config.n_classes,
        in_channels=3,
        dropout=config.dropout,
        dropout2=config.dropout2,
        hidden_dim=config.hidden_dim,
    ).to(device)
    print(f"Model parameters: {model.count_parameters():,}")

    # Class-weighted loss with label smoothing
    class_weights = compute_class_weights(train_ds).to(device)
    print(f"Class weights: {class_weights.cpu().numpy()}")
    criterion = nn.CrossEntropyLoss(weight=class_weights, label_smoothing=0.1)

    # Optimizer + scheduler
    optimizer = torch.optim.AdamW(model.parameters(),
                                  lr=config.lr,
                                  weight_decay=config.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=config.epochs)

    # Checkpoint dir
    os.makedirs(config.checkpoint_dir, exist_ok=True)

    best_val_acc = 0.0
    patience_counter = 0

    for epoch in range(1, config.epochs + 1):
        t0 = time.time()

        # Train
        model.train()
        train_loss = 0.0
        train_correct = 0
        train_total = 0

        for images, labels in train_loader:
            images = images.to(device)
            labels = labels.to(device)

            optimizer.zero_grad()
            outputs = model(images)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item() * images.size(0)
            _, predicted = outputs.max(1)
            train_correct += predicted.eq(labels).sum().item()
            train_total += images.size(0)

        scheduler.step()

        train_loss /= train_total
        train_acc = train_correct / train_total

        # Validate
        model.eval()
        val_loss = 0.0
        val_correct = 0
        val_total = 0

        with torch.no_grad():
            for images, labels in val_loader:
                images = images.to(device)
                labels = labels.to(device)
                outputs = model(images)
                loss = criterion(outputs, labels)
                val_loss += loss.item() * images.size(0)
                _, predicted = outputs.max(1)
                val_correct += predicted.eq(labels).sum().item()
                val_total += images.size(0)

        val_loss /= val_total
        val_acc = val_correct / val_total

        elapsed = time.time() - t0
        lr = optimizer.param_groups[0]['lr']
        print(f"Epoch {epoch:3d}/{config.epochs} | "
              f"train_loss={train_loss:.4f} train_acc={train_acc:.4f} | "
              f"val_loss={val_loss:.4f} val_acc={val_acc:.4f} | "
              f"lr={lr:.6f} | {elapsed:.1f}s")

        # Early stopping
        if val_acc > best_val_acc:
            best_val_acc = val_acc
            patience_counter = 0
            # Save checkpoint
            ckpt = {
                'model_state': model.state_dict(),
                'epoch': epoch,
                'val_acc': val_acc,
                'norm_mean': norm_stats['mean'].tolist() if hasattr(norm_stats['mean'], 'tolist') else norm_stats['mean'],
                'norm_std': norm_stats['std'].tolist() if hasattr(norm_stats['std'], 'tolist') else norm_stats['std'],
                'class_names': GALAXY10_CLASSES,
                'config': {
                    'n_classes': config.n_classes,
                    'in_channels': 3,
                    'dropout': config.dropout,
                    'dropout2': config.dropout2,
                    'hidden_dim': config.hidden_dim,
                    'cutout_size': config.cutout_size,
                },
            }
            ckpt_path = config.checkpoint_path
            torch.save(ckpt, ckpt_path)
            print(f"  -> Saved best model (val_acc={val_acc:.4f})")
        else:
            patience_counter += 1
            if patience_counter >= config.patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(best val_acc={best_val_acc:.4f})")
                break

    print(f"\nBest validation accuracy: {best_val_acc:.4f}")
    print(f"Checkpoint: {config.checkpoint_path}")

    # Evaluate on test set
    print("\n--- Test Set Evaluation ---")
    test_loader = DataLoader(test_ds, batch_size=config.batch_size * 2,
                              shuffle=False, num_workers=config.num_workers)

    # Load best model
    ckpt = torch.load(config.checkpoint_path, map_location=device,
                       weights_only=False)
    model.load_state_dict(ckpt['model_state'])
    model.eval()

    test_correct = 0
    test_total = 0
    per_class_correct = np.zeros(config.n_classes)
    per_class_total = np.zeros(config.n_classes)

    with torch.no_grad():
        for images, labels in test_loader:
            images = images.to(device)
            labels = labels.to(device)
            outputs = model(images)
            _, predicted = outputs.max(1)
            test_correct += predicted.eq(labels).sum().item()
            test_total += images.size(0)

            for c in range(config.n_classes):
                mask = labels == c
                per_class_correct[c] += predicted[mask].eq(labels[mask]).sum().item()
                per_class_total[c] += mask.sum().item()

    test_acc = test_correct / test_total
    print(f"Test accuracy: {test_acc:.4f} ({test_correct}/{test_total})")
    print("\nPer-class accuracy:")
    for c in range(config.n_classes):
        if per_class_total[c] > 0:
            acc = per_class_correct[c] / per_class_total[c]
            print(f"  {c}: {GALAXY10_CLASSES[c]:25s} "
                  f"{acc:.4f} ({int(per_class_correct[c])}/{int(per_class_total[c])})")

    return best_val_acc


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Train Galaxy Morphology CNN')
    parser.add_argument('--data', default=None, help='Path to Galaxy10_DECals.h5')
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--patience', type=int, default=15)
    args = parser.parse_args()

    cfg = TrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        patience=args.patience,
    )
    train(h5_path=args.data, config=cfg)
