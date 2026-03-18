"""Training loop for Star/Galaxy binary classifier."""

import os
import sys
import time
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader

from ..config import TrainConfig
from ..models.cnn import StarFinderCNN
from .dataset import generate_synthetic_stars, load_galaxy_samples, make_train_val_test_split


def train(h5_path=None, config=None):
    """Main training loop with early stopping.

    Args:
        h5_path: path to Galaxy10_DECals.h5 (for galaxy samples)
        config: TrainConfig instance
    """
    if config is None:
        config = TrainConfig()

    if h5_path is None:
        h5_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
            "galaxy_morph", "data", "galaxy10", "Galaxy10_DECals.h5"
        )

    if not os.path.exists(h5_path):
        print(f"ERROR: Galaxy10 dataset not found at {h5_path}")
        print("Download it with:")
        print(f"  wget -P {os.path.dirname(h5_path)} "
              "https://astro.utoronto.ca/~hleung/shared/Galaxy10/Galaxy10_DECals.h5")
        sys.exit(1)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Device: {device}")

    # Generate synthetic stars
    print(f"Generating {config.n_stars} synthetic stars "
          f"(FWHM {config.fwhm_range_min}-{config.fwhm_range_max} px) ...")
    star_images = generate_synthetic_stars(
        n=config.n_stars, size=config.cutout_size,
        fwhm_range=(config.fwhm_range_min, config.fwhm_range_max),
        seed=config.seed)

    # Load galaxy samples
    print(f"Loading {config.n_galaxies} galaxy samples from {h5_path} ...")
    galaxy_images = load_galaxy_samples(
        h5_path, n=config.n_galaxies, size=config.cutout_size,
        seed=config.seed + 1)

    print(f"Stars: {len(star_images)}, Galaxies: {len(galaxy_images)}")

    # Split
    print("Splitting data ...")
    train_ds, val_ds, test_ds = make_train_val_test_split(
        star_images, galaxy_images,
        val_frac=config.val_fraction, test_frac=config.test_fraction,
        seed=config.seed)

    train_loader = DataLoader(train_ds, batch_size=config.batch_size,
                               shuffle=True, num_workers=config.num_workers,
                               pin_memory=True)
    val_loader = DataLoader(val_ds, batch_size=config.batch_size * 2,
                             shuffle=False, num_workers=config.num_workers,
                             pin_memory=True)

    # Model
    model = StarFinderCNN(in_channels=1).to(device)
    print(f"Model parameters: {model.count_parameters():,}")

    # BCE loss
    criterion = nn.BCELoss()

    # Optimizer + scheduler
    optimizer = torch.optim.Adam(model.parameters(), lr=config.lr,
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
            outputs = model(images)  # (B, 1)
            loss = criterion(outputs, labels)
            loss.backward()
            optimizer.step()

            train_loss += loss.item() * images.size(0)
            predicted = (outputs >= 0.5).float()
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
                predicted = (outputs >= 0.5).float()
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
            ckpt = {
                'model_state': model.state_dict(),
                'epoch': epoch,
                'val_acc': val_acc,
                'config': {
                    'in_channels': 1,
                    'cutout_size': config.cutout_size,
                },
            }
            ckpt_path = config.checkpoint_path
            torch.save(ckpt, ckpt_path)
            print(f"  -> Saved best model (val_acc={val_acc:.4f})")
        else:
            patience_counter += 1
            if patience_counter >= config.early_stop_patience:
                print(f"Early stopping at epoch {epoch} "
                      f"(best val_acc={best_val_acc:.4f})")
                break

    print(f"\nBest validation accuracy: {best_val_acc:.4f}")
    print(f"Checkpoint: {config.checkpoint_path}")

    # Evaluate on test set
    print("\n--- Test Set Evaluation ---")
    test_loader = DataLoader(test_ds, batch_size=config.batch_size * 2,
                              shuffle=False, num_workers=config.num_workers)

    ckpt = torch.load(config.checkpoint_path, map_location=device,
                       weights_only=False)
    model.load_state_dict(ckpt['model_state'])
    model.eval()

    test_correct = 0
    test_total = 0
    tp = fp = tn = fn = 0

    with torch.no_grad():
        for images, labels in test_loader:
            images = images.to(device)
            labels = labels.to(device)
            outputs = model(images)
            predicted = (outputs >= 0.5).float()
            test_correct += predicted.eq(labels).sum().item()
            test_total += images.size(0)

            tp += ((predicted == 1) & (labels == 1)).sum().item()
            fp += ((predicted == 1) & (labels == 0)).sum().item()
            tn += ((predicted == 0) & (labels == 0)).sum().item()
            fn += ((predicted == 0) & (labels == 1)).sum().item()

    test_acc = test_correct / test_total
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"Test accuracy: {test_acc:.4f} ({test_correct}/{test_total})")
    print(f"  Star precision: {precision:.4f}, recall: {recall:.4f}, F1: {f1:.4f}")
    print(f"  TP={tp}, FP={fp}, TN={tn}, FN={fn}")

    return best_val_acc


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Train Star/Galaxy Classifier')
    parser.add_argument('--data', default=None, help='Path to Galaxy10_DECals.h5')
    parser.add_argument('--epochs', type=int, default=100)
    parser.add_argument('--batch-size', type=int, default=64)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--patience', type=int, default=15)
    parser.add_argument('--n-stars', type=int, default=8000)
    parser.add_argument('--n-galaxies', type=int, default=8000)
    args = parser.parse_args()

    cfg = TrainConfig(
        epochs=args.epochs,
        batch_size=args.batch_size,
        lr=args.lr,
        early_stop_patience=args.patience,
        n_stars=args.n_stars,
        n_galaxies=args.n_galaxies,
    )
    train(h5_path=args.data, config=cfg)
