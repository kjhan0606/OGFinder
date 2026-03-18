"""Synthetic star/galaxy dataset for training."""

import os
import numpy as np
import torch
from torch.utils.data import Dataset

from ..cutout.extract import _resize_2d


class StarGalaxyDataset(Dataset):
    """Binary classification: star(1) vs galaxy(0)."""

    def __init__(self, images, labels, augment=False):
        """
        Args:
            images: (N, 1, H, W) float32
            labels: (N,) int labels (0=galaxy, 1=star)
            augment: apply random augmentation
        """
        self.images = images
        self.labels = labels
        self.augment = augment

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        img = self.images[idx].copy()  # (1, H, W)

        if self.augment:
            # Random 90-degree rotations
            k = np.random.randint(4)
            img = np.rot90(img, k=k, axes=(1, 2)).copy()
            # Random flips
            if np.random.random() > 0.5:
                img = img[:, :, ::-1].copy()
            if np.random.random() > 0.5:
                img = img[:, ::-1, :].copy()
            # Small Gaussian noise
            if np.random.random() > 0.5:
                noise = np.random.normal(0, 0.02, img.shape).astype(np.float32)
                img = img + noise

        return torch.FloatTensor(img), torch.FloatTensor([self.labels[idx]])


def generate_synthetic_stars(n=8000, size=64, fwhm_range=(1.5, 5.0), seed=42):
    """Generate synthetic star images (Moffat PSF + noise).

    Args:
        n: number of stars to generate
        size: stamp size
        fwhm_range: (min, max) FWHM in pixels
        seed: random seed

    Returns:
        stars: (n, 1, size, size) float32 normalized [0, 1]
    """
    rng = np.random.RandomState(seed)
    stars = np.zeros((n, 1, size, size), dtype=np.float32)
    half = size / 2.0

    y_grid, x_grid = np.mgrid[:size, :size].astype(np.float64)

    for i in range(n):
        fwhm = rng.uniform(fwhm_range[0], fwhm_range[1])
        # Moffat profile: beta=2.5..4.5
        beta = rng.uniform(2.5, 4.5)
        alpha = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))

        # Small random offset from center
        cx = half + rng.uniform(-1.5, 1.5)
        cy = half + rng.uniform(-1.5, 1.5)

        r2 = (x_grid - cx) ** 2 + (y_grid - cy) ** 2
        psf = (1.0 + r2 / alpha ** 2) ** (-beta)

        # Random flux (log-uniform)
        flux = 10.0 ** rng.uniform(2.0, 5.0)
        stamp = psf * flux / psf.sum()

        # Add Poisson-like noise + Gaussian background
        bkg_level = rng.uniform(50, 200)
        bkg_noise = rng.uniform(5, 30)
        stamp = stamp + bkg_level
        stamp = stamp + rng.normal(0, bkg_noise, stamp.shape)

        # Normalize: asinh stretch
        median = np.median(stamp)
        mad = np.median(np.abs(stamp - median))
        sigma = mad * 1.4826 if mad > 0 else 1.0
        stamp = (stamp - median) / sigma
        stamp = np.arcsinh(stamp / 0.1)
        vmin, vmax = stamp.min(), stamp.max()
        if vmax > vmin:
            stamp = (stamp - vmin) / (vmax - vmin)
        else:
            stamp = np.zeros_like(stamp)

        stars[i, 0] = stamp.astype(np.float32)

    return stars


def load_galaxy_samples(h5_path, n=8000, size=64, seed=42):
    """Load galaxy images from Galaxy10 DECaLS H5 -> grayscale.

    Args:
        h5_path: path to Galaxy10_DECals.h5
        n: number of galaxies to sample
        size: output stamp size
        seed: random seed

    Returns:
        galaxies: (n, 1, size, size) float32 normalized [0, 1]
    """
    import h5py

    with h5py.File(h5_path, 'r') as f:
        images_rgb = f['images'][:]  # (N, 256, 256, 3) uint8

    total = len(images_rgb)
    rng = np.random.RandomState(seed)
    indices = rng.choice(total, size=min(n, total), replace=False)

    galaxies = np.zeros((len(indices), 1, size, size), dtype=np.float32)

    for i, idx in enumerate(indices):
        rgb = images_rgb[idx].astype(np.float32)  # (256, 256, 3)
        # Convert to grayscale: 0.299R + 0.587G + 0.114B
        gray = 0.299 * rgb[:, :, 0] + 0.587 * rgb[:, :, 1] + 0.114 * rgb[:, :, 2]
        gray = gray / 255.0  # [0, 1]

        # Resize to target size
        if gray.shape[0] != size or gray.shape[1] != size:
            gray = _resize_2d(gray, size, size)

        # Asinh stretch (already roughly [0,1], apply mild stretch)
        median = np.median(gray)
        mad = np.median(np.abs(gray - median))
        sigma = mad * 1.4826 if mad > 0 else 1.0
        if sigma > 0:
            gray = (gray - median) / sigma
            gray = np.arcsinh(gray / 0.1)
            vmin, vmax = gray.min(), gray.max()
            if vmax > vmin:
                gray = (gray - vmin) / (vmax - vmin)
            else:
                gray = np.zeros_like(gray)

        galaxies[i, 0] = gray.astype(np.float32)

    return galaxies


def make_train_val_test_split(star_images, galaxy_images,
                               val_frac=0.15, test_frac=0.15, seed=42):
    """Combine and split star/galaxy data.

    Args:
        star_images: (N_star, 1, H, W) float32
        galaxy_images: (N_gal, 1, H, W) float32
        val_frac: validation fraction
        test_frac: test fraction
        seed: random seed

    Returns:
        train_ds, val_ds, test_ds
    """
    n_star = len(star_images)
    n_gal = len(galaxy_images)

    # Combine
    images = np.concatenate([star_images, galaxy_images], axis=0)
    labels = np.concatenate([
        np.ones(n_star, dtype=np.int64),   # star=1
        np.zeros(n_gal, dtype=np.int64),   # galaxy=0
    ])

    # Shuffle
    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(labels))
    images = images[perm]
    labels = labels[perm]

    # Split
    n = len(labels)
    n_test = int(n * test_frac)
    n_val = int(n * val_frac)
    n_train = n - n_test - n_val

    train_ds = StarGalaxyDataset(images[:n_train], labels[:n_train], augment=True)
    val_ds = StarGalaxyDataset(images[n_train:n_train + n_val],
                                labels[n_train:n_train + n_val], augment=False)
    test_ds = StarGalaxyDataset(images[n_train + n_val:],
                                 labels[n_train + n_val:], augment=False)

    print(f"  Train: {n_train} (star={int(labels[:n_train].sum())}, "
          f"galaxy={n_train - int(labels[:n_train].sum())})")
    print(f"  Val:   {n_val} (star={int(labels[n_train:n_train+n_val].sum())}, "
          f"galaxy={n_val - int(labels[n_train:n_train+n_val].sum())})")
    print(f"  Test:  {n_test} (star={int(labels[n_train+n_val:].sum())}, "
          f"galaxy={n_test - int(labels[n_train+n_val:].sum())})")

    return train_ds, val_ds, test_ds
