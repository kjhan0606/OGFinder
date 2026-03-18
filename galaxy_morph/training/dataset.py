"""Galaxy10 DECaLS dataset for training."""

import os
import numpy as np
import torch
from torch.utils.data import Dataset

from ..config import GALAXY10_CLASSES


class Galaxy10Dataset(Dataset):
    """Galaxy10 DECaLS dataset: h5 -> RGB 3-channel 64x64.

    The Galaxy10_DECals.h5 file contains:
        images: (N, 256, 256, 3) uint8 RGB
        ans: (N,) int64 labels (0-9)

    Training uses RGB 3-channel to preserve color information
    which is critical for distinguishing spiral vs smooth galaxies.
    Inference on FITS (single-band) duplicates the channel 3x.
    """

    def __init__(self, images, labels, augment=False, size=64,
                 norm_mean=None, norm_std=None):
        """
        Args:
            images: (N, C, H, W) float32, C=3 for RGB
            labels: (N,) int labels
            augment: apply random augmentation
            size: expected spatial size
            norm_mean: (C,) per-channel mean for standardization
            norm_std: (C,) per-channel std for standardization
        """
        self.images = images
        self.labels = labels
        self.augment = augment
        self.size = size
        self.norm_mean = norm_mean
        self.norm_std = norm_std

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        img = self.images[idx].copy()  # (C, H, W) float32

        # Augmentation
        if self.augment:
            # Random 90-degree rotations (rotate last 2 dims)
            k = np.random.randint(4)
            img = np.rot90(img, k=k, axes=(1, 2)).copy()
            # Random flips
            if np.random.random() > 0.5:
                img = img[:, :, ::-1].copy()  # horizontal
            if np.random.random() > 0.5:
                img = img[:, ::-1, :].copy()  # vertical
            # Gaussian noise (small)
            if np.random.random() > 0.5:
                noise = np.random.normal(0, 0.02, img.shape).astype(np.float32)
                img = img + noise

        # Standardize per-channel
        if self.norm_mean is not None and self.norm_std is not None:
            for c in range(img.shape[0]):
                img[c] = (img[c] - self.norm_mean[c]) / (self.norm_std[c] + 1e-8)

        return torch.FloatTensor(img), torch.LongTensor([self.labels[idx]])[0]


def _preprocess_rgb(images_rgb, size=64):
    """Preprocess RGB uint8 images to float32 (N, 3, size, size).

    Args:
        images_rgb: (N, 256, 256, 3) uint8
        size: target resize dimension

    Returns:
        images: (N, 3, size, size) float32 in [0, 1]
    """
    from ..cutout.extract import _resize_2d

    n = len(images_rgb)
    out = np.zeros((n, 3, size, size), dtype=np.float32)

    for i in range(n):
        rgb = images_rgb[i].astype(np.float32) / 255.0  # [0, 1]

        for c in range(3):
            ch = rgb[:, :, c]
            if ch.shape[0] != size or ch.shape[1] != size:
                ch = _resize_2d(ch, size, size)
            out[i, c] = ch

    return out


def make_train_val_test_split(h5_path, size=64, val_frac=0.15, test_frac=0.05,
                               seed=42, asinh_a=0.1):
    """Load Galaxy10 h5 and split into train/val/test with stratification.

    Args:
        h5_path: path to Galaxy10_DECals.h5
        size: cutout resize target
        val_frac: validation fraction
        test_frac: test fraction
        seed: random seed
        asinh_a: unused, kept for API compat

    Returns:
        train_dataset, val_dataset, test_dataset, norm_stats
        norm_stats: dict with 'mean' and 'std' (per-channel, shape (3,))
    """
    import h5py

    print(f"Loading Galaxy10 from {h5_path} ...")
    with h5py.File(h5_path, 'r') as f:
        images_rgb = f['images'][:]    # (N, 256, 256, 3) uint8
        labels = f['ans'][:].astype(int)  # (N,)

    n = len(labels)
    print(f"  Total samples: {n}")
    for i, name in enumerate(GALAXY10_CLASSES):
        count = np.sum(labels == i)
        print(f"  Class {i} ({name}): {count}")

    # Preprocess RGB -> (N, 3, size, size) float32 [0,1]
    print(f"Preprocessing RGB {size}x{size} ...")
    images = _preprocess_rgb(images_rgb, size=size)
    del images_rgb  # free memory

    # Stratified split
    rng = np.random.RandomState(seed)
    train_idx, val_idx, test_idx = [], [], []

    for c in range(10):
        c_idx = np.where(labels == c)[0]
        rng.shuffle(c_idx)
        n_c = len(c_idx)
        n_test = max(1, int(n_c * test_frac))
        n_val = max(1, int(n_c * val_frac))

        test_idx.extend(c_idx[:n_test])
        val_idx.extend(c_idx[n_test:n_test + n_val])
        train_idx.extend(c_idx[n_test + n_val:])

    train_idx = np.array(train_idx)
    val_idx = np.array(val_idx)
    test_idx = np.array(test_idx)

    print(f"  Train: {len(train_idx)}, Val: {len(val_idx)}, Test: {len(test_idx)}")

    # Per-channel normalization stats from training set
    train_images = images[train_idx]
    norm_mean = np.array([train_images[:, c].mean() for c in range(3)], dtype=np.float32)
    norm_std = np.array([train_images[:, c].std() for c in range(3)], dtype=np.float32)
    norm_stats = {'mean': norm_mean, 'std': norm_std}
    print(f"  Norm stats: mean={norm_mean}, std={norm_std}")

    train_ds = Galaxy10Dataset(images[train_idx], labels[train_idx],
                                augment=True, size=size,
                                norm_mean=norm_mean, norm_std=norm_std)
    val_ds = Galaxy10Dataset(images[val_idx], labels[val_idx],
                              augment=False, size=size,
                              norm_mean=norm_mean, norm_std=norm_std)
    test_ds = Galaxy10Dataset(images[test_idx], labels[test_idx],
                               augment=False, size=size,
                               norm_mean=norm_mean, norm_std=norm_std)

    return train_ds, val_ds, test_ds, norm_stats
