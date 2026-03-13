"""PyTorch Dataset for merge/separate pairs from HDF5."""

import numpy as np
import h5py
import torch
from torch.utils.data import Dataset


class MergePairDataset(Dataset):
    """
    Dataset loading features and labels from HDF5.

    The HDF5 file should contain:
    - 'features': (N, 20) float32
    - 'labels': (N,) int64
    """

    def __init__(self, h5_path, indices=None, normalize_stats=None):
        """
        Parameters
        ----------
        h5_path : str
        indices : np.ndarray, optional
            Subset of indices to use (for train/val split).
        normalize_stats : dict, optional
            {'mean': np.ndarray, 'std': np.ndarray} for normalization.
            If None, raw features are returned.
        """
        with h5py.File(h5_path, 'r') as f:
            self.features = f['features'][:].astype(np.float32)
            self.labels = f['labels'][:].astype(np.int64)

        if indices is not None:
            self.features = self.features[indices]
            self.labels = self.labels[indices]

        self.normalize_stats = normalize_stats
        if normalize_stats is not None:
            mean = normalize_stats['mean']
            std = normalize_stats['std']
            std = np.where(std < 1e-8, 1.0, std)
            self.features = (self.features - mean) / std

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        return (torch.from_numpy(self.features[idx]),
                torch.tensor(self.labels[idx], dtype=torch.long))


def make_train_val_split(h5_path, val_fraction=0.2, seed=42):
    """
    Create stratified train/val split.

    Returns
    -------
    train_indices, val_indices : np.ndarray
    normalize_stats : dict with 'mean' and 'std' computed from train set.
    """
    with h5py.File(h5_path, 'r') as f:
        labels = f['labels'][:]
        features = f['features'][:].astype(np.float32)

    rng = np.random.RandomState(seed)
    n = len(labels)

    # Stratified split
    train_idx = []
    val_idx = []
    for c in np.unique(labels):
        c_idx = np.where(labels == c)[0]
        rng.shuffle(c_idx)
        n_val = max(1, int(len(c_idx) * val_fraction))
        val_idx.extend(c_idx[:n_val])
        train_idx.extend(c_idx[n_val:])

    train_idx = np.array(train_idx)
    val_idx = np.array(val_idx)

    # Compute normalization stats from training set only
    train_feats = features[train_idx]
    stats = {
        'mean': train_feats.mean(axis=0),
        'std': train_feats.std(axis=0),
    }

    return train_idx, val_idx, stats
