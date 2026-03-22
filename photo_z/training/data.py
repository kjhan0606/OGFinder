"""SDSS spectro-photometric training data loader for photo-z."""

import numpy as np
import torch
from torch.utils.data import Dataset, Subset


class SDSSPhotoZDataset(Dataset):
    """
    Dataset for photometric redshift training from HDF5 file.

    Expected HDF5 structure:
        /features  - shape (N, n_features), float32
        /redshift  - shape (N,), float32 (spectroscopic redshifts)

    Optionally:
        /feature_names - string array of feature column names
        /feat_mean     - shape (n_features,) normalization mean
        /feat_std      - shape (n_features,) normalization std

    Parameters
    ----------
    features : numpy.ndarray
        Shape (N, n_features).
    redshifts : numpy.ndarray
        Shape (N,).
    """

    def __init__(self, features, redshifts):
        self.features = torch.from_numpy(
            np.asarray(features, dtype=np.float32))
        self.redshifts = torch.from_numpy(
            np.asarray(redshifts, dtype=np.float32))

    def __len__(self):
        return len(self.redshifts)

    def __getitem__(self, idx):
        return self.features[idx], self.redshifts[idx]


def load_sdss_specphoto(path):
    """
    Load SDSS spectro-photometric data from HDF5 file.

    Parameters
    ----------
    path : str
        Path to HDF5 file.

    Returns
    -------
    features : numpy.ndarray
        Shape (N, n_features).
    redshifts : numpy.ndarray
        Shape (N,).
    metadata : dict
        Additional metadata (feature_names, feat_mean, feat_std if available).
    """
    import h5py

    metadata = {}
    with h5py.File(path, 'r') as f:
        features = f['features'][:]
        redshifts = f['redshift'][:]

        if 'feature_names' in f:
            metadata['feature_names'] = [
                s.decode('utf-8') if isinstance(s, bytes) else s
                for s in f['feature_names'][:]
            ]
        if 'feat_mean' in f:
            metadata['feat_mean'] = f['feat_mean'][:]
        if 'feat_std' in f:
            metadata['feat_std'] = f['feat_std'][:]

    return features, redshifts, metadata


def create_datasets(data_path, cfg):
    """
    Create train/val/test datasets from HDF5 data file.

    Parameters
    ----------
    data_path : str
        Path to HDF5 file.
    cfg : TrainConfig
        Training configuration with val_fraction, test_fraction, seed.

    Returns
    -------
    train_ds : Subset
        Training dataset.
    val_ds : Subset
        Validation dataset.
    test_ds : Subset
        Test dataset.
    metadata : dict
        Data metadata including normalization parameters.
    """
    features, redshifts, metadata = load_sdss_specphoto(data_path)

    n = len(redshifts)
    rng = np.random.RandomState(cfg.seed)
    indices = rng.permutation(n)

    n_test = int(n * cfg.test_fraction)
    n_val = int(n * cfg.val_fraction)
    n_train = n - n_val - n_test

    train_idx = indices[:n_train]
    val_idx = indices[n_train:n_train + n_val]
    test_idx = indices[n_train + n_val:]

    # Compute normalization from training set only
    feat_mean = features[train_idx].mean(axis=0)
    feat_std = features[train_idx].std(axis=0)
    feat_std[feat_std < 1e-8] = 1.0

    # Normalize all features
    features_norm = (features - feat_mean) / feat_std

    metadata['feat_mean'] = feat_mean
    metadata['feat_std'] = feat_std
    metadata['n_features'] = features.shape[1]

    # Create full dataset with normalized features
    full_ds = SDSSPhotoZDataset(features_norm, redshifts)

    train_ds = Subset(full_ds, train_idx.tolist())
    val_ds = Subset(full_ds, val_idx.tolist())
    test_ds = Subset(full_ds, test_idx.tolist())

    return train_ds, val_ds, test_ds, metadata
