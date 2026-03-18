"""Synthetic star/galaxy dataset for training.

Stars: Moffat PSF (point sources) with various FWHM/beta/ellipticity
Galaxies: Sersic profiles (extended sources) with various n/r_e/ellipticity

The CNN learns the fundamental distinction: point source vs extended source.
No external data dependency (Galaxy10 h5 not required).
"""

import numpy as np
import torch
from torch.utils.data import Dataset


class StarGalaxyDataset(Dataset):
    """Binary classification: star(1) vs galaxy(0)."""

    def __init__(self, images, labels, augment=False):
        self.images = images
        self.labels = labels
        self.augment = augment

    def __len__(self):
        return len(self.labels)

    def __getitem__(self, idx):
        img = self.images[idx].copy()  # (1, H, W)

        if self.augment:
            k = np.random.randint(4)
            img = np.rot90(img, k=k, axes=(1, 2)).copy()
            if np.random.random() > 0.5:
                img = img[:, :, ::-1].copy()
            if np.random.random() > 0.5:
                img = img[:, ::-1, :].copy()
            # Gaussian noise
            if np.random.random() > 0.5:
                noise = np.random.normal(0, 0.02, img.shape).astype(np.float32)
                img = img + noise
            # Random brightness shift
            if np.random.random() > 0.5:
                shift = np.random.uniform(-0.05, 0.05)
                img = img + shift
            img = np.clip(img, 0, 1)

        return torch.FloatTensor(img), torch.FloatTensor([self.labels[idx]])


def _asinh_normalize(stamp, asinh_a=0.1):
    """Asinh stretch + [0,1] normalization (same as inference pipeline)."""
    median = np.median(stamp)
    mad = np.median(np.abs(stamp - median))
    sigma = mad * 1.4826 if mad > 0 else 1.0
    stamp = (stamp - median) / sigma
    stamp = np.arcsinh(stamp / asinh_a)
    vmin, vmax = stamp.min(), stamp.max()
    if vmax > vmin:
        stamp = (stamp - vmin) / (vmax - vmin)
    else:
        stamp = np.zeros_like(stamp)
    return stamp.astype(np.float32)


def _make_elliptical_r(size, cx, cy, ellip, theta):
    """Compute elliptical radius grid.

    Args:
        size: grid size
        cx, cy: center
        ellip: ellipticity (0=circular, 0.9=very elongated)
        theta: position angle in radians

    Returns:
        r: 2D array of elliptical radii
    """
    y_grid, x_grid = np.mgrid[:size, :size].astype(np.float64)
    dx = x_grid - cx
    dy = y_grid - cy
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    x_rot = dx * cos_t + dy * sin_t
    y_rot = -dx * sin_t + dy * cos_t
    q = 1.0 - ellip  # axis ratio
    q = max(q, 0.1)
    r = np.sqrt(x_rot ** 2 + (y_rot / q) ** 2)
    return r


def generate_synthetic_stars(n=8000, size=64, fwhm_range=(1.5, 5.0), seed=42):
    """Generate synthetic star images (Moffat PSF + noise).

    Includes: variable FWHM, beta, slight ellipticity (tracking errors),
    variable flux/SNR, Poisson-like + Gaussian noise.
    """
    rng = np.random.RandomState(seed)
    stars = np.zeros((n, 1, size, size), dtype=np.float32)
    half = size / 2.0

    for i in range(n):
        fwhm = rng.uniform(fwhm_range[0], fwhm_range[1])
        beta = rng.uniform(2.5, 4.5)
        alpha = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))

        # Small random offset from center
        cx = half + rng.uniform(-2.0, 2.0)
        cy = half + rng.uniform(-2.0, 2.0)

        # Slight ellipticity for stars (tracking errors, optics)
        ellip = rng.uniform(0.0, 0.15)
        theta = rng.uniform(0, np.pi)

        r = _make_elliptical_r(size, cx, cy, ellip, theta)
        psf = (1.0 + (r / alpha) ** 2) ** (-beta)

        # Random flux (log-uniform, wide range)
        flux = 10.0 ** rng.uniform(2.0, 5.0)
        stamp = psf * flux / psf.sum()

        # Poisson noise (shot noise from source)
        stamp_pos = np.maximum(stamp, 0)
        if stamp_pos.max() > 0:
            poisson_noise = rng.poisson(stamp_pos) - stamp_pos
            stamp = stamp + poisson_noise * 0.3  # attenuated

        # Background + readout noise
        bkg_level = rng.uniform(50, 300)
        bkg_noise = rng.uniform(5, 40)
        stamp = stamp + bkg_level
        stamp = stamp + rng.normal(0, bkg_noise, stamp.shape)

        stars[i, 0] = _asinh_normalize(stamp)

    return stars


def generate_synthetic_galaxies(n=8000, size=64, seed=42):
    """Generate synthetic galaxy images (Sersic profiles + noise).

    Mix of galaxy types:
    - Exponential disk (Sersic n=1): 40%
    - De Vaucouleurs (Sersic n=4): 25%
    - General Sersic (n=0.5~6): 35%

    With variable half-light radius, ellipticity, position angle,
    flux, and noise (same noise model as stars for consistency).
    """
    rng = np.random.RandomState(seed)
    galaxies = np.zeros((n, 1, size, size), dtype=np.float32)
    half = size / 2.0

    for i in range(n):
        # Random Sersic index
        roll = rng.random()
        if roll < 0.40:
            sersic_n = rng.uniform(0.8, 1.2)   # exponential disk
        elif roll < 0.65:
            sersic_n = rng.uniform(3.5, 4.5)   # de Vaucouleurs
        else:
            sersic_n = rng.uniform(0.5, 6.0)   # general

        # b_n approximation: b_n ~ 2n - 1/3 + 4/(405n)
        bn = 2.0 * sersic_n - 1.0 / 3.0 + 4.0 / (405.0 * sersic_n)

        # Half-light radius: MUST be larger than typical PSF FWHM
        # to be clearly extended
        r_e = rng.uniform(3.0, 18.0)

        # Ellipticity and orientation
        ellip = rng.uniform(0.0, 0.7)
        theta = rng.uniform(0, np.pi)

        # Center offset
        cx = half + rng.uniform(-2.0, 2.0)
        cy = half + rng.uniform(-2.0, 2.0)

        r = _make_elliptical_r(size, cx, cy, ellip, theta)
        r = np.maximum(r, 0.01)  # avoid division by zero

        # Sersic profile: I(r) = I_e * exp(-b_n * ((r/r_e)^(1/n) - 1))
        profile = np.exp(-bn * ((r / r_e) ** (1.0 / sersic_n) - 1.0))

        # Random flux
        flux = 10.0 ** rng.uniform(2.0, 5.0)
        stamp = profile * flux / profile.sum()

        # Poisson noise
        stamp_pos = np.maximum(stamp, 0)
        if stamp_pos.max() > 0:
            poisson_noise = rng.poisson(np.minimum(stamp_pos, 1e6)) - stamp_pos
            stamp = stamp + poisson_noise * 0.3

        # Background + readout noise (same model as stars)
        bkg_level = rng.uniform(50, 300)
        bkg_noise = rng.uniform(5, 40)
        stamp = stamp + bkg_level
        stamp = stamp + rng.normal(0, bkg_noise, stamp.shape)

        galaxies[i, 0] = _asinh_normalize(stamp)

    return galaxies


def make_train_val_test_split(star_images, galaxy_images,
                               val_frac=0.15, test_frac=0.15, seed=42):
    """Combine and split star/galaxy data."""
    n_star = len(star_images)
    n_gal = len(galaxy_images)

    images = np.concatenate([star_images, galaxy_images], axis=0)
    labels = np.concatenate([
        np.ones(n_star, dtype=np.int64),   # star=1
        np.zeros(n_gal, dtype=np.int64),   # galaxy=0
    ])

    rng = np.random.RandomState(seed)
    perm = rng.permutation(len(labels))
    images = images[perm]
    labels = labels[perm]

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
