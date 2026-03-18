"""Synthetic star/galaxy dataset for training.

Stars: Moffat PSF + diffraction spikes + Airy rings on zero-background
Galaxies: Sersic profile CONVOLVED with clean PSF (no spikes)

Key distinction the CNN learns:
  - PSF structure (spikes, rings) → star
  - Extended smooth profile → galaxy
  - NOT just compactness (which CLASS_STAR already does)
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
            if np.random.random() > 0.5:
                noise = np.random.normal(0, 0.01, img.shape).astype(np.float32)
                img = np.clip(img + noise, 0, 1)

        return torch.FloatTensor(img), torch.FloatTensor([self.labels[idx]])


# ---------------------------------------------------------------------------
#  Normalization (identical to cutout/extract.normalize_cutout)
# ---------------------------------------------------------------------------

def _asinh_normalize(stamp, asinh_a=0.1):
    """Asinh stretch + [0,1]."""
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


# ---------------------------------------------------------------------------
#  PSF components
# ---------------------------------------------------------------------------

def _make_moffat(size, fwhm, beta, cx=None, cy=None):
    """Normalized Moffat PSF core."""
    if cx is None:
        cx = size / 2.0
    if cy is None:
        cy = size / 2.0
    alpha = fwhm / (2.0 * np.sqrt(2.0 ** (1.0 / beta) - 1.0))
    y, x = np.mgrid[:size, :size].astype(np.float64)
    r2 = (x - cx) ** 2 + (y - cy) ** 2
    psf = (1.0 + r2 / alpha ** 2) ** (-beta)
    s = psf.sum()
    if s > 0:
        psf /= s
    return psf


def _add_diffraction_spikes(stamp, cx, cy, n_spikes, spike_angle_offset,
                             spike_length, spike_width, spike_brightness,
                             peak_val):
    """Add diffraction spikes to a star stamp.

    Args:
        stamp: 2D array (modified in-place)
        cx, cy: star center
        n_spikes: number of spikes (4=HST, 6=JWST, 8=ground)
        spike_angle_offset: rotation angle (radians)
        spike_length: spike length in pixels
        spike_width: Gaussian sigma of spike cross-section
        spike_brightness: fraction of peak brightness at r=0 for spike
        peak_val: peak pixel value of the star
    """
    size = stamp.shape[0]
    y_grid, x_grid = np.mgrid[:size, :size].astype(np.float64)

    for k in range(n_spikes):
        angle = k * np.pi / (n_spikes / 2.0) + spike_angle_offset

        # Direction vector
        dx = np.cos(angle)
        dy = np.sin(angle)

        # Distance along spike direction from center
        along = (x_grid - cx) * dx + (y_grid - cy) * dy
        # Perpendicular distance
        perp = np.abs(-(x_grid - cx) * dy + (y_grid - cy) * dx)

        # Spike profile: Gaussian cross-section x exponential falloff along
        mask = along > 0  # only in positive direction (full spike = 2 per pair)
        falloff = np.exp(-along / (spike_length / 3.0))
        cross = np.exp(-0.5 * (perp / spike_width) ** 2)
        spike = spike_brightness * peak_val * falloff * cross * mask

        stamp += spike


def _add_airy_rings(stamp, cx, cy, fwhm, n_rings=3, ring_brightness=0.02,
                     peak_val=1.0):
    """Add faint Airy ring structure around bright stars."""
    size = stamp.shape[0]
    y, x = np.mgrid[:size, :size].astype(np.float64)
    r = np.sqrt((x - cx) ** 2 + (y - cy) ** 2)
    r = np.maximum(r, 0.01)

    # Approximate Airy pattern as sinc-like rings
    # First ring at ~1.22 * lambda/D ≈ 1.6 * FWHM
    ring_spacing = fwhm * 1.6
    for n in range(1, n_rings + 1):
        ring_r = n * ring_spacing
        ring_width = fwhm * 0.3
        ring = np.exp(-0.5 * ((r - ring_r) / ring_width) ** 2)
        brightness = ring_brightness * peak_val / (n ** 1.5)
        stamp += brightness * ring


# ---------------------------------------------------------------------------
#  Star generation
# ---------------------------------------------------------------------------

def generate_synthetic_stars(n=10000, size=64, fwhm_range=(1.5, 5.0), seed=42):
    """Synthetic stars: Moffat PSF + diffraction spikes + Airy rings.

    - Faint stars (SNR < 30): Moffat only (spikes not visible)
    - Moderate stars (30 < SNR < 100): Moffat + faint spikes
    - Bright stars (SNR > 100): Moffat + bright spikes + Airy rings

    Spike types: 4-spike (HST-like), 6-spike (JWST-like), variable
    """
    rng = np.random.RandomState(seed)
    stars = np.zeros((n, 1, size, size), dtype=np.float32)

    for i in range(n):
        fwhm = rng.uniform(fwhm_range[0], fwhm_range[1])
        beta = rng.uniform(2.5, 4.5)
        cx = size / 2.0 + rng.uniform(-1.5, 1.5)
        cy = size / 2.0 + rng.uniform(-1.5, 1.5)

        # Moffat core
        psf = _make_moffat(size, fwhm, beta, cx, cy)

        # SNR range (log-uniform)
        snr = 10.0 ** rng.uniform(0.7, 2.7)
        noise_sigma = rng.uniform(0.005, 0.05)
        peak_flux = snr * noise_sigma

        stamp = psf * peak_flux
        peak_val = stamp.max()

        # Add diffraction spikes for bright stars
        if snr > 25:
            # Choose spike type: 4 (HST), 6 (JWST), or 8 (ground)
            spike_type = rng.choice([4, 4, 4, 6, 6, 8])
            spike_angle = rng.uniform(0, np.pi / spike_type)

            # Spike brightness scales with SNR
            if snr > 100:
                spike_bright = rng.uniform(0.05, 0.20)
                spike_len = rng.uniform(15, 28)
            elif snr > 50:
                spike_bright = rng.uniform(0.02, 0.10)
                spike_len = rng.uniform(10, 20)
            else:
                spike_bright = rng.uniform(0.01, 0.05)
                spike_len = rng.uniform(5, 15)

            spike_width = rng.uniform(0.5, 1.5)

            _add_diffraction_spikes(stamp, cx, cy, spike_type,
                                     spike_angle, spike_len, spike_width,
                                     spike_bright, peak_val)

        # Add Airy rings for very bright stars
        if snr > 80:
            n_rings = rng.randint(1, 4)
            ring_bright = rng.uniform(0.005, 0.03)
            _add_airy_rings(stamp, cx, cy, fwhm, n_rings, ring_bright, peak_val)

        # Gaussian noise (background-subtracted)
        stamp += rng.normal(0, noise_sigma, stamp.shape)

        stars[i, 0] = _asinh_normalize(stamp)

    return stars


# ---------------------------------------------------------------------------
#  Sersic galaxy generation
# ---------------------------------------------------------------------------

def _make_sersic(size, r_e, sersic_n, ellip, theta, cx=None, cy=None):
    """Sersic profile, normalized to sum=1."""
    if cx is None:
        cx = size / 2.0
    if cy is None:
        cy = size / 2.0
    bn = 2.0 * sersic_n - 1.0 / 3.0 + 4.0 / (405.0 * sersic_n)
    y, x = np.mgrid[:size, :size].astype(np.float64)
    dx = x - cx
    dy = y - cy
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    x_rot = dx * cos_t + dy * sin_t
    y_rot = -dx * sin_t + dy * cos_t
    q = max(1.0 - ellip, 0.1)
    r = np.sqrt(x_rot ** 2 + (y_rot / q) ** 2)
    r = np.maximum(r, 0.01)
    profile = np.exp(-bn * ((r / r_e) ** (1.0 / sersic_n) - 1.0))
    s = profile.sum()
    if s > 0:
        profile /= s
    return profile


def _fft_convolve2d(image, kernel):
    """2D convolution via FFT, same-size output."""
    s1 = np.array(image.shape)
    s2 = np.array(kernel.shape)
    shape = s1 + s2 - 1
    result = np.fft.irfft2(np.fft.rfft2(image, shape) *
                           np.fft.rfft2(kernel, shape), shape)
    cy = (s2[0] - 1) // 2
    cx = (s2[1] - 1) // 2
    return result[cy:cy + s1[0], cx:cx + s1[1]]


def generate_synthetic_galaxies(n=10000, size=64, fwhm_range=(1.5, 5.0),
                                 seed=42):
    """Synthetic galaxies: Sersic profile convolved with CLEAN Moffat PSF.

    NO diffraction spikes — galaxies are extended smooth profiles.
    The PSF convolution ensures the same seeing/optics as stars.
    """
    rng = np.random.RandomState(seed)
    galaxies = np.zeros((n, 1, size, size), dtype=np.float32)
    psf_size = 31

    for i in range(n):
        fwhm = rng.uniform(fwhm_range[0], fwhm_range[1])
        beta = rng.uniform(2.5, 4.5)
        psf_kernel = _make_moffat(psf_size, fwhm, beta)

        # Sersic parameters
        roll = rng.random()
        if roll < 0.35:
            sersic_n = rng.uniform(0.8, 1.2)
        elif roll < 0.60:
            sersic_n = rng.uniform(3.5, 4.5)
        else:
            sersic_n = rng.uniform(0.5, 6.0)

        # Half-light radius: 1.5x to 8x FWHM
        r_e = fwhm * rng.uniform(1.5, 8.0)
        ellip = rng.uniform(0.0, 0.7)
        theta = rng.uniform(0, np.pi)
        cx = size / 2.0 + rng.uniform(-1.5, 1.5)
        cy = size / 2.0 + rng.uniform(-1.5, 1.5)

        galaxy = _make_sersic(size, r_e, sersic_n, ellip, theta, cx, cy)
        galaxy_conv = _fft_convolve2d(galaxy, psf_kernel)

        # Same SNR and noise as stars
        snr = 10.0 ** rng.uniform(0.7, 2.7)
        noise_sigma = rng.uniform(0.005, 0.05)
        peak_flux = snr * noise_sigma

        stamp = galaxy_conv * peak_flux
        stamp += rng.normal(0, noise_sigma, stamp.shape)

        galaxies[i, 0] = _asinh_normalize(stamp)

    return galaxies


# ---------------------------------------------------------------------------
#  Train / val / test split
# ---------------------------------------------------------------------------

def make_train_val_test_split(star_images, galaxy_images,
                               val_frac=0.15, test_frac=0.15, seed=42):
    """Combine and split star/galaxy data."""
    n_star = len(star_images)
    n_gal = len(galaxy_images)

    images = np.concatenate([star_images, galaxy_images], axis=0)
    labels = np.concatenate([
        np.ones(n_star, dtype=np.int64),
        np.zeros(n_gal, dtype=np.int64),
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
