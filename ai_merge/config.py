"""Central configuration for AI merge pipeline."""

from dataclasses import dataclass, field
from typing import List
import os


@dataclass
class SimConfig:
    """Simulation parameters for pair generation."""
    # Galaxy profile ranges
    sersic_n_choices: List[float] = field(default_factory=lambda: [0.5, 1.0, 2.0, 4.0])
    re_range: tuple = (2.0, 20.0)       # effective radius in pixels
    flux_range: tuple = (500.0, 50000.0)
    ellip_range: tuple = (0.0, 0.7)

    # Separation (in units of R1+R2)
    sep_ratio_range: tuple = (0.3, 2.5)

    # PSF
    psf_fwhm_choices: List[float] = field(default_factory=lambda: [1.5, 3.0, 5.0])

    # Noise
    noise_std_range: tuple = (1.0, 20.0)

    # Stamp size
    stamp_size: int = 128

    # Counts
    n_merge: int = 25000
    n_separate: int = 25000

    # Merge-specific: secondary feature relative to primary
    merge_flux_frac_range: tuple = (0.05, 0.50)   # 5-50% of primary flux
    merge_size_frac_range: tuple = (0.10, 0.40)   # 10-40% of primary re
    merge_offset_frac: float = 1.5                  # within 1.5 * re of primary

    # SEP detection params
    detect_thresh: float = 1.5
    detect_minarea: int = 5
    deblend_nthresh: int = 32
    deblend_mincont: float = 0.005

    # Mag zeropoint
    mag_zeropoint: float = 25.0

    # Max retries for rejection sampling
    max_retries: int = 50

    # Output
    output_dir: str = "SAOImageDS9/ai_merge/data/simulated"
    output_file: str = "pairs_v1.h5"

    @property
    def output_path(self):
        return os.path.join(self.output_dir, self.output_file)


@dataclass
class TrainConfig:
    """Training hyperparameters."""
    batch_size: int = 256
    lr: float = 1e-3
    weight_decay: float = 1e-4
    epochs: int = 100
    patience: int = 15
    n_classes: int = 2
    val_fraction: float = 0.2
    n_features: int = 20
    hidden_dims: List[int] = field(default_factory=lambda: [128, 64, 32])
    dropout_rates: List[float] = field(default_factory=lambda: [0.3, 0.2, 0.0])
    checkpoint_dir: str = "SAOImageDS9/ai_merge/data/checkpoints"
    checkpoint_file: str = "mlp_best.pt"
    seed: int = 42

    @property
    def checkpoint_path(self):
        return os.path.join(self.checkpoint_dir, self.checkpoint_file)


@dataclass
class CandidateConfig:
    """Candidate pair search parameters."""
    alpha: float = 1.5           # distance factor: d < alpha * (R_i + R_j)
    radius_scale: float = 2.0   # radius = radius_scale * kron_radius * a
    min_flux: float = 0.0       # minimum flux for candidates
