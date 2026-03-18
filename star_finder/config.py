"""Star/Galaxy classification configuration."""

from dataclasses import dataclass


@dataclass
class StarFinderConfig:
    """Inference/cutout configuration."""
    cutout_size: int = 64
    n_channels: int = 1
    batch_size: int = 64
    star_threshold: float = 0.5    # P(star) > threshold -> star


@dataclass
class TrainConfig:
    """Training hyperparameters."""
    epochs: int = 100
    lr: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 64
    val_fraction: float = 0.15
    test_fraction: float = 0.15
    early_stop_patience: int = 15
    seed: int = 42
    cutout_size: int = 64
    checkpoint_file: str = "star_finder_best.pt"
    num_workers: int = 4

    # Synthetic star generation
    n_stars: int = 8000
    n_galaxies: int = 8000
    fwhm_range_min: float = 1.5
    fwhm_range_max: float = 5.0

    @property
    def checkpoint_dir(self):
        import os
        return os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "star_finder", "data", "checkpoints"
        )

    @property
    def checkpoint_path(self):
        import os
        return os.path.join(self.checkpoint_dir, self.checkpoint_file)
