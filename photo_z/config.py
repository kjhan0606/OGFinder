"""Photometric redshift estimation configuration."""

from dataclasses import dataclass


@dataclass
class PhotoZConfig:
    """Inference configuration for PhotoZ MDN."""
    n_features: int = 20
    n_hidden: int = 128
    n_components: int = 5
    dropout: float = 0.1
    batch_size: int = 256
    quality_head: bool = False  # Enable Quality Head v2 for outlier detection


@dataclass
class TrainConfig:
    """Training hyperparameters."""
    epochs: int = 200
    lr: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 256
    val_fraction: float = 0.15
    test_fraction: float = 0.15
    early_stop_patience: int = 20
    seed: int = 42
    checkpoint_file: str = "photo_z_mdn_best.pt"
    num_workers: int = 4
    quality_head: bool = False  # Train with Quality Head v2
    quality_lambda: float = 1.0  # Weight for quality loss

    @property
    def checkpoint_dir(self):
        import os
        return os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "photo_z", "data", "checkpoints"
        )

    @property
    def checkpoint_path(self):
        import os
        return os.path.join(self.checkpoint_dir, self.checkpoint_file)
