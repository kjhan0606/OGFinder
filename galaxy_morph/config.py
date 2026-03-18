"""Galaxy morphology classification configuration."""

from dataclasses import dataclass, field
from typing import List

# Galaxy10 DECaLS 10-class categories
GALAXY10_CLASSES = [
    "Disturbed",
    "Merging",
    "Round Smooth",
    "In-between Smooth",
    "Cigar Smooth",
    "Barred Spiral",
    "Unbarred Tight Spiral",
    "Unbarred Loose Spiral",
    "Edge-on w/o Bulge",
    "Edge-on w/ Bulge",
]

# DS9 marker colors per class
MORPH_COLORS = [
    "red",       # 0: Disturbed
    "magenta",   # 1: Merging
    "cyan",      # 2: Round Smooth
    "blue",      # 3: In-between Smooth
    "green",     # 4: Cigar Smooth
    "yellow",    # 5: Barred Spiral
    "orange",    # 6: Unbarred Tight Spiral
    "white",     # 7: Unbarred Loose Spiral
    "#00ff00",   # 8: Edge-on w/o Bulge
    "#ff69b4",   # 9: Edge-on w/ Bulge
]


@dataclass
class MorphConfig:
    """Inference/cutout configuration."""
    cutout_size: int = 64
    n_classes: int = 10
    min_confidence: float = 0.3
    min_a_image: float = 2.0
    max_class_star: float = 0.8
    asinh_a: float = 0.1


@dataclass
class TrainConfig:
    """Training hyperparameters."""
    batch_size: int = 64
    lr: float = 1e-3
    weight_decay: float = 1e-4
    epochs: int = 100
    patience: int = 15
    dropout: float = 0.3
    dropout2: float = 0.2
    n_classes: int = 10
    cutout_size: int = 64
    val_fraction: float = 0.15
    test_fraction: float = 0.05
    seed: int = 42
    checkpoint_file: str = "cnn_morph_best.pt"
    hidden_dim: int = 64
    num_workers: int = 4

    @property
    def checkpoint_dir(self):
        import os
        return os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "galaxy_morph", "data", "checkpoints"
        )

    @property
    def checkpoint_path(self):
        import os
        return os.path.join(self.checkpoint_dir, self.checkpoint_file)
