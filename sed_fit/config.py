"""SED fitting configuration (Stellar Population Synthesis Emulator)."""

from dataclasses import dataclass
from typing import List

# Physical parameter names for the SPS model
PHYS_PARAM_NAMES = [
    "z",          # redshift
    "log_mass",   # log10(stellar mass / Msun)
    "log_age",    # log10(age / yr)
    "log_Z",      # log10(metallicity / Zsun)
    "Av",         # dust extinction (mag)
    "log_tau",    # log10(e-folding timescale / yr), tau-model SFH
]

# Output parameter names (inverse model predicts these; z is input)
OUTPUT_PARAM_NAMES = [
    "log_mass",
    "log_age",
    "log_Z",
    "Av",
    "log_tau",
]

# Default broadband filter set
DEFAULT_BANDS = [
    "F435W", "F606W", "F775W", "F814W", "F850LP",   # ACS
    "F105W", "F125W", "F140W", "F160W",               # WFC3/IR
]


@dataclass
class SEDConfig:
    """Model architecture and inference configuration."""
    n_phys_params: int = 6          # number of physical input parameters
    n_output_params: int = 5        # inverse model output dims (excl. z)
    n_hidden_emulator: int = 256    # hidden layer width for emulator
    n_hidden_inverse: int = 128     # hidden layer width for inverse model
    n_components: int = 3           # MDN mixture components
    dropout: float = 0.1           # dropout rate
    batch_size: int = 256          # inference batch size


@dataclass
class TrainConfig:
    """Training hyperparameters."""
    epochs: int = 200
    lr: float = 1e-3
    weight_decay: float = 1e-4
    batch_size: int = 256
    val_fraction: float = 0.15
    early_stop_patience: int = 20
    seed: int = 42
    checkpoint_emulator: str = "sps_emulator_best.pt"
    checkpoint_inverse: str = "inverse_sed_best.pt"
    num_workers: int = 4

    @property
    def checkpoint_dir(self):
        import os
        return os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "sed_fit", "data", "checkpoints"
        )

    @property
    def checkpoint_emulator_path(self):
        import os
        return os.path.join(self.checkpoint_dir, self.checkpoint_emulator)

    @property
    def checkpoint_inverse_path(self):
        import os
        return os.path.join(self.checkpoint_dir, self.checkpoint_inverse)
