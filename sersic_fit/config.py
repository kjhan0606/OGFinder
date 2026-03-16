"""Configuration for Sérsic profile fitting."""

from dataclasses import dataclass


@dataclass
class SersicConfig:
    cutout_scale: float = 5.0     # Cutout size = scale * kron_radius * a
    max_sources: int = 200        # Max sources to fit
    n_min: float = 0.2            # Minimum Sérsic index
    n_max: float = 10.0           # Maximum Sérsic index
    re_min: float = 0.5           # Minimum effective radius (pixels)
    mag_zeropoint: float = 25.0
