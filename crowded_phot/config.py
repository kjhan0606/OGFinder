"""Configuration for crowded field photometry."""

from dataclasses import dataclass


@dataclass
class CrowdedPhotConfig:
    group_radius_scale: float = 2.0  # Grouping radius = scale * kron_radius
    max_group_size: int = 20         # Max sources per group
    max_iterations: int = 3          # Iteration cycles
    fit_radius: int = 10             # PSF fitting radius
    detect_thresh: float = 3.0       # Threshold for residual detection
    detect_minarea: int = 3          # Min area for residual detection
    mag_zeropoint: float = 25.0
    max_shift: float = 2.0           # Max centroid shift
