"""Configuration for non-parametric morphometry."""

from dataclasses import dataclass


@dataclass
class MorphometryConfig:
    min_npix: int = 100          # Minimum pixels for morphometry
    cutout_scale: float = 4.0    # Cutout size = scale * kron_radius * a
    petro_eta: float = 0.2       # Petrosian ratio threshold
    petro_rmax: float = 100.0    # Maximum Petrosian search radius (pixels)
    petro_nsteps: int = 50       # Number of radial steps for Petrosian
    asym_center_maxshift: int = 3  # Max pixel shift for asymmetry minimization
