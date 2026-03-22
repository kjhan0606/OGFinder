"""Configuration for 2D Bulge+Disk decomposition."""

from dataclasses import dataclass


@dataclass
class BulgeDiskConfig:
    cutout_scale: float = 5.0      # Cutout size = scale * kron_radius * a
    max_sources: int = 100         # Max sources to decompose
    n_bulge: float = 4.0           # de Vaucouleurs (can be freed)
    n_disk: float = 1.0            # Exponential (fixed)
    free_bulge_n: bool = False     # Allow bulge Sérsic index to vary
    re_min: float = 0.3            # Minimum effective radius (pixels)
    mag_zeropoint: float = 25.0
    pixel_scale: float = 0.263     # arcsec/pixel
