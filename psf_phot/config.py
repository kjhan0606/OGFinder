"""Configuration for PSF photometry."""

from dataclasses import dataclass


@dataclass
class PSFPhotConfig:
    fit_radius: int = 10        # Fitting radius in pixels
    mag_zeropoint: float = 25.0
    max_shift: float = 3.0      # Maximum centroid shift (pixels)
    bg_annulus_inner: float = 15.0  # Background annulus inner radius
    bg_annulus_outer: float = 25.0  # Background annulus outer radius
