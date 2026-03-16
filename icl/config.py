"""ICL (Intra-Cluster Light) detection pipeline configuration."""

from dataclasses import dataclass


@dataclass
class ICLConfig:
    """Configuration for the ICL detection pipeline."""

    # Source masking
    mask_expand_factor: float = 2.5
    bright_star_mag_limit: float = 18.0
    bright_star_radius_scale: float = 10.0
    interp_method: str = 'linear'

    # Background
    bkg_method: str = 'polynomial'
    bkg_poly_order: int = 3
    bkg_sep_mesh: int = 256
    bkg_sigma_clip: float = 3.0

    # SB Profile
    profile_rmin: float = 5.0
    profile_rmax: float = 1000.0
    profile_nsteps: int = 80
    profile_spacing: str = 'log'
    profile_ellipticity: float = 0.0
    profile_pa: float = 0.0
    profile_sector_width: float = 360.0
    profile_sector_pa: float = 0.0

    # Calibration
    mag_zeropoint: float = 25.0
    pixel_scale: float = 0.06

    # ICL measurement
    icl_mu_threshold: float = 26.5
    icl_mu_levels: str = '26.0,27.0,28.0'
    icl_measure_radius: float = 500.0

    n_workers: int = 0
