"""LSBG (Low Surface Brightness Galaxy) detection pipeline configuration."""

from dataclasses import dataclass


@dataclass
class LSBGConfig:
    """Configuration for the LSBG detection pipeline."""

    # Masking
    mask_detect_thresh: float = 1.5
    mask_detect_minarea: int = 5
    mask_expand_factor: float = 3.0
    bright_star_mag_limit: float = 18.0
    bright_star_radius_scale: float = 12.0
    mask_mag_threshold: float = 22.0
    interp_method: str = 'linear'

    # Iterative background
    bkg_method: str = 'sep_large'
    bkg_mesh_size: int = 256
    bkg_poly_order: int = 3
    bkg_sigma_clip: float = 3.0
    bkg_n_iterations: int = 3
    bkg_refine_thresh: float = 2.0

    # Detection (cleaned image)
    detect_thresh: float = 0.8
    detect_minarea: int = 50
    detect_filter_kernel: str = 'gauss5x5'
    deblend_nthresh: int = 32
    deblend_mincont: float = 0.005

    # Photometry
    phot_apertures: str = '5,10,20,40'
    mag_zeropoint: float = 25.0
    pixel_scale: float = 0.06

    # LSBG Filtering
    mu_eff_min: float = 24.0
    mu_eff_max: float = 30.0
    r_eff_min: float = 2.5
    r_eff_max: float = 60.0
    ellipticity_max: float = 0.7
    min_snr: float = 2.0

    n_workers: int = 0
