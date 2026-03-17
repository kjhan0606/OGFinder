"""PSF Deconvolution pipeline configuration."""

from dataclasses import dataclass


@dataclass
class PSFDeconvConfig:
    """Configuration for the PSF deconvolution pipeline."""

    # Star finding
    class_star_thresh: float = 0.8
    max_ellipticity: float = 0.2
    fwhm_sigma: float = 2.0
    min_flux_snr: float = 10.0

    # PSF
    psf_size: int = 51
    psf_method: str = 'median'

    # Extended PSF
    ext_core_mag_min: float = 18.0
    ext_core_mag_max: float = 22.0
    ext_wing_mag_max: float = 16.0
    ext_core_size: int = 51
    ext_wing_size: int = 201
    ext_blend_inner: float = 20.0
    ext_blend_outer: float = 30.0
    ext_saturation_limit: float = 60000.0
    ext_n_core_min: int = 3
    ext_n_wing_min: int = 1

    # Simulation PSF
    sim_telescope: str = 'auto'
    sim_instrument: str = 'auto'
    sim_filter: str = 'auto'
    sim_psf_size: int = 201
    sim_oversample: int = 1
    sim_jitter_sigma: float = 0.007
    sim_focus_offset: float = 0.0

    # Deconvolution
    deconv_method: str = 'rl'
    rl_iterations: int = 30
    wiener_nsr: float = 0.01
    tikhonov_lambda: float = 0.001
    tv_lambda: float = 0.001
    clean_gain: float = 0.1
    clean_niter: int = 1000
    clean_threshold: float = 0.0
    mem_lambda: float = 0.1
    mem_niter: int = 100
