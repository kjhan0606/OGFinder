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
