"""SPS grid generation using a simple analytic SED model.

Generates synthetic broadband photometry from random physical parameters
using a power-law continuum + Calzetti dust extinction + redshift dimming.
No dependency on FSPS or other SPS libraries required.
"""

import os
import argparse
import numpy as np

from ..config import PHYS_PARAM_NAMES, DEFAULT_BANDS


# -----------------------------------------------------------------------
# Approximate effective wavelengths (Angstrom) for common HST filters
# -----------------------------------------------------------------------
FILTER_WAVELENGTHS = {
    # ACS/WFC
    "F435W":  4328.0,
    "F475W":  4746.0,
    "F555W":  5346.0,
    "F606W":  5907.0,
    "F625W":  6311.0,
    "F775W":  7693.0,
    "F814W":  8057.0,
    "F850LP": 9033.0,
    # WFC3/IR
    "F098M": 9864.0,
    "F105W": 10552.0,
    "F110W": 11534.0,
    "F125W": 12486.0,
    "F140W": 13923.0,
    "F160W": 15369.0,
    # WFC3/UVIS
    "F225W": 2372.0,
    "F275W": 2710.0,
    "F336W": 3355.0,
    "F390W": 3923.0,
    # JWST/NIRCam
    "F070W":  7088.0,
    "F090W":  9010.0,
    "F115W": 11540.0,
    "F150W": 15010.0,
    "F200W": 19890.0,
    "F277W": 27620.0,
    "F356W": 35680.0,
    "F444W": 44050.0,
}


def calzetti_extinction(wavelength_angstrom, Av):
    """
    Calzetti (2000) dust extinction curve.

    Parameters
    ----------
    wavelength_angstrom : float or array
        Rest-frame wavelength in Angstrom.
    Av : float or array
        V-band extinction in magnitudes.

    Returns
    -------
    A_lambda : float or array
        Extinction in magnitudes at the given wavelength.
    """
    lam_um = np.asarray(wavelength_angstrom, dtype=np.float64) / 1e4  # um
    Rv = 4.05

    # Piecewise k'(lambda)
    k = np.where(
        lam_um < 0.63,
        2.659 * (-2.156 + 1.509 / lam_um - 0.198 / lam_um**2
                 + 0.011 / lam_um**3) + Rv,
        2.659 * (-1.857 + 1.040 / lam_um) + Rv,
    )

    # Clip to avoid negative extinction at very red wavelengths
    k = np.maximum(k, 0.0)

    return Av * k / Rv


def simple_sed_model(z, log_mass, log_age, log_Z, Av, log_tau,
                     wavelengths_angstrom):
    """
    Simple analytic SED model for grid generation.

    Model components:
    1. Power-law continuum: F_nu ~ lambda^beta, with beta depending
       on age and metallicity (bluer for young/metal-poor populations).
    2. Calzetti dust extinction.
    3. Cosmological dimming and redshift.

    Parameters
    ----------
    z : float
        Redshift.
    log_mass : float
        log10(stellar mass / Msun).
    log_age : float
        log10(age / yr).
    log_Z : float
        log10(metallicity / Zsun).
    Av : float
        V-band dust extinction (mag).
    log_tau : float
        log10(e-folding timescale / yr). Affects UV slope.
    wavelengths_angstrom : array
        Observer-frame effective wavelengths in Angstrom.

    Returns
    -------
    mags : np.ndarray
        Apparent AB magnitudes at each wavelength.
    """
    wavelengths = np.asarray(wavelengths_angstrom, dtype=np.float64)

    # Rest-frame wavelengths
    lam_rest = wavelengths / (1.0 + z)

    # UV slope: younger -> bluer, higher Z -> redder, shorter tau -> bluer
    age_yr = 10.0 ** log_age
    tau_yr = 10.0 ** log_tau
    Z_solar = 10.0 ** log_Z

    # Approximate UV spectral slope beta (F_lambda ~ lambda^beta)
    # Young, metal-poor, short-tau: steep blue (beta ~ -2.5)
    # Old, metal-rich, long-tau: flat/red (beta ~ 0)
    age_factor = np.clip((log_age - 7.0) / 3.0, 0.0, 1.0)   # 0 at 10Myr, 1 at 10Gyr
    Z_factor = np.clip((log_Z + 2.0) / 2.0, 0.0, 1.0)        # 0 at 0.01 Zsun, 1 at Zsun
    tau_factor = np.clip((log_tau - 7.0) / 3.0, 0.0, 1.0)     # 0 at 10Myr, 1 at 10Gyr

    beta = -2.5 + 2.0 * age_factor + 0.5 * Z_factor + 0.3 * tau_factor

    # Reference wavelength (5500 A rest-frame)
    lam_ref = 5500.0

    # Continuum flux (F_lambda ~ lambda^beta, arbitrary normalization)
    f_lambda = (lam_rest / lam_ref) ** beta

    # Apply dust extinction
    A_lambda = calzetti_extinction(lam_rest, Av)
    f_lambda *= 10.0 ** (-0.4 * A_lambda)

    # Luminosity distance (simple Hubble law for demonstration)
    # d_L ~ c*z/H0 for z < ~0.3; for higher z use a rough approximation
    if z < 0.01:
        z_eff = 0.01  # avoid log(0)
    else:
        z_eff = z
    # Rough d_L in Mpc: includes a (1+z) factor for approximate correction
    d_L_approx = (z_eff * 3e5 / 70.0) * (1.0 + 0.5 * z_eff)  # Mpc

    # Distance modulus
    DM = 5.0 * np.log10(np.maximum(d_L_approx, 1e-3)) + 25.0

    # Apparent magnitude calibration:
    # A 1e10 Msun galaxy at 10 pc would have M_AB ~ -20 (V-band).
    # m_AB = M_AB + DM = -20 + DM
    # For log_mass=10 (1e10 Msun): M_abs ~ -20 at reference wavelength
    # For other masses: scale by mass ratio -> delta_M = -2.5*(log_mass - 10)
    # f_lambda provides color variation relative to reference wavelength
    M_abs_ref = -20.0  # absolute mag of 1e10 Msun galaxy at lam_ref
    delta_M_mass = -2.5 * (log_mass - 10.0)  # mass scaling
    delta_M_color = -2.5 * np.log10(np.maximum(f_lambda, 1e-30))  # color term

    mags = M_abs_ref + delta_M_mass + delta_M_color + DM

    return mags


def generate_sps_grid(n_samples, bands, output_path, seed=42):
    """
    Generate a grid of synthetic SPS photometry.

    Parameters
    ----------
    n_samples : int
        Number of random parameter draws.
    bands : list of str
        Filter names (must be keys in FILTER_WAVELENGTHS).
    output_path : str
        Output HDF5 file path.
    seed : int
        Random seed.

    Returns
    -------
    params : np.ndarray
        Shape (n_samples, 6). Physical parameters.
    mags : np.ndarray
        Shape (n_samples, n_bands). Synthetic magnitudes.
    """
    rng = np.random.RandomState(seed)

    # Resolve wavelengths
    wavelengths = np.array([FILTER_WAVELENGTHS[b] for b in bands])
    n_bands = len(bands)

    # Parameter ranges (uniform in these ranges)
    ranges = {
        'z':        (0.01, 6.0),
        'log_mass': (6.0, 12.0),     # log10(Msun)
        'log_age':  (7.0, 10.15),    # log10(yr): 10 Myr to 14 Gyr
        'log_Z':    (-2.0, 0.3),     # log10(Z/Zsun): 0.01 to 2 Zsun
        'Av':       (0.0, 3.0),      # mag
        'log_tau':  (7.0, 10.5),     # log10(yr): 10 Myr to 30 Gyr
    }

    # Draw random parameters
    params = np.zeros((n_samples, 6), dtype=np.float32)
    for i, name in enumerate(PHYS_PARAM_NAMES):
        lo, hi = ranges[name]
        params[:, i] = rng.uniform(lo, hi, size=n_samples).astype(np.float32)

    # Compute magnitudes
    mags = np.zeros((n_samples, n_bands), dtype=np.float32)
    for j in range(n_samples):
        mags[j] = simple_sed_model(
            z=params[j, 0],
            log_mass=params[j, 1],
            log_age=params[j, 2],
            log_Z=params[j, 3],
            Av=params[j, 4],
            log_tau=params[j, 5],
            wavelengths_angstrom=wavelengths,
        )

    # Add photometric noise (0.01-0.1 mag, larger for fainter sources)
    noise_scale = 0.01 + 0.05 * np.clip((mags - 20.0) / 10.0, 0.0, 1.0)
    mags += rng.normal(0, noise_scale).astype(np.float32)

    # Save to HDF5
    import h5py
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with h5py.File(output_path, 'w') as f:
        f.create_dataset('params', data=params)
        f.create_dataset('mags', data=mags)
        f.attrs['param_names'] = PHYS_PARAM_NAMES
        f.attrs['band_names'] = bands
        f.attrs['n_samples'] = n_samples
        f.attrs['n_bands'] = n_bands
        f.attrs['wavelengths'] = wavelengths

    print(f"Generated {n_samples} SPS models with {n_bands} bands")
    print(f"Saved to {output_path}")
    print(f"Parameter ranges:")
    for name in PHYS_PARAM_NAMES:
        lo, hi = ranges[name]
        print(f"  {name:10s}: [{lo:.2f}, {hi:.2f}]")

    return params, mags


def main():
    """CLI entry point for SPS grid generation."""
    parser = argparse.ArgumentParser(
        description='Generate SPS grid for SED fitting emulator training')
    parser.add_argument('--n-samples', type=int, default=100000,
                        help='Number of SPS models to generate (default: 100000)')
    parser.add_argument('--bands', nargs='+', default=None,
                        help='Filter names (default: ACS+WFC3/IR set)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output HDF5 path')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed (default: 42)')
    args = parser.parse_args()

    bands = args.bands if args.bands else DEFAULT_BANDS

    if args.output is None:
        data_dir = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "data"
        )
        output_path = os.path.join(data_dir, "sps_grid.h5")
    else:
        output_path = args.output

    generate_sps_grid(args.n_samples, bands, output_path, seed=args.seed)


if __name__ == '__main__':
    main()
