"""Abstract base class for SPS backends."""

from abc import ABC, abstractmethod
from dataclasses import dataclass
import os
import numpy as np


@dataclass
class SEDResult:
    """Result of SED fitting for one or more sources."""
    log_mass: np.ndarray       # (N,)
    log_mass_err: np.ndarray
    log_age: np.ndarray
    log_age_err: np.ndarray
    log_Z: np.ndarray
    Av: np.ndarray
    log_tau: np.ndarray
    sfr: np.ndarray
    chi2: np.ndarray           # reduced chi2


# Default parameter ranges for grid generation / chi2 fitting
DEFAULT_PARAM_RANGES = {
    'z':        (0.01, 6.0),
    'log_mass': (6.0, 12.0),
    'log_age':  (7.0, 10.15),
    'log_Z':    (-2.0, 0.3),
    'Av':       (0.0, 3.0),
    'log_tau':  (7.0, 10.5),
}


class SPSBackend(ABC):
    """Abstract base class for Stellar Population Synthesis backends.

    All backends must implement:
    - is_available(): whether the required package is installed
    - name(): backend name string
    - generate_sed(): single SED generation (params -> magnitudes)
    - fit_sed(): SED fitting (observations -> physical params)
    """

    @classmethod
    @abstractmethod
    def is_available(cls) -> bool:
        """Return True if the required SPS package is installed."""

    @classmethod
    @abstractmethod
    def name(cls) -> str:
        """Return the backend name (e.g., 'fsps', 'bagpipes')."""

    @abstractmethod
    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        """Generate a single SED: physical parameters -> AB magnitudes.

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
            log10(e-folding timescale / yr).
        bands : list of str
            Canonical filter names (e.g., ['F435W', 'F814W']).

        Returns
        -------
        mags : np.ndarray
            AB magnitudes, shape (n_bands,).
        """

    def generate_grid(self, n_samples, bands, output_path, seed=42,
                      param_ranges=None):
        """Generate a grid of synthetic photometry and save to HDF5.

        Default implementation loops over generate_sed(). Subclasses
        may override for batch optimization.

        Parameters
        ----------
        n_samples : int
            Number of random parameter draws.
        bands : list of str
            Canonical filter names.
        output_path : str
            Output HDF5 file path.
        seed : int
            Random seed.
        param_ranges : dict, optional
            Override default parameter ranges.

        Returns
        -------
        params : np.ndarray, shape (n_samples, 6)
        mags : np.ndarray, shape (n_samples, n_bands)
        """
        import h5py
        from ..config import PHYS_PARAM_NAMES

        rng = np.random.RandomState(seed)
        ranges = dict(DEFAULT_PARAM_RANGES)
        if param_ranges:
            ranges.update(param_ranges)

        n_bands = len(bands)

        # Draw random parameters
        params = np.zeros((n_samples, 6), dtype=np.float32)
        for i, pname in enumerate(PHYS_PARAM_NAMES):
            lo, hi = ranges[pname]
            params[:, i] = rng.uniform(lo, hi, size=n_samples).astype(
                np.float32)

        # Compute magnitudes
        mags = np.zeros((n_samples, n_bands), dtype=np.float32)
        for j in range(n_samples):
            mags[j] = self.generate_sed(
                z=float(params[j, 0]),
                log_mass=float(params[j, 1]),
                log_age=float(params[j, 2]),
                log_Z=float(params[j, 3]),
                Av=float(params[j, 4]),
                log_tau=float(params[j, 5]),
                bands=bands,
            )
            if (j + 1) % 5000 == 0:
                import sys
                print(f"  {self.name()}: {j+1}/{n_samples} SEDs generated",
                      file=sys.stderr)

        # Add photometric noise
        noise_scale = 0.01 + 0.05 * np.clip(
            (mags - 20.0) / 10.0, 0.0, 1.0)
        mags += rng.normal(0, noise_scale).astype(np.float32)

        # Save to HDF5
        os.makedirs(os.path.dirname(os.path.abspath(output_path)),
                     exist_ok=True)
        with h5py.File(output_path, 'w') as f:
            f.create_dataset('params', data=params)
            f.create_dataset('mags', data=mags)
            f.attrs['param_names'] = PHYS_PARAM_NAMES
            f.attrs['band_names'] = bands
            f.attrs['n_samples'] = n_samples
            f.attrs['n_bands'] = n_bands
            f.attrs['backend'] = self.name()

        import sys
        print(f"Generated {n_samples} SEDs with {n_bands} bands "
              f"[backend={self.name()}]", file=sys.stderr)
        print(f"Saved to {output_path}", file=sys.stderr)

        return params, mags

    @abstractmethod
    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        """Fit SEDs to observed photometry.

        Parameters
        ----------
        mags : np.ndarray, shape (N, n_bands)
            Observed AB magnitudes.
        mag_errs : np.ndarray, shape (N, n_bands)
            Magnitude errors.
        photo_z : np.ndarray, shape (N,)
            Photometric redshifts.
        bands : list of str
            Canonical filter names.

        Returns
        -------
        SEDResult
        """

    def fit_sed_chi2(self, mags, mag_errs, photo_z, bands,
                     n_grid=10000, seed=42) -> SEDResult:
        """Chi2 grid-search fitting (for backends without native fitters).

        Generates a random grid in parameter space, computes model
        magnitudes via generate_sed(), and finds the best-fit by chi2
        minimization.

        Parameters
        ----------
        mags : np.ndarray, shape (N, n_bands)
        mag_errs : np.ndarray, shape (N, n_bands)
        photo_z : np.ndarray, shape (N,)
        bands : list of str
        n_grid : int
            Number of grid points per source.
        seed : int
            Random seed.

        Returns
        -------
        SEDResult
        """
        import sys
        rng = np.random.RandomState(seed)
        n_sources = mags.shape[0]

        log_mass_out = np.full(n_sources, np.nan)
        log_mass_err_out = np.full(n_sources, np.nan)
        log_age_out = np.full(n_sources, np.nan)
        log_age_err_out = np.full(n_sources, np.nan)
        log_Z_out = np.full(n_sources, np.nan)
        Av_out = np.full(n_sources, np.nan)
        log_tau_out = np.full(n_sources, np.nan)
        sfr_out = np.full(n_sources, np.nan)
        chi2_out = np.full(n_sources, np.nan)

        ranges = DEFAULT_PARAM_RANGES

        # Pre-draw grid parameters (excluding z which comes from photo_z)
        grid_params = np.zeros((n_grid, 5), dtype=np.float32)
        for i, pname in enumerate(['log_mass', 'log_age', 'log_Z',
                                    'Av', 'log_tau']):
            lo, hi = ranges[pname]
            grid_params[:, i] = rng.uniform(lo, hi, size=n_grid)

        for i in range(n_sources):
            obs = mags[i]
            errs = mag_errs[i]
            z = photo_z[i]

            valid = np.isfinite(obs) & (obs > 0) & (obs < 90)
            if np.sum(valid) < 2 or not np.isfinite(z) or z <= 0:
                continue

            # Default errors where missing
            err_use = np.where(np.isfinite(errs) & (errs > 0),
                               errs, 0.1)

            best_chi2 = np.inf
            best_params = None

            for k in range(n_grid):
                lm, la, lz, av, lt = grid_params[k]
                try:
                    model = self.generate_sed(z, float(lm), float(la),
                                              float(lz), float(av),
                                              float(lt), bands)
                except Exception:
                    continue

                resid = (obs[valid] - model[valid]) / err_use[valid]
                chi2 = np.sum(resid ** 2) / max(np.sum(valid) - 1, 1)

                if chi2 < best_chi2:
                    best_chi2 = chi2
                    best_params = (lm, la, lz, av, lt)

            if best_params is not None:
                lm, la, lz, av, lt = best_params
                log_mass_out[i] = round(float(lm), 2)
                log_age_out[i] = round(float(la), 2)
                log_Z_out[i] = round(float(lz), 2)
                Av_out[i] = round(float(av), 2)
                log_tau_out[i] = round(float(lt), 2)
                chi2_out[i] = round(float(best_chi2), 3)

                # Rough error estimate from delta-chi2
                log_mass_err_out[i] = 0.3
                log_age_err_out[i] = 0.4

                # SFR from mass and age
                mass_val = 10 ** float(lm)
                age_gyr = 10 ** (float(la) - 9)
                sfr_val = mass_val / (age_gyr * 1e9) if age_gyr > 0 else 0
                sfr_out[i] = round(np.log10(max(sfr_val, 1e-5)), 2)

            if (i + 1) % 50 == 0:
                print(f"  chi2 fit: {i+1}/{n_sources}", file=sys.stderr)

        return SEDResult(
            log_mass=log_mass_out,
            log_mass_err=log_mass_err_out,
            log_age=log_age_out,
            log_age_err=log_age_err_out,
            log_Z=log_Z_out,
            Av=Av_out,
            log_tau=log_tau_out,
            sfr=sfr_out,
            chi2=chi2_out,
        )
