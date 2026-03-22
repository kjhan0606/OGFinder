"""CIGALE (pcigale) backend."""

import numpy as np
from .base import SPSBackend, SEDResult


class CIGALEBackend(SPSBackend):
    """SED fitting using CIGALE (Code Investigating GALaxy Emission).

    Requires: pip install pcigale
    Uses subprocess for fitting (pcigale run), direct API for SED generation.
    """

    @classmethod
    def is_available(cls) -> bool:
        try:
            import pcigale
            return True
        except ImportError:
            return False

    @classmethod
    def name(cls) -> str:
        return "cigale"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau,
                     bands) -> np.ndarray:
        from .filters import map_filters
        from .params import canonical_to_cigale

        params = canonical_to_cigale(z, log_mass, log_age, log_Z,
                                      Av, log_tau)
        cig_filters = map_filters(bands, 'cigale')

        from pcigale.sed_modules import get_module

        # Build SED through module pipeline
        sfh = get_module('sfhdelayed')
        sfh_params = {
            'age_main': params['sfhdelayed']['age_main'],
            'tau_main': params['sfhdelayed']['tau_main'],
            'age_burst': 20.0,
            'tau_burst': 50.0,
            'f_burst': 0.0,
            'sfr_A': 1.0,
            'normalise': True,
        }
        sfh.process(sfh_params)
        sed = sfh.sed

        ssp = get_module('bc03')
        ssp_params = {
            'imf': params['bc03']['imf'],
            'metallicity': params['bc03']['metallicity'],
            'separation_age': 10,
        }
        ssp.process(ssp_params, sed)

        dust = get_module('dustatt_calzleit')
        dust_params = {
            'Av_BC': params['dustatt_calzleit']['Av'],
            'Av_ISM': params['dustatt_calzleit']['Av'] * 0.3,
            'slope_BC': 0.0,
            'slope_ISM': 0.0,
        }
        dust.process(dust_params, sed)

        redshift = get_module('redshifting')
        redshift.process({'redshift': z}, sed)

        # Compute photometry
        mags_out = np.zeros(len(bands), dtype=np.float32)
        for j, filt_name in enumerate(cig_filters):
            try:
                flux = sed.compute_fnu(filt_name)
                if flux > 0:
                    mags_out[j] = -2.5 * np.log10(flux) + 8.90
                else:
                    mags_out[j] = 99.0
            except Exception:
                mags_out[j] = 99.0

        # Scale by mass
        mags_out += -2.5 * (log_mass - 10.0)

        return mags_out

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw) -> SEDResult:
        """Fit using CIGALE subprocess pipeline."""
        import sys

        try:
            return self._fit_subprocess(mags, mag_errs, photo_z, bands, **kw)
        except Exception as e:
            print(f"WARNING: CIGALE fit failed ({e}), "
                  f"falling back to chi2", file=sys.stderr)
            return self.fit_sed_chi2(mags, mag_errs, photo_z, bands,
                                     n_grid=kw.get('n_grid', 5000))

    def _fit_subprocess(self, mags, mag_errs, photo_z, bands, **kw):
        """Run pcigale via subprocess for batch fitting."""
        import sys
        import os
        import tempfile
        import subprocess
        from .filters import map_filters

        cig_filters = map_filters(bands, 'cigale')
        n_sources = mags.shape[0]

        # Create temp directory for CIGALE run
        with tempfile.TemporaryDirectory(prefix='cigale_') as tmpdir:
            # Write input catalog
            cat_path = os.path.join(tmpdir, 'input.txt')
            with open(cat_path, 'w') as f:
                cols = ['id', 'redshift']
                for filt in cig_filters:
                    cols.extend([filt, filt + '_err'])
                f.write(' '.join(cols) + '\n')

                for i in range(n_sources):
                    row = [str(i), f"{photo_z[i]:.4f}"]
                    for j in range(len(bands)):
                        flux = 10 ** (-0.4 * (mags[i, j] - 8.90)) * 1e3
                        err = flux * (mag_errs[i, j] if np.isfinite(
                            mag_errs[i, j]) else 0.1) * np.log(10) / 2.5
                        if not np.isfinite(flux) or mags[i, j] > 90:
                            row.extend(['-9999', '-9999'])
                        else:
                            row.extend([f"{flux:.4f}", f"{err:.4f}"])
                    f.write(' '.join(row) + '\n')

            # Write CIGALE config
            cfg_path = os.path.join(tmpdir, 'pcigale.ini')
            with open(cfg_path, 'w') as f:
                f.write(f"data_file = input.txt\n")
                f.write(f"sed_modules = sfhdelayed, bc03, "
                        f"dustatt_calzleit, redshifting\n")
                f.write(f"analysis_method = pdf_analysis\n")

            # Run CIGALE
            result = subprocess.run(
                ['pcigale', 'run'],
                cwd=tmpdir,
                capture_output=True, text=True, timeout=3600,
            )

            if result.returncode != 0:
                raise RuntimeError(f"CIGALE failed: {result.stderr[:200]}")

            # Parse results
            results_path = os.path.join(tmpdir, 'out', 'results.txt')
            return self._parse_cigale_results(results_path, n_sources)

    def _parse_cigale_results(self, results_path, n_sources):
        """Parse CIGALE output results file."""
        log_mass = np.full(n_sources, np.nan)
        log_mass_err = np.full(n_sources, np.nan)
        log_age = np.full(n_sources, np.nan)
        log_age_err = np.full(n_sources, np.nan)
        log_Z = np.full(n_sources, np.nan)
        Av = np.full(n_sources, np.nan)
        log_tau = np.full(n_sources, np.nan)
        sfr = np.full(n_sources, np.nan)
        chi2 = np.full(n_sources, np.nan)

        try:
            with open(results_path, 'r') as f:
                lines = f.readlines()
            if len(lines) < 2:
                raise RuntimeError("Empty results")

            header = lines[0].strip().split()
            for line in lines[1:]:
                vals = line.strip().split()
                if not vals:
                    continue
                idx = int(vals[0])
                if idx >= n_sources:
                    continue

                for j, col in enumerate(header):
                    v = float(vals[j]) if j < len(vals) else np.nan
                    if 'stellar.m_star' in col.lower() and 'err' not in col:
                        log_mass[idx] = np.log10(max(v, 1))
                    elif 'stellar.m_star' in col.lower() and 'err' in col:
                        log_mass_err[idx] = 0.3
                    elif 'sfh.age_main' in col.lower() and 'err' not in col:
                        log_age[idx] = np.log10(max(v * 1e6, 1)) + 0
                    elif 'attenuation.Av' in col.lower():
                        Av[idx] = v
                    elif 'sfr' in col.lower() and 'err' not in col:
                        sfr[idx] = np.log10(max(v, 1e-5))
                    elif 'chi2' in col.lower():
                        chi2[idx] = v
        except Exception:
            pass

        return SEDResult(
            log_mass=log_mass, log_mass_err=log_mass_err,
            log_age=log_age, log_age_err=log_age_err,
            log_Z=log_Z, Av=Av, log_tau=log_tau,
            sfr=sfr, chi2=chi2,
        )
