"""External CIGALE program.

``pcigale run`` is a separate program. It is used only when ``pcigale`` is
already on PATH. The pcigale package is not imported. ``generate_sed`` does
not call it. ``get_backend('auto')`` stays on the in-tree analytic model.
"""

import os
import shutil
import subprocess
import tempfile

import numpy as np

from .base import SEDResult, SPSBackend


class CIGALEBackend(SPSBackend):
    @classmethod
    def is_available(cls) -> bool:
        return shutil.which("pcigale") is not None

    @classmethod
    def name(cls) -> str:
        return "cigale"

    def generate_sed(self, z, log_mass, log_age, log_Z, Av, log_tau, bands):
        del z, log_mass, log_age, log_Z, Av, log_tau, bands
        raise RuntimeError(
            "Backend 'cigale' does not call the pcigale package in-process. "
            "The in-tree SED model is analytic."
        )

    def fit_sed(self, mags, mag_errs, photo_z, bands, **kw):
        del kw
        if not self.is_available():
            raise RuntimeError(
                "the external program 'pcigale' was not found on PATH. "
                "A CIGALE fit is not produced here. "
                "The in-tree SED model is analytic."
            )
        return self._fit_subprocess(mags, mag_errs, photo_z, bands)

    def _fit_subprocess(self, mags, mag_errs, photo_z, bands):
        """Run the user-installed ``pcigale run`` on one catalog."""
        from .filters import map_filters

        cig_filters = map_filters(bands, "cigale")
        n_sources = mags.shape[0]

        with tempfile.TemporaryDirectory(prefix="cigale_") as tmpdir:
            cat_path = os.path.join(tmpdir, "input.txt")
            with open(cat_path, "w", encoding="utf-8") as handle:
                cols = ["id", "redshift"]
                for filt in cig_filters:
                    cols.extend([filt, filt + "_err"])
                handle.write(" ".join(cols) + "\n")
                for i in range(n_sources):
                    row = [str(i), f"{photo_z[i]:.4f}"]
                    for j in range(len(bands)):
                        flux = 10 ** (-0.4 * (mags[i, j] - 8.90)) * 1e3
                        err_mag = mag_errs[i, j] if np.isfinite(mag_errs[i, j]) else 0.1
                        err = flux * err_mag * np.log(10) / 2.5
                        if not np.isfinite(flux) or mags[i, j] > 90:
                            row.extend(["-9999", "-9999"])
                        else:
                            row.extend([f"{flux:.4f}", f"{err:.4f}"])
                    handle.write(" ".join(row) + "\n")

            cfg_path = os.path.join(tmpdir, "pcigale.ini")
            with open(cfg_path, "w", encoding="utf-8") as handle:
                handle.write("data_file = input.txt\n")
                handle.write(
                    "sed_modules = sfhdelayed, bc03, "
                    "dustatt_calzleit, redshifting\n"
                )
                handle.write("analysis_method = pdf_analysis\n")

            result = subprocess.run(
                ["pcigale", "run"],
                cwd=tmpdir,
                capture_output=True,
                text=True,
                timeout=3600,
            )
            if result.returncode != 0:
                raise RuntimeError(
                    "pcigale failed (%s). A CIGALE fit is not produced here."
                    % (result.stderr or result.stdout)[:200]
                )

            results_path = os.path.join(tmpdir, "out", "results.txt")
            if not os.path.isfile(results_path):
                raise RuntimeError(
                    "pcigale did not write out/results.txt. "
                    "A CIGALE fit is not produced here."
                )
            parsed = self._parse_cigale_results(results_path, n_sources)
            if not np.any(np.isfinite(parsed.log_mass)):
                raise RuntimeError(
                    "pcigale did not produce a stellar-mass column. "
                    "A CIGALE fit is not produced here."
                )
            return parsed

    def _parse_cigale_results(self, results_path, n_sources):
        """Read the table written by ``pcigale run``."""
        log_mass = np.full(n_sources, np.nan)
        log_mass_err = np.full(n_sources, np.nan)
        log_age = np.full(n_sources, np.nan)
        log_age_err = np.full(n_sources, np.nan)
        log_Z = np.full(n_sources, np.nan)
        Av = np.full(n_sources, np.nan)
        log_tau = np.full(n_sources, np.nan)
        sfr = np.full(n_sources, np.nan)
        chi2 = np.full(n_sources, np.nan)

        with open(results_path, encoding="utf-8") as handle:
            lines = handle.readlines()
        if len(lines) < 2:
            raise RuntimeError("pcigale results.txt has no data rows")

        header = lines[0].strip().split()
        for line in lines[1:]:
            vals = line.strip().split()
            if not vals:
                continue
            idx = int(vals[0])
            if idx >= n_sources:
                continue
            for j, col in enumerate(header):
                if j >= len(vals):
                    continue
                v = float(vals[j])
                low = col.lower()
                if "stellar.m_star" in low and "err" not in low:
                    log_mass[idx] = np.log10(max(v, 1.0))
                elif "stellar.m_star" in low and "err" in low:
                    log_mass_err[idx] = 0.3
                elif "sfh.age_main" in low and "err" not in low:
                    log_age[idx] = np.log10(max(v * 1e6, 1.0))
                elif "attenuation.av" in low:
                    Av[idx] = v
                elif "sfr" in low and "err" not in low:
                    sfr[idx] = np.log10(max(v, 1e-5))
                elif "chi2" in low:
                    chi2[idx] = v

        return SEDResult(
            log_mass=log_mass,
            log_mass_err=log_mass_err,
            log_age=log_age,
            log_age_err=log_age_err,
            log_Z=log_Z,
            Av=Av,
            log_tau=log_tau,
            sfr=sfr,
            chi2=chi2,
        )
