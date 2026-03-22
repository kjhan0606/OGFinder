#!/usr/bin/env python3
"""
DS9 SED Fitting (AI SPS Emulator + Multi-Backend) CLI.

Uses Inverse MLP+MDN to estimate galaxy physical properties from
multi-band photometry + photo-z. Falls back to empirical color-mass
relation if no checkpoint available.

Supports multiple SPS backends (FSPS, Bagpipes, Prospector, CIGALE,
Dense Basis) for direct SED fitting when --backend is specified.

Usage:
    ds9_sed_fit.py FITSFILE --catalog catalog.tsv
                   --bands g,r,i,z --mag-columns MAG_g,MAG_r,MAG_i,MAG_z
                   [--photo-z-column PHOTO_Z]
                   [--checkpoint-emulator sps_emulator_best.pt]
                   [--checkpoint-inverse inverse_sed_best.pt]
                   [--backend auto|analytic|fsps|bagpipes|prospector|cigale|dense_basis]
                   [--mode fit|list-backends]
                   [--n-workers 4]
"""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def parse_catalog_tsv(path):
    """Parse TSV catalog."""
    with open(path, 'r') as f:
        lines = f.readlines()
    if not lines:
        return None, None

    header = lines[0].strip().split('\t')
    rows = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        rows.append(line.split('\t'))

    return header, rows


def get_column_index(header, name):
    """Find column index by name."""
    for i, h in enumerate(header):
        if h.strip() == name:
            return i
    return -1


def estimate_sed_empirical(mags, mag_errs, photo_z, band_names):
    """Empirical SED fitting from colors (no model needed).

    Uses color-mass, color-age relations calibrated on SDSS.
    """
    n_sources = mags.shape[0]
    n_bands = mags.shape[1]

    log_mass = np.full(n_sources, np.nan)
    log_mass_err = np.full(n_sources, np.nan)
    log_age = np.full(n_sources, np.nan)
    log_age_err = np.full(n_sources, np.nan)
    log_z_met = np.full(n_sources, np.nan)
    av = np.full(n_sources, np.nan)
    sfr = np.full(n_sources, np.nan)

    for i in range(n_sources):
        m = mags[i]
        valid = np.isfinite(m) & (m < 90) & (m > 0)

        if np.sum(valid) < 2:
            continue

        valid_mags = m[valid]
        z = photo_z[i] if np.isfinite(photo_z[i]) else 0.1

        # Color indices
        colors = []
        for j in range(n_bands - 1):
            if valid[j] and valid[j + 1]:
                colors.append(m[j] - m[j + 1])

        if not colors:
            continue

        mean_color = np.mean(colors)
        reddest_color = max(colors) if colors else 0

        # Luminosity distance modulus (simplified)
        dm = 5 * np.log10(max(z, 0.001) * 4285.7) + 25  # Mpc, H0=70

        # Absolute magnitude (use reddest available band)
        m_abs = valid_mags[-1] - dm

        # Mass-to-light from color (Bell+2003 relation)
        # log(M/L) ~ -0.306 + 1.097 * (g-r) for r-band
        ml_ratio = -0.306 + 1.097 * max(0, min(mean_color, 1.5))
        log_lum = -0.4 * (m_abs - 4.65)  # solar luminosities (r-band)
        lm = log_lum + ml_ratio
        lm = max(6.0, min(lm, 13.0))

        # Age from color (redder = older)
        la = 8.5 + 0.8 * mean_color  # Gyr in log
        la = max(7.0, min(la, 10.2))

        # Metallicity from color
        lz = -0.5 + 0.3 * mean_color
        lz = max(-2.0, min(lz, 0.5))

        # Dust extinction from color excess
        av_est = max(0, (mean_color - 0.5) * 0.8)
        av_est = min(av_est, 4.0)

        # SFR from mass and age
        mass_val = 10**lm
        age_val = 10**(la - 9)  # Gyr
        sfr_est = mass_val / (age_val * 1e9) if age_val > 0 else 0

        log_mass[i] = round(lm, 2)
        log_mass_err[i] = round(0.3, 2)
        log_age[i] = round(la, 2)
        log_age_err[i] = round(0.4, 2)
        log_z_met[i] = round(lz, 2)
        av[i] = round(av_est, 2)
        sfr[i] = round(np.log10(max(sfr_est, 1e-5)), 2)

    return log_mass, log_mass_err, log_age, log_age_err, log_z_met, av, sfr


def estimate_sed_mdn(mags, mag_errs, photo_z, band_names,
                     ckpt_emulator, ckpt_inverse):
    """SED fitting with trained MDN model."""
    try:
        from sed_fit.predict import predict_sed
        from sed_fit.config import SEDConfig

        cfg = SEDConfig()
        result = predict_sed(mags, photo_z, ckpt_emulator, ckpt_inverse,
                             len(band_names), cfg)
        return (result['log_mass'], result['log_mass_err'],
                result['log_age'], result['log_age_err'],
                result['log_z'], result['av'], result['sfr'])
    except Exception as e:
        print(f"WARNING: MDN model failed ({e}), falling back to empirical",
              file=sys.stderr)
        return estimate_sed_empirical(mags, mag_errs, photo_z, band_names)


def estimate_sed_backend(mags, mag_errs, photo_z, band_names, backend_name):
    """SED fitting using an SPS backend (direct fitting)."""
    from sed_fit.backends import get_backend

    backend = get_backend(backend_name)
    print(f"Using SPS backend: {backend.name()}", file=sys.stderr)

    result = backend.fit_sed(mags, mag_errs, photo_z, band_names)

    return (result.log_mass, result.log_mass_err,
            result.log_age, result.log_age_err,
            result.log_Z, result.Av, result.sfr)


def mode_list_backends():
    """Print available SPS backends and exit."""
    from sed_fit.backends import list_all, list_available

    all_backends = list_all()
    avail = list_available()

    print("#BACKENDS")
    for name in all_backends:
        status = "available" if name in avail else "not_installed"
        print(f"{name}\t{status}")


def main():
    parser = argparse.ArgumentParser(description='DS9 SED Fitting (AI)')
    parser.add_argument('fits', nargs='?', default=None,
                        help='Path to FITS file')
    parser.add_argument('--catalog', default=None, help='Path to catalog TSV')
    parser.add_argument('--bands', default='g,r,i,z',
                        help='Band names (comma-separated)')
    parser.add_argument('--mag-columns', default='',
                        help='Magnitude column names (comma-separated)')
    parser.add_argument('--photo-z-column', default='PHOTO_Z',
                        help='Photo-z column name')
    parser.add_argument('--checkpoint-emulator', default=None,
                        help='Path to SPS emulator checkpoint')
    parser.add_argument('--checkpoint-inverse', default=None,
                        help='Path to inverse model checkpoint')
    parser.add_argument('--backend', default='auto',
                        help='SPS backend: auto, analytic, fsps, bagpipes, '
                             'prospector, cigale, dense_basis '
                             '(default: auto)')
    parser.add_argument('--mode', default='fit',
                        choices=['fit', 'list-backends'],
                        help='Operation mode (default: fit)')
    parser.add_argument('--n-workers', type=int, default=1,
                        help='Number of workers (unused for SED)')
    args = parser.parse_args()

    # Handle list-backends mode
    if args.mode == 'list-backends':
        mode_list_backends()
        return

    # Fit mode requires FITS + catalog
    if args.fits is None or args.catalog is None:
        print("ERROR: FITS file and --catalog required for fit mode",
              file=sys.stderr)
        sys.exit(1)

    # Parse catalog
    print(f"Loading catalog: {args.catalog}", file=sys.stderr)
    header, rows = parse_catalog_tsv(args.catalog)
    if header is None or not rows:
        print("ERROR: Cannot parse catalog", file=sys.stderr)
        sys.exit(1)

    num_idx = get_column_index(header, 'NUMBER')
    if num_idx < 0:
        print("ERROR: NUMBER column not found", file=sys.stderr)
        sys.exit(1)

    # Determine magnitude columns
    bands = [b.strip() for b in args.bands.split(',')]
    if args.mag_columns:
        mag_col_names = [c.strip() for c in args.mag_columns.split(',')]
    else:
        mag_col_names = []
        for b in bands:
            for pattern in [f'MAG_{b}', f'MAG_AUTO_{b}', f'mag_{b}']:
                idx = get_column_index(header, pattern)
                if idx >= 0:
                    mag_col_names.append(pattern)
                    break
        if not mag_col_names:
            mag_auto = get_column_index(header, 'MAG_AUTO')
            if mag_auto >= 0:
                mag_col_names = ['MAG_AUTO']
                bands = ['auto']

    if not mag_col_names:
        print("ERROR: No magnitude columns found", file=sys.stderr)
        sys.exit(1)

    # Get photo-z column
    pz_idx = get_column_index(header, args.photo_z_column)

    # Get column indices
    mag_indices = []
    for mc in mag_col_names:
        mag_indices.append(get_column_index(header, mc))

    # Check for error columns
    err_indices = []
    for mc in mag_col_names:
        err_name = mc.replace('MAG_', 'MAGERR_')
        err_indices.append(get_column_index(header, err_name))

    # Build arrays
    n_sources = len(rows)
    n_bands = len(mag_col_names)
    mags = np.full((n_sources, n_bands), np.nan)
    mag_errs = np.full((n_sources, n_bands), np.nan)
    photo_z = np.full(n_sources, 0.1)  # default
    numbers = []

    for i, row in enumerate(rows):
        numbers.append(row[num_idx])
        for j, mi in enumerate(mag_indices):
            if mi >= 0 and mi < len(row):
                try:
                    val = float(row[mi])
                    if val < 90 and val > 0:
                        mags[i, j] = val
                except (ValueError, IndexError):
                    pass
            if err_indices[j] >= 0 and err_indices[j] < len(row):
                try:
                    mag_errs[i, j] = float(row[err_indices[j]])
                except (ValueError, IndexError):
                    pass
        if pz_idx >= 0 and pz_idx < len(row):
            try:
                pz = float(row[pz_idx])
                if pz > 0 and pz < 10:
                    photo_z[i] = pz
            except (ValueError, IndexError):
                pass

    print(f"  Sources: {n_sources}, Bands: {n_bands} ({','.join(mag_col_names)})",
          file=sys.stderr)
    if pz_idx >= 0:
        n_pz = np.sum(photo_z > 0.001)
        print(f"  Valid photo-z: {n_pz}", file=sys.stderr)
    else:
        print("  WARNING: No PHOTO_Z column, using default z=0.1",
              file=sys.stderr)

    # Determine fitting method based on backend
    backend = args.backend.lower()
    use_sps_backend = (backend not in ('auto', 'analytic'))

    if use_sps_backend:
        # Direct SPS backend fitting
        print(f"SPS backend: {backend}", file=sys.stderr)
        results = estimate_sed_backend(mags, mag_errs, photo_z, bands,
                                        backend)
    else:
        # AI pipeline (MDN) or empirical fallback
        use_mdn = (args.checkpoint_emulator and
                   os.path.exists(args.checkpoint_emulator) and
                   args.checkpoint_inverse and
                   os.path.exists(args.checkpoint_inverse))

        if use_mdn:
            print(f"Using MDN model", file=sys.stderr)
            results = estimate_sed_mdn(
                mags, mag_errs, photo_z, bands,
                args.checkpoint_emulator, args.checkpoint_inverse)
        else:
            print("Using empirical color-mass relations (no checkpoint)",
                  file=sys.stderr)
            results = estimate_sed_empirical(mags, mag_errs, photo_z, bands)

    log_mass, log_mass_err, log_age, log_age_err, log_z, av_arr, sfr_arr = results

    n_valid = np.sum(np.isfinite(log_mass))
    print(f"  Valid SED fits: {n_valid}/{n_sources}", file=sys.stderr)

    # Output TSV
    print(f"#SED_FIT\tN_FITTED={n_valid}\tN_SOURCES={n_sources}")
    print("NUMBER\tLOG_MASS\tLOG_MASS_ERR\tLOG_AGE\tLOG_AGE_ERR\t"
          "LOG_Z\tAV\tSFR")

    for i in range(n_sources):
        num = numbers[i]

        def fmt(v):
            return f"{v:.2f}" if np.isfinite(v) else "-99"

        print(f"{num}\t{fmt(log_mass[i])}\t{fmt(log_mass_err[i])}\t"
              f"{fmt(log_age[i])}\t{fmt(log_age_err[i])}\t"
              f"{fmt(log_z[i])}\t{fmt(av_arr[i])}\t{fmt(sfr_arr[i])}")


if __name__ == '__main__':
    main()
