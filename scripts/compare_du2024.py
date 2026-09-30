#!/usr/bin/env python3
"""
Compare OGFinder LSBG detection results with Du et al. (2024, RAA 24 055015)
reference catalog for DESI Legacy Survey DR9 brick 2828p530.

Generates publication-quality figures:
  - fig_sky_map.png       : RA-Dec sky distribution
  - fig_mu_reff.png       : mu_eff vs r_eff selection diagram
  - fig_sersic_dist.png   : Sersic n distribution histogram
  - fig_completeness.png  : Completeness & contamination vs mu_eff
  - fig_cutouts.png       : Representative matched/missed LSBG cutouts

Usage: cd OGFinder && python3 scripts/compare_du2024.py
"""

import os
import sys
import csv
import math
import numpy as np

RESULTS_DIR = "results/du2024"
DATA_DIR = "../data"

# ---- Du+2024 reference values ----

# Selection criteria (Du+2024 Table 1)
DU_MU_EFF_MIN = 24.2     # mag/arcsec^2
DU_MU_EFF_MAX = 28.8     # mag/arcsec^2
DU_R_EFF_MIN = 2.5       # arcsec
DU_ELLIPTICITY_MAX = 0.7
DU_TOTAL_LSBGS = 31825   # total catalog size (full survey)

# DESI DR9 pixel scale and ZP
PIXEL_SCALE = 0.262       # arcsec/pixel
ZP = 22.5                 # nanomaggy zeropoint

# Cross-match radius (tightened after R_EFF fix: sep.flux_radius
# centroids align well with Tractor positions)
MATCH_RADIUS_ARCSEC = 5.0

# Du+2024 reported statistics
DU_MEDIAN_SERSIC_N = 1.0  # approximate median
DU_MEDIAN_MU_EFF = 25.5   # approximate median mu_eff,g
DU_MEDIAN_R_EFF = 4.0     # approximate median r_eff (arcsec)


def load_detection_catalog(path):
    """Load OGFinder LSBG detection TSV catalog."""
    sources = []
    with open(path) as f:
        header = None
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if header is None:
                header = line.split('\t')
                continue
            parts = line.split('\t')
            if len(parts) < len(header):
                continue
            src = {}
            for k, v in zip(header, parts):
                try:
                    src[k] = float(v)
                except ValueError:
                    src[k] = v
            sources.append(src)
    return sources


def load_reference_catalog(path):
    """Load Du+2024 reference LSBG CSV catalog."""
    sources = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            src = {}
            for k, v in row.items():
                try:
                    src[k] = float(v)
                except ValueError:
                    src[k] = v
            sources.append(src)
    return sources


def pixel_to_wcs(x_pix, y_pix, fits_path):
    """Convert pixel coordinates to RA, Dec using FITS WCS."""
    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        with fits.open(fits_path) as hdul:
            w = WCS(hdul[0].header)
        coords = w.pixel_to_world(x_pix, y_pix)
        return coords.ra.deg, coords.dec.deg
    except Exception:
        # Fallback: return pixel coords as-is (for plotting only)
        return x_pix, y_pix


def angular_separation(ra1, dec1, ra2, dec2):
    """Angular separation in arcseconds between two positions (degrees)."""
    dra = (ra1 - ra2) * math.cos(math.radians(0.5 * (dec1 + dec2)))
    ddec = dec1 - dec2
    return math.sqrt(dra**2 + ddec**2) * 3600.0


def cross_match(det_sources, ref_sources, radius_arcsec=MATCH_RADIUS_ARCSEC):
    """Cross-match detected sources with reference catalog.

    Returns (matched_det, matched_ref, unmatched_det, unmatched_ref).
    """
    matched_det = []
    matched_ref = []
    matched_ref_idx = set()
    unmatched_det = []

    for det in det_sources:
        det_ra = det.get('ra', det.get('RA', 0))
        det_dec = det.get('dec', det.get('DEC', 0))

        best_dist = float('inf')
        best_idx = -1
        for j, ref in enumerate(ref_sources):
            if j in matched_ref_idx:
                continue
            ref_ra = ref.get('ra', 0)
            ref_dec = ref.get('dec', 0)
            dist = angular_separation(det_ra, det_dec, ref_ra, ref_dec)
            if dist < best_dist:
                best_dist = dist
                best_idx = j

        if best_dist <= radius_arcsec and best_idx >= 0:
            matched_det.append(det)
            matched_ref.append(ref_sources[best_idx])
            matched_ref_idx.add(best_idx)
        else:
            unmatched_det.append(det)

    unmatched_ref = [ref for j, ref in enumerate(ref_sources)
                     if j not in matched_ref_idx]

    return matched_det, matched_ref, unmatched_det, unmatched_ref


def add_wcs_to_detections(det_sources, fits_path):
    """Add RA/Dec to detected sources using FITS WCS."""
    if not os.path.exists(fits_path):
        return

    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        import warnings
        warnings.filterwarnings('ignore', category=fits.verify.VerifyWarning)
        warnings.filterwarnings('ignore', message='.*FITSFixedWarning.*')
        with fits.open(fits_path) as hdul:
            # Try each HDU to find one with a valid 2D WCS
            w = None
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim >= 2:
                    try:
                        wcs_candidate = WCS(hdu.header)
                        if wcs_candidate.naxis >= 2:
                            w = wcs_candidate
                            break
                    except Exception:
                        continue
            if w is None:
                w = WCS(hdul[0].header)
    except Exception:
        return

    for src in det_sources:
        x = src.get('X_IMAGE', 0) - 1  # 1-indexed to 0-indexed
        y = src.get('Y_IMAGE', 0) - 1
        try:
            coord = w.pixel_to_world(x, y)
            src['ra'] = coord.ra.deg
            src['dec'] = coord.dec.deg
        except Exception:
            pass


def plot_sky_map(det_sources, ref_sources, matched_det, unmatched_ref):
    """Figure 1: RA-Dec sky distribution."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    # Reference LSBGs
    ref_ra = [s.get('ra', 0) for s in ref_sources]
    ref_dec = [s.get('dec', 0) for s in ref_sources]
    ax.scatter(ref_ra, ref_dec, s=50, c='red', alpha=0.5, marker='s',
               edgecolors='darkred', linewidth=0.5,
               label=f'Du+2024 reference ({len(ref_sources)})')

    # Our detections
    det_ra = [s.get('ra', 0) for s in det_sources]
    det_dec = [s.get('dec', 0) for s in det_sources]
    ax.scatter(det_ra, det_dec, s=40, c='dodgerblue', alpha=0.6, marker='o',
               edgecolors='navy', linewidth=0.5,
               label=f'OGFinder ({len(det_sources)})')

    # Matched sources
    match_ra = [s.get('ra', 0) for s in matched_det]
    match_dec = [s.get('dec', 0) for s in matched_det]
    if match_ra:
        ax.scatter(match_ra, match_dec, s=80, facecolors='none',
                   edgecolors='green', linewidth=1.5, marker='o',
                   label=f'Matched ({len(matched_det)})')

    ax.set_xlabel('RA (deg)', fontsize=13)
    ax.set_ylabel('Dec (deg)', fontsize=13)
    ax.set_title('Du+2024 LSBG Sky Distribution\n'
                 'DESI DR9 Brick 2828p530', fontsize=12)
    ax.legend(fontsize=9, loc='upper left')
    ax.invert_xaxis()  # RA increases to the left
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    out = os.path.join(RESULTS_DIR, 'fig_sky_map.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_mu_reff(det_sources, ref_sources):
    """Figure 2: mu_eff vs r_eff diagram with Du+2024 selection box."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    # Reference catalog
    ref_mu = [s.get('mu_eff_g', 0) for s in ref_sources]
    ref_r = [s.get('shape_r', 0) for s in ref_sources]
    ax.scatter(ref_r, ref_mu, s=30, c='red', alpha=0.4, marker='s',
               edgecolors='darkred', linewidth=0.3,
               label=f'Du+2024 reference ({len(ref_sources)})')

    # Our detections
    det_mu = [s.get('MU_EFF', s.get('mu_eff', 0)) for s in det_sources]
    det_r = [s.get('R_EFF_ARCSEC', s.get('r_eff_arcsec', 0))
             for s in det_sources]
    ax.scatter(det_r, det_mu, s=25, c='dodgerblue', alpha=0.5, marker='o',
               edgecolors='navy', linewidth=0.3,
               label=f'OGFinder ({len(det_sources)})')

    # Du+2024 selection box
    box_r = [DU_R_EFF_MIN, 20, 20, DU_R_EFF_MIN, DU_R_EFF_MIN]
    box_mu = [DU_MU_EFF_MIN, DU_MU_EFF_MIN, DU_MU_EFF_MAX,
              DU_MU_EFF_MAX, DU_MU_EFF_MIN]
    ax.plot(box_r, box_mu, 'k--', linewidth=1.5, alpha=0.7,
            label=f'Du+2024 selection ({DU_MU_EFF_MIN}<$\\mu_{{eff}}$<{DU_MU_EFF_MAX})')

    # Reference lines
    ax.axhline(DU_MU_EFF_MIN, color='gray', linestyle=':', linewidth=0.8,
               alpha=0.5)
    ax.axhline(DU_MU_EFF_MAX, color='gray', linestyle=':', linewidth=0.8,
               alpha=0.5)
    ax.axvline(DU_R_EFF_MIN, color='gray', linestyle=':', linewidth=0.8,
               alpha=0.5)

    # Constant SB lines (for context)
    r_arr = np.logspace(np.log10(1), np.log10(30), 100)
    for mu_line in [24, 26, 28]:
        ax.plot(r_arr, np.full_like(r_arr, mu_line), ':',
                color='lightgray', linewidth=0.5)

    ax.set_xlabel(r'$r_{eff}$ (arcsec)', fontsize=13)
    ax.set_ylabel(r'$\mu_{eff,g}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_title(r'$\mu_{eff}$ vs $r_{eff}$ Diagram'
                 '\nDu+2024 Selection Region', fontsize=12)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlim(1, 60)
    ax.set_ylim(30, 22)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_mu_reff.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_sersic_distribution(det_sources, ref_sources=None):
    """Figure 3: Sersic n distribution histogram."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    # Our detections
    det_n = [s.get('SERSIC_N', float('nan')) for s in det_sources]
    det_n = [n for n in det_n if np.isfinite(n) and 0 < n < 10]

    bins = np.arange(0, 6.5, 0.5)
    if det_n:
        ax.hist(det_n, bins=bins, color='dodgerblue', alpha=0.7,
                edgecolor='navy', linewidth=0.5,
                label=f'OGFinder (N={len(det_n)}, '
                      f'median={np.median(det_n):.2f})')

    # Du+2024 reference: median Sersic n ~ 1
    ax.axvline(DU_MEDIAN_SERSIC_N, color='red', linestyle='--', linewidth=2,
               label=f'Du+2024 median n={DU_MEDIAN_SERSIC_N}')

    # Typical UDG range
    ax.axvspan(0.5, 1.5, alpha=0.08, color='green',
               label='Typical UDG range (n=0.5-1.5)')

    ax.set_xlabel(r'S\'ersic index $n$', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title(r'S\'ersic Index Distribution of LSBG Candidates', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xlim(0, 6)

    out = os.path.join(RESULTS_DIR, 'fig_sersic_dist.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    if det_n:
        print(f"  Sersic n: median={np.median(det_n):.2f}, "
              f"mean={np.mean(det_n):.2f}, std={np.std(det_n):.2f}")
        print(f"  Range: {np.min(det_n):.2f} - {np.max(det_n):.2f}")


def plot_completeness(det_sources, ref_sources, matched_det, unmatched_det,
                      matched_ref=None):
    """Figure 4: Completeness and contamination vs mu_eff."""
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Bin by mu_eff
    mu_bins = np.arange(23.5, 29.5, 0.5)
    mu_centers = 0.5 * (mu_bins[:-1] + mu_bins[1:])

    # Completeness: fraction of reference sources recovered
    ref_mu = np.array([s.get('mu_eff_g', 0) for s in ref_sources])
    if matched_ref is None:
        matched_ref = []
    matched_ref_mu = np.array([s.get('mu_eff_g', 0) for s in matched_ref])

    completeness = []
    for lo, hi in zip(mu_bins[:-1], mu_bins[1:]):
        n_ref = np.sum((ref_mu >= lo) & (ref_mu < hi))
        n_match = np.sum((matched_ref_mu >= lo) & (matched_ref_mu < hi)) \
            if len(matched_ref_mu) > 0 else 0
        if n_ref > 0:
            completeness.append(n_match / n_ref)
        else:
            completeness.append(np.nan)
    completeness = np.array(completeness)

    valid = np.isfinite(completeness)
    ax1.plot(mu_centers[valid], completeness[valid] * 100, 'o-',
             color='navy', ms=6, linewidth=1.5)
    ax1.axhline(50, color='gray', linestyle=':', linewidth=1,
                label='50% completeness')
    ax1.set_xlabel(r'$\mu_{eff}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax1.set_ylabel('Completeness (%)', fontsize=13)
    ax1.set_title('Recovery Rate vs Surface Brightness', fontsize=12)
    ax1.set_ylim(-5, 105)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Contamination: fraction of our detections NOT in reference
    det_mu = np.array([s.get('MU_EFF', s.get('mu_eff', 0))
                       for s in det_sources])
    unmatch_mu = np.array([s.get('MU_EFF', s.get('mu_eff', 0))
                           for s in unmatched_det])

    contamination = []
    for lo, hi in zip(mu_bins[:-1], mu_bins[1:]):
        n_det = np.sum((det_mu >= lo) & (det_mu < hi))
        n_unmatch = np.sum((unmatch_mu >= lo) & (unmatch_mu < hi))
        if n_det > 0:
            contamination.append(n_unmatch / n_det)
        else:
            contamination.append(np.nan)
    contamination = np.array(contamination)

    valid2 = np.isfinite(contamination)
    ax2.plot(mu_centers[valid2], contamination[valid2] * 100, 's-',
             color='crimson', ms=6, linewidth=1.5)
    ax2.set_xlabel(r'$\mu_{eff}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax2.set_ylabel('Contamination (%)', fontsize=13)
    ax2.set_title('False Positive Rate vs Surface Brightness', fontsize=12)
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.3)

    # Annotate
    ax2.text(0.97, 0.97,
             'Contamination expected:\n'
             'Du+2024 uses XGBoost ML\n'
             'OGFinder uses photometric cuts',
             transform=ax2.transAxes, fontsize=8, va='top', ha='right',
             bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                       edgecolor='gray', alpha=0.7))

    out = os.path.join(RESULTS_DIR, 'fig_completeness.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    # Print summary statistics
    overall_completeness = len(matched_det) / max(len(ref_sources), 1) * 100
    overall_contamination = len(unmatched_det) / max(len(det_sources), 1) * 100
    print(f"  Overall completeness: {overall_completeness:.1f}% "
          f"({len(matched_det)}/{len(ref_sources)})")
    print(f"  Overall contamination: {overall_contamination:.1f}% "
          f"({len(unmatched_det)}/{len(det_sources)})")


def plot_cutouts(det_sources, ref_sources, matched_det, unmatched_ref,
                 fits_path):
    """Figure 5: Representative matched/missed LSBG cutouts (4x4 grid)."""
    import matplotlib.pyplot as plt

    if not os.path.exists(fits_path):
        print("  Cutout figure skipped (FITS file not found)")
        return

    try:
        from astropy.io import fits
        with fits.open(fits_path) as hdul:
            if hdul[0].data is not None:
                image_data = hdul[0].data.astype(np.float64)
            else:
                for ext in hdul[1:]:
                    if ext.data is not None and ext.data.ndim >= 2:
                        image_data = ext.data.astype(np.float64)
                        break
                else:
                    print("  No image data for cutouts")
                    return
        while image_data.ndim > 2:
            image_data = image_data[0]
    except Exception as e:
        print(f"  Cutout figure skipped: {e}")
        return

    ny, nx = image_data.shape
    cutout_size = 80  # pixels (~21 arcsec at 0.262"/pix)

    fig, axes = plt.subplots(4, 4, figsize=(14, 14))
    fig.suptitle('LSBG Candidates: Matched (top 2 rows) & '
                 'Missed by Reference (bottom 2 rows)', fontsize=13)

    # Top 2 rows: matched detections
    show_matched = matched_det[:8]
    for i, src in enumerate(show_matched):
        row, col = divmod(i, 4)
        ax = axes[row][col]
        cx = int(src.get('X_IMAGE', nx // 2)) - 1
        cy = int(src.get('Y_IMAGE', ny // 2)) - 1

        x0 = max(0, cx - cutout_size)
        x1 = min(nx, cx + cutout_size)
        y0 = max(0, cy - cutout_size)
        y1 = min(ny, cy + cutout_size)

        stamp = image_data[y0:y1, x0:x1]
        vmin = np.nanpercentile(stamp, 5)
        vmax = np.nanpercentile(stamp, 99)

        ax.imshow(stamp, origin='lower', cmap='gray_r',
                  vmin=vmin, vmax=vmax)
        mu = src.get('MU_EFF', 0)
        r = src.get('R_EFF_ARCSEC', 0)
        grade = src.get('LSBG_GRADE', '?')
        ax.set_title(f'#{int(src.get("NUMBER", i+1))} '
                     f'$\\mu$={mu:.1f} $r_e$={r:.1f}" ({grade})',
                     fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        # Green border = matched
        for spine in ax.spines.values():
            spine.set_color('green')
            spine.set_linewidth(2)

    # Bottom 2 rows: unmatched detections (excess / missed by reference)
    unmatched_show = [s for s in det_sources
                      if s not in matched_det][:8]
    for i, src in enumerate(unmatched_show):
        row = 2 + i // 4
        col = i % 4
        if row >= 4:
            break
        ax = axes[row][col]
        cx = int(src.get('X_IMAGE', nx // 2)) - 1
        cy = int(src.get('Y_IMAGE', ny // 2)) - 1

        x0 = max(0, cx - cutout_size)
        x1 = min(nx, cx + cutout_size)
        y0 = max(0, cy - cutout_size)
        y1 = min(ny, cy + cutout_size)

        stamp = image_data[y0:y1, x0:x1]
        vmin = np.nanpercentile(stamp, 5)
        vmax = np.nanpercentile(stamp, 99)

        ax.imshow(stamp, origin='lower', cmap='gray_r',
                  vmin=vmin, vmax=vmax)
        mu = src.get('MU_EFF', 0)
        r = src.get('R_EFF_ARCSEC', 0)
        grade = src.get('LSBG_GRADE', '?')
        ax.set_title(f'#{int(src.get("NUMBER", i+1))} '
                     f'$\\mu$={mu:.1f} $r_e$={r:.1f}" ({grade})',
                     fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        # Orange border = unmatched
        for spine in ax.spines.values():
            spine.set_color('orange')
            spine.set_linewidth(2)

    # Clear empty axes
    total_shown = len(show_matched) + len(unmatched_show)
    for idx in range(total_shown, 16):
        row, col = divmod(idx, 4)
        axes[row][col].axis('off')

    out = os.path.join(RESULTS_DIR, 'fig_cutouts.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def print_parameter_comparison():
    """Print parameter comparison table: JWST ICL vs DESI LSBG."""
    print("\n--- Parameter Comparison: JWST (ICL) vs DESI (LSBG) ---")
    print(f"{'Parameter':<25} {'JWST ICL':>12} {'DESI LSBG':>12} {'Reason':>30}")
    print("-" * 82)
    rows = [
        ('pixel-scale', '0.063', '0.262', 'NIRCam vs ground-based'),
        ('mag-zeropoint', 'auto (26.5)', '22.5', 'MJy/sr vs nanomaggies'),
        ('detect-thresh', '5.0 sigma', '1.0 sigma', 'ICL preserve vs LSB detect'),
        ('detect-minarea', '5 pix', '100 pix', 'HST res vs seeing-limited'),
        ('filter-kernel', '-', 'gauss7x7', 'Larger PSF'),
        ('mu-eff-min', '-', '24.0', 'LSBG definition cutoff'),
        ('r-eff-min', '-', '2.5"', 'Du+2024 criterion'),
        ('multiscale', 'No', 'Yes (1,2,4)', 'Extended diffuse sources'),
        ('bkg-mesh-size', '-', '256', 'Large-scale gradients'),
    ]
    for param, jwst, desi, reason in rows:
        print(f"  {param:<23} {jwst:>12} {desi:>12} {reason:>30}")


def main():
    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: Results directory not found: {RESULTS_DIR}")
        print("Run scripts/reproduce_du2024.sh first.")
        sys.exit(1)

    print("=" * 60)
    print("Du et al. (2024, RAA 24 055015) LSBG Comparison")
    print("  Survey: DESI Legacy Survey DR9 (BASS+MzLS)")
    print("  Brick: 2828p530 (RA=282.8, Dec=+53.0)")
    print("  Pixel scale: 0.262\"/pixel, ZP=22.5")
    print("=" * 60)

    # Load detection catalog
    det_file = os.path.join(RESULTS_DIR, 'lsbg_detected.tsv')
    if not os.path.exists(det_file):
        print(f"ERROR: Detection catalog not found: {det_file}")
        sys.exit(1)
    det_sources = load_detection_catalog(det_file)
    print(f"\nLoaded detections: {len(det_sources)} LSBG candidates")

    # Add WCS coordinates to detections
    fits_path = os.path.join(DATA_DIR, 'desi_dr9_2828p530_r.fits')
    add_wcs_to_detections(det_sources, fits_path)

    # Load reference catalog
    ref_file = os.path.join(DATA_DIR, 'du2024_reference_lsbg.csv')
    ref_sources = []
    if os.path.exists(ref_file):
        ref_sources = load_reference_catalog(ref_file)
        print(f"Loaded reference: {len(ref_sources)} Du+2024 LSBGs")
    else:
        print(f"WARNING: Reference catalog not found: {ref_file}")
        print("  Run: python3 scripts/fetch_du2024_data.py")
        print("  Generating figures with detections only...")

    # Cross-match
    matched_det, matched_ref, unmatched_det, unmatched_ref = \
        [], [], det_sources, ref_sources
    if ref_sources and any('ra' in s for s in det_sources):
        print(f"\nCross-matching (radius={MATCH_RADIUS_ARCSEC}\")")
        matched_det, matched_ref, unmatched_det, unmatched_ref = \
            cross_match(det_sources, ref_sources)
        print(f"  Matched: {len(matched_det)}")
        print(f"  Unmatched detections: {len(unmatched_det)}")
        print(f"  Missed references: {len(unmatched_ref)}")

    # Detection statistics
    print("\n--- Detection Statistics ---")
    det_mu = [s.get('MU_EFF', float('nan')) for s in det_sources]
    det_r = [s.get('R_EFF_ARCSEC', float('nan')) for s in det_sources]
    det_mu = [x for x in det_mu if np.isfinite(x)]
    det_r = [x for x in det_r if np.isfinite(x)]

    if det_mu:
        print(f"  mu_eff range: {min(det_mu):.1f} - {max(det_mu):.1f} "
              f"(median {np.median(det_mu):.1f})")
    if det_r:
        print(f"  r_eff range: {min(det_r):.1f} - {max(det_r):.1f} arcsec "
              f"(median {np.median(det_r):.1f})")

    grades = {}
    for s in det_sources:
        g = s.get('LSBG_GRADE', '?')
        grades[g] = grades.get(g, 0) + 1
    if grades:
        print(f"  Grades: {grades}")

    # Generate figures
    print("\n--- Generating Figures ---")

    # Figure 1: Sky map
    if ref_sources and any('ra' in s for s in det_sources):
        print("\nFig 1: Sky map")
        plot_sky_map(det_sources, ref_sources, matched_det, unmatched_ref)

    # Figure 2: mu_eff vs r_eff
    print("\nFig 2: mu_eff vs r_eff")
    plot_mu_reff(det_sources, ref_sources)

    # Figure 3: Sersic n distribution
    print("\nFig 3: Sersic n distribution")
    plot_sersic_distribution(det_sources)

    # Figure 4: Completeness
    if ref_sources and matched_det:
        print("\nFig 4: Completeness & contamination")
        plot_completeness(det_sources, ref_sources,
                          matched_det, unmatched_det, matched_ref)

    # Figure 5: Cutouts
    print("\nFig 5: LSBG cutouts")
    plot_cutouts(det_sources, ref_sources, matched_det, unmatched_ref,
                 fits_path)

    # Parameter comparison table
    print_parameter_comparison()

    print("\n" + "=" * 60)
    print("Done. Figures saved to", RESULTS_DIR)


if __name__ == '__main__':
    main()
