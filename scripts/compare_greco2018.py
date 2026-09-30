#!/usr/bin/env python3
"""
Compare OGFinder LSBG detection results with Greco et al. (2018, ApJ 857 104)
reference catalog on an HSC-SSP test patch.

Generates figures:
  - fig_sky_map.png       : RA-Dec sky distribution (detections vs reference)
  - fig_mu_reff.png       : mu_eff vs r_eff selection diagram
  - fig_sersic_dist.png   : Sersic n distribution histogram
  - fig_completeness.png  : Completeness & contamination vs mu_eff
  - fig_cutouts.png       : Representative matched/missed LSBG cutouts

Usage: cd OGFinder && python3 scripts/compare_greco2018.py
"""

import os
import sys
import csv
import math
import numpy as np

RESULTS_DIR = "results/greco2018"
DATA_DIR = "../data"

# ---- Greco+2018 reference values ----

# Selection criteria (Greco+2018 Section 3)
GRECO_MU_EFF_MIN = 24.3     # mag/arcsec^2 (g-band)
GRECO_MU_EFF_MAX = 28.8
GRECO_R_EFF_MIN = 2.5       # arcsec
GRECO_R_EFF_MAX = 14.0      # arcsec
GRECO_ELLIPTICITY_MAX = 0.7
GRECO_TOTAL_LSBGS = 781     # full catalog

# HSC-SSP pixel scale and ZP
PIXEL_SCALE = 0.168          # arcsec/pixel
ZP = 27.0                   # HSC-SSP coadd zeropoint

# Cross-match radius
MATCH_RADIUS_ARCSEC = 5.0

# Greco+2018 reported statistics
GRECO_MEDIAN_SERSIC_N = 1.0
GRECO_MEDIAN_R_EFF = 5.0     # arcsec (approximate)


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
    """Load Greco+2018 reference LSBG CSV catalog (patch subset)."""
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
            # Normalize VizieR column names
            if 'RAJ2000' in src and 'ra' not in src:
                src['ra'] = src['RAJ2000']
            if 'DEJ2000' in src and 'dec' not in src:
                src['dec'] = src['DEJ2000']
            if 'RAdeg' in src and 'ra' not in src:
                src['ra'] = src['RAdeg']
            if 'DEdeg' in src and 'dec' not in src:
                src['dec'] = src['DEdeg']
            if 'g-r' in src and 'g_r' not in src:
                src['g_r'] = src['g-r']
            if 'g-i' in src and 'g_i' not in src:
                src['g_i'] = src['g-i']
            sources.append(src)
    return sources


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
        det_ra = det.get('RA', det.get('ra', 0))
        det_dec = det.get('DEC', det.get('dec', 0))

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
    """Add RA/Dec to detected sources using FITS WCS (if not already present)."""
    # Check if RA/DEC already present
    if det_sources and ('RA' in det_sources[0] or 'ra' in det_sources[0]):
        ra_key = 'RA' if 'RA' in det_sources[0] else 'ra'
        if isinstance(det_sources[0][ra_key], (int, float)) and det_sources[0][ra_key] != 0:
            return

    if not os.path.exists(fits_path):
        return

    try:
        from astropy.io import fits
        from astropy.wcs import WCS
        import warnings
        warnings.filterwarnings('ignore')
        with fits.open(fits_path) as hdul:
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
            src['RA'] = coord.ra.deg
            src['DEC'] = coord.dec.deg
        except Exception:
            pass


def compute_ref_mu_eff(src):
    """Compute approximate mu_eff from Greco+2018 catalog.

    The catalog stores iSB (i-band central surface brightness = mu_0).
    mu_eff = mu_0 + 2.5*b_n/ln(10) where b_n = sersic_bn(n).
    """
    n = src.get('n', 1.0)
    mu0 = src.get('iSB', 25.0)
    # b_n approximation (Ciotti & Bertin 1999)
    if n > 0.36:
        bn = 2.0 * n - 1.0/3.0 + 4.0/(405.0*n)
    else:
        bn = 0.01
    mu_eff = mu0 + 2.5 * bn / math.log(10)
    return mu_eff


def plot_sky_map(det_sources, ref_sources, matched_det, unmatched_ref):
    """Figure 1: RA-Dec sky distribution."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    # Reference LSBGs
    ref_ra = [s.get('ra', 0) for s in ref_sources]
    ref_dec = [s.get('dec', 0) for s in ref_sources]
    ax.scatter(ref_ra, ref_dec, s=60, c='red', alpha=0.5, marker='s',
               edgecolors='darkred', linewidth=0.5,
               label=f'Greco+2018 reference ({len(ref_sources)})')

    # Our detections
    det_ra = [s.get('RA', s.get('ra', 0)) for s in det_sources]
    det_dec = [s.get('DEC', s.get('dec', 0)) for s in det_sources]
    ax.scatter(det_ra, det_dec, s=40, c='dodgerblue', alpha=0.6, marker='o',
               edgecolors='navy', linewidth=0.5,
               label=f'OGFinder ({len(det_sources)})')

    # Matched sources
    match_ra = [s.get('RA', s.get('ra', 0)) for s in matched_det]
    match_dec = [s.get('DEC', s.get('dec', 0)) for s in matched_det]
    if match_ra:
        ax.scatter(match_ra, match_dec, s=100, facecolors='none',
                   edgecolors='green', linewidth=2, marker='o',
                   label=f'Matched ({len(matched_det)})')

    ax.set_xlabel('RA (deg)', fontsize=13)
    ax.set_ylabel('Dec (deg)', fontsize=13)
    ax.set_title('Greco+2018 LSBG Sky Distribution\n'
                 'HSC-SSP Test Patch', fontsize=12)
    ax.legend(fontsize=9, loc='upper left')
    ax.invert_xaxis()
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal')

    out = os.path.join(RESULTS_DIR, 'fig_sky_map.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_mu_reff(det_sources, ref_sources):
    """Figure 2: mu_eff vs r_eff diagram with Greco+2018 selection box."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    # Reference catalog (compute mu_eff from mu_0 + Sersic correction)
    ref_mu = [compute_ref_mu_eff(s) for s in ref_sources]
    ref_r = [s.get('Reff', 0) for s in ref_sources]
    ax.scatter(ref_r, ref_mu, s=30, c='red', alpha=0.4, marker='s',
               edgecolors='darkred', linewidth=0.3,
               label=f'Greco+2018 reference ({len(ref_sources)})')

    # Our detections
    det_mu = [s.get('MU_EFF', s.get('mu_eff', 0)) for s in det_sources]
    det_r = [s.get('R_EFF_ARCSEC', s.get('r_eff_arcsec', 0))
             for s in det_sources]
    ax.scatter(det_r, det_mu, s=25, c='dodgerblue', alpha=0.5, marker='o',
               edgecolors='navy', linewidth=0.3,
               label=f'OGFinder ({len(det_sources)})')

    # Greco+2018 selection box
    box_r = [GRECO_R_EFF_MIN, GRECO_R_EFF_MAX, GRECO_R_EFF_MAX,
             GRECO_R_EFF_MIN, GRECO_R_EFF_MIN]
    box_mu = [GRECO_MU_EFF_MIN, GRECO_MU_EFF_MIN, GRECO_MU_EFF_MAX,
              GRECO_MU_EFF_MAX, GRECO_MU_EFF_MIN]
    ax.plot(box_r, box_mu, 'k--', linewidth=1.5, alpha=0.7,
            label=f'Greco+2018 selection')

    ax.set_xlabel(r'$r_{eff}$ (arcsec)', fontsize=13)
    ax.set_ylabel(r'$\mu_{eff}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_title(r'$\mu_{eff}$ vs $r_{eff}$ Diagram'
                 '\nGreco+2018 Selection Region', fontsize=12)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlim(1, 30)
    ax.set_ylim(30, 22)
    ax.legend(fontsize=9, loc='upper right')
    ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_mu_reff.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_sersic_distribution(det_sources, ref_sources):
    """Figure 3: Sersic n distribution histogram."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(8, 5))

    bins = np.arange(0, 6.5, 0.5)

    # Reference
    ref_n = [s.get('n', float('nan')) for s in ref_sources]
    ref_n = [n for n in ref_n if np.isfinite(n) and 0 < n < 10]
    if ref_n:
        ax.hist(ref_n, bins=bins, color='red', alpha=0.3,
                edgecolor='darkred', linewidth=0.5,
                label=f'Greco+2018 (N={len(ref_n)}, '
                      f'median={np.median(ref_n):.2f})')

    # Our detections
    det_n = [s.get('SERSIC_N', float('nan')) for s in det_sources]
    det_n = [n for n in det_n if np.isfinite(n) and 0 < n < 10]
    if det_n:
        ax.hist(det_n, bins=bins, color='dodgerblue', alpha=0.5,
                edgecolor='navy', linewidth=0.5,
                label=f'OGFinder (N={len(det_n)}, '
                      f'median={np.median(det_n):.2f})')

    # Typical UDG range
    ax.axvspan(0.5, 1.5, alpha=0.08, color='green',
               label='Typical UDG range (n=0.5-1.5)')

    ax.set_xlabel(r'S\'ersic index $n$', fontsize=13)
    ax.set_ylabel('Count', fontsize=13)
    ax.set_title(r'S\'ersic Index Distribution', fontsize=12)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xlim(0, 6)

    out = os.path.join(RESULTS_DIR, 'fig_sersic_dist.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_completeness(det_sources, ref_sources, matched_det, unmatched_det,
                      matched_ref):
    """Figure 4: Completeness and contamination vs mu_eff."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    mu_bins = np.arange(23.5, 29.5, 0.5)
    mu_centers = 0.5 * (mu_bins[:-1] + mu_bins[1:])

    # Completeness
    ref_mu = np.array([compute_ref_mu_eff(s) for s in ref_sources])
    matched_ref_mu = np.array([compute_ref_mu_eff(s) for s in matched_ref])

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
    ax1.axhline(80, color='green', linestyle=':', linewidth=1, alpha=0.5,
                label='80% (Greco+2018 target)')
    ax1.set_xlabel(r'$\mu_{eff}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax1.set_ylabel('Completeness (%)', fontsize=13)
    ax1.set_title('Recovery Rate vs Surface Brightness', fontsize=12)
    ax1.set_ylim(-5, 105)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)

    # Contamination
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

    out = os.path.join(RESULTS_DIR, 'fig_completeness.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    # Print summary
    overall_completeness = len(matched_det) / max(len(ref_sources), 1) * 100
    overall_contamination = len(unmatched_det) / max(len(det_sources), 1) * 100
    print(f"  Overall completeness: {overall_completeness:.1f}% "
          f"({len(matched_det)}/{len(ref_sources)})")
    print(f"  Overall contamination: {overall_contamination:.1f}% "
          f"({len(unmatched_det)}/{len(det_sources)})")


def plot_cutouts(det_sources, ref_sources, matched_det, unmatched_ref,
                 fits_path):
    """Figure 5: Representative LSBG cutouts."""
    import matplotlib
    matplotlib.use('Agg')
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
    cutout_size = 120  # pixels (~20" at 0.168"/pix)

    n_rows_fig = 4
    n_cols_fig = 4
    fig, axes = plt.subplots(n_rows_fig, n_cols_fig, figsize=(14, 14))
    fig.suptitle('LSBG Candidates: Matched (top 2 rows) & '
                 'Unmatched (bottom 2 rows)', fontsize=13)

    # Top 2 rows: matched detections
    show_matched = matched_det[:8]
    for i, src in enumerate(show_matched):
        row, col = divmod(i, n_cols_fig)
        ax = axes[row][col]
        cx = int(src.get('X_IMAGE', nx // 2)) - 1
        cy = int(src.get('Y_IMAGE', ny // 2)) - 1

        x0, x1 = max(0, cx - cutout_size), min(nx, cx + cutout_size)
        y0, y1 = max(0, cy - cutout_size), min(ny, cy + cutout_size)

        stamp = image_data[y0:y1, x0:x1]
        vmin = np.nanpercentile(stamp, 5)
        vmax = np.nanpercentile(stamp, 99)

        ax.imshow(stamp, origin='lower', cmap='gray_r',
                  vmin=vmin, vmax=vmax)
        mu = src.get('MU_EFF', 0)
        r = src.get('R_EFF_ARCSEC', 0)
        grade = src.get('LSBG_GRADE', '?')
        ax.set_title(f'$\\mu$={mu:.1f} $r_e$={r:.1f}" ({grade})', fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color('green')
            spine.set_linewidth(2)

    # Bottom 2 rows: unmatched detections
    unmatched_show = [s for s in det_sources if s not in matched_det][:8]
    for i, src in enumerate(unmatched_show):
        row = 2 + i // n_cols_fig
        col = i % n_cols_fig
        if row >= n_rows_fig:
            break
        ax = axes[row][col]
        cx = int(src.get('X_IMAGE', nx // 2)) - 1
        cy = int(src.get('Y_IMAGE', ny // 2)) - 1

        x0, x1 = max(0, cx - cutout_size), min(nx, cx + cutout_size)
        y0, y1 = max(0, cy - cutout_size), min(ny, cy + cutout_size)

        stamp = image_data[y0:y1, x0:x1]
        vmin = np.nanpercentile(stamp, 5)
        vmax = np.nanpercentile(stamp, 99)

        ax.imshow(stamp, origin='lower', cmap='gray_r',
                  vmin=vmin, vmax=vmax)
        mu = src.get('MU_EFF', 0)
        r = src.get('R_EFF_ARCSEC', 0)
        grade = src.get('LSBG_GRADE', '?')
        ax.set_title(f'$\\mu$={mu:.1f} $r_e$={r:.1f}" ({grade})', fontsize=8)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_color('orange')
            spine.set_linewidth(2)

    # Clear empty axes
    total_shown = len(show_matched) + len(unmatched_show)
    for idx in range(total_shown, n_rows_fig * n_cols_fig):
        row, col = divmod(idx, n_cols_fig)
        axes[row][col].axis('off')

    out = os.path.join(RESULTS_DIR, 'fig_cutouts.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def plot_parameter_comparison(matched_det, matched_ref):
    """Figure 6: Parameter comparison for matched sources."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not matched_det or not matched_ref:
        print("  Parameter comparison skipped (no matches)")
        return

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # r_eff comparison
    det_r = [s.get('R_EFF_ARCSEC', float('nan')) for s in matched_det]
    ref_r = [s.get('Reff', float('nan')) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_r, ref_r)]
    det_r_v = [d for d, v in zip(det_r, valid) if v]
    ref_r_v = [r for r, v in zip(ref_r, valid) if v]

    if det_r_v:
        ax = axes[0]
        ax.scatter(ref_r_v, det_r_v, s=30, c='dodgerblue', alpha=0.6)
        lim = [0, max(max(ref_r_v), max(det_r_v)) * 1.2]
        ax.plot(lim, lim, 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel('Greco+2018 $r_{eff}$ (arcsec)', fontsize=11)
        ax.set_ylabel('OGFinder $r_{eff}$ (arcsec)', fontsize=11)
        ax.set_title('Effective Radius', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Sersic n comparison
    det_n = [s.get('SERSIC_N', float('nan')) for s in matched_det]
    ref_n = [s.get('n', float('nan')) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_n, ref_n)]
    det_n_v = [d for d, v in zip(det_n, valid) if v]
    ref_n_v = [r for r, v in zip(ref_n, valid) if v]

    if det_n_v:
        ax = axes[1]
        ax.scatter(ref_n_v, det_n_v, s=30, c='dodgerblue', alpha=0.6)
        lim = [0, max(max(ref_n_v), max(det_n_v)) * 1.2]
        ax.plot(lim, lim, 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel(r'Greco+2018 S\'ersic $n$', fontsize=11)
        ax.set_ylabel(r'OGFinder S\'ersic $n$', fontsize=11)
        ax.set_title(r'S\'ersic Index', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Ellipticity comparison
    det_e = [s.get('ELLIPTICITY', float('nan')) for s in matched_det]
    ref_e = [s.get('Ell', float('nan')) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_e, ref_e)]
    det_e_v = [d for d, v in zip(det_e, valid) if v]
    ref_e_v = [r for r, v in zip(ref_e, valid) if v]

    if det_e_v:
        ax = axes[2]
        ax.scatter(ref_e_v, det_e_v, s=30, c='dodgerblue', alpha=0.6)
        ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel('Greco+2018 Ellipticity', fontsize=11)
        ax.set_ylabel('OGFinder Ellipticity', fontsize=11)
        ax.set_title('Ellipticity', fontsize=12)
        ax.set_xlim(0, 0.8)
        ax.set_ylim(0, 0.8)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_param_comparison.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def print_method_comparison():
    """Print method comparison: Greco+2018 SExtractor vs OGFinder SEP."""
    print("\n--- Method Comparison: Greco+2018 vs OGFinder ---")
    print(f"{'Parameter':<25} {'Greco+2018':>15} {'OGFinder':>15}")
    print("-" * 58)
    rows = [
        ('Detection', 'SExtractor', 'SEP'),
        ('Pixel scale', '0.168\"/px', '0.168\"/px'),
        ('DETECT_THRESH', '0.7 sigma', '0.7 sigma'),
        ('DETECT_MINAREA', '100 pix', '100 pix'),
        ('Filter kernel', 'gauss 6px FWHM', 'gauss9x9 (5.9px)'),
        ('BACK_SIZE', '128', '128'),
        ('Selection mu_eff', '>24.3 (g)', '>24.3'),
        ('Selection r_eff', '2.5\"-14\"', '2.5\"-14\"'),
        ('Selection ellip', '<0.7', '<0.7'),
        ('Photometry', 'SExtractor AUTO', 'SEP Kron'),
        ('Sersic fitting', 'imfit (2D)', '1D+2D (scipy)'),
        ('ML classification', 'Random Forest', 'None (cuts only)'),
    ]
    for param, greco, ogf in rows:
        print(f"  {param:<23} {greco:>15} {ogf:>15}")


def main():
    if not os.path.isdir(RESULTS_DIR):
        print(f"ERROR: Results directory not found: {RESULTS_DIR}")
        print("Run scripts/reproduce_greco2018.sh first.")
        sys.exit(1)

    print("=" * 64)
    print("Greco et al. (2018, ApJ 857 104) LSBG Comparison")
    print("  Survey: HSC-SSP Wide")
    print("  Pixel scale: 0.168\"/pixel, ZP=27.0")
    print("  Method: SExtractor (Greco) vs SEP (OGFinder)")
    print("=" * 64)

    # Load detection catalog
    det_file = os.path.join(RESULTS_DIR, 'lsbg_detected.tsv')
    if not os.path.exists(det_file):
        print(f"ERROR: Detection catalog not found: {det_file}")
        sys.exit(1)
    det_sources = load_detection_catalog(det_file)
    print(f"\nLoaded detections: {len(det_sources)} LSBG candidates")

    # Add WCS if needed
    fits_path = os.path.join(DATA_DIR, 'hsc_greco2018_i.fits')
    add_wcs_to_detections(det_sources, fits_path)

    # Load reference catalog (patch subset)
    ref_file = os.path.join(DATA_DIR, 'greco2018_lsbg_patch.csv')
    ref_sources = []
    if os.path.exists(ref_file):
        ref_sources = load_reference_catalog(ref_file)
        print(f"Loaded reference: {len(ref_sources)} Greco+2018 LSBGs in patch")
    else:
        # Try full catalog and filter by image extent
        full_ref = os.path.join(DATA_DIR, 'greco2018_lsbg_catalog.csv')
        if os.path.exists(full_ref):
            ref_sources = load_reference_catalog(full_ref)
            print(f"Loaded full reference: {len(ref_sources)} Greco+2018 LSBGs")
            print("  WARNING: Using full catalog (no patch filter)")
        else:
            print(f"WARNING: Reference catalog not found: {ref_file}")
            print("  Run: python3 scripts/fetch_greco2018_data.py")

    # Cross-match
    matched_det, matched_ref, unmatched_det, unmatched_ref = \
        [], [], det_sources, ref_sources
    if ref_sources and any('RA' in s or 'ra' in s for s in det_sources):
        print(f"\nCross-matching (radius={MATCH_RADIUS_ARCSEC}\")")
        matched_det, matched_ref, unmatched_det, unmatched_ref = \
            cross_match(det_sources, ref_sources)
        print(f"  Matched: {len(matched_det)}")
        print(f"  Unmatched detections: {len(unmatched_det)}")
        print(f"  Missed references: {len(unmatched_ref)}")

        # Print match statistics
        if matched_det:
            dists = []
            for det, ref in zip(matched_det, matched_ref):
                d_ra = det.get('RA', det.get('ra', 0))
                d_dec = det.get('DEC', det.get('dec', 0))
                r_ra = ref.get('ra', 0)
                r_dec = ref.get('dec', 0)
                dists.append(angular_separation(d_ra, d_dec, r_ra, r_dec))
            print(f"  Match distances: median={np.median(dists):.2f}\", "
                  f"max={max(dists):.2f}\"")

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
              f"(median {np.median(det_r):.1f}\")")

    grades = {}
    for s in det_sources:
        g = s.get('LSBG_GRADE', '?')
        grades[g] = grades.get(g, 0) + 1
    if grades:
        print(f"  Grades: {grades}")

    # Generate figures
    print("\n--- Generating Figures ---")

    if ref_sources and any('RA' in s or 'ra' in s for s in det_sources):
        print("\nFig 1: Sky map")
        plot_sky_map(det_sources, ref_sources, matched_det, unmatched_ref)

    print("\nFig 2: mu_eff vs r_eff")
    plot_mu_reff(det_sources, ref_sources)

    print("\nFig 3: Sersic n distribution")
    plot_sersic_distribution(det_sources, ref_sources)

    if ref_sources and matched_det:
        print("\nFig 4: Completeness & contamination")
        plot_completeness(det_sources, ref_sources,
                          matched_det, unmatched_det, matched_ref)

    print("\nFig 5: LSBG cutouts")
    plot_cutouts(det_sources, ref_sources, matched_det, unmatched_ref,
                 fits_path)

    if matched_det:
        print("\nFig 6: Parameter comparison (matched)")
        plot_parameter_comparison(matched_det, matched_ref)

    # Method comparison table
    print_method_comparison()

    print("\n" + "=" * 64)
    print("Done. Figures saved to", RESULTS_DIR)
    print("=" * 64)


if __name__ == '__main__':
    main()
