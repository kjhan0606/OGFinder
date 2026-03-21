#!/usr/bin/env python3
"""
Compare OGFinder LSBG detections with Tanoglidis et al. (2021, ApJS 252 18)
reference catalog on a DESI Legacy Survey test patch.

Generates figures:
  - fig_sky_map.png          : RA-Dec sky distribution
  - fig_mu_reff.png          : mu_eff vs r_eff selection diagram
  - fig_sersic_dist.png      : Sersic n distribution histogram
  - fig_completeness.png     : Completeness & contamination vs mu_eff
  - fig_cutouts.png          : Representative matched/missed LSBG cutouts
  - fig_param_comparison.png : Parameter comparison for matched sources

Usage:
  cd OGFinder && python3 scripts/compare_tanoglidis2021.py          # g-only
  cd OGFinder && python3 scripts/compare_tanoglidis2021.py --riz    # r+i+z coadd
"""

import os
import sys
import csv
import math
import argparse
import numpy as np

RESULTS_DIR = "results/tanoglidis2021"
RESULTS_DIR_RIZ = "results/tanoglidis2021_riz"
DATA_DIR = "../data"

# ---- Tanoglidis+2021 reference values ----
TANO_MU_EFF_MIN = 24.2      # mag/arcsec^2 (g-band)
TANO_MU_EFF_MAX = 28.5
TANO_R_EFF_MIN = 2.5        # arcsec
TANO_R_EFF_MAX = 20.0       # arcsec
TANO_TOTAL_LSBGS = 23790    # full catalog

PIXEL_SCALE = 0.263
ZP = 30.0

MATCH_RADIUS_ARCSEC = 3.0


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
    """Load Tanoglidis+2021 reference CSV (patch subset or full)."""
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
            # Normalize column names
            if 'ra_se' in src and 'ra' not in src:
                src['ra'] = src['ra_se']
            if 'dec_se' in src and 'dec' not in src:
                src['dec'] = src['dec_se']
            sources.append(src)
    return sources


def angular_separation(ra1, dec1, ra2, dec2):
    """Angular separation in arcseconds."""
    dra = (ra1 - ra2) * math.cos(math.radians(0.5 * (dec1 + dec2)))
    ddec = dec1 - dec2
    return math.sqrt(dra**2 + ddec**2) * 3600.0


def cross_match(det_sources, ref_sources, radius_arcsec=MATCH_RADIUS_ARCSEC):
    """Cross-match detected sources with reference catalog."""
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
    """Add RA/Dec to detected sources using FITS WCS."""
    if det_sources and 'RA' in det_sources[0]:
        ra_val = det_sources[0]['RA']
        if isinstance(ra_val, (int, float)) and ra_val != 0:
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
        x = src.get('X_IMAGE', 0) - 1
        y = src.get('Y_IMAGE', 0) - 1
        try:
            coord = w.pixel_to_world(x, y)
            src['RA'] = coord.ra.deg
            src['DEC'] = coord.dec.deg
        except Exception:
            pass


def get_ref_mu_eff(src):
    """Get mu_eff(g) from Tanoglidis+2021 reference source."""
    # Direct SExtractor measurement
    return src.get('mu_mean_g_se', src.get('mu_mean_g', 25.0))


def get_ref_reff(src):
    """Get r_eff(g) in arcsec from reference source."""
    return src.get('flux_radius_g_asec', src.get('r_eff_g', 5.0))


def plot_sky_map(det_sources, ref_sources, matched_det, unmatched_ref):
    """Figure 1: RA-Dec sky distribution."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    ref_ra = [s.get('ra', 0) for s in ref_sources]
    ref_dec = [s.get('dec', 0) for s in ref_sources]
    ax.scatter(ref_ra, ref_dec, s=60, c='red', alpha=0.5, marker='s',
               edgecolors='darkred', linewidth=0.5,
               label=f'Tanoglidis+2021 ({len(ref_sources)})')

    det_ra = [s.get('RA', s.get('ra', 0)) for s in det_sources]
    det_dec = [s.get('DEC', s.get('dec', 0)) for s in det_sources]
    ax.scatter(det_ra, det_dec, s=40, c='dodgerblue', alpha=0.6, marker='o',
               edgecolors='navy', linewidth=0.5,
               label=f'OGFinder ({len(det_sources)})')

    match_ra = [s.get('RA', s.get('ra', 0)) for s in matched_det]
    match_dec = [s.get('DEC', s.get('dec', 0)) for s in matched_det]
    if match_ra:
        ax.scatter(match_ra, match_dec, s=100, facecolors='none',
                   edgecolors='green', linewidth=2, marker='o',
                   label=f'Matched ({len(matched_det)})')

    ax.set_xlabel('RA (deg)', fontsize=13)
    ax.set_ylabel('Dec (deg)', fontsize=13)
    ax.set_title('Tanoglidis+2021 LSBG Sky Distribution\n'
                 'DES DR2 g-band (tile DES0125-0124)', fontsize=12)
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
    """Figure 2: mu_eff vs r_eff diagram."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 1, figsize=(9, 7))

    ref_mu = [get_ref_mu_eff(s) for s in ref_sources]
    ref_r = [get_ref_reff(s) for s in ref_sources]
    ax.scatter(ref_r, ref_mu, s=30, c='red', alpha=0.4, marker='s',
               edgecolors='darkred', linewidth=0.3,
               label=f'Tanoglidis+2021 ({len(ref_sources)})')

    det_mu = [s.get('MU_EFF', s.get('mu_eff', 0)) for s in det_sources]
    det_r = [s.get('R_EFF_ARCSEC', s.get('r_eff_arcsec', 0))
             for s in det_sources]
    ax.scatter(det_r, det_mu, s=25, c='dodgerblue', alpha=0.5, marker='o',
               edgecolors='navy', linewidth=0.3,
               label=f'OGFinder ({len(det_sources)})')

    # Selection box
    box_r = [TANO_R_EFF_MIN, TANO_R_EFF_MAX, TANO_R_EFF_MAX,
             TANO_R_EFF_MIN, TANO_R_EFF_MIN]
    box_mu = [TANO_MU_EFF_MIN, TANO_MU_EFF_MIN, TANO_MU_EFF_MAX,
              TANO_MU_EFF_MAX, TANO_MU_EFF_MIN]
    ax.plot(box_r, box_mu, 'k--', linewidth=1.5, alpha=0.7,
            label='Tanoglidis+2021 selection')

    ax.set_xlabel(r'$r_{eff}$ (arcsec)', fontsize=13)
    ax.set_ylabel(r'$\mu_{eff,g}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax.set_title(r'$\mu_{eff}$ vs $r_{eff}$ Diagram', fontsize=12)
    ax.set_xscale('log')
    ax.invert_yaxis()
    ax.set_xlim(1, 40)
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

    ref_n = [s.get('n', float('nan')) for s in ref_sources]
    ref_n = [n for n in ref_n if np.isfinite(n) and 0 < n < 10]
    if ref_n:
        ax.hist(ref_n, bins=bins, color='red', alpha=0.3,
                edgecolor='darkred', linewidth=0.5,
                label=f'Tanoglidis+2021 (N={len(ref_n)}, '
                      f'med={np.median(ref_n):.2f})')

    det_n = [s.get('SERSIC_N', float('nan')) for s in det_sources]
    det_n = [n for n in det_n if np.isfinite(n) and 0 < n < 10]
    if det_n:
        ax.hist(det_n, bins=bins, color='dodgerblue', alpha=0.5,
                edgecolor='navy', linewidth=0.5,
                label=f'OGFinder (N={len(det_n)}, '
                      f'med={np.median(det_n):.2f})')

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
    ref_mu = np.array([get_ref_mu_eff(s) for s in ref_sources])
    matched_ref_mu = np.array([get_ref_mu_eff(s) for s in matched_ref])

    completeness = []
    for lo, hi in zip(mu_bins[:-1], mu_bins[1:]):
        n_ref = np.sum((ref_mu >= lo) & (ref_mu < hi))
        n_match = np.sum((matched_ref_mu >= lo) & (matched_ref_mu < hi)) \
            if len(matched_ref_mu) > 0 else 0
        completeness.append(n_match / n_ref if n_ref > 0 else np.nan)
    completeness = np.array(completeness)

    valid = np.isfinite(completeness)
    ax1.plot(mu_centers[valid], completeness[valid] * 100, 'o-',
             color='navy', ms=6, linewidth=1.5)
    ax1.axhline(50, color='gray', linestyle=':', linewidth=1,
                label='50% completeness')
    ax1.set_xlabel(r'$\mu_{eff,g}$ (mag arcsec$^{-2}$)', fontsize=13)
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
        contamination.append(n_unmatch / n_det if n_det > 0 else np.nan)
    contamination = np.array(contamination)

    valid2 = np.isfinite(contamination)
    ax2.plot(mu_centers[valid2], contamination[valid2] * 100, 's-',
             color='crimson', ms=6, linewidth=1.5)
    ax2.set_xlabel(r'$\mu_{eff,g}$ (mag arcsec$^{-2}$)', fontsize=13)
    ax2.set_ylabel('Contamination (%)', fontsize=13)
    ax2.set_title('False Positive Rate vs Surface Brightness', fontsize=12)
    ax2.set_ylim(-5, 105)
    ax2.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_completeness.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")

    overall_compl = len(matched_det) / max(len(ref_sources), 1) * 100
    overall_contam = len(unmatched_det) / max(len(det_sources), 1) * 100
    print(f"  Overall completeness: {overall_compl:.1f}% "
          f"({len(matched_det)}/{len(ref_sources)})")
    print(f"  Overall contamination: {overall_contam:.1f}% "
          f"({len(unmatched_det)}/{len(det_sources)})")


def plot_cutouts(det_sources, matched_det, fits_path):
    """Figure 5: LSBG cutouts."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    if not os.path.exists(fits_path):
        print("  Cutout figure skipped (FITS not found)")
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
                    print("  No image data")
                    return
        while image_data.ndim > 2:
            image_data = image_data[0]
    except Exception as e:
        print(f"  Cutout figure skipped: {e}")
        return

    ny, nx = image_data.shape
    cutout_size = 80  # pixels (~21" at 0.262"/px)

    fig, axes = plt.subplots(4, 4, figsize=(14, 14))
    fig.suptitle('LSBG Candidates: Matched (top) & Unmatched (bottom)',
                 fontsize=13)

    # Top 2 rows: matched
    show_matched = matched_det[:8]
    for i, src in enumerate(show_matched):
        row, col = divmod(i, 4)
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

    # Bottom 2 rows: unmatched
    unmatched_show = [s for s in det_sources if s not in matched_det][:8]
    for i, src in enumerate(unmatched_show):
        row = 2 + i // 4
        col = i % 4
        if row >= 4:
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

    total_shown = len(show_matched) + len(unmatched_show)
    for idx in range(total_shown, 16):
        row, col = divmod(idx, 4)
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
    ref_r = [get_ref_reff(s) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_r, ref_r)]
    det_r_v = [d for d, v in zip(det_r, valid) if v]
    ref_r_v = [r for r, v in zip(ref_r, valid) if v]

    if det_r_v:
        ax = axes[0]
        ax.scatter(ref_r_v, det_r_v, s=30, c='dodgerblue', alpha=0.6)
        lim = [0, max(max(ref_r_v), max(det_r_v)) * 1.2]
        ax.plot(lim, lim, 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel('Tanoglidis+2021 $r_{eff}$ (arcsec)', fontsize=11)
        ax.set_ylabel('OGFinder $r_{eff}$ (arcsec)', fontsize=11)
        ax.set_title('Effective Radius', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # mu_eff comparison
    det_mu = [s.get('MU_EFF', float('nan')) for s in matched_det]
    ref_mu = [get_ref_mu_eff(s) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_mu, ref_mu)]
    det_mu_v = [d for d, v in zip(det_mu, valid) if v]
    ref_mu_v = [r for r, v in zip(ref_mu, valid) if v]

    if det_mu_v:
        ax = axes[1]
        ax.scatter(ref_mu_v, det_mu_v, s=30, c='dodgerblue', alpha=0.6)
        lim = [23, 29]
        ax.plot(lim, lim, 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel(r'Tanoglidis+2021 $\mu_{eff,g}$', fontsize=11)
        ax.set_ylabel(r'OGFinder $\mu_{eff,g}$', fontsize=11)
        ax.set_title('Mean Surface Brightness', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    # Sersic n comparison
    det_n = [s.get('SERSIC_N', float('nan')) for s in matched_det]
    ref_n = [s.get('n', float('nan')) for s in matched_ref]
    valid = [(np.isfinite(d) and np.isfinite(r)) for d, r in zip(det_n, ref_n)]
    det_n_v = [d for d, v in zip(det_n, valid) if v]
    ref_n_v = [r for r, v in zip(ref_n, valid) if v]

    if det_n_v:
        ax = axes[2]
        ax.scatter(ref_n_v, det_n_v, s=30, c='dodgerblue', alpha=0.6)
        lim = [0, max(max(ref_n_v), max(det_n_v)) * 1.2]
        ax.plot(lim, lim, 'k--', linewidth=1, alpha=0.5, label='1:1')
        ax.set_xlabel(r'Tanoglidis+2021 S\'ersic $n$', fontsize=11)
        ax.set_ylabel(r'OGFinder S\'ersic $n$', fontsize=11)
        ax.set_title(r'S\'ersic Index', fontsize=12)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3)

    out = os.path.join(RESULTS_DIR, 'fig_param_comparison.png')
    fig.tight_layout()
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved: {out}")


def run_comparison(results_dir, ref_file, fits_path, det_image_label):
    """Run comparison for a single detection run."""
    global RESULTS_DIR
    RESULTS_DIR = results_dir

    if not os.path.isdir(results_dir):
        print(f"ERROR: Results directory not found: {results_dir}")
        return None

    # Load detections
    det_file = os.path.join(results_dir, 'lsbg_detected.tsv')
    if not os.path.exists(det_file):
        print(f"ERROR: Detection catalog not found: {det_file}")
        return None
    det_sources = load_detection_catalog(det_file)
    print(f"\nLoaded detections: {len(det_sources)} LSBG candidates")

    # Add WCS
    add_wcs_to_detections(det_sources, fits_path)

    # Load reference
    ref_sources = []
    if os.path.exists(ref_file):
        ref_sources = load_reference_catalog(ref_file)
        print(f"Loaded reference: {len(ref_sources)} LSBGs")
    else:
        print(f"WARNING: Reference catalog not found: {ref_file}")

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

    # Statistics
    print("\n--- Detection Statistics ---")
    det_mu = [s.get('MU_EFF', float('nan')) for s in det_sources]
    det_r = [s.get('R_EFF_ARCSEC', float('nan')) for s in det_sources]
    det_mu = [x for x in det_mu if np.isfinite(x)]
    det_r = [x for x in det_r if np.isfinite(x)]

    if det_mu:
        print(f"  mu_eff range: {min(det_mu):.1f} - {max(det_mu):.1f} "
              f"(median {np.median(det_mu):.1f})")
    if det_r:
        print(f"  r_eff range: {min(det_r):.1f} - {max(det_r):.1f}\" "
              f"(median {np.median(det_r):.1f}\")")

    grades = {}
    for s in det_sources:
        g = s.get('LSBG_GRADE', '?')
        grades[g] = grades.get(g, 0) + 1
    if grades:
        print(f"  Grades: {grades}")

    # Figures
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
    plot_cutouts(det_sources, matched_det, fits_path)

    if matched_det:
        print("\nFig 6: Parameter comparison")
        plot_parameter_comparison(matched_det, matched_ref)

    # Method comparison
    print("\n--- Method Comparison ---")
    print(f"{'Parameter':<25} {'DES/Tanoglidis':>15} {'OGFinder':>15}")
    print("-" * 58)
    rows = [
        ('Detection', 'SourceExtractor', 'SEP'),
        ('Pixel scale', '0.263\"/px', '0.263\"/px'),
        ('Detection image', 'r+i+z coadd', det_image_label),
        ('Classification', 'SVM', 'Photometric cuts'),
        ('Selection mu_eff', '>24.2 (g)', '>24.2 (g)'),
        ('Selection r_eff', '2.5\"-20\"', '2.5\"-20\"'),
        ('Sersic fitting', 'ngmix/GALFIT', '1D+2D (scipy)'),
    ]
    for param, des, ogf in rows:
        print(f"  {param:<23} {des:>15} {ogf:>15}")

    print("\nDone. Figures saved to", results_dir)

    return {
        'n_det': len(det_sources),
        'n_ref': len(ref_sources),
        'n_matched': len(matched_det),
        'n_missed': len(unmatched_ref),
        'completeness': len(matched_det) / max(len(ref_sources), 1) * 100,
        'contamination': len(unmatched_det) / max(len(det_sources), 1) * 100,
    }


def main():
    parser = argparse.ArgumentParser(
        description='Compare LSBG detections with Tanoglidis+2021')
    parser.add_argument('--riz', action='store_true',
                        help='Compare r+i+z coadd detection results')
    args = parser.parse_args()

    fits_path = os.path.join(DATA_DIR, 'des_tanoglidis2021_g_cutout.fits')

    if args.riz:
        # Use full FoV catalog (32 sources) for r+i+z comparison
        ref_file = os.path.join(DATA_DIR, 'tanoglidis2021_lsbg_fov.csv')

        print("=" * 64)
        print("Tanoglidis+2021 LSBG Comparison: g-only vs r+i+z Coadd")
        print("  Survey: DES DR2 (tile DES0125-0124)")
        print("  Pixel scale: 0.263\"/pixel, ZP=30.0")
        print(f"  Reference: {ref_file}")
        print("=" * 64)

        # g-only results
        print("\n" + "=" * 64)
        print("  [A] g-band only detection")
        print("=" * 64)
        g_stats = run_comparison(
            "results/tanoglidis2021", ref_file, fits_path, 'g-band coadd')

        # r+i+z results
        print("\n" + "=" * 64)
        print("  [B] r+i+z coadd detection + g-band forced photometry")
        print("=" * 64)
        riz_fits = os.path.join(DATA_DIR, 'des_tanoglidis2021_riz_coadd.fits')
        riz_fits_path = riz_fits if os.path.exists(riz_fits) else fits_path
        riz_stats = run_comparison(
            RESULTS_DIR_RIZ, ref_file, riz_fits_path, 'r+i+z coadd')

        # Side-by-side comparison
        print("\n" + "=" * 64)
        print("  g-only vs r+i+z Coadd Comparison")
        print("=" * 64)
        print(f"{'Metric':<30} {'g-only':>12} {'r+i+z':>12}")
        print("-" * 56)

        if g_stats and riz_stats:
            for label, key in [
                ('Reference LSBGs', 'n_ref'),
                ('Detected candidates', 'n_det'),
                ('Matched', 'n_matched'),
                ('Missed', 'n_missed'),
            ]:
                g_val = g_stats[key]
                r_val = riz_stats[key]
                print(f"  {label:<28} {g_val:>12} {r_val:>12}")

            print("-" * 56)
            g_c = g_stats['completeness']
            r_c = riz_stats['completeness']
            print(f"  {'Completeness':<28} {g_c:>11.1f}% {r_c:>11.1f}%")
            g_f = g_stats['contamination']
            r_f = riz_stats['contamination']
            print(f"  {'Contamination':<28} {g_f:>11.1f}% {r_f:>11.1f}%")

            delta = r_c - g_c
            print(f"\n  Completeness improvement: {delta:+.1f}% "
                  f"({g_stats['n_matched']}/{g_stats['n_ref']} -> "
                  f"{riz_stats['n_matched']}/{riz_stats['n_ref']})")

            # Combined completeness analysis
            ref_sources = load_reference_catalog(ref_file)
            g_det = load_detection_catalog(
                os.path.join("results/tanoglidis2021", "lsbg_detected.tsv"))
            riz_det = load_detection_catalog(
                os.path.join(RESULTS_DIR_RIZ, "lsbg_detected.tsv"))
            add_wcs_to_detections(g_det, fits_path)

            n_ref = len(ref_sources)
            g_matched = set()
            riz_matched = set()
            for j, ref in enumerate(ref_sources):
                r_ra = ref.get('ra', 0)
                r_dec = ref.get('dec', 0)
                for d in g_det:
                    d_ra = d.get('RA', d.get('ra', 0))
                    d_dec = d.get('DEC', d.get('dec', 0))
                    if angular_separation(d_ra, d_dec, r_ra, r_dec) < 5.0:
                        g_matched.add(j)
                        break
                for d in riz_det:
                    d_ra = d.get('RA', d.get('ra', 0))
                    d_dec = d.get('DEC', d.get('dec', 0))
                    if angular_separation(d_ra, d_dec, r_ra, r_dec) < 5.0:
                        riz_matched.add(j)
                        break

            combined = g_matched | riz_matched
            g_only_unique = g_matched - riz_matched
            riz_only_unique = riz_matched - g_matched

            print(f"\n  --- Combined Analysis (5\" match) ---")
            print(f"  g-only matches:    {len(g_matched):2d}/{n_ref} "
                  f"({100*len(g_matched)/n_ref:.1f}%)")
            print(f"  r+i+z matches:     {len(riz_matched):2d}/{n_ref} "
                  f"({100*len(riz_matched)/n_ref:.1f}%)")
            print(f"  g-only unique:     {len(g_only_unique):2d} "
                  f"(detected by g but not riz)")
            print(f"  r+i+z unique:      {len(riz_only_unique):2d} "
                  f"(detected by riz but not g)")
            print(f"  Combined (union):  {len(combined):2d}/{n_ref} "
                  f"({100*len(combined)/n_ref:.1f}%)")
            missed_both = set(range(n_ref)) - combined
            if missed_both:
                print(f"  Missed by both:    {len(missed_both)} "
                      f"(refs {sorted(i+1 for i in missed_both)})")
        else:
            if not g_stats:
                print("  g-only results not available.")
                print("  Run: bash scripts/reproduce_tanoglidis2021.sh")
            if not riz_stats:
                print("  r+i+z results not available.")
                print("  Run: bash scripts/reproduce_tanoglidis2021_riz.sh")

        print("\n" + "=" * 64)

    else:
        # Original g-only comparison
        ref_file = os.path.join(DATA_DIR, 'tanoglidis2021_lsbg_patch.csv')

        print("=" * 64)
        print("Tanoglidis et al. (2021, ApJS 252 18) LSBG Comparison")
        print("  Survey: DES DR2 g-band coadd (tile DES0125-0124)")
        print("  Pixel scale: 0.263\"/pixel, ZP=30.0")
        print("  Method: SourceExtractor (DES) vs SEP (OGFinder)")
        print("=" * 64)

        stats = run_comparison(
            RESULTS_DIR, ref_file, fits_path, 'g-band coadd')

        print("\n" + "=" * 64)


if __name__ == '__main__':
    main()
