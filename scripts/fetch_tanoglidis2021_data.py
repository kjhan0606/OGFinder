#!/usr/bin/env python3
"""
Fetch data for Tanoglidis et al. (2021, ApJS 252 18) LSBG reproduction.

Downloads:
  1. DES Y3 LSBG catalog (23,790 sources) from DES data server
  2. DESI Legacy Survey DR10 g-band image (no registration needed)

Usage:
    cd OGFinder && python3 scripts/fetch_tanoglidis2021_data.py

Reference: "Shadows in the Dark: Low-Surface-Brightness Galaxies Discovered
            in the Dark Energy Survey"
Survey: DES Year 3, ~5000 deg^2
Detection: SourceExtractor → SVM classification
Selection: mu_eff(g) > 24.2, 2.5" < r_eff < 20"
"""

import os
import sys
import csv
import math
import urllib.request
import urllib.error
import urllib.parse

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data')

# ---- Output paths ----
CATALOG_OUTPUT = os.path.join(DATADIR, 'tanoglidis2021_lsbg_catalog.csv')
PATCH_CATALOG_OUTPUT = os.path.join(DATADIR, 'tanoglidis2021_lsbg_patch.csv')
FITS_OUTPUT = os.path.join(DATADIR, 'legacy_tanoglidis2021_g.fits')

# ---- Patch parameters ----
PATCH_HALFSIZE_DEG = 0.15  # ~18' half-width -> 36'x36' cutout

# ---- Data URLs ----
DES_CATALOG_URL = (
    'https://desdr-server.ncsa.illinois.edu/despublic/other_files/'
    'y3-lsbg/LSBG_catalog_v2.csv'
)

# Legacy Survey cutout API (no registration)
LEGACY_CUTOUT_URL = 'https://www.legacysurvey.org/viewer/fits-cutout'


def download_file(url, output_path, description="file"):
    """Download a file with progress reporting."""
    if os.path.exists(output_path):
        size_mb = os.path.getsize(output_path) / 1e6
        print(f"  Already exists: {output_path} ({size_mb:.1f} MB)")
        return True

    print(f"  Downloading {description}...")
    print(f"  URL: {url[:120]}...")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'OGFinder/1.0'})
        with urllib.request.urlopen(req, timeout=300) as resp:
            total = resp.headers.get('Content-Length')
            total = int(total) if total else None

            tmp_path = output_path + '.tmp'
            downloaded = 0
            with open(tmp_path, 'wb') as f:
                while True:
                    chunk = resp.read(1024 * 1024)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded += len(chunk)
                    if total:
                        pct = 100.0 * downloaded / total
                        print(f"\r  {downloaded/1e6:.1f} / {total/1e6:.1f} MB "
                              f"({pct:.0f}%)", end='', flush=True)
                    else:
                        print(f"\r  {downloaded/1e6:.1f} MB downloaded",
                              end='', flush=True)
            print()

            os.rename(tmp_path, output_path)
            print(f"  Saved: {output_path} ({downloaded/1e6:.1f} MB)")
            return True

    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        print(f"  ERROR: {e}")
        tmp_path = output_path + '.tmp'
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
        return False


def load_catalog(path):
    """Load DES LSBG catalog CSV.

    Key columns from LSBG_catalog_v2.csv:
      ra_se, dec_se          - SExtractor coordinates
      mu_mean_g_se           - g-band mean SB (SExtractor)
      flux_radius_g_asec     - g-band half-light radius (arcsec, SExtractor)
      ell_se                 - ellipticity (SExtractor)
      n                      - Sersic index (GALFIT-like)
      r_eff_g                - Sersic r_eff in g (arcsec)
      mu_mean_g, mu_0_g      - Sersic mean/central SB
      mag_auto_g_corr        - g-band AUTO mag (extinction-corrected)
    """
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
            # Normalize to common column names
            if 'ra_se' in src:
                src['ra'] = src['ra_se']
            if 'dec_se' in src:
                src['dec'] = src['dec_se']
            sources.append(src)
    return sources


def find_densest_region(sources, halfsize_deg=PATCH_HALFSIZE_DEG):
    """Find the densest square region in the catalog.

    Uses numpy histogram2d for fast density estimation, then refines
    around the densest cell.
    """
    print("\n=== Step 2: Find densest test region ===")
    import numpy as np

    ras = np.array([s['ra'] for s in sources])
    decs = np.array([s['dec'] for s in sources])
    ra_min, ra_max = float(ras.min()), float(ras.max())
    dec_min, dec_max = float(decs.min()), float(decs.max())

    print(f"  Catalog extent: RA [{ra_min:.2f}, {ra_max:.2f}], "
          f"Dec [{dec_min:.2f}, {dec_max:.2f}]")
    print(f"  Total sources: {len(sources)}")

    # Phase 1: Coarse 2D histogram to find dense cells
    cell_size = 2 * halfsize_deg  # bin size = patch size
    ra_bins = np.arange(ra_min, ra_max + cell_size, cell_size)
    dec_bins = np.arange(dec_min, dec_max + cell_size, cell_size)
    hist, _, _ = np.histogram2d(ras, decs, bins=[ra_bins, dec_bins])

    # Find top-10 cells
    flat_idx = np.argsort(hist.ravel())[::-1][:10]
    candidates = []
    for idx in flat_idx:
        ri, di = np.unravel_index(idx, hist.shape)
        ra_c = 0.5 * (ra_bins[ri] + ra_bins[ri + 1])
        dec_c = 0.5 * (dec_bins[di] + dec_bins[di + 1])
        candidates.append((ra_c, dec_c, int(hist[ri, di])))

    # Phase 2: Refine around top candidates with finer grid
    best_ra, best_dec, best_count = 0, 0, 0
    step = halfsize_deg * 0.25
    for cra, cdec, _ in candidates:
        for dra in np.arange(-halfsize_deg, halfsize_deg + step, step):
            for ddec in np.arange(-halfsize_deg, halfsize_deg + step, step):
                ra_t = cra + dra
                dec_t = cdec + ddec
                mask = ((np.abs(ras - ra_t) <= halfsize_deg) &
                        (np.abs(decs - dec_t) <= halfsize_deg))
                count = int(mask.sum())
                if count > best_count:
                    best_count = count
                    best_ra, best_dec = ra_t, dec_t

    print(f"  Best region: RA={best_ra:.4f}, Dec={best_dec:.4f}")
    print(f"  LSBGs in region: {best_count}")
    size_arcmin = 2 * halfsize_deg * 60
    print(f"  Region size: {size_arcmin:.0f}' x {size_arcmin:.0f}' "
          f"({2*halfsize_deg:.2f} deg)")

    # Collect sources in the best patch
    mask = ((np.abs(ras - best_ra) <= halfsize_deg) &
            (np.abs(decs - best_dec) <= halfsize_deg))
    best_sources = [s for s, m in zip(sources, mask) if m]

    if best_count < 5:
        print("  WARNING: Few LSBGs found. Expanding to 0.3 deg...")
        return find_densest_region(sources, halfsize_deg=0.3)

    # Save patch catalog
    if best_sources:
        fieldnames = list(best_sources[0].keys())
        with open(PATCH_CATALOG_OUTPUT, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for src in best_sources:
                fmt = {k: (f"{v:.6f}" if isinstance(v, float) else str(v))
                       for k, v in src.items()}
                writer.writerow(fmt)
        print(f"  Patch catalog: {PATCH_CATALOG_OUTPUT}")

    return best_ra, best_dec, best_count, best_sources


def download_legacy_image(ra_center, dec_center, halfsize_deg):
    """Download DESI Legacy Survey g-band coadd cutout."""
    print("\n=== Step 3: Download Legacy Survey g-band image ===")

    if os.path.exists(FITS_OUTPUT):
        size_mb = os.path.getsize(FITS_OUTPUT) / 1e6
        print(f"  Already exists: {FITS_OUTPUT} ({size_mb:.1f} MB)")
        return True

    # Calculate pixel size for the cutout
    # Legacy Survey pixel scale = 0.262 arcsec/pixel
    pixel_scale = 0.262
    size_arcsec = 2 * halfsize_deg * 3600.0
    size_pixels = int(size_arcsec / pixel_scale)

    # Limit to reasonable cutout size (Legacy Survey limit ~3000-5000 px)
    if size_pixels > 3000:
        print(f"  Requested {size_pixels}px, limiting to 3000px "
              f"({3000*pixel_scale/60:.1f}' x {3000*pixel_scale/60:.1f}')")
        size_pixels = 3000

    params = {
        'ra': f'{ra_center:.6f}',
        'dec': f'{dec_center:.6f}',
        'layer': 'ls-dr10',
        'pixscale': f'{pixel_scale}',
        'size': str(size_pixels),
        'bands': 'g',
    }
    url = LEGACY_CUTOUT_URL + '?' + urllib.parse.urlencode(params)

    print(f"  Cutout: {size_pixels}x{size_pixels} px "
          f"({size_pixels*pixel_scale/60:.1f}' x {size_pixels*pixel_scale/60:.1f}')")
    print(f"  Center: RA={ra_center:.4f}, Dec={dec_center:.4f}")
    print(f"  Layer: ls-dr10, Band: g, Pixscale: {pixel_scale}\"/px")

    return download_file(url, FITS_OUTPUT, "Legacy Survey g-band coadd")


def print_patch_summary(sources):
    """Print statistics of the test patch LSBGs."""
    if not sources:
        return

    import numpy as np

    print("\n--- Test Patch LSBG Statistics ---")
    print(f"  N sources: {len(sources)}")

    ras = [s['ra'] for s in sources]
    decs = [s['dec'] for s in sources]
    print(f"  RA range:  [{min(ras):.4f}, {max(ras):.4f}]")
    print(f"  Dec range: [{min(decs):.4f}, {max(decs):.4f}]")

    # SExtractor measurements
    mu_g = [s.get('mu_mean_g_se', float('nan')) for s in sources]
    mu_g = [x for x in mu_g if not math.isnan(x)]
    if mu_g:
        print(f"  mu_mean_g(SE): [{min(mu_g):.2f}, {max(mu_g):.2f}] "
              f"(median {sorted(mu_g)[len(mu_g)//2]:.2f})")

    r_eff = [s.get('flux_radius_g_asec', float('nan')) for s in sources]
    r_eff = [x for x in r_eff if not math.isnan(x)]
    if r_eff:
        print(f"  r_eff_g(SE):   [{min(r_eff):.2f}, {max(r_eff):.2f}] arcsec "
              f"(median {sorted(r_eff)[len(r_eff)//2]:.2f}\")")

    ell = [s.get('ell_se', float('nan')) for s in sources]
    ell = [x for x in ell if not math.isnan(x)]
    if ell:
        print(f"  Ellipticity:   [{min(ell):.2f}, {max(ell):.2f}]")

    # Sersic measurements
    n_vals = [s.get('n', float('nan')) for s in sources]
    n_vals = [x for x in n_vals if not math.isnan(x) and 0 < x < 10]
    if n_vals:
        print(f"  Sersic n:      [{min(n_vals):.2f}, {max(n_vals):.2f}] "
              f"(median {sorted(n_vals)[len(n_vals)//2]:.2f})")

    mag_g = [s.get('mag_auto_g_corr', float('nan')) for s in sources]
    mag_g = [x for x in mag_g if not math.isnan(x)]
    if mag_g:
        print(f"  g_mag(AUTO):   [{min(mag_g):.2f}, {max(mag_g):.2f}] "
              f"(median {sorted(mag_g)[len(mag_g)//2]:.2f})")


def main():
    os.makedirs(DATADIR, exist_ok=True)

    print("=" * 64)
    print("Tanoglidis et al. (2021, ApJS 252 18) Data Preparation")
    print("  Survey: DES Year 3 (~5000 deg^2)")
    print("  Image: DESI Legacy Survey DR10 (0.262\"/px, ZP=22.5)")
    print("  Selection: mu_eff(g) > 24.2, r_eff 2.5\"-20\"")
    print("=" * 64)

    # Step 1: Download DES LSBG catalog
    print("\n=== Step 1: Download DES Y3 LSBG catalog ===")
    if not download_file(DES_CATALOG_URL, CATALOG_OUTPUT,
                        "DES Y3 LSBG catalog v2"):
        print("\nERROR: Failed to download catalog.")
        sys.exit(1)

    # Load catalog
    print("  Loading catalog...")
    sources = load_catalog(CATALOG_OUTPUT)
    print(f"  Total LSBGs: {len(sources)}")

    # Step 2: Find densest region
    ra_c, dec_c, n_patch, patch_sources = find_densest_region(sources)

    # Print patch summary
    print_patch_summary(patch_sources)

    # Step 3: Download Legacy Survey image
    download_legacy_image(ra_c, dec_c, PATCH_HALFSIZE_DEG)

    print("\n" + "=" * 64)
    print("Data preparation complete.")
    print(f"  Reference catalog: {CATALOG_OUTPUT} ({len(sources)} LSBGs)")
    print(f"  Test patch: RA={ra_c:.4f}, Dec={dec_c:.4f} ({n_patch} LSBGs)")
    print(f"  Patch catalog: {PATCH_CATALOG_OUTPUT}")
    if os.path.exists(FITS_OUTPUT):
        size_mb = os.path.getsize(FITS_OUTPUT) / 1e6
        print(f"  Image: {FITS_OUTPUT} ({size_mb:.1f} MB)")
    print("")
    print("Next: bash scripts/reproduce_tanoglidis2021.sh")
    print("=" * 64)


if __name__ == '__main__':
    main()
