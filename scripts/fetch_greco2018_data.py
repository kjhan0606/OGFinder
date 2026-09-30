#!/usr/bin/env python3
"""
Fetch data for Greco et al. (2018, ApJ 857 104) LSBG reproduction experiment.

Downloads:
  1. VizieR catalog J/ApJ/857/104/table1 (781 LSBGs from HSC-SSP S16A)
  2. HSC-SSP PDR3 coadd image cutout (requires registration)

Usage:
    cd OGFinder && python3 scripts/fetch_greco2018_data.py

Reference: "Illuminating Low Surface Brightness Galaxies with the
            Hyper Suprime-Cam Survey"
Survey: HSC-SSP Wide (S16A, ~200 deg^2)
Pixel scale: 0.168 arcsec/pixel
Depth: ~26 mag (i-band, 5-sigma point source)

Output:
    data/greco2018_lsbg_catalog.csv    - Full LSBG reference catalog (781)
    data/greco2018_lsbg_patch.csv      - Reference LSBGs in test patch
    data/hsc_greco2018_i.fits          - HSC i-band coadd cutout (if credentials)
"""

import os
import sys
import csv
import math
import urllib.request
import urllib.error
import urllib.parse
import json

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data')

# ---- Output paths ----
CATALOG_OUTPUT = os.path.join(DATADIR, 'greco2018_lsbg_catalog.csv')
PATCH_CATALOG_OUTPUT = os.path.join(DATADIR, 'greco2018_lsbg_patch.csv')
FITS_OUTPUT = os.path.join(DATADIR, 'hsc_greco2018_i.fits')

# ---- Test patch parameters ----
# Chosen by densest-region search over the full catalog.
# Will be updated after catalog download.
PATCH_HALFSIZE_DEG = 0.15  # ~18 arcmin half-width -> 36'x36' cutout
PATCH_HALFSIZE_ARCSEC = PATCH_HALFSIZE_DEG * 3600.0

# ---- VizieR catalog URLs ----
# Method 1: Direct CDS FTP (no registration needed)
VIZIER_FTP_URL = 'https://cdsarc.cds.unistra.fr/ftp/J/ApJ/857/104/table1.dat'
VIZIER_README_URL = 'https://cdsarc.cds.unistra.fr/ftp/J/ApJ/857/104/ReadMe'

# Method 2: VizieR TAP (CSV output)
VIZIER_TAP_URL = (
    'http://tapvizier.u-strasbg.fr/TAPVizieR/tap/sync?'
    'REQUEST=doQuery&LANG=ADQL&FORMAT=csv&'
    'QUERY=SELECT+*+FROM+"J/ApJ/857/104/table1"'
)

# ---- HSC-SSP cutout service ----
HSC_CUTOUT_URL = (
    'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout'
)


def download_file(url, output_path, description="file", auth=None):
    """Download a file with progress reporting."""
    if os.path.exists(output_path):
        size_mb = os.path.getsize(output_path) / 1e6
        print(f"  Already exists: {output_path} ({size_mb:.1f} MB)")
        return True

    print(f"  Downloading {description}...")
    print(f"  URL: {url[:120]}...")
    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'OGFinder/1.0'})
        if auth:
            import base64
            credentials = base64.b64encode(
                f"{auth[0]}:{auth[1]}".encode()).decode()
            req.add_header('Authorization', f'Basic {credentials}')

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


def parse_cds_fixed_format(filepath):
    """Parse CDS fixed-format table1.dat using ReadMe column positions.

    Greco+2018 table1.dat columns (from ReadMe):
      Bytes  Format  Units         Label    Description
      1-  3  I3      ---           Seq      Galaxy identifier
      5- 14  F10.6   deg           RAdeg    Right Ascension J2000
     16- 25  F10.6   deg           DEdeg    Declination J2000
     27- 31  F5.2    mag/arcsec2   iSB      i-band central SB
     33- 36  F4.2    mag/arcsec2   e_iSB    uncertainty
     38- 42  F5.2    mag           imag     i-band magnitude
     44- 47  F4.2    mag           e_imag   uncertainty
     49- 53  F5.2    mag           g-r      color
     55- 59  F5.2    mag           g-i      color
     61- 65  F5.2    arcsec        Reff     effective radius
     67- 70  F4.2    arcsec        e_Reff   uncertainty
     72- 75  F4.2    ---           n        Sersic index
     77- 80  F4.2    ---           e_n      uncertainty
     82- 85  F4.2    ---           Ell      ellipticity
     87- 90  F4.2    ---           e_Ell    uncertainty
     92- 95  F4.2    mag           Ag       extinction g
     97-100  F4.2    mag           Ar       extinction r
    102-105  F4.2    mag           Ai       extinction i
    """
    sources = []
    with open(filepath) as f:
        for line in f:
            line = line.rstrip('\n')
            if len(line) < 100:
                continue
            try:
                src = {
                    'Seq': int(line[0:3].strip()),
                    'ra': float(line[4:14].strip()),
                    'dec': float(line[15:25].strip()),
                    'iSB': float(line[26:31].strip()),
                    'e_iSB': float(line[32:36].strip()),
                    'imag': float(line[37:42].strip()),
                    'e_imag': float(line[43:47].strip()),
                    'g_r': float(line[48:53].strip()),
                    'g_i': float(line[54:59].strip()),
                    'Reff': float(line[60:65].strip()),
                    'e_Reff': float(line[66:70].strip()),
                    'n': float(line[71:75].strip()),
                    'e_n': float(line[76:80].strip()),
                    'Ell': float(line[81:85].strip()),
                    'e_Ell': float(line[86:90].strip()),
                    'Ag': float(line[91:95].strip()),
                    'Ar': float(line[96:100].strip()),
                    'Ai': float(line[101:105].strip()),
                }
                sources.append(src)
            except (ValueError, IndexError):
                continue
    return sources


def download_catalog_vizier_tap():
    """Download catalog via VizieR TAP service (CSV format)."""
    print("\n=== Step 1: Download Greco+2018 LSBG catalog ===")

    if os.path.exists(CATALOG_OUTPUT):
        print(f"  Already exists: {CATALOG_OUTPUT}")
        return True

    # Method 1: Try astroquery (best column handling)
    try:
        from astroquery.vizier import Vizier
        print("  Using astroquery.vizier...")
        Vizier.ROW_LIMIT = -1
        result = Vizier(columns=['**']).get_catalogs('J/ApJ/857/104/table1')
        if result:
            table = result[0]
            # Write as CSV
            fieldnames = ['Seq', 'ra', 'dec', 'iSB', 'e_iSB', 'imag',
                         'e_imag', 'g_r', 'g_i', 'Reff', 'e_Reff',
                         'n', 'e_n', 'Ell', 'e_Ell', 'Ag', 'Ar', 'Ai']
            col_map = {
                'Seq': 'Seq', 'ra': 'RAdeg', 'dec': 'DEdeg',
                'iSB': 'iSB', 'e_iSB': 'e_iSB', 'imag': 'imag',
                'e_imag': 'e_imag', 'g_r': 'g-r', 'g_i': 'g-i',
                'Reff': 'Reff', 'e_Reff': 'e_Reff', 'n': 'n', 'e_n': 'e_n',
                'Ell': 'Ell', 'e_Ell': 'e_Ell', 'Ag': 'Ag', 'Ar': 'Ar',
                'Ai': 'Ai'
            }
            # Try alternate column names
            alt_names = {'RAdeg': ['RAJ2000', '_RA', 'RA_ICRS'],
                        'DEdeg': ['DEJ2000', '_DE', 'DE_ICRS']}

            actual_cols = {}
            for our_name, viz_name in col_map.items():
                if viz_name in table.colnames:
                    actual_cols[our_name] = viz_name
                elif viz_name in alt_names:
                    for alt in alt_names[viz_name]:
                        if alt in table.colnames:
                            actual_cols[our_name] = alt
                            break

            with open(CATALOG_OUTPUT, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for row in table:
                    src = {}
                    for our_name in fieldnames:
                        viz_col = actual_cols.get(our_name, our_name)
                        try:
                            val = float(row[viz_col])
                            if our_name in ('Seq',):
                                src[our_name] = str(int(val))
                            else:
                                src[our_name] = f"{val:.6f}"
                        except (KeyError, ValueError, TypeError):
                            src[our_name] = ''
                    writer.writerow(src)

            print(f"  Retrieved {len(table)} sources -> {CATALOG_OUTPUT}")
            return True
    except ImportError:
        print("  astroquery not available, trying CDS FTP...")
    except Exception as e:
        print(f"  astroquery failed: {e}, trying CDS FTP...")

    # Method 2: Direct CDS FTP download (fixed-format)
    tmp_dat = CATALOG_OUTPUT + '.dat'
    if download_file(VIZIER_FTP_URL, tmp_dat, "CDS table1.dat"):
        print("  Parsing CDS fixed-format catalog...")
        sources = parse_cds_fixed_format(tmp_dat)
        if sources:
            fieldnames = list(sources[0].keys())
            with open(CATALOG_OUTPUT, 'w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=fieldnames)
                writer.writeheader()
                for src in sources:
                    fmt = {}
                    for k, v in src.items():
                        if isinstance(v, float):
                            fmt[k] = f"{v:.6f}"
                        else:
                            fmt[k] = str(v)
                    writer.writerow(fmt)
            print(f"  Parsed {len(sources)} sources -> {CATALOG_OUTPUT}")
            os.remove(tmp_dat)
            return True
        os.remove(tmp_dat)

    # Method 3: VizieR TAP
    print("  Trying VizieR TAP service...")
    return download_file(VIZIER_TAP_URL, CATALOG_OUTPUT, "VizieR TAP CSV")


def load_catalog(path):
    """Load catalog CSV and normalize column names."""
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
            # Normalize VizieR column names to our convention
            if 'RAJ2000' in src and 'ra' not in src:
                src['ra'] = src['RAJ2000']
            if 'DEJ2000' in src and 'dec' not in src:
                src['dec'] = src['DEJ2000']
            if 'RAdeg' in src and 'ra' not in src:
                src['ra'] = src['RAdeg']
            if 'DEdeg' in src and 'dec' not in src:
                src['dec'] = src['DEdeg']
            # Normalize color columns
            if 'g-r' in src and 'g_r' not in src:
                src['g_r'] = src['g-r']
            if 'g-i' in src and 'g_i' not in src:
                src['g_i'] = src['g-i']
            sources.append(src)
    return sources


def find_densest_region(sources, halfsize_deg=PATCH_HALFSIZE_DEG):
    """Find the densest square region in the catalog.

    Grid-searches RA, Dec to find the tile with the most LSBGs.

    Returns (ra_center, dec_center, count, sources_in_patch).
    """
    print("\n=== Step 2: Find densest test region ===")

    ras = [s['ra'] for s in sources]
    decs = [s['dec'] for s in sources]
    ra_min, ra_max = min(ras), max(ras)
    dec_min, dec_max = min(decs), max(decs)

    print(f"  Catalog extent: RA [{ra_min:.2f}, {ra_max:.2f}], "
          f"Dec [{dec_min:.2f}, {dec_max:.2f}]")

    # Grid search
    step = halfsize_deg * 0.5  # overlap tiles
    best_ra, best_dec, best_count = 0, 0, 0
    best_sources = []

    ra = ra_min + halfsize_deg
    while ra < ra_max - halfsize_deg:
        dec = dec_min + halfsize_deg
        while dec < dec_max - halfsize_deg:
            patch = [s for s in sources
                     if abs(s['ra'] - ra) <= halfsize_deg
                     and abs(s['dec'] - dec) <= halfsize_deg]
            if len(patch) > best_count:
                best_count = len(patch)
                best_ra, best_dec = ra, dec
                best_sources = patch
            dec += step
        ra += step

    print(f"  Best region: RA={best_ra:.4f}, Dec={best_dec:.4f}")
    print(f"  LSBGs in region: {best_count}")
    print(f"  Region size: {2*halfsize_deg:.2f} x {2*halfsize_deg:.2f} deg "
          f"({2*PATCH_HALFSIZE_ARCSEC/60:.0f}' x {2*PATCH_HALFSIZE_ARCSEC/60:.0f}')")

    if best_count < 3:
        # Expand search area
        print("  WARNING: Few LSBGs found. Expanding search to 0.3 deg...")
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


def download_hsc_image(ra_center, dec_center):
    """Download HSC-SSP i-band coadd cutout.

    Requires HSC-SSP credentials via environment variables:
      HSC_SSP_USER, HSC_SSP_PASS
    or via ~/.netrc entry for hsc-release.mtk.nao.ac.jp
    """
    print("\n=== Step 3: Download HSC-SSP i-band image ===")

    if os.path.exists(FITS_OUTPUT):
        size_mb = os.path.getsize(FITS_OUTPUT) / 1e6
        print(f"  Already exists: {FITS_OUTPUT} ({size_mb:.1f} MB)")
        return True

    # Check for credentials
    user = os.environ.get('HSC_SSP_USER', '')
    passwd = os.environ.get('HSC_SSP_PASS', '')

    if not user or not passwd:
        # Try SSP_PDR_USR/SSP_PDR_PWD (unagi convention)
        user = os.environ.get('SSP_PDR_USR', '')
        passwd = os.environ.get('SSP_PDR_PWD', '')

    if not user or not passwd:
        print("  HSC-SSP credentials not found.")
        print("  To download automatically, set environment variables:")
        print("    export HSC_SSP_USER='your_username'")
        print("    export HSC_SSP_PASS='your_password'")
        print("")
        print("  Register (free) at:")
        print("    https://hsc-release.mtk.nao.ac.jp/datasearch/new_user/new")
        print("")
        print("  Or download manually with curl:")
        sw = int(PATCH_HALFSIZE_ARCSEC)
        print(f"    curl -u 'USER:PASS' \\")
        print(f"      '{HSC_CUTOUT_URL}?ra={ra_center}&dec={dec_center}"
              f"&sw={sw}arcsec&sh={sw}arcsec"
              f"&type=coadd&image=on&mask=off&variance=off"
              f"&filter=HSC-I&rerun=pdr3_wide' \\")
        print(f"      -o {FITS_OUTPUT}")
        print("")
        print("  After downloading, re-run this script or proceed to:")
        print("    bash scripts/reproduce_greco2018.sh")
        return False

    # Build cutout URL
    sw = int(PATCH_HALFSIZE_ARCSEC)
    params = {
        'ra': f'{ra_center:.6f}',
        'dec': f'{dec_center:.6f}',
        'sw': f'{sw}arcsec',
        'sh': f'{sw}arcsec',
        'type': 'coadd',
        'image': 'on',
        'mask': 'off',
        'variance': 'off',
        'filter': 'HSC-I',
        'rerun': 'pdr3_wide',
    }
    url = HSC_CUTOUT_URL + '?' + urllib.parse.urlencode(params)

    return download_file(url, FITS_OUTPUT, "HSC-SSP i-band coadd",
                        auth=(user, passwd))


def print_patch_summary(sources):
    """Print statistics of the test patch LSBGs."""
    if not sources:
        return

    ras = [s['ra'] for s in sources]
    decs = [s['dec'] for s in sources]
    reffs = [s['Reff'] for s in sources]
    ns = [s['n'] for s in sources]
    ells = [s['Ell'] for s in sources]
    isbs = [s['iSB'] for s in sources]
    imags = [s['imag'] for s in sources]

    print("\n--- Test Patch LSBG Statistics ---")
    print(f"  N sources: {len(sources)}")
    print(f"  RA range:  [{min(ras):.4f}, {max(ras):.4f}]")
    print(f"  Dec range: [{min(decs):.4f}, {max(decs):.4f}]")
    print(f"  R_eff:     [{min(reffs):.2f}, {max(reffs):.2f}] arcsec "
          f"(median {sorted(reffs)[len(reffs)//2]:.2f})")
    print(f"  Sersic n:  [{min(ns):.2f}, {max(ns):.2f}] "
          f"(median {sorted(ns)[len(ns)//2]:.2f})")
    print(f"  Ellip:     [{min(ells):.2f}, {max(ells):.2f}]")
    print(f"  iSB (mu_0):[{min(isbs):.2f}, {max(isbs):.2f}] mag/arcsec^2")
    print(f"  i_mag:     [{min(imags):.2f}, {max(imags):.2f}]")

    # Compute approximate mu_eff from Sersic parameters
    # mu_eff = mu_0 + 2.5*log10(e) * b_n for exponential (n=1): ~1.82
    # More generally: mu_eff = mu_0 + 2.5*b_n/ln(10) where b_n ≈ 2n - 1/3
    print("\n  Approximate mu_eff,i distribution:")
    mu_effs = []
    for s in sources:
        n = s['n']
        mu0 = s['iSB']
        # b_n approximation (Ciotti & Bertin 1999)
        bn = 2.0 * n - 1.0/3.0 + 4.0/(405.0*n) if n > 0.36 else 0.01
        mu_eff = mu0 + 2.5 * bn / math.log(10)
        mu_effs.append(mu_eff)
    mu_effs.sort()
    print(f"    Range: [{mu_effs[0]:.2f}, {mu_effs[-1]:.2f}] mag/arcsec^2")
    print(f"    Median: {mu_effs[len(mu_effs)//2]:.2f}")


def main():
    os.makedirs(DATADIR, exist_ok=True)

    print("=" * 64)
    print("Greco et al. (2018, ApJ 857 104) Data Preparation")
    print("  Survey: HSC-SSP Wide (S16A -> PDR3)")
    print("  Pixel scale: 0.168\"/pixel")
    print("  Depth: ~26 mag (i-band, 5-sigma)")
    print("=" * 64)

    # Step 1: Download VizieR catalog
    if not download_catalog_vizier_tap():
        print("\nERROR: Failed to download catalog. Check network.")
        sys.exit(1)

    # Load catalog
    sources = load_catalog(CATALOG_OUTPUT)
    print(f"  Total LSBGs in catalog: {len(sources)}")

    # Step 2: Find densest region
    ra_c, dec_c, n_patch, patch_sources = find_densest_region(sources)

    # Print patch summary
    print_patch_summary(patch_sources)

    # Step 3: Download HSC-SSP image
    download_hsc_image(ra_c, dec_c)

    print("\n" + "=" * 64)
    print("Data preparation complete.")
    print(f"  Reference catalog: {CATALOG_OUTPUT} ({len(sources)} LSBGs)")
    print(f"  Test patch: RA={ra_c:.4f}, Dec={dec_c:.4f} ({n_patch} LSBGs)")
    print(f"  Patch catalog: {PATCH_CATALOG_OUTPUT}")
    if os.path.exists(FITS_OUTPUT):
        print(f"  Image: {FITS_OUTPUT}")
    print("")
    print("Next: bash scripts/reproduce_greco2018.sh")
    print("=" * 64)


if __name__ == '__main__':
    main()
