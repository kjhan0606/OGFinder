#!/usr/bin/env python3
"""
Download SDSS DR18 SpecPhoto data for photo-z training and validation.

Uses SDSS CasJobs / SkyServer SQL query to get galaxies with:
- Spectroscopic redshifts (spec-z)
- ugriz photometry (modelMag)
- Photometric errors
- Morphological parameters

Output: HDF5 file with features + redshifts ready for training.

Usage:
    python3 scripts/fetch_photoz_training.py [--n-max 50000] [--output photo_z/data/sdss_specphoto.h5]
    python3 scripts/fetch_photoz_training.py --quick  # 10K sample for quick test
"""

import os
import sys
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


SDSS_SKYSERVER_URL = "https://skyserver.sdss.org/dr18/SkyServerWS/SearchTools/SqlSearch"


def build_query(n_max=50000, z_min=0.001, z_max=1.0):
    """Build SDSS SQL query for SpecPhoto galaxies."""
    query = f"""
    SELECT TOP {n_max}
        s.z as spec_z, s.zErr as spec_z_err,
        p.ra, p.dec,
        p.modelMag_u, p.modelMag_g, p.modelMag_r, p.modelMag_i, p.modelMag_z,
        p.modelMagErr_u, p.modelMagErr_g, p.modelMagErr_r,
        p.modelMagErr_i, p.modelMagErr_z,
        p.extinction_u, p.extinction_g, p.extinction_r,
        p.extinction_i, p.extinction_z,
        p.petroR50_r, p.petroR90_r,
        p.fracDeV_r
    FROM SpecPhoto AS s
    JOIN PhotoObj AS p ON s.objID = p.objID
    WHERE s.class = 'GALAXY'
      AND s.z BETWEEN {z_min} AND {z_max}
      AND s.zWarning = 0
      AND p.modelMag_r BETWEEN 14.0 AND 22.5
      AND p.modelMag_g BETWEEN 14.0 AND 25.0
      AND p.modelMag_i BETWEEN 14.0 AND 22.5
      AND p.modelMagErr_r < 0.5
      AND p.modelMagErr_g < 0.5
      AND p.modelMagErr_i < 0.5
    ORDER BY NEWID()
    """
    return query.strip()


def fetch_sdss_sql(query, fmt='csv'):
    """Execute SQL query on SDSS SkyServer."""
    import urllib.request
    import urllib.parse

    params = urllib.parse.urlencode({
        'cmd': query,
        'format': fmt,
    })

    url = f"{SDSS_SKYSERVER_URL}?{params}"

    print(f"Querying SDSS SkyServer (this may take a few minutes)...")
    print(f"URL length: {len(url)} chars")

    req = urllib.request.Request(url)
    req.add_header('User-Agent', 'OGFinder/DS9 photo-z training data fetcher')

    try:
        with urllib.request.urlopen(req, timeout=600) as response:
            data = response.read().decode('utf-8')
    except Exception as e:
        print(f"ERROR: SDSS query failed: {e}", file=sys.stderr)
        print("Try with --n-max 10000 for a smaller sample", file=sys.stderr)
        return None

    return data


def parse_csv(csv_data):
    """Parse SDSS CSV response into numpy arrays."""
    lines = csv_data.strip().split('\n')

    # Skip comment lines and find header
    header = None
    data_start = 0
    for i, line in enumerate(lines):
        line = line.strip()
        if line.startswith('#') or line == '':
            continue
        if header is None:
            header = [h.strip() for h in line.split(',')]
            data_start = i + 1
            break

    if header is None:
        print("ERROR: No header found in SDSS response", file=sys.stderr)
        return None, None

    print(f"  Columns: {len(header)}")
    print(f"  Header: {', '.join(header[:10])}...")

    # Parse data rows
    rows = []
    for line in lines[data_start:]:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        fields = line.split(',')
        if len(fields) != len(header):
            continue
        try:
            row = [float(f.strip()) if f.strip() != '' else np.nan
                   for f in fields]
            rows.append(row)
        except ValueError:
            continue

    if not rows:
        print("ERROR: No data rows parsed", file=sys.stderr)
        return None, None

    data = np.array(rows, dtype=np.float64)
    print(f"  Parsed {data.shape[0]} galaxies")
    return header, data


def prepare_features(header, data):
    """Extract features and labels from raw SDSS data.

    Features:
    - 5 dereddened magnitudes (u, g, r, i, z)
    - 5 magnitude errors
    - 5 presence flags (all 1 for SDSS)
    - 4 color indices (u-g, g-r, r-i, i-z)
    - 1 size (petroR50_r)

    Labels:
    - spec_z
    """
    col = {h: i for i, h in enumerate(header)}

    # Dereddened magnitudes
    bands = ['u', 'g', 'r', 'i', 'z']
    mags = np.zeros((data.shape[0], 5))
    mag_errs = np.zeros((data.shape[0], 5))

    for j, b in enumerate(bands):
        mag_col = f'modelMag_{b}'
        err_col = f'modelMagErr_{b}'
        ext_col = f'extinction_{b}'

        mags[:, j] = data[:, col[mag_col]] - data[:, col[ext_col]]
        mag_errs[:, j] = data[:, col[err_col]]

    # Presence flags (all valid for SDSS ugriz)
    flags = np.ones((data.shape[0], 5))

    # Color indices
    colors = np.zeros((data.shape[0], 4))
    for j in range(4):
        colors[:, j] = mags[:, j] - mags[:, j + 1]

    # Size
    if 'petroR50_r' in col:
        size = data[:, col['petroR50_r']].reshape(-1, 1)
    else:
        size = np.ones((data.shape[0], 1))
    size = np.clip(size, 0.1, 100.0)

    # Assemble features: [mags(5), errors(5), flags(5), colors(4), size(1)] = 20
    features = np.hstack([mags, mag_errs, flags, colors, size])

    # Labels
    spec_z = data[:, col['spec_z']]

    # Feature names
    feature_names = (
        [f'mag_{b}' for b in bands] +
        [f'magerr_{b}' for b in bands] +
        [f'flag_{b}' for b in bands] +
        ['color_ug', 'color_gr', 'color_ri', 'color_iz'] +
        ['petroR50_r']
    )

    # Extra columns for analysis
    extra = {
        'mag_r': mags[:, 2],
    }
    for key in ['ra', 'dec', 'spec_z_err', 'fracDeV_r', 'petroR90_r']:
        if key in col:
            extra[key] = data[:, col[key]]

    return features, spec_z, feature_names, extra


def save_hdf5(output_path, features, spec_z, feature_names, extra):
    """Save training data to HDF5."""
    import h5py

    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)

    with h5py.File(output_path, 'w') as f:
        f.create_dataset('features', data=features.astype(np.float32))
        f.create_dataset('redshift', data=spec_z.astype(np.float32))

        # Normalization stats (from full dataset; training code will recompute on train split)
        f.create_dataset('feat_mean', data=features.mean(axis=0).astype(np.float32))
        f.create_dataset('feat_std', data=features.std(axis=0).astype(np.float32))

        # Metadata
        f.attrs['feature_names'] = feature_names
        f.attrs['n_features'] = features.shape[1]
        f.attrs['n_samples'] = features.shape[0]
        f.attrs['bands'] = ['u', 'g', 'r', 'i', 'z']
        f.attrs['z_min'] = float(spec_z.min())
        f.attrs['z_max'] = float(spec_z.max())
        f.attrs['z_median'] = float(np.median(spec_z))

        # Extra columns
        for key, val in extra.items():
            f.create_dataset(f'extra/{key}', data=val.astype(np.float32))

    print(f"\nSaved to {output_path}")
    print(f"  Samples:  {features.shape[0]}")
    print(f"  Features: {features.shape[1]}")
    print(f"  z range:  [{spec_z.min():.4f}, {spec_z.max():.4f}]")
    print(f"  z median: {np.median(spec_z):.4f}")


def main():
    parser = argparse.ArgumentParser(
        description='Download SDSS SpecPhoto data for photo-z training')
    parser.add_argument('--n-max', type=int, default=50000,
                        help='Max number of galaxies (default: 50000)')
    parser.add_argument('--output', type=str, default=None,
                        help='Output HDF5 path')
    parser.add_argument('--quick', action='store_true',
                        help='Quick test with 10K galaxies')
    parser.add_argument('--z-max', type=float, default=1.0,
                        help='Maximum redshift (default: 1.0)')
    args = parser.parse_args()

    if args.quick:
        args.n_max = 10000

    if args.output is None:
        args.output = os.path.join(_project_root, 'photo_z', 'data',
                                   'sdss_specphoto.h5')

    print(f"Fetching {args.n_max} SDSS SpecPhoto galaxies (z < {args.z_max})...")

    # Build and execute query
    query = build_query(n_max=args.n_max, z_max=args.z_max)
    csv_data = fetch_sdss_sql(query)

    if csv_data is None:
        sys.exit(1)

    # Check for error
    if 'error' in csv_data.lower()[:200]:
        print(f"SDSS returned error:\n{csv_data[:500]}", file=sys.stderr)
        sys.exit(1)

    # Parse
    header, data = parse_csv(csv_data)
    if data is None:
        sys.exit(1)

    # Prepare features
    features, spec_z, feature_names, extra = prepare_features(header, data)

    # Quality cuts
    valid = (
        np.all(np.isfinite(features), axis=1) &
        np.isfinite(spec_z) &
        (spec_z > 0.001) &
        (spec_z < args.z_max)
    )
    features = features[valid]
    spec_z = spec_z[valid]
    extra = {k: v[valid] for k, v in extra.items()}
    print(f"  After quality cuts: {features.shape[0]} galaxies")

    # Save
    save_hdf5(args.output, features, spec_z, feature_names, extra)

    # Print summary stats
    print(f"\n=== Dataset Summary ===")
    print(f"  z distribution:")
    for pct in [10, 25, 50, 75, 90]:
        print(f"    {pct}th percentile: z = {np.percentile(spec_z, pct):.4f}")
    print(f"  r-band magnitude:")
    mag_r = extra['mag_r']
    print(f"    range: [{mag_r.min():.1f}, {mag_r.max():.1f}]")
    print(f"    median: {np.median(mag_r):.1f}")


if __name__ == '__main__':
    main()
