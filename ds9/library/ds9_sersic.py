#!/usr/bin/env python3
"""Sérsic profile fitting for each source."""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def _parse_sources(catalog_file, max_sources):
    """Parse catalog file and return list of source info dicts."""
    with open(catalog_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    headers = [h.strip() for h in lines[0].split('\t')]
    col_map = {h: i for i, h in enumerate(headers)}

    sources = []
    for line in lines[1:]:
        if len(sources) >= max_sources:
            break
        fields = line.split('\t')
        try:
            src = {
                'num': int(fields[col_map['NUMBER']].strip()),
                'x': float(fields[col_map['X_IMAGE']].strip()) - 1.0,
                'y': float(fields[col_map['Y_IMAGE']].strip()) - 1.0,
                'a': float(fields[col_map['A_IMAGE']].strip()),
                'b': float(fields[col_map['B_IMAGE']].strip()),
                'theta': float(fields[col_map['THETA_IMAGE']].strip()) * np.pi / 180.0,
                'kron_r': float(fields[col_map['KRON_RADIUS']].strip()) if 'KRON_RADIUS' in col_map else 3.5,
                'flux': float(fields[col_map['FLUX_AUTO']].strip()) if 'FLUX_AUTO' in col_map else 1000.0,
            }
            sources.append(src)
        except (ValueError, IndexError, KeyError):
            continue
    return sources


def _fit_one_source(data, src, cfg):
    """Fit Sérsic profile to a single source. Returns result line string."""
    from sersic_fit.fitter import fit_sersic

    ny, nx = data.shape
    x, y = src['x'], src['y']
    a = src['a']
    kron_r = src['kron_r']

    halfsize = int(cfg.cutout_scale * kron_r * max(a, 3.0))
    halfsize = max(halfsize, 15)
    halfsize = min(halfsize, 150)

    ix, iy = int(x), int(y)
    x0 = max(0, ix - halfsize)
    y0 = max(0, iy - halfsize)
    x1 = min(nx, ix + halfsize + 1)
    y1 = min(ny, iy + halfsize + 1)

    nan_result = {
        'num': src['num'], 'n': np.nan, 're': np.nan, 'Ie': np.nan,
        'ellip': np.nan, 'theta': np.nan, 'chi2': np.nan,
    }

    if x1 - x0 < 10 or y1 - y0 < 10:
        return nan_result

    cutout = data[y0:y1, x0:x1].copy()
    cx = x - x0
    cy = y - y0

    result = fit_sersic(cutout, cx, cy, a, src['b'], src['theta'],
                         src['flux'], cfg)
    if result is None:
        return nan_result

    return {
        'num': src['num'],
        'n': result['n'], 're': result['re'], 'Ie': result['Ie'],
        'ellip': result['ellip'], 'theta': result['theta'],
        'chi2': result['chi2'],
    }


def _worker_fit_source(args):
    """Worker for parallel_map: fit one source from shared data."""
    data_spec, src, cfg = args
    shm = None
    try:
        data, shm = data_spec.attach()
        return _fit_one_source(data, src, cfg)
    except Exception:
        return {
            'num': src.get('num', 0), 'n': np.nan, 're': np.nan,
            'Ie': np.nan, 'ellip': np.nan, 'theta': np.nan, 'chi2': np.nan,
        }
    finally:
        if shm is not None:
            shm.close()


def process_file(fitsfile, catalog_file=None, max_sources=200,
                  mag_zeropoint=25.0, n_workers=0):
    """Process a single FITS file and return TSV string."""
    try:
        import sep
        from astropy.io import fits
    except ImportError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        return ""

    from sersic_fit.config import SersicConfig
    from parallel import parallel_map, resolve_n_workers, SharedArray

    with fits.open(fitsfile) as hdul:
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data.astype(np.float64)
                break
    if data is None:
        print("ERROR: No 2D image data found", file=sys.stderr)
        return ""

    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0
    bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64)
    data_sub = data - bkg

    sources = _parse_sources(catalog_file, max_sources)
    if not sources:
        print("ERROR: No sources in catalog", file=sys.stderr)
        return ""

    cfg = SersicConfig(mag_zeropoint=mag_zeropoint, max_sources=max_sources)
    actual = resolve_n_workers(n_workers)

    print(f"Sérsic fitting {len(sources)} sources...", file=sys.stderr)

    if actual == 1:
        results = [_fit_one_source(data_sub, src, cfg) for src in sources]
    else:
        with SharedArray(data_sub) as spec:
            tasks = [(spec, src, cfg) for src in sources]
            results = parallel_map(_worker_fit_source, tasks,
                                   n_workers=n_workers, label="Sérsic")

    lines = ["NUMBER\tSERSIC_N\tSERSIC_RE\tSERSIC_IE\tSERSIC_ELLIP\tSERSIC_THETA\tSERSIC_CHI2"]
    for r in results:
        def fmt(v):
            return f"{v:.3f}" if np.isfinite(v) else "nan"
        def fmt4(v):
            return f"{v:.4f}" if np.isfinite(v) else "nan"
        def fmt1(v):
            return f"{v:.1f}" if np.isfinite(v) else "nan"
        lines.append(f"{r['num']}\t{fmt(r['n'])}\t{fmt(r['re'])}\t{fmt(r['Ie'])}\t"
                     f"{fmt(r['ellip'])}\t{fmt1(r['theta'])}\t{fmt4(r['chi2'])}")

    print(f"Done: {len(results)} sources fitted", file=sys.stderr)
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(description='Sérsic profile fitting')
    parser.add_argument('fitsfile', nargs='?', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--max-sources', type=int, default=200)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)

    from parallel import add_batch_args
    add_batch_args(parser)

    args = parser.parse_args()

    # Batch mode
    if args.batch is not None:
        from parallel import batch_run
        if not args.batch:
            print("ERROR: --batch requires FITS file list", file=sys.stderr)
            sys.exit(1)
        batch_run(args.batch,
                  lambda f, n_workers=0, **kw: process_file(
                      f, catalog_file=args.catalog,
                      max_sources=args.max_sources,
                      mag_zeropoint=args.mag_zeropoint,
                      n_workers=n_workers),
                  n_workers=args.n_workers,
                  output_dir=args.output_dir)
        return

    if not args.fitsfile:
        print("ERROR: fitsfile required (or use --batch)", file=sys.stderr)
        sys.exit(1)

    output = process_file(args.fitsfile, catalog_file=args.catalog,
                          max_sources=args.max_sources,
                          mag_zeropoint=args.mag_zeropoint,
                          n_workers=args.n_workers)
    if output:
        print(output)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
