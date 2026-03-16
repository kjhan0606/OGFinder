"""Batch processing for multiple FITS files."""

import os
import sys


def add_batch_args(parser):
    """Add --n-workers, --batch, --output-dir arguments to an ArgumentParser."""
    parser.add_argument('--n-workers', type=int, default=0,
                        help='Parallel workers: 0=auto, 1=sequential, N=N workers')
    parser.add_argument('--batch', nargs='*', default=None,
                        help='Batch mode: list of FITS files to process')
    parser.add_argument('--output-dir', default='.',
                        help='Output directory for batch results (default: .)')


def batch_run(fits_files, process_fn, n_workers=0, output_dir='.', **kwargs):
    """Process multiple FITS files sequentially, each with internal parallelism.

    Parameters
    ----------
    fits_files : list of str
        Paths to FITS files.
    process_fn : callable
        Function(fits_path, n_workers=N, **kwargs) → str (TSV output).
    n_workers : int
        Workers for internal parallelism within each file.
    output_dir : str
        Directory for output TSV files.
    **kwargs
        Extra keyword arguments passed to process_fn.

    MPI extension point: fits_files[rank::size] for inter-node parallelism.
    """
    os.makedirs(output_dir, exist_ok=True)

    for i, fpath in enumerate(fits_files, 1):
        basename = os.path.splitext(os.path.basename(fpath))[0]
        outfile = os.path.join(output_dir, f"{basename}_result.tsv")

        print(f"[{i}/{len(fits_files)}] Processing {fpath}...", file=sys.stderr)

        try:
            result = process_fn(fpath, n_workers=n_workers, **kwargs)
            with open(outfile, 'w') as f:
                f.write(result)
            print(f"  -> {outfile}", file=sys.stderr)
        except Exception as e:
            print(f"  ERROR: {e}", file=sys.stderr)
            continue

    print(f"Batch complete: {len(fits_files)} files", file=sys.stderr)
