"""Generate labeled training data by cross-matching detections with a reference catalog.

Usage:
    python3 -m lsbg.training.make_labels \
        --detections results/tanoglidis2021/lsbg_detected.tsv \
        --reference ../data/tanoglidis2021_lsbg_patch.csv \
        --ra-col ra_se --dec-col dec_se \
        --match-radius 3.0 \
        --output lsbg/data/training/labeled.tsv
"""

import argparse
import csv
import math
import os
import sys


def angular_separation(ra1, dec1, ra2, dec2):
    """Angular separation in arcseconds (small-angle approximation)."""
    dra = (ra1 - ra2) * math.cos(math.radians(0.5 * (dec1 + dec2)))
    ddec = dec1 - dec2
    return math.sqrt(dra**2 + ddec**2) * 3600.0


def load_tsv(path):
    """Load TSV file as list of dicts."""
    with open(path, 'r') as f:
        lines = f.read().strip().split('\n')
    if len(lines) < 2:
        return [], []
    headers = lines[0].split('\t')
    rows = []
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) < len(headers):
            continue
        row = {}
        for h, v in zip(headers, fields):
            row[h.strip()] = v.strip()
        rows.append(row)
    return headers, rows


def load_csv(path):
    """Load CSV file as list of dicts."""
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        rows = [row for row in reader]
    return rows


def get_float(row, *keys):
    """Get float value from row, trying multiple key names."""
    for key in keys:
        val = row.get(key)
        if val is not None:
            try:
                v = float(val)
                if math.isfinite(v):
                    return v
            except (ValueError, TypeError):
                continue
    return None


def cross_match(detections, references, ra_col, dec_col, match_radius):
    """Cross-match detections with references by RA/DEC.

    Parameters
    ----------
    detections : list of dict
        Detection catalog rows.
    references : list of dict
        Reference catalog rows.
    ra_col, dec_col : str
        Column names for RA/DEC in the reference catalog.
    match_radius : float
        Match radius in arcseconds.

    Returns
    -------
    labels : list of int
        1 = matched (LSBG), 0 = unmatched (contaminant).
    match_info : list of str
        Match details for each detection.
    """
    labels = []
    match_info = []
    n_matched = 0

    for det in detections:
        det_ra = get_float(det, 'RA', 'ra')
        det_dec = get_float(det, 'DEC', 'dec')

        if det_ra is None or det_dec is None:
            labels.append(0)
            match_info.append('no_coords')
            continue

        best_sep = float('inf')
        for ref in references:
            ref_ra = get_float(ref, ra_col)
            ref_dec = get_float(ref, dec_col)
            if ref_ra is None or ref_dec is None:
                continue
            sep = angular_separation(det_ra, det_dec, ref_ra, ref_dec)
            if sep < best_sep:
                best_sep = sep

        if best_sep <= match_radius:
            labels.append(1)
            match_info.append(f'matched_{best_sep:.2f}')
            n_matched += 1
        else:
            labels.append(0)
            match_info.append(f'unmatched_{best_sep:.2f}')

    return labels, match_info, n_matched


def main():
    parser = argparse.ArgumentParser(
        description='Generate labeled data for LSBG SVM training')
    parser.add_argument('--detections', required=True,
                        help='Detection catalog TSV')
    parser.add_argument('--reference', required=True,
                        help='Reference catalog (CSV or TSV)')
    parser.add_argument('--ra-col', default='ra_se',
                        help='RA column name in reference catalog')
    parser.add_argument('--dec-col', default='dec_se',
                        help='DEC column name in reference catalog')
    parser.add_argument('--match-radius', type=float, default=3.0,
                        help='Cross-match radius in arcseconds')
    parser.add_argument('--output', required=True,
                        help='Output labeled TSV file')
    parser.add_argument('--append', action='store_true',
                        help='Append to existing output file')
    args = parser.parse_args()

    # Load detections
    print(f"Loading detections: {args.detections}", file=sys.stderr)
    det_headers, detections = load_tsv(args.detections)
    print(f"  {len(detections)} sources", file=sys.stderr)

    if len(detections) == 0:
        print("ERROR: No detections loaded", file=sys.stderr)
        sys.exit(1)

    # Load reference catalog
    print(f"Loading reference: {args.reference}", file=sys.stderr)
    if args.reference.endswith('.csv'):
        references = load_csv(args.reference)
    else:
        _, references = load_tsv(args.reference)
    print(f"  {len(references)} reference sources", file=sys.stderr)

    if len(references) == 0:
        print("ERROR: No reference sources loaded", file=sys.stderr)
        sys.exit(1)

    # Cross-match
    print(f"Cross-matching (radius={args.match_radius}\")", file=sys.stderr)
    labels, match_info, n_matched = cross_match(
        detections, references, args.ra_col, args.dec_col, args.match_radius)

    print(f"  Matched: {n_matched} LSBG, "
          f"{len(detections) - n_matched} contaminants", file=sys.stderr)

    # Output
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

    mode = 'a' if args.append else 'w'
    write_header = not args.append or not os.path.exists(args.output) or \
        os.path.getsize(args.output) == 0

    with open(args.output, mode) as f:
        if write_header:
            out_headers = det_headers + ['LABEL', 'MATCH_INFO']
            f.write('\t'.join(out_headers) + '\n')

        for i, det in enumerate(detections):
            vals = [det.get(h.strip(), '') for h in det_headers]
            vals.append(str(labels[i]))
            vals.append(match_info[i])
            f.write('\t'.join(vals) + '\n')

    print(f"Saved: {args.output} ({len(detections)} rows)", file=sys.stderr)


if __name__ == '__main__':
    main()
