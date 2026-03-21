#!/usr/bin/env python3
"""Merge LSBG riz detection catalog with g-band forced photometry.

Applies Tanoglidis+2021 g-band surface brightness selection.
Replaces the inline python3 -c step in reproduce_tanoglidis2021_riz.sh.

Usage:
    python3 scripts/lsbg_merge_forced.py \
        --detection results/tanoglidis2021_riz/lsbg_riz_detected.tsv \
        --forced results/tanoglidis2021_riz/lsbg_forced_g.tsv \
        --mu-eff-min 24.2 --mu-eff-max 28.5 \
        --r-eff-min 2.5 --r-eff-max 20.0 \
        --output results/tanoglidis2021_riz/lsbg_detected.tsv
"""
import argparse
import math
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Merge detection + forced photometry catalogs with SB selection")
    parser.add_argument("--detection", required=True,
                        help="Detection catalog TSV (from sep_detect_riz.py)")
    parser.add_argument("--forced", required=True,
                        help="Forced photometry TSV (from ds9_lsbg.py --mode forced)")
    parser.add_argument("--mu-eff-min", type=float, default=24.2,
                        help="Min g-band mu_eff (mag/arcsec^2)")
    parser.add_argument("--mu-eff-max", type=float, default=28.5,
                        help="Max g-band mu_eff (mag/arcsec^2)")
    parser.add_argument("--r-eff-min", type=float, default=2.5,
                        help="Min effective radius (arcsec)")
    parser.add_argument("--r-eff-max", type=float, default=20.0,
                        help="Max effective radius (arcsec)")
    parser.add_argument("--output", default=None,
                        help="Output TSV file (default: stdout)")
    args = parser.parse_args()

    # Load detection catalog
    with open(args.detection) as f:
        det_lines = f.read().strip().split('\n')
    det_header = det_lines[0].split('\t')
    det_rows = [line.split('\t') for line in det_lines[1:]]

    # Load forced photometry
    with open(args.forced) as f:
        forced_lines = f.read().strip().split('\n')
    forced_header = forced_lines[0].split('\t')
    forced_rows = [line.split('\t') for line in forced_lines[1:]]

    # Build merged header with MU_EFF_G column
    merge_header = det_header + [c for c in forced_header if c != 'NUMBER'] + ['MU_EFF_G']

    out = open(args.output, 'w') if args.output else sys.stdout
    out.write('\t'.join(merge_header) + '\n')

    fr_idx = det_header.index('FLUX_RADIUS_ARCSEC') if 'FLUX_RADIUS_ARCSEC' in det_header else -1
    r_eff_idx = det_header.index('R_EFF_ARCSEC') if 'R_EFF_ARCSEC' in det_header else -1

    n_pass = 0
    for i, det_row in enumerate(det_rows):
        if i >= len(forced_rows):
            break

        forced_dict = dict(zip(forced_header, forced_rows[i]))
        mag_g = float(forced_dict.get('MAG_g', '99'))

        r_eff_moment = float(det_row[r_eff_idx]) if r_eff_idx >= 0 else 5.0
        r_eff_mu = float(det_row[fr_idx]) if fr_idx >= 0 else r_eff_moment
        r_eff_mu = max(r_eff_mu, r_eff_moment * 0.5)
        r_eff_mu = min(r_eff_mu, 20.0)

        if mag_g < 90 and r_eff_mu > 0:
            mu_eff_g = mag_g + 2.5 * math.log10(2 * math.pi * r_eff_mu ** 2)
        else:
            mu_eff_g = 99.0

        # Apply g-band selection
        if mu_eff_g < args.mu_eff_min or mu_eff_g > args.mu_eff_max:
            continue
        if r_eff_moment < args.r_eff_min or r_eff_moment > args.r_eff_max:
            continue

        forced_vals = [forced_dict[c] for c in forced_header if c != 'NUMBER']
        merged = det_row + forced_vals + [f'{mu_eff_g:.3f}']
        out.write('\t'.join(merged) + '\n')
        n_pass += 1

    if args.output:
        out.close()

    print(f'  g-band selection: {n_pass}/{len(det_rows)} passed '
          f'(mu_eff_g {args.mu_eff_min}-{args.mu_eff_max}, '
          f'r_eff {args.r_eff_min}-{args.r_eff_max}")',
          file=sys.stderr)


if __name__ == '__main__':
    main()
