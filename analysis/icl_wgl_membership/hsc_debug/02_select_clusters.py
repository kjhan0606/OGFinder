#!/usr/bin/env python3
"""Select four CAMIRA clusters for the ICL-stack debug.

The catalog columns are id, name, RA, Dec, z_cl, N_mem, and the BCG
spectroscopic redshift (-1 when absent). There is no mass column.

The redshift slice is the 0.05-wide bin, on exact edges from 0.10 to 1.20,
that contains the most rows with a spectroscopic BCG redshift inside the
fractional window. The fractional window is |z_spec - z_cl| / (1 + z_cl) < 0.01,
the scatter Oguri et al. 2018 (arXiv:1701.00818) quote for most of their range.
Richness keeps rows at or above the 90th percentile of that slice
(numpy.percentile, linear). Four rows are then taken in descending richness,
with each pair more than 0.3 degrees apart.
"""
import json
import math
import os
import sys

import numpy as np

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')
CATALOG = os.path.join(ROOT, 'catalog', 'camira_s20a_wide_sm.tbl')
FRAC_MAX = 0.01
SEP_MIN_DEG = 0.3
N_KEEP = 4


def read_catalog(path):
    rows = []
    with open(path, encoding='utf-8') as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 7:
                raise SystemExit('camira row does not have 7 columns: %s' % line.strip())
            rows.append((
                int(parts[0]), parts[1], float(parts[2]), float(parts[3]),
                float(parts[4]), float(parts[5]), float(parts[6])))
    return rows


def angular_sep_deg(a, b):
    ra1, dec1 = math.radians(a[2]), math.radians(a[3])
    ra2, dec2 = math.radians(b[2]), math.radians(b[3])
    cos = (math.sin(dec1) * math.sin(dec2)
           + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos))))


def slice_edges():
    return [round(0.10 + 0.05 * i, 2) for i in range(23)]


def main():
    rows = read_catalog(CATALOG)
    spec = [row for row in rows if row[6] > 0]
    frac = np.array([abs(row[6] - row[4]) / (1.0 + row[4]) for row in spec])
    zcl = np.array([row[4] for row in spec])
    counts = []
    for lo, hi in zip(slice_edges()[:-1], slice_edges()[1:]):
        n = int(np.sum((zcl >= lo) & (zcl < hi) & (frac < FRAC_MAX)))
        counts.append((n, lo, hi))
    n_win, z_lo, z_hi = max(counts)
    window = [
        row for row in spec
        if z_lo <= row[4] < z_hi and abs(row[6] - row[4]) / (1.0 + row[4]) < FRAC_MAX
    ]
    nmem = np.array([row[5] for row in window])
    threshold = float(np.percentile(nmem, 90))
    rich = sorted((row for row in window if row[5] >= threshold), key=lambda row: -row[5])
    chosen = []
    for row in rich:
        if all(angular_sep_deg(row, other) > SEP_MIN_DEG for other in chosen):
            chosen.append(row)
        if len(chosen) == N_KEEP:
            break
    if len(chosen) < N_KEEP:
        raise SystemExit('fewer than %d separated clusters' % N_KEEP)
    seps = [
        angular_sep_deg(chosen[i], chosen[j])
        for i in range(N_KEEP) for j in range(i + 1, N_KEEP)
    ]
    os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
    table = os.path.join(ROOT, 'logs', 'selected_clusters.tsv')
    with open(table, 'w', encoding='utf-8') as handle:
        handle.write('id\tname\tra_deg\tdec_deg\tz_cl\tN_mem\tz_spec_bcg\tfrac_dz\n')
        for row in chosen:
            frac_row = abs(row[6] - row[4]) / (1.0 + row[4])
            handle.write('%d\t%s\t%.6f\t%.6f\t%.4f\t%.3f\t%.5f\t%.6f\n' % (
                row[0], row[1], row[2], row[3], row[4], row[5], row[6], frac_row))
    summary = {
        'n_rows': len(rows),
        'n_spec': len(spec),
        'frac_p50': float(np.percentile(frac, 50)),
        'frac_p68': float(np.percentile(frac, 68)),
        'z_lo': z_lo,
        'z_hi': z_hi,
        'n_window': len(window),
        'nmem_p90': threshold,
        'n_rich': len(rich),
        'min_sep_deg': min(seps),
        'ids': [row[0] for row in chosen],
        'numpy': np.__version__,
        'science_measurement': False,
    }
    with open(os.path.join(ROOT, 'logs', 'selection_summary.json'), 'w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=1)
        handle.write('\n')
    print('selected %s' % ' '.join(str(i) for i in summary['ids']))
    print('window %.2f %.2f n %d Nmem_p90 %.4f' % (z_lo, z_hi, len(window), threshold))


if __name__ == '__main__':
    sys.exit(main())
