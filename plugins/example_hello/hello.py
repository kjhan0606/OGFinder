#!/usr/bin/env python3
"""Example OGFinder plugin driver: read the catalog TSV, add a HELLO column, print the result TSV on stdout.

This is the whole contract of a CLI step: the catalog is passed as --catalog FILE (tab separated, header row);
the program prints a TSV with the key column NUMBER plus the NEW columns (header row + one row per catalog
row) to stdout and exits with 0.  The GUI adds them to the table; the session recorder stores the argv so that the exported
pipeline script runs the very same command on new data."""
import argparse
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--greeting', default='hello')
    ap.add_argument('--min-mag', type=float, default=99.0)
    ap.add_argument('--shout', action='store_true')
    a = ap.parse_args()
    with open(a.catalog) as f:
        lines = [l.rstrip('\n') for l in f if l.strip()]
    head = lines[0].split('\t')
    im = head.index('MAG_AUTO') if 'MAG_AUTO' in head else -1
    text = a.greeting.upper() if a.shout else a.greeting
    inum = head.index('NUMBER') if 'NUMBER' in head else -1
    out = ['NUMBER\tHELLO']
    for k, l in enumerate(lines[1:], 1):
        f = l.split('\t')
        num = f[inum].strip() if inum >= 0 else str(k)
        ok = True
        if im >= 0:
            try:
                ok = float(f[im]) <= a.min_mag
            except ValueError:
                ok = False
        out.append('%s\t%s' % (num, text if ok else ''))
    sys.stdout.write('\n'.join(out) + '\n')


if __name__ == '__main__':
    main()
