#!/usr/bin/env python3
"""Cross-match a source catalog with VizieR catalogs."""

import sys
import os
import argparse
import numpy as np
from urllib.request import urlopen, Request
from urllib.parse import urlencode
import io


# Supported VizieR catalogs
VIZIER_CATALOGS = {
    'GAIA_DR3':     {'table': 'I/355/gaiadr3',  'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
                     'cols': 'RA_ICRS,DE_ICRS,Gmag,BPmag,RPmag,Plx,pmRA,pmDE',
                     'label': 'GAIA_DR3'},
    '2MASS':        {'table': 'II/246/out',      'ra': 'RAJ2000', 'dec': 'DEJ2000',
                     'cols': 'RAJ2000,DEJ2000,Jmag,Hmag,Kmag',
                     'label': '2MASS'},
    'PanSTARRS_DR1':{'table': 'II/349/ps1',      'ra': 'RAJ2000', 'dec': 'DEJ2000',
                     'cols': 'RAJ2000,DEJ2000,gmag,rmag,imag,zmag,ymag',
                     'label': 'PS1'},
    'SDSS_DR17':    {'table': 'V/154/sdss16',    'ra': 'RA_ICRS', 'dec': 'DE_ICRS',
                     'cols': 'RA_ICRS,DE_ICRS,umag,gmag,rmag,imag,zmag',
                     'label': 'SDSS'},
    'ALLWISE':      {'table': 'II/328/allwise',  'ra': 'RAJ2000', 'dec': 'DEJ2000',
                     'cols': 'RAJ2000,DEJ2000,W1mag,W2mag,W3mag,W4mag',
                     'label': 'WISE'},
}


def query_vizier_tap(table, ra_col, dec_col, cols, ra_c, dec_c, radius_deg):
    """Query VizieR TAP service."""
    adql = (f'SELECT {cols} FROM "{table}" '
            f'WHERE 1=CONTAINS(POINT(\'ICRS\',{ra_col},{dec_col}), '
            f'CIRCLE(\'ICRS\',{ra_c:.6f},{dec_c:.6f},{radius_deg:.6f}))')

    url = 'http://tapvizier.cds.unistra.fr/TAPVizieR/tap/sync'
    params = {
        'REQUEST': 'doQuery',
        'LANG': 'ADQL',
        'FORMAT': 'tsv',
        'QUERY': adql,
    }
    full_url = url + '?' + urlencode(params)
    print(f"Querying VizieR: {table}...", file=sys.stderr)

    req = Request(full_url, headers={'User-Agent': 'DS9-SExtractor/1.0'})
    try:
        with urlopen(req, timeout=60) as resp:
            data = resp.read().decode('utf-8', errors='replace')
    except Exception as e:
        print(f"ERROR: VizieR query failed: {e}", file=sys.stderr)
        return None, None, None

    lines = [l.strip() for l in data.split('\n') if l.strip() and not l.startswith('#')]
    if len(lines) < 2:
        return None, None, None

    header = lines[0]
    # Skip type row if present (VizieR TSV has header, types, dashes, data)
    data_lines = []
    for l in lines[1:]:
        if l.startswith('---') or l.startswith('char') or l.startswith('float') or l.startswith('double') or l.startswith('int') or l.startswith('long'):
            continue
        data_lines.append(l)

    # Parse RA/Dec from VizieR results
    hdr_fields = [h.strip() for h in header.split('\t')]
    ra_idx = None
    dec_idx = None
    for i, h in enumerate(hdr_fields):
        if h == ra_col:
            ra_idx = i
        elif h == dec_col:
            dec_idx = i

    if ra_idx is None or dec_idx is None:
        return None, None, None

    viz_ra = []
    viz_dec = []
    viz_rows = []
    for line in data_lines:
        fields = line.split('\t')
        try:
            r = float(fields[ra_idx].strip())
            d = float(fields[dec_idx].strip())
            viz_ra.append(r)
            viz_dec.append(d)
            viz_rows.append(fields)
        except (ValueError, IndexError):
            continue

    return np.array(viz_ra), np.array(viz_dec), viz_rows


def main():
    parser = argparse.ArgumentParser(description='Cross-match with VizieR catalogs')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--vizier-cat', default='II/349',
                        help='VizieR catalog ID or name (e.g., GAIA_DR3, 2MASS, PanSTARRS_DR1, SDSS_DR17, ALLWISE)')
    parser.add_argument('--radius', type=float, default=2.0,
                        help='Match radius in arcseconds (default: 2.0)')
    args = parser.parse_args()

    # Resolve catalog
    cat_info = None
    for name, info in VIZIER_CATALOGS.items():
        if args.vizier_cat == name or args.vizier_cat == info['table']:
            cat_info = info
            break
    if cat_info is None:
        # Try as raw table ID
        cat_info = {'table': args.vizier_cat, 'ra': 'RAJ2000', 'dec': 'DEJ2000',
                     'cols': '*', 'label': args.vizier_cat.replace('/', '_')}

    # Parse input catalog
    with open(args.catalog, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    headers = [h.strip() for h in lines[0].split('\t')]
    col_map = {h: i for i, h in enumerate(headers)}

    if 'ALPHA_J2000' not in col_map or 'DELTA_J2000' not in col_map:
        print("ERROR: Catalog must have ALPHA_J2000 and DELTA_J2000 columns", file=sys.stderr)
        sys.exit(1)

    src_ra = []
    src_dec = []
    src_numbers = []
    for line in lines[1:]:
        fields = line.split('\t')
        try:
            ra = float(fields[col_map['ALPHA_J2000']].strip())
            dec = float(fields[col_map['DELTA_J2000']].strip())
            num = int(fields[col_map['NUMBER']].strip())
            src_ra.append(ra)
            src_dec.append(dec)
            src_numbers.append(num)
        except (ValueError, IndexError):
            continue

    src_ra = np.array(src_ra)
    src_dec = np.array(src_dec)

    if len(src_ra) == 0:
        print("ERROR: No valid sources with RA/Dec", file=sys.stderr)
        sys.exit(1)

    # Field center and radius
    ra_c = np.mean(src_ra)
    dec_c = np.mean(src_dec)
    # Field radius (add margin)
    cos_dec = np.cos(np.radians(dec_c))
    dra = (src_ra - ra_c) * cos_dec
    ddec = src_dec - dec_c
    field_radius = np.max(np.sqrt(dra**2 + ddec**2)) + args.radius / 3600.0

    # Query VizieR
    viz_ra, viz_dec, viz_rows = query_vizier_tap(
        cat_info['table'], cat_info['ra'], cat_info['dec'],
        cat_info['cols'], ra_c, dec_c, field_radius)

    if viz_ra is None or len(viz_ra) == 0:
        print("ERROR: No matches from VizieR", file=sys.stderr)
        sys.exit(1)

    print(f"VizieR returned {len(viz_ra)} sources", file=sys.stderr)

    # Cross-match using angular separation
    match_radius_deg = args.radius / 3600.0
    n_matched = 0

    # Output header
    print(f"NUMBER\tMATCH_DIST\tMATCH_ID")

    for i in range(len(src_ra)):
        # Angular separation
        cos_d = np.cos(np.radians(src_dec[i]))
        dra = (viz_ra - src_ra[i]) * cos_d
        ddec = viz_dec - src_dec[i]
        dist = np.sqrt(dra**2 + ddec**2)
        best_idx = np.argmin(dist)
        best_dist = dist[best_idx]

        if best_dist <= match_radius_deg:
            dist_arcsec = best_dist * 3600.0
            print(f"{src_numbers[i]}\t{dist_arcsec:.3f}\t{best_idx+1}")
            n_matched += 1
        else:
            print(f"{src_numbers[i]}\tnan\t0")

    print(f"Done: {n_matched}/{len(src_ra)} sources matched", file=sys.stderr)


if __name__ == '__main__':
    main()
