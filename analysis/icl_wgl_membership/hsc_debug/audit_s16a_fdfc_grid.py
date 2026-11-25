#!/usr/bin/env python3
"""Pre-registered validation of the S16A FDFC grid (plan 019c8228) on the 66 pilot centres; local only.

Rule fixed in the plan before download: each of the 23 empty pilot randoms (zero S16A sources) must have disc
FDFC fraction < 0.1 and each of the 43 non-empty centres > 0.5. Any mismatch stops use of the grid.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from prepare_s16a_random_candidate_rows import centres
from select_s16a_random_centres import fdfc_fraction, load_fdfc, sha

STATUS = Path('/home/kjhan/BACKUP/Astrafek/s16a_random_point_areas_20261125_01/centre_point_status.tsv')
EMPTY_MAX, FULL_MIN = 0.1, 0.5


def run(output):
    fdfc = load_fdfc()
    if fdfc is None:
        raise SystemExit('FDFC grid receipt absent')
    if output.exists():
        raise FileExistsError(output)
    with STATUS.open() as stream:
        n_src = {r['center_id']: int(r['n_sources']) for r in csv.DictReader(stream, delimiter='\t')}
    lenses = {r['matched_lens_id']: r for r in centres() if r['role'] == 'lens'}
    rows = []
    for r in centres():
        g = r['matched_lens_id']
        f = fdfc_fraction(lenses[g], float(r['ra_deg']), float(r['dec_deg']), float(r['radius_arcsec'])/3600, fdfc['grids'].get(g, {}))
        empty = n_src[r['center_id']] == 0
        rows.append(dict(center_id=r['center_id'], role=r['role'], lens_id=g, n_sources=n_src[r['center_id']],
                         fdfc_fraction=f, expected='empty' if empty else 'nonempty',
                         passes=f < EMPTY_MAX if empty else f > FULL_MIN))
    output.mkdir(parents=True)
    with (output/'centre_fdfc.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    fails = [r for r in rows if not r['passes']]
    summary = dict(schema='astrafek-s16a-fdfc-grid-validation-v1', grid_summary_sha256=fdfc['summary_sha256'],
                   rule=dict(empty_max=EMPTY_MAX, nonempty_min=FULL_MIN), n_centres=len(rows),
                   n_empty=sum(r['expected'] == 'empty' for r in rows), n_fail=len(fails),
                   failures=[{k: r[k] for k in ('center_id', 'lens_id', 'n_sources', 'fdfc_fraction')} for r in fails],
                   grid_validated=not fails, status_sha256=sha(STATUS), code_sha256=sha(Path(__file__)),
                   centre_fdfc_sha256=sha(output/'centre_fdfc.tsv'))
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_fdfc_grid_validation_20261125_01'))
    s = run(parser.parse_args().output)
    print(json.dumps({k: s[k] for k in ('n_centres', 'n_empty', 'n_fail', 'failures', 'grid_validated')}, indent=1))
