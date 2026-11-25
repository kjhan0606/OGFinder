#!/usr/bin/env python3
"""Plan an aggregated S16A full-depth full-colour (FDFC) grid around the 11 pilot lenses; local files only.

No FDFC map was released for S16A (Mandelbaum et al. 2018, sec. 2.2: visits per NSIDE=1024 HEALPix pixel
>= (4,4,4,6,6) in grizy, a fitted i_lim > 25.6, PSF-failure and disconnected-region removal; the HSC release
pages list no window file). Point-level random countinputs do not reproduce it: in the 2026-11-25 inner random
points the (4,4,4,6,6) cut fails 0-100% of points around centres whose S16A sources are all present. The flag
itself, s16a_wide.forced.wl_fulldepth_fullcolor, is the catalogue's FDFC decision per object, so counting
primary objects with the flag true/false per 1-arcmin cell inside the 90-arcmin candidate cones maps the
footprint empirically. Aggregated counts only: no object rows, positions or photometry. Never submits SQL.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from prepare_s16a_bright_object_flag_rows import SUPPLEMENT
from prepare_s16a_random_candidate_rows import LEDGER, N_LENSES, RADIUS_ARCMIN, check_gate, label, lenses
from prepare_s23b_s16a_quality_supplement import CAP, canonical_sha, require, sha, write_json

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_fdfc_grid_plan_2026-11-25_01'
SCHEMA = 'astrafek-s16a-fdfc-grid-plan-v1'
SCHEMA_URL = 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.forced.json'
CELL_ARCMIN = 1.0
LIMIT_PER_LENS = 40000  # pi 90^2 = 25447 cells plus partial edge cells
OUTPUT_COLUMNS = ['lens_id', 'cell_x', 'cell_y', 'n_primary', 'n_fdfc', 'n_fdfc_null', 'n_cells', 'n_objects']
REQUIRED = {'ra': 'double precision', 'dec': 'double precision', 'coord': 'earth', 'detect_is_primary': 'bool',
            'wl_fulldepth_fullcolor': 'bool'}
MIN_CELL_OBJECTS = 20
FDFC_FRACTION = 0.5


def schema_file():
    manifest = json.loads((SUPPLEMENT/'sources.json').read_text())
    entry, = [d for d in manifest['documents'] if d['file'] == 'schema_forced.json']
    file = SUPPLEMENT/entry['file']
    require(entry['url'] == SCHEMA_URL and sha(file) == entry['sha256'], 'forced schema provenance/hash differs')
    return file


def check_schema(file):
    doc = json.loads(file.read_text())
    require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', 'forced'), 'schema identity differs')
    documented = {r[1]: dict(zip(doc['headers'], r)) for r in doc['columns']}
    evidence = []
    for name, kind in REQUIRED.items():
        info = documented.get(name)
        require(info is not None and info['type'] == kind, f'forced column {name} missing or not {kind}')
        evidence.append(dict(column='f.'+name, type=info['type'], description=info['description']))
    return evidence


def cell(expr, origin, scale):
    return f'floor(({expr} - ({origin:.10f}))*{scale:.17g})::integer'


def query(lens):
    ra, dec = float(lens['ra_deg']), float(lens['dec_deg'])
    require(RADIUS_ARCMIN/60 < ra < 360-RADIUS_ARCMIN/60 and abs(dec) < 80, 'cone crosses RA 0/360 or a pole')
    sx, sy = math.cos(math.radians(dec))*60/CELL_ARCMIN, 60/CELL_ARCMIN
    return (f"SELECT\n  {label(lens['lens_id'])} AS lens_id,\n"
            f"  {cell('f.ra', ra, sx)} AS cell_x,\n"
            f"  {cell('f.dec', dec, sy)} AS cell_y,\n"
            f"  count(*) AS n_primary,\n"
            f"  sum(CASE WHEN f.wl_fulldepth_fullcolor THEN 1 ELSE 0 END) AS n_fdfc,\n"
            f"  sum(CASE WHEN f.wl_fulldepth_fullcolor IS NULL THEN 1 ELSE 0 END) AS n_fdfc_null,\n"
            f"  count(*) OVER () AS n_cells,\n"
            f"  sum(count(*)) OVER () AS n_objects\n"
            f"FROM s16a_wide.forced AS f\n"
            f"WHERE coneSearch(f.coord, {ra:.10f}, {dec:.10f}, {RADIUS_ARCMIN*60:.1f})\n"
            f"  AND f.detect_is_primary\n"
            f"GROUP BY 1, 2, 3\nORDER BY 2, 3\nLIMIT {LIMIT_PER_LENS};\n")


def ledger():
    return [(f's16a_fdfc_grid_{i:03d}.sql', lens) for i, lens in enumerate(lenses())]


def constants():
    return dict(radius_arcmin=RADIUS_ARCMIN, cell_arcmin=CELL_ARCMIN, cell_definition='floor((ra-ra0) cos(dec0) 60), floor((dec-dec0) 60)',
                row_limit_per_lens=LIMIT_PER_LENS, min_cell_objects=MIN_CELL_OBJECTS, fdfc_fraction=FDFC_FRACTION)


def inputs(file):
    return {str(p): sha(p) for p in (LEDGER, SUPPLEMENT/'sources.json', file)}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s16a_random_candidate_rows.py', 'prepare_s16a_bright_object_flag_rows.py',
                                                 'prepare_s23b_s16a_quality_supplement.py'))]


def batches(output):
    return [dict(name=name, lens=lens, row_limit=LIMIT_PER_LENS, sql_bytes=(output/name).stat().st_size, sha256=sha(output/name))
            for name, lens in ledger()]


def validate(output):
    """Local exact-plan audit; regenerates every SQL from the pinned ledger. No network."""
    plan = json.loads((output/'plan.json').read_text())
    check_gate(plan)
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['output_columns'] == OUTPUT_COLUMNS and plan['constants'] == constants()
            and plan['cas_release_version'] == 'pdr1' and plan['source_table'] == 's16a_wide.forced'
            and plan['maximum_total_received_bytes'] == CAP, 'plan scope differs')
    file = schema_file()
    require(plan['input_sha256'] == inputs(file) and plan['schema_column_check'] == check_schema(file), 'schema evidence differs')
    for name, lens in ledger():
        require((output/name).read_text() == query(lens), 'SQL scope differs '+name)
    require(plan['batches'] == batches(output) and len(plan['batches']) == N_LENSES, 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], lenses=N_LENSES,
                columns=len(OUTPUT_COLUMNS), batches=N_LENSES, total_sql_bytes=sum(b['sql_bytes'] for b in plan['batches']),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    file = schema_file()
    evidence = check_schema(file)
    output.mkdir(parents=True)
    for name, lens in ledger():
        (output/name).write_text(query(lens))
    items = batches(output)
    cells = math.pi*RADIUS_ARCMIN**2/CELL_ARCMIN**2
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_table='s16a_wide.forced',
        purpose='Per-cell counts of primary S16A forced objects and of those with wl_fulldepth_fullcolor true, on a 1-arcmin grid '
                'inside the 90-arcmin cones of the 11 pilot lenses, to map the FDFC footprint for the random-centre selector.',
        constants=constants(), n_lenses=N_LENSES, n_query_batches=N_LENSES, batches=items,
        output_columns=OUTPUT_COLUMNS, schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_rows=dict(total=round(cells*N_LENSES), note='about 25,400 cells per lens if the whole cone has objects; fewer outside the survey'),
        estimated_csv_bytes=dict(low=int(cells*N_LENSES*40), high=int(cells*N_LENSES*70), hard_stop='the 128 MiB per-plan byte budget'),
        query_form='coneSearch(f.coord, ra, dec, 5400 arcsec) on the documented earth column, detect_is_primary, GROUP BY the cell '
                   'indices, count(*) OVER () AS n_cells and sum(count(*)) OVER () AS n_objects (evaluated before LIMIT; the receiver '
                   'requires rows == n_cells and sum n_primary == n_objects), ORDER BY cell, LIMIT 40000 (a batch at the limit is refused).',
        excluded=['object rows, IDs, positions, photometry, shapes or photo-z', 'non-primary objects', 'other releases'],
        limitations=['Cells are 1 arcmin in a local equirectangular frame (<0.3% distortion at 90 arcmin); FDFC HEALPix pixels are '
                     'about 3.4 arcmin, so boundary cells are mixed.',
                     'Cells with few primary objects (bright-star halos, survey edge) have noisy fractions.'],
        analysis_plan=[f'Cell FDFC fraction = n_fdfc/n_primary; a cell is FDFC if n_primary >= {MIN_CELL_OBJECTS} and the fraction >= '
                       f'{FDFC_FRACTION}; cells with fewer objects or absent cells are not FDFC.',
                       'Disc FDFC area fraction = share of a disc\'s cell-centre lattice (0.25 arcmin sampling) in FDFC cells.',
                       'Validation (pre-registered): each of the 23 empty pilot randoms should have disc FDFC fraction < 0.1 and each '
                       'of the 43 non-empty centres > 0.5; any mismatch is reported and stops use of the grid until reviewed.',
                       'If validated, the selector requires FDFC fraction >= 0.95 over each candidate 3-Mpc disc.',
                       'No selection change, correction or science promotion by this plan.'],
        before_submission=['User approval of this exact plan hash and of the receiver code hashes recorded with it.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(file), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={b['name']: b['sha256'] for b in items},
    )
    plan['query_plan_sha256'] = canonical_sha(plan)
    write_json(output/'plan.json', plan)
    (output/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  plan.json\n')
    return validate(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--verify', action='store_true', help='Verify an existing local plan; no network')
    args = parser.parse_args()
    print(json.dumps(validate(args.output.resolve()) if args.verify else build(args.output), indent=2))
