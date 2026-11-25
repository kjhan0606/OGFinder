#!/usr/bin/env python3
"""Plan a 180-arcmin extension of the S16A random-candidate and FDFC-grid queries for the four short lenses; local only.

With R_c = 2 Mpc and FDFC >= 0.95 (selection 2026-11-25_02) S23B_10600, 3620, 8206 and 8310 get fewer than
MINIMUM random centres inside 90 arcmin. The user chose (2026-11-25) to widen the search cone. Per lens two batches:
- random: the same thinned s16a_wide.random columns as plan cdd04062, in the annulus 90 < d <= 180 arcmin (points
  inside 90 arcmin are already held);
- fdfc: the same aggregated forced-flag cell counts as plan 019c8228, over the full 180-arcmin cone (replaces the
  90-arcmin grid of these lenses, so no cell merging across plans).
Never submits SQL.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import prepare_s16a_fdfc_grid as pf
import prepare_s16a_random_candidate_rows as pr
from prepare_s23b_s16a_quality_supplement import CAP, canonical_sha, require, sha, write_json

DATA = pr.DATA
OUTPUT = DATA/'s16a_wide_extension_plan_2026-11-25_01'
SCHEMA = 'astrafek-s16a-wide-extension-plan-v1'
SELECTION = Path('/home/kjhan/BACKUP/Astrafek/s16a_random_centre_selection_20261125_02/summary.json')
SELECTION_SHA = '2a4dba3ce4f9806a'
LENSES = ('S23B_10600', 'S23B_3620', 'S23B_8206', 'S23B_8310')
INNER_ARCMIN, OUTER_ARCMIN = pr.RADIUS_ARCMIN, 180.0
RANDOM_LIMIT = 45000  # annulus pi (180^2-90^2) x 0.3 x 1.41 = 32,300 at full coverage
FDFC_LIMIT = 110000   # pi 180^2 = 101,788 cells plus partial edge cells


def short_lenses():
    s = json.loads(SELECTION.read_text())
    require(sha(SELECTION).startswith(SELECTION_SHA), 'selection summary differs')
    short = tuple(sorted(g for g, v in s['per_lens'].items() if v['below_minimum']))
    require(short == tuple(sorted(LENSES)), 'short-lens set differs from the selection')
    return [l for l in pr.lenses() if l['lens_id'] in LENSES]


def random_query(lens):
    ra, dec = float(lens['ra_deg']), float(lens['dec_deg'])
    cx, cy, cz = pr.unit_vector(ra, dec)
    columns = ',\n  '.join([f"{pr.label(lens['lens_id'])} AS lens_id", *(f'{e} AS {a}' for a, e in pr.COLUMNS.items()), pr.MATCH_COUNT])
    return (f'SELECT\n  {columns}\nFROM s16a_wide.random AS r\n{pr.where(ra, dec, OUTER_ARCMIN/60)}'
            f'  AND r.cx*({cx:.17g}) + r.cy*({cy:.17g}) + r.cz*({cz:.17g}) < {math.cos(math.radians(INNER_ARCMIN/60)):.17g}\n'
            f'  AND r.iadjust_density < {pr.ADJUST_MAX}\nORDER BY r.object_id\nLIMIT {RANDOM_LIMIT};\n')


def fdfc_query(lens):
    ra, dec = float(lens['ra_deg']), float(lens['dec_deg'])
    require(OUTER_ARCMIN/60 < ra < 360-OUTER_ARCMIN/60 and abs(dec) < 80, 'cone crosses RA 0/360 or a pole')
    sx, sy = math.cos(math.radians(dec))*60/pf.CELL_ARCMIN, 60/pf.CELL_ARCMIN
    return (f"SELECT\n  {pr.label(lens['lens_id'])} AS lens_id,\n"
            f"  {pf.cell('f.ra', ra, sx)} AS cell_x,\n"
            f"  {pf.cell('f.dec', dec, sy)} AS cell_y,\n"
            f"  count(*) AS n_primary,\n"
            f"  sum(CASE WHEN f.wl_fulldepth_fullcolor THEN 1 ELSE 0 END) AS n_fdfc,\n"
            f"  sum(CASE WHEN f.wl_fulldepth_fullcolor IS NULL THEN 1 ELSE 0 END) AS n_fdfc_null,\n"
            f"  count(*) OVER () AS n_cells,\n"
            f"  sum(count(*)) OVER () AS n_objects\n"
            f"FROM s16a_wide.forced AS f\n"
            f"WHERE coneSearch(f.coord, {ra:.10f}, {dec:.10f}, {OUTER_ARCMIN*60:.1f})\n"
            f"  AND f.detect_is_primary\n"
            f"GROUP BY 1, 2, 3\nORDER BY 2, 3\nLIMIT {FDFC_LIMIT};\n")


def ledger():
    out = []
    for i, lens in enumerate(short_lenses()):
        out.append((f's16a_wide_random_{i:03d}.sql', 'random', lens, RANDOM_LIMIT, random_query(lens)))
        out.append((f's16a_wide_fdfc_{i:03d}.sql', 'fdfc', lens, FDFC_LIMIT, fdfc_query(lens)))
    return out


def constants():
    return dict(inner_arcmin=INNER_ARCMIN, outer_arcmin=OUTER_ARCMIN, random_limit=RANDOM_LIMIT, fdfc_limit=FDFC_LIMIT,
                lenses=list(LENSES), random_constants=pr.constants(), fdfc_constants=pf.constants())


def evidence():
    return dict(random=pr.check_schema(pr.schema_file()), forced=pf.check_schema(pf.schema_file()))


def inputs():
    return {str(p): sha(p) for p in (pr.LEDGER, pr.SCHEMA_DIR/'sources.json', pr.schema_file(), pf.SUPPLEMENT/'sources.json',
                                     pf.schema_file(), SELECTION)}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s16a_random_candidate_rows.py', 'prepare_s16a_fdfc_grid.py',
                                                 'prepare_s16a_bright_object_flag_rows.py', 'prepare_s23b_s16a_quality_supplement.py'))]


def batches(output):
    return [dict(name=n, kind=k, lens=l, row_limit=lim, sql_bytes=(output/n).stat().st_size, sha256=sha(output/n))
            for n, k, l, lim, _ in ledger()]


def validate(output):
    """Local exact-plan audit; regenerates every SQL. No network."""
    plan = json.loads((output/'plan.json').read_text())
    pr.check_gate(plan)
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['constants'] == constants() and plan['maximum_total_received_bytes'] == CAP
            and plan['random_output_columns'] == pr.OUTPUT_COLUMNS and plan['fdfc_output_columns'] == pf.OUTPUT_COLUMNS
            and plan['input_sha256'] == inputs() and plan['schema_column_check'] == evidence(), 'plan scope differs')
    for n, _, _, _, sql in ledger():
        require((output/n).read_text() == sql, 'SQL scope differs '+n)
    require(plan['batches'] == batches(output) and len(plan['batches']) == 2*len(LENSES), 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], lenses=len(LENSES),
                batches=len(plan['batches']), total_sql_bytes=sum(b['sql_bytes'] for b in plan['batches']),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    items = ledger()
    output.mkdir(parents=True)
    for n, _, _, _, sql in items:
        (output/n).write_text(sql)
    ring = math.pi*(OUTER_ARCMIN**2-INNER_ARCMIN**2)
    rand_rows = ring*pr.DENSITY_PER_ARCMIN2*pr.ADJUST_MAX*pr.DUPLICATION
    cells = math.pi*OUTER_ARCMIN**2
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_tables=['s16a_wide.random', 's16a_wide.forced'],
        purpose='Widen the random-centre search to 180 arcmin for the four lenses below MINIMUM in selection 2026-11-25_02: thinned '
                'random points in the 90-180 arcmin annulus and the FDFC cell grid over the full 180-arcmin cone.',
        constants=constants(), batches=batches(output), random_output_columns=pr.OUTPUT_COLUMNS, fdfc_output_columns=pf.OUTPUT_COLUMNS,
        random_column_expressions=pr.COLUMNS, schema_column_check=evidence(), maximum_total_received_bytes=CAP,
        estimated_rows=dict(random_per_lens=round(rand_rows), fdfc_cells_per_lens=round(cells),
                            note='full-coverage upper estimates; S23B_8206 and S23B_8310 (92.6 arcmin apart) share most of the area'),
        estimated_csv_bytes=dict(random=int(4*rand_rows*432), fdfc=int(4*cells*41),
                                 method='432 B per random row and 41 B per cell, measured on plans cdd04062 and 019c8228',
                                 hard_stop='the 128 MiB per-plan byte budget'),
        query_form='random: dec/ra box on the (dec, ra) index for the 180-arcmin cone, exact annulus from cx/cy/cz (>= cos 180 arcmin, '
                   '< cos 90 arcmin), iadjust_density thinning, ORDER BY object_id, count(*) OVER () before LIMIT. fdfc: identical to '
                   'plan 019c8228 with coneSearch radius 10800 arcsec. The random dec stripe is 6 degrees, so batches may take longer.',
        analysis_plan=['Random candidates for a lens are all held points (both plans, deduplicated by object_id) within its cone, '
                       'not only rows labelled with that lens (fixes the overlap of the 8206 and 8310 cones).',
                       'Geometry for these four lenses: 2 R <= d <= 180 arcmin - R; all other rules unchanged (R_c = 2 Mpc, '
                       'FDFC >= 0.95 on the 180-arcmin grid, tolerances, round-robin placement, TARGET 8, MINIMUM 5).',
                       'The FDFC grid is checked on the 66 pilot centres with the accepted rule (empty discs < 0.1; source density per '
                       'FDFC area consistent) and must agree cell by cell with the 90-arcmin grid where both exist.'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={b['name']: b['sha256'] for b in batches(output)},
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
