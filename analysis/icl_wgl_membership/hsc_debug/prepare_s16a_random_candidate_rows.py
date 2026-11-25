#!/usr/bin/env python3
"""Plan a thinned S16A random-point query around the 11 pilot lenses; local files only.

The boost audit (2026-11-25) found that the 55 random centres were drawn from a PDR3 pool within
30 arcmin with no control-control separation, so same-lens randoms overlap and 23 of 55 fall
outside the S16A source footprint. s16a_wide.random carries the S16A (Sirius) bright-object,
edge and countinputs flags at 100 points per arcmin^2. A thinned sample (iadjust_density below a
fixed threshold) within 90 arcmin of each lens gives S16A-processed candidate centres for a
non-overlapping random design and S16A seeing/depth/mask ranking over 3-Mpc discs. No flag
filtering, so masked points stay in. Revision 2 follows the independent plan review. Never submits SQL.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from prepare_s23b_s16a_quality_supplement import CAP, canonical_sha, require, sha, write_json

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_random_candidate_rows_plan_2026-11-25_02'
SCHEMA = 'astrafek-s16a-random-candidate-rows-plan-v2'
SCHEMA_DIR = DATA/'s16a_random_schema_2026-11-25_01'
SCHEMA_URL = 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.random.json'
SCHEMA_SHA = '9676a159de2374578a858735451e95051486ca29f246e7253b33adc927a03797'
LEDGER = DATA/'wl_s23b_s16a_shape_qa_pilot_plan_2026-11-17_02'/'center_ledger.tsv'
LEDGER_SHA = '46d9f286c06720a2309ae720f6a50582adacce4a9f379540477055755f819f23'
RADIUS_ARCMIN = 90.0
ADJUST_MAX = 0.003
DENSITY_PER_ARCMIN2 = 100.0
DUPLICATION = 1.41  # patch borders (4200/4000)^2 times ring-tract overlap, from the plan review; not measured
LIMIT_PER_LENS = 27000  # 11 x 27000 x ~442 B stays under the 128 MiB per-plan cap
MARGIN_DEG = 1e-6
BANDS = 'grizy'
N_LENSES = 11


def random_columns(*, thinned):
    columns = {'object_id': 'r.object_id', 'ra': 'r.ra', 'dec': 'r.dec', 'skymap_id': 'r.skymap_id'}
    if thinned:
        columns.update({'tract': 'r.tract', 'patch': 'r.patch', 'iadjust_density': 'r.iadjust_density'})
    for b in BANDS:
        columns[f'{b}countinputs'] = f'r.{b}countinputs'
    for c in ('detect_is_patch_inner', 'detect_is_tract_inner') + (('detect_is_primary',) if thinned else ()):
        columns['i'+c] = 'r.i'+c
    for flag in ('flags_pixel_edge', 'flags_pixel_bright_object_center', 'flags_pixel_bright_object_any'):
        for b in BANDS:
            columns[f'{b}{flag}'] = f'r.{b}{flag}'
    for c in ('flags_pixel_offimage', 'flags_pixel_bad', 'flags_pixel_saturated_center', 'flags_pixel_interpolated_center',
              'flags_pixel_suspect_center', 'flags_pixel_cr_center'):
        columns['i'+c] = 'r.i'+c
    if thinned:
        for c in ('shape_sdss_psf_xx', 'shape_sdss_psf_yy', 'shape_sdss_psf_xy', 'pix_variance'):
            columns['i'+c] = 'r.i'+c
        for b in BANDS:
            columns[f'{b}sky_std'] = f'r.{b}sky_std'
    return columns


COLUMNS = random_columns(thinned=True)
# count(*) OVER () is evaluated before LIMIT: every row carries the number of matching points, so a short
# (cleanly framed) response or a LIMIT cut is detected by rows != n_matched.
OUTPUT_COLUMNS = ['lens_id', *COLUMNS, 'n_matched']
MATCH_COUNT = 'count(*) OVER () AS n_matched'


def kinds(columns):
    integers = tuple(c for c in columns if c in ('object_id', 'skymap_id', 'tract', 'patch') or c.endswith('countinputs'))
    booleans = tuple(c for c in columns if 'detect_is' in c or 'flags_pixel' in c)
    return integers, booleans


INTEGERS, BOOLEANS = kinds(COLUMNS)
FLOATS = tuple(c for c in COLUMNS if c not in INTEGERS and c not in BOOLEANS)


def schema_type(alias, columns):
    integers, booleans = kinds(columns)
    if alias == 'object_id':
        return 'bigint'
    if alias in integers:
        return 'integer'
    if alias in booleans:
        return 'bool'
    return 'double precision'


def unit_vector(ra, dec):
    a, d = math.radians(ra), math.radians(dec)
    return math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d)


def box(ra, dec, radius_deg=RADIUS_ARCMIN/60):
    """Dec and RA ranges that contain the cone; refuses a pole or an RA wrap (none in this sample)."""
    require(abs(dec)+radius_deg < 89.0, 'cone too close to a pole')
    half_ra = math.degrees(math.asin(math.sin(math.radians(radius_deg))/math.cos(math.radians(dec))))+MARGIN_DEG
    lo, hi = ra-half_ra, ra+half_ra
    require(0.0 <= lo and hi < 360.0, 'cone crosses RA 0/360; not supported by this plan')
    return (dec-radius_deg-MARGIN_DEG, dec+radius_deg+MARGIN_DEG), (lo, hi)


def where(ra, dec, radius_deg):
    (dlo, dhi), (alo, ahi) = box(ra, dec, radius_deg)
    cx, cy, cz = unit_vector(ra, dec)
    return (f'WHERE r.dec BETWEEN {dlo:.10f} AND {dhi:.10f}\n'
            f'  AND r.ra BETWEEN {alo:.10f} AND {ahi:.10f}\n'
            f'  AND r.cx*({cx:.17g}) + r.cy*({cy:.17g}) + r.cz*({cz:.17g}) >= {math.cos(math.radians(radius_deg)):.17g}\n')


def label(text):
    require(text.replace('_', '').replace('-', '').isalnum(), 'unexpected identifier '+text)
    return f"'{text}'::text"


def query(lens):
    ra, dec = float(lens['ra_deg']), float(lens['dec_deg'])
    columns = ',\n  '.join([f"{label(lens['lens_id'])} AS lens_id", *(f'{expr} AS {alias}' for alias, expr in COLUMNS.items()), MATCH_COUNT])
    return (f'SELECT\n  {columns}\nFROM s16a_wide.random AS r\n{where(ra, dec, RADIUS_ARCMIN/60)}'
            f'  AND r.iadjust_density < {ADJUST_MAX}\nORDER BY r.object_id\nLIMIT {LIMIT_PER_LENS};\n')


def centres():
    require(sha(LEDGER) == LEDGER_SHA, 'pilot centre ledger differs')
    with LEDGER.open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(len(rows) == 66 and sum(r['role'] == 'lens' for r in rows) == N_LENSES, 'centre ledger differs')
    return rows


def lenses():
    rows = [r for r in centres() if r['role'] == 'lens']
    require(len({r['matched_lens_id'] for r in rows}) == N_LENSES, 'lens set differs')
    return sorted(({'lens_id': r['matched_lens_id'], 'center_id': r['center_id'], 'ra_deg': r['ra_deg'], 'dec_deg': r['dec_deg'],
                    'z_lens': r['z_lens'], 'radius_arcsec': r['radius_arcsec']} for r in rows), key=lambda r: r['lens_id'])


def schema_file():
    manifest = json.loads((SCHEMA_DIR/'sources.json').read_text())
    entry, = [d for d in manifest['documents'] if d['file'] == 'schema_random.json']
    file = SCHEMA_DIR/entry['file']
    require(entry['url'] == SCHEMA_URL and entry['sha256'] == SCHEMA_SHA == sha(file) and manifest['credentials_used'] is False,
            'random schema provenance/hash differs')
    return file


def check_schema(file, columns=COLUMNS):
    doc = json.loads(file.read_text())
    require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', 'random'), 'schema identity differs')
    require('CREATE INDEX ON ... USING btree ("dec", ra)' in doc['indices'], 'documented (dec, ra) index missing')
    documented = {}
    for row in doc['columns']:
        c = dict(zip(doc['headers'], row))
        for band in BANDS:
            documented[c['name'].replace('(g|r|i|z|y)', band)] = c
    for name in ('cx', 'cy', 'cz'):
        require(documented.get(name, {}).get('type') == 'double precision', 'unit-vector column missing: '+name)
    evidence = []
    for alias, expression in columns.items():
        column = expression.split('.', 1)[1]
        info = documented.get(column)
        require(info is not None, 'undocumented column '+expression)
        require(info['type'] == schema_type(alias, columns), f'unexpected schema type {expression}: {info["type"]}')
        evidence.append(dict(alias=alias, expression=expression, table='s16a_wide.random', type=info['type'], unit=info['unit'],
                             description=info['description']))
    return evidence


def expected_rows():
    area = math.pi*RADIUS_ARCMIN**2
    return dict(area_arcmin2=round(area, 1), nominal_full_coverage=round(area*DENSITY_PER_ARCMIN2*ADJUST_MAX, 1),
                with_duplication=round(area*DENSITY_PER_ARCMIN2*ADJUST_MAX*DUPLICATION, 1))


def inputs(file):
    return {str(p): sha(p) for p in (LEDGER, SCHEMA_DIR/'sources.json', file)}


def code_files():
    here = Path(__file__).resolve()
    return [here, here.with_name('prepare_s23b_s16a_quality_supplement.py')]


def ledger():
    return [(f's16a_random_candidate_rows_{i:03d}.sql', lens) for i, lens in enumerate(lenses())]


def constants():
    return dict(radius_arcmin=RADIUS_ARCMIN, iadjust_density_max=ADJUST_MAX, documented_density_per_arcmin2=DENSITY_PER_ARCMIN2,
                assumed_duplication=DUPLICATION, row_limit_per_lens=LIMIT_PER_LENS, box_margin_deg=MARGIN_DEG)


def check_gate(plan):
    require(plan['query_plan_sha256'] == canonical_sha(plan), 'plan hash differs')
    require(plan['status'] == 'prepared_not_submitted' and not plan['query_submitted'] and plan['execution_requires_explicit_user_approval']
            and not plan['science_measurement'] and not plan['calibration_applied'], 'stage gate differs')
    for path, digest in {**plan['input_sha256'], **plan['code_sha256']}.items():
        require(sha(Path(path)) == digest, 'pinned input or code differs: '+path)


def validate(output):
    """Local exact-plan audit; regenerates every SQL from the pinned ledger. No network."""
    plan = json.loads((output/'plan.json').read_text())
    check_gate(plan)
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['output_columns'] == OUTPUT_COLUMNS and plan['column_expressions'] == COLUMNS,
            'requested fields differ')
    require(plan['cas_release_version'] == 'pdr1' and plan['source_table'] == 's16a_wide.random'
            and plan['maximum_total_received_bytes'] == CAP and plan['constants'] == constants(), 'release/budget/constants differ')
    file = schema_file()
    require(plan['input_sha256'] == inputs(file) and plan['schema_column_check'] == check_schema(file), 'schema evidence differs')
    expected = []
    for name, lens in ledger():
        sql = output/name
        require(sql.read_text() == query(lens), 'SQL scope differs '+name)
        expected.append(dict(name=name, lens=lens, row_limit=LIMIT_PER_LENS, expected_rows=expected_rows(),
                             sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    require(plan['batches'] == expected and plan['n_query_batches'] == len(expected) == N_LENSES, 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], lenses=N_LENSES,
                columns=len(OUTPUT_COLUMNS), batches=len(expected), total_sql_bytes=sum(b['sql_bytes'] for b in expected),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    file = schema_file()
    evidence = check_schema(file)
    items = ledger()
    output.mkdir(parents=True)
    batches = []
    for name, lens in items:
        sql = output/name
        sql.write_text(query(lens))
        batches.append(dict(name=name, lens=lens, row_limit=LIMIT_PER_LENS, expected_rows=expected_rows(),
                            sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    rows = expected_rows()['with_duplication']*N_LENSES
    per_row = dict(low=360, high=430)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_table='s16a_wide.random',
        purpose='Thinned S16A random points within 90 arcmin of each of the 11 pilot lenses, with S16A mask, edge, depth and '
                'seeing columns, to choose non-overlapping S16A-processed random centres and to rank them against the lens.',
        constants=constants(), n_lenses=N_LENSES, n_query_batches=len(batches), batches=batches,
        output_columns=OUTPUT_COLUMNS, column_expressions=COLUMNS, schema_column_check=evidence,
        maximum_total_received_bytes=CAP,
        estimated_rows=dict(total=round(rows), note='pi (90 arcmin)^2 x 100 per arcmin^2 x 0.003 x 1.41 per lens if the whole cone '
                            'is processed; the cones of S23B_8206 and S23B_8310 (92.6 arcmin apart) overlap, so some points '
                            'appear in both. The public maskUtils code seeds draws by patch consistently among filters, so '
                            'no per-filter multiplication is expected (plan review; not verified on the server).'),
        estimated_csv_bytes=dict(low=int(rows*per_row['low']), high=int(rows*per_row['high']),
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(OUTPUT_COLUMNS)} columns in CAS CSV style '
                                        '(Python-repr doubles, True/False, CRLF).', hard_stop='the 128 MiB per-plan byte budget'),
        query_form='No coord column exists in s16a_wide.random, so coneSearch is not used: a dec/ra BETWEEN box on the documented '
                   '(dec, ra) index, an exact cone from the cx/cy/cz dot product, iadjust_density thinning, ORDER BY object_id, and '
                   f'LIMIT {LIMIT_PER_LENS} (a batch that returns the limit is treated as truncated and refused). The dec index '
                   'range covers a 3-degree stripe, so a batch may take minutes.',
        not_filtered='No pixel, edge, bright-object, inner or primary flag is used in WHERE; masked points are kept on purpose. '
                     'Rows with NULL iadjust_density (no i-band processing) are dropped by the thinning predicate.',
        excluded=['source objects, shapes, photometry or photo-z', 'images',
                  'random centroids and sky_mean; g/r/z/y PSF moments, pix_variance and adjust_density; pixel flags other than those listed',
                  'S16A deep/udeep, PDR2/PDR3 tables'],
        limitations=['Random points carry no wl_fulldepth_fullcolor flag; candidates can be placed in processed patches but not '
                     'checked against the catalogue full-depth full-colour footprint by this query.',
                     'Random footprints have zero radius, so flags_pixel_bright_object_any on randoms is essentially a one-pixel '
                     'test, not an object-footprint test.',
                     'About 0.3 points per arcmin^2 after deduplication: 58-98 points per 3-Mpc disc, enough for PSF/depth ranking '
                     'and coarse mask fractions (about +/-3-4% at f=0.1), not for areas inside 0.4 Mpc.',
                     'Random PSF moments are SDSS moments in pixels (0.168 arcsec), not the forced adaptive moments in arcsec^2.'],
        analysis_plan=['Deduplicate with idetect_is_patch_inner and idetect_is_tract_inner; a point with NULL i-band inner flags or '
                       'icountinputs = 0 is i-unprocessed and counts as a hole, not as area.',
                       'Geometry (fixed now): candidate centres with 2 R_lens(3 Mpc) <= d <= 90 arcmin - R, i.e. more than 6 proper Mpc '
                       'from the lens and with the whole disc inside the cone; pairwise separation >= 6 proper Mpc; no overlap with any '
                       'other lens disc or candidate set; CAMIRA exclusion as in the existing selector; target at least 8 per lens '
                       'where the footprint allows.',
                       'Masks (fixed after the per-source flag plan): the hole flag and the bands whose edge flags define holes follow '
                       'the pre-registered decision rule of the bright-object flag plan; the same centre-mask rule applies to lenses.',
                       'Ranking: match the distribution of countinputs, isky_std and i PSF size within the disc to the lens disc; no '
                       'hard point-level countinputs cut (15-75% of the catalogue sources would fail (4,4,4,6,6)).',
                       'Never drop or rank centres on S16A source counts, and never reject a lens or random for a masked centre point '
                       '(area normalisation handles it; report the status).',
                       'Still to be fixed in a hashed selector before any count, shape-row, supplement or full-density random query '
                       'for new centres: deterministic greedy placement order, ranking metric and tolerance, minimum per lens, an '
                       'explicit CAMIRA overlap computation (catalogue files and hashes, radius, redshift) for the new positions, and '
                       'the S16A full-depth full-colour footprint (official map, or a separately approved aggregated forced '
                       'wl_fulldepth_fullcolor grid validated on the 23 empty pilot randoms).'],
        before_submission=['User approval of this exact plan hash and of the receiver code hashes recorded with it.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(file), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={name: sha(output/name) for name, _ in items},
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
