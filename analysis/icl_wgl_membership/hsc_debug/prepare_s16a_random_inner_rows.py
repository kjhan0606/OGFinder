#!/usr/bin/env python3
"""Plan a full-density S16A random-point query inside 0.4 proper Mpc (+30 arcsec) of the 66 pilot centres.

The 2026-11-25 boost audit found the 0-0.1 Mpc boost changes by about 0.1 depending on whether label area
dilated from neighbouring sources into source-free cores is counted. Full-density S16A random points
(100 per arcmin^2, no thinning) with S16A (Sirius) bright-object, edge and depth flags measure the
processed and unmasked area of each centre's inner discs directly, independent of source positions,
and show whether each centre point (lens and BCG included) is masked. The 30-arcsec margin covers the
footprint-dilation scale. No flag filtering. Local files only; never submits SQL.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from prepare_s16a_random_candidate_rows import (DENSITY_PER_ARCMIN2, DUPLICATION, LEDGER, MATCH_COUNT, SCHEMA_DIR, centres, check_gate,
                                                check_schema, kinds, label, random_columns, schema_file, where)
from prepare_s23b_s16a_quality_supplement import CAP, canonical_sha, require, sha, write_json

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_random_inner_rows_plan_2026-11-25_01'
SCHEMA = 'astrafek-s16a-random-inner-rows-plan-v1'
INNER_MPC = 0.4
DISC_MPC = 3.0
MARGIN_ARCSEC = 30.0
LIMIT_PER_CENTRE = 6000
N_CENTRES = 66
COLUMNS = random_columns(thinned=False)
# n_matched per centre (window before LIMIT) and n_batch over the whole UNION detect LIMIT cuts and short responses,
# including a dropped trailing centre (the lens centre must return rows, so n_batch is always present).
OUTPUT_COLUMNS = ['center_id', *COLUMNS, 'n_matched', 'n_batch']
INTEGERS, BOOLEANS = kinds(COLUMNS)
FLOATS = tuple(c for c in COLUMNS if c not in INTEGERS and c not in BOOLEANS)


def radius_arcsec(centre):
    """0.4 proper Mpc at the lens redshift (the ledger radius is 3 proper Mpc) plus the footprint margin."""
    return float(centre['radius_arcsec'])*INNER_MPC/DISC_MPC+MARGIN_ARCSEC


def subquery(centre):
    ra, dec = float(centre['ra_deg']), float(centre['dec_deg'])
    columns = ',\n    '.join([f"{label(centre['center_id'])} AS center_id", *(f'{expr} AS {alias}' for alias, expr in COLUMNS.items()), MATCH_COUNT])
    body = where(ra, dec, radius_arcsec(centre)/3600).replace('\n  AND', '\n    AND')
    return (f'(SELECT * FROM (\n  SELECT\n    {columns}\n  FROM s16a_wide.random AS r\n  {body}'
            f'  ORDER BY r.object_id\n  LIMIT {LIMIT_PER_CENTRE}\n) AS c)')


def query(group):
    inner = '\nUNION ALL\n'.join(subquery(c) for c in group)
    return f'SELECT u.*, count(*) OVER () AS n_batch FROM (\n{inner}\n) AS u;\n'


def groups():
    rows = centres()
    by_lens = {}
    for r in rows:
        by_lens.setdefault(r['matched_lens_id'], []).append(r)
    require(len(rows) == N_CENTRES and all(len(v) == 6 for v in by_lens.values()), 'centre groups differ')
    out = []
    for i, lens_id in enumerate(sorted(by_lens)):
        members = sorted(by_lens[lens_id], key=lambda r: (r['role'] != 'lens', r['center_id']))
        out.append((f's16a_random_inner_rows_{i:03d}.sql', lens_id,
                    [dict(center_id=r['center_id'], role=r['role'], ra_deg=r['ra_deg'], dec_deg=r['dec_deg'], z_lens=r['z_lens'],
                          radius_arcsec=r['radius_arcsec'], query_radius_arcsec=radius_arcsec(r)) for r in members]))
    return out


def expected_rows(group):
    nominal = sum(math.pi*(c['query_radius_arcsec']/60)**2*DENSITY_PER_ARCMIN2 for c in group)
    return dict(nominal_full_coverage=round(nominal, 1), with_duplication=round(nominal*DUPLICATION, 1))


def constants():
    return dict(inner_mpc=INNER_MPC, disc_mpc=DISC_MPC, margin_arcsec=MARGIN_ARCSEC, row_limit_per_centre=LIMIT_PER_CENTRE,
                thinning=None, documented_density_per_arcmin2=DENSITY_PER_ARCMIN2, assumed_duplication=DUPLICATION)


def inputs(file):
    return {str(p): sha(p) for p in (LEDGER, SCHEMA_DIR/'sources.json', file)}


def code_files():
    here = Path(__file__).resolve()
    return [here, here.with_name('prepare_s16a_random_candidate_rows.py'), here.with_name('prepare_s23b_s16a_quality_supplement.py')]


def ledger(output):
    out = []
    for name, lens_id, group in groups():
        sql = output/name
        out.append(dict(name=name, lens_id=lens_id, centres=group, row_limit_per_centre=LIMIT_PER_CENTRE,
                        expected_rows=expected_rows(group), sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    return out


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
    require(plan['input_sha256'] == inputs(file) and plan['schema_column_check'] == check_schema(file, COLUMNS), 'schema evidence differs')
    for name, _, group in groups():
        require((output/name).read_text() == query(group), 'SQL scope differs '+name)
    expected = ledger(output)
    require(plan['batches'] == expected and plan['n_query_batches'] == len(expected) == 11
            and sum(len(b['centres']) for b in expected) == N_CENTRES, 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], centres=N_CENTRES,
                columns=len(OUTPUT_COLUMNS), batches=len(expected), total_sql_bytes=sum(b['sql_bytes'] for b in expected),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    file = schema_file()
    evidence = check_schema(file, COLUMNS)
    items = groups()
    output.mkdir(parents=True)
    for name, _, group in items:
        (output/name).write_text(query(group))
    batches = ledger(output)
    rows = sum(b['expected_rows']['with_duplication'] for b in batches)
    per_row = dict(low=240, high=300)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_table='s16a_wide.random',
        purpose='Full-density S16A random points within 0.4 proper Mpc + 30 arcsec of the 11 lens and 55 random pilot centres, to '
                'measure S16A processed and unmasked area in the inner annuli independently of source positions and to record '
                'whether each centre point is masked.',
        constants=constants(), n_centres=N_CENTRES, n_query_batches=len(batches), batches=batches,
        output_columns=OUTPUT_COLUMNS, column_expressions=COLUMNS, schema_column_check=evidence,
        maximum_total_received_bytes=CAP,
        estimated_rows=dict(total=round(rows), note='pi r^2 x 100 per arcmin^2 x 1.41 per centre for r = 0.4 Mpc + 30 arcsec '
                            '(93-112 arcsec) if fully processed; empty pilot randoms may return few or no rows.'),
        estimated_csv_bytes=dict(low=int(rows*per_row['low']), high=int(rows*per_row['high']),
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(OUTPUT_COLUMNS)} columns in CAS CSV style.',
                                 hard_stop='the 128 MiB per-plan byte budget'),
        query_form='One batch per lens group: six parenthesised per-centre subqueries joined by UNION ALL (the form of the executed '
                   'PDR3 control-candidate query), each a dec/ra BETWEEN box, an exact cx/cy/cz cone, ORDER BY object_id, '
                   f'count(*) OVER () AS n_matched and LIMIT {LIMIT_PER_CENTRE}, wrapped in an outer SELECT adding count(*) OVER () '
                   'AS n_batch. The receiver requires rows == n_matched per centre and == n_batch per batch. No thinning.',
        not_filtered='No pixel, edge, bright-object, inner, primary or adjust_density predicate; masked and unprocessed-band points '
                     'are kept on purpose.',
        excluded=['source objects, shapes, photometry or photo-z', 'images', 'random centroids, PSF moments, sky statistics, '
                  'adjust_density and other pixel flags', 'points beyond 0.4 Mpc + 30 arcsec', 'other releases'],
        limitations=['Random footprints have zero radius, so bright_object_any on randoms is essentially a one-pixel test.',
                     'Random points carry no wl_fulldepth_fullcolor flag.',
                     'The 0-0.1 Mpc disc (radius 16-20 arcsec) holds about 22-36 points per centre, so its area carries binomial '
                     'noise where it is partly masked.'],
        analysis_plan=['Deduplicate with idetect_is_patch_inner and idetect_is_tract_inner; points with NULL i-band inner flags, '
                       'icountinputs = 0 or iflags_pixel_edge are i-unprocessed and count as holes.',
                       'Area per centre and annulus (0-0.1, 0.1-0.2, 0.2-0.4 Mpc) = geometric annulus area x (deduplicated points that '
                       'pass the hole rule) / (deduplicated points present); the hole rule follows the bright-object flag plan outcome. '
                       'Counts over the expected 100 per arcmin^2 only flag missing patches (deficit > 3 sigma or skymap_id gaps), '
                       'because they carry 17-21% binomial noise in 0-0.1 Mpc.',
                       'If the 0-0.1 Mpc points are present and pass the hole rule, the full area is kept (the recorded convention) '
                       'and the own-source rule is not used for that centre.',
                       'Centre-point status = the majority flag among deduplicated points within 6 arcsec of the centre; it is '
                       'reported, never used to reject a centre.',
                       'Compare with the source-dilation areas of the 2026-11-25 audit, in particular RC-6065306c1c2aaca36993 and '
                       'RC-1d54db20ff04be924a4f, and recompute the 0-0.4 Mpc boosts with random-point areas as a sensitivity.',
                       'No selection change, correction or science promotion.'],
        before_submission=['User approval of this exact plan hash and of the receiver code hashes recorded with it.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(file), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={b['name']: b['sha256'] for b in batches},
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
