#!/usr/bin/env python3
"""Plan the S16A source-pair and inner random-point queries for the 76 selected random centres; local only.

Selection s16a_random_centre_selection_20261126_01 (summary b12b346d) was approved by the user on 2026-11-26.
Per lens group two batches:
- pairs: the executed pilot pair query (plan _02, corrected shape-weight identifier) for every centre, unchanged,
  wrapped with a per-centre server count before LIMIT and a batch total;
- inner: the executed inner random-point query (plan e8b05d6e) for the same centres (0.4 Mpc + 30 arcsec).
The 11 lens centres already hold both products. Object-ID supplements (Y1 quality, photo-z, flags) need the
returned IDs and are a separate later plan. Never submits SQL.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import prepare_s16a_random_inner_rows as pi
import prepare_s23b_s16a_shape_qa_pilot as pp
from prepare_s16a_random_candidate_rows import check_gate, check_schema, label, schema_file, SCHEMA_DIR
from prepare_s23b_s16a_shape_qa_schema_fix import NEW as NEW_WEIGHT, OLD as OLD_WEIGHT, OUTPUT as PILOT_PLAN
from prepare_s23b_s16a_quality_supplement import CAP, canonical_sha, require, sha, write_json

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_random_source_rows_plan_2026-11-26_01'
SCHEMA = 'astrafek-s16a-random-source-rows-plan-v1'
SELECTION = Path('/home/kjhan/BACKUP/Astrafek/s16a_random_centre_selection_20261126_01')
SELECTION_SHA = 'b12b346d'
FDFC_VALIDATION = Path('/home/kjhan/BACKUP/Astrafek/s16a_fdfc_grid_validation_20261125_01/centre_fdfc.tsv')
N_CENTRES, N_LENSES = 76, 11
PAIR_LIMIT_PER_CENTRE = 20000  # pilot maximum 7,549; a low-FDFC lens scaled to full FDFC stays below 7,000
PAIR_COLUMNS = [*pp.OUTPUT_COLUMNS, 'n_matched', 'n_batch']
INNER_COLUMNS = pi.OUTPUT_COLUMNS
BYTES_PER_PAIR, BYTES_PER_POINT = 190, 245  # pilot receipts: 173 B per pair row + two counts; 241 B per inner point


def selected():
    summary = json.loads((SELECTION/'summary.json').read_text())
    require(sha(SELECTION/'summary.json').startswith(SELECTION_SHA), 'selection summary differs')
    require(sha(SELECTION/'selected_centres.tsv') == summary['selected_centres_sha256'], 'selected centres differ')
    with (SELECTION/'selected_centres.tsv').open(newline='') as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    require(len(rows) == N_CENTRES == summary['n_selected'] and len({r['center_id'] for r in rows}) == N_CENTRES
            and len({r['lens_id'] for r in rows}) == N_LENSES and all(r['center_id'].startswith('RC-') for r in rows), 'selection scope differs')
    return [dict(center_id=r['center_id'], lens_id=r['lens_id'], ra_deg=float(r['ra_deg']), dec_deg=float(r['dec_deg']),
                 z_lens=float(r['z_lens']), radius_arcsec=float(r['radius_arcsec']), fdfc_fraction=float(r['fdfc_fraction']))
            for r in rows]


def pair_subquery(centre):
    body = pp.query(centre)
    require(body.count(OLD_WEIGHT) == 1 and NEW_WEIGHT not in body, 'pilot weight identifier differs')
    body = body.replace(OLD_WEIGHT, NEW_WEIGHT).replace('\n', '\n    ')
    return (f'(SELECT * FROM (\n  SELECT s.*, count(*) OVER () AS n_matched FROM (\n    {body}\n  ) AS s\n'
            f'  ORDER BY s.object_id\n  LIMIT {PAIR_LIMIT_PER_CENTRE}\n) AS p)')


def pair_query(group):
    for c in group:
        label(c['center_id'])
    inner = '\nUNION ALL\n'.join(pair_subquery(c) for c in group)
    return f'SELECT u.*, count(*) OVER () AS n_batch FROM (\n{inner}\n) AS u;\n'


def inner_centre(c):
    return dict(center_id=c['center_id'], role='random', ra_deg=c['ra_deg'], dec_deg=c['dec_deg'], z_lens=c['z_lens'],
                radius_arcsec=c['radius_arcsec'], query_radius_arcsec=pi.radius_arcsec(c))


def expected_pairs():
    """Upper estimate per lens group: lens sources per unit FDFC fraction times the random FDFC fraction."""
    with FDFC_VALIDATION.open(newline='') as stream:
        lens = {r['lens_id']: r for r in csv.DictReader(stream, delimiter='\t') if r['role'] == 'lens'}
    require(len(lens) == N_LENSES, 'FDFC validation lens set differs')
    return {g: int(r['n_sources'])/max(float(r['fdfc_fraction']), 0.2) for g, r in lens.items()}


def groups():
    rows = selected()
    by_lens = {}
    for r in rows:
        by_lens.setdefault(r['lens_id'], []).append(r)
    per = expected_pairs()
    out = []
    for i, g in enumerate(sorted(by_lens)):
        members = sorted(by_lens[g], key=lambda r: r['center_id'])
        pairs = round(sum(per[g]*c['fdfc_fraction'] for c in members))
        points = round(sum(math.pi*(pi.radius_arcsec(c)/60)**2*pi.DENSITY_PER_ARCMIN2*pi.DUPLICATION for c in members))
        out.append((f's16a_random_pairs_{i:03d}.sql', 'pairs', g, members, pairs, pair_query(members)))
        out.append((f's16a_random_inner_{i:03d}.sql', 'inner', g, [inner_centre(c) for c in members], points,
                    pi.query([inner_centre(c) for c in members])))
    return out


def constants():
    return dict(pair_limit_per_centre=PAIR_LIMIT_PER_CENTRE, inner_constants=pi.constants(), pilot_weight_identifier=NEW_WEIGHT,
                bytes_per_pair=BYTES_PER_PAIR, bytes_per_point=BYTES_PER_POINT, n_centres=N_CENTRES, n_lenses=N_LENSES)


def pilot_schema_files():
    return [PILOT_PLAN/n for n in ('schema_column_check.json', 'schema_meas.json', 'schema_meas2.json', 'schema_weaklensing_hsm_regauss.json')]


def inputs():
    return {str(p): sha(p) for p in (SELECTION/'summary.json', SELECTION/'selected_centres.tsv', FDFC_VALIDATION,
                                     PILOT_PLAN/'PLAN_SHA256', SCHEMA_DIR/'sources.json', schema_file(), *pilot_schema_files())}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s16a_random_inner_rows.py', 'prepare_s16a_random_candidate_rows.py',
                                                 'prepare_s23b_s16a_shape_qa_pilot.py', 'prepare_s23b_s16a_shape_qa_schema_fix.py',
                                                 'prepare_s23b_s16a_quality_supplement.py', 'prepare_s23b_s16a_shape_count_plan.py'))]


def batches(output):
    out = []
    for n, kind, g, members, expected, _ in groups():
        limit = PAIR_LIMIT_PER_CENTRE if kind == 'pairs' else pi.LIMIT_PER_CENTRE
        out.append(dict(name=n, kind=kind, lens_id=g, centres=members, row_limit_per_centre=limit, expected_rows=expected,
                        sql_bytes=(output/n).stat().st_size, sha256=sha(output/n)))
    return out


def estimate(bs):
    pairs = sum(b['expected_rows'] for b in bs if b['kind'] == 'pairs')
    points = sum(b['expected_rows'] for b in bs if b['kind'] == 'inner')
    return dict(pairs=pairs, inner_points=points, csv_bytes=pairs*BYTES_PER_PAIR+points*BYTES_PER_POINT)


def validate(output):
    """Local exact-plan audit; regenerates every SQL. No network."""
    plan = json.loads((output/'plan.json').read_text())
    check_gate(plan)
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['constants'] == constants() and plan['maximum_total_received_bytes'] == CAP
            and plan['pair_output_columns'] == PAIR_COLUMNS and plan['inner_output_columns'] == INNER_COLUMNS
            and plan['input_sha256'] == inputs() and plan['inner_schema_column_check'] == check_schema(schema_file(), pi.COLUMNS),
            'plan scope differs')
    for n, _, _, _, _, sql in groups():
        require((output/n).read_text() == sql, 'SQL scope differs '+n)
    bs = batches(output)
    require(plan['batches'] == bs and len(bs) == 2*N_LENSES and sum(len(b['centres']) for b in bs) == 2*N_CENTRES, 'batch ledger differs')
    require(estimate(bs)['csv_bytes'] < CAP, 'estimated bytes exceed the per-plan cap')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], centres=N_CENTRES,
                batches=len(bs), total_sql_bytes=sum(b['sql_bytes'] for b in bs), estimate=estimate(bs),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    items = groups()
    output.mkdir(parents=True)
    for n, _, _, _, _, sql in items:
        (output/n).write_text(sql)
    bs = batches(output)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1',
        source_tables=['s16a_wide.meas2', 's16a_wide.weaklensing_hsm_regauss', 's16a_wide.meas', 's16a_wide.random'],
        purpose='Source pairs (3 proper Mpc) and full-density inner random points (0.4 proper Mpc + 30 arcsec) for the 76 '
                'random centres of the user-approved selection 2026-11-26_01, in the same form as the executed pilot queries.',
        selection_approval='user, 2026-11-26 ("커밋하고 진행" after the selection report)',
        constants=constants(), batches=bs, pair_output_columns=PAIR_COLUMNS, inner_output_columns=INNER_COLUMNS,
        pair_column_expressions={**pp.COLUMN_EXPRESSIONS, 'weight': NEW_WEIGHT},
        inner_column_expressions=pi.COLUMNS, inner_schema_column_check=check_schema(schema_file(), pi.COLUMNS),
        pair_schema_evidence='pilot plan _02 schema_column_check.json and schema JSONs (hashes in input_sha256)',
        maximum_total_received_bytes=CAP,
        estimated=dict(**estimate(bs), method='pairs: lens pilot sources / lens disc FDFC fraction (floor 0.2) x random FDFC '
                       'fraction; points: documented 100/arcmin^2 x duplication 1.41; bytes per row from the pilot receipts',
                       hard_stop='the 128 MiB per-plan byte budget'),
        query_form='pairs: pilot SELECT DISTINCT cone query, unchanged, inside a per-centre count(*) OVER () before ORDER BY '
                   'object_id LIMIT, UNION ALL per lens group with a batch count(*) OVER (). inner: plan e8b05d6e query, unchanged.',
        intake_checks=['Every centre returns rows; per-centre rows equal n_matched and stay below the limit; batch rows equal n_batch.',
                       'Pairs: unique (center_id, object_id), source inside the centre cone, finite shape and weight columns.',
                       'Inner: rows inside the 0.4 Mpc + 30 arcsec cone, unique object_id per centre.'],
        analysis_plan=['Object IDs from the pairs feed one later ID-limited plan (Y1 quality supplement, cmodel photo-z, '
                       'forced/meas flags) under its own approval and byte budget.',
                       'Boosts then use the accepted rule: FDFC areas, all 11 lens groups, expected-count stack, group jackknife.',
                       'Distant randoms (36-166 arcmin for the four wide lenses) are checked for depth and PSF against their lens.'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={b['name']: b['sha256'] for b in bs},
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
