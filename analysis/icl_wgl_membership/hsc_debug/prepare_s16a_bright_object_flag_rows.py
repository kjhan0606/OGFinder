#!/usr/bin/env python3
"""Plan an S16A query of bright-object, edge and full-depth-full-colour flags for the lens-field sources.

The boost audit needs S16A footprint holes from the S16A (Sirius) bright-object mask. Before random
points can define holes, the per-source flags of the selected WL sources must show which flag the
weak-lensing selection already removes (bright_object_center, _any, or neither), and whether
wl_fulldepth_fullcolor is set for them. Reuses the executed ID ledger and batching. Local files only;
never submits SQL.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from prepare_s16a_afterburner_photoz_rows import batches_for, ledger_ids
from prepare_s23b_s16a_quality_supplement import (CHUNK, CAP, OUTPUT as SUPPLEMENT, canonical_sha, require, sha, write_json)
from submit_s23b_s16a_quality_supplement import APPROVED_SHA as SUPPLEMENT_SHA
from s16a_y1_source_contract import object_id

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_bright_object_flag_rows_plan_2026-11-25_03'
SCHEMA = 'astrafek-s16a-bright-object-flag-rows-plan-v3'
SCHEMA_URLS = {'forced': 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.forced.json',
               'meas': 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.meas.json'}
BANDS = 'grizy'
FLAGS = ('flags_pixel_bright_object_center', 'flags_pixel_bright_object_any', 'flags_pixel_edge')
TABLES = {'f': 'forced', 'm': 'meas'}
COLUMNS = {'object_id': 'r.object_id', 'forced_object_id': 'f.object_id', 'wl_fulldepth_fullcolor': 'f.wl_fulldepth_fullcolor'}
for _band in BANDS:
    for _flag in FLAGS:
        COLUMNS[f'{_band}{_flag}'] = f'f.{_band}{_flag}'
# meas carries the i-band flags that the local Table-4 reconstruction uses (plan review).
COLUMNS['meas_object_id'] = 'm.object_id'
for _flag in FLAGS:
    COLUMNS[f'meas_i{_flag}'] = f'm.i{_flag}'
N_OBJECTS = 115993
# Pre-registered decision rule (fixed before any data). T = number of True among the N requested catalogue sources.
APPLIED_MAX = max(3, int(1e-4*N_OBJECTS))
NOT_APPLIED_RANDOM_FRACTION = 0.5
MIN_TESTABLE_EXPECTED = 10*APPLIED_MAX  # 0.5 N f_rand must reach this, else UNTESTABLE
NULL_SOURCE_PLAN = 'cdd04062b79415328e7bbdae08a938a5224bf19daa2717d2d9665f81ed21b356'  # thinned S16A random-candidate plan


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in COLUMNS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{table} AS {alias} ON {alias}.object_id = r.object_id'
                      for alias, table in TABLES.items())
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def schema_files():
    """The forced and meas schemas pinned by the executed quality supplement."""
    manifest = json.loads((SUPPLEMENT/'sources.json').read_text())
    files = {}
    for table, url in SCHEMA_URLS.items():
        entry, = [d for d in manifest['documents'] if d['file'] == f'schema_{table}.json']
        file = SUPPLEMENT/entry['file']
        require(entry['url'] == url and sha(file) == entry['sha256'], table+' schema provenance/hash differs')
        files[table] = file
    return files


def check_schema(files):
    documented = {}
    for table, file in files.items():
        doc = json.loads(file.read_text())
        require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', table), 'schema identity differs')
        documented[table] = {}
        for row in doc['columns']:
            c = dict(zip(doc['headers'], row))
            for band in BANDS:
                documented[table][c['name'].replace('(g|r|i|z|y)', band)] = c
    evidence = []
    for alias, expression in COLUMNS.items():
        prefix, column = expression.split('.')
        if prefix == 'r':
            continue
        table = TABLES[prefix]
        info = documented[table].get(column)
        require(info is not None, 'undocumented column '+expression)
        want = ('bigint', '') if column == 'object_id' else ('bool', '')
        require((info['type'], info['unit']) == want, 'unexpected schema type/unit '+expression)
        evidence.append(dict(alias=alias, expression=expression, table='s16a_wide.'+table,
                             type=info['type'], unit=info['unit'], description=info['description']))
    return evidence


def batches():
    ids, supplement = ledger_ids()
    return ids, [(name.replace('s16a_afterburner_photoz_rows_', 's16a_bright_object_flag_rows_'), group)
                 for name, group in batches_for(ids, supplement)]


def inputs(files):
    paths = [SUPPLEMENT/'quality_supplement_plan.json', SUPPLEMENT/'object_ids.tsv', SUPPLEMENT/'sources.json', *files.values()]
    return {str(p): sha(p) for p in paths}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s16a_afterburner_photoz_rows.py',
                                                'prepare_s23b_s16a_quality_supplement.py', 's16a_y1_source_contract.py'))]


def validate(output):
    """Local exact-plan audit; regenerates every SQL from the verified ID ledger. No network."""
    plan = json.loads((output/'plan.json').read_text())
    require(plan['query_plan_sha256'] == canonical_sha(plan), 'plan hash differs')
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['status'] == 'prepared_not_submitted' and not plan['query_submitted']
            and plan['execution_requires_explicit_user_approval'] and not plan['science_measurement']
            and not plan['calibration_applied'], 'stage gate differs')
    require(plan['output_columns'] == list(COLUMNS) and plan['column_expressions'] == COLUMNS
            and plan['joined_tables'] == {k: 's16a_wide.'+v for k, v in TABLES.items()}, 'requested fields differ')
    require(plan['cas_release_version'] == 'pdr1' and plan['source_catalog'] == 'HSC S16A Wide'
            and plan['maximum_total_received_bytes'] == CAP and plan['batch_size'] == CHUNK, 'release/budget differs')
    for path, digest in {**plan['input_sha256'], **plan['code_sha256']}.items():
        require(sha(Path(path)) == digest, 'pinned input or code differs: '+path)
    for name, digest in plan['output_sha256'].items():
        require(Path(name).name == name and sha(output/name) == digest, 'plan member differs '+name)
    files = schema_files()
    require(plan['input_sha256'] == inputs(files) and plan['schema_column_check'] == check_schema(files), 'schema evidence differs')
    rule = plan['decision_rule']
    require((rule['applied_max_true'], rule['not_applied_random_fraction'], rule['min_testable_expected'], rule['null_source_plan_sha256'])
            == (APPLIED_MAX, NOT_APPLIED_RANDOM_FRACTION, MIN_TESTABLE_EXPECTED, NULL_SOURCE_PLAN), 'decision rule differs')
    ids, groups = batches()
    require(plan['id_scope_plan_sha256'] == SUPPLEMENT_SHA and plan['n_unique_objects'] == plan['expected_output_rows'] == len(ids) == N_OBJECTS,
            'row scope differs')
    expected = []
    for name, group in groups:
        sql = output/name
        require(sql.read_text() == query(group), 'SQL scope differs '+name)
        expected.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                             expected_rows=len(group), sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    require(plan['batches'] == expected and plan['n_query_batches'] == len(expected), 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'],
                unique_objects=len(ids), columns=len(COLUMNS), batches=len(expected),
                total_sql_bytes=sum(b['sql_bytes'] for b in expected), maximum_csv_bytes=CAP,
                remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    files = schema_files()
    evidence = check_schema(files)
    ids, groups = batches()
    output.mkdir(parents=True)
    ledger = []
    for name, group in groups:
        sql = output/name
        sql.write_text(query(group))
        ledger.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                           expected_rows=len(group), sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    per_row = dict(low=150, high=200)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='Establish which S16A mask flags the released weak-lensing catalogue already excludes (bright_object_center/_any, '
                'edge, wl_fulldepth_fullcolor), so that S16A random points can define footprint holes with the same mask.',
        id_scope_plan_sha256=SUPPLEMENT_SHA, id_scope='Exact 115993 object IDs and 24 batch boundaries of the executed quality supplement.',
        n_unique_objects=len(ids), n_new_objects=0, expected_output_rows=len(ids),
        batch_size=CHUNK, n_query_batches=len(ledger), batches=ledger,
        output_columns=list(COLUMNS), column_expressions=COLUMNS,
        joined_tables={k: 's16a_wide.'+v for k, v in TABLES.items()},
        schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(low=len(ids)*per_row['low'], high=len(ids)*per_row['high'],
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(COLUMNS)} columns (3 IDs, 19 booleans of 4-5 bytes, CRLF); the '
                                        'executed 13-column cmodel photo-z query was 15432619 bytes (133 bytes/row).'),
        included='Row-existence IDs, wl_fulldepth_fullcolor and grizy flags_pixel_bright_object_center/_any and flags_pixel_edge from '
                 's16a_wide.forced, and the same three i-band flags from s16a_wide.meas.',
        excluded=['new objects or sky regions', 'object shapes, fluxes and positions', 'other pixel flags',
                  'random points', 'photo-z', 'S17A/S19A/S23B tables'],
        join_policy='VALUES exact signed-bigint IDs, LEFT JOIN by object_id, no filtering. '
                    'Expect exactly one row per requested ID; abort on duplicate, missing or extra IDs or on a joined ID that differs.',
        server_sql_validation='Not submitted. Column names, types and units checked against the public forced schema (pinned by the executed '
                              'quality supplement; byte-identical to the 2026-11-24 public copy). The schemas list no index, but the executed '
                              'quality supplement already LEFT JOINed s16a_wide.forced and s16a_wide.meas by object_id in the same 24 batches.',
        decision_rule=dict(
            population='All 115993 requested IDs; every one is a member of s16a_wide.weaklensing_hsm_regauss (the executed shape '
                       'query used an INNER JOIN), so there is no unselected control group. N counts rows with a non-NULL value.',
            authority='meas i band decides the outcome (the Table-4 contract and the WL shapes use meas); forced gives g, r, z, y and '
                      'a per-object forced-versus-meas agreement count, reported but not used unless meas is AMBIGUOUS.',
            null_source_plan_sha256=NULL_SOURCE_PLAN,
            null_expectation='f_rand = flagged fraction of object_id-unique, inner-deduplicated, i-band non-edge points of the null-source '
                             'plan inside the union of the 43 non-empty source discs (3 proper Mpc); a source-weighted '
                             'sum(n_d f_d)/sum(n_d) over discs is reported beside it.',
            applied_max_true=APPLIED_MAX, not_applied_random_fraction=NOT_APPLIED_RANDOM_FRACTION, min_testable_expected=MIN_TESTABLE_EXPECTED,
            rule_center='meas iflags_pixel_bright_object_center: APPLIED if T <= applied_max_true; else UNTESTABLE if '
                        'not_applied_random_fraction x N x f_rand < min_testable_expected; else NOT APPLIED if T/N >= '
                        'not_applied_random_fraction x f_rand; else AMBIGUOUS. UNTESTABLE and AMBIGUOUS stop the hole definition.',
            rule_any='meas iflags_pixel_bright_object_any: APPLIED if T <= applied_max_true, otherwise NOT APPLIED (an applied cut '
                     'leaves exactly zero; edge-touching sources make the random null inapplicable).',
            rule_edge='Band b edge flags define random holes iff T_b(edge) <= applied_max_true (forced for g, r, z, y; meas for i). '
                      'If meas i edge T > applied_max_true, stop: the Table-4 contract would be inconsistent.',
            rule_fdfc='If T(wl_fulldepth_fullcolor = False) <= applied_max_true, full-depth full-colour is part of the selection and '
                      'must define the random footprint before the random-centre selector is hashed.',
            outcomes=['_center APPLIED and _any NOT APPLIED: the hole mask is the random iflags_pixel_bright_object_center.',
                      'both APPLIED: the catalogue used _any; random points (zero-radius footprints) give the _center geometry '
                      'only, and the footprint-dilation difference is booked as a selection effect.',
                      '_center NOT APPLIED: the database BRIGHT_OBJECT plane is not the catalogue mask; random flags cannot define '
                      'holes and region files would be needed.'],
            also=['Tabulate every flag for every band and table.', 'Cross-tabulate wl_fulldepth_fullcolor with the local forced countinputs.']),
        analysis_plan=['Apply the decision rule above exactly; report T per flag, band and table, f_rand, and the outcome.',
                       'No selection change, correction or science promotion.'],
        before_submission=['User approval of this exact plan hash and of the receiver code hashes recorded with it.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(files), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={name: sha(output/name) for name, _ in groups},
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
