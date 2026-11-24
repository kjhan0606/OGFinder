#!/usr/bin/env python3
"""Plan an S16A query of afterburner fluxes and Ephor_AB row existence; local files only.

Reuses the exact object-ID ledger and 5000-ID batching of the executed quality
supplement, so no new objects are requested. It never submits SQL or asks for
credentials. Row existence is read from the joined tables' object_id columns,
which separates an absent row from a row with NULL values.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from prepare_s23b_s16a_quality_supplement import (CHUNK, CAP, OUTPUT as SUPPLEMENT, canonical_sha, require, sha,
                                                  validate as validate_supplement, write_json)
from submit_s23b_s16a_quality_supplement import APPROVED_SHA as SUPPLEMENT_SHA
from s16a_y1_source_contract import object_id

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_afterburner_photoz_rows_plan_2026-11-24_02'
SOURCES = DATA/'s16a_afterburner_schema_sources_2026-11-24_01'
SCHEMA = 'astrafek-s16a-afterburner-photoz-rows-plan-v1'
BANDS = 'grizy'
N_OBJECTS = 115993
TABLES = {'b': 'afterburner', 'f': 'forced', 'z': 'photoz_ephor_ab'}
COLUMNS = {'object_id': 'r.object_id', 'afterburner_object_id': 'b.object_id'}
for _band in BANDS:
    COLUMNS[_band+'parent_flux_convolved_2_1'] = 'b.'+_band+'parent_flux_convolved_2_1'
for _band in BANDS:
    COLUMNS[_band+'parent_flux_convolved_2_1_flags'] = 'b.'+_band+'parent_flux_convolved_2_1_flags'
for _band in BANDS:
    COLUMNS[_band+'parent_flux_convolved_flag'] = 'b.'+_band+'parent_flux_convolved_flag'
for _band in BANDS:
    # Set when the local PSF is wider than the 1.09 arcsec target, so no flux is measured.
    COLUMNS[_band+'parent_flux_convolved_2_deconv'] = 'b.'+_band+'parent_flux_convolved_2_deconv'
for _band in BANDS:
    COLUMNS[_band+'parent_flux_convolved_seeing'] = 'b.'+_band+'parent_flux_convolved_seeing'
for _band in BANDS:
    COLUMNS[_band+'countinputs'] = 'f.'+_band+'countinputs'
COLUMNS.update(photoz_object_id='z.object_id', photoz_best='z.photoz_best', photoz_risk_best='z.photoz_risk_best')
TYPES = {'object_id': 'bigint', 'flux_convolved_2_1': 'double precision', 'flags': 'bool', 'flag': 'bool',
         'deconv': 'bool', 'seeing': 'real', 'countinputs': 'smallint', 'photoz_best': 'real', 'photoz_risk_best': 'real'}


def expected_type(column):
    return next(t for suffix, t in TYPES.items() if column.endswith(suffix))


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in COLUMNS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{table} AS {alias} ON {alias}.object_id = r.object_id'
                      for alias, table in TABLES.items())
    # LEFT JOIN keeps every requested ID; a NULL joined object_id marks an absent row.
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def schema_documents(sources):
    manifest = json.loads((sources/'sources.json').read_text())
    entry, = manifest['documents']
    require(entry['url'] == 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.afterburner.json'
            and sha(sources/entry['file']) == entry['sha256'], 'afterburner schema provenance/hash differs')
    files = {'afterburner': sources/entry['file'], 'forced': SUPPLEMENT/'schema_forced.json',
             'photoz_ephor_ab': SUPPLEMENT/'schema_photoz_ephor_ab.json'}
    return files


def check_schemas(files):
    docs = {}
    for table, file in files.items():
        doc = json.loads(file.read_text())
        require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', table), 'schema identity differs')
        columns = {}
        for row in doc['columns']:
            c = dict(zip(doc['headers'], row))
            for band in BANDS:
                columns[c['name'].replace('(g|r|i|z|y)', band)] = c
        docs[table] = columns
    evidence = []
    for alias, expression in COLUMNS.items():
        prefix, column = expression.split('.')
        if prefix == 'r':
            continue
        info = docs[TABLES[prefix]].get(column)
        require(info is not None, 'undocumented column '+expression)
        require(info['type'] == expected_type(column), 'unexpected schema type '+expression)
        evidence.append(dict(alias=alias, expression=expression, table='s16a_wide.'+TABLES[prefix],
                             type=info['type'], unit=info['unit'], description=info['description']))
    return evidence


def ledger_ids():
    checked = validate_supplement(SUPPLEMENT)
    require(checked['plan_sha256'] == SUPPLEMENT_SHA, 'ID ledger is not the executed supplement plan')
    with (SUPPLEMENT/'object_ids.tsv').open() as stream:
        ids = [object_id(r['object_id']) for r in csv.DictReader(stream, delimiter='\t')]
    require(len(ids) == len(set(ids)) == N_OBJECTS and ids == sorted(ids, key=int), 'ID ledger differs')
    supplement = json.loads((SUPPLEMENT/'quality_supplement_plan.json').read_text())
    return ids, supplement


def batches_for(ids, supplement):
    batches = []
    for i, start in enumerate(range(0, len(ids), CHUNK)):
        group = ids[start:start+CHUNK]
        old = supplement['batches'][i]
        require((old['first_object_id'], old['last_object_id'], old['expected_rows']) == (group[0], group[-1], len(group)),
                'batch boundaries differ from the executed supplement')
        batches.append((f's16a_afterburner_photoz_rows_{i:03d}.sql', group))
    require(len(batches) == supplement['n_query_batches'], 'batch count differs')
    return batches


def inputs(files):
    paths = [SUPPLEMENT/'quality_supplement_plan.json', SUPPLEMENT/'object_ids.tsv', SOURCES/'sources.json', *files.values()]
    return {str(p): sha(p) for p in paths}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s23b_s16a_quality_supplement.py', 's16a_y1_source_contract.py'))]


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
    for path, digest in plan['input_sha256'].items():
        require(sha(Path(path)) == digest, 'pinned input differs')
    for path, digest in plan['code_sha256'].items():
        require(sha(Path(path)) == digest, 'code differs')
    for name, digest in plan['output_sha256'].items():
        require(Path(name).name == name and sha(output/name) == digest, 'plan member differs '+name)
    files = schema_documents(SOURCES)
    require(plan['input_sha256'] == inputs(files) and plan['schema_column_check'] == check_schemas(files), 'schema evidence differs')
    ids, supplement = ledger_ids()
    require(plan['id_scope_plan_sha256'] == SUPPLEMENT_SHA and plan['n_unique_objects'] == plan['expected_output_rows'] == len(ids),
            'row scope differs')
    expected = []
    for name, group in batches_for(ids, supplement):
        file = output/name
        require(file.read_text() == query(group), 'SQL scope differs '+name)
        expected.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                             expected_rows=len(group), sql_bytes=file.stat().st_size, sha256=sha(file)))
    require(plan['batches'] == expected and plan['n_query_batches'] == len(expected), 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'],
                unique_objects=len(ids), columns=len(COLUMNS), batches=len(expected),
                total_sql_bytes=sum(b['sql_bytes'] for b in expected), maximum_csv_bytes=CAP,
                remote_query_submitted=False, science_measurement=False)


def build(output):
    output = Path(output).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    files = schema_documents(SOURCES)
    evidence = check_schemas(files)
    ids, supplement = ledger_ids()
    batches = batches_for(ids, supplement)
    output.mkdir(parents=True)
    ledger = []
    for name, group in batches:
        file = output/name
        file.write_text(query(group))
        ledger.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                           expected_rows=len(group), sql_bytes=file.stat().st_size, sha256=sha(file)))
    per_row = dict(low=280, high=560)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='Test whether a missing Ephor_AB row coincides with an absent or NULL/NaN afterburner '
                'parent_flux_convolved_2_1 flux in any band, and whether that flux was skipped because the local '
                'seeing exceeded the 1.09-arcsec target; compare per-band coadd input counts.',
        id_scope_plan_sha256=SUPPLEMENT_SHA, id_scope='Exact 115993 object IDs and 24 batch boundaries of the executed quality supplement.',
        n_unique_objects=len(ids), n_new_objects=0, expected_output_rows=len(ids),
        batch_size=CHUNK, n_query_batches=len(ledger), batches=ledger,
        output_columns=list(COLUMNS), column_expressions=COLUMNS,
        joined_tables={k: 's16a_wide.'+v for k, v in TABLES.items()},
        schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(low=len(ids)*per_row['low'], high=len(ids)*per_row['high'],
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(COLUMNS)} columns; '
                                        'the executed 31-column supplement was 45335621 bytes (391 bytes/row).'),
        included='Afterburner and Ephor_AB row-existence IDs; grizy parent_flux_convolved_2_1 with its two failure flags, '
                 'the 1.09-arcsec deconvolution flag and the afterburner seeing at the object; forced grizy countinputs; Ephor_AB best/risk.',
        excluded=['new objects or sky regions', 'flux errors, magnitudes and other afterburner apertures',
                  'shapes, positions, images, PSF images or models', 'full photo-z PDFs or other photo-z codes', 'S17A/S19A/S23B tables'],
        join_policy='VALUES exact signed-bigint IDs, LEFT JOIN each table by object_id, no filtering. '
                    'Expect exactly one row per requested ID; abort on duplicate, missing or extra IDs or on a joined ID that differs.',
        server_sql_validation='Not submitted. Column names and types checked against public schema documents only. '
                              'The afterburner schema lists no index, so the join cost on CAS is unknown.',
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        after_receipt=['Join to the unchanged object ledger; keep NULL and NaN distinct.',
                       'Report row absence and value missingness per band for objects with and without Ephor_AB best/risk.',
                       'Compare with the earlier supplement photo-z; no selection change, correction or science promotion.'],
        input_sha256=inputs(files), code_sha256={str(f): sha(f) for f in code_files()},
        output_sha256={name: sha(output/name) for name, _ in batches},
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
