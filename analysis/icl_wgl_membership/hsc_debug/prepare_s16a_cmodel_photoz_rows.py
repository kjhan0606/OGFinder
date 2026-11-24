#!/usr/bin/env python3
"""Plan an S16A query of photo-z codes that do not need afterburner fluxes; local files only.

photoz_ephor uses cmodel_{dev,exp} fluxes (schema description), so it should exist for
sources whose Ephor_AB row is missing. photoz_mizuki is an independent template code whose
input photometry the schema does not state; its row existence is itself a result.
Reuses the exact ID ledger and batching of the executed quality supplement. Never submits SQL.
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
OUTPUT = DATA/'s16a_cmodel_photoz_rows_plan_2026-11-24_01'
SOURCES = DATA/'s16a_photoz_schema_sources_2026-11-24_01'
SCHEMA = 'astrafek-s16a-cmodel-photoz-rows-plan-v1'
TABLES = {'e': 'photoz_ephor', 'm': 'photoz_mizuki'}
PREFIX = {'e': 'ephor', 'm': 'mizuki'}
FIELDS = {'e': ('object_id', 'photoz_best', 'photoz_risk_best', 'photoz_std_best'),
          'm': ('object_id', 'photoz_best', 'photoz_risk_best', 'photoz_std_best',
                'reduced_chisq', 'prob_star', 'prob_gal', 'prob_qso')}
TYPES = {'object_id': 'bigint'}
COLUMNS = {'object_id': 'r.object_id'}
for _alias in TABLES:
    for _field in FIELDS[_alias]:
        COLUMNS[PREFIX[_alias]+'_'+_field] = _alias+'.'+_field
N_OBJECTS = 115993


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in COLUMNS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{table} AS {alias} ON {alias}.object_id = r.object_id'
                      for alias, table in TABLES.items())
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def schema_files():
    manifest = json.loads((SOURCES/'sources.json').read_text())
    files = {}
    for entry in manifest['documents']:
        table = entry['file'].removeprefix('schema_').removesuffix('.json')
        require(entry['url'] == f'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.{table}.json'
                and sha(SOURCES/entry['file']) == entry['sha256'], 'schema provenance/hash differs')
        files[table] = SOURCES/entry['file']
    require(set(files) == set(TABLES.values()), 'schema set differs')
    return files


def check_schemas(files):
    evidence = []
    for alias, table in TABLES.items():
        doc = json.loads(files[table].read_text())
        require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', table), 'schema identity differs')
        require(any('object_id' in i for i in doc.get('indices', [])), 'no object_id index documented for '+table)
        columns = {dict(zip(doc['headers'], r))['name']: dict(zip(doc['headers'], r)) for r in doc['columns']}
        for field in FIELDS[alias]:
            require(field in columns and columns[field]['type'] == TYPES.get(field, 'real'), f'unexpected schema for {table}.{field}')
            evidence.append(dict(alias=PREFIX[alias]+'_'+field, expression=alias+'.'+field, table='s16a_wide.'+table,
                                 type=columns[field]['type'], description=columns[field]['description']))
        evidence.append(dict(alias=PREFIX[alias]+'_table_description', table='s16a_wide.'+table,
                             description=' '.join(doc['desc'].split())))
    return evidence


def batches():
    ids, supplement = ledger_ids()
    return ids, [(name.replace('s16a_afterburner_photoz_rows_', 's16a_cmodel_photoz_rows_'), group)
                 for name, group in batches_for(ids, supplement)]


def inputs(files):
    paths = [SUPPLEMENT/'quality_supplement_plan.json', SUPPLEMENT/'object_ids.tsv', SOURCES/'sources.json', *files.values()]
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
    require(plan['input_sha256'] == inputs(files) and plan['schema_column_check'] == check_schemas(files), 'schema evidence differs')
    ids, groups = batches()
    require(plan['id_scope_plan_sha256'] == SUPPLEMENT_SHA and plan['n_unique_objects'] == plan['expected_output_rows'] == len(ids) == N_OBJECTS,
            'row scope differs')
    expected = []
    for name, group in groups:
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
    files = schema_files()
    evidence = check_schemas(files)
    ids, groups = batches()
    output.mkdir(parents=True)
    ledger = []
    for name, group in groups:
        file = output/name
        file.write_text(query(group))
        ledger.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                           expected_rows=len(group), sql_bytes=file.stat().st_size, sha256=sha(file)))
    per_row = dict(low=120, high=240)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='Compare redshifts of sources without Ephor_AB (afterburner) photo-z against same-cell sources with it, '
                'using photo-z codes that do not require afterburner fluxes.',
        id_scope_plan_sha256=SUPPLEMENT_SHA, id_scope='Exact 115993 object IDs and 24 batch boundaries of the executed quality supplement.',
        n_unique_objects=len(ids), n_new_objects=0, expected_output_rows=len(ids),
        batch_size=CHUNK, n_query_batches=len(ledger), batches=ledger,
        output_columns=list(COLUMNS), column_expressions=COLUMNS,
        joined_tables={k: 's16a_wide.'+v for k, v in TABLES.items()},
        schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(low=len(ids)*per_row['low'], high=len(ids)*per_row['high'],
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(COLUMNS)} columns; the executed '
                                        '35-column afterburner query was 38752722 bytes (334 bytes/row).'),
        included='Row-existence IDs and photoz_best, photoz_risk_best, photoz_std_best from photoz_ephor (cmodel inputs) and photoz_mizuki, '
                 'plus Mizuki fit diagnostics reduced_chisq and prob_star/prob_gal/prob_qso.',
        excluded=['new objects or sky regions', 'full photo-z PDFs and other point estimates', 'Mizuki physical parameters (mass, SFR, rest-frame magnitudes)',
                  'other photo-z codes (DEmP, FRANKEN-Z, MLZ, NNPZ)', 'S17A/S19A/S23B tables'],
        join_policy='VALUES exact signed-bigint IDs, LEFT JOIN each table by object_id, no filtering. '
                    'Expect exactly one row per requested ID; abort on duplicate, missing or extra IDs or on a joined ID that differs.',
        server_sql_validation='Not submitted. Column names and types checked against public schema documents; both tables document a unique object_id index.',
        analysis_plan=['Calibrate Ephor(cmodel) against Ephor_AB z_best on sources that have both, as a function of the local '
                       'afterburner r-band FWHM. Sources with Ephor_AB have r FWHM <= 1.092 arcsec, missing ones 1.092-1.171 arcsec, '
                       'so use the near-threshold finite sources (r FWHM 1.0-1.092 arcsec) and report the trend with FWHM; '
                       'do not apply an all-source calibration to the missing sources.',
                       'Compare Ephor(cmodel) and Mizuki z_best of photo-z-missing vs finite sources with the same-cell stratified permutation test, by cell and by cell x i-mag.',
                       'Per lens, compare the lensing-efficiency-weighted mean of Ephor(cmodel) z for all selected sources against Ephor_AB-available sources only.',
                       'No selection change, correction or science promotion.'],
        before_submission=['User approval of this exact plan hash.',
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
