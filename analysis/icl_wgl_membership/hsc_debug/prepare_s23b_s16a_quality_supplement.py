#!/usr/bin/env python3
"""Build an ID-limited S16A quality/photo-z plan using local files only.

Never submits SQL, opens a network connection or requests credentials.
Schema JSONs must be pre-fetched public documentation, not CAS query results.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import shutil

from fetch_ancillary_catalogues import catalogue_reader
from submit_s23b_s16a_shape_qa_pilot import validate_plan, verify_complete, CORRECTED_PLAN
from s16a_y1_source_contract import EXPRESSIONS, TABLES, REQUIRED, contract, object_id, select_row

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'wl_s23b_s16a_quality_supplement_plan_2026-11-18_01'
SOURCES = DATA/'s16a_y1_input_contract_sources_2026-11-18_01'
RECEIPT_SHA = '373b1bca841d733bd4269421dd00fb2b6afb89295cc2f39dbab3d8424bd45c6c'
CHUNK = 5000
CAP = 128 * 1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(file):
    return hashlib.sha256(Path(file).read_bytes()).hexdigest()


def canonical_sha(plan):
    value = dict(plan)
    value.pop('query_plan_sha256', None)
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'),
                                     ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def write_json(file, value):
    file.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in EXPRESSIONS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{table} AS {alias} ON {alias}.object_id = r.object_id'
                      for alias, table in TABLES.items())
    # LEFT JOIN deliberately retains missing photometry/photo-z as NULL. No cuts,
    # footprint expansion, new shape table, DISTINCT, COALESCE or LIMIT.
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def check_schemas(sources):
    manifest = json.loads((sources/'sources.json').read_text())
    docs = {}
    for table in TABLES.values():
        filename = 'schema_'+table+'.json'
        entry = next(x for x in manifest['documents'] if x['file'] == filename)
        url = 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.'+table+'.json'
        require(entry['url'] == url and sha(sources/filename) == entry['sha256'], 'schema provenance/hash differs')
        doc = json.loads((sources/filename).read_text())
        require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', table), 'schema identity differs')
        # Expand documented five-band macro, including a_(g|r|i|z|y).
        columns = {}
        for row in doc['columns']:
            c = dict(zip(doc['headers'], row))
            for band in 'grizy':
                columns[c['name'].replace('(g|r|i|z|y)', band)] = c
        require(columns['object_id']['type'] == 'bigint', 'join ID is not bigint')
        docs[table] = columns
    evidence = []
    for alias, expression in EXPRESSIONS.items():
        a, column = expression.split('.')
        if a == 'r':
            continue
        require(column in docs[TABLES[a]], 'unknown documented column '+expression)
        info = docs[TABLES[a]][column]
        expected = 'bool' if alias in ('is_primary', *contract_flags()) else None
        require(info['type'] == expected if expected else info['type'] in ('real', 'double precision'),
                'unexpected schema type '+expression)
        evidence.append(dict(alias=alias, expression=expression, table=TABLES[a],
                             type=info['type'], unit=info['unit']))
    return evidence


def contract_flags():
    from s16a_y1_source_contract import BAD_FLAGS
    return BAD_FLAGS


def load_pilot(directory):
    plan = validate_plan(directory)
    receipt_file = directory/'results/download_summary.json'
    require(sha(receipt_file) == RECEIPT_SHA, 'completed pilot receipt differs')
    state = json.loads((directory/'results/shape_qa_status.json').read_text())
    require(state['status'] == 'complete', 'pilot is incomplete')
    verify_complete(plan, directory/'results', state)
    rows, inputs = [], {str(receipt_file): sha(receipt_file),
                        str(directory/'shape_qa_plan.json'): sha(directory/'shape_qa_plan.json')}
    for batch in plan['batches']:
        file = directory/'results'/(Path(batch['name']).stem+'.csv')
        inputs[str(file)] = sha(file)
        with file.open() as stream:
            rows.extend(dict(r) for r in catalogue_reader(stream))
    require(len(rows) == 169535, 'pilot pair count differs')
    ids = sorted({object_id(r['object_id']) for r in rows}, key=int)
    require(len(ids) == 115993, 'pilot unique object count differs')
    require(len({(r['center_id'], r['object_id']) for r in rows}) == len(rows), 'duplicate input pair')
    return plan, rows, ids, inputs


def validate(output):
    """Local exact-plan audit; regenerates every SQL and its complete ID list."""
    plan = json.loads((output/'quality_supplement_plan.json').read_text())
    require(plan['query_plan_sha256'] == canonical_sha(plan), 'plan hash differs')
    require((output/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    for name, digest in plan['output_sha256'].items():
        require(Path(name).name == name and sha(output/name) == digest, 'plan member differs '+name)
    for path, digest in plan['input_sha256'].items():
        require(sha(Path(path)) == digest, 'pinned input differs')
    for path, digest in plan['code_sha256'].items():
        require(sha(Path(path)) == digest, 'code differs')
    require(plan['status'] == 'prepared_not_submitted' and not plan['query_submitted']
            and plan['execution_requires_explicit_user_approval'] and not plan['science_measurement'], 'stage gate differs')
    require(plan['output_columns'] == list(EXPRESSIONS) and plan['column_expressions'] == EXPRESSIONS,
            'requested fields differ')
    require(plan['cas_release_version'] == 'pdr1' and plan['source_catalog'] == 'HSC S16A Wide'
            and plan['maximum_total_received_bytes'] == CAP
            and plan['batch_size'] == CHUNK, 'release/budget differs')
    source_plan, rows, ids, inputs = load_pilot(Path(plan['pilot_directory']))
    require(plan['pilot_executed_plan_sha256'] == source_plan['query_plan_sha256']
            and plan['input_sha256'] == inputs, 'pilot inputs differ')
    require(plan['n_unique_objects'] == len(ids) and plan['expected_output_rows'] == len(ids)
            and plan['n_original_source_pairs'] == len(rows), 'row budget differs')
    with (output/'object_ids.tsv').open() as f:
        require(list(csv.DictReader(f, delimiter='\t')) == [{'object_id': i} for i in ids], 'ID ledger differs')
    with (output/'center_object_pairs.tsv').open() as f:
        require(list(csv.DictReader(f, delimiter='\t')) == [dict(center_id=r['center_id'], object_id=r['object_id']) for r in rows], 'pair ledger differs')
    expected_batches = []
    for i, start in enumerate(range(0, len(ids), CHUNK)):
        group = ids[start:start+CHUNK]
        file = output/f's16a_quality_supplement_{i:03d}.sql'
        require(file.read_text() == query(group), 'SQL scope differs')
        expected_batches.append(dict(name=file.name, first_object_id=group[0], last_object_id=group[-1],
                                     expected_rows=len(group), sql_bytes=file.stat().st_size, sha256=sha(file)))
    require(plan['batches'] == expected_batches and plan['n_query_batches'] == len(expected_batches), 'batch ledger differs')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'],
                unique_objects=len(ids), columns=len(EXPRESSIONS), batches=len(expected_batches),
                total_sql_bytes=sum(b['sql_bytes'] for b in expected_batches),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(output, sources):
    output, sources = Path(output).resolve(), Path(sources).resolve()
    require(not output.exists(), 'refusing to overwrite plan')
    evidence = check_schemas(sources)
    pilot, rows, ids, inputs = load_pilot(CORRECTED_PLAN)
    # A readiness audit, not a spurious retained-count estimate from missing cuts.
    diagnostic = Counter()
    available = set(rows[0])
    for row in rows:
        diagnostic[select_row(row)['table4_diagnostic']] += 1
    output.mkdir(parents=True, exist_ok=False)
    for f in sources.iterdir():
        if f.is_file():
            shutil.copy2(f, output/f.name)
    write_json(output/'schema_column_check.json', evidence)
    write_json(output/'y1_input_contract.json', contract())
    for filename, fields, records in (
        ('object_ids.tsv', ['object_id'], [dict(object_id=i) for i in ids]),
        ('center_object_pairs.tsv', ['center_id', 'object_id'], rows),
    ):
        with (output/filename).open('w') as f:
            writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore', lineterminator='\n')
            writer.writeheader(); writer.writerows(records)
    write_json(output/'existing_input_readiness.json', dict(
        source_pairs=len(rows), unique_objects=len(ids),
        available_table4_inputs=sorted(available.intersection(REQUIRED)),
        absent_table4_inputs=sorted(set(REQUIRED)-available),
        table4_diagnostic_counts=dict(diagnostic), physical_retained_fraction=None,
        point_photoz_available=False, calibrated_background_selection_available=False,
        new_query_submitted=False, science_measurement=False))
    batches = []
    for i, start in enumerate(range(0, len(ids), CHUNK)):
        group = ids[start:start+CHUNK]
        file = output/f's16a_quality_supplement_{i:03d}.sql'; file.write_text(query(group))
        batches.append(dict(name=file.name, first_object_id=group[0], last_object_id=group[-1],
                            expected_rows=len(group), sql_bytes=file.stat().st_size, sha256=sha(file)))
    plan = dict(
        schema='astrafek-s16a-quality-supplement-plan-v1', status='prepared_not_submitted',
        execution_requires_explicit_user_approval=True, query_submitted=False,
        new_catalog_rows_received=False, science_measurement=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        pilot_directory=str(CORRECTED_PLAN), pilot_executed_plan_sha256=pilot['query_plan_sha256'],
        n_unique_objects=len(ids), n_original_source_pairs=len(rows), n_lenses=11,
        n_centers_in_original_ledger=66, n_positive_centers_previously_queried=43,
        n_new_centers=0, n_new_objects=0, expected_output_rows=len(ids),
        batch_size=CHUNK, n_query_batches=len(batches), batches=batches,
        output_columns=list(EXPRESSIONS), column_expressions=EXPRESSIONS,
        maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(low=len(ids)*128+len(batches)*2000,
                                 high=len(ids)*512+len(batches)*2000,
                                 method='128-512 bytes/object plus 2000 bytes/batch header; estimate, not observed'),
        included='Exact pilot object IDs; raw selection flags, unforced grizy cmodel photometry, i extinction, raw sigma/resolution and Ephor_AB best/risk.',
        excluded=['new sky centers or new objects', 'source positions or repeated shapes/weights/m/c',
                  'images/PSFs/PSF moments', 'full photo-z PDFs', 'COSMOS calibration files', 'S19A/Y3 catalogs'],
        join_policy='VALUES exact signed-bigint IDs with LEFT JOIN by object_id; no filtering. Expect one row per requested ID; abort on duplicate/missing/extra IDs.',
        server_sql_validation='Not submitted or tested on CAS; documentation-only identifier/type verification.',
        before_submission=['Separate user approval of this exact scope/hash.',
                           'Implement and test a serial resumable receiver, pin plan and SQL, enforce cumulative cap including retries, and verify ID/column equality.',
                           'First approved batch is the SQL compatibility check; stop on error/ambiguity, no automatic scope expansion.'],
        after_receipt=['Retain NULLs; exact one-row-per-ID join to the unchanged original pair ledger.',
                       'Run Table-4 diagnostics and report attrition/tails/radial support, including known-zero controls.',
                       'Point-z/risk completeness only; no new redshift cuts, shear, mass, boost or science promotion.'],
        input_sha256=inputs,
        code_sha256={str(f):sha(f) for f in (Path(__file__).resolve(),
            Path(__file__).with_name('s16a_y1_source_contract.py').resolve())},
        output_sha256={f.name:sha(f) for f in sorted(output.iterdir()) if f.is_file()},
    )
    plan['query_plan_sha256'] = canonical_sha(plan)
    write_json(output/'quality_supplement_plan.json', plan)
    (output/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  quality_supplement_plan.json\n')
    return validate(output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=OUTPUT)
    parser.add_argument('--schemas', type=Path, default=SOURCES)
    parser.add_argument('--verify', action='store_true', help='Verify an existing local plan; no network')
    args = parser.parse_args()
    print(json.dumps(validate(args.output.resolve()) if args.verify else build(args.output, args.schemas), indent=2))
