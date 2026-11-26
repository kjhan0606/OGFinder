#!/usr/bin/env python3
"""Plan an ID-limited S16A query of local-background and deblending columns for the standard-boost sources; local only.

Purpose: candidate (a) of the central |e| deficit (BCG/ICL light). No sky-level column exists in s16a_wide; the
documented proxies are the forced variance at the object position and input count, the parent flux and the parent
and child second moments measured with the child-matched Gaussian weight (light from neighbours and the BCG/ICL),
and the meas deblending flags. IDs: exactly the ledger of the received P(z)-summary plan (be6bd1a9), reused unchanged.
Form: FROM VALUES ids LEFT JOIN meas, forced; 5,000 IDs per batch, exactly one row per ID. Never submits SQL.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from prepare_s16a_photoz_interval_rows import OUTPUT as PZ_PLAN, items, ledger_ids
from prepare_s16a_random_candidate_rows import check_gate
from prepare_s23b_s16a_quality_supplement import CAP, CHUNK, canonical_sha, require, sha, write_json
from s16a_y1_source_contract import object_id

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_local_background_rows_plan_2026-11-27_01'
SOURCES = DATA/'hsc_schema_tables_2026-11-24_01'
PZ_APPROVED = 'be6bd1a95bee49638a5d5337fa4c63f3cc93f8ac77d40a37e3fafe4a558c60cb'
SCHEMA = 'astrafek-s16a-local-background-rows-plan-v1'
TABLES = {'m': 'meas', 'f': 'forced'}
COLUMNS = {
    'object_id': ('r', 'object_id', None),
    'meas_object_id': ('m', 'object_id', 'bigint'),
    'forced_object_id': ('f', 'object_id', 'bigint'),
    'i_variance': ('f', 'ivariance', 'real'),
    'g_variance': ('f', 'gvariance', 'real'),
    'i_countinputs': ('f', 'icountinputs', 'smallint'),
    'g_countinputs': ('f', 'gcountinputs', 'smallint'),
    'parent_id': ('m', 'parent_id', 'bigint'),
    'deblend_nchild': ('m', 'deblend_nchild', 'integer'),
    'i_blendedness_abs_flux_child': ('m', 'iblendedness_abs_flux_child', 'double precision'),
    'i_blendedness_abs_flux_parent': ('m', 'iblendedness_abs_flux_parent', 'double precision'),
    **{f'i_blendedness_abs_shape_{w}_{k}': ('m', f'iblendedness_abs_shape_{w}_{k}', 'real') for w in ('child', 'parent') for k in ('11', '22', '12')},
    **{f'i_deblend_{k}': ('m', f'ideblend_{k}', 'bool') for k in ('has_stray_flux', 'too_many_peaks', 'parent_too_big', 'masked')},
}
EXPRESSIONS = {alias: f'{t}.{c}' for alias, (t, c, _) in COLUMNS.items()}
JOINED_IDS = ('meas_object_id', 'forced_object_id')
BYTES_PER_ROW = 300  # three 19-digit IDs, parent_id, two doubles, eight reals, small ints and four booleans


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in EXPRESSIONS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{t} AS {a} ON {a}.object_id = r.object_id' for a, t in TABLES.items())
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def schema_files():
    manifest = json.loads((SOURCES/'sources.json').read_text())
    out = {}
    for table in TABLES.values():
        name = f'pdr1.s16a_wide.{table}.json'
        entry = next(d for d in manifest['documents'] if d['file'] == name)
        require(entry['url'] == 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/'+name and sha(SOURCES/name) == entry['sha256'],
                'schema provenance/hash differs '+name)
        out[table] = SOURCES/name
    return out


def check_schema(files):
    docs = {}
    for table, file in files.items():
        doc = json.loads(file.read_text())
        require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', table), 'schema identity differs')
        cols = {}
        for row in doc['columns']:
            c = dict(zip(doc['headers'], row))
            for band in 'grizy':
                cols[c['name'].replace('(g|r|i|z|y)', band)] = c
        docs[table] = cols
    evidence = []
    for alias, (t, c, typ) in COLUMNS.items():
        if t == 'r':
            continue
        col = docs[TABLES[t]].get(c)
        require(col is not None and col['type'] == typ, f'unexpected schema for {t}.{c}')
        evidence.append(dict(alias=alias, table=TABLES[t], column=c, type=col['type'], unit=col['unit'], description=col['description']))
    return evidence


def pz_ledger():
    plan = json.loads((PZ_PLAN/'plan.json').read_text())
    require(plan['query_plan_sha256'] == PZ_APPROVED and canonical_sha(plan) == PZ_APPROVED, 'P(z) plan differs from the approved one')
    ledger = PZ_PLAN/'object_ids.tsv'
    require(sha(ledger) == plan['ledger_sha256'], 'P(z) ledger differs')
    return ledger, plan


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s16a_photoz_interval_rows.py', 'prepare_s23b_s16a_quality_supplement.py',
                                                 's16a_y1_source_contract.py'))]


def batch_name(name):
    return name.replace('photoz_interval_', 'local_background_')


def validate(directory=OUTPUT):
    plan = json.loads((directory/'plan.json').read_text())
    check_gate(plan)
    require((directory/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['cas_release_version'] == 'pdr1' and plan['chunk'] == CHUNK
            and plan['maximum_total_received_bytes'] == CAP and plan['output_columns'] == list(EXPRESSIONS)
            and plan['column_expressions'] == EXPRESSIONS and plan['schema_column_check'] == check_schema(schema_files()), 'plan scope differs')
    ledger, _ = pz_ledger()
    require(sha(ledger) == plan['ledger_sha256'], 'ID ledger differs')
    ids = ledger_ids(ledger)
    require(len(ids) == plan['n_ids'], 'ID count differs')
    expected = []
    for name, group in items(ids):
        name = batch_name(name)
        require((directory/name).read_text() == query(group), 'SQL scope differs '+name)
        expected.append(dict(name=name, n_ids=len(group), first_id=group[0], last_id=group[-1],
                             sql_bytes=(directory/name).stat().st_size, sha256=sha(directory/name)))
    require(plan['batches'] == expected, 'batch ledger differs')
    require(len(ids)*BYTES_PER_ROW < CAP, 'estimated bytes exceed the per-plan cap')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], ids=len(ids), batches=len(expected),
                columns=len(EXPRESSIONS), total_sql_bytes=sum(b['sql_bytes'] for b in expected), estimated_csv_bytes=len(ids)*BYTES_PER_ROW,
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(directory=OUTPUT):
    require(not directory.exists(), 'refusing to overwrite plan')
    files = schema_files()
    evidence = check_schema(files)
    ledger, pz_plan = pz_ledger()
    ids = ledger_ids(ledger)
    directory.mkdir(parents=True)
    batches = []
    for name, group in items(ids):
        name = batch_name(name)
        (directory/name).write_text(query(group))
        batches.append(dict(name=name, n_ids=len(group), first_id=group[0], last_id=group[-1],
                            sql_bytes=(directory/name).stat().st_size, sha256=sha(directory/name)))
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='Candidate (a) of the S16A central |e| deficit of blue primary sources at 30-60 arcsec: local background and '
                'neighbour/BCG light at the source position (no sky-level column exists), lens versus randoms.',
        n_ids=len(ids), ledger_sha256=sha(ledger), ledger_counts=pz_plan['ledger_counts'], ledger_source=str(ledger),
        id_rule='The ID ledger of the received P(z)-summary plan be6bd1a9, unchanged (selected, Ephor_AB-available sources of '
                'the 11 lens discs and the 76 new random discs).',
        chunk=CHUNK, batches=batches, output_columns=list(EXPRESSIONS), column_expressions=EXPRESSIONS,
        joined_tables={a: 's16a_wide.'+t for a, t in TABLES.items()}, schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(bytes=len(ids)*BYTES_PER_ROW, method=f'{BYTES_PER_ROW} bytes per row upper estimate'),
        not_requested=['images or pixel data', 'other bands than g and i', 'new objects or sky regions', 'shape-catalogue columns'],
        intake_checks=['Exactly one row per requested ID per batch; meas_object_id and forced_object_id equal the requested ID.'],
        analysis_plan=['Lens - random (group jackknife) in angular and proper-radius bins of: log i/g variance at matched i and '
                       'raw_sigma, countinputs, parent/child flux ratio, parent - child second moments and the induced '
                       'ellipticity of the extra light (tangential to the centre or not), deblending flags.',
                       'Random |e| versus each proxy and the proxy-matched lens - random |e| of blue primary sources at theta < 60 arcsec.',
                       'No selection change or science promotion without a separate user decision.'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256={str(p): sha(p) for p in (SOURCES/'sources.json', *files.values(), PZ_PLAN/'plan.json', ledger)},
        code_sha256={str(f): sha(f) for f in code_files()},
    )
    plan['query_plan_sha256'] = canonical_sha(plan)
    write_json(directory/'plan.json', plan)
    (directory/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  plan.json\n')
    return validate(directory)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Verify the existing local plan; no network')
    args = parser.parse_args()
    print(json.dumps(validate() if args.verify else build(), indent=2))
