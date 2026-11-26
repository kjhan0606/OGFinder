#!/usr/bin/env python3
"""Plan an ID-limited S16A Ephor_AB query of photo-z PDF summaries for the standard-boost sources; local only.

Scope: the selected, Ephor_AB-available sources (Y1 diagnostic + numeric calibration + finite photoz_best and
risk) of the 11 lens centres (pilot pairs) and of the 76 designed random centres (receipt c1dfe5a4 and its ID
supplements). No background cut is applied to the IDs, so any cut can be re-made from the received rows.
Columns: the PDF summaries of s16a_wide.photoz_ephor_ab not yet received - photoz_mc (one Monte-Carlo draw
from the PDF), mean, median, mode, std_best, conf_best and the 68/95 per cent interval bounds. Full PDFs are not
CAS columns and are not requested. Form: FROM VALUES ids LEFT JOIN, 5,000 IDs per batch, exactly one row per ID.
Never submits SQL.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from prepare_s16a_random_candidate_rows import check_gate
from prepare_s23b_s16a_quality_supplement import CAP, CHUNK, canonical_sha, require, sha, write_json
from s16a_y1_source_contract import object_id

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OUTPUT = DATA/'s16a_photoz_interval_rows_plan_2026-11-26_01'
SOURCES = DATA/'hsc_schema_tables_2026-11-24_01'
SCHEMA_FILE = 'pdr1.s16a_wide.photoz_ephor_ab.json'
SCHEMA = 'astrafek-s16a-photoz-interval-rows-plan-v1'
TABLE = 'photoz_ephor_ab'
FIELDS = ('object_id', 'photoz_mc', 'photoz_mean', 'photoz_median', 'photoz_mode', 'photoz_std_best', 'photoz_conf_best',
          'photoz_err68_min', 'photoz_err68_max', 'photoz_err95_min', 'photoz_err95_max')
COLUMNS = {'object_id': 'r.object_id', **{'ephor_ab_'+f: 'e.'+f for f in FIELDS}}
JOINED_IDS = ('ephor_ab_object_id',)
BYTES_PER_ROW = 150  # 19-digit ID twice and ten reals of up to 11 characters


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in COLUMNS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    return (f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n'
            f'LEFT JOIN s16a_wide.{TABLE} AS e ON e.object_id = r.object_id\nORDER BY r.object_id;\n')


def schema_file():
    manifest = json.loads((SOURCES/'sources.json').read_text())
    entry = next(d for d in manifest['documents'] if d['file'] == SCHEMA_FILE)
    require(entry['url'] == 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/'+SCHEMA_FILE and sha(SOURCES/SCHEMA_FILE) == entry['sha256'],
            'schema provenance/hash differs')
    return SOURCES/SCHEMA_FILE


def check_schema(file):
    doc = json.loads(file.read_text())
    require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', TABLE), 'schema identity differs')
    cols = {dict(zip(doc['headers'], r))['name']: dict(zip(doc['headers'], r)) for r in doc['columns']}
    evidence = []
    for f in FIELDS:
        require(f in cols and cols[f]['type'] == ('bigint' if f == 'object_id' else 'real'), 'unexpected schema for '+f)
        evidence.append(dict(alias='ephor_ab_'+f, type=cols[f]['type'], description=cols[f]['description']))
    require(cols['object_id']['constraint'] == 'primary', 'object_id is not the documented primary key')
    return evidence


def derive_ids():
    """Heavy: loads the verified receipts exactly as the standard boost does."""
    import audit_s16a_boost_mask_field_tests as bt
    from audit_s16a_random_extension_boosts import load_new, load_new_centres
    centres, selection_sha = load_new_centres()
    new, counts, provenance, _ = load_new(centres)
    pilot_centres, pilot, pilot_prov, _ = bt.load_centres()
    rows = []
    for c in pilot_centres:
        if c['role'] == 'lens':
            e = pilot[c['center_id']]
            rows += [(o, c['center_id']) for o, m in zip(e['ids'], e['missing']) if not m]
    for cid, e in new.items():
        rows += [(o, cid) for o, m in zip(e['ids'], e['missing']) if not m]
    return rows, dict(selection_summary_sha256=selection_sha, new_receipts=provenance, pilot_receipts=pilot_prov, counts=counts)


def write_ledger(path, rows):
    ids = sorted({int(o) for o, _ in rows})
    lens = {int(o) for o, c in rows if c.startswith('LC-')}
    with path.open('w', newline='') as stream:
        w = csv.writer(stream, delimiter='\t', lineterminator='\n')
        w.writerow(['object_id', 'in_lens_disc'])
        w.writerows([i, i in lens] for i in ids)
    return dict(n_ids=len(ids), n_lens_disc_ids=len(lens), n_pairs=len(rows))


def ledger_ids(path):
    with path.open() as stream:
        r = list(csv.DictReader(stream, delimiter='\t'))
    ids = [x['object_id'] for x in r]
    require(ids == sorted(ids, key=int) and len(set(ids)) == len(ids), 'ledger not sorted unique')
    return ids


def items(ids):
    return [(f'photoz_interval_{i:03d}.sql', ids[s:s+CHUNK]) for i, s in enumerate(range(0, len(ids), CHUNK))]


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s23b_s16a_quality_supplement.py', 's16a_y1_source_contract.py',
                                                 'audit_s16a_random_extension_boosts.py', 'audit_s16a_boost_mask_field_tests.py'))]


def validate(directory=OUTPUT):
    plan = json.loads((directory/'plan.json').read_text())
    check_gate(plan)
    require((directory/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['cas_release_version'] == 'pdr1' and plan['chunk'] == CHUNK
            and plan['maximum_total_received_bytes'] == CAP and plan['output_columns'] == list(COLUMNS)
            and plan['column_expressions'] == COLUMNS and plan['schema_column_check'] == check_schema(schema_file()), 'plan scope differs')
    ledger = directory/'object_ids.tsv'
    require(sha(ledger) == plan['ledger_sha256'], 'ID ledger differs')
    ids = ledger_ids(ledger)
    require(len(ids) == plan['n_ids'], 'ID count differs')
    expected = []
    for name, group in items(ids):
        require((directory/name).read_text() == query(group), 'SQL scope differs '+name)
        expected.append(dict(name=name, n_ids=len(group), first_id=group[0], last_id=group[-1],
                             sql_bytes=(directory/name).stat().st_size, sha256=sha(directory/name)))
    require(plan['batches'] == expected, 'batch ledger differs')
    require(len(ids)*BYTES_PER_ROW < CAP, 'estimated bytes exceed the per-plan cap')
    return dict(status='local_plan_verified_not_submitted', plan_sha256=plan['query_plan_sha256'], ids=len(ids),
                lens_disc_ids=plan['ledger_counts']['n_lens_disc_ids'], batches=len(expected),
                total_sql_bytes=sum(b['sql_bytes'] for b in expected), estimated_csv_bytes=len(ids)*BYTES_PER_ROW,
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(directory=OUTPUT):
    require(not directory.exists(), 'refusing to overwrite plan')
    evidence = check_schema(schema_file())
    rows, provenance = derive_ids()
    directory.mkdir(parents=True)
    counts = write_ledger(directory/'object_ids.tsv', rows)
    ids = ledger_ids(directory/'object_ids.tsv')
    batches = []
    for name, group in items(ids):
        (directory/name).write_text(query(group))
        batches.append(dict(name=name, n_ids=len(group), first_id=group[0], last_id=group[-1],
                            sql_bytes=(directory/name).stat().st_size, sha256=sha(directory/name)))
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='P(z)-based dilution check of the S16A standard boost (background cut photoz_best - z_lens > 0.2): the fraction of '
                'background-selected sources whose PDF draws or interval bounds fall at or below the lens redshift, lens versus randoms.',
        n_ids=len(ids), ledger_counts=counts, ledger_sha256=sha(directory/'object_ids.tsv'), id_provenance=provenance,
        id_rule='Unique object_id of selected, Ephor_AB-available sources in the 11 lens discs (pilot) and the 76 new random discs.',
        chunk=CHUNK, batches=batches, output_columns=list(COLUMNS), column_expressions=COLUMNS,
        joined_table='s16a_wide.'+TABLE, schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(bytes=len(ids)*BYTES_PER_ROW, method=f'{BYTES_PER_ROW} bytes per row upper estimate'),
        not_requested=['full photo-z PDFs (not CAS columns)', 'other photo-z codes', 'new objects or sky regions'],
        intake_checks=['Exactly one row per requested ID per batch; joined ephor_ab_object_id equals the requested ID.'],
        analysis_plan=['Per centre and annulus: fraction of background-cut sources with photoz_mc <= z_lens and with photoz_err95_min <= z_lens; '
                       'lens/random ratio and the implied member dilution of the standard boost.',
                       'Standard boost recomputed with photoz_err95_min > z_lens + 0.1 and with photoz_mc > z_lens + 0.2 as source cuts.',
                       'No selection change or science promotion without a separate user decision.'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256={str(p): sha(p) for p in (SOURCES/'sources.json', schema_file())},
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
