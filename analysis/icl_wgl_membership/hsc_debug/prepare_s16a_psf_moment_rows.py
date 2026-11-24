#!/usr/bin/env python3
"""Plan an S16A query of PSF-model adaptive moments for the lens-field sources; local files only.

The seeing-transfer test compared lens sources on the afterburner seeing scale with COSMOS
calibration objects on a PSF-moment proxy (2.3548 sqrt of the mean of shape_sdss_psf_11/_22,
S17A companion files). s16a_wide.forced carries the same shape_sdss_psf columns (arcsec^2), so
the same 115993 sources already holding afterburner seeing give a pointwise cross-calibration
and put the lenses on the COSMOS proxy scale. Reuses the executed ID ledger and batching.
Never submits SQL.
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
OUTPUT = DATA/'s16a_psf_moment_rows_plan_2026-11-24_01'
SCHEMA = 'astrafek-s16a-psf-moment-rows-plan-v1'
SCHEMA_URL = 'https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr1.s16a_wide.forced.json'
BANDS = 'grizy'
MOMENTS = ('11', '22', '12')
TABLES = {'f': 'forced'}
COLUMNS = {'object_id': 'r.object_id', 'forced_object_id': 'f.object_id'}
for _band in BANDS:
    for _m in MOMENTS:
        COLUMNS[f'{_band}shape_sdss_psf_{_m}'] = f'f.{_band}shape_sdss_psf_{_m}'
N_OBJECTS = 115993


def query(ids):
    ids = [object_id(i) for i in ids]
    require(ids and len(ids) == len(set(ids)), 'empty or duplicate SQL ID list')
    columns = ',\n  '.join(f'{expr} AS {alias}' for alias, expr in COLUMNS.items())
    values = ',\n  '.join(f'({i}::bigint)' for i in ids)
    joins = '\n'.join(f'LEFT JOIN s16a_wide.{table} AS {alias} ON {alias}.object_id = r.object_id'
                      for alias, table in TABLES.items())
    return f'SELECT\n  {columns}\nFROM (VALUES\n  {values}\n) AS r(object_id)\n{joins}\nORDER BY r.object_id;\n'


def schema_file():
    """The forced schema pinned by the executed quality supplement (same bytes as the 2026-11-24 public copy)."""
    manifest = json.loads((SUPPLEMENT/'sources.json').read_text())
    entry, = [d for d in manifest['documents'] if d['file'] == 'schema_forced.json']
    file = SUPPLEMENT/entry['file']
    require(entry['url'] == SCHEMA_URL and sha(file) == entry['sha256'], 'forced schema provenance/hash differs')
    return file


def check_schema(file):
    doc = json.loads(file.read_text())
    require((doc['db'], doc['schema'], doc['table']) == ('pdr1', 's16a_wide', 'forced'), 'schema identity differs')
    columns = {}
    for row in doc['columns']:
        c = dict(zip(doc['headers'], row))
        for band in BANDS:
            columns[c['name'].replace('(g|r|i|z|y)', band)] = c
    evidence = []
    for alias, expression in COLUMNS.items():
        prefix, column = expression.split('.')
        if prefix == 'r':
            continue
        info = columns.get(column)
        require(info is not None, 'undocumented column '+expression)
        want = ('bigint', '') if column == 'object_id' else ('real', 'arcsec^2')
        require((info['type'], info['unit']) == want, 'unexpected schema type/unit '+expression)
        evidence.append(dict(alias=alias, expression=expression, table='s16a_wide.forced',
                             type=info['type'], unit=info['unit'], description=info['description']))
    return evidence


def batches():
    ids, supplement = ledger_ids()
    return ids, [(name.replace('s16a_afterburner_photoz_rows_', 's16a_psf_moment_rows_'), group)
                 for name, group in batches_for(ids, supplement)]


def inputs(file):
    paths = [SUPPLEMENT/'quality_supplement_plan.json', SUPPLEMENT/'object_ids.tsv', SUPPLEMENT/'sources.json', file]
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
    file = schema_file()
    require(plan['input_sha256'] == inputs(file) and plan['schema_column_check'] == check_schema(file), 'schema evidence differs')
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
    file = schema_file()
    evidence = check_schema(file)
    ids, groups = batches()
    output.mkdir(parents=True)
    ledger = []
    for name, group in groups:
        sql = output/name
        sql.write_text(query(group))
        ledger.append(dict(name=name, first_object_id=group[0], last_object_id=group[-1],
                           expected_rows=len(group), sql_bytes=sql.stat().st_size, sha256=sha(sql)))
    per_row = dict(low=150, high=300)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', source_catalog='HSC S16A Wide',
        purpose='Cross-calibrate the afterburner seeing (parent_flux_convolved_seeing) against the PSF-model adaptive-moment '
                'FWHM used as the COSMOS seeing proxy, on the same sources, and place the lens sources on the COSMOS proxy scale.',
        id_scope_plan_sha256=SUPPLEMENT_SHA, id_scope='Exact 115993 object IDs and 24 batch boundaries of the executed quality supplement.',
        n_unique_objects=len(ids), n_new_objects=0, expected_output_rows=len(ids),
        batch_size=CHUNK, n_query_batches=len(ledger), batches=ledger,
        output_columns=list(COLUMNS), column_expressions=COLUMNS,
        joined_tables={k: 's16a_wide.'+v for k, v in TABLES.items()},
        schema_column_check=evidence, maximum_total_received_bytes=CAP,
        estimated_csv_bytes=dict(low=len(ids)*per_row['low'], high=len(ids)*per_row['high'],
                                 method=f'{per_row["low"]}-{per_row["high"]} bytes/row for {len(COLUMNS)} columns (2 IDs, 15 reals); the executed '
                                        '13-column cmodel photo-z query was 15432619 bytes (133 bytes/row).'),
        included='Row-existence ID and grizy shape_sdss_psf_11/_22/_12 (adaptive moments of the PSF model at the object, arcsec^2) from s16a_wide.forced.',
        excluded=['new objects or sky regions', 'object shapes, fluxes, positions and flags', 'HSM PSF moments (meas2)',
                  'PSF images or models', 'photo-z', 'S17A/S19A/S23B tables'],
        join_policy='VALUES exact signed-bigint IDs, LEFT JOIN by object_id, no filtering. '
                    'Expect exactly one row per requested ID; abort on duplicate, missing or extra IDs or on a joined ID that differs.',
        server_sql_validation='Not submitted. Column names, types and units checked against the public forced schema (pinned by the executed '
                              'quality supplement; byte-identical to the 2026-11-24 public copy). The schema lists no index, but the executed '
                              'afterburner query already LEFT JOINed s16a_wide.forced by object_id in the same 24 batches and completed.',
        analysis_plan=['Per source and band: afterburner FWHM = seeing x 2.3548 x 0.168 arcsec; trace FWHM = 2.3548 sqrt((I11+I22)/2) '
                       '(the COSMOS proxy); determinant FWHM = 2.3548 (I11 I22 - I12^2)^(1/4). Report the median ratio and scatter in FWHM bins, '
                       'and which definition the afterburner seeing follows.',
                       'Map the 1.092-arcsec deconvolution threshold (Ephor_AB missing) onto the trace-FWHM scale.',
                       'Recompute each lens\'s eff-weighted r FWHM on the trace scale and its two-sided support relative to the COSMOS '
                       'weighted q01/q99 proxy; report which lenses change support group. S16A vs S17A pipeline PSF models may differ; '
                       'state this as a residual.',
                       'No selection change, correction or science promotion.'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        input_sha256=inputs(file), code_sha256={str(f): sha(f) for f in code_files()},
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
