#!/usr/bin/env python3
"""Plan the ID-limited S16A supplements for the sources of the 76 random centres; local only.

The IDs are the 373,449 unique source objects of the verified pair receipt (plan c1dfe5a4, approved and
received 2026-11-26). Three kinds, each the executed pilot query form unchanged (FROM VALUES ... LEFT JOIN):
- quality: the Y1 input contract of the quality supplement (plan 44b9d4df), 31 columns;
- cmodel: photoz_ephor and photoz_mizuki (plan 35265ee4);
- flags: forced/meas bright-object, edge and FDFC flags (plan e73c0f90).
Afterburner fluxes (plan 951a81fb) are a missingness diagnostic and are not requested. The sorted IDs are cut
into three contiguous parts, one plan each, so every plan stays under the 128 MiB per-plan budget.
Never submits SQL.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import prepare_s16a_bright_object_flag_rows as pflag
import prepare_s16a_cmodel_photoz_rows as pcm
import prepare_s23b_s16a_quality_supplement as pq
from prepare_s16a_random_candidate_rows import check_gate
from prepare_s16a_random_source_rows import OUTPUT as SOURCE_PLAN
from prepare_s23b_s16a_quality_supplement import CAP, CHUNK, canonical_sha, require, sha, write_json
from s16a_y1_source_contract import EXPRESSIONS

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SCHEMA = 'astrafek-s16a-random-id-supplement-plan-v1'
SOURCE_SHA = 'c1dfe5a4860ceb830281e1444b90b438a99e93702726490518fa32a90073b211'
N_OBJECTS, PARTS = 373449, 3
KINDS = {'quality': (pq.query, list(EXPRESSIONS), 391), 'cmodel': (pcm.query, list(pcm.COLUMNS), 133),
         'flags': (pflag.query, list(pflag.COLUMNS), 168)}  # bytes per row from the pilot receipts


def output(part):
    require(part in range(1, PARTS+1), 'invalid part')
    return DATA/f's16a_random_id_supplement_plan_2026-11-26_{part:02d}'


def receipt():
    results = SOURCE_PLAN/'results'
    status = json.loads((results/'status.json').read_text())
    summary = json.loads((results/'download_summary.json').read_text())
    require(status['status'] == 'complete' and status['plan_sha256'] == SOURCE_SHA == summary['approved_scope_plan_sha256']
            and sha(results/'download_summary.json') == status['summary_sha256'], 'source-row receipt incomplete or changed')
    return results, summary


def all_ids():
    results, summary = receipt()
    ids = set()
    for name, digest in sorted(summary['output_sha256'].items()):
        if not name.startswith('s16a_random_pairs_'):
            continue
        file = results/name
        require(sha(file) == digest, 'pair CSV changed: '+name)
        with file.open() as stream:
            stream.readline()
            for line in stream:
                ids.add(int(line.split(',', 2)[1]))
    require(len(ids) == N_OBJECTS, 'unique source count differs')
    return sorted(ids)


def part_ids(part, ids=None):
    ids = all_ids() if ids is None else ids
    size = -(-len(ids)//PARTS)
    return [str(i) for i in ids[(part-1)*size:part*size]]


def ledger(part, ids=None):
    chosen = part_ids(part, ids)
    out = []
    for kind, (query, _, _) in KINDS.items():
        for i, start in enumerate(range(0, len(chosen), CHUNK)):
            group = chosen[start:start+CHUNK]
            out.append((f'{kind}_{i:03d}.sql', kind, group, query(group)))
    return chosen, out


def evidence():
    return dict(quality=pq.check_schemas(pq.SOURCES), cmodel=pcm.check_schemas(pcm.schema_files()),
                flags=pflag.check_schema(pflag.schema_files()))


def inputs():
    results, _ = receipt()
    files = [results/'status.json', results/'download_summary.json', SOURCE_PLAN/'PLAN_SHA256', pq.SOURCES/'sources.json',
             pcm.SOURCES/'sources.json', *pcm.schema_files().values(), *pflag.schema_files().values()]
    return {str(p): sha(p) for p in files}


def code_files():
    here = Path(__file__).resolve()
    return [here, *(here.with_name(n) for n in ('prepare_s23b_s16a_quality_supplement.py', 'prepare_s16a_cmodel_photoz_rows.py',
                                                 'prepare_s16a_bright_object_flag_rows.py', 'prepare_s16a_afterburner_photoz_rows.py',
                                                 's16a_y1_source_contract.py', 'prepare_s16a_random_source_rows.py'))]


def batches(directory, items):
    return [dict(name=n, kind=k, n_ids=len(g), first_id=g[0], last_id=g[-1], sql_bytes=(directory/n).stat().st_size,
                 sha256=sha(directory/n)) for n, k, g, _ in items]


def estimate(n):
    return dict(rows_per_kind=n, csv_bytes=sum(n*b for _, _, b in KINDS.values()))


def validate(directory, part):
    plan = json.loads((directory/'plan.json').read_text())
    check_gate(plan)
    require((directory/'PLAN_SHA256').read_text().split()[0] == plan['query_plan_sha256'], 'sidecar differs')
    require(plan['schema'] == SCHEMA and plan['part'] == part and plan['parts'] == PARTS and plan['chunk'] == CHUNK
            and plan['maximum_total_received_bytes'] == CAP and plan['source_plan_sha256'] == SOURCE_SHA
            and plan['output_columns'] == {k: c for k, (_, c, _) in KINDS.items()}
            and plan['input_sha256'] == inputs() and plan['schema_column_check'] == evidence(), 'plan scope differs')
    chosen, items = ledger(part)
    for n, _, _, sql in items:
        require((directory/n).read_text() == sql, 'SQL scope differs '+n)
    require(plan['batches'] == batches(directory, items) and plan['n_ids'] == len(chosen)
            and plan['id_range'] == [chosen[0], chosen[-1]], 'batch ledger differs')
    require(estimate(len(chosen))['csv_bytes'] < CAP, 'estimated bytes exceed the per-plan cap')
    return dict(status='local_plan_verified_not_submitted', part=part, plan_sha256=plan['query_plan_sha256'], ids=len(chosen),
                batches=len(items), total_sql_bytes=sum(b['sql_bytes'] for b in plan['batches']), estimate=estimate(len(chosen)),
                maximum_csv_bytes=CAP, remote_query_submitted=False, science_measurement=False)


def build(part):
    directory = output(part)
    require(not directory.exists(), 'refusing to overwrite plan')
    chosen, items = ledger(part)
    directory.mkdir(parents=True)
    for n, _, _, sql in items:
        (directory/n).write_text(sql)
    plan = dict(
        schema=SCHEMA, status='prepared_not_submitted', execution_requires_explicit_user_approval=True,
        query_submitted=False, new_catalog_rows_received=False, science_measurement=False, calibration_applied=False,
        cas_release_version='pdr1', part=part, parts=PARTS, chunk=CHUNK, source_plan_sha256=SOURCE_SHA,
        purpose=f'Part {part}/{PARTS} of the ID-limited quality, cmodel photo-z and flag supplements for the sources of the 76 '
                'random centres, in the executed pilot query forms; the user approved making these plans on 2026-11-26.',
        n_ids=len(chosen), id_range=[chosen[0], chosen[-1]], id_rule='sorted unique object_id of the verified pair receipt, '
        f'contiguous part {part} of {PARTS} (ceil division), {CHUNK} IDs per batch',
        output_columns={k: c for k, (_, c, _) in KINDS.items()}, schema_column_check=evidence(),
        maximum_total_received_bytes=CAP,
        estimated=dict(**estimate(len(chosen)), method='bytes per row of the pilot receipts (quality 391, cmodel 133, flags 168)'),
        not_requested=['afterburner fluxes (diagnostic only; plan 951a81fb)'],
        intake_checks=['Exactly one row per requested ID in every batch, joined IDs equal the requested ID (pilot validators).'],
        before_submission=['User approval of this exact plan hash.',
                           'The first batch is the SQL compatibility check; stop on any error, no automatic scope change.'],
        batches=batches(directory, items), input_sha256=inputs(), code_sha256={str(f): sha(f) for f in code_files()},
    )
    plan['query_plan_sha256'] = canonical_sha(plan)
    write_json(directory/'plan.json', plan)
    (directory/'PLAN_SHA256').write_text(plan['query_plan_sha256']+'  plan.json\n')
    return validate(directory, part)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--part', type=int, required=True)
    parser.add_argument('--verify', action='store_true', help='Verify an existing local plan; no network')
    args = parser.parse_args()
    print(json.dumps(validate(output(args.part), args.part) if args.verify else build(args.part), indent=2))
