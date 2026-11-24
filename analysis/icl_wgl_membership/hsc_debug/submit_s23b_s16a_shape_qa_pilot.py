#!/usr/bin/env python3
"""Receive only the exact approved S16A diagnostic pilot, serially and resumably.

Reuses the established HSC CAS client. Never changes the SQL, expands the
sample, writes credentials, or performs a scientific lensing measurement.
The byte cap covers downloaded CSV payload, including interrupted attempts.
"""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
import getpass
import json
import math
from pathlib import Path
import platform
import re
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from prepare_s23b_s16a_shape_count_plan import canonical_sha
from prepare_s23b_s16a_shape_qa_pilot import (
    OUTPUT as DEFAULT_PLAN, OUTPUT_COLUMNS, COLUMN_EXPRESSIONS, query,
    AUDIT, AUDIT_SHA, TARGET_MANIFEST, pilot_lenses, make_centers,
)
from audit_s23b_s16a_wl_random_support import read_tsv, require, angular_separation_arcmin
from prepare_s23b_s16a_shape_qa_schema_fix import same_scope, OLD as OLD_WEIGHT, NEW as NEW_WEIGHT, OUTPUT as CORRECTED_PLAN

APPROVED_SHA = '7af70c2be81f83e61bc6f74caf98b66a6309ae1b0caed0b9b7d1a32f42bd3f10'
IDENTIFIER_FIX_SHA = '1e3ac9451e519761e0f773915de79728d5443808efedc74155963a79c3ae4b30'
MAX_BYTES = 134217728
CHUNK = 1024*1024
SCHEMA = 'astrafek-s16a-shape-qa-pilot-receipt-v1'
FLOAT_COLUMNS = OUTPUT_COLUMNS[2:]


def validate_plan(directory):
    directory = Path(directory).resolve()
    plan = json.loads((directory/'shape_qa_plan.json').read_text())
    executed_hash = plan.get('query_plan_sha256')
    require(canonical_sha(plan) == executed_hash and executed_hash in (APPROVED_SHA, IDENTIFIER_FIX_SHA), 'plan differs from exact user-approved scope')
    require((directory/'PLAN_SHA256').read_text().split()[0] == executed_hash, 'plan hash sidecar differs')
    corrected = executed_hash == IDENTIFIER_FIX_SHA
    expected_expressions = dict(COLUMN_EXPRESSIONS)
    if corrected:
        original = validate_plan(DEFAULT_PLAN)
        same_scope(original, plan)
        amendment = plan['identifier_correction']
        require(amendment['predecessor_plan_sha256'] == APPROVED_SHA
                and amendment['semantic_data_scope_changed'] is False
                and amendment['original_expression'] == OLD_WEIGHT
                and amendment['corrected_expression'] == NEW_WEIGHT, 'unreviewed identifier amendment')
        require(sha(Path(__file__).with_name('prepare_s23b_s16a_shape_qa_schema_fix.py')) == plan['schema_fix_script_sha256'], 'schema-fix generator differs')
        predecessor = json.loads((DEFAULT_PLAN/'results/shape_qa_status.json').read_text())
        require(predecessor['confirmed_bytes_received'] == predecessor['budget_bytes_charged'] == 0
                and predecessor['query_submitted'] is False
                and all('job_id' not in j for j in predecessor['jobs'].values()), 'predecessor has downloads/jobs; review cumulative scope')
        expected_expressions['weight'] = NEW_WEIGHT
    require(plan['schema'] == 'astrafek-s16a-diagnostic-shape-qa-pilot-plan-v1'
            and plan['status'] == 'prepared_not_submitted'
            and plan['execution_requires_explicit_user_approval'] is True, 'invalid plan approval gate')
    for flag in ('query_submitted','network_used','source_rows_downloaded','images_returned','psfs_returned','photoz_returned','science_measurement'):
        require(plan[flag] is False, 'immutable plan scope flags differ')
    require(plan['cas_release_version'] == 'pdr1' and plan['source_catalog'] == 'HSC S16A Wide', 'release differs')
    require(plan['output_columns'] == OUTPUT_COLUMNS and plan['column_expressions'] == expected_expressions
            and len(OUTPUT_COLUMNS) == 12, 'output columns differ')
    require(plan['maximum_total_received_bytes'] == MAX_BYTES
            and plan['expected_source_pairs_from_prior_counts'] == 169535
            and plan['n_lenses'] == 11 and plan['n_centers_in_ledger'] == 66
            and plan['n_query_batches'] == 5 and plan['batch_size'] == 10
            and plan['n_positive_centers_to_query'] == 43
            and plan['n_known_zero_centers_not_queried'] == 23, 'approved row/byte/center budget differs')
    for path, digest in plan['input_sha256'].items():
        require(sha(Path(path)) == digest, 'pinned plan input differs')
    require(sha(AUDIT/'summary.json') == AUDIT_SHA, 'support audit differs')
    audit = json.loads((AUDIT/'summary.json').read_text())
    for name in ('per_lens_support.tsv','per_random_center.tsv'):
        require(sha(AUDIT/name) == audit['output_sha256'][name], 'support ledger differs')
    require(sha(Path(__file__).with_name('prepare_s23b_s16a_shape_qa_pilot.py')) == plan['script_sha256'], 'plan builder differs')
    targets = {str(r['id']):r for r in json.loads(TARGET_MANIFEST.read_text())['clusters']}
    expected_centers = make_centers(pilot_lenses(read_tsv(AUDIT/'per_lens_support.tsv')),
                                    read_tsv(AUDIT/'per_random_center.tsv'), targets)
    require(plan['centers'] == expected_centers, 'center identities, counts or coordinates differ')
    for name, digest in plan['output_sha256'].items():
        require(Path(name).name == name and sha(directory/name) == digest, 'plan member differs')
    positive = [r for r in expected_centers if r['query_planned']]
    require(len(plan['batches']) == 5 and sum(r['expected_raw_pairs'] for r in positive) == 169535, 'batch total differs')
    for i, batch in enumerate(plan['batches']):
        group = positive[10*i:10*(i+1)]
        name = f's16a_shape_qa_{i:03d}.sql'
        file = directory/name
        sql = '\nUNION ALL\n'.join(query(r) for r in group)+';\n'
        if corrected:
            sql = sql.replace(OLD_WEIGHT, NEW_WEIGHT)
        require(batch['name'] == name and batch['center_ids'] == [r['center_id'] for r in group]
                and batch['expected_raw_pairs'] == sum(r['expected_raw_pairs'] for r in group)
                and file.read_text() == sql and sha(file) == batch['sha256']
                and file.stat().st_size == batch['sql_bytes'], 'SQL scope differs')
    return plan


def numeric(value):
    # Keep missing values visible to QA. They are never repaired to zero.
    if value.strip().lower() in ('', 'null', 'none', 'nan'):
        return float('nan')
    return float(value)


def validate_csv(path, centers):
    expected = {r['center_id']:r for r in centers}
    require(len(expected) == len(centers) and all(r['expected_raw_pairs'] > 0 for r in centers), 'invalid expected CSV centers')
    counts, seen = Counter(), set()
    nonfinite = Counter()
    maximum_offset_excess_arcsec = 0.
    with Path(path).open(newline='', encoding='utf-8-sig') as f:
        reader = catalogue_reader(f)
        require(reader.fieldnames == OUTPUT_COLUMNS, 'CAS columns differ from approved twelve fields')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), 'malformed CAS row')
            ident, obj = row['center_id'], row['object_id']
            require(ident in expected and re.fullmatch(r'[0-9]+', obj) is not None and int(obj) > 0, 'unexpected center or object ID')
            key = (ident, int(obj))
            require(key not in seen, 'duplicate (center, object) pair')
            seen.add(key); counts[ident] += 1
            require(counts[ident] <= expected[ident]['expected_raw_pairs'], 'center exceeds approved expected row count')
            values = {k:numeric(row[k]) for k in FLOAT_COLUMNS}
            ra, dec = values['source_ra'], values['source_dec']
            require(math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90, 'invalid raw sky coordinate; units not accepted')
            c = expected[ident]
            sep = angular_separation_arcmin(c['ra_deg'], c['dec_deg'], ra, dec)*60
            excess = sep-c['radius_arcsec']
            require(excess <= .01, 'raw coordinates do not match degree-valued query cone; no conversion attempted')
            maximum_offset_excess_arcsec = max(maximum_offset_excess_arcsec, excess)
            nonfinite.update(k for k,v in values.items() if not math.isfinite(v))
    require(dict(counts) == {k:r['expected_raw_pairs'] for k,r in expected.items()}, 'center counts differ from approved prior count ledger')
    return dict(rows=sum(counts.values()), per_center_counts=dict(counts),
                nonfinite_by_column=dict(nonfinite), coordinates_pass_degree_cones=True,
                maximum_positive_cone_excess_arcsec=maximum_offset_excess_arcsec)


def budgeted_download(response, destination, state, checkpoint):
    remaining = MAX_BYTES-state['budget_bytes_charged']
    require(remaining > 0, 'approved total CSV byte budget exhausted')
    length = response.headers.get('Content-Length')
    if length is not None:
        length = int(length)
        require(0 <= length <= remaining, 'CAS payload would exceed remaining approved byte budget')
    received = 0
    with destination.open('xb') as stream:
        while True:
            if length is not None and received == length:
                break
            remaining = MAX_BYTES-state['budget_bytes_charged']
            require(remaining > 0, 'approved total CSV byte budget exhausted; partial file preserved')
            amount = min(CHUNK, remaining, length-received if length is not None else CHUNK)
            # Reserve before I/O: interruption can overcharge, never undercount.
            state['budget_bytes_charged'] += amount
            save(checkpoint, state)
            block = response.read(amount)
            require(len(block) <= amount, 'response exceeded requested block size')
            state['budget_bytes_charged'] -= amount-len(block)
            state['confirmed_bytes_received'] += len(block)
            if block:
                state['source_rows_downloaded'] = True
            save(checkpoint, state)
            if not block:
                break
            stream.write(block)
            received += len(block)
    require(length is None or received == length, 'truncated CAS payload; partial file preserved')
    return received


def fetch_batch(client, batch, sql, centers, results, state, checkpoint, sleep=time.sleep):
    name = Path(batch['name']).stem
    item = state['jobs'].setdefault(name, dict(status='pending'))
    target = results/(name+'.csv')
    if item['status'] == 'complete':
        require(target.is_file() and sha(target) == item['sha256'], 'completed CSV missing or changed')
        validate_csv(target, centers)
        return
    require(item['status'] not in ('error','downloaded_invalid'), 'prior CAS/data error requires review')
    require(not (item['status'] == 'submitting' and 'job_id' not in item), 'submission outcome unknown; inspect CAS before any retry')
    if 'job_id' not in item:
        item['status'] = 'submitting'; save(checkpoint, state)
        try:
            submitted = client.json('submit', catalog_job=dict(sql=sql, out_format='csv',
                include_metainfo_to_body=False, release_version='pdr1'), nomail=True, skip_syntax_check=False)
        except urllib.error.HTTPError as error:
            if error.code in (400,406):
                message = error.read(8192).decode('utf-8', errors='replace')
                for secret in client.credential.values():
                    if secret:
                        message = message.replace(secret, '[redacted]').replace(secret.replace('/', '\\/'), '[redacted]')
                item['redacted_schema_error'] = message[:2000]
                item['http_code'] = error.code
            if error.code in (400,401,403,406,429):
                item['status'] = 'pending'; save(checkpoint, state)
            raise
        item.update(status='submitted', job_id=str(submitted['id']))
        state['query_submitted'] = True
        save(checkpoint, state)
    while True:
        result = client.json('status', id=item['job_id'])
        if result.get('status') == 'done':
            break
        if result.get('status') == 'error':
            message = str(result.get('error','CAS query failed'))
            for secret in client.credential.values():
                if secret: message = message.replace(secret, '[redacted]')
            item.update(status='error', query_error=message[:2000]); save(checkpoint, state)
            raise ValueError('CAS query failed; redacted error preserved in checkpoint')
        print(f"{name}: {result.get('status','unknown')}", flush=True)
        sleep(45)
    partial = results/(name+'.csv.part')
    if partial.exists():
        partial.rename(results/(name+f'.csv.incomplete-{time.time_ns()}'))
    with client.post('download', id=item['job_id']) as response:
        budgeted_download(response, partial, state, checkpoint)
    try:
        metadata = validate_csv(partial, centers)
    except (ValueError, TypeError, OverflowError):
        item.update(status='downloaded_invalid', partial_file=partial.name, sha256=sha(partial))
        save(checkpoint, state)
        raise
    partial.replace(target)
    item.update(status='complete', bytes=target.stat().st_size, sha256=sha(target), **metadata)
    save(checkpoint, state)
    print(f"{name}: verified {metadata['rows']} rows, {item['bytes']} bytes", flush=True)


def verify_complete(plan, results, state):
    centers = {r['center_id']:r for r in plan['centers']}
    rows, payload_bytes = 0, 0
    files = {}
    require(set(state['jobs']) == {Path(b['name']).stem for b in plan['batches']}, 'completed job set differs')
    for batch in plan['batches']:
        name = Path(batch['name']).stem; item = state['jobs'][name]
        file = results/(name+'.csv')
        require(item['status'] == 'complete' and sha(file) == item['sha256'] and file.stat().st_size == item['bytes'], 'completed download differs')
        valid = validate_csv(file, [centers[i] for i in batch['center_ids']])
        require(valid['rows'] == batch['expected_raw_pairs'] == item['rows'], 'final row count differs')
        rows += valid['rows']; payload_bytes += item['bytes']; files[file.name] = sha(file)
    require(rows == 169535 and payload_bytes <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= MAX_BYTES, 'row or byte budget differs')
    return dict(schema=SCHEMA, status='complete_diagnostic_shape_download',
        approved_scope_plan_sha256=APPROVED_SHA, executed_plan_sha256=plan['query_plan_sha256'],
        n_csv_files=5, n_queried_centers=43,
        n_lenses=11, n_random_centers_in_ledger=55, known_zero_centers_not_queried=23,
        output_columns=OUTPUT_COLUMNS, n_source_pair_rows=rows, verified_csv_bytes=payload_bytes,
        confirmed_bytes_received=state['confirmed_bytes_received'], budget_bytes_charged=state['budget_bytes_charged'],
        maximum_total_received_bytes=MAX_BYTES, output_sha256=files,
        credentials_saved=False, source_rows_downloaded=True, images_returned=False,
        psfs_returned=False, photoz_returned=False, science_measurement=False,
        submitter_sha256=sha(Path(__file__)),
        limitations=['Payload cap includes interrupted attempts; reservations surviving interruption are conservatively charged.',
            'Shape/calibration and radial support QA only; no redshift or complete selection/PSF fields were received.'])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    mode = p.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--submit', action='store_true')
    p.add_argument('--plan-dir', type=Path, default=CORRECTED_PLAN)
    p.add_argument('--user', default='kjhan0606')
    p.add_argument('--approved-plan-sha256')
    args = p.parse_args()
    plan_dir = args.plan_dir.resolve(); plan = validate_plan(plan_dir)
    print(json.dumps(dict(plan_sha256=plan['query_plan_sha256'], approved_scope_plan_sha256=APPROVED_SHA, centers=43, columns=OUTPUT_COLUMNS,
        expected_rows=169535, maximum_total_received_bytes=MAX_BYTES, batches=5,
        network_used=False if args.plan else 'only explicitly approved plan')), flush=True)
    if args.plan: return 0
    require(args.approved_plan_sha256 == APPROVED_SHA, '--submit requires the exact approved plan hash')
    require(plan['query_plan_sha256'] == IDENTIFIER_FIX_SHA,
            'original plan has the rejected weight identifier; use the corrected _02 plan')
    require(platform.node().lower().startswith('lageunha'), 'run the established CPU workflow on LagEunha')
    results = plan_dir/'results'; checkpoint = results/'shape_qa_status.json'
    if results.exists():
        require(checkpoint.is_file(), 'refusing existing results without a checkpoint')
    else:
        results.mkdir(exist_ok=False)
    with (results/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if checkpoint.exists():
            state = json.loads(checkpoint.read_text())
            require(state['schema'] == SCHEMA and state['plan_sha256'] == plan['query_plan_sha256']
                    and state.get('approved_scope_plan_sha256', APPROVED_SHA) == APPROVED_SHA, 'checkpoint scope differs')
            require(0 <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= MAX_BYTES, 'checkpoint byte budget invalid')
        else:
            state = dict(schema=SCHEMA, status='pending', plan_sha256=plan['query_plan_sha256'], approved_scope_plan_sha256=APPROVED_SHA,
                jobs={}, query_submitted=False, confirmed_bytes_received=0, budget_bytes_charged=0,
                credentials_saved=False, source_rows_downloaded=False, science_measurement=False)
            save(checkpoint, state)
        if state['status'] == 'complete':
            verify_complete(plan, results, state)
            require(sha(results/'download_summary.json') == state['summary_sha256'], 'download summary differs')
            print('Already complete; locally verified without network or credentials.')
            return 0
        require(sys.stdin.isatty(), 'interactive terminal required for hidden credential input')
        password = getpass.getpass('HSC password (hidden): ')
        require(bool(password), 'empty password')
        client = Cas(args.user, password)
        state.update(status='running', last_error_type=None)
        save(checkpoint, state)
        try:
            centers = {r['center_id']:r for r in plan['centers']}
            for batch in plan['batches']:
                state['current_query'] = Path(batch['name']).stem; save(checkpoint, state)
                fetch_batch(client, batch, (plan_dir/batch['name']).read_text(),
                    [centers[i] for i in batch['center_ids']], results, state, checkpoint)
            summary = verify_complete(plan, results, state)
            file = results/'download_summary.json'
            if file.exists():
                require(json.loads(file.read_text()) == summary, 'existing download summary differs')
            else:
                save(file, summary)
            state.update(status='complete', current_query=None, source_rows_downloaded=True,
                         summary_sha256=sha(file), total_rows=169535)
            save(checkpoint, state)
        except BaseException as error:
            state.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'stopped_on_error',
                         last_error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear(); password = None
    print('Complete: 5 CSVs, 169535 source-pair rows; science_measurement=false.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f'HSC CAS HTTP {error.code}; no response body or credentials printed.', file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print('Stopped; checkpoint and partial data preserved.', file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print('S16A pilot: '+(str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
