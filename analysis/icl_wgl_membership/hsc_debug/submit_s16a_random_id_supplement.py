#!/usr/bin/env python3
"""Receive only an approved part of the S16A ID supplements for the 76 random centres' sources, serially.

Same safeguards as the other S16A receivers: exact plan hash per part, pinned plan directory, cumulative byte
budget per plan, SQL hash re-check at use, no stored credentials (--use-env-auth or one hidden prompt), transport
framing and final line break. Every batch must return exactly one row per requested ID (the pilot validators).
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import json
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, save, sha
from submit_s23b_s16a_shape_qa_pilot import budgeted_download, MAX_BYTES
from submit_s16a_afterburner_photoz_rows import batch_ids, batch_text, code_hashes, credentials
from cas_transfer_checks import redact, require_line_end, require_poll_budget, transport
import submit_s16a_bright_object_flag_rows as sflag
import submit_s16a_cmodel_photoz_rows as scm
import submit_s23b_s16a_quality_supplement as sq
from prepare_s16a_random_id_supplement import PARTS, output, validate
from prepare_s23b_s16a_quality_supplement import CAP, require

APPROVED_SHA = {1: 'aeda229e886582c75074c3c58e19bbb727a364fc5dff3fbc44463893cd801fe7',  # user approved 2026-11-26
                2: 'bce80af866047312eeb1d8327a3d6378cc5366ff45a10a13575296f22f73f022',
                3: '7d7ed612bea5678ec115e9d8687a6ce997fd28826c177dbd97b9ef3843ee2cca'}
SCHEMA = 'astrafek-s16a-random-id-supplement-receipt-v1'
POLL_SECONDS = 20
VALIDATORS = dict(quality=sq.validate_csv, cmodel=scm.validate_csv, flags=sflag.validate_csv)


def load_plan(directory, part):
    checked = validate(directory, part)
    require(APPROVED_SHA.get(part) is not None and checked['plan_sha256'] == APPROVED_SHA[part], 'not the exact approved ID supplement plan')
    require(MAX_BYTES == CAP == 134217728, 'shared cumulative download cap differs')
    return json.loads((directory/'plan.json').read_text())


def ids_of(directory, batch):
    ids = batch_ids(directory, batch)
    require(len(ids) == batch['n_ids'] and ids[0] == batch['first_id'] and ids[-1] == batch['last_id'], 'batch IDs differ from the ledger')
    return ids


def fetch_batch(client, batch, sql, ids, results, state, checkpoint, sleep=time.sleep):
    name = Path(batch['name']).stem
    check = lambda file, _: VALIDATORS[batch['kind']](file, ids)
    item = state['jobs'].setdefault(name, dict(status='pending'))
    target = results/(name+'.csv')
    if item['status'] == 'complete':
        require(target.is_file() and sha(target) == item['sha256'] and target.stat().st_size == item['bytes'], 'completed CSV missing or changed')
        check(target, batch)
        return
    require(item['status'] not in ('error', 'downloaded_invalid'), 'prior CAS/data error requires review')
    require(not (item['status'] == 'submitting' and 'job_id' not in item), 'submission outcome unknown; inspect CAS before retry')
    if 'job_id' not in item:
        item['status'] = 'submitting'; save(checkpoint, state)
        try:
            submitted = client.json('submit', catalog_job=dict(sql=sql, out_format='csv',
                include_metainfo_to_body=False, release_version='pdr1'), nomail=True, skip_syntax_check=False)
        except urllib.error.HTTPError as error:
            if error.code in (400, 406):
                item['redacted_schema_error'] = redact(error.read(8192).decode('utf-8', errors='replace'), client)
                item['http_code'] = error.code
            if error.code in (400, 401, 403, 406, 429):
                item['status'] = 'pending'; save(checkpoint, state)
            raise
        item.update(status='submitted', job_id=str(submitted['id']))
        state['query_submitted'] = True
        save(checkpoint, state)
    waited, started = 0, time.monotonic()
    while True:
        require_poll_budget(max(waited, time.monotonic()-started))
        status = client.json('status', id=item['job_id'])
        if status.get('status') == 'done':
            break
        if status.get('status') == 'error':
            item.update(status='error', query_error=redact(str(status.get('error', 'CAS error')), client))
            save(checkpoint, state)
            raise ValueError('CAS query failed; redacted diagnostic saved')
        print(f"{name}: {status.get('status', 'unknown')}", flush=True)
        sleep(POLL_SECONDS)
        waited += POLL_SECONDS
    partial = results/(name+'.csv.part')
    if partial.exists():
        partial.rename(results/(name+f'.csv.incomplete-{time.time_ns()}'))
    with client.post('download', id=item['job_id']) as response:
        item['transport'] = transport(response); save(checkpoint, state)
        budgeted_download(response, partial, state, checkpoint)
    try:
        require_line_end(partial)
        metadata = check(partial, batch)
    except (ValueError, TypeError, OverflowError, csv.Error):
        item.update(status='downloaded_invalid', partial_file=partial.name, sha256=sha(partial))
        save(checkpoint, state)
        raise
    partial.replace(target)
    item.update(status='complete', bytes=target.stat().st_size, sha256=sha(target), **metadata)
    save(checkpoint, state)
    print(f"part {state['part']} {name}: verified {metadata['rows']} rows, {item['bytes']} bytes", flush=True)


def verify_complete(plan, directory, state):
    payload_bytes, files, per_batch = 0, {}, {}
    results = directory/'results'
    require(set(state['jobs']) == {Path(b['name']).stem for b in plan['batches']}, 'job set differs')
    for batch in plan['batches']:
        name = Path(batch['name']).stem; item = state['jobs'][name]
        file = results/(name+'.csv')
        require(item['status'] == 'complete' and sha(file) == item['sha256'] and file.stat().st_size == item['bytes'], 'completed CSV differs')
        valid = VALIDATORS[batch['kind']](file, ids_of(directory, batch))
        require(valid['rows'] == item['rows'] == batch['n_ids'], 'batch count differs')
        payload_bytes += item['bytes']; files[file.name] = item['sha256']
        per_batch[name] = dict(kind=batch['kind'], rows=valid['rows'])
    for kind in VALIDATORS:
        require(sum(v['rows'] for v in per_batch.values() if v['kind'] == kind) == plan['n_ids'], 'ID coverage differs: '+kind)
    require(payload_bytes <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP, 'byte budget differs')
    return dict(schema=SCHEMA, status='complete_random_id_supplement_download', part=plan['part'],
                approved_scope_plan_sha256=APPROVED_SHA[plan['part']], n_ids=plan['n_ids'],
                n_csv_files=len(files), per_batch=per_batch, verified_csv_bytes=payload_bytes,
                confirmed_bytes_received=state['confirmed_bytes_received'], budget_bytes_charged=state['budget_bytes_charged'],
                maximum_total_received_bytes=CAP, output_sha256=files, code_sha256=state['code_sha256'], credentials_saved=False,
                source_objects_returned=True, science_measurement=False, calibration_applied=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true', help='Verify the local plan only; no network')
    mode.add_argument('--submit', action='store_true')
    parser.add_argument('--part', type=int, required=True, choices=range(1, PARTS+1))
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--user', default='kjhan0606')
    parser.add_argument('--use-env-auth', action='store_true', help='Use HSC_USER/HSC_PASSWORD already set in this shell')
    args = parser.parse_args()
    part = args.part; directory = output(part).resolve(); plan = load_plan(directory, part)
    print(json.dumps(dict(part=part, plan_sha256=APPROVED_SHA[part], ids=plan['n_ids'], batches=len(plan['batches']), maximum_csv_bytes=CAP,
                          mode='local_plan_check' if args.plan else 'approved_transfer')), flush=True)
    if args.plan:
        return 0
    require(APPROVED_SHA[part] is not None and args.approved_plan_sha256 == APPROVED_SHA[part],
            'submission requires explicit user approval of the exact plan hash')
    require(platform.node().lower().startswith('lageunha'), 'run the established CPU workflow on LagEunha')
    results = directory/'results'; checkpoint = results/'status.json'
    if results.exists():
        require(checkpoint.is_file(), 'existing results have no checkpoint')
    else:
        results.mkdir(exist_ok=False)
    with (results/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        code = code_hashes()
        if checkpoint.exists():
            state = json.loads(checkpoint.read_text())
            require(state['schema'] == SCHEMA and state['plan_sha256'] == APPROVED_SHA[part] and state['part'] == part, 'checkpoint scope differs')
            require(0 <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP, 'checkpoint budget differs')
            require(state['code_sha256'] == code, 'receiver dependencies changed; review before resume')
            require(set(state['jobs']) <= {Path(b['name']).stem for b in plan['batches']}, 'checkpoint has unplanned jobs')
        else:
            state = dict(schema=SCHEMA, status='pending', part=part, plan_sha256=APPROVED_SHA[part], code_sha256=code,
                         jobs={}, query_submitted=False, confirmed_bytes_received=0, budget_bytes_charged=0,
                         credentials_saved=False, source_rows_downloaded=False, science_measurement=False)
            save(checkpoint, state)
        if state['status'] == 'complete':
            verify_complete(plan, directory, state)
            require(sha(results/'download_summary.json') == state['summary_sha256'], 'receipt summary differs')
            print('Already complete; verified locally without credentials or network.')
            return 0
        user, password = credentials(args)
        client = Cas(user, password); password = None
        state.update(status='running', last_error_type=None); save(checkpoint, state)
        try:
            for batch in plan['batches']:
                state['current_query'] = Path(batch['name']).stem; save(checkpoint, state)
                fetch_batch(client, batch, batch_text(directory, batch), ids_of(directory, batch), results, state, checkpoint)
            summary = verify_complete(plan, directory, state)
            file = results/'download_summary.json'
            if file.exists():
                require(json.loads(file.read_text()) == summary, 'existing receipt summary differs')
            else:
                save(file, summary)
            state.update(status='complete', current_query=None, summary_sha256=sha(file))
            save(checkpoint, state)
        except BaseException as error:
            state.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'stopped_on_error',
                         last_error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
    print(f'Complete: {len(plan["batches"])} CSVs; science_measurement=false.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f'HSC CAS HTTP {error.code}; no credentials or response body printed.', file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print('Stopped; checkpoint and partial files preserved.', file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        print('S16A random ID supplement: '+(str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
