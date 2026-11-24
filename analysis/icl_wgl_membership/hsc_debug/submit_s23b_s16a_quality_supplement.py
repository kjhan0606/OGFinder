#!/usr/bin/env python3
"""Receive only the approved 31-column S16A pilot-ID supplement, serially.

Uses the established CAS client and conservative download-budget implementation.
No credentials are stored, no SQL is regenerated, no sample is expanded.
"""
from __future__ import annotations

import argparse
from collections import Counter
import fcntl
import getpass
import json
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from submit_s23b_s16a_shape_qa_pilot import budgeted_download, MAX_BYTES
from prepare_s23b_s16a_quality_supplement import OUTPUT, CAP, validate
from s16a_y1_source_contract import EXPRESSIONS, BAD_FLAGS, finite, boolean, object_id

APPROVED_SHA = '44b9d4df4904dac295f08a426abb20aa536ae85fc2c91239ac599f39deb412a1'
SCHEMA = 'astrafek-s16a-quality-supplement-receipt-v1'
CODE_FILES = ('submit_s23b_s16a_quality_supplement.py', 'submit_s23b_s16a_shape_qa_pilot.py',
              'prepare_s23b_s16a_quality_supplement.py', 's16a_y1_source_contract.py',
              'fetch_ancillary_catalogues.py')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def load_plan(directory):
    checked = validate(directory)
    require(checked['plan_sha256'] == APPROVED_SHA, 'not the exact approved supplement plan')
    require(MAX_BYTES == CAP == 134217728, 'shared cumulative download cap differs')
    return json.loads((directory/'quality_supplement_plan.json').read_text())


def batch_ids(directory, batch):
    # SQL was already checked against the original pilot by load_plan().
    import re
    return re.findall(r'\(([1-9][0-9]*)::bigint\)', (directory/batch['name']).read_text())


def validate_csv(file, expected_ids):
    expected = set(expected_ids)
    require(len(expected) == len(expected_ids) > 0, 'invalid expected object IDs')
    seen, missing, bad_booleans = set(), Counter(), Counter()
    flags = ('is_primary', *BAD_FLAGS)
    with Path(file).open(newline='', encoding='utf-8-sig') as stream:
        reader = catalogue_reader(stream)
        require(reader.fieldnames == list(EXPRESSIONS), 'CAS columns differ from approved 31 fields')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), 'malformed CAS CSV row')
            ident = object_id(row['object_id'])
            require(ident in expected and ident not in seen, 'unexpected or duplicate object_id')
            seen.add(ident)
            for key, value in row.items():
                if key == 'object_id':
                    continue
                if key in flags:
                    if boolean(value) is None:
                        missing[key] += 1
                        if value.strip().lower() not in ('', 'null', 'none', 'nan'):
                            bad_booleans[key] += 1
                else:
                    # NULL/nonfinite cells remain visible; never filter or repair.
                    if finite(value) is None:
                        missing[key] += 1
    require(seen == expected, 'missing requested object_id')
    require(not bad_booleans, 'unrecognized non-null CAS boolean representation')
    return dict(rows=len(seen), missing_or_nonfinite_by_column=dict(missing), exact_ids_verified=True)


def redact(message, client):
    for secret in client.credential.values():
        if secret:
            message = message.replace(secret, '[redacted]').replace(secret.replace('/', '\\/'), '[redacted]')
    return message[:2000]


def fetch_batch(client, batch, sql, ids, results, state, checkpoint, sleep=time.sleep):
    name = Path(batch['name']).stem
    item = state['jobs'].setdefault(name, dict(status='pending'))
    target = results/(name+'.csv')
    if item['status'] == 'complete':
        require(target.is_file() and sha(target) == item['sha256'] and target.stat().st_size == item['bytes'], 'completed CSV missing or changed')
        validate_csv(target, ids)
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
    while True:
        status = client.json('status', id=item['job_id'])
        if status.get('status') == 'done':
            break
        if status.get('status') == 'error':
            item.update(status='error', query_error=redact(str(status.get('error', 'CAS error')), client))
            save(checkpoint, state)
            raise ValueError('CAS query failed; redacted diagnostic saved')
        print(f"{name}: {status.get('status', 'unknown')}", flush=True)
        sleep(45)
    partial = results/(name+'.csv.part')
    if partial.exists():
        partial.rename(results/(name+f'.csv.incomplete-{time.time_ns()}'))
    with client.post('download', id=item['job_id']) as response:
        budgeted_download(response, partial, state, checkpoint)
    try:
        metadata = validate_csv(partial, ids)
    except (ValueError, TypeError, OverflowError):
        item.update(status='downloaded_invalid', partial_file=partial.name, sha256=sha(partial))
        save(checkpoint, state)
        raise
    partial.replace(target)
    item.update(status='complete', bytes=target.stat().st_size, sha256=sha(target), **metadata)
    save(checkpoint, state)
    print(f"{name}: verified {metadata['rows']} objects, {item['bytes']} bytes", flush=True)


def verify_complete(plan, directory, state):
    rows, payload_bytes, files, missing = 0, 0, {}, Counter()
    results = directory/'results'
    require(set(state['jobs']) == {Path(b['name']).stem for b in plan['batches']}, 'job set differs')
    for batch in plan['batches']:
        name = Path(batch['name']).stem; item = state['jobs'][name]
        file = results/(name+'.csv')
        require(item['status'] == 'complete' and sha(file) == item['sha256'] and file.stat().st_size == item['bytes'], 'completed CSV differs')
        valid = validate_csv(file, batch_ids(directory, batch))
        require(valid['rows'] == batch['expected_rows'] == item['rows'], 'batch count differs')
        rows += valid['rows']; payload_bytes += item['bytes']; files[file.name] = sha(file)
        missing.update(valid['missing_or_nonfinite_by_column'])
    require(rows == 115993 and payload_bytes <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP, 'row/byte budget differs')
    return dict(schema=SCHEMA, status='complete_quality_supplement_download',
        approved_scope_plan_sha256=APPROVED_SHA, n_csv_files=24, n_unique_objects=rows,
        output_columns=list(EXPRESSIONS), verified_csv_bytes=payload_bytes,
        confirmed_bytes_received=state['confirmed_bytes_received'], budget_bytes_charged=state['budget_bytes_charged'],
        maximum_total_received_bytes=CAP, missing_or_nonfinite_by_column=dict(missing), output_sha256=files,
        code_sha256=state['code_sha256'], credentials_saved=False,
        photoz_point_estimates_returned=True, full_photoz_pdfs_returned=False,
        images_returned=False, psfs_returned=False, new_objects_returned=False,
        science_measurement=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--submit', action='store_true')
    parser.add_argument('--plan-dir', type=Path, default=OUTPUT)
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--user', default='kjhan0606')
    args = parser.parse_args()
    directory = args.plan_dir.resolve(); plan = load_plan(directory)
    print(json.dumps(dict(plan_sha256=APPROVED_SHA, objects=115993, columns=31, batches=24,
                         maximum_csv_bytes=CAP, mode='local_plan_check' if args.plan else 'approved_transfer')), flush=True)
    if args.plan:
        return 0
    require(args.approved_plan_sha256 == APPROVED_SHA, 'submission requires the exact approved plan hash')
    require(platform.node().lower().startswith('lageunha'), 'run the established CPU workflow on LagEunha')
    results = directory/'results'; checkpoint = results/'quality_supplement_status.json'
    if results.exists():
        require(checkpoint.is_file(), 'existing results have no checkpoint')
    else:
        results.mkdir(exist_ok=False)
    with (results/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        code = {name: sha(Path(__file__).with_name(name)) for name in CODE_FILES}
        if checkpoint.exists():
            state = json.loads(checkpoint.read_text())
            require(state['schema'] == SCHEMA and state['plan_sha256'] == APPROVED_SHA, 'checkpoint scope differs')
            require(0 <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP, 'checkpoint budget differs')
            require(state['code_sha256'] == code, 'receiver dependencies changed; review before resume')
            require(set(state['jobs']) <= {Path(b['name']).stem for b in plan['batches']}, 'checkpoint has unplanned jobs')
        else:
            state = dict(schema=SCHEMA, status='pending', plan_sha256=APPROVED_SHA, code_sha256=code,
                         jobs={}, query_submitted=False, confirmed_bytes_received=0, budget_bytes_charged=0,
                         credentials_saved=False, source_rows_downloaded=False, science_measurement=False)
            save(checkpoint, state)
        if state['status'] == 'complete':
            verify_complete(plan, directory, state)
            require(sha(results/'download_summary.json') == state['summary_sha256'], 'receipt summary differs')
            print('Already complete; verified locally without credentials or network.')
            return 0
        require(sys.stdin.isatty(), 'interactive terminal required for hidden credential input')
        password = getpass.getpass('HSC password (hidden): ')
        require(bool(password), 'empty credential')
        client = Cas(args.user, password)
        state.update(status='running', last_error_type=None); save(checkpoint, state)
        try:
            for batch in plan['batches']:
                state['current_query'] = Path(batch['name']).stem; save(checkpoint, state)
                fetch_batch(client, batch, (directory/batch['name']).read_text(), batch_ids(directory, batch), results, state, checkpoint)
            summary = verify_complete(plan, directory, state)
            file = results/'download_summary.json'
            if file.exists():
                require(json.loads(file.read_text()) == summary, 'existing receipt summary differs')
            else:
                save(file, summary)
            state.update(status='complete', current_query=None, summary_sha256=sha(file), total_rows=115993)
            save(checkpoint, state)
        except BaseException as error:
            state.update(status='interrupted' if isinstance(error, KeyboardInterrupt) else 'stopped_on_error',
                         last_error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear(); password = None
    print('Complete: 24 CSVs, 115993 objects; science_measurement=false.')
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
        print('S16A supplement: '+(str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
