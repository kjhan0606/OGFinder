#!/usr/bin/env python3
"""Receive only the approved aggregated S16A FDFC grid query, serially.

Same safeguards as the other S16A receivers: exact plan hash, pinned plan directory, cumulative byte
budget, SQL hash re-check at use, no stored credentials (--use-env-auth or one hidden prompt). Each
batch must hold one lens, unique cells inside the planned cone, counts that match the server window
totals, and fewer rows than the planned limit (a batch at the limit is treated as truncated).
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import fcntl
import json
import math
from pathlib import Path
import platform
import re
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from submit_s23b_s16a_shape_qa_pilot import budgeted_download, MAX_BYTES
from submit_s16a_afterburner_photoz_rows import batch_text, code_hashes, credentials
from cas_transfer_checks import redact, require_line_end, require_poll_budget, transport
from prepare_s16a_fdfc_grid import CELL_ARCMIN, OUTPUT, OUTPUT_COLUMNS, validate
from prepare_s16a_random_candidate_rows import N_LENSES, RADIUS_ARCMIN
from prepare_s23b_s16a_quality_supplement import CAP, require

APPROVED_SHA = '019c8228ba2f39264763bfaa9defd2f7e588f02b72d1a93d376d42db41f0189a'  # user approved 2026-11-25
SCHEMA = 'astrafek-s16a-fdfc-grid-receipt-v1'
POLL_SECONDS = 20


def load_plan(directory):
    checked = validate(directory)
    require(APPROVED_SHA is not None and checked['plan_sha256'] == APPROVED_SHA, 'not the exact approved FDFC grid plan')
    require(MAX_BYTES == CAP == 134217728, 'shared cumulative download cap differs')
    return json.loads((directory/'plan.json').read_text())


def validate_csv(file, batch):
    """One lens per batch, unique cells inside the cone's cell range, counts consistent, fewer rows than the limit."""
    lens = batch['lens']
    half = math.ceil(RADIUS_ARCMIN/CELL_ARCMIN)+1
    cells, objects, fdfc, nulls, n_cells, n_objects = set(), 0, 0, 0, set(), set()
    with Path(file).open(newline='', encoding='utf-8-sig') as stream:
        reader = catalogue_reader(stream)
        require(reader.fieldnames == OUTPUT_COLUMNS, f'CAS columns differ from the approved {len(OUTPUT_COLUMNS)} fields')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), 'malformed CAS CSV row')
            require(row['lens_id'] == lens['lens_id'], 'row belongs to another lens')
            values = {}
            for key in OUTPUT_COLUMNS[1:]:
                text = row[key].strip()
                require(re.fullmatch(r'-?[0-9]+', text) is not None, 'non-integer value in '+key)
                values[key] = int(text)
            x, y = values['cell_x'], values['cell_y']
            require(abs(x) <= half and abs(y) <= half, 'cell outside the planned cone')
            require((x, y) not in cells, 'duplicate cell')
            cells.add((x, y))
            require(values['n_primary'] >= 1 and 0 <= values['n_fdfc'] and 0 <= values['n_fdfc_null']
                    and values['n_fdfc']+values['n_fdfc_null'] <= values['n_primary'], 'inconsistent cell counts')
            objects += values['n_primary']; fdfc += values['n_fdfc']; nulls += values['n_fdfc_null']
            n_cells.add(values['n_cells']); n_objects.add(values['n_objects'])
    require(0 < len(cells) < batch['row_limit'], 'empty batch or row limit reached (possible truncation)')
    require(n_cells == {len(cells)} and n_objects == {objects}, 'rows differ from the server totals (truncated or cut response)')
    return dict(rows=len(cells), n_primary=objects, n_fdfc=fdfc, n_fdfc_null=nulls)


def fetch_batch(client, batch, sql, results, state, checkpoint, sleep=time.sleep):
    name = Path(batch['name']).stem
    item = state['jobs'].setdefault(name, dict(status='pending'))
    target = results/(name+'.csv')
    if item['status'] == 'complete':
        require(target.is_file() and sha(target) == item['sha256'] and target.stat().st_size == item['bytes'], 'completed CSV missing or changed')
        validate_csv(target, batch)
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
        metadata = validate_csv(partial, batch)
    except (ValueError, TypeError, OverflowError, csv.Error):
        item.update(status='downloaded_invalid', partial_file=partial.name, sha256=sha(partial))
        save(checkpoint, state)
        raise
    partial.replace(target)
    item.update(status='complete', bytes=target.stat().st_size, sha256=sha(target), **metadata)
    save(checkpoint, state)
    print(f"{name}: verified {metadata['rows']} cells, {item['bytes']} bytes", flush=True)


def verify_complete(plan, directory, state):
    rows, payload_bytes, files = 0, 0, {}
    totals, per_lens = Counter(), {}
    results = directory/'results'
    require(set(state['jobs']) == {Path(b['name']).stem for b in plan['batches']}, 'job set differs')
    for batch in plan['batches']:
        name = Path(batch['name']).stem; item = state['jobs'][name]
        file = results/(name+'.csv')
        require(item['status'] == 'complete' and sha(file) == item['sha256'] and file.stat().st_size == item['bytes'], 'completed CSV differs')
        valid = validate_csv(file, batch)
        require(valid['rows'] == item['rows'], 'batch count differs')
        rows += valid['rows']; payload_bytes += item['bytes']; files[file.name] = item['sha256']; per_lens[batch['lens']['lens_id']] = valid['rows']
        totals.update({k: valid[k] for k in ('n_primary', 'n_fdfc', 'n_fdfc_null')})
    require(len(per_lens) == N_LENSES and payload_bytes <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP,
            'lens/byte budget differs')
    return dict(schema=SCHEMA, status='complete_fdfc_grid_download',
        approved_scope_plan_sha256=APPROVED_SHA, n_csv_files=len(files), n_cells=rows, cells_per_lens=per_lens,
        output_columns=OUTPUT_COLUMNS, verified_csv_bytes=payload_bytes,
        confirmed_bytes_received=state['confirmed_bytes_received'], budget_bytes_charged=state['budget_bytes_charged'],
        maximum_total_received_bytes=CAP, object_totals=dict(totals), output_sha256=files, code_sha256=state['code_sha256'], credentials_saved=False,
        object_rows_returned=False, science_measurement=False, calibration_applied=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true', help='Verify the local plan only; no network')
    mode.add_argument('--submit', action='store_true')
    parser.add_argument('--plan-dir', type=Path, default=OUTPUT)
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--user', default='kjhan0606')
    parser.add_argument('--use-env-auth', action='store_true', help='Use HSC_USER/HSC_PASSWORD already set in this shell')
    args = parser.parse_args()
    directory = args.plan_dir.resolve(); plan = load_plan(directory)
    print(json.dumps(dict(plan_sha256=APPROVED_SHA, lenses=N_LENSES, columns=len(OUTPUT_COLUMNS), batches=len(plan['batches']),
                          maximum_csv_bytes=CAP, mode='local_plan_check' if args.plan else 'approved_transfer')), flush=True)
    if args.plan:
        return 0
    require(APPROVED_SHA is not None and args.approved_plan_sha256 == APPROVED_SHA,
            'submission requires explicit user approval of the exact plan hash')
    require(directory == OUTPUT.resolve(), 'run only from the pinned plan directory; a copy would reset the byte budget')
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
        user, password = credentials(args)
        client = Cas(user, password); password = None
        state.update(status='running', last_error_type=None); save(checkpoint, state)
        try:
            for batch in plan['batches']:
                state['current_query'] = Path(batch['name']).stem; save(checkpoint, state)
                fetch_batch(client, batch, batch_text(directory, batch), results, state, checkpoint)
            summary = verify_complete(plan, directory, state)
            file = results/'download_summary.json'
            if file.exists():
                require(json.loads(file.read_text()) == summary, 'existing receipt summary differs')
            else:
                save(file, summary)
            state.update(status='complete', current_query=None, summary_sha256=sha(file), total_rows=summary['n_cells'])
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
        print('S16A FDFC grid: '+(str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
