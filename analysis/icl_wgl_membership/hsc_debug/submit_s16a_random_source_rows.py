#!/usr/bin/env python3
"""Receive only the approved S16A source-pair and inner random-point rows for the 76 random centres, serially.

Same safeguards as the other S16A receivers: exact plan hash, pinned plan directory, cumulative byte budget,
SQL hash re-check at use, no stored credentials (--use-env-auth or one hidden prompt), transport framing and
final line break, server window totals before LIMIT (per centre and per batch), and a centre at its limit is
refused. Every centre must return rows. Pair sources must lie in their 3-Mpc cone (degree-valued ira/idec,
0.01 arcsec as in the pilot receiver); inner points in their 0.4 Mpc + 30 arcsec cone.
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
from submit_s23b_s16a_shape_qa_pilot import FLOAT_COLUMNS, MAX_BYTES, budgeted_download, numeric
from submit_s16a_afterburner_photoz_rows import batch_text, code_hashes, credentials
from audit_s23b_s16a_wl_random_support import angular_separation_arcmin
from cas_transfer_checks import redact, require_line_end, require_poll_budget, transport
from prepare_s16a_random_candidate_rows import unit_vector
import prepare_s16a_random_inner_rows as pi
from prepare_s16a_random_source_rows import INNER_COLUMNS, N_CENTRES, OUTPUT, PAIR_COLUMNS, validate
from prepare_s23b_s16a_quality_supplement import CAP, require
from s16a_y1_source_contract import boolean

APPROVED_SHA = 'c1dfe5a4860ceb830281e1444b90b438a99e93702726490518fa32a90073b211'  # user approved 2026-11-26
SCHEMA = 'astrafek-s16a-random-source-rows-receipt-v1'
POLL_SECONDS = 20
CONE_TOLERANCE = 1e-12
PAIR_EXCESS_ARCSEC = 0.01


def load_plan(directory):
    checked = validate(directory)
    require(APPROVED_SHA is not None and checked['plan_sha256'] == APPROVED_SHA, 'not the exact approved random source-row plan')
    require(MAX_BYTES == CAP == 134217728, 'shared cumulative download cap differs')
    return json.loads((directory/'plan.json').read_text())


def counts(row, matched, totals):
    require(all(re.fullmatch(r'[0-9]+', row[k].strip()) for k in ('n_matched', 'n_batch')), 'invalid match counts')
    matched[row['center_id']].add(int(row['n_matched'])); totals.add(int(row['n_batch']))


def finish(batch, seen, matched, totals):
    require(all(0 < len(v) < batch['row_limit_per_centre'] for v in seen.values()), 'empty centre or per-centre row limit reached')
    require(all(matched[c] == {len(v)} for c, v in seen.items()), 'rows of a centre differ from its server match count')
    require(totals == {sum(len(v) for v in seen.values())}, 'rows differ from the batch match count (a centre may be missing)')


def validate_pairs(file, batch):
    centres = {c['center_id']: c for c in batch['centres']}
    seen, matched, totals = {c: set() for c in centres}, {c: set() for c in centres}, set()
    nonfinite, excess = Counter(), -math.inf
    with Path(file).open(newline='', encoding='utf-8-sig') as stream:
        reader = catalogue_reader(stream)
        require(reader.fieldnames == PAIR_COLUMNS, 'CAS columns differ from the approved pair fields')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), 'malformed CAS CSV row')
            cid = row['center_id']
            require(cid in centres, 'row belongs to another centre')
            counts(row, matched, totals)
            require(re.fullmatch(r'[1-9][0-9]*', row['object_id'].strip()) is not None, 'invalid source object_id')
            ident = int(row['object_id'])
            require(ident not in seen[cid], 'duplicate (centre, object) pair')
            seen[cid].add(ident)
            values = {k: numeric(row[k]) for k in FLOAT_COLUMNS}
            ra, dec = values['source_ra'], values['source_dec']
            require(math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90, 'invalid raw sky coordinate')
            c = centres[cid]
            e = angular_separation_arcmin(c['ra_deg'], c['dec_deg'], ra, dec)*60-c['radius_arcsec']
            require(e <= PAIR_EXCESS_ARCSEC, 'source outside the degree-valued query cone; no conversion attempted')
            excess = max(excess, e)
            nonfinite.update(k for k, v in values.items() if not math.isfinite(v))
    finish(batch, seen, matched, totals)
    return dict(rows=sum(len(v) for v in seen.values()), rows_by_centre={k: len(v) for k, v in seen.items()},
                nonfinite_by_column=dict(nonfinite), maximum_cone_excess_arcsec=excess)


def validate_inner(file, batch):
    cones = {c['center_id']: (unit_vector(float(c['ra_deg']), float(c['dec_deg'])), math.cos(math.radians(c['query_radius_arcsec']/3600)))
             for c in batch['centres']}
    seen, matched, totals = {c: set() for c in cones}, {c: set() for c in cones}, set()
    null, nonfinite, true = Counter(), Counter(), Counter()
    with Path(file).open(newline='', encoding='utf-8-sig') as stream:
        reader = catalogue_reader(stream)
        require(reader.fieldnames == INNER_COLUMNS, 'CAS columns differ from the approved inner random fields')
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()), 'malformed CAS CSV row')
            cid = row['center_id']
            require(cid in cones, 'row belongs to another centre')
            counts(row, matched, totals)
            require(re.fullmatch(r'[1-9][0-9]*', row['object_id'].strip()) is not None, 'invalid random object_id')
            ident = int(row['object_id'])
            require(ident not in seen[cid], 'duplicate random object_id within a centre')
            seen[cid].add(ident)
            ra, dec = float(row['ra']), float(row['dec'])
            require(math.isfinite(ra) and math.isfinite(dec), 'non-finite random position')
            centre, cos_r = cones[cid]
            require(sum(a*b for a, b in zip(unit_vector(ra, dec), centre)) >= cos_r-CONE_TOLERANCE, 'random point outside the planned cone')
            for key in INNER_COLUMNS[2:-2]:
                value = row[key].strip()
                if value == '':
                    null[key] += 1
                elif key in pi.BOOLEANS:
                    flag = boolean(value)
                    require(flag is not None, 'unrecognized CAS boolean representation')
                    true[key] += flag
                elif key in pi.INTEGERS:
                    require(re.fullmatch(r'-?[0-9]+', value) is not None, 'non-integer value in '+key)
                else:
                    nonfinite[key] += not math.isfinite(float(value))
    finish(batch, seen, matched, totals)
    return dict(rows=sum(len(v) for v in seen.values()), rows_by_centre={k: len(v) for k, v in seen.items()},
                null_by_column=dict(null), nonfinite_by_column={k: v for k, v in nonfinite.items() if v},
                true_by_column={k: v for k, v in true.items() if v})


VALIDATORS = dict(pairs=validate_pairs, inner=validate_inner)


def fetch_batch(client, batch, sql, results, state, checkpoint, sleep=time.sleep):
    name = Path(batch['name']).stem
    check = VALIDATORS[batch['kind']]
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
    print(f"{name}: verified {metadata['rows']} rows, {item['bytes']} bytes", flush=True)


def verify_complete(plan, directory, state):
    payload_bytes, files, per_batch = 0, {}, {}
    results = directory/'results'
    require(set(state['jobs']) == {Path(b['name']).stem for b in plan['batches']}, 'job set differs')
    for batch in plan['batches']:
        name = Path(batch['name']).stem; item = state['jobs'][name]
        file = results/(name+'.csv')
        require(item['status'] == 'complete' and sha(file) == item['sha256'] and file.stat().st_size == item['bytes'], 'completed CSV differs')
        valid = VALIDATORS[batch['kind']](file, batch)
        require(valid['rows'] == item['rows'], 'batch count differs')
        payload_bytes += item['bytes']; files[file.name] = item['sha256']
        per_batch[name] = dict(kind=batch['kind'], lens_id=batch['lens_id'], rows=valid['rows'], rows_by_centre=valid['rows_by_centre'])
    for kind in VALIDATORS:
        require(sum(len(v['rows_by_centre']) for v in per_batch.values() if v['kind'] == kind) == N_CENTRES, 'centre set differs: '+kind)
    require(payload_bytes <= state['confirmed_bytes_received'] <= state['budget_bytes_charged'] <= CAP, 'byte budget differs')
    return dict(schema=SCHEMA, status='complete_random_source_rows_download', approved_scope_plan_sha256=APPROVED_SHA,
                n_csv_files=len(files), per_batch=per_batch, verified_csv_bytes=payload_bytes,
                confirmed_bytes_received=state['confirmed_bytes_received'], budget_bytes_charged=state['budget_bytes_charged'],
                maximum_total_received_bytes=CAP, output_sha256=files, code_sha256=state['code_sha256'], credentials_saved=False,
                source_objects_returned=True, science_measurement=False, calibration_applied=False)


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
    print(json.dumps(dict(plan_sha256=APPROVED_SHA, centres=N_CENTRES, batches=len(plan['batches']), maximum_csv_bytes=CAP,
                          mode='local_plan_check' if args.plan else 'approved_transfer')), flush=True)
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
        print('S16A random source rows: '+(str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        sys.exit(1)
