#!/usr/bin/env python3
"""Receive the three pinned Ephor_AB/COSMOS tables after explicit approval."""
from __future__ import annotations

import argparse
import base64
import fcntl
import getpass
import hashlib
import json
import os
from pathlib import Path
import sys
from urllib.error import HTTPError, URLError
from urllib.request import Request, build_opener

from astropy.io import fits

from fetch_ancillary_catalogues import NoRedirect, save, sha

DIRECTORY = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01')
APPROVAL_HASH = '73a28a83e9a2579b290ab68347702e0f2d28adfa9de162bf37e97df360727975'
CHECKPOINT = 'receipt.json'


def load_plan(path):
    plan = json.loads(path.read_text())
    payload = {key: value for key, value in plan.items() if key != 'plan_sha256'}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if digest != APPROVAL_HASH or plan['plan_sha256'] != APPROVAL_HASH:
        raise ValueError('plan differs from the pinned three-file scope')
    return plan


def inspect_fits(path):
    tables = []
    with fits.open(path, mode='readonly', memmap=True) as hdus:
        hdus.verify('exception')
        for hdu in hdus:
            if 'CHECKSUM' in hdu.header and hdu.verify_checksum() != 1:
                raise ValueError('FITS checksum mismatch')
            if 'DATASUM' in hdu.header and hdu.verify_datasum() != 1:
                raise ValueError('FITS data checksum mismatch')
            if isinstance(hdu, fits.BinTableHDU) and hdu.data is not None and len(hdu.data):
                tables.append(dict(hdu=hdu.name, rows=len(hdu.data), columns=list(hdu.columns.names)))
    if not tables:
        raise ValueError('expected a nonempty binary FITS table')
    return tables


def verify_state(plan, destination, state):
    if state['plan_sha256'] != plan['plan_sha256']:
        raise ValueError('receipt belongs to a different plan')
    sizes = {entry['filename']: entry['content_length'] for entry in plan['files']}
    if (not set(state['files']) <= sizes.keys()
            or state['reserved_bytes'] != sum(sizes[name] for name in state['files'])
            or state['reserved_bytes'] > plan['maximum_total_payload_bytes']):
        raise ValueError('receipt exceeds or differs from the byte plan')
    for name, item in state['files'].items():
        path = destination/name
        if item['status'] != 'complete':
            raise ValueError(f'{name}: previous transfer did not complete. Review the preserved partial file and byte reservation before proceeding.')
        if not path.is_file() or path.stat().st_size != sizes[name] or sha(path) != item['sha256']:
            raise ValueError(f'completed FITS changed: {name}')


def receive_file(entry, plan, destination, state, open_response, authorization):
    name, size = entry['filename'], entry['content_length']
    target = destination/name
    partial = target.with_suffix('.fits.part')
    if name in state['files'] or target.exists() or partial.exists():
        raise ValueError(f'file already attempted or present: {name}')
    if state['reserved_bytes'] + size > plan['maximum_total_payload_bytes']:
        raise ValueError('payload budget exhausted')

    # Reserve before GET. A failed transfer does not automatically release the
    # reservation, because transport buffering can precede the first saved byte.
    state['reserved_bytes'] += size
    state['files'][name] = {'status': 'in_progress', 'expected_bytes': size}
    save(destination/CHECKPOINT, state)
    request = Request(entry['url'], method='GET', headers={'Authorization': authorization,
        'If-Match': entry['etag'], 'Accept-Encoding': 'identity'})
    try:
        with open_response(request, timeout=60) as response:
            headers = response.headers
            if (response.status != 200 or headers.get('ETag') != entry['etag']
                    or headers.get('Content-Length') != str(size)
                    or headers.get('Content-Encoding', 'identity') != 'identity'):
                raise ValueError('server response differs from planned ETag or length')
            remaining = size
            with partial.open('xb') as stream:
                while remaining:
                    block = response.read(min(1024**2, remaining))
                    if not block or len(block) > remaining:
                        raise ValueError('truncated or invalid response')
                    stream.write(block)
                    remaining -= len(block)
                stream.flush()
                os.fsync(stream.fileno())
        tables = inspect_fits(partial)
        digest = sha(partial)
        partial.rename(target)
        state['files'][name].update(status='complete', bytes=size, sha256=digest, tables=tables)
        save(destination/CHECKPOINT, state)
    except BaseException:
        state['files'][name]['status'] = 'stopped'
        state['files'][name]['saved_partial_bytes'] = partial.stat().st_size if partial.exists() else 0
        save(destination/CHECKPOINT, state)
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--receive', action='store_true')
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--user', default='kjhan0606')
    args = parser.parse_args()
    plan = load_plan(DIRECTORY/'plan.json')
    if args.plan:
        print(json.dumps(plan, indent=2))
        return
    if args.approved_plan_sha256 != APPROVAL_HASH:
        parser.error('receipt requires explicit user approval of the pinned plan')

    destination = Path(plan['destination'])
    try:
        destination.mkdir(exist_ok=False)
    except FileExistsError:
        if not (destination/CHECKPOINT).is_file():
            raise ValueError('existing results directory has no receipt')
    with (destination/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        receipt = destination/CHECKPOINT
        if receipt.exists():
            state = json.loads(receipt.read_text())
        else:
            state = dict(plan_sha256=APPROVAL_HASH, reserved_bytes=0, files={},
                         science_measurement=False, credentials_saved=False)
            save(receipt, state)
        verify_state(plan, destination, state)
        pending = [entry for entry in plan['files'] if entry['filename'] not in state['files']]
        if pending and not sys.stdin.isatty():
            raise ValueError('hidden credential entry requires an interactive terminal')
        if pending:
            password = getpass.getpass('HSC archive credential (hidden): ')
            if not password:
                raise ValueError('empty credential')
            authorization = 'Basic ' + base64.b64encode((args.user+':'+password).encode()).decode()
            del password
            opener = build_opener(NoRedirect())
            for entry in pending:
                receive_file(entry, plan, destination, state, opener.open, authorization)
        verify_state(plan, destination, state)
        print(f'Verified {len(state["files"])} FITS tables. Receipt: {receipt}')


if __name__ == '__main__':
    try:
        main()
    except HTTPError as exc:
        sys.exit(f'HSC archive HTTP {exc.code}. Transfer stopped.')
    except URLError:
        sys.exit('HSC archive connection failed. Transfer stopped.')
