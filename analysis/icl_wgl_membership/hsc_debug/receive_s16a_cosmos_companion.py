#!/usr/bin/env python3
"""Receive the two pinned S17A Wide-COSMOS companion photometry tables after approval."""
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
from urllib.request import build_opener

from fetch_ancillary_catalogues import NoRedirect, save
from receive_s16a_ephor_ab_cosmos import CHECKPOINT, receive_file, verify_state

DIRECTORY = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_cosmos_companion_plan_2026-11-24_01')
APPROVAL_HASH = '29e3cb624a420b44244ce793674b1561fe1c8041b2d3cce01ed64f5fb4fa8788'


def load_plan(path, approval=APPROVAL_HASH):
    plan = json.loads(path.read_text())
    payload = {key: value for key, value in plan.items() if key != 'plan_sha256'}
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':')).encode()).hexdigest()
    if digest != approval or plan.get('plan_sha256') != approval:
        raise ValueError('plan differs from the pinned two-file companion scope')
    if sum(entry['content_length'] for entry in plan['files']) != plan['maximum_total_payload_bytes']:
        raise ValueError('payload cap differs from the planned file sizes')
    return plan


def credential(args):
    if args.use_env_auth:
        user, password = os.environ.get('HSC_USER') or args.user, os.environ.get('HSC_PASSWORD')
        if not password:
            raise ValueError('HSC_PASSWORD is not set in this shell; no prompt or file fallback')
    else:
        if not sys.stdin.isatty():
            raise ValueError('hidden credential entry requires an interactive terminal')
        user, password = args.user, getpass.getpass('HSC archive credential (hidden): ')
        if not password:
            raise ValueError('empty credential')
    authorization = 'Basic ' + base64.b64encode((user+':'+password).encode()).decode()
    del password
    return authorization


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--plan', action='store_true')
    mode.add_argument('--receive', action='store_true')
    parser.add_argument('--approved-plan-sha256')
    parser.add_argument('--user', default='kjhan0606')
    parser.add_argument('--use-env-auth', action='store_true', help='read HSC_PASSWORD (and HSC_USER) already set in this shell')
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
        if pending:
            authorization = credential(args)
            opener = build_opener(NoRedirect())
            for entry in pending:
                print(f"Receiving {entry['filename']} ({entry['content_length']} bytes)", flush=True)
                receive_file(entry, plan, destination, state, opener.open, authorization)
            del authorization
        verify_state(plan, destination, state)
        print(f'Verified {len(state["files"])} FITS tables. Receipt: {receipt}')


if __name__ == '__main__':
    try:
        main()
    except HTTPError as exc:
        sys.exit(f'HSC archive HTTP {exc.code}. Transfer stopped.')
    except URLError:
        sys.exit('HSC archive connection failed. Transfer stopped.')
