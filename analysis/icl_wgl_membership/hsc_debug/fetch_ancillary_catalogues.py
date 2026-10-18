#!/usr/bin/env python3
"""Fetch prepared HSC CAS catalogues with one hidden password prompt.

Uses the documented hscReleaseQuery API protocol. No image download, automatic
scientific acceptance, email, credential persistence, or concurrent CAS jobs.
"""
import argparse
import csv
import getpass
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
import urllib.error
import urllib.request

API = 'https://hsc-release.mtk.nao.ac.jp/datasearch/api/catalog_jobs/'
DEFAULT = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        raise ValueError('CAS redirect refused; credentials were not forwarded')


class Cas:
    def __init__(self, user, password):
        self.credential = dict(account_name=user, password=password)
        self.opener = urllib.request.build_opener(NoRedirect())

    def post(self, action, **fields):
        data = dict(fields, credential=self.credential, clientVersion=20190514.1)
        request = urllib.request.Request(API+action, json.dumps(data).encode(),
                                         {'Content-Type': 'application/json'})
        return self.opener.open(request, timeout=120)

    def json(self, action, **fields):
        with self.post(action, **fields) as response:
            return json.load(response)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def save(path, state):
    temporary = path.with_suffix('.tmp')
    with temporary.open('w') as stream:
        json.dump(state, stream, indent=2, allow_nan=False)
        stream.flush(); os.fsync(stream.fileno())
    temporary.replace(path)


def plan(requests, qa_sql, batch):
    if not requests or len({r['cluster_id'] for r in requests}) != len(requests):
        raise ValueError('Request list must have unique cluster IDs')
    jobs = [('patch_qa', qa_sql)]
    for kind, field in [('extinction', 'extinction_sql'), ('random', 'random_candidates_sql')]:
        for start in range(0, len(requests), batch):
            queries = [r[field].strip().removesuffix(';') for r in requests[start:start+batch]]
            if any(not q.upper().startswith('SELECT ') or ';' in q for q in queries):
                raise ValueError('Only prepared single SELECT statements are accepted')
            jobs.append((f'{kind}_{start//batch:03d}', '\nUNION ALL\n'.join('('+q+')' for q in queries)+';'))
    return jobs


def catalogue_reader(stream):
    """HSC CAS prefixes the first CSV column with '# ', even without metadata."""
    reader = csv.DictReader(stream)
    names = [name.strip() for name in (reader.fieldnames or [])]
    if names:
        names[0] = names[0].removeprefix('#').strip()
    if not names or any(not name for name in names) or len(set(names)) != len(names):
        raise ValueError('Missing or duplicate CSV column names')
    reader.fieldnames = names
    return reader


def validate_csv(path, name):
    with path.open(newline='', encoding='utf-8-sig') as stream:
        reader = catalogue_reader(stream)
        names = reader.fieldnames or []
        required = {'target_cluster_id', 'object_id', 'ra', 'dec'} if name != 'patch_qa' else {'skymap_id', 'ra', 'dec'}
        if name.startswith('extinction_'): required.add('a_i')
        if not names or not required.issubset(names):
            raise ValueError('Unexpected CSV columns; preserve download for schema review')
        count = 0
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError('Malformed CSV record')
            count += 1
        if name == 'patch_qa' and count == 0:
            raise ValueError('Empty patch QA table')
    return count


def fetch_one(client, name, sql, work, state, checkpoint, sleep=time.sleep):
    item = state['jobs'].setdefault(name, {'status': 'pending'})
    target = work/(name+'.csv')
    if item['status'] == 'complete':
        if not target.is_file() or sha(target) != item['sha256']:
            raise ValueError('Saved result missing or changed: '+name)
        print(name+': verified existing result', flush=True)
        return
    if item['status'] == 'submitting' and 'job_id' not in item:
        raise ValueError('Submission outcome unknown for '+name+'; check CAS job list before resubmission')
    if item['status'] == 'error':
        raise ValueError('CAS query failed for '+name+'; review query_error and schema before retry')
    if 'job_id' not in item:
        item['status'] = 'submitting'; save(checkpoint, state)
        try:
            result = client.json('submit', catalog_job=dict(sql=sql, out_format='csv',
                include_metainfo_to_body=False, release_version='pdr3'), nomail=True, skip_syntax_check=False)
        except urllib.error.HTTPError as error:
            if error.code in (400, 401, 403, 406, 429):
                item['status'] = 'pending'; save(checkpoint, state)
            raise
        item.update(job_id=result['id'], status='submitted'); save(checkpoint, state)
    while True:
        result = client.json('status', id=item['job_id'])
        if result['status'] == 'done': break
        if result['status'] == 'error':
            error = str(result.get('error', 'query failed'))
            for secret in client.credential.values():
                if secret: error = error.replace(secret, '[redacted]')
            item.update(status='error', query_error=error[:2000]); save(checkpoint,state)
            raise ValueError(name+': '+error[:2000])
        print(name+': '+str(result['status']), flush=True)
        sleep(60)
    partial = work/(name+'.csv.part')
    # A previous partial download is retained, not mistaken for a complete table.
    if partial.exists(): partial.rename(work/(name+f'.csv.incomplete-{time.time_ns()}'))
    with client.post('download', id=item['job_id']) as response, partial.open('xb') as stream:
        for block in iter(lambda: response.read(1024*1024), b''):
            stream.write(block)
    rows = validate_csv(partial, name)
    partial.replace(target)
    item.update(status='complete', rows=rows, bytes=target.stat().st_size, sha256=sha(target))
    save(checkpoint,state)
    print(f'{name}: saved {rows} rows', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--user', default='kjhan0606')
    p.add_argument('--prepared', type=Path, default=DEFAULT/'ancillary_2026-10-13_01')
    p.add_argument('--work', type=Path, default=DEFAULT/'catalogues_2026-10-13_01')
    p.add_argument('--batch-size', type=int, default=10)
    p.add_argument('--plan', action='store_true', help='show counts; no network or password')
    args = p.parse_args()
    if not 1 <= args.batch_size <= 25: p.error('batch-size must be 1..25')
    req = args.prepared/'hsc_catalogue_requests.json'
    qa = args.prepared/'patch_qa.sql'
    jobs = plan(json.loads(req.read_text()), qa.read_text(), args.batch_size)
    fingerprint = hashlib.sha256(json.dumps(jobs).encode()).hexdigest()
    print(f'{len(jobs)} sequential catalogue queries; existing verified results will be reused.', flush=True)
    if args.plan: return 0
    if not platform.node().lower().startswith('lageunha'): p.error('Run on lageunha CPU')
    if not sys.stdin.isatty(): p.error('Run in an interactive terminal for the one hidden password prompt')
    args.work.mkdir(parents=True, exist_ok=True)
    # Prevent two terminals submitting duplicate work into the same checkpoint.
    import fcntl
    with (args.work/'run.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        checkpoint = args.work/'catalogue_status.json'
        state = json.loads(checkpoint.read_text()) if checkpoint.exists() else dict(
            schema=1, status='pending', plan_sha256=fingerprint, jobs={},
            credentials_saved=False, science_measurement=False)
        if state['plan_sha256'] != fingerprint: raise ValueError('Query plan changed; use a new output directory')
        password = getpass.getpass('HSC password (hidden; one prompt): ')
        if not password: raise ValueError('Empty password')
        client = Cas(args.user, password)
        state['status'] = 'running'; save(checkpoint, state)
        try:
            for name, sql in jobs:
                state['current_query'] = name; save(checkpoint,state)
                fetch_one(client, name, sql, args.work, state, checkpoint)
            state['status'] = 'complete'; save(checkpoint,state)
        except BaseException as error:
            state['status'] = 'interrupted' if isinstance(error, KeyboardInterrupt) else 'stopped_on_error'
            state['error_type'] = type(error).__name__
            save(checkpoint,state)
            raise
        finally:
            client.credential.clear(); password = None
    print('Catalogue downloads complete. Matching, extinction acceptance and control-image downloads remain separate.')
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f'HSC CAS HTTP {error.code}; no credentials or response body printed. Saved jobs/results remain resumable.', file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print('\nStopped; submitted CAS job and completed files retained for resume.', file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        # Application errors are curated; do not serialize request objects or tracebacks.
        message = str(error) if isinstance(error, ValueError) else type(error).__name__
        print('fetch_ancillary_catalogues: '+message, file=sys.stderr)
        sys.exit(1)
