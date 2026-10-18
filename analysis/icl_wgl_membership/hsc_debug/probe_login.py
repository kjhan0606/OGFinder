#!/usr/bin/env python3
"""Prompt in the terminal, then run the existing one-cutout header probe.

No credentials are placed in arguments, files, logs, or the parent shell.
The password is read literally, including slash and shell metacharacters.
"""
import getpass
import argparse
import importlib.util
import io
import json
import os
import re
from pathlib import Path
import sys
import urllib.parse
import urllib.request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--global-sky', action='store_true',
                        help='Read the authenticated manual, then probe one coadd/bg with all HDU headers')
    parser.add_argument('--record-name', default='cutout_global',
                        help='New log basename for a global probe; previous records are never replaced')
    args = parser.parse_args()
    if not re.fullmatch(r'[A-Za-z0-9_-]+', args.record_name):
        parser.error('record-name must be a simple log basename')
    if not sys.stdin.isatty():
        print('Run this command in an interactive SSH terminal.', file=sys.stderr)
        return 1
    user = input('HSC account name: ').strip()
    password = getpass.getpass('HSC password (hidden): ')
    if not user or not password:
        print('Both account name and password are required.', file=sys.stderr)
        return 1
    os.environ['HSC_USER'] = user
    os.environ['HSC_PASSWORD'] = password
    path = Path(__file__).with_name('07_fetch_cutouts.py')
    spec = importlib.util.spec_from_file_location('hsc_header_probe', path)
    probe = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(probe)
    try:
        if args.global_sky:
            return global_probe(probe, user, password, args.record_name)
        return probe.main() or 0
    finally:
        os.environ.pop('HSC_USER', None)
        os.environ.pop('HSC_PASSWORD', None)


def global_probe(probe, user, password, record_name='cutout_global'):
    """Keep the manual and headers, never a bulk download or a science claim."""
    from astropy.io import fits
    logs = Path(probe.ROOT) / 'logs'
    logs.mkdir(parents=True, exist_ok=True)
    for suffix in ('_manual.html', '_headers.txt', '_status.json'):
        if (logs / (record_name + suffix)).exists():
            raise ValueError('Previous record exists; use a new --record-name')
    manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
    manager.add_password(None, 'https://hsc-release.mtk.nao.ac.jp/', user, password)
    opener = urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(manager))
    manual_url = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/manual.html'
    with opener.open(manual_url, timeout=60) as response:
        manual = response.read(2000001)
    if len(manual) > 2000000 or b'coadd/bg' not in manual:
        raise ValueError('Manual does not document coadd/bg; no new cutout was requested')
    # Exclusive files preserve every previous authentication/header probe.
    with (logs / (record_name + '_manual.html')).open('xb') as handle:
        handle.write(manual)
    ident, ra, dec = probe.first_selected()
    query = urllib.parse.urlencode(dict(ra=ra, dec=dec, sw='20arcsec', sh='20arcsec',
                                        type='coadd/bg', image='on', mask='on', variance='on',
                                        filter='HSC-I', rerun='pdr3_wide'))
    code, content_type, body, limit = probe.probe(probe.BASE + query, user, password)
    status = {'probe_id': ident, 'requested_product': 'coadd/bg', 'http_status': int(code),
              'credential_env': True, 'content_type': content_type,
              'manual_source': manual_url, 'manual_file': 'logs/' + record_name + '_manual.html',
              'bulk_download': False, 'frames_tsv_written': False, 'science_measurement': False}
    if code == 200 and body.startswith(b'SIMPLE') and len(body) <= limit:
        with fits.open(io.BytesIO(body)) as hdul:
            with (logs / (record_name + '_headers.txt')).open('x', encoding='utf-8') as handle:
                for i, hdu in enumerate(hdul):
                    handle.write('HDU %d\n' % i)
                    handle.write(''.join(card.image.rstrip() + '\n' for card in hdu.header.cards))
                    handle.write('\n')
        status.update(status='all_headers_recorded', header_file='logs/' + record_name + '_headers.txt',
                      product_verified=False, note='Source mapping and inherited primary calibration still need review')
    else:
        status.update(status='refused', product_verified=False)
    with (logs / (record_name + '_status.json')).open('x', encoding='utf-8') as handle:
        json.dump(status, handle, indent=1)
        handle.write('\n')
    print('global probe %s http %s' % (status['status'], code))
    return 0 if code == 200 and status['status'] == 'all_headers_recorded' else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (EOFError, KeyboardInterrupt):
        print('\nProbe cancelled.', file=sys.stderr)
        sys.exit(1)
    except (ValueError, OSError) as exc:
        print('Probe stopped: %s' % exc, file=sys.stderr)
        sys.exit(1)
