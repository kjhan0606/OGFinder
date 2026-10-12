#!/usr/bin/env python3
"""Probe one documented HSC PDR3 cutout. This script does not guess a URL.

The endpoint and the query keys are the PDR3 cutout client in unagi
(config.Rerun img_url and DEFAULT_CUTOUT_CENTER):

    https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout
    ra dec sw sh type=coadd image mask variance filter=HSC-I rerun=pdr3_wide

Optional basic authentication reads HSC_USER and HSC_PASSWORD from the
environment. The values are not written to the log. A response that is not a
FITS image stops the download. A FITS response is reduced to the header of
the first 2-D HDU and is not copied into cutouts/. The bulk download stays
off until that header is identified as the global-sky coadd. calexp, the
local-sky coadd, is not an ICL input.
"""
import json
import os
import sys
import urllib.error
import urllib.request

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')
BASE = 'https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout?'


def first_selected():
    path = os.path.join(ROOT, 'logs', 'selected_clusters.tsv')
    with open(path, encoding='utf-8') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        line = handle.readline().rstrip('\n').split('\t')
    return line[header.index('id')], line[header.index('ra_deg')], line[header.index('dec_deg')]


def write_status(payload):
    path = os.path.join(ROOT, 'logs', 'cutout_fetch_status.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=1)
        handle.write('\n')
    print('cutout fetch %s http %s' % (payload['status'], payload.get('http_status')))


def opener_for(user, password):
    if user and password:
        manager = urllib.request.HTTPPasswordMgrWithDefaultRealm()
        manager.add_password(None, BASE, user, password)
        return urllib.request.build_opener(urllib.request.HTTPBasicAuthHandler(manager))
    return urllib.request.build_opener()


def probe(url, user, password, limit=40000000):
    request = urllib.request.Request(url)
    opener = opener_for(user, password)
    try:
        with opener.open(request, timeout=90) as response:
            body = response.read(limit + 1)
            return response.status, response.headers.get('Content-Type', ''), body, limit
    except urllib.error.HTTPError as exc:
        return exc.code, exc.headers.get('Content-Type', ''), exc.read(4096), limit


def header_text(body):
    from astropy.io import fits
    import io
    with fits.open(io.BytesIO(body)) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, 'ndim', 0) == 2:
                return ''.join(card.image.rstrip() + '\n' for card in hdu.header.cards)
    raise SystemExit('probe FITS has no 2-D image')


def main():
    cluster_id, ra, dec = first_selected()
    query = '&'.join([
        'ra=%s' % ra,
        'dec=%s' % dec,
        'sw=20arcsec',
        'sh=20arcsec',
        'type=coadd',
        'image=on',
        'mask=off',
        'variance=off',
        'filter=HSC-I',
        'rerun=pdr3_wide',
    ])
    url = BASE + query
    user = os.environ.get('HSC_USER') or ''
    password = os.environ.get('HSC_PASSWORD') or ''
    code, content_type, body, limit = probe(url, user, password)
    fits_ok = body.startswith(b'SIMPLE') and len(body) <= limit
    status = {
        'endpoint': BASE,
        'query_keys': [
            'ra', 'dec', 'sw', 'sh', 'type', 'image', 'mask', 'variance', 'filter', 'rerun',
        ],
        'probe_id': cluster_id,
        'probe_size': '20arcsec',
        'http_status': int(code),
        'content_type': content_type,
        'fits': bool(fits_ok),
        'credential_env': bool(user and password),
        'downloaded': False,
        'bulk_download': False,
        'frames_tsv_written': False,
        'science_measurement': False,
    }
    if code == 401:
        status['status'] = 'unauthorized'
        status['reason'] = 'the cutout service returned 401 Authorization Required'
        write_status(status)
        return 0
    if code != 200 or not fits_ok:
        status['status'] = 'refused'
        status['reason'] = 'the cutout response is not a FITS image'
        write_status(status)
        return 0
    logs = os.path.join(ROOT, 'logs')
    header_path = os.path.join(logs, 'cutout_probe_header.txt')
    with open(header_path, 'w', encoding='utf-8') as handle:
        handle.write(header_text(body))
    status['status'] = 'header_recorded'
    status['reason'] = (
        'the probe header is in logs/cutout_probe_header.txt; '
        'bulk download stays off until that header is the global-sky coadd'
    )
    status['header_file'] = 'logs/cutout_probe_header.txt'
    write_status(status)
    return 0


if __name__ == '__main__':
    sys.exit(main() or 0)
