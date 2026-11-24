#!/usr/bin/env python3
"""HEAD-only inventory of two documented COSMOS companions; never GET FITS."""
import argparse
import base64
from datetime import datetime, timezone
import getpass
import hashlib
import json
import os
from pathlib import Path
import sys
from urllib.error import HTTPError, URLError
from urllib.request import Request, HTTPRedirectHandler, build_opener

DOCUMENT='https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/'
BASE='https://hsc-release.mtk.nao.ac.jp/archive/filetree/cosmos_photoz_catalog_reweighted_to_s16a_shape_catalog/'
FILES=('target_wide_s17a_9812.fits','target_wide_s17a_9813.fits')


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):
        return None


def metadata(name, opener, authorization=None):
    headers={'User-Agent':'Astrafek-COSMOS-metadata/1','Accept-Encoding':'identity'}
    if authorization:
        headers['Authorization']=authorization
    req=Request(BASE+name,method='HEAD',headers=headers)
    try:
        response=opener.open(req,timeout=30)
    except HTTPError as error:
        response=error
    except URLError:
        return dict(filename=name,url=BASE+name,status='connection_error',file_size_bytes=None)
    with response:
        status=response.code
        length=response.headers.get('Content-Length')
        # An error page's length must never become a proposed file size.
        size=int(length) if status==200 and length and length.isdecimal() else None
        return dict(filename=name,url=BASE+name,http_status=status,
                    file_size_bytes=size,etag=response.headers.get('ETag') if status==200 else None,
                    last_modified=response.headers.get('Last-Modified') if status==200 else None,
                    content_type=response.headers.get('Content-Type'),
                    authentication_required=status==401)


def run(output,use_env_auth=False,prompt_user=None):
    if output.exists():
        raise FileExistsError(output)
    authorization=None
    if prompt_user:
        if not sys.stdin.isatty():
            raise ValueError('hidden credential entry requires an interactive terminal')
        password=getpass.getpass('HSC archive credential (hidden; HEAD only): ')
        if not password:
            raise ValueError('empty credential')
        authorization='Basic '+base64.b64encode((prompt_user+':'+password).encode()).decode()
        del password
    elif use_env_auth:
        user,password=os.environ.get('HSC_USER'),os.environ.get('HSC_PASSWORD')
        if not user or not password:
            raise ValueError('HSC_USER and HSC_PASSWORD must already be set; no prompt or file fallback')
        authorization='Basic '+base64.b64encode((user+':'+password).encode()).decode()
        del password
    opener=build_opener(NoRedirect())
    rows=[metadata(name,opener,authorization) for name in FILES]
    report=dict(documentation_url=DOCUMENT,checked_utc=datetime.now(timezone.utc).isoformat(),
        method='HEAD',fits_body_requests=0,credentials_saved=False,
        authenticated=authorization is not None,files=rows,
        science_measurement=False,download_approved=False,
        size_plan_ready=all(r.get('file_size_bytes') is not None for r in rows),
        limitations=['No FITS schema or payload inspected. HEAD sizes do not verify file integrity.',
                    'A subsequent FITS receipt requires a separate explicit scope approval.'],
        runner_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    output.parent.mkdir(parents=True,exist_ok=True)
    with output.open('x') as f:
        json.dump(report,f,indent=2,allow_nan=False);f.write('\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--use-env-auth',action='store_true')
    parser.add_argument('--prompt-auth',action='store_true',help='one hidden credential prompt in an interactive terminal')
    parser.add_argument('--user',default='kjhan0606')
    args=parser.parse_args()
    run(args.output.resolve(),args.use_env_auth,args.user if args.prompt_auth else None)
