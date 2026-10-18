#!/usr/bin/env python3
"""Anonymous Gaia DR3 position crossmatch; evidence only, no PSF approval."""
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import platform
import time
import urllib.request
import urllib.error
import uuid
from astropy.table import Table

ENDPOINT='https://gea.esac.esa.int/tap-server/tap/sync'
QUERY='''SELECT t.candidate_id, g.source_id, g.ra, g.dec, g.ref_epoch,
g.parallax_over_error, g.pmra, g.pmra_error, g.pmdec, g.pmdec_error,
g.ruwe, g.duplicated_source, g.phot_g_mean_mag, g.phot_variable_flag,
g.in_qso_candidates, g.in_galaxy_candidates, g.classprob_dsc_combmod_star,
DISTANCE(POINT('ICRS',t.ra,t.dec),POINT('ICRS',g.ra,g.dec))*3600 AS separation_arcsec
FROM tap_upload.targets AS t JOIN gaiadr3.gaia_source AS g
ON 1=CONTAINS(POINT('ICRS',g.ra,g.dec),CIRCLE('ICRS',t.ra,t.dec,2.0/3600.0))'''


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--requests',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--async-query',action='store_true')
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    req=json.loads(a.requests.read_text())['candidates'];a.work.mkdir(parents=True,exist_ok=False)
    table=Table(dict(candidate_id=list(range(len(req))),ra=[r['ra_deg'] for r in req],dec=[r['dec_deg'] for r in req]))
    table.write(a.work/'targets.vot',format='votable');(a.work/'query.adql').write_text(QUERY)
    (a.work/'candidates.json').write_text(json.dumps(req,indent=2))
    boundary='----GaiaTargets'+uuid.uuid4().hex;parts=[]
    fields=dict(REQUEST='doQuery',LANG='ADQL',FORMAT='csv',MAXREC='20000',UPLOAD='targets,param:table1',QUERY=QUERY)
    if a.async_query:fields['PHASE']='RUN'
    for key,value in fields.items():
        parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="{key}"\r\n\r\n{value}\r\n'.encode())
    parts.append(f'--{boundary}\r\nContent-Disposition: form-data; name="table1"; filename="targets.vot"\r\nContent-Type: application/x-votable+xml\r\n\r\n'.encode()+(a.work/'targets.vot').read_bytes()+b'\r\n')
    parts.append(f'--{boundary}--\r\n'.encode())
    status=dict(endpoint=ENDPOINT,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),n_candidates=len(req),
        request_sha256=hashlib.sha256(a.requests.read_bytes()).hexdigest(),query_sha256=hashlib.sha256(QUERY.encode()).hexdigest(),
        science_measurement=False,psf_approved=False,epoch_propagation_applied=False,match_radius_arcsec=2.,status='requested')
    try:
        endpoint=ENDPOINT.replace('/sync','/async') if a.async_query else ENDPOINT
        status['endpoint']=endpoint
        request=urllib.request.Request(endpoint,b''.join(parts),{'Content-Type':'multipart/form-data; boundary='+boundary})
        with urllib.request.urlopen(request,timeout=180) as response:
            raw=response.read(20*1024*1024+1);job_url=response.geturl().rstrip('/')
        if a.async_query:
            if not job_url.startswith(endpoint+'/'):raise ValueError('Missing async job location')
            status.update(job_url=job_url,status='submitted')
            (a.work/'status.json').write_text(json.dumps(status,indent=2)+'\n')
            for attempt in range(11):
                with urllib.request.urlopen(job_url+'/phase',timeout=60) as response:phase=response.read().decode().strip()
                print('Gaia job phase: '+phase,flush=True)
                if phase=='COMPLETED':break
                if phase in ('ERROR','ABORTED'):
                    with urllib.request.urlopen(job_url+'/error',timeout=60) as response:(a.work/'tap_error.txt').write_bytes(response.read(100000))
                    raise ValueError('Gaia asynchronous query '+phase)
                if attempt==10:raise TimeoutError('Gaia job pending; saved job URL for follow-up')
                time.sleep(60)
            with urllib.request.urlopen(job_url+'/results/result',timeout=180) as response:raw=response.read(20*1024*1024+1)
        (a.work/'gaia_matches.csv').write_bytes(raw)
        if len(raw)>20*1024*1024:raise ValueError('Response exceeds limit')
        reader=csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
        if not {'candidate_id','source_id','separation_arcsec'}.issubset(reader.fieldnames or []):raise ValueError('Unexpected TAP schema; raw response preserved')
        rows=list(reader)
        if len(rows)>=20000:raise ValueError('Possible MAXREC truncation')
        for row in rows:
            if None in row or any(v is None for v in row.values()):raise ValueError('Malformed CSV row')
            if not 0<=int(row['candidate_id'])<len(req):raise ValueError('Unknown candidate ID')
            if not 0<=float(row['separation_arcsec'])<=2.00001:raise ValueError('Match outside search radius')
        status.update(status='complete',n_match_rows=len(rows),response_sha256=hashlib.sha256(raw).hexdigest())
    except Exception as exc:
        if isinstance(exc,urllib.error.HTTPError):(a.work/'http_error.txt').write_bytes(exc.read(100000))
        status.update(status='failed',error_type=type(exc).__name__,error=str(exc));raise
    finally:
        (a.work/'status.json').write_text(json.dumps(status,indent=2)+'\n')
    print(f'Gaia anonymous crossmatch: {len(req)} candidates, {len(rows)} rows; epoch and PSF review pending')


if __name__=='__main__':main()
