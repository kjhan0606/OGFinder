#!/usr/bin/env python3
"""Compare audited member PSF-core vs unconvolved injections on common support.

No downloads or image reads. The comparison is conditional on both variants
passing profile selection, with failures counted explicitly. Ratios have their
own convolved/unconvolved nominal annular input denominators; recovered flux
increments are also compared in physical surface-brightness units. Intervals
resample target sky cells, not individual injection rows.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import numpy as np

PLAN_SHA='03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d'
OLD_SUMMARY_SHA='e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77'
BASE=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
OLD=BASE/'s23b_camira_member_outer_injection_recovery_2026-11-07_01'
TARGET=BASE/'s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json'
TARGET_SHA='e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44'


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def read(p):
    with Path(p).open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def key(r):return (r['target_id'],r['filter'],r['member_index'],r['source_row'],float(r['amplitude_at_50_kpc_njy_arcsec2']))


def paired_rows(old,new):
    a={key(r):r for r in old};b={key(r):r for r in new}
    if len(a)!=len(old) or len(b)!=len(new) or set(a)!=set(b):
        raise ValueError('injection member/amplitude rows do not match exactly')
    out=[]
    for k in sorted(a):
        x,y=a[k],b[k]
        if y.get('psf_convolution_applied')!='True' or not y.get('member_psf_sha256'):
            raise ValueError('PSF row lacks convolved flag or PSF provenance')
        if abs(float(x['baseline_delta_sb_njy_arcsec2'])-float(y['baseline_delta_sb_njy_arcsec2']))>2e-5:
            raise ValueError('baseline replay changed between variants')
        supported=x['injected_status']=='supported' and y['injected_status']=='supported'
        row=dict(target_id=k[0],filter=k[1],member_index=k[2],source_row=k[3],amplitude=k[4],
                 unconvolved_status=x['injected_status'],psf_status=y['injected_status'],both_supported=supported,
                 unconvolved_response=None,psf_response=None,delta_response=None,
                 unconvolved_increment=None,psf_increment=None,delta_increment=None)
        if supported:
            for label,column in [('response','recovered_increment_over_nominal_input'),('increment','recovered_delta_increment_njy_arcsec2')]:
                v0,v1=float(x[column]),float(y[column])
                if not math.isfinite(v0) or not math.isfinite(v1):raise ValueError('nonfinite recovered value')
                row['unconvolved_'+label]=v0;row['psf_'+label]=v1;row['delta_'+label]=v1-v0
        out.append(row)
    return out


def median_block_interval(rows,metric,frames,replicas,seed):
    groups={}
    for r in rows:
        f=frames[r['target_id']]
        cell=(math.floor(((float(f['ra_deg'])-84.056247280625)%360)/5),
              math.floor((float(f['dec_deg'])+90)/5))
        groups.setdefault(cell,[]).append(float(r[metric]))
    cells=sorted(groups); values=np.array([r[metric] for r in rows],dtype=float)
    if len(cells)<5:return dict(n_cells=len(cells),median=float(np.median(values)) if len(values) else None,ci_low=None,ci_high=None,status='insufficient_sky_cells')
    rng=np.random.default_rng(seed);draws=[]
    for _ in range(replicas):
        sample=np.concatenate([groups[cells[int(i)]] for i in rng.integers(0,len(cells),len(cells))])
        draws.append(float(np.median(sample)))
    low,high=np.quantile(draws,[.025,.975])
    return dict(n_cells=len(cells),median=float(np.median(values)),ci_low=float(low),ci_high=float(high),status='complete')


def run(psf_dir,out,replicas):
    if out.exists():raise FileExistsError(out)
    if sha(OLD/'summary.json')!=OLD_SUMMARY_SHA or sha(TARGET)!=TARGET_SHA:raise ValueError('pinned baseline changed')
    oldsum=json.loads((OLD/'summary.json').read_text());newsum=json.loads((psf_dir/'summary.json').read_text())
    for d,s in ((OLD,oldsum),(psf_dir,newsum)):
        if s.get('science_measurement') is not False or s.get('n_frames_complete')!=427 or s.get('n_frames_failed')!=0:
            raise ValueError('incomplete or promoted injection run')
        for filename in ('run_manifest.json','injection_responses.tsv','selected_members.tsv'):
            if sha(d/filename)!=s['output_sha256'][filename]:raise ValueError('injection product hash mismatch')
    method=newsum['method'];psf=method.get('psf_input') or {}
    if (newsum.get('psf_convolution_applied') is not True or psf.get('plan_sha256')!=PLAN_SHA
        or method.get('source_poisson_added') is not False or method.get('variance_updated') is not False):
        raise ValueError('PSF request or noise scope differs from approved calibration')
    if sha(Path(psf['audit']))!=psf['audit_sha256']:raise ValueError('PSF audit changed')
    audit=json.loads(Path(psf['audit']).read_text())
    if audit.get('status')!='verified' or audit.get('approved_plan_sha256')!=PLAN_SHA or audit.get('n_returned')!=427:
        raise ValueError('position-specific PSF receipt not verified')
    rows=paired_rows(read(OLD/'injection_responses.tsv'),read(psf_dir/'injection_responses.tsv'))
    if len(rows)!=1281:raise ValueError('expected 427 x 3 injection rows')
    frames={str(r['id']):r for r in json.loads(TARGET.read_text())['clusters']}
    summaries=[]
    for bi,band in enumerate(('i','i2')):
        for ai,amp in enumerate((.5,1.,2.)):
            subset=[r for r in rows if r['filter']==band and r['amplitude']==amp]
            both=[r for r in subset if r['both_supported']]
            for metric in ('delta_response','delta_increment'):
                result=median_block_interval(both,metric,frames,replicas,20261117+bi*10+ai)
                summaries.append(dict(filter=band,amplitude=amp,metric=metric,n_total=len(subset),n_both=len(both),
                    status_pairs={str(k):v for k,v in Counter((r['unconvolved_status'],r['psf_status']) for r in subset).items()},**result))
    out.mkdir(parents=True,exist_ok=False)
    with (out/'paired_ledger.tsv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    result=dict(status='complete_exploratory_psf_core_paired_comparison',science_measurement=False,network_used=False,images_read=False,
        old_summary_sha256=sha(OLD/'summary.json'),psf_summary_sha256=sha(psf_dir/'summary.json'),
        approved_plan_sha256=PLAN_SHA,psf_audit_sha256=sha(Path(psf['audit'])),target_manifest_sha256=sha(TARGET),
        script_sha256=sha(Path(__file__)),output_sha256={'paired_ledger.tsv':sha(out/'paired_ledger.tsv')},
        bootstrap_replicas=replicas,seed=20261117,comparisons=summaries,
        limitations=['Paired sample includes only profiles supported in both injections; no completeness correction.',
          'Local model-PSF core excludes extended wings.',
          'Ratios use their respective nominal inputs; physical flux increments are reported separately.',
          'Source Poisson, coadd weights, sky transfer and covariance remain uncalibrated.',
          'Marginal intervals across amplitudes are exploratory, not multiplicity-adjusted detection tests.'])
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--psf-dir',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--replicas',type=int,default=5000);a=p.parse_args()
    if a.replicas<100:p.error('at least 100 bootstrap replicates required')
    run(a.psf_dir.resolve(),a.output.resolve(),a.replicas)
