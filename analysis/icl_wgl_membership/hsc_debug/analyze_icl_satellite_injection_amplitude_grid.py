#!/usr/bin/env python3
"""Paired amplitude-linearity and mask-response audit for ICL/satellite injections."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

FIELDS = {
    'satellite_unit_response': ('transfer_response', 'satellite_response_to_satellite_unit_injection'),
    'satellite_envelopes_to_icl': ('transfer_response', 'icl_response_to_all_satellite_envelopes'),
    'icl_diffuse_unit_response': ('transfer_response', 'icl_response_to_diffuse_unit_injection'),
    'diffuse_to_satellite_response': ('transfer_response', 'satellite_response_to_smooth_diffuse'),
}


def get(d, path):
    for k in path:
        d = d[k]
    return d


def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1<<20),b''): h.update(block)
    return h.hexdigest()


def load(path, expected_amp=None):
    out={}
    for f in sorted(path.glob('cluster_*.json')):
        d=json.loads(f.read_text())
        if d.get('status')!='complete' or d.get('science_measurement') is not False:
            raise ValueError(f'unexpected result status/scope {f}')
        amp=d['injection']['satellite_amplitude_at_50_per_galaxy_njy_arcsec2']
        if expected_amp is not None and not np.isclose(amp,expected_amp):
            raise ValueError(f'{f}: expected amplitude {expected_amp}, got {amp}')
        key=(str(d['cluster_id']),d['filter'])
        if key in out: raise ValueError(f'duplicate frame {key}')
        out[key]=(f,d)
    return out


def qsum(vals):
    x=np.asarray([v for v in vals if v is not None and np.isfinite(float(v))],float)
    if not len(x):
        return {'n':0,'median':None,'q16':None,'q84':None,'mean':None}
    return {'n':len(x),'median':float(np.median(x)),'q16':float(np.quantile(x,.16)),
            'q84':float(np.quantile(x,.84)),'mean':float(np.mean(x))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--amp-0p5',type=Path,required=True)
    p.add_argument('--amp-1p0',type=Path,required=True)
    p.add_argument('--amp-2p0',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--expected-frames',type=int,default=57,
                   help='required paired-frame count (defaults to the original 57-frame validation sample)')
    p.add_argument('--bootstrap',type=int,default=10000)
    a=p.parse_args()
    if a.expected_frames < 1 or a.bootstrap < 1:
        raise ValueError('expected frame and bootstrap counts must be positive')
    runs={.5:load(a.amp_0p5,.5),1.:load(a.amp_1p0,1.),2.:load(a.amp_2p0,2.)}
    keys=set(runs[.5])
    if any(set(v)!=keys for v in runs.values()) or len(keys)!=a.expected_frames:
        raise ValueError(f'expected identical paired {a.expected_frames}-frame sets at all amplitudes')
    keys=sorted(keys)
    records=[]
    for key in keys:
        samples={amp:runs[amp][key][1] for amp in runs}
        ref=samples[1.]
        for amp,d in samples.items():
            for f in ('filter','z','icl_profile_bin_kpc','satellite_profile_bin_kpc','cluster_radius_selection_kpc'):
                if d.get(f)!=ref.get(f): raise ValueError(f'{f} mismatch for {key}')
            inj=d['injection']; inj_ref=ref['injection']
            for f in ('diffuse_amplitude_at_100_njy_arcsec2','diffuse_scale_kpc',
                      'satellite_scale_kpc','n_injected_satellite_centers','psf_file_count',
                      'psf_model','unresolved_bcg_policy','source_detection_and_masks_rerun'):
                if inj.get(f)!=inj_ref.get(f): raise ValueError(f'injection config mismatch {f} for {key}')
            imghash=lambda x: next((v for n,v in x['input_sha256'].items() if n.endswith(f'/images/cluster_{key[0]}.fits')),None)
            psfhash=lambda x: sorted(v for n,v in x['input_sha256'].items() if '/psf_' in n and n.endswith('_HSC-I.fits'))
            if imghash(d)!=imghash(ref) or psfhash(d)!=psfhash(ref):
                raise ValueError(f'image or position PSF content changed for {key}')
            if d['redetection_by_scenario']['baseline']!=ref['redetection_by_scenario']['baseline']:
                raise ValueError(f'no-injection baseline re-detection changed for {key}')
        row={'cluster_id':key[0],'filter':key[1],'z':ref['z'],
             'image_sha256':imghash(ref),'n_position_psfs':len(psfhash(ref))}
        for amp,d in samples.items():
            tag=str(amp).replace('.','p')
            for metric,path in FIELDS.items(): row[f'{metric}_{tag}']=get(d,path)
            base=d['redetection_by_scenario']['baseline']
            sat=d['redetection_by_scenario']['satellite_only']
            both=d['redetection_by_scenario']['both']
            row[f'detections_satellite_only_minus_baseline_{tag}']=sat['n_detections']-base['n_detections']
            row[f'mask_pixels_satellite_only_minus_baseline_{tag}']=sat['icl_mask_pixels']-base['icl_mask_pixels']
            row[f'detections_both_minus_baseline_{tag}']=both['n_detections']-base['n_detections']
            row[f'mask_pixels_both_minus_baseline_{tag}']=both['icl_mask_pixels']-base['icl_mask_pixels']
        records.append(row)
    if a.output_dir.exists(): raise FileExistsError(a.output_dir)
    a.output_dir.mkdir(parents=True)
    pair=a.output_dir/'paired_amplitude_frames.tsv'
    with pair.open('w',newline='') as f:
        w=csv.DictWriter(f,list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
    rng=np.random.default_rng(20261105)
    summary={'schema':'icl-satellite-amplitude-grid-audit-v1','science_measurement':False,
             'status':'complete','interpretation':'injection/recovery transfer sensitivity; not detection or completeness correction',
             'n_paired_frames':len(records),'n_cluster_ids':len({x['cluster_id'] for x in records}),
             'frames_by_filter':{band:sum(r['filter']==band for r in records) for band in sorted({r['filter'] for r in records})},
             'injected_satellite_amplitudes_njy_arcsec2':[.5,1.,2.],
             'fixed_diffuse_amplitude_njy_arcsec2':1.,'image_and_position_psf_hashes_identical':True,
             'no_injection_baseline_detection_masks_identical':True,'metrics_by_filter':{}}
    for band in [*sorted({r['filter'] for r in records}),'all']:
        rr=[r for r in records if band=='all' or r['filter']==band]
        summary['metrics_by_filter'][band]={}
        for metric in FIELDS:
            summary['metrics_by_filter'][band][metric]={}
            for amp in (.5,1.,2.):
                tag=str(amp).replace('.','p')
                summary['metrics_by_filter'][band][metric][str(amp)]=qsum([r[f'{metric}_{tag}'] for r in rr])
            summary['metrics_by_filter'][band][metric]['paired_delta_vs_1']={}
            for amp in (.5,2.):
                tag=str(amp).replace('.','p')
                matched={}
                for cid in sorted({r['cluster_id'] for r in rr}):
                    values=[]
                    for row in rr:
                        if row['cluster_id']!=cid:
                            continue
                        high=row[f'{metric}_{tag}']; reference=row[f'{metric}_1p0']
                        if (high is not None and reference is not None
                                and np.isfinite(float(high)) and np.isfinite(float(reference))):
                            values.append(float(high)-float(reference))
                    if values:
                        matched[cid]=np.asarray(values,float)
                ids=sorted(matched)
                point=np.concatenate([matched[cid] for cid in ids]) if ids else np.asarray([],float)
                boots=[]
                if ids:
                    for _ in range(a.bootstrap):
                        draws=rng.choice(ids,len(ids),replace=True)
                        boots.append(float(np.median(np.concatenate([matched[c] for c in draws]))))
                summary['metrics_by_filter'][band][metric]['paired_delta_vs_1'][str(amp)]={
                    'n_paired_frames':int(len(point)),
                    'n_paired_clusters':int(len(ids)),
                    'median':float(np.median(point)) if len(point) else None,
                    'cluster_bootstrap_95pct_ci':([float(np.quantile(boots,.025)),float(np.quantile(boots,.975))]
                                                   if boots else [None,None])}
        summary['metrics_by_filter'][band]['mask_and_detection_changes']={}
        for amp in (.5,1.,2.):
            tag=str(amp).replace('.','p')
            summary['metrics_by_filter'][band]['mask_and_detection_changes'][str(amp)]={}
            for metric in ('detections_satellite_only_minus_baseline','mask_pixels_satellite_only_minus_baseline',
                           'detections_both_minus_baseline','mask_pixels_both_minus_baseline'):
                vals=[r[f'{metric}_{tag}'] for r in rr]
                summary['metrics_by_filter'][band]['mask_and_detection_changes'][str(amp)][metric]={
                    **qsum(vals),'n_frames_changed':sum(v!=0 for v in vals)}
    summary['input_json_sha256']={str(amp):{f'cluster_{cid}_{band}.json':sha(runs[amp][(cid,band)][0]) for cid,band in keys} for amp in runs}
    (a.output_dir/'summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    lines=['# Satellite-injection amplitude response audit','',
           f'This paired offline analysis compares 0.5, 1, and 2 nJy/arcsec² satellite injections on the same {len(records)} HSC frames, images, and position-dependent PSFs. The diffuse injection is fixed at 1 nJy/arcsec². Detection and masks were regenerated in each scenario. All products retain `science_measurement=false`.', '',
           '| Filter | Amplitude | Unit satellite recovery median (q16–q84) | Frames with changed detection count (satellite only / both) |', '|---|---:|---:|---:|']
    for band in sorted({r['filter'] for r in records}):
        v=summary['metrics_by_filter'][band]
        for amp in ('0.5','1.0','2.0'):
            unit=v['satellite_unit_response'][amp]; m=v['mask_and_detection_changes'][amp]
            lines.append(f"| {band} | {amp} | {unit['median']:.3f} ({unit['q16']:.3f}–{unit['q84']:.3f}) | {m['detections_satellite_only_minus_baseline']['n_frames_changed']} / {m['detections_both_minus_baseline']['n_frames_changed']} |")
    lines += ['', 'The unit-response coefficients are normalized to the injected amplitude. Their frame-to-frame spread is descriptive; paired cluster bootstrap intervals are in `summary.json`. Changes in the regenerated masks can alter the measurement even when the source image and PSFs are identical. This experiment assesses transfer nonlinearity and mask sensitivity only; it does not infer a real satellite population or correct completeness.', '']
    (a.output_dir/'assessment.md').write_text('\n'.join(lines))
    print(json.dumps({'output_dir':str(a.output_dir),'paired_frames':len(records),
                      'summary_sha256':sha(a.output_dir/'summary.json'),
                      'table_sha256':sha(pair),'assessment_sha256':sha(a.output_dir/'assessment.md')},indent=2))

if __name__=='__main__': main()
