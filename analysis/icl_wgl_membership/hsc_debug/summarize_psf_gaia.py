#!/usr/bin/env python3
"""Conservative Gaia evidence labels for review, never a certified PSF selection."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path


def number(row,key):
    try:
        value=float(row.get(key,''));return value if math.isfinite(value) else None
    except (TypeError,ValueError):return None


def classify(rows):
    if not rows:return 'no_gaia_match',[]
    if len(rows)!=1:return 'ambiguous_position',[]
    row=rows[0];sep=number(row,'separation_arcsec')
    if sep is None or not 0<=sep<=1.:return 'offset_requires_epoch_review',[]
    ruwe=number(row,'ruwe');prob=number(row,'classprob_dsc_combmod_star')
    evidence=[]
    parallax=number(row,'parallax_over_error')
    if parallax is not None and parallax>=5:evidence.append('positive_parallax_snr_ge5')
    for axis in ('pmra','pmdec'):
        pm=number(row,axis);err=number(row,axis+'_error')
        if pm is not None and err is not None and err>0 and abs(pm)/err>=5:evidence.append(axis+'_component_snr_ge5')
    flags=all(str(row.get(k,'')).lower() in ('false','0') for k in ('duplicated_source','in_qso_candidates','in_galaxy_candidates'))
    good=ruwe is not None and 0<ruwe<=1.4 and prob is not None and .99<=prob<=1 and flags
    if good and evidence:return 'stellar_evidence_requires_epoch_and_psf_review',evidence
    return 'matched_but_stellar_evidence_insufficient',evidence


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    status=json.loads((a.work/'status.json').read_text());assert status['status']=='complete'
    raw=(a.work/'gaia_matches.csv').read_bytes();assert hashlib.sha256(raw).hexdigest()==status['response_sha256']
    req=json.loads((a.work/'candidates.json').read_text());group=collections.defaultdict(list)
    with (a.work/'gaia_matches.csv').open(newline='') as f:
        for row in csv.DictReader(f):group[int(row['candidate_id'])].append(row)
    records=[]
    for i,r in enumerate(req):
        label,evidence=classify(group[i]);records.append(dict(candidate_id=i,cluster_id=r['cluster_id'],number=r['number'],
            status=label,source_ids=[x['source_id'] for x in group[i]],separations_arcsec=[float(x['separation_arcsec']) for x in group[i]],
            evidence=evidence,epoch_review_required=True,psf_approved=False))
    source_counts=collections.Counter(source for r in records for source in r['source_ids'])
    summary=dict(n_candidates=len(req),status_counts=dict(collections.Counter(r['status'] for r in records)),
        n_unique_gaia_sources=len(source_counts),n_sources_shared_by_candidates=sum(n>1 for n in source_counts.values()),
        n_clusters_with_stellar_evidence=len({r['cluster_id'] for r in records if r['status'].startswith('stellar_evidence')}),
        science_measurement=False,psf_approved=False,records=records,
        thresholds=dict(unique_within_arcsec=2,accepted_position_arcsec=1,ruwe_max=1.4,star_probability_min=.99,astrometric_component_snr_min=5),
        caveats=['Engineering prioritization, not a calibrated stellar purity/completeness model.',
        'Coordinates at Gaia reference epoch; HSC coadd effective epoch not determined.',
        'No-match does not imply galaxy; motion and Gaia depth can cause missing matches.',
        'PSF isolation, variability, saturation and wide wings still require image review.'])
    with (a.work/'evidence_summary.json').open('x') as f:json.dump(summary,f,indent=2)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('records','caveats')},indent=2))


if __name__=='__main__':main()
