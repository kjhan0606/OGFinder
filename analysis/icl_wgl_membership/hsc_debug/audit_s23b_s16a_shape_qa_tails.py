#!/usr/bin/env python3
"""Flag post-intake tails in already received fields; never apply source cuts.

These are post hoc diagnostics, not Y1 selection criteria or per-object shear
validity rules. The derived sigma_e is not the unreceived raw HSM sigma field.
"""
from collections import Counter, defaultdict
import json
import math
from pathlib import Path

from fetch_ancillary_catalogues import catalogue_reader, sha
from submit_s23b_s16a_shape_qa_pilot import numeric
from prepare_s23b_s16a_shape_qa_schema_fix import OUTPUT as PLAN

OUTPUT=PLAN.parent/'wl_s23b_s16a_shape_qa_pilot_evidence_2026-11-17_01'


def flags(m,e1,e2,sigma):
    return [name for name,active in (
        ('m_le_minus_one',math.isfinite(m) and m<=-1),
        ('ellipticity_norm_ge_two',math.isfinite(e1) and math.isfinite(e2) and math.hypot(e1,e2)>=2),
        ('derived_sigma_e_gt_0p4',math.isfinite(sigma) and sigma>.4)) if active]


def run():
    results=PLAN/'results'
    receipt=json.loads((results/'download_summary.json').read_text())
    counts=Counter();objects=defaultdict(set);by_center=defaultdict(Counter)
    flagged=[];rows=0
    for name,digest in receipt['output_sha256'].items():
        file=results/name
        if sha(file)!=digest:raise ValueError('received CSV changed')
        with file.open(newline='',encoding='utf-8-sig') as f:
            for row in catalogue_reader(f):
                rows+=1
                m,e1,e2,sigma=(numeric(row[k]) for k in ('m','e1','e2','sigma_e'))
                labels=flags(m,e1,e2,sigma)
                if labels:
                    flagged.append(dict(center_id=row['center_id'],object_id=row['object_id'],
                        diagnostic_flags=','.join(labels),m=m,e_norm=math.hypot(e1,e2),derived_sigma_e=sigma))
                for flag in labels:
                    counts[flag]+=1;objects[flag].add(row['object_id']);by_center[row['center_id']][flag]+=1
    if rows!=169535:raise ValueError('tail audit pair count differs')
    output=OUTPUT/'tail_check';output.mkdir(parents=True,exist_ok=False)
    from audit_s23b_s16a_wl_random_support import write_table
    if flagged:write_table(output/'flagged_source_pairs.tsv',flagged)
    result=dict(schema='astrafek-s16a-pilot-post-intake-tail-diagnostics-v1',
        source_pair_rows=rows,flagged_pair_union=len(flagged),
        flags={name:dict(source_pairs=counts[name],unique_objects=len(objects[name]))
            for name in ('m_le_minus_one','ellipticity_norm_ge_two','derived_sigma_e_gt_0p4')},
        by_center={k:dict(v) for k,v in sorted(by_center.items())},
        source_csv_sha256=receipt['output_sha256'],script_sha256=sha(Path(__file__)),
        source_qa_summary_sha256=sha(PLAN.parent/'wl_s23b_s16a_shape_qa_pilot_intake_2026-11-17_01/summary.json'),
        used_as_source_selection=False,source_rows_deleted_or_repaired=False,
        network_used=False,science_measurement=False,
        limitations=['Post hoc descriptive tails; no threshold is declared a complete Y1 selection.',
            'Per-object m<=-1 flags a tail for review, not an automatic invalidity rule for the ensemble m correction.',
            'The exported derived sigma_e differs from the raw HSM sigma used in source-selection recipes; .4 is only a diagnostic marker here.',
            'No shear, physical validity, or detection significance is inferred.'])
    (output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('by_center','source_csv_sha256')},indent=2))


if __name__=='__main__':run()
