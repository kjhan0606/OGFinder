#!/usr/bin/env python3
"""Describe missing point-z inputs in the existing S16A pilot; no calibration."""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

from audit_s16a_photoz_radial_support import PLAN, OBJECT_QA, load_inputs, sha, catalogue_reader, write_table
from s16a_y1_source_contract import finite


def distribution(values):
    x = np.asarray([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    return dict(n=len(x), median=float(np.median(x)) if len(x) else None,
                p16=float(np.percentile(x, 16)) if len(x) else None,
                p84=float(np.percentile(x, 84)) if len(x) else None)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    receipt = json.loads((PLAN/'results'/'download_summary.json').read_text())
    properties = {}
    for name, digest in receipt['output_sha256'].items():
        path = PLAN/'results'/name
        if sha(path) != digest:
            raise ValueError('supplement digest mismatch')
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                obj = row['object_id']
                if obj in properties:
                    raise ValueError('duplicate object')
                mag, extinction = finite(row['i_mag']), finite(row['a_i'])
                flux, error = finite(row['i_flux']), finite(row['i_flux_err'])
                properties[obj] = (mag-extinction if mag is not None and extinction is not None else None,
                                   flux/error if flux is not None and error is not None and error > 0 else None)
    if not selected <= properties.keys() or not available <= properties.keys():
        raise ValueError('selected or available object absent')
    c_map = {c['center_id']:c for c in centers}
    groups = defaultdict(list)
    unique_weights = {}
    for p in pairs:
        oid = p['object_id']
        if oid not in selected:
            continue
        weight = float(p['weight'])
        if not np.isfinite(weight) or weight <= 0:
            raise ValueError('invalid shape weight')
        if oid in unique_weights and not np.isclose(unique_weights[oid], weight, rtol=1e-12):
            raise ValueError('object shape weight changes between centers')
        unique_weights[oid] = weight
        groups[p['center_id']].append((oid, weight))
    if set(unique_weights) != selected:
        raise ValueError('selected object ledger differs')
    rows = []
    for c in centers:
        data = groups[c['center_id']]
        total = sum(w for _,w in data)
        missing = [(o,w) for o,w in data if o not in available]
        rows.append(dict(center_id=c['center_id'], lens_id=c['matched_lens_id'], role=c['role'],
                         selected_pairs=len(data), missing_pairs=len(missing),
                         missing_pair_fraction=len(missing)/len(data) if data else None,
                         missing_shape_weight_fraction=sum(w for _,w in missing)/total if total else None))
    unique = {}
    for label, ids in [('finite',selected & available), ('missing',selected-available)]:
        unique[label] = dict(n=len(ids), i_mag_extinction_corrected=distribution(properties[o][0] for o in ids),
                             i_snr=distribution(properties[o][1] for o in ids),
                             shape_weight=distribution(unique_weights[o] for o in ids))
    role = {}
    for label in ('lens','random'):
        members = [r for r in rows if r['role']==label]
        entries = [x for c in centers if c['role']==label for x in groups[c['center_id']]]
        weight = sum(w for _,w in entries)
        role[label] = dict(selected_pairs=len(entries), missing_pairs=sum(o not in available for o,_ in entries),
                           missing_shape_weight_fraction=sum(w for o,w in entries if o not in available)/weight if weight else None,
                           nonempty_centers=sum(bool(r['selected_pairs']) for r in members))
    output.mkdir(parents=True)
    table = output/'center_missingness.tsv'
    write_table(table, rows)
    report = dict(science_measurement=False, calibration_applied=False, source_object_counts=unique,
                  pair_counts_by_role=role, selected_objects=len(selected),
                  limitations=['Missing photo-z values have no lensing efficiency, so missing shape weight is not missing lensing-estimator weight.',
                    'Magnitude and S/N contrasts are descriptive and do not prove missingness causes or transportability.',
                    'COSMOS lacks fields needed to replay all local quality, mask and missingness criteria; no calibration shift is identifiable.',
                    'Source reuse makes pairs dependent; no significance test or correction is applied.'],
                  provenance=provenance, runner_sha256=sha(Path(__file__)), output_sha256={table.name:sha(table)})
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',required=True,type=Path)
    args=parser.parse_args()
    run(args.output.resolve())
