#!/usr/bin/env python3
"""Measure the catalog ensemble response under pilot shape weights only."""
import argparse
from collections import defaultdict
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wl_y1 import y1_shape_response
from audit_s16a_photoz_radial_support import OBJECT_QA, PLAN, load_inputs, sha, write_table


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    roles = {c['center_id']: c['role'] for c in centers}
    groups = defaultdict(list)
    for pair in pairs:
        obj = pair['object_id']
        if obj not in selected:
            continue
        group = roles[pair['center_id']]
        values = tuple(float(pair[key]) for key in ('weight', 'rms_e', 'm'))
        groups[group, 'table4'].append(values)
        if obj in available:
            groups[group, 'finite_photoz'].append(values)
    rows = [dict(role=role, sample=sample, **y1_shape_response(*zip(*values, strict=True)),
                 individual_m_le_minus_one=sum(row[2] <= -1 for row in values))
            for (role, sample), values in sorted(groups.items())]
    output.mkdir(parents=True, exist_ok=False)
    table = output/'shape_response.tsv'
    write_table(table, rows)
    record = dict(**provenance, results=rows, weight_definition='catalog shape weight per source-center pair',
                  science_measurement=False, delta_sigma_computed=False,
                  limitations='No shear, distance efficiency, redshift responsivity or selection correction applied.',
                  output_sha256={table.name:sha(table)}, runner_sha256=sha(Path(__file__)),
                  kernel_sha256=sha(Path(__file__).resolve().parents[3]/'ogfkit/wl_y1.py'))
    with (output/'summary.json').open('x') as stream:
        json.dump(record, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(rows, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
