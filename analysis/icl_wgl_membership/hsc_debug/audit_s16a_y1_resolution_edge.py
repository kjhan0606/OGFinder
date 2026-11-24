#!/usr/bin/env python3
"""Measure shape-weighted resolution-cut edge density in the fixed S16A pilot."""
import argparse
from collections import defaultdict
import json
from pathlib import Path

import numpy as np

from audit_s16a_photoz_radial_support import (
    FlatLambdaCDM, OBJECT_QA, PLAN, RADIAL_EDGES, catalogue_reader, load_inputs,
    radial_bins, separations_radians, sha, write_table,
)

EDGE = .3
WIDTHS = (.01, .02, .05)


def load_resolution():
    root = PLAN/'results'
    receipt = json.loads((root/'download_summary.json').read_text())
    resolutions = {}
    for name, digest in receipt['output_sha256'].items():
        path = root/name
        if sha(path) != digest:
            raise ValueError(f'changed supplement file {name}')
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                ident = row['object_id']
                if ident in resolutions:
                    raise ValueError('duplicate resolution entry')
                resolutions[ident] = row['resolution']
    return resolutions


def metrics(rows, role, sample, annulus):
    weights = np.array([r[0] for r in rows])
    resolution = np.array([r[1] for r in rows])
    if (not np.all(np.isfinite(weights) & (weights > 0))
            or not np.all(np.isfinite(resolution) & (resolution >= EDGE))):
        raise ValueError('selected sample violates the expected weight or R2 domain')
    total = float(weights.sum())
    for width in WIDTHS:
        near = resolution < round(EDGE + width, 2)
        weighted = float(weights[near].sum())
        density = weighted / total / width if total else None
        yield dict(role=role, sample=sample, annulus=annulus, width=width,
                   n_pairs=len(rows), n_edge_pairs=int(near.sum()), shape_weight=total,
                   edge_shape_weight=weighted, edge_density=density)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, prior, provenance = load_inputs(PLAN, OBJECT_QA)
    resolutions = load_resolution()
    if set(resolutions) != {r['object_id'] for r in pairs}:
        raise ValueError('resolution ledger does not match the pilot object IDs')
    by_center = defaultdict(list)
    for row in pairs:
        if row['object_id'] in selected:
            by_center[row['center_id']].append(row)
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    summaries = defaultdict(list)
    expected = {(r['center_id'], int(r['radial_bin'])): int(r['selected_count']) for r in prior}
    availability = PLAN.parent/'s16a_photoz_radial_support_2026-11-18_01'/'summary.json'
    if sha(availability) != '89e0f8728fa8d1097ac77d107e0401f05892a20c2bb07dd7eeec7a7366c88a75':
        raise ValueError('prior annular availability check changed')
    expected_roles = json.loads(availability.read_text())['role_pair_diagnostics']
    for center in centers:
        ident, role = center['center_id'], center['role']
        rows = by_center[ident]
        radius = separations_radians([float(r['source_ra']) for r in rows],
            [float(r['source_dec']) for r in rows], center['ra_deg'], center['dec_deg'])
        bins = radial_bins(radius * cosmo.angular_diameter_distance(center['z_lens']).value)
        for k, n in enumerate(np.bincount(bins, minlength=len(RADIAL_EDGES)-1)):
            if n != expected[ident, k]:
                raise ValueError('annular support differs from previous pilot QA')
        for row, k in zip(rows, bins, strict=True):
            obj = row['object_id']
            entry = (float(row['weight']), float(resolutions[obj]))
            for annulus in ('all', f'{RADIAL_EDGES[k]:g}–{RADIAL_EDGES[k+1]:g}'):
                summaries[role, 'table4', annulus].append(entry)
                if obj in available:
                    summaries[role, 'finite_photoz', annulus].append(entry)
    for role, values in expected_roles.items():
        if (len(summaries[role, 'table4', 'all']) != values['selected_pairs']
                or len(summaries[role, 'finite_photoz', 'all']) != values['finite_photoz_pairs']):
            raise ValueError('overall pair counts differ from prior availability diagnostic')
    table = [row for (role, sample, annulus), values in sorted(summaries.items())
             for row in metrics(values, role, sample, annulus)]
    output.mkdir(parents=True, exist_ok=False)
    path = output/'resolution_edge_density.tsv'
    write_table(path, table)
    result = dict(science_measurement=False, **provenance, edge=EDGE, widths=WIDTHS,
                  overall=[r for r in table if r['annulus']=='all'], n_summary_rows=len(table),
                  weight_definition='Catalog shape weight per source-center pair',
                  limitation='One-sided finite-bin density sensitivity. Finite photo-z selects available inputs only. Final lensing pair weights and calibrated m_sel are pending.',
                  radial_availability_summary_sha256=sha(availability), output_sha256={path.name:sha(path)},
                  runner_sha256=sha(Path(__file__)))
    (output/'summary.json').write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['overall'], indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
