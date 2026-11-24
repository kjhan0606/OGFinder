#!/usr/bin/env python3
"""Test point-z gap sensitivity before applying calibrated Y1 source cuts."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
import sys

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit.wlcluster import C_KMS, G_MPC
from audit_s16a_photoz_radial_support import (
    FlatLambdaCDM, OBJECT_QA, PLAN, RADIAL_EDGES, catalogue_reader, load_inputs, sha, write_table,
)
from audit_s16a_y1_reuse_shape_noise import GEOMETRY, GEOMETRY_SHA

GAPS = (0., .1, .2)


def finite_point_redshifts(available):
    results = PLAN/'results'
    receipt = json.loads((results/'download_summary.json').read_text())
    result = {}
    for name, digest in receipt['output_sha256'].items():
        path = results/name
        if sha(path) != digest:
            raise ValueError('point-z input changed')
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                obj = row['object_id']
                if obj in available:
                    result[obj] = float(row['photoz_best'])
    if result.keys() != available or not all(np.isfinite(z) and z >= 0 for z in result.values()):
        raise ValueError('point-z availability ledger differs or contains invalid redshifts')
    return result


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    center_lookup = {c['center_id']:c for c in centers}
    point_z = finite_point_redshifts(available)
    geometry_file = GEOMETRY/'summary.json'
    if sha(geometry_file) != GEOMETRY_SHA:
        raise ValueError('changed source geometry summary')
    geometry = json.loads(geometry_file.read_text())
    table = GEOMETRY/'pair_geometry.tsv'
    if any(geometry[k] != v for k,v in provenance.items()) or sha(table) != geometry['output_sha256'][table.name]:
        raise ValueError('source geometry lineage differs')

    groups = defaultdict(list)
    with table.open(newline='') as stream:
        for pair, geo in zip(pairs, csv.DictReader(stream, delimiter='\t'), strict=True):
            obj, cid = pair['object_id'], pair['center_id']
            if (geo['object_id'], geo['center_id']) != (obj,cid):
                raise ValueError('pair identities differ')
            if obj in selected and obj in point_z:
                groups[cid].append((int(geo['radial_bin']), point_z[obj], float(pair['weight'])))

    cosmology = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    conversion = 4*np.pi*G_MPC/C_KMS**2
    n_bins = len(RADIAL_EDGES)-1
    counts, weights_by_role, results = {}, defaultdict(list), []
    for cid, center in center_lookup.items():
        values = np.asarray(groups[cid], dtype=float).reshape(-1,3)
        bins, z, shape_weight = values.T
        lens_z = float(center['z_lens'])
        chi_lens = cosmology.comoving_distance(lens_z).value
        chi_source = cosmology.comoving_distance(z).value
        beta = np.zeros_like(z)
        behind = z > lens_z
        beta[behind] = conversion * chi_lens/(1+lens_z)*(1-chi_lens/chi_source[behind])
        w_pair = shape_weight*beta**2
        if not np.all(np.isfinite(w_pair) & (w_pair >= 0)):
            raise ValueError('invalid diagnostic lensing pair weight')
        for gap in GAPS:
            keep = z > lens_z+gap
            weights_by_role[center['role'],gap].extend(w_pair[keep])
            for k in range(n_bins):
                in_bin = bins == k
                weights = w_pair[in_bin & keep]
                total = float(weights.sum())
                n_eff = total**2/float(weights @ weights) if total else 0.
                counts[cid,k,gap] = len(weights)
                results.append(dict(center_id=cid, lens_id=center['matched_lens_id'], role=center['role'],
                                    annulus=k, gap=gap, available_pairs=int(in_bin.sum()),
                                    candidate_pairs=len(weights), pair_weight_sum=total, pair_neff=n_eff))

    if sum(r['available_pairs'] for r in results if r['gap'] == 0) != 161036:
        raise ValueError('point-z availability no longer matches verified pilot')
    lens_groups = defaultdict(list)
    for c in centers:
        lens_groups[c['matched_lens_id']].append(c['center_id'])
    membership = {str(gap):[sorted(lens_id for lens_id,cids in lens_groups.items() if all(counts[cid,k,gap] for cid in cids))
                            for k in range(n_bins)] for gap in GAPS}
    common = {gap:sorted(set.intersection(*map(set,bins[1:]))) for gap,bins in membership.items()}
    totals = []
    for (role,gap), weights in sorted(weights_by_role.items()):
        w = np.asarray(weights)
        base = sum(weights_by_role[role,0.])
        totals.append(dict(role=role, gap=gap, n_pairs=len(w), pair_weight_fraction=float(w.sum()/base),
                           pair_neff=float(w.sum()**2/(w@w))))
    output.mkdir(parents=True, exist_ok=False)
    target = output/'center_annulus_pointz_support.tsv'
    write_table(target, results)
    summary = dict(**provenance, geometry_sha256=GEOMETRY_SHA, radial_edges_proper_mpc=RADIAL_EDGES.tolist(),
                   role_totals=totals, per_annulus_common_lens_ids=membership, common_0p1_to_3=common,
                   science_measurement=False, measured_shapes_used=False,
                   definition='Table 4 pairs with finite best/risk and z_best strictly above lens_z + gap. No risk threshold. Pair weights are catalog shape weight times squared proper inverse critical density.',
                   limitations='Point estimate sensitivity only. Pair N_eff does not account for source reuse. Calibrated background selection is pending.',
                   output_sha256={target.name:sha(target)}, runner_sha256=sha(Path(__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(totals,indent=2))
    print(json.dumps(common,indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
