#!/usr/bin/env python3
"""Join received Ephor_AB point estimates to the S16A-matched COSMOS weights.

See the dsigma HSC example for the ID/S17a_objid join. This stage checks
availability and calibration weights before any lens-dependent source cut.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.table import Table

from fetch_ancillary_catalogues import sha
from receive_s16a_ephor_ab_cosmos import CHECKPOINT, DIRECTORY, load_plan, verify_state

POINT_FIELDS = {'photo_z': 'photoz_best', 'photo_z_risk': 'photoz_risk_best',
                'photo_z_low68': 'photoz_err68_min', 'photo_z_std': 'photoz_std_best'}
CALIBRATION_FIELDS = {'reference_z': 'cosmos_photoz', 'shape_weight': 'weight_source', 'som_weight': 'som_weight'}


def columns(table):
    lookup = {name.casefold(): name for name in table.colnames}
    if len(lookup) != len(table.colnames):
        raise ValueError('ambiguous case-insensitive FITS column names')
    return lookup


def integer_ids(column):
    ident = np.ma.asarray(column)
    if ident.ndim != 1 or ident.dtype.kind not in 'iu' or np.any(np.ma.getmaskarray(ident)):
        raise ValueError('join IDs must be exact unmasked integer scalars')
    if (ident <= 0).any() or (ident > np.iinfo(np.int64).max).any():
        raise ValueError('ID outside positive int64 domain')
    return np.asarray(ident, dtype=np.int64)


def numeric(column):
    value = np.ma.asarray(column)
    if value.ndim != 1 or value.dtype.kind not in 'iuf':
        raise ValueError('calibration column must be numeric and scalar')
    return value.astype(float).filled(np.nan)


def require_unique(ids, label):
    if np.unique(ids).size != len(ids):
        raise ValueError(f'duplicate IDs in {label}')


def join_tables(point_tables, reference):
    if len(point_tables) != 2:
        raise ValueError('expected two Ephor_AB point tables')
    point_columns = [columns(table) for table in point_tables]
    ref_columns = columns(reference)
    ids = np.concatenate([integer_ids(t[c['id']]) for t, c in zip(point_tables, point_columns, strict=True)])
    ref_ids = integer_ids(reference[ref_columns['s17a_objid']])
    require_unique(ids, 'both Ephor_AB tracts')
    require_unique(ref_ids, 'COSMOS reference table')
    _, point_index, ref_index = np.intersect1d(ids, ref_ids, assume_unique=True, return_indices=True)
    if not len(point_index):
        raise ValueError('no shared S17A IDs')
    availability = {key: [name in c for c in point_columns] for key, name in POINT_FIELDS.items()}
    if not all(availability['photo_z']):
        raise ValueError('photoz_best must be present in both point tables')
    joined = {'s17a_objid': ref_ids[ref_index], 'tract': np.where(point_index < len(point_tables[0]), 9812, 9813)}
    for key, name in POINT_FIELDS.items():
        joined[key] = np.concatenate([numeric(t[c[name]]) if name in c else np.full(len(t), np.nan)
                                      for t, c in zip(point_tables, point_columns, strict=True)])[point_index]
    for key, name in CALIBRATION_FIELDS.items():
        joined[key] = numeric(reference[ref_columns[name]])[ref_index]

    z_fields = ('photo_z', 'reference_z')
    w_fields = ('shape_weight', 'som_weight')
    valid = {key: np.isfinite(joined[key]) & (joined[key] >= 0 if key in z_fields else joined[key] > 0)
             for key in (*z_fields, *w_fields)}
    usable = np.logical_and.reduce(list(valid.values()))
    combined = joined['shape_weight'][usable]*joined['som_weight'][usable]
    if not combined.size or not np.all(np.isfinite(combined) & (combined > 0)):
        raise ValueError('no valid calibration support or invalid combined weight')
    normalized = combined / combined.max()
    ref_shape = numeric(reference[ref_columns['weight_source']])
    ref_som = numeric(reference[ref_columns['som_weight']])
    ref_weight_ok = np.isfinite(ref_shape) & np.isfinite(ref_som) & (ref_shape > 0) & (ref_som > 0)
    ref_weight = np.zeros(len(reference))
    ref_weight[ref_weight_ok] = ref_shape[ref_weight_ok] * ref_som[ref_weight_ok]
    if not np.all(np.isfinite(ref_weight)):
        raise ValueError('invalid combined weight in COSMOS reference table')
    ref_weight /= ref_weight.max()
    ref_weight_total = ref_weight.sum()
    joined['eligible_for_selection_review'] = usable
    missing = np.setdiff1d(ref_ids, ids, assume_unique=True)
    summary = dict(point_rows=len(ids), reference_rows=len(ref_ids), joined_rows=len(point_index),
                   point_only_rows=len(ids)-len(point_index), reference_only_rows=len(missing),
                   eligible_rows=int(usable.sum()), invalid_by_column={k:int((~v).sum()) for k,v in valid.items()},
                   finite_optional_columns={k:int(np.isfinite(joined[k]).sum()) for k in POINT_FIELDS if k != 'photo_z'},
                   joined_reference_weight_fraction=float(ref_weight[ref_index].sum()/ref_weight_total),
                   eligible_reference_weight_fraction=float(ref_weight[ref_index[usable]].sum()/ref_weight_total),
                   point_columns_present_by_tract=availability,
                   joined_by_tract={str(tract):int((joined['tract']==tract).sum()) for tract in (9812,9813)},
                   effective_count_shape_som=float(normalized.sum()**2 / (normalized @ normalized)),
                   science_measurement=False, lens_selection_applied=False, photoz_calibration_applied=False)
    return joined, missing, summary


def write_tsv(path, table):
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(table)
        writer.writerows(zip(*table.values(), strict=True))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    plan = load_plan(DIRECTORY/'plan.json')
    destination = Path(plan['destination'])
    receipt = destination/CHECKPOINT
    state = json.loads(receipt.read_text())
    verify_state(plan, destination, state)
    entries = plan['files']
    if set(state['files']) != {e['filename'] for e in entries}:
        raise ValueError('requires all three verified FITS files')
    tables = [Table.read(destination/e['filename'], format='fits') for e in entries]
    joined, missing, summary = join_tables(tables[:2], tables[2])
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output/'cosmos_join.tsv', joined)
    write_tsv(output/'reference_only_ids.tsv', {'s17a_objid':missing})
    summary.update(plan_sha256=plan['plan_sha256'], receipt_sha256=sha(receipt),
                   input_sha256={e['filename']:state['files'][e['filename']]['sha256'] for e in entries},
                   output_sha256={p.name:sha(p) for p in output.iterdir()}, runner_sha256=sha(Path(__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(summary, stream, indent=2, allow_nan=False)
        stream.write('\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
