#!/usr/bin/env python3
"""Compare received S17A COSMOS companion columns with the S16A Table-4 inputs; no selection.

Each required alias is matched by its literal name and by the S16A database column
it was queried from. Every match still needs semantic verification. The reference
table's whichfile/idx_in_file pointers are checked against the companion rows.
"""
import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import numpy as np
from astropy.io import fits

from prepare_s23b_s16a_shape_qa_pilot import COLUMN_EXPRESSIONS
from receive_s16a_cosmos_companion import DIRECTORY, load_plan
from receive_s16a_ephor_ab_cosmos import CHECKPOINT, verify_state
from s16a_y1_source_contract import REQUIRED

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
SUPPLEMENT_SQL = DATA/'wl_s23b_s16a_quality_supplement_plan_2026-11-18_01/s16a_quality_supplement_000.sql'
INTAKE = DATA/'s16a_ephor_ab_cosmos_intake_2026-11-19_01/summary.json'
INTAKE_SHA = '4b0c08e933c6e9d68209eb06b3fb773bc9d4cbc7c89c20ffc651a0c64fd4e84b'
REFERENCE = DATA/'s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01/results/Afterburner_reweighted_COSMOS_photoz_FDFC.fits'
TABLES = {'m': 's16a_wide.meas', 'a': 's16a_wide.meas2', 'f': 's16a_wide.forced'}


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 24), b''):
            digest.update(block)
    return digest.hexdigest()


def source_columns():
    """Alias -> S16A database column, read from the executed supplement SQL and shape query map."""
    mapping = {}
    for line in SUPPLEMENT_SQL.read_text().splitlines():
        match = re.fullmatch(r'\s*([a-z])\.([a-z0-9_]+) AS ([a-z0-9_]+),?', line)
        if match and match.group(1) in TABLES:
            mapping[match.group(3)] = dict(table=TABLES[match.group(1)], column=match.group(2))
    for alias in ('e1', 'e2'):
        prefix, column = COLUMN_EXPRESSIONS[alias].split('.')
        mapping[alias] = dict(table='s16a_wide.meas2', column=column)
    missing = [k for k in REQUIRED if k not in mapping]
    if missing:
        raise ValueError('no S16A source column for: '+', '.join(missing))
    return mapping


def companion_tables(plan, destination):
    state = json.loads((destination/CHECKPOINT).read_text())
    verify_state(plan, destination, state)
    if set(state['files']) != {e['filename'] for e in plan['files']}:
        raise ValueError('companion receipt is incomplete')
    tables = {}
    for name in sorted(state['files']):
        with fits.open(destination/name, memmap=True) as hdus:
            found = [h for h in hdus if isinstance(h, fits.BinTableHDU) and h.header.get('NAXIS2', 0)]
            if len(found) != 1:
                raise ValueError(f'expected one nonempty table in {name}')
            table = found[0]
            tables[name] = dict(hdu=table.name, rows=int(table.header['NAXIS2']),
                                columns=[dict(name=c.name, format=c.format, unit=c.unit) for c in table.columns])
    return state, tables


def match_fields(mapping, tables):
    columns = {c['name'].casefold(): (name, c['name']) for name, t in tables.items() for c in t['columns']}
    rows = []
    for alias in REQUIRED:
        source = mapping[alias]['column'].casefold()
        hit = columns.get(source) or columns.get(alias.casefold())
        candidates = sorted({orig for key, (_, orig) in columns.items()
                             if key != source and (source in key or (len(alias) > 3 and alias.casefold() in key))})
        status = ('s16a_column_name_present' if source in columns else
                  'alias_name_present' if alias.casefold() in columns else
                  'substring_candidates_only' if candidates else 'absent')
        rows.append(dict(required_field=alias, s16a_table=mapping[alias]['table'], s16a_column=mapping[alias]['column'],
                         status=status, matched_column=hit[1] if hit else '', matched_file=hit[0] if hit else '',
                         substring_candidates=';'.join(candidates[:12]),
                         note='name match only; semantics, units and release definition unverified' if status != 'absent' else ''))
    return rows


def pointer_check(tables, destination):
    if sha(INTAKE) != INTAKE_SHA:
        raise ValueError('reference intake summary changed')
    intake = json.loads(INTAKE.read_text())
    if sha(REFERENCE) != intake['input_sha256'][REFERENCE.name]:
        raise ValueError('reference FITS changed')
    with fits.open(REFERENCE, memmap=True) as hdus:
        ref = [h for h in hdus if isinstance(h, fits.BinTableHDU)][0].data
        names = {n.casefold(): n for n in ref.columns.names}
        which = np.char.strip(np.asarray(ref[names['whichfile']]).astype(str))
        index = np.asarray(ref[names['idx_in_file']]).astype(np.int64)
        objid = np.asarray(ref[names['s17a_objid']]).astype(np.int64)
    result = dict(reference_rows=int(len(which)), whichfile_values={str(k): int(v) for k, v in
                  zip(*np.unique(which, return_counts=True))}, files={})
    for name, table in tables.items():
        stem = name.casefold().removesuffix('.fits')
        chosen = np.array([stem in w.casefold() or w.casefold() in stem for w in which]) if len(which) else np.array([], bool)
        idx = index[chosen]
        in_range = (idx >= 0) & (idx < table['rows'])
        entry = dict(reference_rows_pointing_here=int(chosen.sum()), index_in_range=int(in_range.sum()),
                     index_out_of_range=int((~in_range).sum()), id_column=None, id_equal=None, id_differs=None)
        id_names = [c['name'] for c in table['columns'] if c['name'].casefold() in ('object_id', 'id', 's17a_objid', 'objid')]
        if id_names and in_range.any():
            with fits.open(destination/name, memmap=True) as hdus:
                data = [h for h in hdus if isinstance(h, fits.BinTableHDU) and h.header.get('NAXIS2', 0)][0].data
                ids = np.asarray(data[id_names[0]][idx[in_range]]).astype(np.int64)
            equal = ids == objid[chosen][in_range]
            entry.update(id_column=id_names[0], id_equal=int(equal.sum()), id_differs=int((~equal).sum()))
        result['files'][name] = entry
    covered = sum(f['reference_rows_pointing_here'] for f in result['files'].values())
    result['reference_rows_not_pointing_to_received_files'] = result['reference_rows']-covered
    return result


def run(output):
    if output.exists():
        raise FileExistsError(output)
    plan = load_plan(DIRECTORY/'plan.json')
    destination = Path(plan['destination'])
    state, tables = companion_tables(plan, destination)
    mapping = source_columns()
    fields = match_fields(mapping, tables)
    pointers = pointer_check(tables, destination)
    output.mkdir(parents=True)
    field_table = output/'required_field_matches.tsv'
    with field_table.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(fields[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(fields)
    schema = output/'companion_schema.json'
    schema.write_text(json.dumps(tables, indent=2)+'\n')
    counts = {s: sum(r['status'] == s for r in fields) for s in
              ('s16a_column_name_present', 'alias_name_present', 'substring_candidates_only', 'absent')}
    report = dict(science_measurement=False, calibration_applied=False, selection_replayed=False,
        plan_sha256=plan['plan_sha256'], receipt_sha256=sha(destination/CHECKPOINT),
        companion_sha256={name: item['sha256'] for name, item in state['files'].items()},
        supplement_sql_sha256=sha(SUPPLEMENT_SQL), reference_sha256=sha(REFERENCE),
        n_required_fields=len(fields), match_counts=counts, pointer_check=pointers,
        limitations=['A column-name match does not establish the same measurement, units, flag semantics or release algorithm (S17A vs S16A).',
                     'No Table-4 cut, accepted-pixel or bright-star mask is evaluated here.',
                     'Pointer checks verify row addressing only, not photometric agreement with S16A.'],
        runner_sha256=sha(Path(__file__)), output_sha256={p.name: sha(p) for p in (field_table, schema)})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(match_counts=counts, pointer_check=pointers), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
