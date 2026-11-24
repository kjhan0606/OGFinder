#!/usr/bin/env python3
"""Inventory received FITS selection inputs without guessing field equivalence."""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from astropy.io import fits
from s16a_y1_source_contract import REQUIRED, select_row

INTAKE=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_ephor_ab_cosmos_intake_2026-11-19_01')
FILES=INTAKE.parent/'s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01/results'
EXPECTED='4b0c08e933c6e9d68209eb06b3fb773bc9d4cbc7c89c20ffc651a0c64fd4e84b'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary_path=INTAKE/'summary.json'
    if sha(summary_path)!=EXPECTED:
        raise ValueError('intake summary changed')
    intake=json.loads(summary_path.read_text())
    schemas={}
    for name,digest in intake['input_sha256'].items():
        path=FILES/name
        if sha(path)!=digest:
            raise ValueError('received FITS hash mismatch: '+name)
        with fits.open(path,memmap=True) as hdus:
            tables=[h for h in hdus if isinstance(h,fits.BinTableHDU)]
            if len(tables)!=1:
                raise ValueError('expected one FITS table')
            table=tables[0]
            schemas[name]=dict(rows=table.header['NAXIS2'],
                columns=[dict(name=c.name,format=c.format,unit=c.unit) for c in table.columns])
    # Inventory literal names only. Do not map weight_source to raw shape errors,
    # or infer selection flags from membership in a reference catalog.
    names={c['name'].casefold() for s in schemas.values() for c in s['columns']}
    fields=[dict(required_field=k,literal_column_present=k.casefold() in names,
                 status='requires_semantic_verification' if k.casefold() in names else 'absent') for k in REQUIRED]
    diagnostic=select_row({})
    output.mkdir(parents=True)
    table=output/'required_selection_fields.tsv'
    with table.open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(fields[0]),delimiter='\t')
        writer.writeheader();writer.writerows(fields)
    schema=output/'fits_schema.json'
    schema.write_text(json.dumps(schemas,indent=2)+'\n')
    report=dict(science_measurement=False,calibration_applied=False,
        intake_summary_sha256=EXPECTED,input_sha256=intake['input_sha256'],
        n_required_fields=len(REQUIRED),n_literal_names_present=sum(r['literal_column_present'] for r in fields),
        empty_input_diagnostic=diagnostic['table4_diagnostic'],
        n_unknown_criteria=sum(v=='unknown' for v in diagnostic['criteria'].values()),
        photoz_inputs_present=all(k in names for k in ('photoz_best','photoz_risk_best')),
        reference_inputs_present=all(k in names for k in ('cosmos_photoz','weight_source','som_weight','s17a_objid')),
        limitations=['Column presence is necessary but does not prove units, release compatibility or selection equivalence.',
            'No photometry or shape quality selection can be replayed from the received files alone.',
            'The all-unknown diagnostic describes missing input capability, not failure of every reference galaxy.',
            'Accepted pixel and Sirius masks need separate geometry evidence.'],
        runner_sha256=sha(Path(__file__)),output_sha256={p.name:sha(p) for p in (table,schema)})
    (output/'summary.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
