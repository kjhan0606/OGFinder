#!/usr/bin/env python3
"""Relate absent Ephor_AB rows to null afterburner fluxes in the received COSMOS tables.

S17A Wide-COSMOS only. It describes how Ephor_AB rows exist for these tracts and
does not establish the cause of missing S16A Wide photo-z in the lens fields.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from receive_s16a_cosmos_companion import DIRECTORY, load_plan
from receive_s16a_ephor_ab_cosmos import CHECKPOINT, verify_state

DATA = Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
INTAKE = DATA/'s16a_ephor_ab_cosmos_intake_2026-11-19_01/summary.json'
INTAKE_SHA = '4b0c08e933c6e9d68209eb06b3fb773bc9d4cbc7c89c20ffc651a0c64fd4e84b'
EPHOR = DATA/'s16a_y1_ephor_ab_cosmos_plan_2026-11-18_01/results'
BANDS = 'grizy'


def sha(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as stream:
        for block in iter(lambda: stream.read(1 << 24), b''):
            digest.update(block)
    return digest.hexdigest()


def table(path):
    with fits.open(path, memmap=True) as hdus:
        found = [h for h in hdus if isinstance(h, fits.BinTableHDU) and h.header.get('NAXIS2', 0)]
        if len(found) != 1:
            raise ValueError(f'expected one nonempty table: {path.name}')
        return found[0].data.copy()


def unavailable(data, column):
    return np.asarray(data[column+'_isnull'], bool) | ~np.isfinite(np.asarray(data[column], float))


def run(output):
    if output.exists():
        raise FileExistsError(output)
    plan = load_plan(DIRECTORY/'plan.json')
    destination = Path(plan['destination'])
    state = json.loads((destination/CHECKPOINT).read_text())
    verify_state(plan, destination, state)
    if sha(INTAKE) != INTAKE_SHA:
        raise ValueError('Ephor intake summary changed')
    intake = json.loads(INTAKE.read_text())
    tracts, totals = {}, dict(companion=0, absent=0, absent_afterburner_unavailable=0,
                              present=0, present_afterburner_unavailable=0)
    for tract in ('9812', '9813'):
        cname, ename = f'target_wide_s17a_{tract}.fits', f'pdf-s17a_wide-{tract}.cat.fits'
        if sha(EPHOR/ename) != intake['input_sha256'][ename]:
            raise ValueError('Ephor table changed: '+ename)
        companion, ephor = table(destination/cname), table(EPHOR/ename)
        cid, eid = np.asarray(companion['object_id'], np.int64), np.asarray(ephor['ID'], np.int64)
        if len(np.unique(cid)) != len(cid) or len(np.unique(eid)) != len(eid):
            raise ValueError('duplicate object IDs')
        present = np.isin(cid, eid)
        after = np.zeros(len(cid), bool)
        per_band = {}
        for band in BANDS:
            bad = unavailable(companion, band+'parent_flux_convolved_2_1')
            per_band[band] = dict(absent=int(bad[~present].sum()), present=int(bad[present].sum()))
            after |= bad
        cmodel = np.zeros(len(cid), bool)
        for band in BANDS:
            cmodel |= unavailable(companion, band+'cmodel_flux')
        nonfinite_best = ~np.isfinite(np.asarray(ephor['PHOTOZ_BEST'], float)) | ~np.isfinite(np.asarray(ephor['PHOTOZ_RISK_BEST'], float))
        tracts[tract] = dict(companion_rows=int(len(cid)), ephor_rows=int(len(eid)),
            ephor_ids_absent_from_companion=int((~np.isin(eid, cid)).sum()),
            companion_absent_from_ephor=int((~present).sum()),
            absent_with_afterburner_unavailable=int(after[~present].sum()),
            present_with_afterburner_unavailable=int(after[present].sum()),
            absent_with_cmodel_unavailable=int(cmodel[~present].sum()),
            present_with_cmodel_unavailable=int(cmodel[present].sum()),
            ephor_rows_nonfinite_best_or_risk=int(nonfinite_best.sum()),
            afterburner_unavailable_by_band=per_band)
        totals['companion'] += len(cid)
        totals['absent'] += int((~present).sum())
        totals['absent_afterburner_unavailable'] += int(after[~present].sum())
        totals['present'] += int(present.sum())
        totals['present_afterburner_unavailable'] += int(after[present].sum())
    exact = (totals['absent'] == totals['absent_afterburner_unavailable'] and totals['present_afterburner_unavailable'] == 0)
    report = dict(science_measurement=False, calibration_applied=False,
        afterburner_unavailable_iff_absent_ephor_row=exact, totals=totals, tracts=tracts,
        definition='Afterburner unavailable = any of g,r,i,z,y parent_flux_convolved_2_1 flagged _isnull or non-finite; '
                   'absent = companion object_id with no row in the received Ephor_AB point table.',
        limitations=['S17A Wide-COSMOS tracts 9812/9813 only; transfer to S16A Wide lens fields is an inference, not a measurement.',
                     'Received Ephor tables only contain rows; a row with non-finite best/risk would be a different missingness mode.',
                     'Does not identify why afterburner fluxes are unavailable (e.g. coadd coverage) or test sky geometry.'],
        plan_sha256=plan['plan_sha256'], companion_sha256={k: v['sha256'] for k, v in state['files'].items()},
        ephor_sha256={k: v for k, v in intake['input_sha256'].items() if k.startswith('pdf-')},
        runner_sha256=sha(Path(__file__)))
    output.mkdir(parents=True)
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(exact=exact, totals=totals), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
