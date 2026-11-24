#!/usr/bin/env python3
"""Relate S16A Ephor_AB row absence to afterburner flux availability and coadd input counts.

Reads only the approved afterburner/photo-z receipt and the earlier verified supplement;
no network. Descriptive: no selection change, correction or science measurement.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path

import numpy as np

from audit_s16a_missingness_geometry import CELL_DEG, MIN_CELL
from audit_s16a_photoz_radial_support import PLAN as SUPPLEMENT, OBJECT_QA, load_inputs
from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_s16a_afterburner_photoz_rows import BANDS, OUTPUT as DIRECTORY
from submit_s16a_afterburner_photoz_rows import load_plan, verify_complete

AB_STATES = ('afterburner_row_absent', 'flux_null_any_band', 'flux_nonfinite_any_band', 'flux_finite_all_bands')
Z_STATES = ('photoz_row_absent', 'photoz_row_value_unavailable', 'photoz_row_finite')
FLAGS = ('parent_flux_convolved_2_1_flags', 'parent_flux_convolved_flag', 'parent_flux_convolved_2_deconv')


def number(value):
    value = value.strip()
    return None if value == '' else float(value)


def classify(row):
    if row['afterburner_object_id'].strip() == '':
        ab = 'afterburner_row_absent'
    else:
        fluxes = [number(row[b+'parent_flux_convolved_2_1']) for b in BANDS]
        ab = ('flux_null_any_band' if any(f is None for f in fluxes) else
              'flux_nonfinite_any_band' if not all(math.isfinite(f) for f in fluxes) else 'flux_finite_all_bands')
    if row['photoz_object_id'].strip() == '':
        z = 'photoz_row_absent'
    else:
        best, risk = number(row['photoz_best']), number(row['photoz_risk_best'])
        z = 'photoz_row_finite' if best is not None and risk is not None and math.isfinite(best) and math.isfinite(risk) \
            else 'photoz_row_value_unavailable'
    return ab, z


def band_state(row, band):
    if row['afterburner_object_id'].strip() == '':
        return 'row_absent'
    value = number(row[band+'parent_flux_convolved_2_1'])
    return 'null' if value is None else 'nonfinite' if not math.isfinite(value) else 'nonpositive' if value <= 0 else 'positive'


def contingency(rows, ids):
    table = Counter(classify(rows[i]) for i in ids)
    return [dict(afterburner=a, **{z: table[(a, z)] for z in Z_STATES}) for a in AB_STATES]


def quantiles(values):
    values = np.asarray([v for v in values if v is not None], float)
    if not len(values):
        return dict(n=0)
    q = np.quantile(values, [0.05, 0.25, 0.5, 0.75, 0.95])
    return dict(n=int(len(values)), q05=float(q[0]), q25=float(q[1]), median=float(q[2]), q75=float(q[3]), q95=float(q[4]))


def write_tsv(path, rows):
    with path.open('x', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader(); writer.writerows(rows)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    plan = load_plan(DIRECTORY)
    results = DIRECTORY/'results'
    state = json.loads((results/'status.json').read_text())
    receipt = verify_complete(plan, DIRECTORY, state)
    if state['status'] != 'complete' or sha(results/'download_summary.json') != state['summary_sha256']:
        raise ValueError('afterburner receipt is incomplete or changed')
    rows = {}
    for name in receipt['output_sha256']:
        with (results/name).open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                rows[row['object_id']] = row
    centers, pairs, selected, available, _, provenance = load_inputs(SUPPLEMENT, OBJECT_QA)
    position = {}
    for p in pairs:
        position.setdefault(p['object_id'], (float(p['source_ra']), float(p['source_dec'])))
    if set(rows) != set(position):
        raise ValueError('receipt IDs differ from the pilot ledger')
    everyone = sorted(rows, key=int)
    chosen = [o for o in everyone if o in selected]
    missing = {o for o in everyone if o not in available}

    new_finite = {o for o in everyone if classify(rows[o])[1] == 'photoz_row_finite'}
    consistency = dict(earlier_finite=len(available), new_finite=len(new_finite),
                       finite_now_missing_before=len(new_finite-available), missing_now_finite_before=len(available-new_finite))

    absent = np.array([classify(rows[o])[1] == 'photoz_row_absent' for o in chosen])
    unavailable = np.array([classify(rows[o])[0] != 'flux_finite_all_bands' for o in chosen])
    iff = dict(selected_objects=len(chosen), photoz_row_absent=int(absent.sum()),
               afterburner_unavailable=int(unavailable.sum()),
               absent_and_unavailable=int((absent & unavailable).sum()),
               absent_but_afterburner_finite=int((absent & ~unavailable).sum()),
               present_but_afterburner_unavailable=int((~absent & unavailable).sum()),
               exact=bool(np.array_equal(absent, unavailable)))

    band_rows = []
    for group, ids in (('selected_photoz_missing', [o for o in chosen if o in missing]),
                       ('selected_photoz_finite', [o for o in chosen if o not in missing])):
        for band in BANDS:
            counts = Counter(band_state(rows[o], band) for o in ids)
            flags, unmeasured_deconv = Counter(), 0
            for o in ids:
                for suffix in FLAGS:
                    if rows[o][band+suffix].strip().lower() in ('true', 't'):
                        flags[suffix] += 1
                if band_state(rows[o], band) in ('null', 'nonfinite') and \
                        rows[o][band+'parent_flux_convolved_2_deconv'].strip().lower() in ('true', 't'):
                    unmeasured_deconv += 1
            band_rows.append(dict(group=group, band=band, objects=len(ids),
                                  **{s: counts[s] for s in ('row_absent', 'null', 'nonfinite', 'nonpositive', 'positive')},
                                  flags_2_1_true=flags['parent_flux_convolved_2_1_flags'],
                                  convolved_flag_true=flags['parent_flux_convolved_flag'],
                                  deconv_2_true=flags['parent_flux_convolved_2_deconv'],
                                  unavailable_with_deconv_2=unmeasured_deconv))

    def cause(row):
        found = []
        for band in BANDS:
            if band_state(row, band) in ('row_absent', 'null', 'nonfinite'):
                seeing = number(row[band+'parent_flux_convolved_seeing'])
                found.append(band+(':deconv' if row[band+'parent_flux_convolved_2_deconv'].strip().lower() in ('true', 't')
                                   else ':seeing_unavailable' if seeing is None or not math.isfinite(seeing) else ':other'))
        return '+'.join(found) or 'none'
    causes = {scope: dict(Counter(cause(rows[o]) for o in ids if classify(rows[o])[1] == 'photoz_row_absent'))
              for scope, ids in (('selected', chosen), ('all_pilot', everyone))}
    threshold = {}
    for band in BANDS:
        flagged, unflagged = [], []
        for o in everyone:
            seeing = number(rows[o][band+'parent_flux_convolved_seeing'])
            if seeing is not None and math.isfinite(seeing):
                (flagged if rows[o][band+'parent_flux_convolved_2_deconv'].strip().lower() in ('true', 't') else unflagged).append(seeing)
        threshold[band] = dict(deconv_true=len(flagged), min_seeing_deconv_true=min(flagged) if flagged else None,
                               max_seeing_deconv_false=max(unflagged) if unflagged else None,
                               separated=bool(flagged) and min(flagged) > max(unflagged))
    # Inference, not a schema statement: 6.5 px FWHM = 1.09 arcsec at 0.168 arcsec/px, as a Gaussian sigma in px.
    sigma_target = 6.5/(2*math.sqrt(2*math.log(2)))

    cells = defaultdict(list)
    for o in chosen:
        ra, dec = position[o]
        cells[(int(np.floor(ra/CELL_DEG)), int(np.floor(dec/CELL_DEG)))].append(o)
    mixed = [o for v in cells.values() if len(v) >= MIN_CELL and 0.1 <= np.mean([x in missing for x in v]) <= 0.9 for o in v]
    count_rows = []
    for scope, ids in (('all_selected', chosen), ('mixed_cells', mixed)):
        for quantity in ('countinputs', 'parent_flux_convolved_seeing'):
            for band in BANDS:
                for label, group in (('missing', [o for o in ids if o in missing]), ('finite', [o for o in ids if o not in missing])):
                    values = [number(rows[o][band+quantity]) for o in group]
                    count_rows.append(dict(scope=scope, quantity=quantity, band=band, photoz=label,
                                           null=sum(v is None for v in values),
                                           nonfinite=sum(v is not None and not math.isfinite(v) for v in values),
                                           **quantiles([v for v in values if v is not None and math.isfinite(v)])))

    output.mkdir(parents=True)
    tables = {'contingency_selected.tsv': contingency(rows, chosen), 'contingency_all_pilot.tsv': contingency(rows, everyone),
              'band_status.tsv': band_rows, 'band_inputs_seeing.tsv': count_rows}
    for name, table in tables.items():
        write_tsv(output/name, table)
    report = dict(science_measurement=False, calibration_applied=False, selection_changed=False,
        photoz_consistency_with_earlier_supplement=consistency, selected_row_absence_vs_afterburner=iff,
        absent_row_causes=causes, deconv_seeing_threshold_by_band=threshold,
        inferred_seeing_threshold=dict(value=sigma_target, interpretation='Gaussian sigma in pixels of a 6.5 px (1.09 arcsec) FWHM; '
                                       'inferred from the observed flag boundary, the schema gives no unit'),
        mixed_cell_definition=f'{CELL_DEG} deg RA/Dec cells with >= {MIN_CELL} selected objects and missing fraction 0.1-0.9',
        mixed_cell_objects=len(mixed),
        definitions=dict(afterburner_unavailable='afterburner row absent, or any grizy parent_flux_convolved_2_1 NULL or non-finite',
                         photoz_missing='no finite photoz_best and photoz_risk_best in the earlier supplement'),
        limitations=['countinputs is from s16a_wide.forced and counts coadd inputs at the object centre only.',
                     'parent_flux_convolved_seeing is the afterburner PSF size at the object; the 1.09-arcsec target and '
                     'the deconv flag follow the schema description, not a re-measurement.',
                     'Coincidence of row absence and afterburner unavailability does not identify why the afterburner failed.'],
        plan_sha256=plan['query_plan_sha256'], receipt_sha256=sha(results/'download_summary.json'),
        supplement_provenance=provenance, runner_sha256=sha(Path(__file__)),
        output_sha256={name: sha(output/name) for name in tables})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps(dict(consistency=consistency, iff=iff, causes=causes, threshold=threshold), indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
