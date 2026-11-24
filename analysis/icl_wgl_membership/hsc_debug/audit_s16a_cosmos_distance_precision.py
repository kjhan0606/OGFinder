#!/usr/bin/env python3
"""Compare the saved photo-z responses with direct FLRW distances."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

from audit_s16a_cosmos_lens_response import INTAKE, sha, write_table
from ogfkit.wl_photoz import point_redshift_response


class DirectCosmology:
    def __init__(self):
        self.model = FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0)

    def comoving(self, z):
        return np.asarray(self.model.comoving_distance(z).value)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    intake_path = INTAKE / 'summary.json'
    intake = json.loads(intake_path.read_text())
    table = INTAKE / 'cosmos_join.tsv'
    if sha(table) != intake['output_sha256'][table.name]:
        raise ValueError('Changed COSMOS table')
    with table.open() as stream:
        rows = [r for r in csv.DictReader(stream, delimiter='\t')
                if r['eligible_for_selection_review'] == 'True']
    if len(rows) != intake['eligible_rows']:
        raise ValueError('Eligibility count differs')
    photo, reference, shape, som, risk = np.asarray([
        [float(r[k]) for k in ('photo_z', 'reference_z', 'shape_weight', 'som_weight', 'photo_z_risk')]
        for r in rows]).T
    prior_dir = INTAKE.parent/'s16a_cosmos_lens_response_2026-11-19_01'
    prior = json.loads((prior_dir/'summary.json').read_text())
    response_table = prior_dir/'lens_response.tsv'
    if sha(response_table) != prior['output_sha256'][response_table.name]:
        raise ValueError('Changed response table')
    if prior['cosmos_intake_sha256'] != sha(intake_path):
        raise ValueError('Input lineage differs')
    with response_table.open() as stream:
        cases = [r for r in csv.DictReader(stream, delimiter='\t') if r['cosmos_tract'] == 'combined']
    if len(cases) != 33:
        raise ValueError('Expected 11 lenses and 3 gaps')
    model = DirectCosmology()
    result = []
    for case in cases:
        zl, gap = float(case['lens_z']), float(case['gap'])
        keep = np.isfinite(risk) & (photo > zl + gap)
        direct = point_redshift_response(zl, photo[keep], reference[keep], shape[keep], som[keep], model)
        old = float(case['delta_sigma_multiplier'])
        result.append(dict(lens_id=case['lens_id'], lens_z=zl, gap=gap,
                           saved_multiplier=old, direct_multiplier=direct['delta_sigma_multiplier'],
                           fractional_difference=direct['delta_sigma_multiplier']/old-1,
                           n_calibration=direct['used_pairs']))
    output.mkdir(parents=True)
    target = output/'distance_precision.tsv'
    write_table(target, result)
    report = dict(science_measurement=False, calibration_applied=False,
                  intake_summary_sha256=sha(intake_path), response_summary_sha256=sha(prior_dir/'summary.json'),
                  runner_sha256=sha(Path(__file__)), output_sha256={target.name:sha(target)},
                  n_cases=len(result), max_abs_fractional_difference=max(abs(r['fractional_difference']) for r in result),
                  definition='Compare saved 512-point trapezoid distance grid with Astropy FLRW distances '
                             'for the same calibration rows and estimator.',
                  limitation='Numerical distance precision only; does not address sample selection, '
                             'reference redshift error, cosmic variance, or covariance.')
    (output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    run(parser.parse_args().output.resolve())
