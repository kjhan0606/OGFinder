#!/usr/bin/env python3
"""Leave one COSMOS object out of each redshift-bin reweighted response."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, sha, write_table
from ogfkit.wlcluster import Cosmology
from ogfkit.wl_photoz import point_redshift_response

BASE = INTAKE.parent
SUPPORT = BASE/'s16a_cosmos_bin_support_2026-11-21_01'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    intake_file = INTAKE/'summary.json'
    support_file = SUPPORT/'summary.json'
    intake, support = (json.loads(p.read_text()) for p in (intake_file, support_file))
    for directory, payload in ((INTAKE, intake), (SUPPORT, support)):
        for filename, digest in payload['output_sha256'].items():
            if sha(directory/filename) != digest:
                raise ValueError('Input output digest mismatch')
    if support['intake_sha256'] != sha(intake_file):
        raise ValueError('Intake lineage differs')
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        data = [r for r in csv.DictReader(stream, delimiter='\t') if r['eligible_for_selection_review'] == 'True']
    if len(data) != intake['eligible_rows']:
        raise ValueError('Unexpected COSMOS eligible count')
    ids = np.asarray([r['s17a_objid'] for r in data])
    photo, reference, shape, som = np.asarray([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in data]).T
    with (SUPPORT/'photoz_bin_support.tsv').open() as stream:
        bins = list(csv.DictReader(stream, delimiter='\t'))
    with (SUPPORT/'lens_photoz_reweight.tsv').open() as stream:
        lenses = list(csv.DictReader(stream, delimiter='\t'))
    edges = np.asarray(support['edges'])
    cosmology = Cosmology()
    results = []
    for lens in lenses:
        zl = float(lens['lens_z'])
        keep = photo > zl + 0.2
        p, ref, sh, so, ident = photo[keep], reference[keep], shape[keep], som[keep], ids[keep]
        chi_lens = float(cosmology.comoving(zl))
        cp, cr = cosmology.comoving(p), cosmology.comoving(ref)
        est = 1-chi_lens/cp
        true = np.zeros_like(ref)
        behind = ref > zl
        true[behind] = 1-chi_lens/cr[behind]
        weight = sh*so*est**2
        ratio = true/est
        index = np.searchsorted(edges, p, side='right')-1
        count = len(edges)-1
        mass = np.bincount(index, weights=weight, minlength=count)
        numerator = np.bincount(index, weights=weight*ratio, minlength=count)
        observed = np.asarray([float(r['observed_weight_fraction']) for r in bins if r['lens_id']==lens['lens_id']])
        if len(observed) != count or np.any((observed > 0) & (mass <= 0)):
            raise ValueError('Incomplete bin support')
        averages = np.divide(numerator, mass, out=np.zeros_like(mass), where=mass>0)
        response = float(observed@averages)
        if not np.isclose(1/response, float(lens['reweighted']), rtol=1e-12):
            raise ValueError('Saved multiplier not reproduced')
        remaining = mass[index]-weight
        evaluable = (remaining > mass[index]*1e-12) & (observed[index]>0) & (weight>0)
        new_response = np.full(len(weight), np.nan)
        k = np.flatnonzero(evaluable)
        new_response[k] = response + observed[index[k]] * (
            (numerator[index[k]] - weight[k]*ratio[k]) / remaining[k] - averages[index[k]])
        change = (1/new_response[k])/(1/response)-1
        worst = k[np.argmax(np.abs(change))]
        # Recompute the most influential deletion using the production estimator.
        mask = np.arange(len(p)) != worst
        adjusted_mass = np.bincount(index[mask], weights=weight[mask], minlength=count)
        factor = np.divide(observed, adjusted_mass, out=np.zeros_like(observed), where=adjusted_mass>0)
        brute = point_redshift_response(zl,p[mask],ref[mask],sh[mask],so[mask]*factor[index[mask]],cosmology)
        if not np.isclose(brute['delta_sigma_multiplier'], 1/new_response[worst], rtol=1e-11):
            raise ValueError('Deletion formula differs from direct recalculation')
        results.append(dict(lens_id=lens['lens_id'], lens_z=zl, n_evaluable=len(k),
                            n_unsupported_deletions=int(np.sum((remaining <= mass[index]*1e-12) & (observed[index]>0) & (weight>0))),
                            baseline_multiplier=1/response, max_abs_fractional_shift=float(np.max(np.abs(change))),
                            signed_fractional_shift=float(1/new_response[worst]*response-1),
                            influential_object_id=ident[worst], photo_z=p[worst], reference_z=ref[worst],
                            bin_z_min=edges[index[worst]], bin_z_max=edges[index[worst]+1]))
    output.mkdir(parents=True)
    table=output/'lens_max_object_influence.tsv'
    write_table(table,results)
    summary=dict(science_measurement=False, calibration_applied=False,
                 input_intake_sha256=sha(intake_file), input_support_sha256=sha(support_file),
                 runner_sha256=sha(Path(__file__)), output_sha256={table.name:sha(table)},
                 max_abs_fractional_shift=max(r['max_abs_fractional_shift'] for r in results),
                 n_lenses=len(results), n_unsupported_deletions=sum(r['n_unsupported_deletions'] for r in results),
                 limitation='Single-object influence conditional on the observed histogram and 0.1 bins. '
                            'Not an uncertainty estimate; correlated objects, reference errors and cosmic variance remain.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
