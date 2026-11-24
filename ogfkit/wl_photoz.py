"""Photo-z response for a weak-shear, point-redshift DeltaSigma estimator.

See Nakajima et al. (2012), section 7.1.  Supply the calibration sample after
the same estimated-redshift and shape selection used in the actual pairs.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from ogfkit.wlcluster import C_KMS, G_MPC, Cosmology


def point_redshift_response(lens_z, photo_z, reference_z, shape_weight, som_weight, cosmology):
    """Return 1+b_z and its reciprocal for one fixed lens redshift.

    Photo-z foregrounds have zero estimator weight. Reference foregrounds
    with background photo-z stay in the denominator and dilute the response.
    The reference redshift is assumed to supply the calibration truth proxy.
    """
    z_est, z_ref, shape, som = (np.asarray(v, dtype=float) for v in
                                 (photo_z, reference_z, shape_weight, som_weight))
    arrays = (z_est, z_ref, shape, som)
    if z_est.ndim != 1 or any(a.shape != z_est.shape for a in arrays):
        raise ValueError('redshift and weight vectors must have the same one-dimensional shape')
    if not np.isfinite(lens_z) or lens_z <= 0:
        raise ValueError('expected a finite positive lens redshift')
    if any(not np.all(np.isfinite(a) & (a >= 0)) for a in arrays):
        raise ValueError('redshifts and weights must be finite and nonnegative')

    lens_chi = float(cosmology.comoving(lens_z))
    if not np.isfinite(lens_chi) or lens_chi <= 0:
        raise ValueError('lens redshift outside the distance domain')
    prefactor = 4*np.pi*G_MPC / C_KMS**2 * lens_chi / (1+lens_z)

    def inverse_sigma(source_z):
        source_chi = cosmology.comoving(source_z)
        if not np.all(np.isfinite(source_chi)):
            raise ValueError('source redshift outside the distance domain')
        beta = np.zeros_like(source_chi)
        behind = source_z > lens_z
        beta[behind] = prefactor * (1-lens_chi/source_chi[behind])
        return beta

    estimated = inverse_sigma(z_est)
    reference = inverse_sigma(z_ref)
    weights = shape * som * estimated**2
    weight_sum = float(weights.sum())
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        raise ValueError('no positive calibration weight')
    used = weights > 0
    normalized = weights[used]/weight_sum
    response = float(normalized @ (reference[used]/estimated[used]))
    if not np.isfinite(response) or response <= 0:
        raise ValueError('reference response is zero or invalid')
    return dict(response=response, bias=response-1, delta_sigma_multiplier=1/response,
                used_pairs=int(used.sum()), weight_sum=weight_sum,
                neff_calibration=float(1/(normalized @ normalized)),
                reference_foreground_weight_fraction=float(normalized[reference[used] == 0].sum()))


def stack_photoz_response(responses, observed_pair_weight_sums):
    """Combine lens responses using actual measurement denominators in one bin."""
    response = np.asarray(responses, dtype=float)
    weights = np.asarray(observed_pair_weight_sums, dtype=float)
    if (response.ndim != 1 or weights.shape != response.shape
            or not np.all(np.isfinite(weights) & (weights >= 0))):
        raise ValueError('expected matching vectors and finite nonnegative stack weights')
    used = weights > 0
    total = weights[used].sum()
    if not np.isfinite(total) or total <= 0:
        raise ValueError('no positive observed pair weight')
    if not np.all(np.isfinite(response[used]) & (response[used] > 0)):
        raise ValueError('every contributing lens needs a positive calibrated response')
    combined = float((weights[used]/total) @ response[used])
    return dict(response=combined, bias=combined-1, delta_sigma_multiplier=1/combined)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('table', type=Path, help='Selected per-lens TSV with photo_z, reference_z, shape_weight, som_weight')
    parser.add_argument('--lens-z', type=float, required=True)
    parser.add_argument('--h0', type=float, required=True)
    parser.add_argument('--omega-m', type=float, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not (np.isfinite(args.h0) and args.h0 > 0 and 0 < args.omega_m <= 1):
        parser.error('expected positive H0 and 0 < Omega_m <= 1')
    if args.output.exists():
        raise FileExistsError(args.output)
    names = ('photo_z', 'reference_z', 'shape_weight', 'som_weight')
    raw = args.table.read_bytes()
    reader = csv.DictReader(raw.decode('utf-8-sig').splitlines(), delimiter='\t')
    if not set(names) <= set(reader.fieldnames or []):
        parser.error('required calibration columns are missing')
    data = list(reader)
    columns = [[float(row[name]) for row in data] for name in names]
    result = point_redshift_response(args.lens_z, *columns, Cosmology(H0=args.h0, Om=args.omega_m))
    report = dict(**result, science_measurement=False, lens_z=args.lens_z, h0=args.h0, omega_m=args.omega_m,
                  surface_density_convention='proper', input_path=str(args.table.resolve()),
                  input_sha256=hashlib.sha256(raw).hexdigest(), input_rows=len(data))
    with args.output.open('x') as stream:
        json.dump(report, stream, indent=2, allow_nan=False)
        stream.write('\n')


if __name__ == '__main__':
    main()
