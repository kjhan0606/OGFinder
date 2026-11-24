#!/usr/bin/env python3
"""Save proper radii and source-frame rotations for the verified S16A pairs."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from ogfkit import wl_y1
from audit_s16a_photoz_radial_support import FlatLambdaCDM, OBJECT_QA, PLAN, load_inputs, radial_bins, sha, separations_radians
from s19a_shape_stack import pair_geometry


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, prior, provenance = load_inputs(PLAN, OBJECT_QA)
    ledger = {row['center_id']:row for row in centers}
    center_rows = [ledger[p['center_id']] for p in pairs]
    coords = [np.array([float(r[k]) for r in rows]) for rows,k in
              [(center_rows,'ra_deg'), (center_rows,'dec_deg'), (pairs,'source_ra'), (pairs,'source_dec')]]
    lens_z = np.array([float(row['z_lens']) for row in center_rows])
    angle = wl_y1.source_frame_angle(*coords)
    lens = SkyCoord(coords[0]*u.deg, coords[1]*u.deg)
    source = SkyCoord(coords[2]*u.deg, coords[3]*u.deg)
    reference_angle = np.pi/2-source.position_angle(lens).rad

    def spin_difference(left, right):
        return np.abs(np.angle(np.exp(2j*(left-right))))/2

    residual = spin_difference(angle, reference_angle)
    if residual.max() > 1e-9:
        raise ValueError('source-frame angle differs from Astropy')
    sep = separations_radians(coords[2], coords[3], coords[0], coords[1])
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    distance = cosmo.angular_diameter_distance(lens_z).value
    radius = sep*distance
    reference_radius = source.separation(lens).rad*distance
    if not np.allclose(radius, reference_radius, rtol=1e-9, atol=1e-10):
        raise ValueError('proper radius differs from Astropy')
    legacy_radius, legacy_angle = pair_geometry(dict(zip(('lens_ra_deg','lens_dec_deg','source_ra_deg','source_dec_deg','lens_z'),
                                                         [*coords,lens_z], strict=True)),
                                               SimpleNamespace(angular_diameter=lambda z: cosmo.angular_diameter_distance(z).value))
    if not np.allclose(radius, legacy_radius, rtol=1e-9, atol=1e-10):
        raise ValueError('previous and new radius differ')
    bins = radial_bins(radius)
    counts = Counter((row['center_id'], int(k)) for row,k in zip(pairs,bins,strict=True) if row['object_id'] in selected)
    if any(counts[row['center_id'],int(row['radial_bin'])] != int(row['selected_count']) for row in prior):
        raise ValueError('selected annular counts changed')
    output.mkdir(parents=True, exist_ok=False)
    table = output/'pair_geometry.tsv'
    with table.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(('center_id','object_id','radius_proper_mpc','radial_bin','source_angle_rad','table4','finite_photoz'))
        for row,r,k,phi in zip(pairs,radius,bins,angle,strict=True):
            obj = row['object_id']
            writer.writerow((row['center_id'],obj,r,k,phi,obj in selected,obj in available))
    difference = np.rad2deg(spin_difference(angle,legacy_angle))*3600
    summary = dict(**provenance, pairs=len(pairs), checked_center_annuli=len(prior),
                   max_astropy_angle_error_rad=float(residual.max()), max_astropy_radius_error_mpc=float(np.max(np.abs(radius-reference_radius))),
                   legacy_lens_frame_difference_arcsec=dict(zip(('median','q95','max'),map(float,np.quantile(difference,[.5,.95,1])),strict=True)),
                   cosmology=dict(h0=70,omega_m=.3,tcmb0=0),
                   science_measurement=False, observed_shapes_rotated=False, shear_computed=False,
                   output_sha256={table.name:sha(table)}, runner_sha256=sha(Path(__file__)), kernel_sha256=sha(Path(wl_y1.__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(summary,indent=2))


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    run(p.parse_args().output.resolve())
