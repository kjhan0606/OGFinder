#!/usr/bin/env python3
"""Locate missing S16A Ephor_AB best/risk inputs on the sky; no calibration.

Reads the verified quality-supplement receipt and the fixed 66-center pilot.
Reports per tract/patch missingness (decoded from object_id), nearest-neighbour
clustering, same-cell property contrasts and a sky map. Descriptive only.
"""
import argparse
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_s16a_photoz_radial_support import PLAN, OBJECT_QA, load_inputs, catalogue_reader, write_table

UPSTREAM = Path('/home/kjhan/BACKUP/Astrafek/s16a_selection_missingness_20261121_01/summary.json')
UPSTREAM_SHA = '755e1dbf6ff721da772d6765f96434efa318eb3f64d7d672141272f3a2cb4d9b'
CELL_DEG = 0.05
MIN_CELL = 20


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def number(row, key):
    try:
        value = float(row[key])
    except (TypeError, ValueError):
        return np.nan
    return value if np.isfinite(value) else np.nan


def quantiles(values):
    x = values[np.isfinite(values)]
    return [float(v) for v in np.percentile(x, [16, 50, 84])] if len(x) else None


def run(output):
    if output.exists():
        raise FileExistsError(output)
    if sha(UPSTREAM) != UPSTREAM_SHA:
        raise ValueError('upstream missingness summary changed')
    upstream = json.loads(UPSTREAM.read_text())
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    if provenance != upstream['provenance']:
        raise ValueError('pilot lineage differs from upstream missingness')
    receipt = json.loads((PLAN/'results'/'download_summary.json').read_text())
    props = {}
    for name, digest in receipt['output_sha256'].items():
        path = PLAN/'results'/name
        if sha(path) != digest:
            raise ValueError('supplement digest mismatch')
        with path.open(newline='', encoding='utf-8-sig') as stream:
            for row in catalogue_reader(stream):
                props[row['object_id']] = row
    lens_of = {c['center_id']: c['matched_lens_id'] for c in centers}
    objects = {}
    for p in pairs:
        oid = p['object_id']
        if oid not in selected:
            continue
        entry = objects.setdefault(oid, dict(ra=float(p['source_ra']), dec=float(p['source_dec']),
                                             weight=float(p['weight']), lenses=set()))
        if abs(entry['ra']-float(p['source_ra'])) > 1e-9 or abs(entry['dec']-float(p['source_dec'])) > 1e-9:
            raise ValueError('object position changes between centers')
        entry['lenses'].add(lens_of[p['center_id']])
    if set(objects) != selected:
        raise ValueError('selected ledger differs')
    ids = sorted(objects)
    missing = np.array([o not in available for o in ids])
    if int(missing.sum()) != upstream['source_object_counts']['missing']['n'] or len(ids) != upstream['selected_objects']:
        raise ValueError('counts differ from upstream missingness')

    # HSC coadd object_id: tract in bits >= 42, patch code in bits 29-41 (x = code>>8, y = (code>>3)&31).
    code = {o: (int(o) >> 42, (int(o) >> 29) & 0x1FFF) for o in ids}
    units = defaultdict(lambda: dict(n=0, missing=0, lenses=set()))
    for o, m in zip(ids, missing):
        u = units[code[o]]
        u['n'] += 1
        u['missing'] += int(m)
        u['lenses'] |= objects[o]['lenses']
    patch_rows = [dict(tract=t, patch_code=c, patch_x=c >> 8, patch_y=(c >> 3) & 31, selected_objects=u['n'],
                       missing_objects=u['missing'], missing_fraction=u['missing']/u['n'],
                       lens_groups=','.join(sorted(u['lenses'])))
                  for (t, c), u in sorted(units.items())]
    fractions = np.array([r['missing_fraction'] for r in patch_rows])

    by_lens = []
    for lens in sorted({lens_of[c['center_id']] for c in centers}):
        sub = [i for i, o in enumerate(ids) if lens in objects[o]['lenses']]
        m = missing[sub]
        ra = np.array([objects[ids[i]]['ra'] for i in sub])
        dec = np.array([objects[ids[i]]['dec'] for i in sub])
        row = dict(lens_id=lens, unique_objects=len(sub), missing_objects=int(m.sum()),
                   missing_fraction=float(m.mean()) if len(sub) else None,
                   nn_missing_given_missing=None, nn_missing_given_finite=None)
        if m.sum() and len(sub) > 2:
            xy = np.c_[ra*np.cos(np.radians(dec.mean())), dec]
            nn = cKDTree(xy).query(xy, k=2)[1][:, 1]
            row['nn_missing_given_missing'] = float(m[nn][m].mean())
            row['nn_missing_given_finite'] = float(m[nn][~m].mean())
        by_lens.append(row)

    cells = defaultdict(list)
    for i, o in enumerate(ids):
        cells[(int(np.floor(objects[o]['ra']/CELL_DEG)), int(np.floor(objects[o]['dec']/CELL_DEG)))].append(i)
    mixed = [i for v in cells.values() if len(v) >= MIN_CELL and 0.1 <= missing[v].mean() <= 0.9 for i in v]
    mixed = np.array(mixed, dtype=int)
    mm = missing[mixed]
    features = {'i_mag_extinction_corrected': lambda r: number(r, 'i_mag')-number(r, 'a_i'),
                'extendedness': lambda r: number(r, 'extendedness'),
                'blendedness': lambda r: number(r, 'blendedness'),
                'resolution': lambda r: number(r, 'resolution')}
    for band in 'grizy':
        features[band+'_cmodel_snr'] = (lambda b: lambda r: number(r, b+'_flux')/number(r, b+'_flux_err'))(band)
    contrast = {}
    for name, fn in features.items():
        values = np.array([fn(props[ids[i]]) for i in mixed])
        contrast[name] = dict(finite=quantiles(values[~mm]), missing=quantiles(values[mm]))
    weights = np.array([objects[ids[i]]['weight'] for i in mixed])
    contrast['shape_weight'] = dict(finite=quantiles(weights[~mm]), missing=quantiles(weights[mm]))
    cmodel_nan = {b: dict(finite=float(np.mean([np.isnan(number(props[o], b+'_flux')) for o, m in zip(ids, missing) if not m])),
                          missing=float(np.mean([np.isnan(number(props[o], b+'_flux')) for o, m in zip(ids, missing) if m])))
                  for b in 'grizy'}

    output.mkdir(parents=True)
    patch_table = output/'patch_missingness.tsv'
    lens_table = output/'lens_spatial_clustering.tsv'
    write_table(patch_table, patch_rows)
    write_table(lens_table, by_lens)
    figure = output/'missing_sky.png'
    plot(figure, objects, ids, missing, centers, [r['lens_id'] for r in sorted(by_lens, key=lambda r: -r['missing_objects'])[:2]])
    report = dict(science_measurement=False, calibration_applied=False,
        selected_objects=len(ids), missing_objects=int(missing.sum()),
        tract_patch_units=len(patch_rows), units_without_missing=int((fractions == 0).sum()),
        units_fully_missing=int((fractions == 1).sum()),
        units_partially_missing=int(((fractions > 0) & (fractions < 1)).sum()),
        largest_lens_share=dict(lens_id=max(by_lens, key=lambda r: r['missing_objects'])['lens_id'],
                                fraction_of_missing=max(r['missing_objects'] for r in by_lens)/int(missing.sum())),
        mixed_cell_definition=f'{CELL_DEG} deg RA/Dec cells with >= {MIN_CELL} objects and missing fraction 0.1-0.9',
        mixed_cell_objects=int(len(mixed)), mixed_cell_missing=int(mm.sum()),
        same_cell_property_quantiles_p16_p50_p84=contrast, cmodel_flux_nan_fraction=cmodel_nan,
        interpretation='Missing best/risk values are spatially coherent (axis-aligned stripes and blocks inside patches), '
                       'not whole missing patches; same-cell object properties are indistinguishable.',
        limitations=['Tract/patch decoding assumes the HSC coadd object_id bit layout; patch boundaries were not drawn from a skymap.',
                     'Ephor_AB uses afterburner parent_flux_convolved fluxes that are not in the local supplement; the cause of '
                     'missing rows (input-exposure coverage, afterburner failure or table absence) is not identified.',
                     'Same-cell contrasts are descriptive quantiles without uncertainties; sources are reused across centers.',
                     'No correction, reweighting or lensing-efficiency estimate is derived.'],
        provenance=provenance, upstream_summary_sha256=UPSTREAM_SHA,
        runner_sha256=sha(Path(__file__)),
        output_sha256={p.name: sha(p) for p in (patch_table, lens_table, figure)})
    (output/'summary.json').write_text(json.dumps(report, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: report[k] for k in ('selected_objects', 'missing_objects', 'tract_patch_units',
                      'units_without_missing', 'units_fully_missing', 'units_partially_missing', 'largest_lens_share')}, indent=2))


def plot(path, objects, ids, missing, centers, lenses):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(lenses), figsize=(7.5*len(lenses), 7))
    for ax, lens in zip(np.atleast_1d(axes), lenses):
        sub = [i for i, o in enumerate(ids) if lens in objects[o]['lenses']]
        ra = np.array([objects[ids[i]]['ra'] for i in sub])
        dec = np.array([objects[ids[i]]['dec'] for i in sub])
        m = missing[sub]
        ax.scatter(ra[~m], dec[~m], s=0.3, c='0.7', label=f'finite best/risk ({(~m).sum()})', rasterized=True)
        ax.scatter(ra[m], dec[m], s=0.6, c='crimson', label=f'missing best/risk ({m.sum()})', rasterized=True)
        for c in centers:
            if c['matched_lens_id'] == lens:
                ax.plot(c['ra_deg'], c['dec_deg'], 'k*' if c['role'] == 'lens' else 'bx', ms=12)
        ax.invert_xaxis()
        ax.set_aspect(1/np.cos(np.radians(dec.mean())))
        ax.set_xlabel('RA [deg]')
        ax.set_ylabel('Dec [deg]')
        ax.set_title(f'{lens}: lens (star), fixed randoms (x)')
        ax.legend(markerscale=15, loc='lower left', fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=110, metadata={'Software': None})
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    run(parser.parse_args().output.resolve())
