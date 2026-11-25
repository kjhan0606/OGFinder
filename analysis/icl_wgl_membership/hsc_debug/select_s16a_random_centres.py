#!/usr/bin/env python3
"""Deterministic S16A random-centre selector for the 11 pilot lenses (local only; never contacts CAS).

Candidates are the deduplicated thinned S16A random points of the approved candidate query. Rules (fixed here):
- geometry: 2 R <= d(lens) <= 90 arcmin - R, R = the lens 3-proper-Mpc radius; a disc may not overlap any pilot
  lens disc (separation > R + R_lens) or any already accepted disc of any lens (> R + R_other; within one lens
  this is the 6-Mpc pairwise separation);
- CAMIRA: zero overlap of the 3-Mpc disc with any S20A or S23B CAMIRA cluster's 3-proper-Mpc aperture at the
  cluster redshift (the WL random-coordinate screen rule, both catalogues unioned);
- footprint: the FDFC grid (disc FDFC fraction >= 0.95) is required. Without it the selector only runs with
  --provisional-without-fdfc, writes 'provisional_not_for_queries' and uses an i-processed point-coverage proxy;
- masks: never rejects on the centre point or on mask fraction (reported); never uses S16A source counts;
- ranking: per-disc medians of icountinputs, i PSF FWHM (SDSS moments) and isky_std versus the lens disc;
  score = sum of squared differences scaled by the tolerances; candidates outside any tolerance are dropped;
- placement: round-robin over lenses in ID order, each taking its best compatible candidate (score, then
  sha256(lens|object_id)), until TARGET per lens or exhaustion; a lens with fewer than MINIMUM is reported.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree

from audit_s16a_random_point_areas import CAND, hole, i_null, inner_ok, load_rows, truthy, unit
from audit_s23b_icl_random_control_candidates import read_camira
from prepare_s16a_random_candidate_rows import ADJUST_MAX, DENSITY_PER_ARCMIN2, RADIUS_ARCMIN, centres

DATA = Path('/gpfs/kjhan/Astrafek/data')
CAMIRA = {DATA/'catalog/camira_s20a_wide_sm.tbl': 'd929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa',
          DATA/'catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl': '5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71'}
FDFC_GRID = DATA/'sample_800_2026-10-13/s16a_fdfc_grid_plan_2026-11-25_01'
H0, OM0 = 70.0, 0.3
CAMIRA_MPC = 3.0
TARGET, MINIMUM = 8, 5
FDFC_MIN = 0.95
COVERAGE_MIN = 0.8  # provisional proxy: deduplicated i-processed points / expected (0.3 per arcmin^2)
TOL = dict(icountinputs=2.0, psf_fwhm_ratio=1.2, sky_ratio=1.3)
PIXEL = 0.168
FWHM = 2*math.sqrt(2*math.log(2))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sep_deg(u, v):
    return np.degrees(2*np.arcsin(np.clip(np.linalg.norm(u-v, axis=-1)/2, 0, 1)))


def chord(deg):
    return 2*math.sin(math.radians(deg)/2)


def load_points():
    _, rows = load_rows(CAND)
    seen, pts = set(), []
    for r in rows:
        if r['object_id'] in seen or i_null(r) or not inner_ok(r):
            continue
        seen.add(r['object_id'])
        xx, yy = float(r['ishape_sdss_psf_xx'] or 'nan'), float(r['ishape_sdss_psf_yy'] or 'nan')
        pts.append(dict(lens_id=r['lens_id'], object_id=r['object_id'], ra=float(r['ra']), dec=float(r['dec']),
                        icountinputs=int(r['icountinputs']), hole=hole(r), sky=float(r['isky_std'] or 'nan'),
                        psf=FWHM*PIXEL*math.sqrt((xx+yy)/2) if xx > 0 and yy > 0 else math.nan))
    return pts


def disc_stats(idx, pts, area):
    ci = np.array([pts[i]['icountinputs'] for i in idx], float)
    psf = np.array([pts[i]['psf'] for i in idx]); sky = np.array([pts[i]['sky'] for i in idx])
    return dict(n_points=len(idx), coverage=len(idx)/(DENSITY_PER_ARCMIN2*ADJUST_MAX*area),
                hole_fraction=float(np.mean([pts[i]['hole'] for i in idx])) if len(idx) else None,
                icountinputs=float(np.median(ci)) if len(idx) else math.nan,
                psf_fwhm=float(np.nanmedian(psf)) if np.isfinite(psf).any() else math.nan,
                sky=float(np.nanmedian(sky)) if np.isfinite(sky).any() else math.nan)


def fdfc_fraction(lens, ra, dec, radius_deg, grid):
    """Share of a 0.25-arcmin lattice over the disc lying in FDFC cells of the lens grid."""
    ra0, dec0 = float(lens['ra_deg']), float(lens['dec_deg'])
    step = 0.25/60
    n = int(math.ceil(radius_deg/step))
    ax = np.arange(-n, n+1)*step
    x, y = np.meshgrid(ax, ax)
    keep = x**2+y**2 <= radius_deg**2
    d = dec+y[keep]; a = ra+x[keep]/math.cos(math.radians(dec))
    cx = np.floor((a-ra0)*math.cos(math.radians(dec0))*60).astype(int); cy = np.floor((d-dec0)*60).astype(int)
    return float(np.mean([grid.get((i, j), False) for i, j in zip(cx, cy)]))


def load_fdfc():
    """FDFC cells per lens from a completed, receipt-verified grid download; None if absent."""
    summary = FDFC_GRID/'results'/'download_summary.json'
    if not summary.is_file():
        return None
    s = json.loads(summary.read_text())
    from prepare_s16a_fdfc_grid import FDFC_FRACTION, MIN_CELL_OBJECTS
    grids = {}
    for name, digest in s['output_sha256'].items():
        path = FDFC_GRID/'results'/name
        if sha(path) != digest:
            raise ValueError('FDFC grid CSV changed '+name)
        with path.open() as stream:
            reader = csv.reader(stream)
            head = [h.lstrip('# ').strip() for h in next(reader)]
            for r in reader:
                r = dict(zip(head, r))
                n, f = int(r['n_primary']), int(r['n_fdfc'])
                grids.setdefault(r['lens_id'], {})[(int(r['cell_x']), int(r['cell_y']))] = n >= MIN_CELL_OBJECTS and f/n >= FDFC_FRACTION
    return dict(summary_sha256=sha(summary), grids=grids)


def run(output, provisional):
    fdfc = load_fdfc()
    if fdfc is None and not provisional:
        raise SystemExit('FDFC grid not available: the bright-object rule requires FDFC before the selector is hashed '
                         '(use --provisional-without-fdfc for a non-binding feasibility run)')
    if output.exists():
        raise FileExistsError(output)
    for path, digest in CAMIRA.items():
        if sha(path) != digest:
            raise ValueError('CAMIRA catalogue changed '+str(path))
    cosmo = FlatLambdaCDM(H0=H0, Om0=OM0)
    camira = [c for path in CAMIRA for c in read_camira(path)]
    cu = unit([c['ra_deg'] for c in camira], [c['dec_deg'] for c in camira])
    crad = np.degrees(CAMIRA_MPC/cosmo.angular_diameter_distance(np.array([c['z'] for c in camira])).value)
    ctree = cKDTree(cu)
    ledger = centres()
    lenses = {r['matched_lens_id']: r for r in ledger if r['role'] == 'lens'}
    lu = {g: unit([float(r['ra_deg'])], [float(r['dec_deg'])])[0] for g, r in lenses.items()}
    lr = {g: float(r['radius_arcsec'])/3600 for g, r in lenses.items()}
    pts = load_points()
    pu = unit([p['ra'] for p in pts], [p['dec'] for p in pts])
    ptree = cKDTree(pu)
    pools, lens_stats, rejected = {}, {}, {}
    for g in sorted(lenses):
        R = lr[g]
        area = math.pi*(R*60)**2
        lens_stats[g] = disc_stats(ptree.query_ball_point(lu[g], chord(R)), pts, area)
        counts = dict(geometry=0, lens_overlap=0, camira=0, footprint=0, tolerance=0, accepted=0)
        pool = []
        for i, p in enumerate(pts):
            if p['lens_id'] != g:
                continue
            d = float(sep_deg(pu[i], lu[g]))
            if not (2*R <= d <= RADIUS_ARCMIN/60-R):
                counts['geometry'] += 1; continue
            if any(float(sep_deg(pu[i], lu[h])) <= R+lr[h] for h in lenses if h != g):
                counts['lens_overlap'] += 1; continue
            near = np.array(ctree.query_ball_point(pu[i], chord(R+crad.max())), int)
            if len(near) and np.any(sep_deg(pu[i], cu[near]) <= R+crad[near]):
                counts['camira'] += 1; continue
            st = disc_stats(ptree.query_ball_point(pu[i], chord(R)), pts, area)
            if fdfc is not None:
                st['fdfc_fraction'] = fdfc_fraction(lenses[g], p['ra'], p['dec'], R, fdfc['grids'].get(g, {}))
                ok = st['fdfc_fraction'] >= FDFC_MIN
            else:
                ok = st['coverage'] >= COVERAGE_MIN
            if not ok:
                counts['footprint'] += 1; continue
            L = lens_stats[g]
            dci = st['icountinputs']-L['icountinputs']
            rp = st['psf_fwhm']/L['psf_fwhm']; rs = st['sky']/L['sky']
            if not (abs(dci) <= TOL['icountinputs'] and 1/TOL['psf_fwhm_ratio'] <= rp <= TOL['psf_fwhm_ratio']
                    and 1/TOL['sky_ratio'] <= rs <= TOL['sky_ratio']):
                counts['tolerance'] += 1; continue
            score = (dci/TOL['icountinputs'])**2+(math.log(rp)/math.log(TOL['psf_fwhm_ratio']))**2+(math.log(rs)/math.log(TOL['sky_ratio']))**2
            tie = hashlib.sha256(f"{g}|{p['object_id']}".encode()).hexdigest()
            pool.append(dict(index=i, score=score, tie=tie, **st))
            counts['accepted'] += 1
        pools[g] = sorted(pool, key=lambda c: (c['score'], c['tie']))
        rejected[g] = counts
    chosen, placed = {g: [] for g in lenses}, []
    active = True
    while active:
        active = False
        for g in sorted(lenses):
            if len(chosen[g]) >= TARGET:
                continue
            while pools[g]:
                c = pools[g].pop(0)
                if all(float(sep_deg(pu[c['index']], pu[j])) > lr[g]+r for j, r in placed):
                    chosen[g].append(c); placed.append((c['index'], lr[g])); active = True
                    break
    rows = []
    for g in sorted(lenses):
        for rank, c in enumerate(chosen[g]):
            p = pts[c['index']]
            rows.append(dict(center_id='RC-'+hashlib.sha256(f"s16a-wl-random|{g}|{p['object_id']}".encode()).hexdigest()[:20],
                             lens_id=g, rank=rank, random_object_id=p['object_id'], ra_deg=repr(p['ra']), dec_deg=repr(p['dec']),
                             z_lens=lenses[g]['z_lens'], radius_arcsec=lenses[g]['radius_arcsec'],
                             offset_arcmin=float(sep_deg(pu[c['index']], lu[g]))*60, score=c['score'],
                             **{k: c[k] for k in ('n_points', 'coverage', 'hole_fraction', 'icountinputs', 'psf_fwhm', 'sky')},
                             fdfc_fraction=c.get('fdfc_fraction')))
    output.mkdir(parents=True)
    with (output/'selected_centres.tsv').open('w', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    summary = dict(
        schema='astrafek-s16a-random-centre-selection-v1',
        status='provisional_not_for_queries' if fdfc is None else 'selected_pending_user_approval',
        footprint='i-processed coverage proxy (FDFC grid absent)' if fdfc is None else 'FDFC grid',
        fdfc_grid_summary_sha256=None if fdfc is None else fdfc['summary_sha256'],
        constants=dict(h0=H0, om0=OM0, camira_mpc=CAMIRA_MPC, target=TARGET, minimum=MINIMUM, fdfc_min=FDFC_MIN,
                       coverage_min=COVERAGE_MIN, tolerances=TOL, cone_arcmin=RADIUS_ARCMIN),
        per_lens={g: dict(selected=len(chosen[g]), below_minimum=len(chosen[g]) < MINIMUM, candidate_filter=rejected[g],
                          lens_disc=lens_stats[g]) for g in sorted(lenses)},
        n_selected=len(rows), camira_sha256={str(k): v for k, v in CAMIRA.items()},
        candidate_receipt_sha256=sha(CAND/'results'/'download_summary.json'),
        selector_sha256=sha(Path(__file__)), selected_centres_sha256=sha(output/'selected_centres.tsv'),
        science_ready=False)
    (output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    return summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--provisional-without-fdfc', action='store_true')
    args = parser.parse_args()
    s = run(args.output, args.provisional_without_fdfc)
    print(json.dumps({g: (v['selected'], v['candidate_filter']) for g, v in s['per_lens'].items()}, indent=1))
    print(s['status'], s['n_selected'])
