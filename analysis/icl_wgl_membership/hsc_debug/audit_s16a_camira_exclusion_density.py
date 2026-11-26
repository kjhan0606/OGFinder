#!/usr/bin/env python3
"""Does the 2-Mpc CAMIRA exclusion place S16A random discs in underdense sky?  Local only, no new queries.

Density tracer: FDFC objects per FDFC cell (1 arcmin cells of the accepted grids; the four wide lenses use their
checked 180-arcmin grid), averaged over a disc of the lens's 3-Mpc radius by a disc-kernel convolution.
Candidate disc centres are every cell centre allowed by the selector geometry (2R <= d <= cone - R, no overlap
with another lens disc) with cell-based FDFC fraction >= 0.95. For each candidate the selector's CAMIRA
clearance (S20A and S23B catalogues, min over clusters of (sep - R) D_A(z_c)) is computed. Per lens, densities
are normalised by the mean over all candidates. Reported: CAMIRA-pass and -fail means, the 76 selected centres,
the pilot random centres, the lens disc and its 2-3 Mpc annulus; and the tracer against the WL source density of
the random centres (proxy check). Candidate discs overlap, so per-lens means are not independent samples; the
spread over the 11 lenses is the uncertainty used.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.signal import fftconvolve
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

import select_s16a_random_centres as sel
from audit_s23b_icl_random_control_candidates import read_camira

BOOSTS = Path('/home/kjhan/BACKUP/Astrafek/s16a_random_extension_boosts_20261126_02')
SELECTION = Path('/home/kjhan/BACKUP/Astrafek/s16a_random_centre_selection_20261126_01')
CAMIRA_MPC = 2.0
OUTER_MPC = (2.0, 3.0)


def raw_cells(path_dir, summary_name='download_summary.json', kind=None):
    s = json.loads((path_dir/'results'/summary_name).read_text())
    out = defaultdict(dict)
    for name, digest in s['output_sha256'].items():
        if kind and s['per_batch'][Path(name).stem]['kind'] != kind:
            continue
        path = path_dir/'results'/name
        if sel.sha(path) != digest:
            raise ValueError('grid CSV changed '+name)
        with path.open() as stream:
            reader = csv.reader(stream)
            head = [h.lstrip('# ').strip() for h in next(reader)]
            for r in reader:
                r = dict(zip(head, r))
                out[r['lens_id']][(int(r['cell_x']), int(r['cell_y']))] = (int(r['n_primary']), int(r['n_fdfc']))
    return out, sel.sha(path_dir/'results'/summary_name)


def disc_kernel(radius_cells):
    n = int(math.ceil(radius_cells))
    ax = np.arange(-n, n+1)
    x, y = np.meshgrid(ax, ax)
    return (x**2+y**2 <= radius_cells**2).astype(float)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fdfc = sel.load_fdfc()
    acc = json.loads(sel.ACCEPTANCE.read_text())
    if not (acc['accepted'] and acc['grid_summary_sha256'] == fdfc['summary_sha256']):
        raise ValueError('FDFC grid not accepted')
    wide = sel.load_wide()
    sel.merge_wide_fdfc(fdfc['grids'], wide)  # raises on any disagreement
    base, base_sha = raw_cells(sel.FDFC_GRID)
    wide_cells, wide_sha = raw_cells(sel.WIDE, kind='fdfc')
    cells = {**base, **wide_cells}
    for path, digest in sel.CAMIRA.items():
        if sel.sha(path) != digest:
            raise ValueError('CAMIRA catalogue changed')
    cosmo = FlatLambdaCDM(H0=sel.H0, Om0=sel.OM0)
    camira = [c for path in sel.CAMIRA for c in read_camira(path)]
    cu = sel.unit([c['ra_deg'] for c in camira], [c['dec_deg'] for c in camira])
    da_c = cosmo.angular_diameter_distance(np.array([c['z'] for c in camira])).value
    ctree = cKDTree(cu)
    lenses = {r['matched_lens_id']: r for r in sel.centres() if r['role'] == 'lens'}
    pilot_rand = defaultdict(list)
    for r in sel.centres():
        if r['role'] == 'random':
            pilot_rand[r['matched_lens_id']].append(r)
    lu = {g: sel.unit([float(r['ra_deg'])], [float(r['dec_deg'])])[0] for g, r in lenses.items()}
    lr = {g: float(r['radius_arcsec'])/3600 for g, r in lenses.items()}
    cone = {g: sel.WIDE_ARCMIN if g in wide['lenses'] else sel.RADIUS_ARCMIN for g in lenses}
    with (SELECTION/'selected_centres.tsv').open() as stream:
        chosen = list(csv.DictReader(stream, delimiter='\t'))
    with (BOOSTS/'random_depth_check.tsv').open() as stream:
        wl = {r['center_id']: r for r in csv.DictReader(stream, delimiter='\t')}

    per_lens, centre_rows = [], []
    for g in sorted(lenses):
        L = lenses[g]; ra0, dec0 = float(L['ra_deg']), float(L['dec_deg'])
        keys = np.array(list(cells[g]))
        off = keys.min(0); shape = keys.max(0)-off+1
        n_obj, is_fd = np.zeros(shape), np.zeros(shape)
        for (x, y), (n, f) in cells[g].items():
            if sel.fdfc_cell(n, f):
                n_obj[x-off[0], y-off[1]] = f; is_fd[x-off[0], y-off[1]] = 1
        Rm = lr[g]*60
        k = disc_kernel(Rm)
        s_obj = fftconvolve(n_obj, k, mode='same'); s_fd = fftconvolve(is_fd, k, mode='same')
        frac = s_fd/k.sum()
        dens = np.where(s_fd > 0.5, s_obj/np.maximum(s_fd, 1), np.nan)

        def at(ra, dec):
            i = int(math.floor((ra-ra0)*math.cos(math.radians(dec0))*60))-off[0]
            j = int(math.floor((dec-dec0)*60))-off[1]
            return (float(dens[i, j]), float(frac[i, j])) if 0 <= i < shape[0] and 0 <= j < shape[1] else (math.nan, 0.0)

        xs, ys = np.meshgrid(np.arange(shape[0])+off[0]+0.5, np.arange(shape[1])+off[1]+0.5, indexing='ij')
        dec = dec0+ys/60; ra = ra0+xs/60/math.cos(math.radians(dec0))
        cu_cell = sel.unit(ra.ravel(), dec.ravel())
        d = sel.sep_deg(cu_cell, lu[g])
        ok = (d >= 2*lr[g]) & (d <= cone[g]/60-lr[g]) & (frac.ravel() >= sel.FDFC_MIN) & np.isfinite(dens.ravel())
        for h in lenses:
            if h != g:
                ok &= sel.sep_deg(cu_cell, lu[h]) > lr[g]+lr[h]
        idx = np.flatnonzero(ok)
        reach = sel.chord(lr[g]+np.degrees(CAMIRA_MPC/da_c.min()))
        clear = np.full(len(idx), np.inf)
        for n, i in enumerate(idx):
            near = np.array(ctree.query_ball_point(cu_cell[i], reach), int)
            if len(near):
                clear[n] = float(np.min(np.radians(np.maximum(sel.sep_deg(cu_cell[i], cu[near])-lr[g], 0))*da_c[near]))
        D = dens.ravel()[idx]; ref = float(D.mean())
        passing = clear > CAMIRA_MPC
        # lens disc and 2-3 Mpc annulus densities from the cells themselves
        rr = np.hypot(xs, ys)/60  # degrees from the lens to cell centres (tangent approximation)
        mpc = np.radians(rr)*float(cosmo.angular_diameter_distance(float(L['z_lens'])).value)
        fd = is_fd.astype(bool)
        disc = fd & (mpc <= 3.0); ann = fd & (mpc > OUTER_MPC[0]) & (mpc <= OUTER_MPC[1])
        lens_disc = float(n_obj[disc].sum()/disc.sum())/ref if disc.sum() else None
        lens_ann = float(n_obj[ann].sum()/ann.sum())/ref if ann.sum() else None
        sel_d = [at(float(c['ra_deg']), float(c['dec_deg']))[0]/ref for c in chosen if c['lens_id'] == g]
        pil = [(r, *at(float(r['ra_deg']), float(r['dec_deg']))) for r in pilot_rand[g]]
        pil_d = [v/ref for r, v, f in pil if f >= sel.FDFC_MIN]
        for c in chosen:
            if c['lens_id'] == g:
                v, f = at(float(c['ra_deg']), float(c['dec_deg']))
                w = wl.get(c['center_id'], {})
                centre_rows.append(dict(center_id=c['center_id'], lens_id=g, set='new', tracer_density=v, tracer_ratio=v/ref,
                                        cell_fdfc_fraction=f, wl_density=w.get('density_per_arcmin2'), camira_clearance_mpc=c['camira_clearance_mpc']))
        for r, v, f in pil:
            w = wl.get(r['center_id'], {})
            centre_rows.append(dict(center_id=r['center_id'], lens_id=g, set='pilot', tracer_density=v, tracer_ratio=v/ref,
                                    cell_fdfc_fraction=f, wl_density=w.get('density_per_arcmin2'), camira_clearance_mpc=None))
        per_lens.append(dict(lens_id=g, cone_arcmin=cone[g], n_candidates=len(idx), n_camira_pass=int(passing.sum()),
                             mean_density_all=ref, pass_over_all=float(D[passing].mean()/ref) if passing.any() else None,
                             fail_over_all=float(D[~passing].mean()/ref) if (~passing).any() else None,
                             pass_over_fail=float(D[passing].mean()/D[~passing].mean()) if passing.any() and (~passing).any() else None,
                             selected_over_all=float(np.mean(sel_d)), pilot_fdfc_randoms_over_all=float(np.mean(pil_d)) if pil_d else None,
                             n_pilot_fdfc_randoms=len(pil_d), lens_disc_over_all=lens_disc, lens_2_3mpc_over_all=lens_ann))

    def spread(key):
        v = np.array([r[key] for r in per_lens if r[key] is not None], float)
        return dict(n=len(v), mean=float(v.mean()), se=float(v.std(ddof=1)/math.sqrt(len(v))), min=float(v.min()), max=float(v.max()))
    proxy = [(float(r['tracer_density']), float(r['wl_density'])) for r in centre_rows
             if r['wl_density'] not in (None, '', 'None') and math.isfinite(r['tracer_density']) and r['cell_fdfc_fraction'] >= sel.FDFC_MIN]
    summary = dict(schema='astrafek-s16a-camira-exclusion-density-v1', science_measurement=False, selection_change=False,
                   tracer='FDFC objects per FDFC 1-arcmin cell, disc mean', camira_mpc=CAMIRA_MPC,
                   fdfc_grid_summary_sha256=base_sha, wide_summary_sha256=wide_sha,
                   across_lenses={k: spread(k) for k in ('pass_over_all', 'fail_over_all', 'pass_over_fail', 'selected_over_all',
                                                          'pilot_fdfc_randoms_over_all', 'lens_disc_over_all', 'lens_2_3mpc_over_all')},
                   proxy_check=dict(n=len(proxy), spearman=spearmanr(*zip(*proxy))._asdict() if len(proxy) > 2 else None),
                   limitations=['Tracer counts all FDFC primary objects (stars and galaxies, all magnitudes), not WL sources.',
                                'Cell-based FDFC fraction (1-arcmin cells) approximates the selector 0.25-arcmin lattice.',
                                'Lens annulus densities use cell centres; cluster members inflate the lens disc.'],
                   runner_sha256=sel.sha(Path(__file__)))
    output.mkdir(parents=True)
    for fname, tab in (('per_lens.tsv', per_lens), ('centres.tsv', centre_rows)):
        with (output/fname).open('w', newline='') as stream:
            w = csv.DictWriter(stream, fieldnames=list(tab[0]), delimiter='\t', lineterminator='\n')
            w.writeheader(); w.writerows(tab)
    summary['output_sha256'] = {n: sel.sha(output/n) for n in ('per_lens.tsv', 'centres.tsv')}
    (output/'summary.json').write_text(json.dumps(summary, indent=2, default=float)+'\n')
    return per_lens, summary


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('/home/kjhan/BACKUP/Astrafek/s16a_camira_exclusion_density_20261126_01'))
    per_lens, summary = run(parser.parse_args().output)
    for r in per_lens:
        print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
    print(json.dumps({k: summary[k] for k in ('across_lenses', 'proxy_check')}, indent=1, default=float))
