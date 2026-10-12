"""Synthetic self-consistency of scaled ICL stacking and a weak-lensing mass.

The clusters, redshifts, masses and catalogue cuts below are inputs of this run.
history.json records those inputs and the numbers this process produces.

From the OGFinder root:

    python analysis/icl_wgl_membership/run_analysis.py
"""
import json
import math
import subprocess
import sys
from datetime import date
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from icl.members import (  # noqa: E402
    background_for_lensing, classify_members, diffuse_stack_masked, fraction_kept,
    object_mask)
from icl.wgl import r_delta_mpc, scaled_radius, stack_clusters, theta_delta_arcsec  # noqa: E402
from ogfkit.wlcluster import (  # noqa: E402
    ARCSEC, Cosmology, HaloModel, bin_edges, components_from_tangential, fit_mass,
    reduced_shear, select_sources)

HERE = Path(__file__).resolve().parent
S0 = 2.5
U_SCALE = 0.45
DELTA = 200.0
U_MAX = 2.2
N_SIDE = 89
DZ = 0.03
RS_WIDTH = 0.10
MAG_LIMIT = 22.5
MIN_SATELLITES = 2
P_MIN = 0.8
BLOB_COUNTS = 20000.0
LENS_Z = 0.30
SOURCE_Z = 1.0
LENS_MASS = 1.0e15
LENS_C = 4.0
RMIN = 0.3
RMAX = 3.0
N_BACKGROUND = 500
N_STARS = 180
N_SATELLITES = 40
N_PHOTO_BACKGROUND = 30


def sequence_color(mag):
    mag = np.asarray(mag, float)
    return 1.10 - 0.03 * (mag - 20.0)


def catalogue(z_cluster):
    def row(name, z, z_phot, mag, color, star, blob, pz):
        return dict(name=name, z=z, z_phot=z_phot, mag=mag, color=color, star=star,
                    blob=blob, pz=pz)

    return [
        row('bcg', z_cluster, z_cluster, 17.0, sequence_color(17.0), False, False, None),
        row('spec_red_1', z_cluster + 0.008, np.nan, 19.0, sequence_color(19.0), False, True, None),
        row('spec_red_2', z_cluster + 0.012, np.nan, 20.2, sequence_color(20.2), False, True, None),
        row('spec_blue', z_cluster - 0.010, np.nan, 23.2, sequence_color(23.2) - 0.40, False, True, None),
        row('spec_background_red', z_cluster + 0.55, z_cluster, 21.0, sequence_color(21.0), False, True, None),
        row('photo_satellite', np.nan, z_cluster + 0.005, 21.2, sequence_color(21.2) - 0.35, False, True, None),
        row('photo_background_red', np.nan, z_cluster + 0.70, 21.4, sequence_color(21.4), False, True, None),
        row('red_sequence', np.nan, z_cluster + 0.15, 20.6, sequence_color(20.6), False, True, None),
        row('unknown_blue', np.nan, z_cluster + 0.15, 21.1, sequence_color(21.1) - 0.45, False, True, None),
        row('star', z_cluster + 0.80, z_cluster + 0.80, 18.4, sequence_color(18.4), True, True, None),
        row('pz_satellite', np.nan, np.nan, 22.0, sequence_color(22.0) - 0.40, False, True, 'cluster'),
        row('faint_red', np.nan, z_cluster + 0.15, 23.6, sequence_color(23.6), False, True, None),
    ]


def classify_catalogue(rows, z_cluster):
    z = np.array([row['z'] for row in rows], float)
    z_phot = np.array([row['z_phot'] for row in rows], float)
    color = np.array([row['color'] for row in rows], float)
    mag = np.array([row['mag'] for row in rows], float)
    star = np.array([row['star'] for row in rows], bool)
    grid = np.linspace(0.0, 2.0, 401)
    pz = np.full((len(rows), grid.size), np.nan)
    for i, row in enumerate(rows):
        if row['pz'] == 'cluster':
            pz[i] = np.exp(-0.5 * ((grid - z_cluster) / 0.008) ** 2)
    return classify_members(
        z, z_cluster, DZ, bcg_index=0, star=star, z_phot=z_phot,
        cluster_window=(z_cluster - 0.02, z_cluster + 0.02),
        background_window=(z_cluster + 0.40, z_cluster + 1.20),
        pz=pz, z_grid=grid, p_min=P_MIN, color=color, mag=mag,
        rs_width=RS_WIDTH, mag_limit=MAG_LIMIT, min_satellites=MIN_SATELLITES)


def paint(spec, rows, cosmo):
    theta = theta_delta_arcsec(spec['mass'], spec['z'], cosmo, DELTA)
    pixel_scale = theta / spec['pixels_per_radius']
    blob_rows = [row for row in rows if row['blob']]
    for k, row in enumerate(blob_rows):
        row['dx'] = 14 + k * 4
        row['dy'] = 0
    centers = [(row['dx'], row['dy']) for row in blob_rows]
    for i, left in enumerate(centers):
        for right in centers[i + 1:]:
            if max(abs(left[0] - right[0]), abs(left[1] - right[1])) < 4:
                raise RuntimeError('resolved sources overlap')
    reach = int(math.ceil(U_MAX * spec['pixels_per_radius'])) + 6
    half = max(reach, max(row['dx'] for row in blob_rows) + 3)
    n = 2 * half + 1
    yy, xx = np.indices((n, n))
    theta_pix = pixel_scale * np.hypot(xx - half, yy - half)
    u = scaled_radius(theta_pix, spec['mass'], spec['z'], cosmo, DELTA)
    image = S0 * np.exp(-u / U_SCALE) * pixel_scale ** 2
    seg = np.zeros((n, n), int)
    for i, row in enumerate(rows):
        row['counts'] = 0.0
        if not row['blob']:
            continue
        block = (np.abs(xx - (half + row['dx'])) <= 1) & (np.abs(yy - (half + row['dy'])) <= 1)
        image[block] += BLOB_COUNTS
        seg[block] = i + 1
        row['counts'] = float(BLOB_COUNTS * np.count_nonzero(block))
    return dict(image=image, seg=seg, center=(half, half), pixel_scale=float(pixel_scale),
                theta_delta_arcsec=float(theta))


def radial_mean(stack):
    u = stack['u']
    uy, ux = np.meshgrid(u, u, indexing='ij')
    radius = np.hypot(ux, uy)
    edges = np.linspace(0.0, float(u.max()), 8)
    rows = []
    for i in range(edges.size - 1):
        model = S0 * np.exp(-radius / U_SCALE)
        sel = ((radius >= edges[i]) & (radius < edges[i + 1]) & (stack['weight'] > 0)
               & np.isfinite(stack['sb']))
        if not np.any(sel):
            continue
        rows.append(dict(
            u_min=float(edges[i]), u_max=float(edges[i + 1]), n=int(np.count_nonzero(sel)),
            sb_mean=float(np.mean(stack['sb'][sel])),
            model_mean=float(np.mean(model[sel]))))
    return rows


def draw_sources(cosmo, n, seed):
    rng = np.random.default_rng(seed)
    unit = rng.random(n)
    radius = np.sqrt(RMIN ** 2 + unit * (RMAX ** 2 - RMIN ** 2))
    phi = rng.random(n) * 2.0 * math.pi
    grid = np.geomspace(RMIN * 0.8, RMAX * 1.2, 80)
    sigma, excess = HaloModel(cosmo, LENS_Z, grid, two_halo=False).profiles(LENS_MASS, LENS_C)
    sigma_crit = np.full(n, float(cosmo.sigma_crit(LENS_Z, SOURCE_Z)))
    shear = reduced_shear(np.interp(radius, grid, sigma), np.interp(radius, grid, excess),
                          sigma_crit, reduced=True)
    g1, g2 = components_from_tangential(phi, shear, np.zeros(n))
    distance = float(cosmo.angular_diameter(LENS_Z))
    angle = radius / distance
    return (angle * np.cos(phi) / ARCSEC, angle * np.sin(phi) / ARCSEC, g1, g2)


def fit_lensing(cosmo, x, y, g1, g2, redshift, labels):
    selected = select_sources(
        x, y, g1, g2, redshift, np.ones(len(x)), LENS_Z, cosmo, z_gap=0.1,
        rmin=RMIN, rmax=RMAX, member=labels)
    edges = bin_edges(RMIN, RMAX, 6)
    summary, rows = fit_mass(
        selected, edges, cosmo, LENS_Z, sigma_e=1e-3, concentration='fixed',
        c_fixed=LENS_C, two_halo=False, samples=0, reduced=True, power='linear',
        lss_cov=False, min_per_bin=8)
    profile = [dict(R_Mpc=row['R'], n=row['n'], gt=row['gt'], gt_err=row['gt_err'],
                    gt_model=row.get('gt_model'), deltasigma=row['deltasigma'],
                    boost=row['boost']) for row in rows]
    kept = dict(success=summary['success'], method=summary['method'],
                M200c_Msun=summary.get('M200c_Msun'), c200=summary.get('c200'),
                chi2=summary.get('chi2'), n_src=summary['n_src'], n_bins=summary.get('n_bins'),
                dof=summary.get('dof'), covariance=summary.get('covariance'),
                power=summary.get('power'), two_halo=summary.get('two_halo'),
                profile=profile)
    return kept


def json_ready(value):
    if isinstance(value, dict):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_ready(item) for item in value]
    if isinstance(value, np.ndarray):
        return json_ready(value.tolist())
    if isinstance(value, np.floating):
        value = float(value)
    if isinstance(value, np.integer):
        return int(value)
    if isinstance(value, np.bool_):
        return bool(value)
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    return value


def parent_commit():
    try:
        return subprocess.check_output(
            ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def main():
    cosmo = Cosmology()
    expected = ['bcg', 'satellite', 'satellite', 'satellite', 'background', 'satellite',
                'background', 'satellite', 'unknown', 'star', 'satellite', 'unknown']
    clusters = [
        dict(name='A', z=0.25, mass=1.2e15, pixels_per_radius=40),
        dict(name='B', z=0.62, mass=3.0e14, pixels_per_radius=26),
    ]
    records = []
    open_records = []
    cluster_rows = []
    for spec in clusters:
        rows = catalogue(spec['z'])
        classified = classify_catalogue(rows, spec['z'])
        labels = classified['labels']
        if list(labels) != expected:
            raise RuntimeError('membership labels for %s are %s' % (spec['name'], list(labels)))
        if not classified['red_sequence']['ok']:
            raise RuntimeError('red sequence did not fit for %s' % spec['name'])
        painted = paint(spec, rows, cosmo)
        diffuse = object_mask(painted['seg'], labels, 'diffuse')
        records.append(dict(
            image=painted['image'], center=painted['center'], pixel_scale=painted['pixel_scale'],
            mass=spec['mass'], z=spec['z'], mask=diffuse))
        open_records.append(dict(
            image=painted['image'], center=painted['center'], pixel_scale=painted['pixel_scale'],
            mass=spec['mass'], z=spec['z']))
        kept = fraction_kept(labels)
        flux_kept = float(sum(row['counts'] for row, flag in zip(rows, kept) if flag))
        flux_out = float(sum(row['counts'] for row, flag in zip(rows, kept) if not flag))
        masked = diffuse_stack_masked(labels)
        cluster_rows.append(dict(
            name=spec['name'], z=spec['z'], mass_Msun=spec['mass'],
            pixels_per_radius=spec['pixels_per_radius'],
            pixel_scale_arcsec=painted['pixel_scale'],
            theta_delta_arcsec=painted['theta_delta_arcsec'],
            r_delta_mpc=float(r_delta_mpc(spec['mass'], spec['z'], cosmo, DELTA)),
            red_sequence=classified['red_sequence'],
            labels=[dict(name=row['name'], label=str(label), counts=row['counts'])
                    for row, label in zip(rows, labels)],
            n_diffuse_masked=int(np.count_nonzero(masked)),
            fraction_counts_kept=flux_kept,
            fraction_counts_rejected=flux_out))

    stacked = stack_clusters(records, cosmo, DELTA, U_MAX, N_SIDE)
    opened = stack_clusters(open_records, cosmo, DELTA, U_MAX, N_SIDE)
    uy, ux = np.meshgrid(stacked['u'], stacked['u'], indexing='ij')
    radius = np.hypot(ux, uy)
    model = S0 * np.exp(-radius / U_SCALE)
    comparable = (stacked['weight'] > 0) & np.isfinite(stacked['sb']) & (model > 1e-3)
    residual = np.abs(stacked['sb'][comparable] - model[comparable]) / model[comparable]
    max_rel = float(np.max(residual))
    if max_rel >= 0.02:
        raise RuntimeError('scaled stack residual %s' % max_rel)
    open_peak = float(np.nanmax(opened['sb']))
    masked_peak = float(np.nanmax(stacked['sb']))
    if open_peak <= 5.0 * S0:
        raise RuntimeError('unmasked stack did not keep the resolved sources')

    x, y, g1, g2 = draw_sources(cosmo, N_BACKGROUND, 17)
    xs, ys, _g1s, _g2s = draw_sources(cosmo, N_STARS, 19)
    xt, yt, _g1t, _g2t = draw_sources(cosmo, N_SATELLITES, 23)
    x = np.concatenate([x, xs, xt])
    y = np.concatenate([y, ys, yt])
    g1 = np.concatenate([g1, np.zeros(N_STARS + N_SATELLITES)])
    g2 = np.concatenate([g2, np.zeros(N_STARS + N_SATELLITES)])
    redshift = np.concatenate([
        np.full(N_BACKGROUND, SOURCE_Z), np.full(N_STARS, SOURCE_Z),
        np.full(N_SATELLITES, LENS_Z + 0.02)])
    star = np.concatenate([
        np.zeros(N_BACKGROUND, bool), np.ones(N_STARS, bool), np.zeros(N_SATELLITES, bool)])
    labels = classify_members(redshift, LENS_Z, 0.05, star=star)['labels']
    labelled = fit_lensing(cosmo, x, y, g1, g2, redshift, labels)
    unlabelled = fit_lensing(cosmo, x, y, g1, g2, redshift, None)
    if not labelled['success'] or not unlabelled['success']:
        raise RuntimeError('mass fit did not succeed')
    log_ratio = math.log10(labelled['M200c_Msun'] / LENS_MASS)
    if abs(log_ratio) >= 0.02:
        raise RuntimeError('labelled mass ratio %s dex' % log_ratio)
    dilution = math.log10(unlabelled['M200c_Msun'] / labelled['M200c_Msun'])
    if dilution >= -0.05:
        raise RuntimeError('stars with background redshifts did not move the mass (%s dex)' % dilution)
    if unlabelled['n_src'] - labelled['n_src'] != N_STARS:
        raise RuntimeError('star rows were not the only extra lensing sources')

    photo_z = np.concatenate([
        np.full(N_PHOTO_BACKGROUND, 1.0), np.full(10, LENS_Z)])
    photo_labels = classify_members(
        np.full(photo_z.size, np.nan), LENS_Z, 0.05, z_phot=photo_z,
        cluster_window=(LENS_Z - 0.02, LENS_Z + 0.02),
        background_window=(0.50, 1.50))['labels']
    photo_keep = background_for_lensing(photo_labels, np.full(photo_z.size, np.nan), LENS_Z, 0.1)
    if int(np.count_nonzero(photo_keep)) != 0:
        raise RuntimeError('photometric labels entered the lensing sample without a redshift')
    n_photo_background = int(np.count_nonzero(photo_labels[:N_PHOTO_BACKGROUND] == 'background'))
    if n_photo_background != N_PHOTO_BACKGROUND:
        raise RuntimeError('photometric background window was not applied')

    history = dict(
        product='Astrafex Desktop',
        analysis='scaled ICL stack, ICL light fraction, weak-lensing mass, spectroscopic and photometric membership',
        date=date.today().isoformat(),
        kind='synthetic self-consistency',
        statement=(
            'The numbers in this file are the output of run_analysis.py on the synthetic '
            'clusters defined in that script. They record this run. A paper cites the '
            'repository path below and the commit that contains this directory.'),
        cosmology=dict(H0=cosmo.H0, Om=cosmo.Om, ns=cosmo.ns, sigma8=cosmo.sigma8,
                       note='Caller settings of Cosmology(), the same defaults as the weak-lensing library.'),
        delta=DELTA,
        membership=dict(
            spectroscopic_half_width=DZ,
            spectroscopic_boundary='|z - z_cluster| = dz is outside the cluster',
            photoz_cluster_window='[z_cluster - 0.02, z_cluster + 0.02)',
            photoz_background_window='[z_cluster + 0.40, z_cluster + 1.20)',
            pz_p_min=P_MIN,
            red_sequence_width=RS_WIDTH,
            red_sequence_magnitude_limit=MAG_LIMIT,
            red_sequence_min_satellites=MIN_SATELLITES,
            red_sequence_form='colour = a + b (mag - m0), fit to spectroscopic satellites only'),
        icl=dict(
            surface_brightness='S0 exp(-u / u_scale), counts per arcsec^2',
            S0=S0, u_scale=U_SCALE, u_max=U_MAX, n_side=N_SIDE,
            max_relative_residual=max_rel,
            masked_peak_sb=masked_peak,
            unmasked_peak_sb=open_peak,
            radial_profile=radial_mean(stacked),
            clusters=cluster_rows),
        weak_lensing=dict(
            z_lens=LENS_Z, z_source=SOURCE_Z, input_M200c_Msun=LENS_MASS,
            input_c200=LENS_C, rmin_Mpc=RMIN, rmax_Mpc=RMAX, z_gap=0.1,
            samples=0, two_halo=False, power='linear', concentration='fixed',
            n_background=N_BACKGROUND, n_stars=N_STARS, n_satellites=N_SATELLITES,
            labelled=labelled, unlabelled=unlabelled,
            log10_labelled_over_input=log_ratio,
            log10_unlabelled_over_labelled=dilution,
            n_photometric_background_without_spectroscopic_redshift=n_photo_background,
            photometric_background_kept_for_lensing=int(np.count_nonzero(photo_keep))),
        data_availability=dict(
            repository='https://github.com/kjhan0606/OGFinder',
            path='analysis/icl_wgl_membership/',
            script='analysis/icl_wgl_membership/run_analysis.py',
            history='analysis/icl_wgl_membership/history.json',
            modules=['icl/wgl.py', 'icl/members.py', 'ogfkit/wlcluster.py'],
            parent_commit=parent_commit(),
            parent_commit_note=(
                'HEAD when this file was written. The citation commit is the commit '
                'that adds this directory.')))
    path = HERE / 'history.json'
    path.write_text(json.dumps(json_ready(history), indent=2) + '\n', encoding='utf-8')
    print('wrote %s' % path)
    print('icl max relative residual %.6f' % max_rel)
    print('icl masked peak %.6f unmasked peak %.6f' % (masked_peak, open_peak))
    print('lensing labelled M200c %.6e log10 ratio %+.5f' % (labelled['M200c_Msun'], log_ratio))
    print('lensing unlabelled M200c %.6e log10 shift %+.5f' % (unlabelled['M200c_Msun'], dilution))
    print('lensing n_src labelled %d unlabelled %d' % (labelled['n_src'], unlabelled['n_src']))


if __name__ == '__main__':
    main()
