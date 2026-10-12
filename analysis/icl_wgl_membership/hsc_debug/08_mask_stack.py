#!/usr/bin/env python3
"""Mask the BCG and CAMIRA members, then average the frames on a kpc grid.

frames.tsv columns are id, filename, sky, product. sky is global or
global-sky. product is empty, deepCoadd, or coadd/bg. local-sky and calexp
are refused. The image product is copied into the log and is not renamed.

Membership comes from logs/members.tsv. One HSC-I frame does not assign
membership. SEP supplies the pixel footprint of a catalog position: a
detection within 1 arcsec is masked with an ellipse of 2.5 Kron radii, and
a position without a detection is masked with a circle of 2 arcsec. Other
detections are written to the detection table and stay in the image.

The stack calls icl.wgl.warp_metric and mean_stack. It does not call
stack_clusters, does not read a mass, and does not use the library product
name bcg+diffuse. The stack product is bcg-and-members-masked. The survey
checklist is left empty and science_measurement stays false.
"""
import json
import math
import os
import sys

os.environ.setdefault('OMP_NUM_THREADS', '2')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '2')
os.environ.setdefault('MKL_NUM_THREADS', '2')

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')
OGF = os.environ.get('ASTRAFEX_OGF_ROOT', '/home/kjhan/BACKUP/Astrafek/OGFinder')
KPC_MAX = 500.0
N_SIDE = 201
H0 = 70.0
OM = 0.3
MATCH_ARCSEC = 1.0
FALLBACK_ARCSEC = 2.0
KRON_FACTOR = 2.5
SEP_THRESH = 1.5
ALLOWED_SKY = ('global', 'global-sky')
ALLOWED_PRODUCT = ('', 'deepCoadd', 'coadd/bg')


def write_status(payload, code):
    path = os.path.join(ROOT, 'logs', 'mask_stack_status.json')
    os.makedirs(os.path.dirname(path), exist_ok=True)
    payload['science_measurement'] = False
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=1)
        handle.write('\n')
    print('mask stack %s' % payload['status'])
    return code


def read_tsv(path):
    with open(path, encoding='utf-8') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
    return header, rows


def pixel_scale_arcsec(header):
    if 'CD1_1' in header and 'CD2_1' in header:
        return 3600.0 * math.hypot(float(header['CD1_1']), float(header['CD2_1']))
    if 'CDELT1' in header:
        return 3600.0 * abs(float(header['CDELT1']))
    raise SystemExit('pixel scale is missing from the cutout header')


def world_xy(header, ra, dec):
    from astropy.coordinates import SkyCoord
    from astropy.wcs import WCS
    import astropy.units as u
    x, y = WCS(header).world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
    return float(x), float(y)


def first_image(path):
    from astropy.io import fits
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if hdu.data is not None and getattr(hdu.data, 'ndim', 0) == 2:
                return np_image(hdu.data), hdu.header
    raise SystemExit('no 2-D image in %s' % path)


def np_image(data):
    import numpy as np
    return np.array(data, dtype=float, copy=True)


def _paint(mask, x, y, a, b, theta, circle):
    import numpy as np
    a = max(float(a), 1.0)
    b = max(float(b), 1.0)
    radius = int(math.ceil(max(a, b))) + 1
    x0 = int(math.floor(x))
    y0 = int(math.floor(y))
    y1 = max(0, y0 - radius)
    y2 = min(mask.shape[0], y0 + radius + 1)
    x1 = max(0, x0 - radius)
    x2 = min(mask.shape[1], x0 + radius + 1)
    if y1 >= y2 or x1 >= x2:
        return
    yy, xx = np.mgrid[y1:y2, x1:x2]
    dx = xx - x
    dy = yy - y
    if circle:
        mask[y1:y2, x1:x2] |= (dx * dx + dy * dy) <= a * a
        return
    ct = math.cos(theta)
    st = math.sin(theta)
    xr = dx * ct + dy * st
    yr = -dx * st + dy * ct
    mask[y1:y2, x1:x2] |= (xr / a) ** 2 + (yr / b) ** 2 <= 1.0


def build_mask(image, pixel_scale, positions):
    """Return a True-excluded mask and a detection table for one frame."""
    import numpy as np
    import sep
    mask = np.zeros(image.shape, dtype=bool)
    finite = np.isfinite(image)
    if not np.any(finite):
        raise SystemExit('image has no finite pixels')
    work = np.ascontiguousarray(np.where(finite, image, np.median(image[finite])), dtype=np.float64)
    background = sep.Background(work)
    residual = work - np.array(background)
    objects = sep.extract(residual, SEP_THRESH, err=background.globalrms)
    n_det = 0 if objects is None else len(objects)
    detections = []
    used = np.zeros(n_det, dtype=bool)
    match_pix = max(MATCH_ARCSEC / pixel_scale, 1.5)
    fallback = max(FALLBACK_ARCSEC / pixel_scale, 1.5)
    if n_det:
        kron, flag = sep.kron_radius(
            work, objects['x'], objects['y'], objects['a'], objects['b'], objects['theta'], 6.0)
    else:
        kron = flag = None
    n_ellipse = 0
    n_circle = 0
    for px, py, role in positions:
        chosen = None
        if n_det:
            dist = np.hypot(objects['x'] - px, objects['y'] - py)
            nearest = int(np.argmin(dist))
            if dist[nearest] <= match_pix:
                chosen = nearest
                used[nearest] = True
        if chosen is None:
            _paint(mask, px, py, fallback, fallback, 0.0, circle=True)
            n_circle += 1
            continue
        obj = objects[chosen]
        if flag[chosen] == 0 and np.isfinite(kron[chosen]) and kron[chosen] > 0:
            scale = KRON_FACTOR * float(kron[chosen])
        else:
            scale = 3.0
        _paint(mask, float(obj['x']), float(obj['y']), scale * float(obj['a']), scale * float(obj['b']), float(obj['theta']), circle=False)
        n_ellipse += 1
    for i in range(n_det):
        detections.append({
            'x': float(objects['x'][i]),
            'y': float(objects['y'][i]),
            'a': float(objects['a'][i]),
            'b': float(objects['b'][i]),
            'matched': bool(used[i]),
        })
    info = {
        'n_detections': n_det,
        'n_ellipse': n_ellipse,
        'n_circle': n_circle,
        'n_unmatched_detections': int(n_det - np.sum(used)) if n_det else 0,
    }
    return mask, detections, info


def bunit_refused(header):
    raw = str(header.get('BUNIT', '')).strip()
    value = raw.lower()
    if not value:
        return False, ''
    if 'njy' in value or value in ('jy', 'mag', 'mag/arcsec2', 'mag/arcsec^2'):
        return True, raw
    return False, raw


def load_inputs():
    frames_path = os.path.join(ROOT, 'cutouts', 'frames.tsv')
    members_path = os.path.join(ROOT, 'logs', 'members.tsv')
    selected_path = os.path.join(ROOT, 'logs', 'selected_clusters.tsv')
    if not os.path.isfile(members_path) or not os.path.isfile(frames_path):
        return None
    sheader, selected = read_tsv(selected_path)
    findex = {name: sheader.index(name) for name in sheader}
    clusters = {}
    for row in selected:
        clusters[row[findex['id']]] = {
            'ra': float(row[findex['ra_deg']]),
            'dec': float(row[findex['dec_deg']]),
            'z': float(row[findex['z_spec_bcg']]),
        }
    mheader, member_rows = read_tsv(members_path)
    mindex = {name: mheader.index(name) for name in mheader}
    members = {}
    for row in member_rows:
        members.setdefault(row[mindex['id']], []).append(
            (float(row[mindex['ra_deg']]), float(row[mindex['dec_deg']]))
        )
    fheader, frame_rows = read_tsv(frames_path)
    need = ['id', 'filename', 'sky', 'product']
    if fheader[:4] != need:
        raise SystemExit('frames.tsv header is not id filename sky product')
    frames = []
    for row in frame_rows:
        if len(row) < 4:
            raise SystemExit('frames.tsv row is short')
        sky = row[2].strip()
        product = row[3].strip()
        if sky not in ALLOWED_SKY or product not in ALLOWED_PRODUCT:
            raise SystemExit('frame %s sky=%s product=%s is not an ICL input' % (row[0], sky, product or '(empty)'))
        frames.append({'id': row[0], 'filename': row[1], 'sky': sky, 'product': product})
    return clusters, members, frames


def annular_profile(stacked):
    import numpy as np
    axis = np.asarray(stacked['kpc'], float)
    yy, xx = np.meshgrid(axis, axis, indexing='ij')
    radius = np.hypot(xx, yy)
    sb = np.asarray(stacked['sb'], float)
    weight = np.asarray(stacked['weight'], float)
    edges = np.linspace(0.0, float(np.max(np.abs(axis))), 21)
    lines = ['r_kpc_min\tr_kpc_max\tn_pix\tmean_sb\tweight_sum']
    for i in range(len(edges) - 1):
        if i == len(edges) - 2:
            sel = (radius >= edges[i]) & (radius <= edges[i + 1]) & (weight > 0) & np.isfinite(sb)
        else:
            sel = (radius >= edges[i]) & (radius < edges[i + 1]) & (weight > 0) & np.isfinite(sb)
        mean = float(np.mean(sb[sel])) if np.any(sel) else float('nan')
        lines.append('%.6g\t%.6g\t%d\t%.8g\t%.6g' % (
            edges[i], edges[i + 1], int(np.sum(sel)), mean, float(np.sum(weight[sel])) if np.any(sel) else 0.0))
    return '\n'.join(lines) + '\n'


def run_stack(clusters, members, frames):
    import numpy as np
    sys.path.insert(0, OGF)
    from icl.wgl import mean_stack, warp_metric
    from ogfkit.imageio import save_fits
    from ogfkit.wlcluster import Cosmology
    cosmo = Cosmology(H0=H0, Om=OM)
    warped = []
    detection_lines = ['id\tx\ty\ta\tb\tmatched']
    per_cluster = []
    for frame in frames:
        cluster = clusters.get(frame['id'])
        if cluster is None:
            raise SystemExit('frame %s is not in the selected list' % frame['id'])
        path = frame['filename']
        if not os.path.isabs(path):
            path = os.path.join(ROOT, 'cutouts', path)
        if not os.path.isfile(path):
            raise SystemExit('cutout was not found: %s' % path)
        image, header = first_image(path)
        refused, bunit = bunit_refused(header)
        if refused:
            raise SystemExit('frame %s BUNIT %s is not counts per pixel' % (frame['id'], bunit))
        scale = pixel_scale_arcsec(header)
        sky_positions = [(ra, dec, 'member') for ra, dec in members.get(frame['id'], [])]
        sky_positions.append((cluster['ra'], cluster['dec'], 'bcg'))
        pixels = [world_xy(header, ra, dec) + (role,) for ra, dec, role in sky_positions]
        mask, detections, info = build_mask(image, scale, pixels)
        center = world_xy(header, cluster['ra'], cluster['dec'])
        warped_frame = warp_metric(
            image, center, scale, cluster['z'], cosmo,
            mask=mask, kpc_max=KPC_MAX, n_side=N_SIDE)
        warped_frame['id'] = frame['id']
        warped.append(warped_frame)
        for det in detections:
            detection_lines.append('%s\t%.3f\t%.3f\t%.3f\t%.3f\t%s' % (
                frame['id'], det['x'], det['y'], det['a'], det['b'],
                'yes' if det['matched'] else 'no'))
        per_cluster.append({
            'id': frame['id'],
            'z_spec_bcg': cluster['z'],
            'pixel_scale_arcsec': scale,
            'center_xy': [center[0], center[1]],
            'bunit': bunit,
            'sky': frame['sky'],
            'image_product': frame['product'],
            **info,
        })
    stacked = mean_stack(warped)
    work = os.path.join(ROOT, 'work')
    os.makedirs(work, exist_ok=True)
    extra = {
        'PRODUCT': 'bcg-and-members-masked',
        'SCIENCE': False,
        'COORDSYS': 'kpc',
        'KPCMAX': KPC_MAX,
        'H0': H0,
        'OMEGAM': OM,
    }
    save_fits(os.path.join(work, 'masked_kpc_stack.fits'), stacked['sb'], extra=extra)
    save_fits(os.path.join(work, 'masked_kpc_weight.fits'), stacked['weight'], extra=extra)
    profile = annular_profile(stacked)
    with open(os.path.join(work, 'masked_kpc_profile.tsv'), 'w', encoding='utf-8') as handle:
        handle.write('# product bcg-and-members-masked\n')
        handle.write('# science_measurement false\n')
        handle.write('# coordinate kpc H0 %s Om %s kpc_max %s n_side %s\n' % (H0, OM, KPC_MAX, N_SIDE))
        handle.write(profile)
    with open(os.path.join(ROOT, 'logs', 'mask_detections.tsv'), 'w', encoding='utf-8') as handle:
        handle.write('\n'.join(detection_lines) + '\n')
    center = stacked['weight'].shape[0] // 2
    payload = {
        'product': 'bcg-and-members-masked',
        'product_note': (
            'BCG and CAMIRA member footprints are excluded. '
            'The library product bcg+diffuse is the mean that keeps the BCG. '
            'This file is the masked mean.'),
        'coordinate': 'kpc',
        'kpc_max': KPC_MAX,
        'n_side': N_SIDE,
        'h0': H0,
        'omega_m': OM,
        'mass_used': False,
        'n_clusters': len(warped),
        'ids': stacked['ids'],
        'central_weight': float(stacked['weight'][center, center]),
        'mask_rule': {
            'membership': 'CAMIRA member file rows, no weight cut',
            'bcg': 'cluster catalog coordinate',
            'sep_detection': 'SEP ellipse of 2.5 Kron radii within 1 arcsec',
            'fallback': 'circle of 2 arcsec when SEP has no match',
            'other_detections': 'recorded and left in the image',
            'true_mask_pixel': 'excluded from the mean',
        },
        'clusters': per_cluster,
        'radial_sample_changes': bool(stacked['radial']['sample_changes_with_radius']),
        'science_measurement': False,
        'stack_executed': True,
        'status': 'written',
    }
    with open(os.path.join(work, 'masked_stack.json'), 'w', encoding='utf-8') as handle:
        json.dump(payload, handle, indent=1)
        handle.write('\n')
    nan_center = payload['central_weight']
    if not np.isfinite(nan_center):
        raise SystemExit('central weight is not finite')
    return payload


def main():
    try:
        loaded = load_inputs()
    except SystemExit as exc:
        return write_status({
            'status': 'refused',
            'reason': str(exc),
            'stack_executed': False,
            'product': 'bcg-and-members-masked',
        }, 1)
    if loaded is None:
        return write_status({
            'status': 'stopped',
            'reason': 'cutouts/frames.tsv is absent' if os.path.isfile(os.path.join(ROOT, 'logs', 'members.tsv')) else 'logs/members.tsv is absent',
            'stack_executed': False,
            'product': 'bcg-and-members-masked',
            'coordinate': 'kpc',
            'mass_used': False,
        }, 0)
    clusters, members, frames = loaded
    missing = [frame['id'] for frame in frames if frame['id'] not in clusters]
    if missing:
        return write_status({
            'status': 'refused',
            'reason': 'frames.tsv has ids outside the selected list',
            'ids': missing,
            'stack_executed': False,
            'product': 'bcg-and-members-masked',
        }, 1)
    if len(frames) != len(clusters):
        return write_status({
            'status': 'refused',
            'reason': 'frames.tsv does not list every selected cluster',
            'n_frames': len(frames),
            'n_selected': len(clusters),
            'stack_executed': False,
            'product': 'bcg-and-members-masked',
        }, 1)
    try:
        payload = run_stack(clusters, members, frames)
    except SystemExit as exc:
        return write_status({
            'status': 'refused',
            'reason': str(exc),
            'stack_executed': False,
            'product': 'bcg-and-members-masked',
        }, 1)
    summary = dict(payload)
    summary.pop('clusters', None)
    summary['n_cluster_rows'] = len(payload['clusters'])
    return write_status(summary, 0)


if __name__ == '__main__':
    sys.exit(main())
