"""Cluster / lensing tools against synthetic truth (and a real HUDF crop when the FITS files exist).  Numbers are printed with -s."""
def _import_ogfmeas():
    """In-tree measurements. Finds ogfmeas from this file so a script does not need PYTHONPATH."""
    import pathlib
    import sys
    for parent in pathlib.Path(__file__).resolve().parents:
        if (parent / "ogfmeas" / "__init__.py").is_file():
            folder = str(parent)
            if folder not in sys.path:
                sys.path.insert(0, folder)
            break
    import ogfmeas
    return ogfmeas.measurement_library()

import json
import math
import os
import subprocess
import sys

import numpy as np
import pytest

from ogfkit import cluster as cl, clustersynth as cs, imageio, tsvio

HERE = os.path.dirname(os.path.abspath(__file__))
PLUGIN = os.path.dirname(HERE)
CLI = os.path.join(PLUGIN, 'cluster.py')
FITS = os.environ.get('OGF_TEST_FITS', '/workspace/fits')


def say(*a):
    sys.stderr.write(' '.join(str(x) for x in a) + '\n')


def write_cat(path, cols, data):
    n = len(next(iter(data.values())))
    with open(path, 'w') as f:
        f.write('\t'.join(['NUMBER'] + cols) + '\n')
        for i in range(n):
            f.write('\t'.join([str(i + 1)] + ['%.5g' % data[c][i] if not isinstance(data[c][i], str) else data[c][i] for c in cols]) + '\n')


def run_cli(args):
    p = subprocess.run([sys.executable, CLI] + args, capture_output=True, text=True)
    assert p.returncode == 0, p.stderr[-800:]
    lines = [l for l in p.stdout.splitlines() if l]
    cols = lines[0].split('\t')
    rows = [dict(zip(cols, l.split('\t'))) for l in lines[1:]]
    return cols, rows, p.stderr


def fnum(rows, c):
    return np.array([tsvio.fnum(r.get(c)) for r in rows])


def test_red_sequence_recovery_many_seeds():
    A = []
    for seed in range(1, 41):
        d = cs.cluster_catalog(seed)
        r = np.hypot(d['x'] - d['centre'][0], d['y'] - d['centre'][1])
        s = r < 400
        rs = cl.red_sequence(d['color'][s], d['mag'][s], d['color_err'][s], mag_range=(19, 24), m0=21)
        assert rs['ok']
        A.append([rs['a'], rs['b'], rs['scatter'], rs['a_err'], rs['b_err']])
    A = np.array(A)
    t = np.array([1.2, -0.04, 0.05])
    bias = A[:, :3].mean(0) - t
    rms = np.sqrt(((A[:, :3] - t) ** 2).mean(0))
    pa, pb = ((A[:, 0] - 1.2) / A[:, 3]).std(), ((A[:, 1] + 0.04) / A[:, 4]).std()
    say('RS recovery (40 clusters of 150 members + field): mean colour(21) %.4f (true 1.2), slope %.4f (-0.04), scatter %.4f (0.05); rms errors %.4f %.4f %.4f; pull std %.2f %.2f'
        % (A[:, 0].mean(), A[:, 1].mean(), A[:, 2].mean(), rms[0], rms[1], rms[2], pa, pb))
    assert abs(bias[0]) < 0.01 and abs(bias[1]) < 0.006 and abs(bias[2]) < 0.012
    assert pa < 1.6 and pb < 1.6                                   # quoted errors are slightly optimistic (see docs)


def test_fixed_slope_and_mixture_beats_nothing():
    d = cs.cluster_catalog(5)
    r = np.hypot(d['x'] - d['centre'][0], d['y'] - d['centre'][1])
    s = r < 400
    rs0 = cl.red_sequence(d['color'][s], d['mag'][s], d['color_err'][s], mag_range=(19, 24), m0=21, fit_slope=False)
    assert rs0['ok'] and rs0['b'] == 0.0 and rs0['b_err'] == 0.0
    assert abs(rs0['a'] - (1.2 - 0.04 * (np.median(d['mag'][s & (d['mag'] > 19) & (d['mag'] < 24)]) - 21)) - 0.0) < 0.2


def test_rs_too_few_objects():
    rs = cl.red_sequence([1.0] * 3, [20.0, 21.0, 22.0])
    assert not rs['ok'] and rs['note']


def test_members_cli_purity_completeness(tmp_path):
    d = cs.cluster_catalog(3)
    n = len(d['x'])
    cat = tmp_path / 'c.tsv'
    write_cat(str(cat), ['X_IMAGE', 'Y_IMAGE', 'MAG_R', 'MAG_G', 'MAGERR_R', 'MAGERR_G', 'Z_SPEC'],
              dict(X_IMAGE=d['x'] + 1, Y_IMAGE=d['y'] + 1, MAG_R=d['mag'], MAG_G=d['color'] + d['mag'], MAGERR_R=d['color_err'] * 0.7, MAGERR_G=d['color_err'] * 0.7, Z_SPEC=d['z']))
    base = ['--task', 'members', '--catalog', str(cat), '--work', str(tmp_path / 'w'), '--meta-out', str(tmp_path / 'meta.json'), '--mag-column', 'MAG_R', '--blue-column', 'MAG_G',
            '--red-column', 'MAG_R', '--blue-err-column', 'MAGERR_G', '--red-err-column', 'MAGERR_R', '--mag-min', '19.5', '--mag-max', '24']
    for label, extra in (('colour + radial', []), ('colour + radial + redshift', ['--z-column', 'Z_SPEC'])):
        cols, rows, err = run_cli(base + extra)
        assert cols == ['NUMBER', 'CL_RS_RES', 'CL_RS_PULL', 'CL_RS', 'CL_R', 'CL_PMEM', 'CL_MEMBER']
        mem = fnum(rows, 'CL_MEMBER') > 0
        truth = d['member'] & (d['mag'] >= 19.5) & (d['mag'] <= 24)
        tp = (mem & truth).sum()
        pur, comp = tp / max(mem.sum(), 1), tp / truth.sum()
        pm = fnum(rows, 'CL_PMEM')
        rsflag = (fnum(rows, 'CL_RS') > 0)
        say('members (%s, p_mem >= 0.5): %d selected, purity %.3f, completeness %.3f (true members %d in the magnitude range); all %d red-sequence objects: sum p_mem %.1f vs %d true members among them'
            % (label, mem.sum(), pur, comp, truth.sum(), rsflag.sum(), np.nansum(pm[rsflag]), (rsflag & d['member']).sum()))
        assert pur > 0.75 and comp > 0.65
        # the radial probabilities are calibrated: expected true members among the colour-selected objects agree with the truth to ~ 20 %
        tr = (rsflag & d['member']).sum()
        assert abs(np.nansum(pm[rsflag]) - tr) < 0.15 * tr + 3
    meta = json.load(open(tmp_path / 'meta.json'))['cluster']
    cx, cy = meta['centre']
    say('auto centre (%.0f, %.0f) vs true (%.0f, %.0f): %.1f px; RS a %.3f b %.4f scatter %.3f' % (cx, cy, d['centre'][0], d['centre'][1], math.hypot(cx - d['centre'][0], cy - d['centre'][1]), meta['rs_a'], meta['rs_b'], meta['rs_scatter']))
    assert math.hypot(cx - d['centre'][0], cy - d['centre'][1]) < 40
    assert os.path.exists(tmp_path / 'w' / 'cluster_rs.png')


def test_density_map_is_calibrated_poisson():
    shape = (2000, 2000)
    rng = np.random.default_rng(1)
    stds, means, mx = [], [], []
    for k in range(12):
        n = 2500
        x, y = rng.uniform(0, 2000, n), rng.uniform(0, 2000, n)
        dm = cl.density_map(x, y, shape, 60.0, binsize=5)
        g = dm['good']
        stds.append(dm['sig'][g].std()); means.append(dm['sig'][g].mean()); mx.append(dm['sig'][g].max())
    say('density significance, 12 uniform Poisson fields (2500 objects, sigma 60 px): map mean %.3f, std %.3f (expect 0, 1), max %.2f' % (np.mean(means), np.mean(stds), np.max(mx)))
    assert abs(np.mean(stds) - 1.0) < 0.12 and abs(np.mean(means)) < 0.1 and max(mx) < 5.0


def test_density_peak_recovery_and_masked_area():
    d = cs.cluster_catalog(2, shape=(3000, 3000), centre=(1500.0, 1200.0), n_field=2500, n_mem=250)
    valid = np.ones((3000, 3000), bool)
    valid[:, 2400:] = False                                  # a chip gap / edge
    ok = d['x'] < 2400
    dm = cl.density_map(d['x'][ok], d['y'][ok], (3000, 3000), 90.0, binsize=6, valid=valid)
    pk = cl.find_peaks(dm, 4.0)
    err = math.hypot(pk[0]['x'] - 1500, pk[0]['y'] - 1200)
    edge = dm['sig'][:, 2400 // 6 + 2:].max() if dm['sig'][:, 2400 // 6 + 2:].size else 0
    say('density peak at (%.0f, %.0f), %.1f px from the true centre, %.1f sigma; %d peaks > 4 sigma; masked columns have sig %.2f' % (pk[0]['x'], pk[0]['y'], err, pk[0]['sig'], len(pk), edge))
    assert err < 30 and pk[0]['sig'] > 8
    assert edge == 0.0
    # background rate recovered from the clipped mean (true: field only, 2500 over 3000^2 minus the masked part -> per 6 px cell)
    expect = 2500 * 36 / 9e6
    say('background per cell %.5f vs field-only expectation %.5f' % (dm['lam'], expect))
    assert abs(dm['lam'] / expect - 1) < 0.15


def _arc_catalog(seed, tmp, **kw):
    sep = _import_ogfmeas()
    from scipy import ndimage as ndi
    img, truth = cs.arc_field(seed, **kw)
    bk = sep.Background(np.ascontiguousarray(img))
    sub = img - bk.back()
    objs, seg = sep.extract(sub, 1.5, err=bk.globalrms, minarea=8, deblend_cont=1.0, segmentation_map=True)
    fits_path = os.path.join(tmp, 'arc%d.fits' % seed)
    imageio.save_fits(fits_path, sub.astype(np.float32))
    cat_path = os.path.join(tmp, 'arc%d.tsv' % seed)
    write_cat(cat_path, ['X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE'], dict(X_IMAGE=objs['x'] + 1, Y_IMAGE=objs['y'] + 1, A_IMAGE=objs['a'], B_IMAGE=objs['b']))
    # truth object -> detection index via the segmentation map along the centre line
    match = []
    for t in truth:
        rx, ry = cs.truth_ridge(t, (350.0, 350.0))
        lab = [seg[int(round(y)), int(round(x))] for x, y in zip(rx, ry) if 0 <= x < seg.shape[1] and 0 <= y < seg.shape[0]]
        lab = [l for l in lab if l > 0]
        match.append(-1 if not lab else int(np.bincount(lab).argmax()) - 1)
    return fits_path, cat_path, truth, match


def _arc_stats(seeds, tmp, **kw):
    rec = dict(arc=[0, 0], other=[0, 0]); errL, errR, errA, errS = [], [], [], []
    fp_kinds = {}
    for seed in seeds:
        f, c, truth, match = _arc_catalog(seed, tmp, **kw)
        cols, rows, _ = run_cli(['--task', 'arcs', '--catalog', c, '--image', f, '--work', os.path.join(tmp, 'w%d' % seed), '--center-mode', 'manual', '--center-x', '351', '--center-y', '351'])
        flag = fnum(rows, 'ARC_FLAG'); L = fnum(rows, 'ARC_LEN'); R = fnum(rows, 'ARC_RCURV'); al = fnum(rows, 'ARC_ALIGN'); sp = fnum(rows, 'ARC_SPAN')
        for t, j in zip(truth, match):
            hit = j >= 0 and flag[j] == 1
            key = 'arc' if t['kind'] == 'arc' else 'other'
            rec[key][1] += 1
            rec[key][0] += int(hit)
            if t['kind'] != 'arc' and hit:
                fp_kinds[t['kind']] = fp_kinds.get(t['kind'], 0) + 1
            if t['kind'] == 'arc' and j >= 0 and np.isfinite(L[j]) and np.isfinite(R[j]):
                errL.append(L[j] / t['length'] - 1); errR.append(R[j] / t['R'] - 1); errA.append(al[j]); errS.append(sp[j] - t['span_deg'])
    return rec, fp_kinds, np.array(errL), np.array(errR), np.array(errA), np.array(errS)


def test_arcs_bright(tmp_path):
    rec, fp, eL, eR, eA, eS = _arc_stats(range(1, 7), str(tmp_path))
    say('arcs, amp 3-8 sigma (6 fields, %d truth arcs): recovered as flagged candidates %d/%d; false candidates %d/%d (%s); length error median %+.3f (rms %.3f), curvature radius error median %+.3f (rms %.3f), '
        'alignment max %.1f deg, span error median %+.1f deg' % (rec['arc'][1], rec['arc'][0], rec['arc'][1], rec['other'][0], rec['other'][1], fp, np.median(eL), eL.std(), np.median(eR), eR.std(), eA.max(), np.median(eS)))
    assert rec['arc'][0] / rec['arc'][1] > 0.8
    assert rec['other'][0] == 0
    assert abs(np.median(eL)) < 0.04 and abs(np.median(eR)) < 0.04


def test_arcs_faint(tmp_path):
    rec, fp, eL, eR, eA, eS = _arc_stats(range(11, 17), str(tmp_path), amp_range=(1.2, 2.5))
    say('arcs, amp 1.2-2.5 sigma (6 fields): recovered %d/%d; false candidates %d/%d (%s); length error median %+.3f (rms %.3f), R error median %+.3f (rms %.3f)'
        % (rec['arc'][0], rec['arc'][1], rec['other'][0], rec['other'][1], fp, np.median(eL) if len(eL) else float('nan'), eL.std() if len(eL) else float('nan'), np.median(eR) if len(eR) else float('nan'), eR.std() if len(eR) else float('nan')))
    assert rec['other'][0] <= 1


def test_arc_alignment_discriminates_without_centre_flag(tmp_path):
    """Radial lines are rejected by the alignment cut and straight tangential lines by the curvature cut."""
    f, c, truth, match = _arc_catalog(1, str(tmp_path))
    args = ['--task', 'arcs', '--catalog', c, '--image', f, '--work', str(tmp_path / 'w'), '--center-mode', 'manual', '--center-x', '351', '--center-y', '351']
    _, rows, _ = run_cli(args)
    flag = fnum(rows, 'ARC_FLAG'); al = fnum(rows, 'ARC_ALIGN')
    rad = [al[j] for t, j in zip(truth, match) if t['kind'] == 'radial' and j >= 0]
    st = [flag[j] for t, j in zip(truth, match) if t['kind'] == 'straight_tangential' and j >= 0]
    say('radial lines: alignment angles %s (flagged %d); tangential straight lines flagged %s with require-curved, %s without' % (np.round(rad, 1), sum(flag[j] for t, j in zip(truth, match) if t['kind'] == 'radial' and j >= 0), st, ''))
    assert min(rad) > 80 and sum(st) == 0
    _, rows2, _ = run_cli(args + ['--arc-allow-straight'])
    flag2 = fnum(rows2, 'ARC_FLAG')
    st2 = [flag2[j] for t, j in zip(truth, match) if t['kind'] == 'straight_tangential' and j >= 0]
    say('straight tangential lines flagged without the curvature requirement: %d/%d' % (sum(st2), len(st2)))
    assert sum(st2) >= len(st2) - 1


def test_cli_deterministic_and_density_columns(tmp_path):
    d = cs.cluster_catalog(4, n_field=800, n_mem=80)
    cat = tmp_path / 'c.tsv'
    write_cat(str(cat), ['X_IMAGE', 'Y_IMAGE', 'MAG_R'], dict(X_IMAGE=d['x'] + 1, Y_IMAGE=d['y'] + 1, MAG_R=d['mag']))
    args = ['--task', 'density', '--catalog', str(cat), '--work', str(tmp_path / 'w'), '--mag-column', 'MAG_R', '--dens-sigma', '100']
    c1, r1, _ = run_cli(args)
    c2, r2, _ = run_cli(args)
    assert c1 == ['NUMBER', 'CL_SIGMA'] and r1 == r2
    assert os.path.exists(tmp_path / 'w' / 'cluster_density.fits') and os.path.exists(tmp_path / 'w' / 'cluster_peaks.tsv')
    s = fnum(r1, 'CL_SIGMA')
    say('density CLI: members have mean CL_SIGMA %.1f vs field %.1f' % (np.nanmean(s[d['member']]), np.nanmean(s[~d['member']])))
    assert np.nanmean(s[d['member']]) > 3 * max(np.nanmean(s[~d['member']]), 0.5)


def test_missing_columns_reported(tmp_path):
    cat = tmp_path / 'c.tsv'
    cat.write_text('NUMBER\tX_IMAGE\tY_IMAGE\n1\t1\t1\n')
    p = subprocess.run([sys.executable, CLI, '--task', 'members', '--catalog', str(cat), '--work', str(tmp_path)], capture_output=True, text=True)
    assert p.returncode != 0 and 'missing column' in (p.stderr + p.stdout)


# ---------------------------------------------------------------- real HUDF crop (three WFC3 bands)

HUDF = [os.path.join(FITS, 'hudf_f105w.fits'), os.path.join(FITS, 'hudf_f125w.fits'), os.path.join(FITS, 'hudf_f160w.fits')]


@pytest.mark.skipif(not all(os.path.exists(p) for p in HUDF), reason='HUDF FITS files not available')
def test_hudf_multiband_null_and_injection(tmp_path):
    sep = _import_ogfmeas()
    from scipy import ndimage as ndi
    sl = (slice(1000, 2800), slice(1000, 2800))
    imgs = [imageio.load_image(p)[0][sl].astype(np.float64) for p in HUDF]
    det = imgs[2]
    bk = sep.Background(np.ascontiguousarray(det), bw=64, bh=64)
    sub = det - bk.back()
    objs, seg0 = sep.extract(sub, 2.0, err=bk.globalrms, minarea=10, segmentation_map=True)
    occupied = ndi.binary_dilation(seg0 > 0, iterations=4)
    zp = 25.94
    mags = []
    for im, zpb in ((imgs[0], 26.27), (imgs[2], 25.94)):
        b = sep.Background(np.ascontiguousarray(im), bw=64, bh=64)
        fl, fe, _ = sep.sum_circle(im - b.back(), objs['x'], objs['y'], 4.0, err=b.globalrms)
        mags.append((zpb - 2.5 * np.log10(np.clip(fl, 1e-9, None)), np.abs(1.0857 * fe / np.clip(fl, 1e-9, None)), fl / np.clip(fe, 1e-9, None)))
    snr = np.minimum(mags[0][2], mags[1][2])
    ok = snr > 8
    n = len(objs)
    cat = tmp_path / 'hudf.tsv'
    write_cat(str(cat), ['X_IMAGE', 'Y_IMAGE', 'MAG_160', 'MAG_105', 'E160', 'E105', 'A_IMAGE', 'B_IMAGE'],
              dict(X_IMAGE=objs['x'][ok] + 1, Y_IMAGE=objs['y'][ok] + 1, MAG_160=mags[1][0][ok], MAG_105=mags[0][0][ok], E160=mags[1][1][ok], E105=mags[0][1][ok], A_IMAGE=objs['a'][ok], B_IMAGE=objs['b'][ok]))
    nobj = int(ok.sum())
    crop = str(tmp_path / 'crop160.fits')
    imageio.save_fits(crop, det.astype(np.float32))
    # null test: a blank-ish field has no cluster -> the density map should be consistent with noise-like fluctuations of the galaxy clustering, no huge peak
    cols, rows, err = run_cli(['--task', 'density', '--catalog', str(cat), '--work', str(tmp_path / 'wd'), '--mag-column', 'MAG_160', '--dens-sigma', '60', '--dens-bin', '6', '--image', crop])
    s = fnum(rows, 'CL_SIGMA')
    summ = json.load(open(tmp_path / 'wd' / 'cluster_summary.json'))['density']
    say('HUDF crop 1800x1800 px, %d objects (S/N>8 in F105W and F160W): density map max %.1f sigma, %d peaks >3 sigma (a Poisson field of this size would give ~0-1; galaxy clustering raises it)' % (nobj, summ['max_sigma'], len(summ['peaks'])))
    assert summ['max_sigma'] < 8.0
    # colour-magnitude diagram of the real field: the fit must run and report the sequence it finds without crashing; members are NOT expected
    cols, rows, err = run_cli(['--task', 'members', '--catalog', str(cat), '--work', str(tmp_path / 'wm'), '--mag-column', 'MAG_160', '--blue-column', 'MAG_105', '--red-column', 'MAG_160',
                               '--blue-err-column', 'E105', '--red-err-column', 'E160', '--mag-min', '22', '--mag-max', '27', '--fit-radius', '700'])
    sm = json.load(open(tmp_path / 'wm' / 'cluster_summary.json'))
    say('HUDF radial model:', {k: (round(v, 6) if isinstance(v, float) else v) for k, v in sm['radial_profile'].items()}, 'members', sm.get('n_members'), 'centre', sm['centre'], sm['centre_how'])
    say('HUDF CMD (F105W-F160W vs F160W): fit %s; colour(m0) %.3f slope %.3f scatter %.3f; %d objects within the sequence, %d with radial p_mem>0 (field: expected contamination %.1f)'
        % (sm['red_sequence']['ok'], sm['red_sequence']['a'], sm['red_sequence']['b'], sm['red_sequence']['scatter'], sm.get('n_rs', 0), int(np.sum(fnum(rows, 'CL_PMEM') > 0)), sm.get('expected_contamination', float('nan'))))
    # injection: tangential arcs added to the real F160W crop, measured on the real noise and real neighbours
    rng = np.random.default_rng(7)
    inj = np.zeros_like(det)
    centre = (900.0, 900.0)
    truths = []
    sd = imageio.robust_sigma(sub)
    def ridge_pts(R, phi, span):
        return [(centre[0] + R * math.cos(phi + f * span), centre[1] + R * math.sin(phi + f * span)) for f in np.linspace(-0.5, 0.5, 60)]
    ridges = []
    for _try in range(20000):
        if len(truths) >= 10:
            break
        R = rng.uniform(120, 420); phi = rng.uniform(0, 2 * math.pi)
        span = math.radians(rng.uniform(35, 70))
        pts = ridge_pts(R, phi, span)
        if not all(40 < x < 1760 and 40 < y < 1760 for x, y in pts):
            continue
        if sum(occupied[int(round(y)), int(round(x))] for x, y in pts) > 2:         # inject on empty sky (completeness-test protocol)
            continue
        if any(np.min(np.hypot(np.array([q[0] for q in pts])[:, None] - np.array([q[0] for q in o])[None, :], np.array([q[1] for q in pts])[:, None] - np.array([q[1] for q in o])[None, :])) < 40 for o in ridges):
            continue
        ridges.append(pts)
        mx, my = centre[0] + R * math.cos(phi), centre[1] + R * math.sin(phi)
        cs.add_arc(inj, centre[0], centre[1], R, phi, span, 5.0, 8.0 * sd * 1.0)
        truths.append((mx, my, R, math.degrees(span), R * span))
    injected = sub + __import__('scipy.ndimage', fromlist=['x']).gaussian_filter(inj, 1.5)
    path = str(tmp_path / 'inj.fits')
    imageio.save_fits(path, injected.astype(np.float32))
    b2 = sep.Background(np.ascontiguousarray(injected), bw=64, bh=64)
    o2, seg = sep.extract(injected - b2.back(), 3.0, err=b2.globalrms, minarea=15, deblend_cont=1.0, segmentation_map=True)
    c2 = tmp_path / 'inj.tsv'
    write_cat(str(c2), ['X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE'], dict(X_IMAGE=o2['x'] + 1, Y_IMAGE=o2['y'] + 1, A_IMAGE=o2['a'], B_IMAGE=o2['b']))
    cols, rows, err = run_cli(['--task', 'arcs', '--catalog', str(c2), '--image', path, '--work', str(tmp_path / 'wa'), '--center-mode', 'manual', '--center-x', '901', '--center-y', '901'])
    flag = fnum(rows, 'ARC_FLAG'); L = fnum(rows, 'ARC_LEN'); LW = fnum(rows, 'ARC_LW'); SN = fnum(rows, 'ARC_SNR'); CUR = fnum(rows, 'ARC_RCURV'); AL = fnum(rows, 'ARC_ALIGN')
    found, errs = 0, []
    for (mx, my, R, spd, length) in truths:
        ph = math.atan2(my - centre[1], mx - centre[0])
        pts = [(centre[0] + R * math.cos(ph + f * math.radians(spd)), centre[1] + R * math.sin(ph + f * math.radians(spd))) for f in np.linspace(-0.5, 0.5, 40)]
        lab = [seg[int(round(y)), int(round(x))] for x, y in pts]
        lab = [l for l in lab if l > 0]
        if lab:
            j = int(np.bincount(lab).argmax()) - 1
            if flag[j] == 1:
                found += 1; errs.append(L[j] / length - 1)
    nfp = int(np.nansum(flag)) - found
    say('HUDF + 10 injected arcs (amp 8 sigma): flagged %d/10; length error median %+.3f; other flagged objects among %d detections: %d (real galaxies/noise)' % (found, np.median(errs) if errs else float('nan'), len(o2), nfp))
    assert found >= 5 and nfp <= 3
