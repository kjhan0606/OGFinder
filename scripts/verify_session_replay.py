#!/usr/bin/env python3
"""
Verify that the "Save Session as Python Script" export reproduces a ds9 GUI session.

  verify_session_replay.py [--workdir DIR] [--python /path/to/venv/python] [--display :77]
                           [--skip-hudf] [--skip-m51] [--skip-pipeline] [--keep]

Run it with a python that has numpy + astropy (the OGFinder venv).  Tests:

  T1  HUDF (F160W+F105W+F125W): scripted GUI session under Xvfb (ds9 -source
      scripts/verify_gui_session.tcl; the Tcl calls the procs the menus call) ->
      Save Session as Python Script -> run the script with --mode replay from a plain
      shell (empty HOME, fresh cwd, env -i) -> compare catalog / masks / manifest.
  T2  m51 (single image): same, plus Sersic, morphometry and an ICL chain (Auto Mask,
      masked image, background, BCG centre, profile, measure) -> replay compare.
  T3  determinism: replay m51 with --override '*:--n-workers=1' and '=3' and rerun twice;
      outputs must equal the GUI ones (thread-count independence).
  T4  pipeline mode on NEW data from the HUDF session script: m51 (1 band), two HUDF
      cut-outs (3 bands and 2 bands), a 1-band cut-out; batch (--jobs 3), --resume
      (everything cached; a changed parameter re-runs), summary table.

"Identical" means: catalog TSV text byte-for-byte equal AND every column equal as numbers
with tolerance 0; mask / ICL FITS: pixel arrays bitwise equal (headers compared separately);
manifest: every non-volatile output hash equals the hash recomputed from the file, and equals
the hash of the corresponding GUI file.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
FITS = '/workspace/fits'
results = []          # (test, check, ok, detail)


def rec(test, check, ok, detail=''):
    results.append((test, check, bool(ok), detail))
    print('  [%s] %s: %s%s' % ('PASS' if ok else 'FAIL', test, check, ('  -- ' + detail) if detail else ''), flush=True)
    return ok


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def fits_data_hash(path):
    from astropy.io import fits
    outer = hashlib.sha256()
    with fits.open(path, memmap=False) as hd:
        for h in hd:
            d = h.data
            if d is None:
                continue
            outer.update(hashlib.sha256(np.ascontiguousarray(d).byteswap().tobytes()
                                        if d.dtype.byteorder in ('<', '=') and sys.byteorder == 'little'
                                        else np.ascontiguousarray(d).tobytes()).digest())
    return outer.hexdigest()


def read_tsv(p):
    with open(p, encoding='utf-8', newline='') as f:
        txt = f.read()
    lines = [l for l in txt.split('\n') if l.strip() != '']
    hdr = [h.strip() for h in lines[0].split('\t')]
    rows = [[c.strip() for c in l.split('\t')] for l in lines[1:]]
    return txt, hdr, rows


def compare_catalogs(test, tag, gui, rep):
    t1, h1, r1 = read_tsv(gui)
    t2, h2, r2 = read_tsv(rep)
    ok = rec(test, tag + ': shape', (h1 == h2 and len(r1) == len(r2)),
             'GUI %d rows x %d cols, replay %d rows x %d cols' % (len(r1), len(h1), len(r2), len(h2)))
    if not ok:
        return False
    rec(test, tag + ': text byte-identical', t1 == t2, 'sha256 %s' % hashlib.sha256(t1.encode()).hexdigest()[:16])
    bad = []
    n_num = n_txt = 0
    for j, name in enumerate(h1):
        a = [r[j] if j < len(r) else '' for r in r1]
        b = [r[j] if j < len(r) else '' for r in r2]
        try:
            fa = np.array([float(x) for x in a])
            fb = np.array([float(x) for x in b])
            n_num += 1
            same = np.array_equal(fa, fb, equal_nan=True)
        except ValueError:
            n_txt += 1
            same = (a == b)
        if not same:
            bad.append(name)
    rec(test, tag + ': every column equal (tol 0)', not bad,
        '%d numeric + %d text columns compared%s' % (n_num, n_txt, ('; DIFFER: ' + ','.join(bad)) if bad else ''))
    return not bad


def compare_fits(test, tag, gui, rep):
    from astropy.io import fits
    with fits.open(gui, memmap=False) as a, fits.open(rep, memmap=False) as b:
        da, db = a[0].data, b[0].data
        same = (da is not None and db is not None and da.shape == db.shape and da.dtype == db.dtype
                and np.array_equal(da, db, equal_nan=True))
        hdiff = sorted(k for k in set(a[0].header) | set(b[0].header)
                       if k not in ('', 'COMMENT', 'HISTORY') and a[0].header.get(k) != b[0].header.get(k))
        extra = ''
        if da is not None:
            extra = 'shape %s %s' % (da.shape, da.dtype)
    rec(test, tag + ': pixel data bitwise equal', same, extra + ('; header keys differing: %s' % hdiff if hdiff else '; headers equal'))
    return same


def run(cmd, env=None, cwd=None, log=None, timeout=3600):
    t0 = time.time()
    p = subprocess.run(cmd, env=env, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=timeout)
    out = p.stdout.decode('utf-8', 'replace')
    if log:
        with open(log, 'w') as f:
            f.write(out)
    return p.returncode, out, time.time() - t0


def ensure_xvfb(display):
    num = display.lstrip(':')
    if os.path.exists('/tmp/.X11-unix/X%s' % num):
        return None
    p = subprocess.Popen(['Xvfb', display, '-screen', '0', '1400x1000x24'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    time.sleep(2)
    return p


def run_gui(kind, files, outdir, home, args, extra_env=None):
    """scripted GUI session; returns (rc, seconds, log)"""
    os.makedirs(home, exist_ok=True)
    env = dict(os.environ)
    env.update(HOME=home, DISPLAY=args.display, OGF_VERIFY_KIND=kind, OGF_VERIFY_OUT=outdir,
               OGF_VERIFY_FILES=' '.join(files), OGF_VERIFY_FITS='x', OGFINDER_PYTHON=args.python)
    env.pop('OGF_SESSION_PYTHON', None)
    cmd = [os.path.join(ROOT, 'bin', 'ds9')] + files + ['-geometry', '1300x950', '-source',
                                                       os.path.join(HERE, 'verify_gui_session.tcl')]
    os.makedirs(outdir, exist_ok=True)
    rc, out, dt = run(cmd, env=env, log=os.path.join(outdir, 'gui_stdout.log'), timeout=3000)
    return rc, dt, out


def run_script(script, argv, outdir, home, args, cwd, extra=None):
    """plain-shell run: clean environment, empty HOME, fresh cwd"""
    os.makedirs(home, exist_ok=True)
    os.makedirs(cwd, exist_ok=True)
    env = {'PATH': '/usr/bin:/bin', 'HOME': home, 'LANG': 'C.UTF-8'}
    cmd = [args.python, script, '--python', args.python, '--outdir', outdir] + argv
    return run(cmd, env=env, cwd=cwd, log=outdir + '.log', timeout=3000)


def check_manifest(test, field_dir, gui_files):
    """gui_files: {relative output path or basename -> GUI file path} for hash comparison"""
    mp = os.path.join(field_dir, 'manifest.json')
    if not rec(test, 'manifest written', os.path.exists(mp)):
        return
    man = json.load(open(mp))
    bad = 0
    n = 0
    for o in man['outputs']:
        if o.get('volatile'):
            continue
        p = os.path.join(field_dir, o['path'])
        n += 1
        if not os.path.exists(p) or sha256_file(p) != o['sha256']:
            bad += 1
    rec(test, 'manifest sha256 of %d outputs equals recomputed hash' % n, bad == 0 and n > 0, '%d mismatches' % bad)
    byp = {o['path']: o for o in man['outputs']}
    eq_file = eq_data = tot = 0
    diffs = []
    for rel, gp in gui_files.items():
        o = byp.get(rel)
        if o is None or not os.path.exists(gp):
            diffs.append(rel + ' (missing)')
            continue
        tot += 1
        if o['sha256'] == sha256_file(gp):
            eq_file += 1
        else:
            diffs.append(rel + ' (file hash)')
        if rel.endswith('.fits'):
            if o.get('data_sha256') == fits_data_hash(gp):
                eq_data += 1
            else:
                diffs.append(rel + ' (data hash)')
    nfits = sum(1 for r in gui_files if r.endswith('.fits'))
    rec(test, 'manifest hashes equal the GUI files', not diffs,
        '%d files: %d whole-file sha256 equal, %d/%d FITS data hashes equal%s' % (
            tot, eq_file, eq_data, nfits, ('; DIFFER: ' + ', '.join(diffs)) if diffs else ''))
    t = man.get('tools', {})
    rec(test, 'manifest has tool versions (ds9 git/sextract hash/python pkgs)',
        bool(t.get('ogfinder_git_head')) and 'sha256' in t.get('ds9_sextract', {}) and 'numpy' in t.get('ogfinder_python', {}),
        'git %s, ds9_sextract sha256 %s, numpy %s, astropy %s, sep %s' % (
            t.get('ogfinder_git_head', '')[:10], t.get('ds9_sextract', {}).get('sha256', '')[:10],
            t.get('ogfinder_python', {}).get('numpy'), t.get('ogfinder_python', {}).get('astropy'),
            t.get('ogfinder_python', {}).get('sep')))
    steps = man.get('steps', [])
    rec(test, 'manifest has step timings', bool(steps) and all('seconds' in s for s in steps if s['status'] == 'ok' and s.get('class') != 'config'),
        '%d steps, total %.1f s' % (len(steps), man.get('total_seconds', 0)))
    ins = man.get('inputs', [])
    rec(test, 'manifest input hashes', all(i.get('sha256') == sha256_file(i['path']) for i in ins),
        ', '.join('%s %s' % (i['key'], i['sha256'][:10]) for i in ins))
    return man


def replay_checks(test, kind, gui_out, gui_home, rep_out, field, files):
    field_dir = os.path.join(rep_out, field)
    work = os.path.join(field_dir, 'work')
    outp = os.path.join(field_dir, 'outputs')
    gui_ds9 = os.path.join(gui_home, '.ds9')
    gui_files = {}
    # catalogs saved by the GUI session
    for cat in sorted(f for f in os.listdir(gui_out) if f.startswith(('catalog_gui', 'icl_profile_gui', 'icl_measure_gui'))):
        if not os.path.exists(os.path.join(outp, cat)):
            rec(test, 'catalog ' + cat, False, 'not produced by the replay')
            continue
        compare_catalogs(test, 'catalog ' + cat, os.path.join(gui_out, cat), os.path.join(outp, cat))
        gui_files['outputs/' + cat] = os.path.join(gui_out, cat)
    # data products (masks, ICL) in the GUI's ~/.ds9 vs replay work dir
    prods = sorted(f for f in os.listdir(gui_ds9)
                   if f.endswith('.fits') and f.startswith(('mask_', 'icl_')) and not f.endswith('.tmp'))
    for f in prods:
        rp = os.path.join(work, f)
        if not os.path.exists(rp):
            rec(test, 'product ' + f, False, 'not produced by the replay')
            continue
        compare_fits(test, 'product ' + f, os.path.join(gui_ds9, f), rp)
        gui_files['work/' + f] = os.path.join(gui_ds9, f)
    for f in sorted(os.listdir(gui_ds9)):
        if f.startswith('icl_profile') and f.endswith('.tsv') and os.path.exists(os.path.join(work, f)):
            a = open(os.path.join(gui_ds9, f)).read()
            b = open(os.path.join(work, f)).read()
            rec(test, 'product ' + f, a == b, 'TSV text equal' if a == b else 'TSV text differs')
    man = check_manifest(test, field_dir, gui_files)
    # per-step stdout / catalog checkpoints written by the script itself
    if man:
        st = [s for s in man['steps'] if s['status'] == 'ok']
        sm = [s for s in st if s.get('stdout_matches_gui') is False]
        cm = [s for s in st if s.get('catalog_matches_gui') is False]
        n_so = sum(1 for s in st if s.get('stdout_matches_gui') is True)
        n_cm = sum(1 for s in st if s.get('catalog_matches_gui') is True)
        rec(test, 'script-internal checkpoints (stdout crc32 / catalog crc32 vs GUI)', not sm and not cm,
            '%d tool-stdout and %d catalog checkpoints identical; %d stdout / %d catalog mismatches' % (n_so, n_cm, len(sm), len(cm)))
    return man


def test_replay(name, kind, files, args, base):
    print('\n=== %s: GUI session + replay ===' % name, flush=True)
    gui_out = os.path.join(base, name + '_gui')
    gui_home = os.path.join(base, name + '_guihome')
    shutil.rmtree(gui_out, ignore_errors=True)
    shutil.rmtree(gui_home, ignore_errors=True)
    rc, dt, out = run_gui(kind, files, gui_out, gui_home, args)
    tcl_err = [l for l in out.splitlines() if 'rror' in l and 'VERIFY' not in l and 'filtered' not in l.lower()]
    rec(name, 'GUI session completed (exit 0, VERIFY: DONE)', rc == 0 and 'VERIFY: DONE' in out, '%.1f s' % dt)
    rec(name, 'no Tcl errors on stderr/stdout of the GUI run', not tcl_err, '; '.join(tcl_err[:3]))
    if rc != 0:
        return None
    script = os.path.join(gui_out, 'ogfinder_session.py')
    rec(name, 'session script exported', os.path.exists(script))
    rep_out = os.path.join(base, name + '_replay')
    shutil.rmtree(rep_out, ignore_errors=True)
    argv = ['--mode', 'replay', '--data-dir', FITS] if kind == 'hudf' else ['--mode', 'replay'] + files
    if kind == 'hudf':
        # only the three HUDF files, exactly as in the GUI session
        d = os.path.join(base, name + '_inputs')
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        for f in files:
            os.symlink(f, os.path.join(d, os.path.basename(f)))
        argv = ['--mode', 'replay', '--data-dir', d, '--field-name', 'hudf']
        field = 'hudf'
    else:
        field = 'm51'
    rc, out, dt2 = run_script(script, argv, rep_out, os.path.join(base, name + '_emptyhome'), args,
                              os.path.join(base, name + '_cwd'))
    rec(name, 'replay script exit code 0 (plain shell, empty HOME, fresh cwd)', rc == 0, '%.1f s' % dt2)
    if rc != 0:
        print(out[-3000:])
        return None
    # nothing from the GUI session leaked in: replay must not have touched the emptyhome .ds9
    dd = os.path.join(base, name + '_emptyhome', '.ds9')
    nfiles = sum(len(f) for _, _, f in os.walk(dd)) if os.path.exists(dd) else 0
    rec(name, 'replay created no files in ~/.ds9 (all products go to --outdir)', nfiles == 0, '%d files' % nfiles)
    man = replay_checks(name, kind, gui_out, gui_home, rep_out, field, files)
    return dict(script=script, gui_out=gui_out, gui_home=gui_home, rep_out=rep_out, gui_s=dt, replay_s=dt2, field=field)


def test_determinism(args, base, ctx):
    name = 'T3-determinism'
    print('\n=== %s ===' % name, flush=True)
    script = ctx['script']
    gui_out, gui_home = ctx['gui_out'], ctx['gui_home']
    wk = lambda n: sum((['--override', '%s:--n-workers=%d' % (pat, n)] for pat in ('mask.*', 'analysis.*', 'icl.*')), [])
    for label, extra in (('n-workers=1', wk(1)),
                         ('n-workers=3', wk(3)),
                         ('rerun (default workers)', [])):
        out = os.path.join(base, 'det_' + label.split()[0].replace('=', ''))
        shutil.rmtree(out, ignore_errors=True)
        rc, o, dt = run_script(script, ['--mode', 'replay', FITS + '/m51.fits'] + extra, out,
                               os.path.join(base, 'det_home'), args, os.path.join(base, 'det_cwd'))
        if not rec(name, label + ': script ran', rc == 0, '%.1f s' % dt):
            continue
        fd = os.path.join(out, 'm51')
        for cat in ('catalog_gui.tsv', 'catalog_gui_final.tsv'):
            a = open(os.path.join(gui_out, cat), encoding='utf-8').read()
            b = open(os.path.join(fd, 'outputs', cat), encoding='utf-8').read()
            rec(name, '%s: %s byte-identical to the GUI catalog' % (label, cat), a == b)
        same = True
        for f in ('mask_m51.fits', 'mask_m51_bool.fits', 'icl_background_m51.fits', 'icl_bgsub_m51.fits'):
            gp = os.path.join(gui_home, '.ds9', f)
            from astropy.io import fits
            with fits.open(gp, memmap=False) as A, fits.open(os.path.join(fd, 'work', f), memmap=False) as B:
                same &= np.array_equal(A[0].data, B[0].data, equal_nan=True)
        rec(name, '%s: mask + ICL FITS pixel data bitwise equal to the GUI' % label, same)


def make_cutouts(base, args):
    from astropy.io import fits
    d = os.path.join(base, 'cutouts')
    shutil.rmtree(d, ignore_errors=True)
    os.makedirs(d)
    spec = {'A': (400, 400, 1400), 'B': (1500, 600, 1200), 'C': (2000, 2000, 1000)}   # name: x0,y0,size (0-based)
    out = {}
    for tag, (x0, y0, n) in spec.items():
        for band in ('f160w', 'f105w', 'f125w'):
            with fits.open(os.path.join(FITS, 'hudf_%s.fits' % band), memmap=True) as h:
                data = np.array(h[0].data[y0:y0 + n, x0:x0 + n])
                hd = h[0].header.copy()
            hd['CRPIX1'] = hd['CRPIX1'] - x0
            hd['CRPIX2'] = hd['CRPIX2'] - y0
            p = os.path.join(d, 'cut%s_%s.fits' % (tag, band))
            fits.PrimaryHDU(data.astype(np.float32), hd).writeto(p, overwrite=True)
            out[(tag, band)] = p
    return out


def test_pipeline(args, base, ctx):
    name = 'T4-pipeline'
    print('\n=== %s: pipeline mode on new data ===' % name, flush=True)
    script = ctx['script']          # recorded on HUDF 3 bands (T1)
    cuts = make_cutouts(base, args)
    lst = os.path.join(base, 'fields.txt')
    with open(lst, 'w') as f:
        f.write('# NAME  files (band names come from the FITS FILTER keyword)\n')
        f.write('cutA3 %s %s %s\n' % (cuts[('A', 'f105w')], cuts[('A', 'f125w')], cuts[('A', 'f160w')]))
        f.write('cutB2 %s %s\n' % (cuts[('B', 'f105w')], cuts[('B', 'f160w')]))
        f.write('cutC1 %s\n' % cuts[('C', 'f160w')])
        f.write('m51   %s\n' % (FITS + '/m51.fits'))
    out = os.path.join(base, 'pipe_out')
    shutil.rmtree(out, ignore_errors=True)
    rc, o, dt = run_script(script, ['--input-list', lst, '--jobs', '3'], out, os.path.join(base, 'pipe_home'), args,
                           os.path.join(base, 'pipe_cwd'))
    rec(name, 'batch run (4 fields, --jobs 3) exit code 0', rc == 0, '%.1f s' % dt)
    if rc != 0:
        print(o[-3000:])
    summ = os.path.join(out, 'summary.tsv')
    rec(name, 'summary.tsv written', os.path.exists(summ))
    rows = {}
    if os.path.exists(summ):
        L = [l.split('\t') for l in open(summ).read().splitlines()]
        hdr = L[0]
        rows = {r[0]: dict(zip(hdr, r)) for r in L[1:]}
        print('  summary.tsv:')
        for l in open(summ).read().splitlines():
            print('    ' + l.replace('\t', ' | '))
    for fld, nb in (('cutA3', 3), ('cutB2', 2), ('cutC1', 1), ('m51', 1)):
        r = rows.get(fld, {})
        ok = r.get('status') == 'ok' and int(r.get('objects', 0) or 0) > 0
        mf = r.get('masked_fraction', '')
        rec(name, '%s (%d band%s): completed, objects > 0, mask made' % (fld, nb, 's' if nb > 1 else ''), ok and mf != '',
            '%s objects, masked fraction %s, %s s' % (r.get('objects'), mf, r.get('seconds')))
        mpath = os.path.join(out, fld, 'manifest.json')
        if os.path.exists(mpath):
            man = json.load(open(mpath))
            sk = [(s['step'], s.get('reason', '')) for s in man['steps'] if s['status'] == 'skipped']
            fl = [s for s in man['steps'] if s['status'] == 'failed']
            cols = r.get('columns')
            exp_mb = nb >= 2
            have_mb = any(s['step'] == 'bands.measure' and s['status'] == 'ok' for s in man['steps'])
            rec(name, '%s: band-dependent steps %s' % (fld, 'ran' if exp_mb else 'skipped with a message'),
                have_mb == exp_mb and not fl,
                'columns=%s; skipped: %s' % (cols, '; '.join('%s [%s]' % (a, b[:60]) for a, b in sk if a in ('bands.measure', 'mask.copy_to_bands'))))
            if nb >= 2:
                cat = os.path.join(out, fld, 'catalog_final.tsv')
                hdr = open(cat).readline().rstrip('\n').split('\t')
                bands = sorted(b for b in ('F105W', 'F125W', 'F160W') if 'MAG_' + b in hdr)
                rec(name, '%s: MAG_<band> columns for the %d supplied bands' % (fld, nb), len(bands) == nb, ','.join(bands))
    # sensible outputs: magnitudes in a plausible range for HUDF cut-outs
    cat = os.path.join(out, 'cutA3', 'catalog_final.tsv')
    if os.path.exists(cat):
        _, h, r = read_tsv(cat)
        i = h.index('MAG_AUTO')
        mags = np.array([float(x[i]) for x in r])
        mags = mags[mags < 90]
        rec(name, 'cutA3: MAG_AUTO range sensible (header ZP used)', 18 < np.median(mags) < 30 and len(mags) > 100,
            'median %.2f, 5-95%% %.2f..%.2f, N=%d' % (np.median(mags), *np.percentile(mags, [5, 95]), len(mags)))
    # --resume: everything cached
    rc, o, dt2 = run_script(script, ['--input-list', lst, '--jobs', '3', '--resume'], out, os.path.join(base, 'pipe_home'), args,
                            os.path.join(base, 'pipe_cwd'))
    ncached = sum(int(json.load(open(os.path.join(out, f, 'manifest.json')))['step_counts']['cached']) for f in rows)
    nran = sum(int(json.load(open(os.path.join(out, f, 'manifest.json')))['step_counts']['ran']) for f in rows)
    rec(name, '--resume rerun skips completed steps', rc == 0 and nran == 0 and ncached > 0,
        '%d steps cached, %d re-run, %.1f s (first run %.1f s)' % (ncached, nran, dt2, dt))
    # changed parameter -> re-run
    rc, o, dt3 = run_script(script, ['--input-list', lst, '--jobs', '3', '--resume', '--override',
                                    'mask.auto:--detect-thresh=4.0'], out, os.path.join(base, 'pipe_home'), args,
                            os.path.join(base, 'pipe_cwd'))
    m = json.load(open(os.path.join(out, 'm51', 'manifest.json')))
    ran = [s['step'] for s in m['steps'] if s['status'] == 'ok']
    cached = [s['step'] for s in m['steps'] if s['status'] == 'cached']
    rec(name, '--resume with a changed parameter re-runs the changed step (and what follows)', rc == 0 and 'mask.auto' in ran and not [s for s in m['steps'] if s['status'] == 'failed'],
        '(design: a changed step restarts the field from the beginning because steps edit masks in place) m51: cached=%s ran=%s (%.1f s)' % (','.join(cached), ','.join(ran), dt3))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--workdir', default='/workspace/ogf_verify')
    ap.add_argument('--python', default=sys.executable)
    ap.add_argument('--display', default=':77')
    ap.add_argument('--skip-hudf', action='store_true')
    ap.add_argument('--skip-m51', action='store_true')
    ap.add_argument('--skip-pipeline', action='store_true')
    a = ap.parse_args()
    a.python = os.path.abspath(a.python) if os.sep in a.python else a.python
    base = os.path.abspath(a.workdir)
    os.makedirs(base, exist_ok=True)
    xv = ensure_xvfb(a.display)
    t0 = time.time()
    hud = m51 = None
    timings = {}
    try:
        if not a.skip_hudf:
            hud = test_replay('T1-hudf', 'hudf', [FITS + '/hudf_f160w.fits', FITS + '/hudf_f105w.fits', FITS + '/hudf_f125w.fits'], a, base)
            if hud:
                timings['T1 GUI session'] = hud['gui_s']; timings['T1 replay'] = hud['replay_s']
        if not a.skip_m51:
            m51 = test_replay('T2-m51', 'm51', [FITS + '/m51.fits'], a, base)
            if m51:
                timings['T2 GUI session'] = m51['gui_s']; timings['T2 replay'] = m51['replay_s']
                test_determinism(a, base, m51)
        if not a.skip_pipeline and hud:
            t1 = time.time()
            test_pipeline(a, base, hud)
            timings['T4 pipeline (batch+resume)'] = time.time() - t1
    finally:
        if xv:
            xv.terminate()
    print('\n' + '=' * 78)
    print('SUMMARY')
    nfail = 0
    for t in sorted(set(r[0] for r in results)):
        rs = [r for r in results if r[0] == t]
        f = [r for r in rs if not r[2]]
        nfail += len(f)
        print('  %-16s %3d checks, %d failed  %s' % (t, len(rs), len(f), 'PASS' if not f else 'FAIL'))
        for r in f:
            print('      FAILED: %s -- %s' % (r[1], r[3]))
    print('  timings (s): ' + ', '.join('%s %.1f' % kv for kv in timings.items()))
    print('  total %.1f s;  %s' % (time.time() - t0, 'ALL CHECKS PASSED' if nfail == 0 else '%d CHECK(S) FAILED' % nfail))
    return 1 if nfail else 0


if __name__ == '__main__':
    sys.exit(main())
