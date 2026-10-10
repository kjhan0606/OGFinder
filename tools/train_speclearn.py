#!/usr/bin/env python3
"""Train the in-tree spectrum models on local SDSS stars and DESI galaxies.

Writes npz files under --out (default /tmp/speclearn_train).  Does not modify
the survey archive and does not write a supernova template library.

The star rate is how often a held-out SDSS spectrum receives the temperature
letter of its subclass.  The galaxy rate is how often the neighbour vote
agrees with galaxy_analysis on spectra left out of both the autoencoder and
the neighbour library.  That galaxy rate repeats this line diagnostic.  It is
not a FastSpecFit or SDSS subclass accuracy.
"""
import argparse
import hashlib
import json
import os
import sys
import time

os.environ.setdefault('OMP_NUM_THREADS', '1')
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('MKL_NUM_THREADS', '1')
os.environ.setdefault('NUMEXPR_NUM_THREADS', '1')

import numpy as np

ROOT = '/gpfs/kjhan/Astro_Foundation'
SDSS_MANIFEST = os.path.join(ROOT, 'data/raw/spectra/sdss_dr17/manifest.jsonl')
DESI_DIR = os.path.join(ROOT, 'data/raw/spectra/desi_dr1_sparcl')


def _holdout(key):
    digest = hashlib.sha1(str(key).encode('utf-8')).hexdigest()
    return int(digest[:8], 16) % 5 == 0


def select_stars(cap, seed):
    from ogfkit.speclearn import temperature_letter, LETTERS
    rng = np.random.default_rng(seed)
    buckets = {L: [] for L in LETTERS}
    seen = {L: 0 for L in LETTERS}
    if not os.path.isfile(SDSS_MANIFEST):
        raise FileNotFoundError('SDSS manifest was not found: %s' % SDSS_MANIFEST)
    with open(SDSS_MANIFEST) as fh:
        for line in fh:
            rec = json.loads(line)
            cat = rec.get('catalog') or {}
            if int(cat.get('zWarning') or 0) != 0:
                continue
            if (cat.get('class') or '').upper() != 'STAR':
                continue
            try:
                z = float(cat['z'])
            except (TypeError, ValueError, KeyError):
                continue
            if not np.isfinite(z) or abs(z) >= 0.002:
                continue
            letter = temperature_letter(cat.get('subclass'))
            if letter is None:
                continue
            seen[letter] += 1
            item = (os.path.join(ROOT, rec['path']), z, letter, str(cat.get('specObjID') or rec.get('specobjid') or ''))
            bucket = buckets[letter]
            if len(bucket) < cap:
                bucket.append(item)
            else:
                j = int(rng.integers(0, seen[letter]))
                if j < cap:
                    bucket[j] = item
    jobs = []
    for letter in LETTERS:
        for path, z, L, sid in buckets[letter]:
            jobs.append((path, z, L, sid, 'test' if _holdout(sid) else 'train'))
    return jobs, seen


def _read_star(item):
    path, z, letter, sid, split = item
    try:
        from ogfkit.spectra import read_spectrum
        from ogfkit.speclearn import on_grid, STAR_WAVE
        s = read_spectrum(path)
        err = s.get('err')
        x = on_grid(s['wave'], s['flux'], err, z, STAR_WAVE)
        row = dict(letter=letter, sid=sid, split=split, ok=x is not None, error=None, x=None, z=z)
        if x is None:
            return row
        if split == 'train':
            row['x'] = x
        else:
            row['wave'] = np.asarray(s['wave'], float)
            row['flux'] = np.asarray(s['flux'], float)
            row['err'] = None if err is None else np.asarray(err, float)
        return row
    except Exception as ex:
        return dict(letter=letter, sid=sid, split=split, ok=False, error=str(ex)[:240], x=None, z=None)


def train_stars(jobs, out_dir):
    import multiprocessing as mp
    from ogfkit.speclearn import fit_star, classify_star, save_model, LETTERS
    t0 = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(2) as pool:
        rows = pool.map(_read_star, jobs, chunksize=8)
    train_x, train_y, test = [], [], []
    n_drop = n_err = 0
    for row in rows:
        if row['error']:
            n_err += 1
            continue
        if not row['ok']:
            n_drop += 1
            continue
        if row['split'] == 'train':
            train_x.append(row['x'])
            train_y.append(row['letter'])
        else:
            test.append(row)
    print('stars read %d train %d test %d drop %d err %d in %.1fs' % (
        len(rows), len(train_y), len(test), n_drop, n_err, time.time() - t0), flush=True)
    if len(train_y) < 4:
        raise RuntimeError('star training kept %d spectra' % len(train_y))
    model = fit_star(np.vstack(train_x), train_y, n_comp=30, seed=1)
    path = os.path.join(out_dir, 'star.npz')
    save_model(path, model)
    per = {L: dict(train=int(np.sum(np.array(train_y) == L)), test=0, correct=0) for L in LETTERS}
    for row in test:
        an = classify_star(row['wave'], row['flux'], row['err'], row['z'], model)
        L = row['letter']
        per[L]['test'] += 1
        per[L]['correct'] += int(an['type'] == L)
    n_test = sum(v['test'] for v in per.values())
    n_correct = sum(v['correct'] for v in per.values())
    return dict(path=path, n_test=n_test, n_correct=n_correct, per_letter=per,
                n_train=len(train_y), n_drop=n_drop, n_error=n_err, note=model['note'])


def _desi_spectra(path, limit, rng):
    from ogfkit.speclearn import on_grid, GAL_WAVE
    out = []
    with np.load(path, allow_pickle=True) as blob:
        offs = blob['offsets']
        wave, flux, ivar, mask = blob['wavelength'], blob['flux'], blob['ivar'], blob['mask']
        meta = blob['metadata_json']
        order = rng.permutation(len(meta))
        for i in order:
            if len(out) >= limit:
                break
            raw = meta[int(i)]
            if isinstance(raw, bytes):
                raw = raw.decode('utf-8')
            info = json.loads(raw)
            if (info.get('spectype') or '').upper() != 'GALAXY':
                continue
            if int(info.get('redshift_warning') or 0) != 0:
                continue
            try:
                z = float(info['redshift'])
            except (TypeError, ValueError, KeyError):
                continue
            if not np.isfinite(z) or not (0.02 < z < 0.5):
                continue
            a, b = int(offs[int(i)]), int(offs[int(i) + 1])
            w = np.array(wave[a:b], dtype=float, copy=True)
            f = np.array(flux[a:b], dtype=float, copy=True)
            iv = np.array(ivar[a:b], dtype=float, copy=True)
            m = np.array(mask[a:b], copy=True)
            bad = (m != 0) | ~(iv > 0) | ~np.isfinite(f) | ~np.isfinite(w)
            err = np.full(f.shape, np.nan)
            err[~bad] = 1.0 / np.sqrt(iv[~bad])
            f[bad] = np.nan
            x = on_grid(w, f, err, z, GAL_WAVE)
            if x is None:
                continue
            tid = str(info.get('targetid') or info.get('specid') or '%s:%d' % (os.path.basename(path), int(i)))
            out.append(dict(w=w, f=f, err=err, z=z, tid=tid, x=x))
    return out


def collect_galaxies(n_unlab, n_lib, n_hold, per_file, seed):
    if not os.path.isdir(DESI_DIR):
        raise FileNotFoundError('DESI directory was not found: %s' % DESI_DIR)
    rng = np.random.default_rng(seed)
    paths = [ent.path for ent in os.scandir(DESI_DIR) if ent.name.startswith('galaxy-') and ent.name.endswith('.npz')]
    order = rng.permutation(len(paths))
    unlab, lib, hold = [], [], []
    opened = 0
    for j in order:
        if len(unlab) >= n_unlab and len(lib) >= n_lib and len(hold) >= n_hold:
            break
        opened += 1
        try:
            specs = _desi_spectra(paths[int(j)], per_file, rng)
        except Exception as ex:
            print('desi skip', os.path.basename(paths[int(j)]), ex, flush=True)
            continue
        for spec in specs:
            test = _holdout(spec['tid'])
            if test and len(hold) < n_hold:
                hold.append(spec)
            elif (not test) and len(lib) < n_lib:
                lib.append(spec)
            elif len(unlab) < n_unlab:
                unlab.append(spec['x'])
        if opened % 25 == 0:
            print('desi files %d unlab %d lib %d hold %d' % (opened, len(unlab), len(lib), len(hold)), flush=True)
    print('desi files %d unlab %d lib %d hold %d' % (opened, len(unlab), len(lib), len(hold)), flush=True)
    return unlab, lib, hold


_GAL_JOBS = []


def _label_galaxy(i):
    from ogfkit.specscience import galaxy_analysis
    spec = _GAL_JOBS[i]
    try:
        an = galaxy_analysis(spec['w'], spec['f'], spec['err'], frame='vacuum', z_known=float(spec['z']), zmax=1.0)
        return an.get('type') or 'unknown'
    except Exception as ex:
        return 'unknown'


def train_galaxies(unlab, lib, hold, out_dir):
    import multiprocessing as mp
    from ogfkit.speclearn import fit_galaxy, classify_galaxy, save_model
    global _GAL_JOBS
    _GAL_JOBS = list(lib) + list(hold)
    n_lib = len(lib)
    t0 = time.time()
    ctx = mp.get_context('fork')
    with ctx.Pool(2) as pool:
        types = pool.map(_label_galaxy, list(range(len(_GAL_JOBS))), chunksize=4)
    print('galaxy labels %d in %.1fs' % (len(types), time.time() - t0), flush=True)
    lib_x, lib_y = [], []
    hist = {}
    for spec, typ in zip(lib, types[:n_lib]):
        hist[typ] = hist.get(typ, 0) + 1
        if typ and typ != 'unknown':
            lib_x.append(spec['x'])
            lib_y.append(typ)
    hold_rows = []
    for spec, typ in zip(hold, types[n_lib:]):
        hist[typ] = hist.get(typ, 0) + 1
        if typ and typ != 'unknown':
            hold_rows.append((spec, typ))
    # Spectra whose line diagnostic is unknown still teach the autoencoder.
    X = list(unlab) + lib_x
    labels = [''] * len(unlab) + lib_y
    # Keep a copy of the unknown library rows as extra unlabeled examples.
    for spec, typ in zip(lib, types[:n_lib]):
        if not typ or typ == 'unknown':
            X.append(spec['x'])
            labels.append('')
    if len(lib_y) < 4:
        raise RuntimeError('galaxy neighbour library kept %d labelled spectra' % len(lib_y))
    model = fit_galaxy(np.vstack(X), labels, latent=16, hidden=96, epochs=40, batch=128, lr=0.02, k=9, seed=1)
    path = os.path.join(out_dir, 'galaxy.npz')
    save_model(path, model)
    n_ok = 0
    per = {}
    for spec, typ in hold_rows:
        an = classify_galaxy(spec['w'], spec['f'], spec['err'], spec['z'], model)
        per.setdefault(typ, dict(n=0, correct=0))
        per[typ]['n'] += 1
        per[typ]['correct'] += int(an['type'] == typ)
        n_ok += int(an['type'] == typ)
    return dict(path=path, n_unlabeled=len(unlab), n_library=len(lib_y), n_holdout=len(hold_rows),
                n_correct=n_ok, label_counts=hist, per_type=per, loss=model['loss'], note=model['note'])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--out', default='/tmp/speclearn_train')
    ap.add_argument('--star-cap', type=int, default=300)
    ap.add_argument('--unlabeled', type=int, default=1200)
    ap.add_argument('--library', type=int, default=500)
    ap.add_argument('--holdout', type=int, default=160)
    ap.add_argument('--per-file', type=int, default=16)
    ap.add_argument('--seed', type=int, default=20261010)
    a = ap.parse_args(argv)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()
    jobs, seen = select_stars(a.star_cap, a.seed)
    print('star candidates', {L: seen[L] for L in seen}, 'kept', len(jobs), flush=True)
    star = train_stars(jobs, a.out)
    unlab, lib, hold = collect_galaxies(a.unlabeled, a.library, a.holdout, a.per_file, a.seed)
    galaxy = train_galaxies(unlab, lib, hold, a.out)
    report = dict(star=star, galaxy=galaxy, seconds=time.time() - t0)
    # Grid vectors are not part of the written report.
    dest = os.path.join(a.out, 'report.json')
    with open(dest, 'w') as fh:
        json.dump(report, fh, indent=1)
    st = star['n_correct'] / star['n_test'] if star['n_test'] else float('nan')
    gt = galaxy['n_correct'] / galaxy['n_holdout'] if galaxy['n_holdout'] else float('nan')
    print('STAR %d/%d %.4f' % (star['n_correct'], star['n_test'], st), flush=True)
    print('GALAXY %d/%d %.4f' % (galaxy['n_correct'], galaxy['n_holdout'], gt), flush=True)
    print('report', dest, flush=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
