"""Learned spectrum labels that stay inside the published methods.

Stars: a temperature letter (O B A F G K M) from principal components and a
softmax, trained on catalog letters.  Galaxies and AGN: an autoencoder trained
without labels, then a nearest-neighbour vote among spectra whose line
diagnostic already has a type.  Supernovae: continuum-removed cross-correlation
against a directory of template spectra.  No template library is shipped.

The star letter repeats the catalog's temperature class.  The galaxy neighbour
repeats the line diagnostic of its neighbours.  Neither is an MK standard or a
FastSpecFit label.
"""
import os

import numpy as np
from scipy.spatial import cKDTree

from ogfkit.spectra import continuum

C_KMS = 299792.458
STAR_WAVE = np.geomspace(3900.0, 8800.0, 640)
GAL_WAVE = np.geomspace(4200.0, 6900.0, 480)
SN_WAVE = np.geomspace(4500.0, 7000.0, 400)
LETTERS = ('O', 'B', 'A', 'F', 'G', 'K', 'M')
MIN_COVER = 0.75
SN_MIN_R = 0.40
NOTE_STAR = 'Temperature letter from catalog subclasses. Not an MK luminosity class.'
NOTE_GALAXY = 'Nearest neighbours in an autoencoder latent space. Names are the line diagnostic of those neighbours.'
NOTE_SN = 'Cross-correlation against the template directory. No template library is included.'


def temperature_letter(subclass):
    """First MK temperature letter, or None for white dwarfs, CVs, and blank labels."""
    s = (subclass or '').strip().upper()
    if not s or s.startswith('WD') or s.startswith('CV'):
        return None
    return s[0] if s[0] in LETTERS else None


def _finite_pairs(wave, flux, err):
    w = np.asarray(wave, float).ravel()
    f = np.asarray(flux, float).ravel()
    good = np.isfinite(w) & np.isfinite(f) & (w > 0)
    if err is not None:
        e = np.asarray(err, float).ravel()
        if e.shape == f.shape:
            good &= np.isfinite(e) & (e > 0)
    w, f = w[good], f[good]
    if w.size < 20:
        return w, f
    order = np.argsort(w, kind='mergesort')
    w, f = w[order], f[order]
    # interp needs strictly increasing samples
    keep = np.concatenate([[True], np.diff(w) > 0])
    return w[keep], f[keep]


def on_grid(wave, flux, err, z, grid):
    """Rest-frame flux on ``grid``, divided by its median absolute value.

    Returns None when too little of the grid is covered.
    """
    z = float(z)
    if not np.isfinite(z) or z <= -0.5:
        return None
    w, f = _finite_pairs(wave, flux, err)
    if w.size < 20:
        return None
    wr = w / (1.0 + z)
    out = np.interp(grid, wr, f, left=np.nan, right=np.nan)
    cover = float(np.mean(np.isfinite(out)))
    if cover < MIN_COVER:
        return None
    med = float(np.nanmedian(np.abs(out)))
    if not np.isfinite(med) or med <= 0:
        return None
    out = out / med
    # A narrow spike (a residual cosmic ray or a very strong line) otherwise
    # dominates the autoencoder update.  Continuum stays near 1.
    out = np.clip(out, -5.0, 25.0)
    # Missing edges stay at 0 so the vector is finite.  Coverage already passed.
    return np.where(np.isfinite(out), out, 0.0)


def _softmax(z):
    z = z - np.max(z, axis=1, keepdims=True)
    e = np.exp(z)
    return e / np.sum(e, axis=1, keepdims=True)


def _fit_softmax(X, y, n_class, n_iter=600, lr=0.35, l2=1e-3, seed=1):
    rng = np.random.default_rng(seed)
    n, d = X.shape
    W = rng.normal(0.0, 0.01, size=(d, n_class))
    b = np.zeros(n_class)
    Y = np.eye(n_class)[y]
    counts = np.bincount(y, minlength=n_class).astype(float)
    class_w = counts.sum() / (np.maximum(counts, 1.0) * n_class)
    sw = class_w[y][:, None]
    for _ in range(n_iter):
        P = _softmax(X @ W + b)
        G = (P - Y) * sw / n
        W -= lr * (X.T @ G + l2 * W)
        b -= lr * G.sum(axis=0)
    return W, b


def fit_star(X, letters, n_comp=30, seed=1):
    """PCA plus a softmax.  ``letters`` are temperature letters, one per row of X."""
    X = np.asarray(X, float)
    letters = np.asarray(list(letters))
    if X.ndim != 2 or X.shape[0] < 4:
        raise ValueError('star training needs at least 4 spectra')
    classes = np.array([c for c in LETTERS if np.any(letters == c)])
    if classes.size < 2:
        raise ValueError('star training needs at least two temperature letters')
    y = np.array([int(np.where(classes == L)[0][0]) for L in letters])
    mu = X.mean(axis=0)
    _, _, vt = np.linalg.svd(X - mu, full_matrices=False)
    n_comp = int(min(n_comp, vt.shape[0], X.shape[0] - 1))
    comp = vt[:n_comp]
    Z = (X - mu) @ comp.T
    zmu = Z.mean(axis=0)
    zsd = Z.std(axis=0)
    zsd = np.where(zsd > 0, zsd, 1.0)
    Zn = (Z - zmu) / zsd
    W, b = _fit_softmax(Zn, y, classes.size, seed=seed)
    return dict(kind='star', wave=STAR_WAVE.copy(), mean=mu, components=comp, z_mean=zmu, z_scale=zsd,
                weight=W, bias=b, classes=classes, note=NOTE_STAR)


def _star_proba(x, model):
    z = ((x - model['mean']) @ model['components'].T - model['z_mean']) / model['z_scale']
    p = _softmax(z.reshape(1, -1) @ model['weight'] + model['bias'])[0]
    i = int(np.argmax(p))
    return str(model['classes'][i]), float(p[i])


def _xavier(rng, fan_in, fan_out):
    return rng.normal(0.0, 1.0 / np.sqrt(fan_in), size=(fan_in, fan_out))


def _ae_forward(X, w):
    h1 = np.tanh(X @ w['W1'] + w['b1'])
    lat = np.tanh(h1 @ w['W2'] + w['b2'])
    h3 = np.tanh(lat @ w['W3'] + w['b3'])
    y = h3 @ w['W4'] + w['b4']
    return h1, lat, h3, y


def _apply(w, key, grad, lr, limit=1.0):
    step = lr * grad
    peak = float(np.max(np.abs(step))) if step.size else 0.0
    if np.isfinite(peak) and peak > limit:
        step = step * (limit / peak)
    w[key] = w[key] - step


def _ae_step(X, w, lr):
    h1, lat, h3, y = _ae_forward(X, w)
    n = X.shape[0]
    dy = (y - X) / n
    _apply(w, 'W4', h3.T @ dy, lr)
    _apply(w, 'b4', dy.sum(axis=0), lr)
    dh3 = (dy @ w['W4'].T) * (1.0 - h3 * h3)
    _apply(w, 'W3', lat.T @ dh3, lr)
    _apply(w, 'b3', dh3.sum(axis=0), lr)
    dlat = (dh3 @ w['W3'].T) * (1.0 - lat * lat)
    _apply(w, 'W2', h1.T @ dlat, lr)
    _apply(w, 'b2', dlat.sum(axis=0), lr)
    dh1 = (dlat @ w['W2'].T) * (1.0 - h1 * h1)
    _apply(w, 'W1', X.T @ dh1, lr)
    _apply(w, 'b1', dh1.sum(axis=0), lr)
    loss = float(np.mean((y - X) ** 2))
    if not np.isfinite(loss):
        raise RuntimeError('galaxy autoencoder loss is not finite')
    return loss


def fit_galaxy(X, labels, latent=16, hidden=96, epochs=25, batch=128, lr=0.08, k=9, seed=1):
    """Autoencoder on ``X``.  ``labels`` name the rows that enter the neighbour library.

    A row with an empty label is used only to train the autoencoder.
    """
    X = np.asarray(X, float)
    labels = np.asarray(list(labels), dtype=object)
    if X.shape[0] < 8:
        raise ValueError('galaxy training needs at least 8 spectra')
    rng = np.random.default_rng(seed)
    n, d = X.shape
    h, m = int(hidden), int(latent)
    w = dict(W1=_xavier(rng, d, h), b1=np.zeros(h), W2=_xavier(rng, h, m), b2=np.zeros(m),
             W3=_xavier(rng, m, h), b3=np.zeros(h), W4=_xavier(rng, h, d), b4=np.zeros(d))
    order = np.arange(n)
    last = None
    for _ in range(int(epochs)):
        rng.shuffle(order)
        for i in range(0, n, batch):
            last = _ae_step(X[order[i:i + batch]], w, lr)
    named = np.array([bool(s) for s in labels])
    if named.sum() < 4:
        raise ValueError('galaxy neighbour library needs at least 4 labelled spectra')
    _, lat, _, _ = _ae_forward(X[named], w)
    kk = int(min(k, int(named.sum())))
    w.update(kind='galaxy', wave=GAL_WAVE.copy(), latent=lat, labels=labels[named].astype(str),
             k=kk, loss=last, note=NOTE_GALAXY)
    return w


def _encode(x, model):
    return _ae_forward(x.reshape(1, -1), model)[1][0]


def _neighbour_vote(x, model):
    tree = cKDTree(model['latent'])
    kk = int(model['k'])
    dist, idx = tree.query(_encode(x, model), k=kk)
    idx = np.atleast_1d(idx)
    dist = np.atleast_1d(dist)
    names = [str(model['labels'][i]) for i in idx]
    counts = {}
    for name in names:
        counts[name] = counts.get(name, 0) + 1
    best = max(counts, key=lambda n: (counts[n], -names.index(n)))
    return best, counts[best] / float(len(names)), float(np.mean(dist)), kk


def save_model(path, model):
    os.makedirs(os.path.dirname(os.path.abspath(path)) or '.', exist_ok=True)
    kind = model['kind']
    if kind == 'star':
        np.savez_compressed(path, kind=np.array('star'), wave=model['wave'], mean=model['mean'],
                            components=model['components'], z_mean=model['z_mean'], z_scale=model['z_scale'],
                            weight=model['weight'], bias=model['bias'], classes=np.array(model['classes']),
                            note=np.array(model['note']))
    elif kind == 'galaxy':
        np.savez_compressed(path, kind=np.array('galaxy'), wave=model['wave'],
                            W1=model['W1'], b1=model['b1'], W2=model['W2'], b2=model['b2'],
                            W3=model['W3'], b3=model['b3'], W4=model['W4'], b4=model['b4'],
                            latent=model['latent'], labels=np.array(model['labels']),
                            k=np.array(model['k']), note=np.array(model['note']))
    else:
        raise ValueError('unknown model kind %r' % kind)


def load_model(path):
    if not os.path.isfile(path):
        raise FileNotFoundError('model was not found: %s' % path)
    z = np.load(path, allow_pickle=False)
    kind = str(z['kind'])
    if kind == 'star':
        return dict(kind='star', wave=z['wave'], mean=z['mean'], components=z['components'],
                    z_mean=z['z_mean'], z_scale=z['z_scale'], weight=z['weight'], bias=z['bias'],
                    classes=z['classes'], note=str(z['note']))
    if kind == 'galaxy':
        return dict(kind='galaxy', wave=z['wave'], W1=z['W1'], b1=z['b1'], W2=z['W2'], b2=z['b2'],
                    W3=z['W3'], b3=z['b3'], W4=z['W4'], b4=z['b4'], latent=z['latent'],
                    labels=z['labels'], k=int(z['k']), note=str(z['note']))
    raise ValueError('unknown model kind %r' % kind)


def classify_star(wave, flux, err, z, model):
    x = on_grid(wave, flux, err, 0.0 if z is None else z, model['wave'])
    if x is None:
        return dict(type='unknown', p=None, quality=0, note=model['note'])
    letter, p = _star_proba(x, model)
    return dict(type=letter, p=p, quality=1 if p >= 0.45 else 0, note=model['note'])


def classify_galaxy(wave, flux, err, z, model):
    x = on_grid(wave, flux, err, z, model['wave'])
    if x is None or z is None or not np.isfinite(z):
        return dict(type='unknown', p=None, k=0, distance=None, quality=0, note=model['note'])
    name, frac, dist, kk = _neighbour_vote(x, model)
    return dict(type=name, p=frac, k=kk, distance=dist, quality=1 if frac >= 0.5 else 0, note=model['note'])


def _continuum_residual(wave, flux):
    cont, _ = continuum(wave, flux, width=101)
    r = np.asarray(flux, float) - cont
    r = np.where(np.isfinite(r), r, 0.0)
    s = float(np.sqrt(np.mean(r * r)))
    if s <= 0:
        return None
    return r / s


def load_templates(directory):
    """Text spectra in ``directory``.  The type is the ``# type`` header, else the file prefix."""
    if not directory or not os.path.isdir(directory):
        raise FileNotFoundError('template directory was not found: %s' % directory)
    out = []
    for name in sorted(os.listdir(directory)):
        if not name.lower().endswith(('.txt', '.dat')):
            continue
        path = os.path.join(directory, name)
        kind = name.split('_', 1)[0].split('.', 1)[0]
        rows = []
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if not s:
                    continue
                if s.startswith('#'):
                    parts = s[1:].split()
                    if len(parts) >= 2 and parts[0].lower() == 'type':
                        kind = parts[1]
                    continue
                bits = s.split()
                if len(bits) >= 2:
                    rows.append((float(bits[0]), float(bits[1])))
        if len(rows) < 20:
            continue
        arr = np.asarray(rows, float)
        out.append(dict(type=kind, wave=arr[:, 0], flux=arr[:, 1], file=name))
    if not out:
        raise FileNotFoundError('template directory was not found: %s' % directory)
    return out


def match_supernova(wave, flux, err, z, templates, v_min=-4000.0, v_max=4000.0, v_step=500.0):
    """Highest continuum-removed correlation against ``templates`` over a velocity lag.

    A host redshift is required.  Line velocities are only the lag around that redshift.
    """
    if z is None or not np.isfinite(z) or z <= -0.5:
        return dict(type='unknown', score=None, v_lag=None, template=None, quality=0, note=NOTE_SN)
    w, f = _finite_pairs(wave, flux, err)
    if w.size < 20:
        return dict(type='uncertain', score=None, v_lag=None, template=None, quality=0, note=NOTE_SN)
    wr = w / (1.0 + float(z))
    grid = SN_WAVE
    obs = np.interp(grid, wr, f, left=np.nan, right=np.nan)
    if float(np.mean(np.isfinite(obs))) < MIN_COVER:
        return dict(type='uncertain', score=None, v_lag=None, template=None, quality=0, note=NOTE_SN)
    obs = np.where(np.isfinite(obs), obs, np.nanmedian(obs[np.isfinite(obs)]))
    obs_r = _continuum_residual(grid, obs)
    if obs_r is None:
        return dict(type='uncertain', score=None, v_lag=None, template=None, quality=0, note=NOTE_SN)
    best = (-1.0, None, None, None)
    lags = np.arange(v_min, v_max + 0.5 * v_step, v_step)
    for tpl in templates:
        tw, tf = _finite_pairs(tpl['wave'], tpl['flux'], None)
        if tw.size < 20:
            continue
        for v in lags:
            scale = 1.0 - float(v) / C_KMS
            shifted = np.interp(grid, tw * scale, tf, left=np.nan, right=np.nan)
            if float(np.mean(np.isfinite(shifted))) < MIN_COVER:
                continue
            med = float(np.nanmedian(shifted))
            shifted = np.where(np.isfinite(shifted), shifted, med)
            tr = _continuum_residual(grid, shifted)
            if tr is None:
                continue
            score = float(np.mean(obs_r * tr))
            if score > best[0]:
                best = (score, str(tpl['type']), float(v), tpl.get('file'))
    score, kind, v_lag, fname = best
    if kind is None or score < SN_MIN_R:
        return dict(type='uncertain', score=None if kind is None else score, v_lag=v_lag,
                    template=fname, quality=0, note=NOTE_SN)
    return dict(type=kind, score=score, v_lag=v_lag, template=fname, quality=1, note=NOTE_SN)


def model_note(model):
    return model.get('note') or ''
