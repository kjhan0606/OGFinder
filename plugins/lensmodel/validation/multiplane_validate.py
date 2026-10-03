"""Validation of the multi-source-plane fit on synthetic double-source-plane lenses (SIE + shear, quads in both planes)."""
import os, sys, json, math, time
import numpy as np
ROOT = os.environ.get('OGF_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..')))
sys.path.insert(0, ROOT)
from ogfkit import lensmodel as L, lenssynth as S, lensextra as X


def make(rng, zl=0.222, z1=0.609, z2=2.035, noise=0.005):
    w2 = X.dist_weight(zl, z2, z1)
    for _ in range(100):
        d = S.random_lens(rng, quad=True, theta_E=rng.uniform(1.0, 1.8))
        m = L.build_sie_shear(d['params'])
        g1 = d['images']; b1 = d['source']
        # second plane: weight w2 > 1; source near centre of its caustic
        sm = X.ScaledLens(m, w2)
        for _t in range(50):
            r = 0.12 * rng.uniform(0.3, 1) * d['params']['theta_E'] * w2
            a = rng.uniform(0, 2 * math.pi)
            b2 = (r * math.cos(a), r * math.sin(a))
            g2 = sm.solve_images(b2[0], b2[1], half=3 * d['params']['theta_E'] * w2, n=301)
            if len(g2) == 4:
                break
        else:
            continue
        xy1 = np.array([[i['x'], i['y']] for i in g1]); xy2 = np.array([[i['x'], i['y']] for i in g2])
        return d['params'], w2, xy1 + rng.normal(0, noise, xy1.shape), xy2 + rng.normal(0, noise, xy2.shape)
    raise RuntimeError


def angd(a, b):
    return (a - b + 90) % 180 - 90


def run(n=24, noise=0.005, seed=5):
    rng = np.random.default_rng(seed)
    zl, z1, z2 = 0.222, 0.609, 2.035
    rows = []
    t0 = time.time()
    for i in range(n):
        p, w2, a, b = make(rng, zl, z1, z2, noise)
        c = (0.0, 0.0)
        out = {}
        out['one'] = X.fit_multiplane([dict(xy=a, z=z1)], c, sigma=noise, n_starts=6)
        out['two_alone'] = X.fit_multiplane([dict(xy=b, z=z2)], c, sigma=noise, n_starts=6)
        out['joint_fixed'] = X.fit_multiplane([dict(xy=a, z=z1), dict(xy=b, z=z2)], c, sigma=noise, z_l=zl, n_starts=6)
        out['joint_free'] = X.fit_multiplane([dict(xy=a, z=z1), dict(xy=b, z=None)], c, sigma=noise, n_starts=8)
        out['same_plane'] = X.fit_multiplane([dict(xy=np.vstack([a, b]), z=z1)], c, sigma=noise, n_starts=6)   # wrongly merge both planes
        row = dict(truth=p, w2=w2)
        for k, r in out.items():
            pp = r['params']
            row[k] = dict(dthE=(pp['theta_E'] - p['theta_E']) / p['theta_E'] if k != 'two_alone' else (pp['theta_E'] / w2 - p['theta_E']) / p['theta_E'],
                          dq=pp['q'] - p['q'], dphi=angd(pp['phi'], p['phi']), dgam=pp['gamma'] - p['gamma'], chi2=r['chi2'], dof=r['dof'])
            if k == 'joint_free':
                row[k]['w2'] = r['weights'][1]; row[k]['w2err'] = r['errors'].get('w1')
        rows.append(row)
        print(i, round(time.time() - t0), flush=True)
    return rows


def summarize(rows):
    s = {}
    for k in ('one', 'two_alone', 'joint_fixed', 'joint_free', 'same_plane'):
        v = lambda f: np.array([r[k][f] for r in rows])
        s[k] = dict(thE_sig=float(np.std(v('dthE'))), thE_bias=float(np.mean(v('dthE'))), q_sig=float(np.std(v('dq'))), phi_sig=float(np.std(v('dphi'))),
                    gam_sig=float(np.std(v('dgam'))), chi2_med=float(np.median(v('chi2'))), dof=int(np.median(v('dof'))),
                    rms_dthE_max=float(np.max(np.abs(v('dthE')))))
    w = np.array([r['joint_free']['w2'] for r in rows]); wt = np.array([r['w2'] for r in rows])
    s['w_ratio'] = dict(median=float(np.median(w / wt)), sig=float(np.std(w / wt)), true_w2_median=float(np.median(wt)))
    return s


if __name__ == '__main__':
    out = sys.argv[1]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 24
    rows = run(n, seed=int(sys.argv[3]) if len(sys.argv) > 3 else 5)
    s = summarize(rows)
    json.dump(dict(n=len(rows), noise_arcsec=0.005, zl=0.222, z1=0.609, z2=2.035, summary=s, rows=rows), open(out, 'w'), indent=1, default=float)
    print(json.dumps(s, indent=0)[:2500])
