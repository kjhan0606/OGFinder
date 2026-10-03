"""Strong-lens extensions: multiple source planes and pixelated (regularised, linear) source reconstruction.

* ``dist_weight`` / ``ScaledLens``: a source at redshift z_s sees the deflection scaled by  w = (D_ls/D_s)(z_s) / (D_ls/D_s)(z_ref).
* ``fit_multiplane``: image-plane fit of SIE + shear to several sets of multiple images of sources at *different* redshifts that
  share one lens.  Each set has its own (nuisance) source position; the weight of a set is fixed from its redshift or fitted.
* ``SourceInversion``: pixelated source on a regular grid, bilinear lensing operator + PSF, Gaussian noise, Tikhonov (zeroth /
  gradient / curvature) regularisation.  The regularisation strength is chosen by maximising the Bayesian evidence; the same evidence
  can be maximised over the lens parameters (``fit_lens_pixels``).  numpy / scipy only.
"""
import math
import numpy as np
from . import lensmodel as LM


# ------------------------------------------------------------------------------------------------------------------- multiple source planes
def dist_weight(z_l, z_s, z_ref, H0=70.0, Om=0.3):
    """deflection scaling of a source at z_s relative to a source at z_ref: (D_ls/D_s)(z_s) / (D_ls/D_s)(z_ref)"""
    _, Ds, Dls = LM.distances(z_l, z_s, H0, Om)
    _, Dr, Dlr = LM.distances(z_l, z_ref, H0, Om)
    return (Dls / Ds) / (Dlr / Dr)


class ScaledLens(LM.LensModel):
    """the same mass model, but deflection angles multiplied by w (source plane at a different redshift)"""

    def __init__(self, model, w):
        super().__init__(model.comps)
        self.w = float(w)

    def alpha(self, x, y):
        ax, ay = super().alpha(x, y)
        return self.w * ax, self.w * ay

    def psi(self, x, y):
        return self.w * super().psi(x, y)


_PN = ['theta_E', 'q', 'phi', 'gamma', 'phi_g']


def fit_multiplane(groups, centre, sigma=0.005, z_l=None, z_ref=None, free=None, start=None, fixed=None, n_starts=10, seed=0, H0=70.0, Om=0.3):
    """Joint image-plane fit of an SIE + shear lens to several image sets.

    groups: list of dict(xy=Nx2 image positions (arcsec), z=source redshift or None (weight fitted), sigma=optional per-group position error)
    The first group defines the reference plane (weight 1, theta_E refers to it).  With redshifts given (and z_l) the weights are fixed from
    the cosmology (flat LCDM, H0, Om); a group with z=None gets a fitted weight w_k (> 0).
    -> dict(params, weights, sources, chi2, dof, rms_arcsec, errors, success)
    """
    from scipy.optimize import least_squares
    free = list(free or _PN)
    fixed = dict(fixed or {})
    G = len(groups)
    if start is None and G > 1:                       # start from the single-plane fit of the first set (the other planes are then a refinement)
        r0 = fit_multiplane(groups[:1], centre, sigma=sigma, free=free, fixed=fixed, n_starts=n_starts, seed=seed)
        start = {k: float(r0['params'][k]) for k in _PN}
    obs = [np.asarray(g['xy'], float).reshape(-1, 2) for g in groups]
    sig = [float(g.get('sigma', sigma)) for g in groups]
    wfix = [1.0] + [None] * (G - 1)
    for k in range(1, G):
        z = groups[k].get('z')
        if z is not None and z_l is not None:
            wfix[k] = dist_weight(z_l, z, z_ref if z_ref is not None else groups[0]['z'], H0, Om)
    wfree = [k for k in range(1, G) if wfix[k] is None]
    x0, y0 = centre
    d0 = dict(theta_E=float(np.mean(np.hypot(obs[0][:, 0] - x0, obs[0][:, 1] - y0))), q=0.85, phi=0.0, gamma=0.0, phi_g=0.0)
    d0.update(start or {})
    d0.update(fixed)
    lo = dict(theta_E=0.05, q=0.2, phi=-360.0, gamma=0.0, phi_g=-360.0)
    hi = dict(theta_E=30.0, q=1.0, phi=360.0, gamma=0.5, phi_g=360.0)

    def unpack(v):
        p = dict(d0)
        for i, n in enumerate(free):
            p[n] = v[i]
        j = len(free)
        w = list(wfix)
        for k in wfree:
            w[k] = v[j]
            j += 1
        src = [(v[j + 2 * k], v[j + 2 * k + 1]) for k in range(G)]
        return p, w, src

    def model_of(p):
        return LM.build_sie_shear(dict(p, x0=x0, y0=y0))

    def resid(v):
        p, w, src = unpack(v)
        m = model_of(p)
        out = []
        for k in range(G):
            sm = ScaledLens(m, w[k])
            pred = LM._local_images(sm, obs[k], src[k][0], src[k][1])
            r = (obs[k] - pred) / sig[k]
            r[~np.isfinite(r)] = 100.0
            out.append(r.ravel())
        return np.concatenate(out)

    rng = np.random.default_rng(seed)
    best = None
    trials = []
    n_starts = max(n_starts, 20) if G > 1 else n_starts
    for s in range(n_starts):
        if s == 0 or G == 1:
            v0 = [d0[n] * (1 + (0.05 * rng.standard_normal() if s else 0.0)) if n in ('theta_E', 'q') else d0[n] + (rng.uniform(-5, 5) if s and n in ('phi', 'phi_g') else (0.02 * rng.random() if s and n == 'gamma' else 0.0))
                  for n in free]
        else:                                          # several planes: the first set alone is degenerate, so draw the shape parameters widely
            draw = dict(theta_E=d0['theta_E'] * (1 + 0.05 * rng.standard_normal()), q=rng.uniform(0.5, 1.0), phi=rng.uniform(0, 180), gamma=rng.uniform(0, 0.1), phi_g=rng.uniform(0, 180))
            v0 = [draw[n] for n in free]
        v0 = [min(max(a, lo[n] + 1e-6), hi[n] - 1e-6) for a, n in zip(v0, free)]
        w0 = [wfix[0]]
        for k in wfree:
            v0.append(1.0 + (0.1 * rng.standard_normal() if s else 0.0))
        p0 = dict(d0)
        p0.update({n: v0[i] for i, n in enumerate(free)})
        m0 = model_of(p0)
        for k in range(G):
            wk = wfix[k] if wfix[k] is not None else 1.0
            bx, by = ScaledLens(m0, wk).beta(obs[k][:, 0], obs[k][:, 1])
            v0 += [float(np.mean(bx)), float(np.mean(by))]
        lb = [lo[n] for n in free] + [0.05] * len(wfree) + [-1e3] * (2 * G)
        ub = [hi[n] for n in free] + [20.0] * len(wfree) + [1e3] * (2 * G)
        try:
            r = least_squares(resid, v0, bounds=(lb, ub), x_scale='jac', max_nfev=60)      # short run from every start ...
        except Exception:
            continue
        trials.append((r.cost, v0, lb, ub, r))
    trials.sort(key=lambda t: t[0])
    for cost, v0, lb, ub, r0 in trials[:3]:                                                  # ... the best three are polished
        try:
            r = least_squares(resid, r0.x, bounds=(lb, ub), x_scale='jac', max_nfev=400)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    r = best
    p, w, src = unpack(r.x)
    N = sum(2 * len(o) for o in obs)
    npar = len(r.x)
    dof = N - npar
    chi2 = float(2 * r.cost)
    J = r.jac
    cov = np.linalg.pinv(J.T @ J) * max(1.0, chi2 / dof if dof > 0 else 1.0)
    err = np.sqrt(np.clip(np.diag(cov), 0, None))
    names = free + ['w%d' % k for k in wfree]
    errors = {n: float(err[i]) for i, n in enumerate(names)}
    rms = float(np.sqrt(np.mean((resid(r.x) * np.repeat(sig, [2 * len(o) for o in obs])) ** 2) * 2))   # arcsec, radial (both coordinates)
    return dict(params=p, weights=w, sources=src, chi2=chi2, dof=dof, rms_arcsec=rms, errors=errors, success=bool(r.success))


# ------------------------------------------------------------------------------------------------------------------- pixelated source
def psf_matrix(kernel, shape):
    """sparse (Nimg x Nimg) convolution with a normalised kernel (zero padding)"""
    import scipy.sparse as sp
    ny, nx = shape
    k = np.asarray(kernel, float)
    k = k / k.sum()
    hy, hx = k.shape[0] // 2, k.shape[1] // 2
    yy, xx = np.mgrid[0:ny, 0:nx]
    idx = np.arange(ny * nx).reshape(ny, nx)
    rows, cols, vals = [], [], []
    for dy in range(-hy, hy + 1):
        for dx in range(-hx, hx + 1):
            v = k[dy + hy, dx + hx]
            if v == 0:
                continue
            sy, sx = yy + dy, xx + dx
            ok = (sy >= 0) & (sy < ny) & (sx >= 0) & (sx < nx)
            rows.append(idx[ok])
            cols.append(idx[sy[ok], sx[ok]])
            vals.append(np.full(int(ok.sum()), v))
    return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ny * nx, ny * nx))


def gaussian_kernel(sigma_pix, size=None):
    s = max(float(sigma_pix), 1e-3)
    h = size // 2 if size else int(math.ceil(4 * s))
    x = np.arange(-h, h + 1)
    k = np.exp(-0.5 * (x[:, None] ** 2 + x[None, :] ** 2) / s ** 2)
    return k / k.sum()


def reg_matrix(n, kind='gradient'):
    """regularisation matrix R (n^2 x n^2, sparse):  s^T R s = |D s|^2"""
    import scipy.sparse as sp
    I = sp.identity(n, format='csr')
    if kind == 'zero':
        return sp.identity(n * n, format='csr')
    if kind == 'gradient':
        D = sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n))
        Dx = sp.kron(I, D)
        Dy = sp.kron(D, I)
        return (Dx.T @ Dx + Dy.T @ Dy).tocsr()
    if kind == 'curvature':
        D = sp.diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n))
        Dxx = sp.kron(I, D)
        Dyy = sp.kron(D, I)
        return (Dxx.T @ Dxx + Dyy.T @ Dyy).tocsr()
    raise ValueError(kind)


def lens_operator(model, shape, pixscale, origin_xy, centre_xy, grid, oversample=1):
    """sparse bilinear operator L (Nimg x n^2): image pixel -> source-plane pixel weights (average over oversample^2 sub-pixel rays).
    grid: dict(n, x0, y0, pix) (arcsec of the centre of source pixel [0,0])"""
    import scipy.sparse as sp
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    n = grid['n']
    o = int(oversample)
    ii = np.arange(ny * nx)
    acc = None
    for sy in range(o):
        for sx in range(o):
            px = xx + (sx + 0.5) / o - 0.5
            py = yy + (sy + 0.5) / o - 0.5
            tx = (px - origin_xy[0]) * pixscale + centre_xy[0]
            ty = (py - origin_xy[1]) * pixscale + centre_xy[1]
            bx, by = model.beta(tx, ty)
            fx = ((bx - grid['x0']) / grid['pix']).ravel()
            fy = ((by - grid['y0']) / grid['pix']).ravel()
            ix = np.floor(fx).astype(int)
            iy = np.floor(fy).astype(int)
            ax = fx - ix
            ay = fy - iy
            rows, cols, vals = [], [], []
            for dy, dx, w in ((0, 0, (1 - ax) * (1 - ay)), (0, 1, ax * (1 - ay)), (1, 0, (1 - ax) * ay), (1, 1, ax * ay)):
                jx, jy = ix + dx, iy + dy
                ok = (jx >= 0) & (jx < n) & (jy >= 0) & (jy < n) & np.isfinite(w)
                rows.append(ii[ok])
                cols.append(jy[ok] * n + jx[ok])
                vals.append(w[ok] / (o * o))
            M = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(ny * nx, n * n))
            acc = M if acc is None else acc + M
    return acc


def auto_grid(model, mask, shape, pixscale, origin_xy, centre_xy, n=40, margin=1.25):
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    tx = (xx - origin_xy[0]) * pixscale + centre_xy[0]
    ty = (yy - origin_xy[1]) * pixscale + centre_xy[1]
    bx, by = model.beta(tx[mask], ty[mask])
    cx, cy = float(np.median(bx)), float(np.median(by))
    half = margin * max(np.percentile(np.abs(bx - cx), 99), np.percentile(np.abs(by - cy), 99))
    pix = 2 * half / (n - 1)
    return dict(n=n, x0=cx - half, y0=cy - half, pix=pix)


class SourceInversion:
    """Linear source inversion for a given lens model.

    data (image), sigma (scalar or map), mask (bool, pixels used), kernel (PSF stamp or None).  ``solve(lam)`` / ``best()`` give the regularised
    source; ``evidence(lam)`` = -1/2 chi2 - 1/2 lam s^T R s - 1/2 ln det A + 1/2 ln det(lam R)  (constants dropped; comparable between lens models for
    the same data / noise / mask / grid size / regularisation kind)."""

    def __init__(self, model, data, sigma, mask, pixscale, origin_xy, centre_xy, kernel=None, n=40, reg='gradient', grid=None, oversample=2):
        from scipy.linalg import cho_factor
        self.shape = data.shape
        self.mask = np.asarray(mask, bool) & np.isfinite(data)
        self.grid = grid or auto_grid(model, self.mask, data.shape, pixscale, origin_xy, centre_xy, n)
        n = self.grid['n']
        L = lens_operator(model, data.shape, pixscale, origin_xy, centre_xy, self.grid, oversample)
        F = L if kernel is None else psf_matrix(kernel, data.shape) @ L
        sel = np.flatnonzero(self.mask.ravel())
        self.F = F[sel]
        self.d = np.asarray(data, float).ravel()[sel]
        sg = np.broadcast_to(np.asarray(sigma, float), data.shape).ravel()[sel]
        self.w = 1.0 / sg ** 2
        Fw = self.F.multiply(self.w[:, None]).tocsr()
        self.FtWF = (self.F.T @ Fw).toarray()
        self.FtWd = np.asarray(Fw.T @ self.d).ravel()
        self.R = reg_matrix(n, reg).toarray()
        self.Reps = self.R + 1e-6 * np.trace(self.R) / n ** 2 * np.eye(n * n)
        self.ldR = 2 * np.log(np.diag(np.linalg.cholesky(self.Reps))).sum()
        self.ndata = len(sel)
        self.sel = sel
        self._cache = {}

    def solve(self, lam):
        A = self.FtWF + lam * self.Reps
        c = np.linalg.cholesky(A)
        from scipy.linalg import cho_solve
        s = cho_solve((c, True), self.FtWd)
        ldA = 2 * np.log(np.diag(c)).sum()
        res = self.d - self.F @ s
        chi2 = float((self.w * res ** 2).sum())
        pen = float(lam * s @ self.Reps @ s)
        n = self.grid['n']
        ev = -0.5 * chi2 - 0.5 * pen - 0.5 * ldA + 0.5 * (n * n * math.log(lam) + self.ldR)
        return dict(s=s, chi2=chi2, penalty=pen, evidence=float(ev), lam=lam, cho=c)

    def best(self, lo=-6.0, hi=6.0, iters=24, lam0=None, width=1.0):
        """maximise the evidence over log10(lam) (golden section on [lo, hi] after a coarse scan)"""
        scale = float(np.trace(self.FtWF) / max(np.trace(self.Reps), 1e-30))
        if lam0 is not None:                         # local search around a previous optimum (used inside the lens optimiser)
            a, b = math.log10(lam0 / scale) - width, math.log10(lam0 / scale) + width
        else:
            g = np.linspace(lo, hi, 13)
            ev = [self.solve(scale * 10 ** x)['evidence'] for x in g]
            k = int(np.argmax(ev))
            a, b = g[max(k - 1, 0)], g[min(k + 1, len(g) - 1)]
        gr = (math.sqrt(5) - 1) / 2
        c, d = b - gr * (b - a), a + gr * (b - a)
        fc, fd = self.solve(scale * 10 ** c)['evidence'], self.solve(scale * 10 ** d)['evidence']
        for _ in range(iters):
            if fc > fd:
                b, d, fd = d, c, fc
                c = b - gr * (b - a)
                fc = self.solve(scale * 10 ** c)['evidence']
            else:
                a, c, fc = c, d, fd
                d = a + gr * (b - a)
                fd = self.solve(scale * 10 ** d)['evidence']
        r = self.solve(scale * 10 ** (0.5 * (a + b)))
        return r

    def n_eff(self, r):
        """effective number of fitted parameters  tr[(FtWF + lam R)^-1 FtWF]"""
        from scipy.linalg import cho_solve
        return float(np.trace(cho_solve((r['cho'], True), self.FtWF)))

    def images(self, r):
        n = self.grid['n']
        src = r['s'].reshape(n, n)
        model = np.zeros(self.shape).ravel()
        model[self.sel] = self.F @ r['s']
        return src, model.reshape(self.shape)


def fit_lens_pixels(data, sigma, mask, pixscale, origin_xy, centre_xy, start, free=None, kernel=None, n=28, reg='gradient', fixed_lam=None, maxiter=400, xtol=1e-4, oversample=1):
    """maximise the source-marginalised evidence over the SIE + shear parameters (Nelder-Mead).  -> dict(params, evidence, nfev)"""
    from scipy.optimize import minimize
    free = list(free or _PN)
    p0 = dict(start)
    x0, y0 = centre_xy
    state = dict(n=0, grid=None, lam=None)

    def build(v):
        p = dict(p0)
        p.update({k: v[i] for i, k in enumerate(free)})
        return p

    def neg(v):
        p = build(v)
        if not (0.05 < p['theta_E'] < 30 and 0.2 < p['q'] <= 1.0 and 0 <= p.get('gamma', 0) < 0.5):
            return 1e30
        m = LM.build_sie_shear(dict(p, x0=x0, y0=y0))
        try:
            S = SourceInversion(m, data, sigma, mask, pixscale, origin_xy, centre_xy, kernel=kernel, n=n, reg=reg, oversample=oversample)
            r = (S.best(iters=6, lam0=state['lam'], width=0.7) if state['lam'] else S.best(iters=10)) if fixed_lam is None else S.solve(fixed_lam)
            state['lam'] = r['lam']
        except np.linalg.LinAlgError:
            return 1e30
        state['n'] += 1
        return -r['evidence']

    v0 = np.array([p0[k] for k in free], float)
    step = {'theta_E': 0.03, 'q': 0.04, 'phi': 4.0, 'gamma': 0.02, 'phi_g': 10.0}
    simplex = np.vstack([v0] + [v0 + np.eye(len(v0))[i] * step[k] for i, k in enumerate(free)])
    r = minimize(neg, v0, method='Nelder-Mead', options=dict(initial_simplex=simplex, maxiter=maxiter, xatol=xtol, fatol=1e-2))
    return dict(params=build(r.x), evidence=float(-r.fun), nfev=state['n'], success=bool(r.success))


def arc_mask(data, sigma, origin_xy, rmax_pix, nsig=4.0, smooth=1.5, dilate=3):
    """pixels used for the inversion: smoothed data above nsig sigma_smooth (arcs), dilated, inside a circle of radius rmax_pix about the lens"""
    from scipy.ndimage import gaussian_filter, binary_dilation
    d = np.nan_to_num(np.asarray(data, float))
    sm = gaussian_filter(d, smooth)
    sg = float(np.median(sigma)) if np.ndim(sigma) else float(sigma)
    ssm = sg / (2 * math.sqrt(math.pi) * smooth)               # white noise smoothed with a Gaussian of sigma = smooth
    m = sm > nsig * ssm
    m = binary_dilation(m, iterations=int(dilate))
    yy, xx = np.mgrid[0:d.shape[0], 0:d.shape[1]]
    return m & (np.hypot(xx - origin_xy[0], yy - origin_xy[1]) < rmax_pix)
