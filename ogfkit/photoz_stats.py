"""Photo-z distribution tools: PIT / calibration scores, percentile accuracy, sample representativeness (numpy / scipy only).

Predictive distributions (one per object) are described by `Predictive`:
    gauss      zphot, sigma
    split      zphot (median), lo, hi: asymmetric Gaussian from the 16/50/84 percentiles (sigma_lo = p50 - p16, sigma_hi = p84 - p50)
    mixture    pi, mu, sigma  (N x K): Gaussian mixture (the MDN of photo_z)
    grid       zgrid (Nz,), pdf (N x Nz): tabulated p(z)
All are optionally truncated to z >= zmin (default 0) and renormalised.

Point metrics use dz = (zphot - zspec) / (1 + zspec): bias, sigma_NMAD, outlier fraction (|dz| > 0.15), percentiles of |dz|.
Calibration of the PDFs: PIT_i = F_i(zspec_i) ~ U(0,1) if calibrated; KS / Cramer-von Mises tests, histogram, credible-interval coverage,
CRPS and log score, optimal scale of the widths.  Representativeness: 1-D and multi-dimensional comparison of the spectroscopic (training /
calibration) sample with the target sample in feature space, k-NN density-ratio weights (Lima et al. 2008) and weighted accuracy.
"""
import math

import numpy as np
from scipy import special, stats

SQ2 = math.sqrt(2.0)


def ndtr(x):
    return special.ndtr(x)


def _col(a, n=None):
    a = np.asarray(a, float)
    return a


class Predictive:
    def __init__(self, kind, zmin=0.0, **kw):
        self.kind = kind
        self.zmin = zmin
        if kind == 'gauss':
            self.mu = np.asarray(kw['zphot'], float)[:, None]
            self.sig = np.maximum(np.asarray(kw['sigma'], float), 1e-6)[:, None]
            self.pi = np.ones_like(self.mu)
        elif kind == 'mixture':
            self.pi = np.asarray(kw['pi'], float)
            self.pi = self.pi / self.pi.sum(axis=1, keepdims=True)
            self.mu = np.asarray(kw['mu'], float)
            self.sig = np.maximum(np.asarray(kw['sigma'], float), 1e-6)
        elif kind == 'split':
            self.med = np.asarray(kw['zphot'], float)
            self.slo = np.maximum(np.asarray(kw['lo'], float), 1e-6)
            self.shi = np.maximum(np.asarray(kw['hi'], float), 1e-6)
        elif kind == 'grid':
            self.zg = np.asarray(kw['zgrid'], float)
            p = np.asarray(kw['pdf'], float)
            p = np.clip(p, 0, None)
            norm = np.trapezoid(p, self.zg, axis=1)
            self.p = p / np.where(norm > 0, norm, 1.0)[:, None]
            self.cg = np.concatenate([np.zeros((len(p), 1)), np.cumsum(0.5 * (self.p[:, 1:] + self.p[:, :-1]) * np.diff(self.zg), axis=1)], axis=1)
        else:
            raise ValueError('unknown predictive kind %s' % kind)
        self.n = len(self.mu) if kind in ('gauss', 'mixture') else (len(self.med) if kind == 'split' else len(self.p))

    # raw (untruncated) cdf / pdf at z (N,)
    def _cdf_raw(self, z):
        z = np.asarray(z, float)
        if self.kind in ('gauss', 'mixture'):
            return (self.pi * ndtr((z[:, None] - self.mu) / self.sig)).sum(axis=1)
        if self.kind == 'split':
            # asymmetric Gaussian with continuous density: weights proportional to the two sigmas
            w = self.slo + self.shi
            lo = z < self.med
            f_lo = 2 * self.slo / w * ndtr((z - self.med) / self.slo)
            f_hi = self.slo / w * 1.0 * 2 * 0.5 + 2 * self.shi / w * (ndtr((z - self.med) / self.shi) - 0.5)
            return np.where(lo, f_lo, f_hi)
        out = np.empty(len(z))
        for i in range(len(z)):
            out[i] = np.interp(z[i], self.zg, self.cg[i], left=0.0, right=self.cg[i, -1])
        return out

    def _pdf_raw(self, z):
        z = np.asarray(z, float)
        if self.kind in ('gauss', 'mixture'):
            return (self.pi * np.exp(-0.5 * ((z[:, None] - self.mu) / self.sig) ** 2) / (self.sig * math.sqrt(2 * math.pi))).sum(axis=1)
        if self.kind == 'split':
            w = self.slo + self.shi
            s = np.where(z < self.med, self.slo, self.shi)
            return 2.0 / w * np.exp(-0.5 * ((z - self.med) / s) ** 2) / math.sqrt(2 * math.pi)
        out = np.empty(len(z))
        for i in range(len(z)):
            out[i] = np.interp(z[i], self.zg, self.p[i], left=0.0, right=0.0)
        return out

    def _f0(self):
        return self._cdf_raw(np.full(self.n, self.zmin)) if self.zmin is not None else np.zeros(self.n)

    def cdf(self, z):
        f0 = self._f0()
        return np.clip((self._cdf_raw(z) - f0) / np.maximum(1.0 - f0, 1e-12), 0.0, 1.0)

    def pdf(self, z):
        f0 = self._f0()
        z = np.asarray(z, float)
        return np.where(z >= (self.zmin if self.zmin is not None else -np.inf), self._pdf_raw(z) / np.maximum(1.0 - f0, 1e-12), 0.0)

    def grid(self, zgrid):
        """(N, Nz) pdf on a common grid (rows renormalised on the grid)."""
        zgrid = np.asarray(zgrid, float)
        f0 = self._f0()
        P = np.empty((self.n, len(zgrid)))
        for j, z in enumerate(zgrid):
            P[:, j] = self._pdf_raw(np.full(self.n, z))
        P = np.where(zgrid[None, :] >= (self.zmin if self.zmin is not None else -np.inf), P, 0.0)
        norm = np.trapezoid(P, zgrid, axis=1)
        return P / np.where(norm > 0, norm, 1.0)[:, None]

    def quantile(self, q, zgrid=None):
        zgrid = np.linspace(-0.2 if self.zmin is None else self.zmin, 4.0, 4001) if zgrid is None else zgrid
        P = self.grid(zgrid)
        C = np.concatenate([np.zeros((self.n, 1)), np.cumsum(0.5 * (P[:, 1:] + P[:, :-1]) * np.diff(zgrid), axis=1)], axis=1)
        C = C / np.maximum(C[:, -1:], 1e-12)
        out = np.empty(self.n)
        qv = float(q)
        idx = (C < qv).sum(axis=1).clip(1, len(zgrid) - 1)
        c0 = C[np.arange(self.n), idx - 1]
        c1 = C[np.arange(self.n), idx]
        t = np.where(c1 > c0, (qv - c0) / np.where(c1 > c0, c1 - c0, 1), 0.0)
        return zgrid[idx - 1] + t * (zgrid[idx] - zgrid[idx - 1])

    def mean_std(self, zgrid=None):
        if zgrid is None and self.kind in ('gauss', 'mixture'):          # analytic (untruncated)
            m = (self.pi * self.mu).sum(axis=1)
            v = (self.pi * (self.sig ** 2 + self.mu ** 2)).sum(axis=1) - m ** 2
            return m, np.sqrt(np.maximum(v, 0))
        if zgrid is None and self.kind == 'split':
            return self.med, 0.5 * (self.slo + self.shi)
        zgrid = np.linspace(-0.2 if self.zmin is None else self.zmin, 4.0, 2001) if zgrid is None else zgrid
        P = self.grid(zgrid)
        m = np.trapezoid(P * zgrid, zgrid, axis=1)
        v = np.trapezoid(P * zgrid ** 2, zgrid, axis=1) - m ** 2
        return m, np.sqrt(np.maximum(v, 0))

    def scaled(self, s):
        """Copy with all widths multiplied by s (gauss / mixture / split)."""
        c = Predictive.__new__(Predictive)
        c.__dict__.update(self.__dict__)
        if self.kind in ('gauss', 'mixture'):
            c.sig = self.sig * s
        elif self.kind == 'split':
            c.slo = self.slo * s
            c.shi = self.shi * s
        else:
            raise ValueError('scaling needs a parametric predictive')
        return c


# ---------------------------------------------------------------------------------------------------------- point accuracy
def delta_z(zphot, zspec):
    return (np.asarray(zphot, float) - np.asarray(zspec, float)) / (1.0 + np.asarray(zspec, float))


def nmad(dz):
    dz = np.asarray(dz, float)
    dz = dz[np.isfinite(dz)]
    return float(1.4826 * np.median(np.abs(dz - np.median(dz)))) if len(dz) else float('nan')


def point_metrics(zphot, zspec, outlier=0.15, weights=None):
    """bias, median bias, sigma_NMAD, outlier fractions, percentiles of |dz|.  `weights` gives weighted versions."""
    dz = delta_z(zphot, zspec)
    ok = np.isfinite(dz)
    dz = dz[ok]
    w = np.ones(len(dz)) if weights is None else np.asarray(weights, float)[ok]
    ad = np.abs(dz)

    def wquant(x, q):
        o = np.argsort(x)
        c = np.cumsum(w[o])
        return float(np.interp(q * c[-1], c, x[o]))

    out = dict(n=int(len(dz)))
    if not len(dz):
        return out
    wm = wquant(dz, 0.5)
    out.update(bias_mean=float(np.average(dz, weights=w)), bias_median=wm, sigma_nmad=float(1.4826 * wquant(np.abs(dz - wm), 0.5)),
               sigma_std=float(math.sqrt(np.average((dz - np.average(dz, weights=w)) ** 2, weights=w))),
               outlier_frac=float(np.average(ad > outlier, weights=w)), outlier_cut=outlier,
               outlier_3sig_frac=float(np.average(ad > 3 * 1.4826 * wquant(np.abs(dz - wm), 0.5), weights=w)),
               absdz_p68=wquant(ad, 0.68), absdz_p90=wquant(ad, 0.90), absdz_p95=wquant(ad, 0.95), absdz_p99=wquant(ad, 0.99))
    return out


def bootstrap_point(zphot, zspec, n=200, seed=1, outlier=0.15):
    rng = np.random.default_rng(seed)
    N = len(zphot)
    keys = ('bias_mean', 'sigma_nmad', 'outlier_frac', 'absdz_p95')
    arr = {k: [] for k in keys}
    for _ in range(n):
        j = rng.integers(0, N, N)
        m = point_metrics(np.asarray(zphot)[j], np.asarray(zspec)[j], outlier)
        for k in keys:
            arr[k].append(m[k])
    return {k + '_err': float(np.std(arr[k], ddof=1)) for k in keys}


def binned_metrics(values, zphot, zspec, edges, pred=None, outlier=0.15, min_n=30):
    """Point (and, with `pred`, PDF) metrics in bins of `values` (e.g. zspec, zphot, magnitude)."""
    out = []
    pit = pred_pit = None
    if pred is not None:
        pit = pred.cdf(zspec)
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (values >= lo) & (values < hi) & np.isfinite(values)
        r = dict(lo=float(lo), hi=float(hi), n=int(s.sum()))
        if s.sum() >= min_n:
            r.update({k: v for k, v in point_metrics(zphot[s], zspec[s], outlier).items() if k in ('bias_mean', 'sigma_nmad', 'outlier_frac', 'absdz_p68', 'absdz_p95')})
            if pit is not None:
                r['pit_ks_p'] = float(stats.kstest(pit[s], 'uniform').pvalue)
                r['cover68'] = float(np.mean(np.abs(pit[s] - 0.5) <= 0.34))
                r['pit_var'] = float(np.var(pit[s]))
        out.append(r)
    return out


# --------------------------------------------------------------------------------------------------------- PDF calibration
LEVELS = (0.5, 0.68, 0.8, 0.9, 0.95, 0.99)


def wilson(k, n, z=1.0):
    if n == 0:
        return 0.0, 1.0
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return max(0.0, (c - h) / d), min(1.0, (c + h) / d)


def crps_gauss_mix(pred, y):
    """Closed-form CRPS of a Gaussian mixture (untruncated) at y."""
    def A(m, s2):
        s = np.sqrt(s2)
        return m * (2 * ndtr(m / s) - 1) + 2 * s * np.exp(-0.5 * (m / s) ** 2) / math.sqrt(2 * math.pi)
    y = np.asarray(y, float)[:, None]
    t1 = (pred.pi * A(y - pred.mu, pred.sig ** 2)).sum(axis=1)
    K = pred.mu.shape[1]
    t2 = np.zeros(pred.n)
    for i in range(K):
        for j in range(K):
            t2 += pred.pi[:, i] * pred.pi[:, j] * A(pred.mu[:, i] - pred.mu[:, j], pred.sig[:, i] ** 2 + pred.sig[:, j] ** 2)[:]
    return t1 - 0.5 * t2


def crps(pred, y, zgrid=None):
    """CRPS per object.  Exact for gauss / mixture without truncation effect (zmin tiny mass), numerical on a grid otherwise."""
    y = np.asarray(y, float)
    if pred.kind in ('gauss', 'mixture') and (pred.zmin is None or np.all(pred._f0() < 1e-3)):
        return crps_gauss_mix(pred, y)
    zg = np.linspace(0.0 if pred.zmin is None else pred.zmin, 4.0, 2001) if zgrid is None else zgrid
    P = pred.grid(zg)
    C = np.concatenate([np.zeros((pred.n, 1)), np.cumsum(0.5 * (P[:, 1:] + P[:, :-1]) * np.diff(zg), axis=1)], axis=1)
    H = (zg[None, :] >= y[:, None]).astype(float)
    return np.trapezoid((C - H) ** 2, zg, axis=1)


def log_score(pred, y, floor=1e-12):
    return -np.log(np.maximum(pred.pdf(np.asarray(y, float)), floor))


def pit_report(pit, nbins=10):
    pit = np.asarray(pit, float)
    pit = pit[np.isfinite(pit)]
    n = len(pit)
    h, e = np.histogram(pit, bins=nbins, range=(0, 1))
    exp = n / nbins
    chi2 = float(((h - exp) ** 2 / exp).sum())
    ks = stats.kstest(pit, 'uniform')
    cvm = stats.cramervonmises(pit, 'uniform')
    var = float(np.var(pit))
    mean = float(pit.mean())
    shape = 'consistent with uniform'
    if ks.pvalue < 0.01 or cvm.pvalue < 0.01:
        if abs(mean - 0.5) > 0.03:
            shape = 'tilted (biased: truth systematically %s the PDFs)' % ('above' if mean > 0.5 else 'below')
        elif var > 1 / 12 * 1.1:
            shape = 'U-shaped (overconfident: PDFs too narrow or outliers)'
        elif var < 1 / 12 * 0.9:
            shape = 'hump (underconfident: PDFs too wide)'
        else:
            shape = 'non-uniform'
    cover = []
    for L in LEVELS:
        k = int(np.sum(np.abs(pit - 0.5) <= L / 2))
        lo, hi = wilson(k, n)
        cover.append(dict(level=L, observed=k / n, lo=lo, hi=hi))
    return dict(n=n, hist=[int(x) for x in h], hist_expected=float(exp), chi2_10bins=chi2, chi2_p=float(stats.chi2.sf(chi2, nbins - 1)),
                ks_stat=float(ks.statistic), ks_p=float(ks.pvalue), cvm_stat=float(cvm.statistic), cvm_p=float(cvm.pvalue), mean=mean, var=var, var_uniform=1 / 12,
                frac_extreme=float(np.mean((pit < 0.001) | (pit > 0.999))), shape=shape, coverage=cover)


def optimal_scale(pred, zspec, lo=0.3, hi=3.0):
    """Width scale s minimising the mean log score (gauss / mixture / split); and the score gain."""
    from scipy.optimize import minimize_scalar
    f = lambda s: float(np.mean(log_score(pred.scaled(s), zspec)))
    r = minimize_scalar(f, bounds=(lo, hi), method='bounded', options=dict(xatol=1e-3))
    return dict(scale=float(r.x), logscore_at_scale=float(r.fun), logscore_original=f(1.0))


def recalibrate_pit(pit_train, pit_apply=None):
    """Probability-integral-transform recalibration: map PIT through the empirical CDF of the training PITs (uniform by construction)."""
    t = np.sort(np.asarray(pit_train, float))
    a = np.asarray(pit_train if pit_apply is None else pit_apply, float)
    return (np.searchsorted(t, a, side='right') + 0.5 * (np.searchsorted(t, a, side='left') - np.searchsorted(t, a, side='right'))) / len(t)


def stacked_nz(pred, zspec, zgrid, bins):
    """Stacked p(z) of the sample vs the histogram of the spectroscopic redshifts: (z centres, stacked n(z), spec histogram, KS distance of the CDFs)."""
    P = pred.grid(zgrid)
    stack = P.sum(axis=0)
    stack_cdf = np.cumsum(stack) * (zgrid[1] - zgrid[0])
    stack_cdf = stack_cdf / stack_cdf[-1]
    zs = np.sort(np.asarray(zspec, float))
    spec_cdf = np.searchsorted(zs, zgrid, side='right') / len(zs)
    h, e = np.histogram(zspec, bins=bins)
    cen = 0.5 * (e[1:] + e[:-1])
    sn = np.array([np.trapezoid(stack[(zgrid >= a) & (zgrid <= b)], zgrid[(zgrid >= a) & (zgrid <= b)]) for a, b in zip(e[:-1], e[1:])])
    return dict(z=cen, stacked=sn, spec=h.astype(float), ks_distance=float(np.max(np.abs(stack_cdf - spec_cdf))),
                mean_stacked=float(np.trapezoid(stack * zgrid, zgrid) / np.trapezoid(stack, zgrid)), mean_spec=float(np.mean(zspec)))


def calibration_summary(pred, zspec, nboot=0, seed=1):
    pit = pred.cdf(zspec)
    out = dict(pit=pit_report(pit))
    cr = crps(pred, zspec)
    ls = log_score(pred, zspec)
    out['crps_mean'] = float(np.mean(cr))
    out['crps_median'] = float(np.median(cr))
    out['logscore_mean'] = float(np.mean(ls))
    m, sd = pred.mean_std()
    zsc = (np.asarray(zspec) - m) / np.maximum(sd, 1e-6)
    out['zscore'] = dict(mean=float(np.mean(zsc)), std=float(np.std(zsc)), nmad=nmad(zsc), frac_gt3=float(np.mean(np.abs(zsc) > 3)))
    return out, pit, cr, ls, zsc


# ----------------------------------------------------------------------------------------------------- representativeness
def robust_scale(X, ref=None):
    ref = X if ref is None else ref
    med = np.nanmedian(ref, axis=0)
    mad = 1.4826 * np.nanmedian(np.abs(ref - med), axis=0)
    mad = np.where(mad > 0, mad, np.nanstd(ref, axis=0))
    mad = np.where(mad > 0, mad, 1.0)
    return (X - med) / mad


def feature_table(Xs, Xt, names):
    """1-D comparison spectroscopic (s) vs target (t): KS, standardised mean difference, quantiles."""
    rows = []
    for j, nme in enumerate(names):
        a = Xs[:, j][np.isfinite(Xs[:, j])]
        b = Xt[:, j][np.isfinite(Xt[:, j])]
        ks = stats.ks_2samp(a, b)
        sd = math.sqrt(0.5 * (a.var() + b.var())) or 1.0
        rows.append(dict(feature=nme, n_spec=int(len(a)), n_target=int(len(b)), ks=float(ks.statistic), ks_p=float(ks.pvalue), std_mean_diff=float((a.mean() - b.mean()) / sd),
                         median_spec=float(np.median(a)), median_target=float(np.median(b)), p5_spec=float(np.percentile(a, 5)), p95_spec=float(np.percentile(a, 95)),
                         p5_target=float(np.percentile(b, 5)), p95_target=float(np.percentile(b, 95)),
                         frac_target_outside_spec_range=float(np.mean((b < a.min()) | (b > a.max())))))
    return rows


def knn_classifier_test(Xs, Xt, k=20, nperm=200, seed=1, max_n=20000):
    """Two-sample test in feature space: leave-one-out k-NN probability of 'spec' among the k neighbours in the pooled sample; AUC (0.5 = indistinguishable)
    with a permutation p-value.  Samples are balanced by sub-sampling to the smaller one (max_n each)."""
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(seed)
    n = min(len(Xs), len(Xt), max_n)
    A = Xs[rng.choice(len(Xs), n, replace=False)]
    B = Xt[rng.choice(len(Xt), n, replace=False)]
    X = np.vstack([A, B])
    y = np.r_[np.ones(n), np.zeros(n)]
    nb = cKDTree(X).query(X, k=k + 1)[1][:, 1:]

    def auc(lab):
        score = lab[nb].mean(axis=1)
        r = stats.rankdata(score)
        n1 = lab.sum()
        n0 = len(lab) - n1
        return float((r[lab == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))
    a0 = auc(y)
    perm = np.array([auc(rng.permutation(y)) for _ in range(nperm)])
    p = float((1 + np.sum(perm >= a0)) / (nperm + 1))
    return dict(auc=a0, p=p, perm_mean=float(perm.mean()), perm_std=float(perm.std()), n_per_sample=int(n), k=k)


def knn_weights(Xs, Xt, k=5):
    """Lima et al. (2008) weights of the spectroscopic objects: w = (k / N_target) / (n_spec(<R) / N_spec), R = distance to the k-th target neighbour.
    Weights are normalised to mean 1."""
    from scipy.spatial import cKDTree
    tt = cKDTree(Xt)
    ts = cKDTree(Xs)
    R = tt.query(Xs, k=k)[0]
    R = R[:, -1] if k > 1 else R
    ns = np.array([len(ts.query_ball_point(x, r)) for x, r in zip(Xs, R)], float)
    w = (k / len(Xt)) / (ns / len(Xs))
    return w / w.mean()


def effective_fraction(w):
    w = np.asarray(w, float)
    return float(w.sum() ** 2 / (len(w) * (w ** 2).sum()))


def coverage_distance(Xs, Xt, q=95.0):
    """Distance of every target object to its nearest spec object, in units of the q-th percentile of the leave-one-out nearest-neighbour distances within the
    spec sample; > 1 means 'outside the support of the spec sample'.  Returns (ratio per target object, fraction outside, expected fraction 1 - q/100)."""
    from scipy.spatial import cKDTree
    ts = cKDTree(Xs)
    dss = ts.query(Xs, k=2)[0][:, 1]
    dts = ts.query(Xt, k=1)[0]
    ref = np.percentile(dss, q)
    ratio = dts / ref
    return ratio, float(np.mean(ratio > 1)), 1 - q / 100.0
