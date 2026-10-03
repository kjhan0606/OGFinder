"""Unified photometric error model for aperture photometry: one place for the error terms that several tools used to estimate separately.

Terms (variances add):
  sky        pixel noise summed over the N pixels of the aperture, inflated by the pixel correlation of the image,
             sigma_N = sigma_pix * alpha * N^beta (blank-aperture law, ogfkit.noise; uncorrelated: alpha 1, beta 0.5)
  poisson    max(flux, 0) / gain for the source (electrons; gain 0 = ignore)
  sky_est    error of the local sky level subtracted from all N aperture pixels: N * sigma_mean, sigma_mean = std of the sky estimator
             over the annulus (nsky pixels).  Uncorrelated white noise and a mean estimator give sigma_pix / sqrt(nsky) (the classical DAOPHOT
             term N^2 sigma^2 / nsky); correlated noise and the median / mode estimators are covered by an empirical `sky_factor`
             (kappa = std of the annulus estimate * sqrt(nsky) / sigma_pix, measured on blank annuli, `annulus_factor`)
  sky_sys    systematic sky/background-model error per pixel (counts, e.g. noisemodel mesh error ~0.2 counts), fully correlated over the aperture: N * sky_sys
  contam     uncertainty of the flux of neighbours falling in the aperture (`neighbour_contamination`: sum_j c_ij sigma_fj in quadrature
             plus a relative PSF-model uncertainty of the contamination itself); the contamination itself is a bias that can be subtracted
  apcorr     error of the aperture correction: total = f_ap / frac, var = (sigma_ap / frac)^2 + (total sigma_frac / frac)^2

Pure numpy functions, JSON-serialisable.  Conventions: x, y 0-based pixel coordinates, flux in counts above sky, apertures circular radius r px with
exact / sub-pixel pixel overlap (sep.sum_circle subpix=5).  Everything is a model: validate against blank apertures and injections of the image
(plugins/noisemodel/validation/photerr_validate.py).
"""
import math

import numpy as np

from . import noise as nz

MAG = 2.5 / math.log(10.0)          # 1.0857


class NoiseModel:
    """Pixel rms + correlation law (+ optional empirical sky-estimator factors)."""

    def __init__(self, sigma_pix, alpha=1.0, beta=0.5, ok=False, n_range=None, sky_factors=None, local=None):
        self.sigma_pix = float(sigma_pix)
        self.alpha = float(alpha)
        self.beta = float(beta)
        self.ok = bool(ok)
        self.n_range = n_range
        self.sky_factors = dict(sky_factors or {})
        self.local = local          # dict(ann, mode, alpha, beta): law of the local-sky-subtracted blank-aperture flux, sigma = sigma_pix alpha N^beta

    @classmethod
    def white(cls, sigma_pix):
        return cls(sigma_pix, 1.0, 0.5, False)

    @classmethod
    def from_meta(cls, d):
        """From the dict noisemodel writes (`--meta-out`, key `noisemodel`: sigma_pix, alpha, beta) or from `to_dict`."""
        d = d.get('noisemodel', d)
        return cls(d['sigma_pix'], d.get('alpha', 1.0), d.get('beta', 0.5), d.get('ok', 'alpha' in d), d.get('n_range'), d.get('sky_factors'), d.get('local'))

    @classmethod
    def measure(cls, sub, mask, radii=(1.5, 2, 3, 4, 6, 8, 12), n_aper=500, seed=1):
        """Blank-aperture noise law of a background-subtracted image (mask True = excluded)."""
        an = nz.aperture_noise(sub, mask, radii, n_aper=n_aper, seed=seed)
        law = nz.fit_noise_law(an['N'], an['sigma'], an['sigma_pix'])
        return cls(an['sigma_pix'], law['alpha'], law['beta'], law['ok'], [min(an['N']), max(an['N'])] if an['N'] else None)

    def to_dict(self):
        return dict(sigma_pix=self.sigma_pix, alpha=self.alpha, beta=self.beta, ok=self.ok, n_range=self.n_range, sky_factors=self.sky_factors, local=self.local)

    @property
    def law(self):
        return dict(alpha=self.alpha, beta=self.beta, ok=self.ok)

    def corr(self, N):
        """sigma_N / (sigma_pix sqrt(N)) (1 for white noise); the law is used as fitted, even outside n_range."""
        N = np.asarray(N, float)
        if not self.ok:
            return np.ones_like(N)
        return self.alpha * np.maximum(N, 1e-9) ** (self.beta - 0.5)

    def add_local(self, img, mask, radii, ann, mode='median', n=500, seed=1):
        """Measure the local-sky-subtracted aperture noise law for this annulus / sky mode (stored in `.local`)."""
        t = local_sky_noise(img, mask, radii, ann, mode, n, seed)
        if len(t['N']) < 3:
            return self
        law = nz.fit_noise_law(t['N'], t['sigma'], self.sigma_pix)
        self.local = dict(ann=[float(ann[0]), float(ann[1])], mode=mode, alpha=law['alpha'], beta=law['beta'], ok=law['ok'], N=t['N'], sigma=t['sigma'])
        return self

    def local_sigma(self, N):
        """Std of the local-sky-subtracted aperture sum (includes the sky-estimation error); requires add_local."""
        N = np.asarray(N, float)
        return self.sigma_pix * self.local['alpha'] * np.maximum(N, 1e-9) ** self.local['beta']

    def sigma_sum(self, N, sigma_pix=None):
        sp = self.sigma_pix if sigma_pix is None else sigma_pix
        return sp * np.sqrt(np.asarray(N, float)) * self.corr(N)


def local_sky_noise(img, mask, radii, ann, mode='median', n=500, seed=1):
    """Std of blank-aperture fluxes *after subtracting the local annulus sky* (the photometry actually done), vs aperture area.  Correlated noise with power
    on scales larger than the annulus is removed by the local sky, so this can be well below the global-sky law; it contains the sky-estimation error.
    -> dict(N, sigma, nsky, n_used)"""
    rng = np.random.default_rng(seed)
    out = dict(N=[], sigma=[], nsky=[], n_used=[])
    pos = nz.blank_positions(mask, ann[1] + 1, n, rng)
    if len(pos) < 30:
        return out
    for r in radii:
        ph = aperture_photometry(img, pos, r, ann, mode, mask=mask)
        f = ph['flux']
        f = f[np.isfinite(f)]
        if len(f) < 30:
            continue
        out['N'].append(math.pi * r * r); out['sigma'].append(float(nz.clipped_std(f, 4.0))); out['nsky'].append(float(np.nanmedian(ph['nsky']))); out['n_used'].append(int(len(f)))
    return out


def annulus_factor(sub, mask, r_in, r_out, mode='median', n=400, seed=1, sigma_pix=None):
    """kappa = std over blank positions of the annulus sky estimate * sqrt(nsky) / sigma_pix (mode 'mean' | 'median').  White noise: 1 (mean), ~1.25 (median);
    correlated noise: larger.  Returns (kappa, nsky, n_used)."""
    sep = nz._sep()
    rng = np.random.default_rng(seed)
    pos = nz.blank_positions(mask, r_out, n, rng)
    if len(pos) < 30:
        return float('nan'), 0, len(pos)
    sp = sigma_pix if sigma_pix is not None else nz.robust_std(sub[~mask][:2000000] if (~mask).sum() > 2000000 else sub[~mask])
    iy, ix = np.mgrid[:int(2 * r_out + 3), :int(2 * r_out + 3)]
    h = int(r_out + 1)
    rr = np.hypot(ix - h, iy - h)
    ann = (rr >= r_in) & (rr <= r_out)
    nsky = int(ann.sum())
    est = []
    for x, y in pos:
        cx, cy = int(round(x)), int(round(y))
        cut = sub[cy - h:cy + h + 1, cx - h:cx + h + 1]
        if cut.shape != ann.shape:
            continue
        v = cut[ann & ~mask[cy - h:cy + h + 1, cx - h:cx + h + 1]]
        if len(v) < 10:
            continue
        est.append(np.mean(v) if mode == 'mean' else np.median(v))
    if len(est) < 30:
        return float('nan'), nsky, len(est)
    return float(nz.clipped_std(np.array(est), 4.0) * math.sqrt(nsky) / sp), nsky, len(est)


def sky_est_sigma(nsky, sigma_pix, model=None, sky_factor=None):
    """Std of the sky level estimated from an annulus of `nsky` pixels."""
    nsky = np.maximum(np.asarray(nsky, float), 1.0)
    if sky_factor is not None:
        k = sky_factor
    elif model is not None and model.ok:
        k = model.corr(nsky)                                   # compact-region law as an upper estimate of an annulus
    else:
        k = 1.0
    return sigma_pix * k / np.sqrt(nsky)


def aperture_error(flux, N, sigma_pix, model=None, gain=None, nsky=0, sky_factor=None, sky_sys=0.0, contam_err=0.0, frac=1.0, frac_err=0.0, local=False):
    """Error budget of an aperture flux (counts above sky).  Arguments broadcast.  nsky = 0: no sky-estimation term (global / model sky).
    -> dict(sky, poisson, sky_est, sky_sys, contam, ap [error of the aperture flux], total [error of flux/frac incl. aperture-correction error],
            flux_total, var_ap)"""
    flux = np.asarray(flux, float)
    N = np.asarray(N, float)
    sp = np.asarray(sigma_pix, float)
    if local and model is not None and model.local and model.local.get('ok'):
        sky = model.local_sigma(N) * sp / model.sigma_pix
        nsky = 0
    elif model is not None and model.ok:
        sky = sp * np.sqrt(N) * model.corr(N)
    else:
        sky = sp * np.sqrt(N)
    pois = np.sqrt(np.clip(flux, 0, None) / gain) if gain else np.zeros_like(sky)
    nsky_a = np.asarray(nsky, float)
    se = np.where(nsky_a > 0, N * sky_est_sigma(nsky_a, sp, model, sky_factor), 0.0)
    ss = N * np.asarray(sky_sys, float)
    ce = np.asarray(contam_err, float)
    var_ap = sky ** 2 + pois ** 2 + se ** 2 + ss ** 2 + ce ** 2
    ap = np.sqrt(var_ap)
    frac = np.asarray(frac, float)
    ftot = flux / frac
    tot = np.sqrt(var_ap / frac ** 2 + (ftot * np.asarray(frac_err, float) / frac) ** 2)
    return dict(sky=sky, poisson=pois, sky_est=se, sky_sys=ss + 0 * sky, contam=ce + 0 * sky, ap=ap, total=tot, flux_total=ftot, var_ap=var_ap)


def mag_error(flux, err):
    flux = np.asarray(flux, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(flux > 0, MAG * np.asarray(err, float) / flux, np.nan)


# ----------------------------------------------------------------------------------------------------------- PSF enclosed fractions
def _psf_stamp(psf, x, y, dx, dy):
    if callable(psf):
        return np.asarray(psf(x, y, dx, dy), float)
    p = np.asarray(psf, float)
    if abs(dx) < 1e-9 and abs(dy) < 1e-9:
        return p
    from scipy.ndimage import shift as ndshift
    return ndshift(p, (dy, dx), order=3, mode='constant')


def enclosed(stamp, cx, cy, r, r_in=0.0, r_out=None):
    """Flux of a stamp inside a circle (r) or annulus (r_in..r_out) centred at stamp pixel (cx, cy) (0-based, may be fractional)."""
    sep = nz._sep()
    s = np.ascontiguousarray(stamp, np.float64)
    if r_out is None:
        v, _, _ = sep.sum_circle(s, np.array([cx]), np.array([cy]), r, subpix=5)
    else:
        v, _, _ = sep.sum_circann(s, np.array([cx]), np.array([cy]), r_in, r_out, subpix=5)
    return float(v[0])


def aperture_fraction(psf, r, sky_annulus=None, x=0.0, y=0.0):
    """Net fraction of a point source's flux measured in an aperture of radius r when the sky is subtracted using an annulus (r_in, r_out) that
    contains the source's own wings: frac = f_ap - N / A_ann * f_ann (sky_annulus None: plain enclosed fraction).  The stamp is the unit-sum PSF
    (array, centred on the middle pixel, or callable(x, y, dx, dy)).  Flux outside the stamp is lost: use an extended PSF (psfext) for large radii."""
    st = _psf_stamp(psf, x, y, 0.0, 0.0)
    h = (st.shape[0] - 1) / 2.0
    f = enclosed(st, h, h, r)
    if sky_annulus is None:
        return f
    r_in, r_out = sky_annulus
    a = enclosed(st, h, h, 0, r_in, r_out)
    A = math.pi * (r_out ** 2 - r_in ** 2)
    return f - math.pi * r * r / A * a


def aperture_fraction_scatter(psf_model, xy, r, sky_annulus=None):
    """Mean and rms of the net aperture fraction of a spatially varying PSFModel (`.stamp(x, y)`) over positions xy -> (mean, rms)."""
    fr = np.array([aperture_fraction(lambda x, y, dx, dy: psf_model.stamp(x, y, dx, dy), r, sky_annulus, x, y) for x, y in xy])
    return float(fr.mean()), float(fr.std())


def neighbour_contamination(xy, flux, r, psf, sky_annulus=None, flux_err=None, psf_rel_err=0.1, targets=None, rmax=None):
    """Flux that the PSF wings of the *other* sources put into each target aperture (and, with sky_annulus, the compensating sky over-subtraction):
    contam_i = sum_{j != i} flux_j * [f_ap(j -> i) - N_i / A_ann * f_ann(j -> i)].  The neighbour flux must be the total (not aperture) flux.
    -> dict(contam [bias in counts, positive = target too bright], err [quadrature of c_ij sigma_fj and psf_rel_err * |contam|], n_neigh, dominant_fraction)"""
    from scipy.spatial import cKDTree
    xy = np.asarray(xy, float).reshape(-1, 2)
    flux = np.asarray(flux, float)
    fe = np.zeros_like(flux) if flux_err is None else np.asarray(flux_err, float)
    nt = len(xy)
    idx = range(nt) if targets is None else targets
    st0 = _psf_stamp(psf, xy[0, 0], xy[0, 1], 0, 0) if nt else np.zeros((3, 3))
    half = (st0.shape[0] - 1) / 2.0
    rmax = (half + r) if rmax is None else rmax
    tree = cKDTree(xy)
    out = dict(contam=np.zeros(nt), err=np.zeros(nt), n_neigh=np.zeros(nt, int), max_term=np.zeros(nt))
    for i in idx:
        var = 0.0
        tot = 0.0
        mx = 0.0
        k = 0
        for j in tree.query_ball_point(xy[i], rmax):
            if j == i or not np.isfinite(flux[j]):
                continue
            ix, iy = int(round(xy[j, 0])), int(round(xy[j, 1]))
            st = _psf_stamp(psf, xy[j, 0], xy[j, 1], xy[j, 0] - ix, xy[j, 1] - iy)
            h = (st.shape[0] - 1) / 2.0
            cx, cy = xy[i, 0] - ix + h, xy[i, 1] - iy + h
            c = enclosed(st, cx, cy, r)
            if sky_annulus is not None:
                r_in, r_out = sky_annulus
                c -= math.pi * r * r / (math.pi * (r_out ** 2 - r_in ** 2)) * enclosed(st, cx, cy, 0, r_in, r_out)
            tot += flux[j] * c
            var += (c * fe[j]) ** 2
            mx = max(mx, abs(flux[j] * c))
            k += 1
        out['contam'][i] = tot
        out['err'][i] = math.sqrt(var + (psf_rel_err * tot) ** 2)
        out['n_neigh'][i] = k
        out['max_term'][i] = mx
    return out


# ----------------------------------------------------------------------------------------------------------- measuring with it
def aperture_photometry(data, xy, r, sky_annulus=(12.0, 18.0), sky_mode='median', model=None, gain=None, bkg=None, mask=None):
    """Circular-aperture photometry with local annulus sky and the full error budget (without neighbour / aperture-correction terms, add those
    with `aperture_error` arguments).  data: image; xy: 0-based; mask True = bad (excluded from the sky and aperture).
    -> dict(flux, sky, nsky, N, sigma_pix [local], err dict from aperture_error)"""
    from photutils.aperture import CircularAperture, CircularAnnulus
    d = np.asarray(data, float)
    xy = np.asarray(xy, float).reshape(-1, 2)
    bad = np.zeros(d.shape, bool) if mask is None else np.asarray(mask, bool)
    dd = np.where(bad | ~np.isfinite(d), np.nan, d)
    ann = CircularAnnulus(xy, *sky_annulus).to_mask(method='center')
    ap = CircularAperture(xy, r).to_mask(method='exact')
    ann = ann if isinstance(ann, list) else [ann]
    ap = ap if isinstance(ap, list) else [ap]
    n = len(xy)
    flux = np.full(n, np.nan); sky = np.full(n, np.nan); nsky = np.zeros(n); N = np.full(n, np.nan); sig = np.full(n, np.nan)
    for i in range(n):
        c = ann[i].cutout(dd, fill_value=np.nan)
        if c is None:
            continue
        v = c[(ann[i].data > 0) & np.isfinite(c)]
        if len(v) < 10:
            continue
        if sky_mode == 'mean':
            s = float(np.mean(v))
        elif sky_mode == 'mode':
            m = float(np.median(v)); mu = float(np.mean(v)); s = 3 * m - 2 * mu
        else:
            s = float(np.median(v))
        sig[i] = nz.robust_std(v)
        a = ap[i].cutout(dd, fill_value=np.nan)
        if a is None:
            continue
        w = ap[i].data
        ok = np.isfinite(a)
        N[i] = float((w * ok).sum())
        flux[i] = float(((a - s) * w)[ok].sum())
        sky[i] = s
        nsky[i] = len(v)
    sp = np.where(np.isfinite(sig), sig, model.sigma_pix if model else np.nan)
    return dict(flux=flux, sky=sky, nsky=nsky, N=N, sigma_pix=sp)
