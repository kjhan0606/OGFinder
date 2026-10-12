"""Cluster weak lensing: tangential shear, excess surface density, and an NFW mass.

Coordinates. Positions are arcsec on a local tangent plane, x to the right and y up.
For a sky catalogue the plane is x = -(RA - RA_c) cos(Dec_c) and y = Dec - Dec_c, in arcsec.
Shear components g1, g2 are in that same plane: g1 > 0 stretches along +x.
The tangential and cross components about the cluster centre (angle phi from +x toward +y) are

    g_t = -g1 cos 2phi - g2 sin 2phi
    g_x =  g1 sin 2phi - g2 cos 2phi

so a positive g_t is a tangential alignment. A multiplicative bias m and additive terms
c1, c2 are removed first: g <- (g - c) / (1 + m).

The mass model is a spherical NFW halo (Bartelmann 1996; Wright & Brainerd 2000) defined
by M_200c and c_200c. The reduced shear of each source is g = gamma / (1 - kappa), with
gamma = DeltaSigma / Sigma_crit and kappa = Sigma / Sigma_crit. Sigma_crit uses a flat
LambdaCDM angular-diameter distance. Optional pieces, all off or on by the caller:

* concentration from Duffy et al. 2008, Table 1, NFW, full sample, Delta=200 rho_crit,
  z = 0-2: c = 5.71 (M h / 2e12 Msun)^{-0.084} (1+z)^{-0.47}, with h = H0/100.
  That fit was made on WMAP5 simulations.
* a two-halo term. The library default is a linear Eisenstein & Hu 1998 zero-baryon
  P(k). The menu default is Takahashi et al. 2012 Halofit (w = -1) on that linear
  spectrum, with a Tinker et al. 2010 bias. Not a Boltzmann code. Full Halofit
  includes the one-halo power of the field, so the two-halo term can double-count
  small-scale matter.
* miscentring: a centred fraction plus a Rayleigh offset (Johnston et al. 2007),
  azimuthally averaged. Both can be fitted with Gaussian priors supplied by the caller.
* a boost from a random catalogue of the same selection, normalised to 1 in the
  outermost radial bin. Members dilute the shear, so the correction multiplies g_t.

Magnification, when the caller supplies a count slope s = d log10 N(<m) / dm, multiplies
the model shear by 1 + (5s - 2) kappa. s = 0 leaves it unchanged. s = 0.4 is a null test.
The fit covariance is the shape noise unless the caller adds the Gaussian large-scale-structure
variance of the annulus-averaged tangential shear, or an azimuthal jackknife that replaces both.
A selection response for a resolution cut and an absolute-ellipticity cut is separate from the
per-object shear response. An internal m and c from a Gaussian or Sersic drawing is not a
survey calibration and is not applied to the mass unless the caller copies m across.

A stack of cluster profiles lives in stack_shear_profiles. Each cluster is a shear
catalogue: positions, two shear components, a redshift, and a weight. The function
bins each catalogue with select_sources, averages those profiles, and can fit the
average. It does not accept an image, and it does not measure shapes on a coadd.
""",
import math

import numpy as np
from scipy.integrate import quad
from scipy.ndimage import map_coordinates
from scipy.optimize import least_squares
from scipy.signal import fftconvolve
from scipy.special import j0, j1

def _trapz(y, x, axis=0):
    fn = getattr(np, 'trapezoid', None) or np.trapz
    return fn(y, x, axis=axis)


C_KMS = 299792.458
ARCSEC = math.pi / 648000.0
G_MPC = 4.30091e-9  # Mpc (km/s)^2 / Msun
TCMB = 2.725
DELTA_C = 1.686

# Duffy et al. 2008, Table 1, NFW, full sample, 200 rho_crit, redshift range 0-2.
DUFFY_A = 5.71
DUFFY_B = -0.084
DUFFY_C = -0.47
DUFFY_PIVOT = 2.0e12  # h^-1 Msun

METACAL_SCOPE = (
    'measure_stamps uses one round Gaussian PSF and a per-object diagonal response. '
    'A spatially varying PSF, an ensemble response beyond the resolution and ellipticity cuts, '
    'noise correlation, and detection or blending bias remain survey-calibration tasks. '
    'Completing this step does not complete a survey calibration.')
DELTASIGMA_NOTE = (
    'deltasigma is the mean tangential shear divided by the mean of Sigma_crit^{-1}, '
    'written in Msun/pc^2 so the column matches a single-cluster profile. '
    'The mass model uses the reduced shear of each source.')
DENSITY_NOTE = (
    'A change in source number density combines member contamination, obscuration, '
    'detection loss, and magnification. The boost factor does not separate those contributions.')
POPULATION_NOTE = (
    'The ICL weights and the lensing weights differ, so the two averages do not represent the same population. '
    'Shirasaki & Takada 2018 MNRAS 478, 4277 combine profiles drawn from one sample.')
WEIGHT_NOTE = (
    'sigma_e/sqrt(sum w) applies when w is an inverse-variance weight and weight 1 has variance sigma_e^2. '
    'These weights were not certified on that scale, so the formula is omitted.')
LSS_SKIP_NOTE = (
    'The large-scale-structure term was left off. It is added only for one source plane, '
    'and it is not a survey covariance.')
LSS_USED_NOTE = (
    'The covariance is the shape noise plus the Gaussian large-scale-structure variance of the '
    'annulus-averaged tangential shear, for one source plane. It is not a survey covariance.')


def rho_crit(z, H0, Om):
    """Critical density at redshift z, Msun/Mpc^3."""
    H = H0 * math.sqrt(Om * (1 + z) ** 3 + (1 - Om))
    return 3 * H * H / (8 * math.pi * G_MPC)


def _chi_grid(H0, Om, zmax=8.0, n=512):
    zs = np.linspace(0.0, zmax, n)
    Ol = 1.0 - Om
    E = np.sqrt(Om * (1 + zs) ** 3 + Ol)
    # cumulative trapezoid of dz/E, times c/H0
    dz = np.diff(zs)
    integ = np.concatenate([[0.0], np.cumsum(0.5 * (1 / E[1:] + 1 / E[:-1]) * dz)])
    return zs, C_KMS / H0 * integ


class Cosmology:
    """Flat LambdaCDM distances and a cached linear spectrum. Numpy and scipy only."""

    def __init__(self, H0=70.0, Om=0.3, ns=0.96, sigma8=0.8):
        self.H0 = float(H0)
        self.Om = float(Om)
        self.ns = float(ns)
        self.sigma8 = float(sigma8)
        self.h = self.H0 / 100.0
        self._zs, self._chi = _chi_grid(self.H0, self.Om)
        self._growth = {}
        self._field = None

    def comoving(self, z):
        z = np.asarray(z, float)
        return np.interp(z, self._zs, self._chi, left=0.0, right=np.nan)

    def angular_diameter(self, z):
        z = np.asarray(z, float)
        return self.comoving(z) / (1.0 + z)

    def angular_diameter_z1z2(self, z1, z2):
        """D_A(z1, z2) in Mpc. Negative redshifts or z2 <= z1 give nan."""
        z1 = np.asarray(z1, float)
        z2 = np.asarray(z2, float)
        chi1 = self.comoving(z1)
        chi2 = self.comoving(z2)
        out = (chi2 - chi1) / (1.0 + z2)
        return np.where((z2 > z1) & (z1 >= 0) & np.isfinite(out), out, np.nan)

    def sigma_crit(self, z_l, z_s):
        """Critical surface density, Msun/Mpc^2. nan when the source is not behind the lens."""
        Dl = self.angular_diameter(z_l)
        Ds = self.angular_diameter(z_s)
        Dls = self.angular_diameter_z1z2(z_l, z_s)
        geom = Ds / (Dl * Dls)
        pref = C_KMS ** 2 / (4 * math.pi * G_MPC)
        return pref * geom

    def growth(self, z):
        """Linear growth normalised to 1 at z=0."""
        z = float(z)
        if z in self._growth:
            return self._growth[z]
        Om = self.Om
        Ol = 1.0 - Om

        def E(a):
            return math.sqrt(Om / a ** 3 + Ol)

        def integrand(a):
            ea = E(a)
            return 1.0 / (a ** 3 * ea ** 3)

        def raw(a):
            I, _ = quad(integrand, 1e-4, a, epsabs=1e-6, limit=200)
            return 2.5 * Om * E(a) * I

        d = raw(1.0 / (1.0 + z)) / raw(1.0)
        self._growth[z] = float(d)
        return self._growth[z]

    def duffy_c(self, M200, z):
        """c_200c from Duffy et al. 2008 Table 1 (NFW, full, 200c, z=0-2)."""
        Mh = float(M200) * self.h
        return DUFFY_A * (Mh / DUFFY_PIVOT) ** DUFFY_B * (1.0 + z) ** DUFFY_C

    def linear_field(self):
        if self._field is None:
            self._field = LinearField(self)
        return self._field


def transfer_eh98(k_mpc, H0, Om):
    """Eisenstein & Hu 1998 zero-baryon transfer function, their equation (29).

    k_mpc is in 1/Mpc. T -> 1 as k -> 0. Theta uses T_CMB = 2.725 K.
    """
    k = np.asarray(k_mpc, float)
    h = H0 / 100.0
    theta = TCMB / 2.7
    gamma = Om * h
    q = (k / h) * theta ** 2 / gamma
    L0 = np.log(2.0 * math.e + 1.8 * q)
    C0 = 14.2 + 731.0 / (1.0 + 62.5 * q)
    T = L0 / (L0 + C0 * q * q)
    return np.where(k <= 0, 1.0, T)


def _tophat(x):
    x = np.maximum(np.asarray(x, float), 1e-6)
    return 3.0 * (np.sin(x) - x * np.cos(x)) / x ** 3


class LinearField:
    """z=0 linear P(k), xi(r), and sigma(R). Distances in Mpc, k in 1/Mpc."""

    def __init__(self, cosmo):
        self.cosmo = cosmo
        k = np.geomspace(1e-4, 30.0, 768)
        T = transfer_eh98(k, cosmo.H0, cosmo.Om)
        shape = k ** cosmo.ns * T ** 2
        R8 = 8.0 / cosmo.h
        sig2 = self._sigma2(k, shape, R8)
        amp = (cosmo.sigma8 ** 2) / sig2
        self._amp = float(amp)
        self.k = k
        self.P = amp * shape
        self._halofit_cache = {}
        r = np.geomspace(0.02, 400.0, 240)
        kr = k[:, None] * r[None, :]
        integ = (k ** 2 * self.P)[:, None] * np.sinc(kr / math.pi)
        self.r = r
        self._xi = _trapz(integ, k, axis=0) / (2 * math.pi ** 2)
        r8 = 8.0 / cosmo.h
        self._sig_r = np.unique(np.sort(np.append(np.geomspace(0.05, 80.0, 80), r8)))
        self._sig = np.array([math.sqrt(max(self._sigma2(k, self.P, R), 0.0)) for R in self._sig_r])

    @staticmethod
    def _sigma2(k, P, R):
        W = _tophat(k * R)
        return float(_trapz(k ** 2 * P * W ** 2, k) / (2 * math.pi ** 2))

    def xi(self, r):
        return np.interp(np.asarray(r, float), self.r, self._xi, left=self._xi[0], right=0.0)

    def sigma(self, R):
        return np.interp(np.asarray(R, float), self._sig_r, self._sig, left=self._sig[0], right=self._sig[-1])

    def sigma8(self):
        """sigma(R=8 Mpc/h) from the same integral that fixed the amplitude."""
        return math.sqrt(max(self._sigma2(self.k, self.P, 8.0 / self.cosmo.h), 0.0))

    def linear_pk(self, k, z):
        """Eisenstein & Hu linear P(k, z), Mpc^3. Growth is D(z)^2 and is not in Halofit twice."""
        k = np.asarray(k, float)
        T = transfer_eh98(k, self.cosmo.H0, self.cosmo.Om)
        D2 = self.cosmo.growth(float(z)) ** 2
        return self._amp * np.power(np.maximum(k, 0.0), self.cosmo.ns) * T ** 2 * D2

    def _gauss_sigma2(self, k, P, R):
        """Gaussian-filtered variance, Smith et al. 2003: int Delta^2(k) exp(-k^2 R^2) d ln k."""
        integ = k ** 2 * P * np.exp(-(k * R) ** 2)
        return float(_trapz(integ, k) / (2 * math.pi ** 2))

    def halofit_params(self, z):
        """Takahashi et al. 2012 appendix, w = -1, so the Omega_w (1+w) terms are zero.

        Returns k_sigma, n_eff and C. n_eff = -3 - d ln sigma^2 / d ln R at sigma = 1.
        C = - d^2 ln sigma^2 / d ln R^2 at the same point.
        """
        z = float(z)
        cached = self._halofit_cache.get(('par', round(z, 6)))
        if cached is not None:
            return cached
        k = np.geomspace(1e-4, 200.0, 1024)
        P = self.linear_pk(k, z)

        def sig2(R):
            return self._gauss_sigma2(k, P, R)

        lo, hi = 1e-3, 200.0
        slo, shi = sig2(lo), sig2(hi)
        if not (slo > 1.0 > shi):
            raise ValueError('Halofit could not bracket sigma(R) = 1')
        for _ in range(60):
            mid = math.sqrt(lo * hi)
            if sig2(mid) > 1.0:
                lo = mid
            else:
                hi = mid
        R = math.sqrt(lo * hi)
        eps = 0.05
        lp = math.log(max(sig2(R * math.exp(eps)), 1e-30))
        lm = math.log(max(sig2(R * math.exp(-eps)), 1e-30))
        d1 = (lp - lm) / (2.0 * eps)
        d2 = (lp + lm) / (eps * eps)
        n_eff = -3.0 - d1
        C = -d2
        out = (1.0 / R, float(n_eff), float(C))
        self._halofit_cache[('par', round(z, 6))] = out
        return out

    def halofit_pk(self, z, field_one_halo=True):
        """Nonlinear P(k) from Takahashi et al. 2012. The array is Mpc^3 on a fixed k grid.

        Growth is already inside the linear spectrum at z, so a later projection must not
        multiply by D(z)^2 again. field_one_halo False keeps only the quasi-linear piece.
        """
        z = float(z)
        key = ('pk', round(z, 6), bool(field_one_halo))
        cached = self._halofit_cache.get(key)
        if cached is not None:
            return cached
        k_sigma, n_eff, Ccurv = self.halofit_params(z)
        ne = n_eff
        ne2, ne3, ne4 = ne * ne, ne ** 3, ne ** 4
        Omz = self._omega_m(z)
        # w = -1: the published Omega_w (1+w) terms vanish.
        log_a = 1.5222 + 2.8553 * ne + 2.3706 * ne2 + 0.9903 * ne3 + 0.2250 * ne4 - 0.6038 * Ccurv
        log_b = -0.5642 + 0.5864 * ne + 0.5716 * ne2 - 1.5474 * Ccurv
        log_c = 0.3698 + 2.0404 * ne + 0.8161 * ne2 + 0.5869 * Ccurv
        gamma_n = 0.1971 - 0.0843 * ne + 0.8460 * Ccurv
        alpha_n = abs(6.0835 + 1.3373 * ne - 0.1959 * ne2 - 5.5274 * Ccurv)
        beta_n = 2.0379 - 0.7354 * ne + 0.3157 * ne2 + 1.2490 * ne3 + 0.3980 * ne4 - 0.1682 * Ccurv
        log_nu = 5.2105 + 3.6902 * ne
        a_n = 10 ** log_a
        b_n = 10 ** log_b
        c_n = 10 ** log_c
        nu_n = 10 ** log_nu
        f1 = Omz ** -0.0307
        f2 = Omz ** -0.0585
        f3 = Omz ** 0.0743
        k = np.geomspace(1e-4, 100.0, 640)
        P_L = self.linear_pk(k, z)
        delta_L = k ** 3 * P_L / (2 * math.pi ** 2)
        y = k / k_sigma
        damp = np.exp(-0.25 * y - 0.125 * y * y)
        delta_Q = delta_L * np.power(1.0 + delta_L, beta_n) / (1.0 + alpha_n * delta_L) * damp
        delta = delta_Q
        if field_one_halo:
            y_safe = np.maximum(y, 1e-8)
            numer = a_n * np.power(y_safe, 3.0 * f1)
            denom = 1.0 + b_n * np.power(y_safe, f2) + np.power(c_n * f3 * y_safe, 3.0 - gamma_n)
            delta_H = numer / denom
            # mu_n = 0 in Takahashi et al. 2012, so only nu_n / y^2 remains.
            delta_H = delta_H / (1.0 + nu_n / y_safe ** 2)
            delta = delta_Q + delta_H
        P = np.where(k > 0, delta * (2 * math.pi ** 2) / k ** 3, 0.0)
        P = np.where(np.isfinite(P) & (P > 0), P, 0.0)
        out = (k, P, dict(k_sigma=k_sigma, n_eff=n_eff, C=Ccurv))
        self._halofit_cache[key] = out
        return out

    def _omega_m(self, z):
        Ez2 = self.cosmo.Om * (1 + z) ** 3 + (1 - self.cosmo.Om)
        return self.cosmo.Om * (1 + z) ** 3 / Ez2


def tinker_bias(nu, Delta):
    """Tinker et al. 2010, equation (6) and Table 2. Delta is the overdensity in units of the mean density.

    nu = delta_c / sigma(M). The table is calibrated for 200 <= Delta <= 3200; values outside that
    range are clamped. y = log10(Delta).
    """
    D = float(np.clip(Delta, 200.0, 3200.0))
    y = math.log10(D)
    y4 = (4.0 / y) ** 4
    damp = math.exp(-y4)
    A = 1.0 + 0.24 * y * damp
    a = 0.44 * y - 0.88
    B = 0.183
    b = 1.5
    C = 0.019 + 0.107 * y + 0.19 * damp
    c = 2.4
    nu = np.asarray(nu, float)
    return 1.0 - A * nu ** a / (nu ** a + DELTA_C ** a) + B * nu ** b + C * nu ** c


def halo_bias(cosmo, M200, z):
    """Tinker bias of a halo of mass M200 (Msun, 200 rho_crit) at redshift z."""
    field = cosmo.linear_field()
    rho_m0 = cosmo.Om * rho_crit(0.0, cosmo.H0, cosmo.Om)
    R = (3.0 * M200 / (4.0 * math.pi * rho_m0)) ** (1.0 / 3.0)
    sig = float(field.sigma(R)) * cosmo.growth(z)
    nu = DELTA_C / max(sig, 1e-4)
    Ez2 = cosmo.Om * (1 + z) ** 3 + (1 - cosmo.Om)
    Omz = cosmo.Om * (1 + z) ** 3 / Ez2
    Delta_mean = 200.0 / Omz
    return float(tinker_bias(nu, Delta_mean))


def nfw_scales(M200, c, z, cosmo):
    """Return (rho_s, r_s, R200) for an NFW halo. M200 is M_200c in Msun, radii in Mpc."""
    rc = rho_crit(z, cosmo.H0, cosmo.Om)
    R200 = (3.0 * M200 / (4.0 * math.pi * 200.0 * rc)) ** (1.0 / 3.0)
    rs = R200 / c
    dc = (200.0 / 3.0) * c ** 3 / (math.log(1.0 + c) - c / (1.0 + c))
    return dc * rc, rs, R200


def _nfw_fh(x):
    """Dimensionless projected NFW factors f and h. Sigma = 2 rho_s r_s f, Sigma_bar = 4 rho_s r_s h / x^2."""
    x = np.maximum(np.asarray(x, float), 1e-6)
    f = np.empty_like(x)
    h = np.empty_like(x)
    lo = x < 0.999
    hi = x > 1.001
    mid = ~lo & ~hi
    if np.any(lo):
        xl = x[lo]
        sl = np.sqrt(1.0 - xl * xl)
        F = np.arccosh(np.clip(1.0 / xl, 1.0, None)) / sl
        f[lo] = (1.0 - F) / (xl * xl - 1.0)
        h[lo] = np.log(xl / 2.0) + F
    if np.any(hi):
        xh = x[hi]
        sh = np.sqrt(xh * xh - 1.0)
        F = np.arccos(np.clip(1.0 / xh, -1.0, 1.0)) / sh
        f[hi] = (1.0 - F) / (xh * xh - 1.0)
        h[hi] = np.log(xh / 2.0) + F
    if np.any(mid):
        f[mid] = 1.0 / 3.0
        h[mid] = math.log(0.5) + 1.0
    return f, h


def nfw_sigma(R, rho_s, r_s):
    """Projected surface density and its excess, Msun/Mpc^2, at proper radius R (Mpc)."""
    R = np.asarray(R, float)
    x = np.maximum(R / r_s, 1e-6)
    f, h = _nfw_fh(x)
    Sigma = 2.0 * rho_s * r_s * f
    Sigma_bar = 4.0 * rho_s * r_s * h / (x * x)
    return Sigma, Sigma_bar - Sigma


def mass_enclosed(R, rho_s, r_s):
    """NFW mass inside a 3-D radius R (Mpc), Msun."""
    x = np.maximum(np.asarray(R, float) / r_s, 1e-8)
    return 4.0 * math.pi * rho_s * r_s ** 3 * (np.log(1.0 + x) - x / (1.0 + x))


def _mean_interior(R, Sigma):
    """Azimuthally averaged excess: (2/R^2) int_0^R Sigma R' dR' - Sigma."""
    R = np.asarray(R, float)
    Sigma = np.asarray(Sigma, float)
    order = np.argsort(R)
    Rs, Ss = R[order], Sigma[order]
    Rp = np.concatenate([[0.0], Rs])
    Sp = np.concatenate([[Ss[0]], Ss])
    integ = Sp * Rp
    cum = np.cumsum(0.5 * (integ[1:] + integ[:-1]) * np.diff(Rp))
    bar = 2.0 * cum / np.maximum(Rs, 1e-12) ** 2
    out = np.empty_like(R)
    out[order] = bar - Ss
    return out


def miscenter_sigma(R, sigma_of_R, sigma_off, n_ang=24, n_off=12):
    """Surface density about the observed centre for a Rayleigh offset of scale sigma_off (Mpc)."""
    R = np.asarray(R, float)
    if sigma_off <= 0:
        return sigma_of_R(R)
    u = (np.arange(n_off) + 0.5) / n_off
    roff = sigma_off * np.sqrt(-2.0 * np.log(np.clip(1.0 - u, 1e-8, 1.0)))
    ang = np.linspace(0.0, 2.0 * math.pi, n_ang, endpoint=False)
    acc = np.zeros_like(R, float)
    for ro in roff:
        ring = np.zeros_like(R, float)
        for a in ang:
            Rt = np.sqrt(R * R + ro * ro + 2.0 * R * ro * math.cos(a))
            ring += sigma_of_R(Rt)
        acc += ring / n_ang
    return acc / n_off


def tangential(x, y, g1, g2, xc=0.0, yc=0.0):
    """Tangential and cross shear. See the module docstring for the sign."""
    phi = np.arctan2(np.asarray(y, float) - yc, np.asarray(x, float) - xc)
    c = np.cos(2.0 * phi)
    s = np.sin(2.0 * phi)
    g1 = np.asarray(g1, float)
    g2 = np.asarray(g2, float)
    gt = -g1 * c - g2 * s
    gx = g1 * s - g2 * c
    return gt, gx


def components_from_tangential(phi, gt, gx=0.0):
    """Inverse of tangential() at position angle phi."""
    c = np.cos(2.0 * phi)
    s = np.sin(2.0 * phi)
    g1 = -np.asarray(gt, float) * c + np.asarray(gx, float) * s
    g2 = -np.asarray(gt, float) * s - np.asarray(gx, float) * c
    return g1, g2


def _xi_from_power(k, P, r):
    """xi(r) = int k^2 P(k) sinc(kr) dk / (2 pi^2). k may be logarithmic; the integral uses a linear grid."""
    r = np.asarray(r, float)
    positive = r[r > 0]
    rmin = float(np.min(positive)) if len(positive) else 0.2
    rmax = float(np.max(r)) if len(r) else 300.0
    dk = (2.0 * math.pi / max(rmax, 1.0)) / 8.0
    kmax = min(80.0, max(30.0, 6.0 / max(rmin, 0.05)))
    n = int(min(16000, max(2000, kmax / max(dk, 1e-4))))
    kk = np.linspace(dk, kmax, n)
    logk = np.log(k)
    Pk = np.exp(np.interp(np.log(kk), logk, np.log(np.maximum(P, 1e-40)), left=-80.0, right=-80.0))
    out = np.empty(len(r))
    block = 24
    pref = kk ** 2 * Pk
    for i0 in range(0, len(r), block):
        rr = r[i0:i0 + block]
        kr = kk[:, None] * rr[None, :]
        integ = pref[:, None] * np.sinc(kr / math.pi)
        out[i0:i0 + block] = _trapz(integ, kk, axis=0) / (2 * math.pi ** 2)
    return out


def _matter_projection(cosmo, z, R_grid, power='linear'):
    """Matter Sigma and DeltaSigma (b = 1) at proper radii R_grid, Msun/Mpc^2, at the lens redshift.

    power 'linear' grows the z=0 correlation by D(z)^2. power 'halofit' uses Takahashi et al. 2012
    at the lens redshift, and that spectrum already contains the growth.
    """
    field = cosmo.linear_field()
    rho_m0 = cosmo.Om * rho_crit(0.0, cosmo.H0, cosmo.Om)
    u = np.linspace(0.0, 300.0, 181)
    if power == 'linear':
        grow = cosmo.growth(z) ** 2

        def xi_at(r3):
            return field.xi(r3)
    elif power == 'halofit':
        grow = 1.0
        k_nl, P_nl, _info = field.halofit_pk(z, field_one_halo=True)
        # Radii the line of sight will ask for, so the transform is done once.
        rc = np.asarray(R_grid, float) * (1.0 + z)
        r_need = np.sqrt(rc[:, None] ** 2 + u[None, :] ** 2).ravel()
        r_tab = np.geomspace(max(float(np.min(r_need)) * 0.8, 0.01), 420.0, 160)
        xi_tab = _xi_from_power(k_nl, P_nl, r_tab)

        def xi_at(r3):
            return np.interp(r3, r_tab, xi_tab, left=xi_tab[0], right=0.0)
    else:
        raise ValueError('power must be linear or halofit')
    Sigma = np.empty(len(R_grid))
    for i, Rp in enumerate(R_grid):
        Rc = Rp * (1.0 + z)
        r3 = np.sqrt(Rc * Rc + u * u)
        xi = xi_at(r3)
        # 2 * int_0^umax, the line of sight is symmetric
        I = 2.0 * float(_trapz(xi, u))
        Sigma[i] = rho_m0 * (1.0 + z) ** 2 * grow * I
    return Sigma, _mean_interior(R_grid, Sigma)


class HaloModel:
    """NFW, optional two-halo term and miscentring, evaluated on a radial grid.

    power is 'linear' (the library default) or 'halofit'. The menu asks for Halofit.
    profiles() reads f_cen and sigma_off at call time, so a fit can change them.
    """

    def __init__(self, cosmo, z, R, two_halo=True, miscenter=False, f_cen=1.0, sigma_off=0.1, power='linear'):
        self.cosmo = cosmo
        self.z = float(z)
        self.R = np.asarray(R, float)
        self.two_halo = bool(two_halo)
        self.miscenter = bool(miscenter)
        self.f_cen = float(f_cen)
        self.sigma_off = float(sigma_off)
        self.power = power if self.two_halo else 'linear'
        if self.power not in ('linear', 'halofit'):
            raise ValueError('power must be linear or halofit')
        if self.two_halo:
            self._Sm, self._dSm = _matter_projection(cosmo, self.z, self.R, power=self.power)
        else:
            self._Sm = np.zeros_like(self.R)
            self._dSm = np.zeros_like(self.R)

    def profiles(self, M200, c):
        """Return Sigma, DeltaSigma in Msun/Mpc^2 on self.R."""
        rho_s, r_s, _ = nfw_scales(M200, c, self.z, self.cosmo)
        Sig, dSig = nfw_sigma(self.R, rho_s, r_s)
        if self.miscenter and self.f_cen < 1.0 and self.sigma_off > 0:
            Sig_m = miscenter_sigma(self.R, lambda r, rs=rho_s, rr=r_s: nfw_sigma(r, rs, rr)[0], self.sigma_off)
            dSig_m = _mean_interior(self.R, Sig_m)
            f = self.f_cen
            Sig = f * Sig + (1.0 - f) * Sig_m
            dSig = f * dSig + (1.0 - f) * dSig_m
        if self.two_halo:
            b = halo_bias(self.cosmo, M200, self.z)
            Sig = Sig + b * self._Sm
            dSig = dSig + b * self._dSm
        return Sig, dSig


def reduced_shear(Sigma, DeltaSigma, sigma_crit, reduced=True):
    """Per-source reduced shear. sigma_crit is Msun/Mpc^2."""
    sc = np.asarray(sigma_crit, float)
    kappa = np.asarray(Sigma, float) / sc
    gamma = np.asarray(DeltaSigma, float) / sc
    if not reduced:
        return gamma
    denom = 1.0 - kappa
    g = np.where(denom > 0.05, gamma / denom, np.nan)
    return g


def correct_shear(g1, g2, m=0.0, c1=0.0, c2=0.0):
    fac = 1.0 + float(m)
    if fac <= 0:
        raise ValueError('multiplicative bias leaves a non-positive response')
    return (np.asarray(g1, float) - c1) / fac, (np.asarray(g2, float) - c2) / fac


def effective_source_z(cosmo, z_lens, z, n, z_gap=0.1):
    """Source redshift whose Sigma_crit^{-1} matches the n(z)-weighted mean behind the lens."""
    z = np.asarray(z, float)
    n = np.asarray(n, float)
    m = (z > float(z_lens) + float(z_gap)) & (n > 0) & np.isfinite(z) & np.isfinite(n)
    if not np.any(m):
        raise ValueError('n(z) has no source behind the lens')
    zz, nn = z[m], n[m]
    inv = 1.0 / cosmo.sigma_crit(z_lens, zz)
    target = float(np.sum(nn * inv) / np.sum(nn))
    grid = np.linspace(float(z_lens) + float(z_gap) + 1e-4, min(7.5, float(np.max(zz)) + 0.5), 4000)
    invg = 1.0 / cosmo.sigma_crit(z_lens, grid)
    return float(grid[np.argmin(np.abs(invg - target))])


def boost_factor(R_src, R_rand, edges):
    """n_src/n_rand divided by the same ratio in the outermost bin that has randoms. 1 if there is no random catalogue."""
    if R_rand is None or len(R_rand) == 0:
        return np.ones(len(edges) - 1)
    ns, _ = np.histogram(R_src, edges)
    nr, _ = np.histogram(R_rand, edges)
    with np.errstate(divide='ignore', invalid='ignore'):
        raw = ns / np.maximum(nr, 0)
    raw = np.where(nr > 0, raw, np.nan)
    finite = np.where(np.isfinite(raw) & (raw > 0))[0]
    if len(finite) == 0:
        return np.ones(len(edges) - 1)
    norm = raw[finite[-1]]
    out = raw / norm
    return np.where(np.isfinite(out), out, 1.0)


def bin_edges(rmin, rmax, n):
    return np.geomspace(rmin, rmax, int(n) + 1)


def _weighted_mean(values, w):
    sw = np.sum(w)
    if sw <= 0:
        return np.nan, np.nan, 0.0
    mean = np.sum(w * values) / sw
    neff = sw ** 2 / np.sum(w * w)
    return float(mean), float(neff), float(sw)


def select_sources(x_arcsec, y_arcsec, g1, g2, z, weight, z_lens, cosmo, z_gap=0.1, z_shift=0.0,
                   rmin=0.3, rmax=3.0, gmax=1.5, m=0.0, c1=0.0, c2=0.0, member=None):
    """Sources behind the lens, inside the radial window, with a finite corrected shear.

    x_arcsec, y_arcsec are relative to the cluster centre. Returns a dict of arrays.
    member, when given, is the membership label of each row. Only rows labelled
    background are kept. Satellites, the BCG, stars and unknowns drop out.
    """
    g1c, g2c = correct_shear(g1, g2, m, c1, c2)
    zs = np.asarray(z, float) + float(z_shift)
    R_ang = np.hypot(x_arcsec, y_arcsec) * ARCSEC  # radians
    Dl = float(cosmo.angular_diameter(z_lens))
    R = R_ang * Dl
    sc = cosmo.sigma_crit(z_lens, zs)
    w = np.asarray(weight, float)
    w = np.where(np.isfinite(w) & (w > 0), w, 0.0)
    gt, gx = tangential(x_arcsec, y_arcsec, g1c, g2c)
    phi = np.arctan2(np.asarray(y_arcsec, float), np.asarray(x_arcsec, float))
    ok = (np.isfinite(gt) & np.isfinite(gx) & np.isfinite(sc) & (sc > 0)
          & (zs > z_lens + z_gap) & (R >= rmin) & (R <= rmax)
          & (np.abs(g1c) < gmax) & (np.abs(g2c) < gmax) & (w > 0))
    if member is not None:
        member = np.asarray(member)
        if member.shape != zs.shape:
            raise ValueError('member labels do not match the catalogue')
        ok = ok & (member == 'background')
    return dict(R=R[ok], gt=gt[ok], gx=gx[ok], w=w[ok], sigma_crit=sc[ok],
                z=zs[ok], index=np.where(ok)[0], Dl=Dl, phi=phi[ok])


def _profile_from_selection(sel, edges, sigma_e, random_R=None):
    R, gt, gx, w, sc = sel['R'], sel['gt'], sel['gx'], sel['w'], sel['sigma_crit']
    bins = np.digitize(R, edges) - 1
    B = boost_factor(R, random_R, edges)
    rows = []
    for i in range(len(edges) - 1):
        m = bins == i
        if not np.any(m):
            continue
        # Boost multiplies the observed shear: members dilute g_t.
        gtb = gt[m] * B[i]
        gxb = gx[m] * B[i]
        mt, neff, sw = _weighted_mean(gtb, w[m])
        mx, _, _ = _weighted_mean(gxb, w[m])
        # sigma_e is the rms of one component at weight 1. var(mean) = sigma_e^2 / sum(w).
        err = sigma_e / math.sqrt(sw) if sw > 0 else np.nan
        inv = _weighted_mean(1.0 / sc[m], w[m])[0]
        dS = mt / inv if inv and np.isfinite(inv) and inv > 0 else np.nan
        dSerr = err / inv if inv and np.isfinite(inv) and inv > 0 else np.nan
        Rm = (edges[i] * edges[i + 1]) ** 0.5
        Dl = sel.get('Dl')
        rows.append(dict(bin=i, R=Rm, n=int(m.sum()), neff=neff,
                         gt=mt, gx=mx, gt_err=err, gx_err=err, boost=float(B[i]),
                         deltasigma=dS / 1e12, deltasigma_err=dSerr / 1e12,
                         sigma_crit_inv=inv, sw=sw,
                         R_arcsec=(float(Rm / Dl / ARCSEC) if Dl else None)))
    return rows, bins, B


def _bin_mean(values, bins, rows, weights):
    out = []
    for row in rows:
        m = bins == row['bin']
        gm = values[m]
        ww = weights[m]
        good = np.isfinite(gm) & (ww > 0)
        if not np.any(good):
            out.append(np.nan)
        else:
            out.append(_weighted_mean(gm[good], ww[good])[0])
    return np.array(out, float)


def _model_profiles_on_sources(sel, model, M200, c):
    Sig, dSig = model.profiles(M200, c)
    S = np.interp(sel['R'], model.R, Sig)
    dS = np.interp(sel['R'], model.R, dSig)
    return S, dS


def _model_gt(sel, bins, rows, model, M200, c, reduced):
    """Weighted mean of the model reduced shear in each kept bin."""
    S, dS = _model_profiles_on_sources(sel, model, M200, c)
    g = reduced_shear(S, dS, sel['sigma_crit'], reduced=reduced)
    return _bin_mean(g, bins, rows, sel['w'])


def _model_kappa(sel, bins, rows, model, M200, c):
    """Weighted mean of kappa = Sigma / Sigma_crit in each kept bin."""
    S, _dS = _model_profiles_on_sources(sel, model, M200, c)
    return _bin_mean(S / sel['sigma_crit'], bins, rows, sel['w'])


def _annulus_j2(ell, theta1, theta2):
    """Area average of J_2(ell * theta) on an annulus. theta is in radians.

    The integral of theta J_2 is -2 J_0(u) - u J_1(u), with u = ell * theta.
    """
    ell = np.asarray(ell, float)
    area = float(theta2) ** 2 - float(theta1) ** 2
    if area <= 0:
        return np.zeros_like(ell)

    def _F(u):
        return -2.0 * j0(u) - u * j1(u)

    u1 = ell * float(theta1)
    u2 = ell * float(theta2)
    return 2.0 * (_F(u2) - _F(u1)) / (area * np.maximum(ell, 1e-6) ** 2)


def convergence_power(cosmo, z_source, ell, power='linear'):
    """Limber P_kappa(ell) for one source plane. Dimensionless."""
    field = cosmo.linear_field()
    chi_s = float(cosmo.comoving(z_source))
    ell = np.asarray(ell, float)
    if not np.isfinite(chi_s) or chi_s <= 0:
        return np.zeros_like(ell)
    chi = np.linspace(chi_s * 0.02, chi_s * 0.98, 32)
    z = np.interp(chi, cosmo._chi, cosmo._zs)
    a = 1.0 / (1.0 + z)
    pref = 1.5 * cosmo.Om * (cosmo.H0 / C_KMS) ** 2
    q = pref * (chi / a) * (chi_s - chi) / chi_s
    Pk = np.empty((len(chi), len(ell)))
    for i, (zz, cc) in enumerate(zip(z, chi)):
        kk = np.maximum(ell / cc, 1e-8)
        if power == 'halofit':
            k_grid, P_grid, _info = field.halofit_pk(float(zz))
            lnP = np.interp(np.log(kk), np.log(k_grid), np.log(np.maximum(P_grid, 1e-40)),
                            left=np.log(max(P_grid[0], 1e-40)), right=np.log(max(P_grid[-1], 1e-40)))
            Pk[i] = np.exp(lnP)
        else:
            Pk[i] = np.interp(kk, field.k, field.P, left=field.P[0], right=0.0) * cosmo.growth(float(zz)) ** 2
    weight = (q * q) / (chi * chi)
    return _trapz(weight[:, None] * Pk, chi, axis=0)


def lss_bin_covariance(theta1, theta2, cosmo, z_source, power='linear'):
    """Gaussian covariance of the annulus-averaged tangential shear, shape noise not included."""
    theta1 = np.asarray(theta1, float)
    theta2 = np.asarray(theta2, float)
    ell = np.geomspace(8.0, 4.0e4, 96)
    Pk = convergence_power(cosmo, z_source, ell, power=power)
    nb = len(theta1)
    J = np.vstack([_annulus_j2(ell, theta1[i], theta2[i]) for i in range(nb)])
    weight = ell * Pk / (2.0 * math.pi)
    C = np.empty((nb, nb))
    for i in range(nb):
        for j in range(i, nb):
            cij = float(_trapz(weight * J[i] * J[j], ell))
            C[i, j] = C[j, i] = cij
    return C


def _subset_selection(sel, mask):
    n = len(sel['R'])
    out = {}
    for key, value in sel.items():
        if isinstance(value, np.ndarray) and value.shape == (n,):
            out[key] = value[mask]
        else:
            out[key] = value
    return out


def jackknife_covariance(sel, edges, rows, n_jack, sigma_e, random_R):
    """Leave-one-wedge-out covariance of the binned tangential shear. None if it cannot be formed."""
    if n_jack < 2 or 'phi' not in sel:
        return None
    phi = np.asarray(sel['phi'], float)
    wedge = np.clip(np.floor((phi + math.pi) / (2.0 * math.pi) * n_jack).astype(int), 0, n_jack - 1)
    bin_ids = [r['bin'] for r in rows]
    reps = []
    for leave in range(int(n_jack)):
        keep = wedge != leave
        if int(np.sum(keep)) < 4:
            continue
        prof, _bins, _B = _profile_from_selection(_subset_selection(sel, keep), edges, sigma_e, random_R)
        by = {p['bin']: p['gt'] for p in prof}
        vec = []
        ok = True
        for b in bin_ids:
            if b not in by or not np.isfinite(by[b]):
                ok = False
                break
            vec.append(by[b])
        if ok:
            reps.append(vec)
    if len(reps) < max(3, int(n_jack) // 2):
        return None
    sample = np.asarray(reps, float)
    n = len(sample)
    delta = sample - sample.mean(axis=0)
    return (n - 1.0) / n * (delta.T @ delta)


def _cholesky_cov(C):
    C = 0.5 * (np.asarray(C, float) + np.asarray(C, float).T)
    try:
        return np.linalg.cholesky(C), 0.0
    except np.linalg.LinAlgError:
        jitter = 1e-8 * float(np.mean(np.abs(np.diag(C)))) + 1e-14
        return np.linalg.cholesky(C + np.eye(C.shape[0]) * jitter), jitter


def fit_mass(sel, edges, cosmo, z_lens, sigma_e=0.26, concentration='duffy', c_fixed=4.0,
             two_halo=True, reduced=True, miscenter=False, f_cen=0.8, sigma_off=0.1,
             random_R=None, samples=200, seed=1, min_per_bin=8, power='linear',
             lss_cov=False, n_jackknife=0, count_slope=0.0, fit_miscenter=False,
             f_cen_sigma=0.0, sigma_off_sigma=0.0):
    """Fit M_200c to the binned tangential shear. Returns a JSON-ready dict and the profile rows.

    concentration is 'duffy', 'free', or 'fixed'. samples is the Metropolis length after the
    least-squares start; 0 reports only the least-squares point. power is 'linear' or 'halofit'.
    lss_cov adds the Gaussian large-scale-structure covariance. n_jackknife > 1 replaces that
    matrix with an azimuthal jackknife when the wedges can be formed. count_slope is
    d log10 N(<m) / dm. fit_miscenter frees the centred fraction and the offset when the
    corresponding prior width is positive. The library defaults keep a shape-noise, linear fit.
    """
    rows, bins, _B = _profile_from_selection(sel, edges, sigma_e, random_R)
    rows = [r for r in rows if r['n'] >= min_per_bin and np.isfinite(r['gt_err']) and r['gt_err'] > 0]
    count_slope = float(count_slope)
    # A slope of 0 means the correction is off. 5*0-2 would otherwise be -2 and
    # would multiply the model by 1-2 kappa.
    if count_slope == 0.0:
        mag = 0.0
    else:
        mag = 5.0 * count_slope - 2.0
        if abs(mag) < 1e-12:
            mag = 0.0
    use_mis = bool(miscenter) or bool(fit_miscenter)
    free_f = bool(fit_miscenter) and float(f_cen_sigma) > 0
    free_s = bool(fit_miscenter) and float(sigma_off_sigma) > 0
    if power not in ('linear', 'halofit'):
        raise ValueError('power must be linear or halofit')
    cov_name = 'shape'
    cov_note = ''
    if two_halo and power == 'halofit':
        two_halo_limit = (
            'The two-halo term is Takahashi et al. 2012 Halofit on an Eisenstein & Hu 1998 zero-baryon spectrum, '
            'with a Tinker et al. 2010 bias. The line of sight is frozen at the lens out to 300 Mpc comoving. '
            'Full Halofit includes the one-halo power of the field, so small-scale matter can be counted twice. '
            'It is not a Boltzmann code.')
    elif two_halo:
        two_halo_limit = (
            'The two-halo term is a linear Eisenstein & Hu 1998 zero-baryon spectrum with a Tinker et al. 2010 bias. '
            'It is not a Boltzmann code and not Halofit.')
    else:
        two_halo_limit = 'The two-halo term is off.'
    if free_f or free_s:
        mis_limit = (
            'The centred fraction and the Rayleigh offset are fitted inside the prior widths supplied by the caller. '
            'Those priors are inputs, not a calibration from data.')
    elif use_mis:
        mis_limit = 'Miscentring is a centred fraction plus a Rayleigh offset set by the caller.'
    else:
        mis_limit = 'Miscentring is off.'
    if mag == 0.0 and count_slope == 0.0:
        mag_limit = 'No magnification correction (the count slope is 0).'
    elif mag == 0.0:
        mag_limit = 'The count slope is 0.4, so the magnification factor (5s-2) is zero.'
    else:
        mag_limit = (
            'Magnification uses the caller count slope s = d log10 N(<m)/dm. '
            'The model shear is multiplied by 1+(5s-2) kappa. Members dilute through the boost.')
    limits = [
        'Spherical NFW. Reduced shear is gamma/(1-kappa) at each source when reduced is on.',
        'The uncertainty is the shape noise in the bins. Large-scale structure sample variance is not included.',
        two_halo_limit,
        'Duffy et al. 2008 Table 1 (NFW, full sample, 200 rho_crit, z=0-2) was calibrated on WMAP5 simulations.',
        mis_limit,
        mag_limit,
        'No triaxial halo and no baryonic contraction of the dark-matter profile.',
    ]
    base = dict(z_lens=z_lens, H0=cosmo.H0, Om=cosmo.Om, ns=cosmo.ns, sigma8=cosmo.sigma8,
                n_src=int(len(sel['R'])), concentration=concentration, two_halo=bool(two_halo),
                reduced=bool(reduced), miscenter=use_mis, f_cen=float(f_cen), sigma_off=float(sigma_off),
                power=power, count_slope=count_slope, fit_miscenter=bool(fit_miscenter),
                f_cen_sigma=float(f_cen_sigma), sigma_off_sigma=float(sigma_off_sigma),
                covariance='shape', n_jackknife=int(n_jackknife),
                radius_unit='proper Mpc', deltasigma_unit='Msun/pc^2', limits=limits, success=False)
    if len(rows) < 2:
        base['note'] = 'fewer than two radial bins have enough sources'
        return base, rows

    Rgrid = np.geomspace(max(edges[0] * 0.5, 1e-3), edges[-1] * 1.5, 48)
    model = HaloModel(cosmo, z_lens, Rgrid, two_halo=two_halo, miscenter=use_mis,
                      f_cen=f_cen, sigma_off=sigma_off, power=power)
    obs = np.array([r['gt'] for r in rows])
    err = np.array([r['gt_err'] for r in rows])
    C = None
    if int(n_jackknife) >= 2:
        try:
            Cj = jackknife_covariance(sel, edges, rows, int(n_jackknife), sigma_e, random_R)
        except (ValueError, FloatingPointError):
            Cj = None
        if Cj is None:
            cov_note = 'jackknife was requested but not enough wedges had sources, so it was not used'
        else:
            C = Cj
            cov_name = 'jackknife'
    if C is None and lss_cov:
        try:
            Dl = float(sel['Dl'])
            th1 = np.array([edges[r['bin']] / Dl for r in rows])
            th2 = np.array([edges[r['bin'] + 1] / Dl for r in rows])
            zmed = float(np.median(sel['z']))
            Clss = lss_bin_covariance(th1, th2, cosmo, zmed, power=power)
            C = Clss + np.diag(err ** 2)
            cov_name = 'shape+lss'
        except (ValueError, FloatingPointError) as exc:
            cov_note = 'large-scale-structure covariance was not used (%s)' % exc
            C = None
            cov_name = 'shape'
    if C is not None:
        try:
            chol, _jitter = _cholesky_cov(C)
        except np.linalg.LinAlgError:
            cov_note = 'the covariance was not positive definite, so shape noise was used'
            C = None
            cov_name = 'shape'
    if C is not None:
        def whiten(v):
            return np.linalg.solve(chol, v)
        marginal = np.sqrt(np.clip(np.diag(0.5 * (C + C.T)), 0, None))
        for r, s in zip(rows, marginal):
            r['gt_err_shape'] = r['gt_err']
            r['gt_err'] = float(s)
    else:
        def whiten(v):
            return np.asarray(v, float) / err
    if cov_name == 'jackknife':
        limits[1] = 'The covariance is an azimuthal jackknife of the binned tangential shear. It replaces the shape-noise and large-scale-structure matrix.'
    elif cov_name == 'shape+lss':
        limits[1] = 'The covariance is the shape noise plus the Gaussian large-scale-structure variance of the annulus-averaged tangential shear, for one source plane at the median source redshift. It is not a survey covariance.'
    base['limits'] = limits
    base['covariance'] = cov_name
    if cov_note:
        base['covariance_note'] = cov_note
    free_c = concentration == 'free'

    def c_of(logM, logc=None):
        M = 10 ** logM
        if concentration == 'fixed':
            return float(c_fixed)
        if concentration == 'free':
            return 10 ** float(logc)
        return float(cosmo.duffy_c(M, z_lens))

    def decode(p):
        i = 0
        logM = float(p[0])
        i += 1
        logc = None
        if free_c:
            logc = float(p[i])
            i += 1
        f = float(f_cen)
        soff = float(sigma_off)
        if free_f:
            f = float(p[i])
            i += 1
        if free_s:
            soff = float(p[i])
        return logM, logc, f, soff

    def predict(p):
        logM, logc, f, soff = decode(p)
        model.f_cen = f
        model.sigma_off = soff
        M = 10 ** logM
        cc = c_of(logM, logc)
        S, dS = _model_profiles_on_sources(sel, model, M, cc)
        g = reduced_shear(S, dS, sel['sigma_crit'], reduced=reduced)
        gt = _bin_mean(g, bins, rows, sel['w'])
        if mag != 0.0:
            kap = _bin_mean(S / sel['sigma_crit'], bins, rows, sel['w'])
            gt = gt * (1.0 + mag * kap)
        return gt

    def prior_residual(p):
        _logM, _logc, f, soff = decode(p)
        extra = []
        if free_f:
            extra.append((f - float(f_cen)) / float(f_cen_sigma))
        if free_s:
            extra.append((soff - float(sigma_off)) / float(sigma_off_sigma))
        return np.asarray(extra, float)

    def residual(p):
        data = whiten(obs - predict(p))
        if free_f or free_s:
            return np.concatenate([data, prior_residual(p)])
        return data

    p0 = [14.8]
    lo = [12.5]
    hi = [16.0]
    if free_c:
        p0.append(math.log10(4.0))
        lo.append(math.log10(1.5))
        hi.append(math.log10(12.0))
    if free_f:
        p0.append(float(np.clip(f_cen, 0.05, 1.0)))
        lo.append(0.05)
        hi.append(1.0)
    if free_s:
        p0.append(float(np.clip(sigma_off, 1e-3, 2.0)))
        lo.append(1e-3)
        hi.append(2.0)
    n_extra = int(free_f) + int(free_s)
    fit = least_squares(residual, p0, bounds=(lo, hi), xtol=1e-6, ftol=1e-6, max_nfev=80 + 40 * n_extra)
    p_best = fit.x.copy()
    jac = fit.jac
    try:
        jtj = jac.T @ jac
        cov = np.linalg.inv(jtj)
        curv = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        cov = None
        curv = np.full(len(p_best), 0.05)

    def loglike(p):
        if np.any(p < np.array(lo)) or np.any(p > np.array(hi)):
            return -1e300
        r = residual(p)
        if not np.all(np.isfinite(r)):
            return -1e300
        return -0.5 * float(np.sum(r * r))

    chain = None
    if samples and samples > 0:
        rng = np.random.default_rng(seed)
        step = np.maximum(curv, 0.01)
        x = p_best.copy()
        lp = loglike(x)
        burn = max(40, samples // 3)
        acc = 0
        kept = []
        total = burn + int(samples)
        for i in range(total):
            prop = x + rng.normal(0.0, step)
            lp2 = loglike(prop)
            if math.log(rng.random()) < lp2 - lp:
                x, lp = prop, lp2
                acc += 1
            if i == burn // 2:
                rate = acc / max(i, 1)
                if rate < 0.15:
                    step *= 0.5
                elif rate > 0.5:
                    step *= 1.5
                acc = 0
            if i >= burn:
                kept.append(x.copy())
        chain = np.array(kept)

    def pack(p):
        logM, logc, f, soff = decode(p)
        return 10 ** logM, c_of(logM, logc), f, soff

    if chain is not None and len(chain):
        packed = [pack(rowp) for rowp in chain]
        Ms = np.array([t[0] for t in packed])
        cs = np.array([t[1] for t in packed])
        fs = np.array([t[2] for t in packed])
        ss = np.array([t[3] for t in packed])
        M_med = float(np.median(Ms))
        c_med = float(np.median(cs))
        f_med = float(np.median(fs))
        s_med = float(np.median(ss))
        M_p16, M_p84 = [float(v) for v in np.percentile(Ms, [16, 84])]
        c_p16, c_p84 = [float(v) for v in np.percentile(cs, [16, 84])]
        method = 'metropolis'
        n_keep = int(len(chain))
    else:
        M_med, c_med, f_med, s_med = pack(p_best)
        dlog = float(curv[0]) if len(curv) else 0.0
        M_p16 = 10 ** (math.log10(M_med) - dlog)
        M_p84 = 10 ** (math.log10(M_med) + dlog)
        if free_c and len(curv) > 1:
            c_p16 = 10 ** (math.log10(c_med) - float(curv[1]))
            c_p84 = 10 ** (math.log10(c_med) + float(curv[1]))
        else:
            c_p16 = c_p84 = c_med
        method = 'least_squares'
        n_keep = 0

    model.f_cen = f_med
    model.sigma_off = s_med
    S_fin, dS_fin = _model_profiles_on_sources(sel, model, M_med, c_med)
    g_fin = reduced_shear(S_fin, dS_fin, sel['sigma_crit'], reduced=reduced)
    gmodel = _bin_mean(g_fin, bins, rows, sel['w'])
    if mag != 0.0:
        kap_fin = _bin_mean(S_fin / sel['sigma_crit'], bins, rows, sel['w'])
        gmodel = gmodel * (1.0 + mag * kap_fin)
    chi_t = float(np.sum(whiten(obs - gmodel) ** 2))
    chi_x = float(np.sum((np.array([r['gx'] for r in rows]) / err) ** 2))
    npar = 1 + int(free_c) + int(free_f) + int(free_s)
    _Sig, dSig = model.profiles(M_med, c_med)
    for r, gm in zip(rows, gmodel):
        r['gt_model'] = float(gm) if np.isfinite(gm) else None
        r['deltasigma_model'] = float(np.interp(r['R'], model.R, dSig) / 1e12)
    span = np.array(hi) - np.array(lo)
    at_bound = bool(np.any(np.abs(p_best - np.array(lo)) <= 1e-4 * np.maximum(span, 1.0))
                    or np.any(np.abs(np.array(hi) - p_best) <= 1e-4 * np.maximum(span, 1.0)))
    base.update(success=bool(fit.success) or method == 'metropolis', at_bound=at_bound, method=method,
                M200c_Msun=M_med, M200c_p16=M_p16, M200c_p84=M_p84,
                c200=c_med, c200_p16=c_p16, c200_p84=c_p84,
                f_cen=float(f_med), sigma_off=float(s_med),
                M200c_least_squares=float(10 ** p_best[0]),
                chi2=chi_t, chi2_cross=chi_x, dof=int(len(rows) - npar),
                n_bins=len(rows), n_samples=n_keep, nfev=int(fit.nfev))
    if at_bound:
        base['note'] = 'a fitted parameter is at a bound'
    return base, rows


def _plain(value):
    """JSON-ready copy. Non-finite floats become None."""
    if isinstance(value, dict):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(item) for item in value]
    if isinstance(value, np.ndarray):
        return _plain(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        number = float(value)
        return number if np.isfinite(number) else None
    return value


def r200c_mpc(mass, z, cosmo):
    """Proper R_200c in Mpc. M_200c is a solar mass."""
    mass = float(mass)
    z = float(z)
    if not np.isfinite(mass) or mass <= 0:
        raise ValueError('mass must be positive')
    if not np.isfinite(z) or z < 0:
        raise ValueError('redshift must be non-negative')
    rho = rho_crit(z, cosmo.H0, cosmo.Om)
    return (3.0 * mass / (4.0 * math.pi * 200.0 * rho)) ** (1.0 / 3.0)


def photometric_geometry(cosmo, z_lens, z_grid, probability, z_gap=0.1,
                         cluster_window=None, background_window=None):
    """Label window and lensing efficiency from a sample n(z) or per-source P(z).

    `probability` matches `z_grid`, or has shape (n_source, n_z). The mean of
    Sigma_crit^{-1} is the probability-weighted mean behind the lens.
    sigma_crit_inv_of_mean_z evaluates Sigma_crit at the probability-weighted mean
    redshift. The two numbers stay separate. No spectroscopic redshift column is filled.
    """
    z_grid = np.asarray(z_grid, float).reshape(-1)
    prob = np.asarray(probability, float)
    if prob.ndim == 1:
        if prob.shape != z_grid.shape:
            raise ValueError('n(z) length does not match z_grid')
        prob = prob.reshape(1, -1)
    elif prob.ndim != 2 or prob.shape[1] != z_grid.size:
        raise ValueError('P(z) must have one column per redshift')
    behind = (z_grid > float(z_lens) + float(z_gap)) & np.isfinite(z_grid)
    if not np.any(behind):
        raise ValueError('n(z) has no source behind the lens')
    zz = z_grid[behind]
    pp = np.where(np.isfinite(prob[:, behind]), np.clip(prob[:, behind], 0.0, None), 0.0)
    weight = pp.sum(axis=1)
    usable = weight > 0
    if not np.any(usable):
        raise ValueError('n(z) has no source behind the lens')
    inv = 1.0 / np.asarray(cosmo.sigma_crit(z_lens, zz), float)
    finite = np.isfinite(inv) & (inv > 0)
    if not np.any(finite):
        raise ValueError('Sigma_crit is not finite behind the lens')
    pp = pp[:, finite]
    zz = zz[finite]
    inv = inv[finite]
    weight = pp.sum(axis=1)
    usable = weight > 0
    per_inv = (pp[usable] * inv).sum(axis=1) / weight[usable]
    per_z = (pp[usable] * zz).sum(axis=1) / weight[usable]
    src_w = weight[usable]
    mean_inv = float(np.sum(src_w * per_inv) / np.sum(src_w))
    mean_z = float(np.sum(src_w * per_z) / np.sum(src_w))
    inv_at_mean = float(1.0 / np.asarray(cosmo.sigma_crit(z_lens, mean_z), float))
    return dict(
        label_window=dict(
            cluster=None if cluster_window is None else [float(cluster_window[0]), float(cluster_window[1])],
            background=None if background_window is None else [float(background_window[0]), float(background_window[1])]),
        mean_sigma_crit_inv=mean_inv,
        sigma_crit_inv_of_mean_z=inv_at_mean,
        mean_z=mean_z,
        n_sources=int(np.count_nonzero(usable)),
        note=('The mean of Sigma_crit^{-1} and Sigma_crit^{-1} at the mean redshift are different quantities. '
              'This result does not fill a spectroscopic redshift column.'))


def _column(record, name, n=None, dtype=float):
    if name not in record:
        raise ValueError('cluster record needs %s' % name)
    values = np.asarray(record[name], dtype).reshape(-1)
    if n is not None and values.size != n:
        raise ValueError('%s length does not match the catalogue' % name)
    return values


def _reject_image_input(clusters):
    if isinstance(clusters, np.ndarray):
        raise ValueError('a stacked image is not a shear catalogue')
    if isinstance(clusters, (list, tuple)) and clusters and isinstance(clusters[0], np.ndarray):
        raise ValueError('a stacked image is not a shear catalogue')
    if isinstance(clusters, dict) and 'image' in clusters and 'g1' not in clusters:
        raise ValueError('a stacked image is not a shear catalogue')


def _filter_cluster_rows(record):
    """Drop quality failures and bad cross-matches before select_sources.

    A finite spectroscopic redshift is not a quality flag. A photometric point
    estimate is not accepted in place of the spectroscopic column.
    """
    if not isinstance(record, dict):
        raise ValueError('a stacked image is not a shear catalogue')
    if 'image' in record and 'g1' not in record:
        raise ValueError('a stacked image is not a shear catalogue')
    if 'z_phot' in record and 'z' not in record:
        raise ValueError('a photometric point redshift is not a spectroscopic column')
    x = _column(record, 'x_arcsec')
    n = int(x.size)
    y = _column(record, 'y_arcsec', n)
    g1 = _column(record, 'g1', n)
    g2 = _column(record, 'g2', n)
    z = _column(record, 'z', n) if 'z' in record else np.full(n, np.nan)
    w = _column(record, 'weight', n)
    keep = np.ones(n, bool)
    quality_supplied = 'spec_quality' in record
    if quality_supplied:
        quality = np.asarray(record['spec_quality']).reshape(-1)
        if quality.size != n:
            raise ValueError('spec_quality length does not match the catalogue')
        if quality.dtype == bool or quality.dtype == np.bool_:
            keep &= quality.astype(bool)
        else:
            quality = np.asarray(quality, float)
            minimum = float(record.get('min_spec_quality', 1))
            keep &= np.isfinite(quality) & (quality >= minimum)
    if 'crossmatch_ok' in record:
        matched = np.asarray(record['crossmatch_ok'], bool).reshape(-1)
        if matched.size != n:
            raise ValueError('crossmatch_ok length does not match the catalogue')
        keep &= matched
    shape = np.isfinite(g1) & np.isfinite(g2) & np.isfinite(x) & np.isfinite(y)
    if quality_supplied:
        intersection = int(np.count_nonzero(keep & shape & np.isfinite(z)))
    else:
        intersection = None
    member = None
    if record.get('member') is not None:
        member = np.asarray(record['member']).reshape(-1)
        if member.size != n:
            raise ValueError('member labels do not match the catalogue')
        member = member[keep]
    return dict(
        x=x[keep], y=y[keep], g1=g1[keep], g2=g2[keep], z=z[keep], w=w[keep], member=member,
        n_in=n, n_kept=int(np.count_nonzero(keep)), n_spec_intersection=intersection,
        spec_quality_supplied=quality_supplied)


def _weights_certified(selections, inverse_variance_weights):
    if inverse_variance_weights is True:
        return True
    if inverse_variance_weights is False:
        return False
    for sel in selections:
        weight = np.asarray(sel['w'], float)
        if weight.size and not np.allclose(weight, 1.0):
            return False
    return True


def _annotate_shape_error(rows, sigma_e, certified, boost_uncertainty):
    for row in rows:
        if not certified or not (row['sw'] > 0):
            row['gt_err'] = None
            row['gx_err'] = None
            row['deltasigma_err'] = None
            row['error_formula'] = 'omitted'
            continue
        boost = float(row['boost'])
        err = float(sigma_e) / math.sqrt(float(row['sw'])) * boost
        if boost_uncertainty is not None and np.isfinite(boost_uncertainty) and boost != 0.0:
            raw = float(row['gt']) / boost
            err = math.sqrt(err * err + (raw * float(boost_uncertainty)) ** 2)
        row['gt_err'] = float(err)
        row['gx_err'] = float(err)
        inv = row['sigma_crit_inv']
        row['deltasigma_err'] = float(err / inv / 1e12) if inv and np.isfinite(inv) and inv > 0 else None
        row['error_formula'] = 'sigma_e/sqrt(sum w), with the boost factor'
    return rows


def _mag_factor(count_slope):
    count_slope = float(count_slope)
    if count_slope == 0.0:
        return 0.0
    mag = 5.0 * count_slope - 2.0
    return 0.0 if abs(mag) < 1e-12 else mag


def _concentration_of(mode, mass, z, cosmo, c_fixed, logc):
    if mode == 'fixed':
        return float(c_fixed)
    if mode == 'free':
        return 10.0 ** float(logc)
    if mode == 'duffy':
        return float(cosmo.duffy_c(mass, z))
    raise ValueError('concentration must be duffy, free, or fixed')


def _model_on_bins(state, mass, concentration, reduced, mag):
    shear = _model_gt(state['sel'], state['bins'], state['rows'], state['model'], mass, concentration, reduced)
    if mag == 0.0:
        return {row['bin']: float(value) for row, value in zip(state['rows'], shear)}
    kappa = _model_kappa(state['sel'], state['bins'], state['rows'], state['model'], mass, concentration)
    out = {}
    for row, value, kap in zip(state['rows'], shear, kappa):
        out[row['bin']] = float(value * (1.0 + mag * kap)) if np.isfinite(value) and np.isfinite(kap) else np.nan
    return out


def _average_model(states, obs_rows, mass, mode, c_fixed, logc, cosmo, reduced, mag):
    predicted = []
    for state in states:
        conc = _concentration_of(mode, mass, state['z_lens'], cosmo, c_fixed, logc)
        state['_cache'] = _model_on_bins(state, mass, conc, reduced, mag)
    for row in obs_rows:
        values = []
        weights = []
        for index, weight in row['contributors']:
            value = states[index]['_cache'].get(row['bin'], np.nan)
            if np.isfinite(value) and weight > 0:
                values.append(value)
                weights.append(weight)
        if not values:
            predicted.append(np.nan)
        else:
            weights = np.asarray(weights, float)
            predicted.append(float(np.sum(weights * np.asarray(values, float)) / np.sum(weights)))
    return np.asarray(predicted, float)


def _fit_average(states, obs_rows, cosmo, concentration, c_fixed, reduced, count_slope,
                 samples, seed, covariance=None):
    """Least-squares NFW on the stacked tangential shear. samples>0 adds a short Metropolis."""
    if concentration not in ('duffy', 'free', 'fixed'):
        raise ValueError('concentration must be duffy, free, or fixed')
    mag = _mag_factor(count_slope)
    free_c = concentration == 'free'
    obs = np.array([row['gt'] for row in obs_rows], float)
    if len(obs_rows) < 2 or not np.isfinite(obs).all():
        return dict(success=False, note='fewer than two radial bins have a tangential shear')
    certified = all(row.get('gt_err') not in (None, 0) and np.isfinite(row['gt_err']) and row['gt_err'] > 0
                    for row in obs_rows)
    if certified:
        err = np.array([row['gt_err'] for row in obs_rows], float)
    else:
        err = np.ones(len(obs_rows), float)
    chol = None
    if covariance is not None:
        cov = np.asarray(covariance, float)
        if cov.shape != (len(obs_rows), len(obs_rows)):
            raise ValueError('covariance shape does not match the stacked bins')
        chol, _jitter = _cholesky_cov(cov)

    def residual(p):
        logc = float(p[1]) if free_c else None
        model = _average_model(states, obs_rows, 10.0 ** float(p[0]), concentration, c_fixed, logc, cosmo, reduced, mag)
        if not np.isfinite(model).all():
            return np.full(obs.shape, 1e3)
        delta = obs - model
        if chol is not None:
            return np.linalg.solve(chol, delta)
        return delta / err

    p0 = [14.8]
    lo = [12.5]
    hi = [16.0]
    if free_c:
        p0.append(math.log10(4.0))
        lo.append(math.log10(1.5))
        hi.append(math.log10(12.0))
    fit = least_squares(residual, p0, bounds=(lo, hi), xtol=1e-6, ftol=1e-6, max_nfev=80)
    try:
        cov = np.linalg.inv(fit.jac.T @ fit.jac)
        curv = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        curv = np.full(len(fit.x), 0.05)

    def pack(p):
        logc = float(p[1]) if free_c else None
        mass = 10.0 ** float(p[0])
        # One concentration number for the record: the first lens, or the fixed value.
        conc = _concentration_of(concentration, mass, states[0]['z_lens'], cosmo, c_fixed, logc)
        return mass, conc

    chain_mass = None
    if samples and int(samples) > 0:
        rng = np.random.default_rng(seed)
        step = np.maximum(curv, 0.01)
        x = fit.x.copy()
        lp = -0.5 * float(np.sum(residual(x) ** 2))
        burn = max(20, int(samples) // 3)
        kept = []
        for i in range(burn + int(samples)):
            prop = x + rng.normal(0.0, step)
            if np.any(prop < np.asarray(lo)) or np.any(prop > np.asarray(hi)):
                continue
            trial = residual(prop)
            if not np.isfinite(trial).all():
                continue
            lp2 = -0.5 * float(np.sum(trial ** 2))
            if math.log(rng.random()) < lp2 - lp:
                x, lp = prop, lp2
            if i >= burn:
                kept.append(x.copy())
        if kept:
            packed = [pack(row) for row in kept]
            chain_mass = np.array([item[0] for item in packed])
            chain_conc = np.array([item[1] for item in packed])
    mass, conc = pack(fit.x)
    if chain_mass is not None:
        mass = float(np.median(chain_mass))
        conc = float(np.median(chain_conc))
        p16, p84 = [float(v) for v in np.percentile(chain_mass, [16, 84])]
        method = 'metropolis'
        n_keep = int(chain_mass.size)
    else:
        dlog = float(curv[0]) if len(curv) else 0.0
        p16 = 10 ** (math.log10(mass) - dlog)
        p84 = 10 ** (math.log10(mass) + dlog)
        method = 'least_squares'
        n_keep = 0
    model = _average_model(states, obs_rows, mass, concentration, c_fixed,
                           math.log10(conc) if free_c else None, cosmo, reduced, mag)
    final = [math.log10(mass)]
    if free_c:
        final.append(math.log10(conc))
    chi2 = float(np.sum(residual(np.asarray(final, float)) ** 2))
    npar = 1 + int(free_c)
    span = np.asarray(hi) - np.asarray(lo)
    at_bound = bool(np.any(np.abs(fit.x - np.asarray(lo)) <= 1e-4 * np.maximum(span, 1.0))
                    or np.any(np.abs(np.asarray(hi) - fit.x) <= 1e-4 * np.maximum(span, 1.0)))
    for row, pred in zip(obs_rows, model):
        row['gt_model'] = float(pred) if np.isfinite(pred) else None
    return dict(success=bool(fit.success) or method == 'metropolis', at_bound=at_bound, method=method,
                M200c_Msun=float(mass), M200c_p16=float(p16), M200c_p84=float(p84),
                c200=float(conc), chi2=chi2, dof=int(len(obs_rows) - npar),
                n_bins=len(obs_rows), n_samples=n_keep, concentration=concentration,
                error_formula='sigma_e/sqrt(sum w)' if certified else 'omitted')


def _halo(state, cosmo, two_halo, power, miscenter, f_cen, sigma_off):
    radii = np.asarray(state['sel']['R'], float)
    grid = np.geomspace(max(float(np.min(radii)) * 0.5, 1e-3), float(np.max(radii)) * 1.5, 48)
    return HaloModel(cosmo, state['z_lens'], grid, two_halo=two_halo, miscenter=miscenter,
                     f_cen=f_cen, sigma_off=sigma_off, power=power)


def _refit_variant(states, obs_rows, cosmo, settings, miscenter=False, f_cen=0.8, sigma_off=0.1, rows=None):
    """One extra least-squares fit. The caller passes the rows that define the radial window."""
    saved = [(state['model'], state['rows']) for state in states]
    used = obs_rows if rows is None else rows
    if len(used) < 2:
        for state, (model, old_rows) in zip(states, saved):
            state['model'] = model
            state['rows'] = old_rows
        return None
    for state in states:
        state['model'] = _halo(state, cosmo, settings['two_halo'], settings['power'], miscenter, f_cen, sigma_off)
        state['rows'] = [row for row in state['rows'] if any(row['bin'] == item['bin'] for item in used)]
    # A cluster row list that no longer matches `bins` still works: _model_gt indexes by bin.
    try:
        fitted = _fit_average(states, used, cosmo, settings['concentration'], settings['c_fixed'],
                              settings['reduced'], settings['count_slope'], 0, 1)
    finally:
        for state, (model, old_rows) in zip(states, saved):
            state['model'] = model
            state['rows'] = old_rows
    if not fitted.get('success'):
        return None
    return fitted.get('M200c_Msun')


def _sensitivity_table(states, obs_rows, cosmo, settings, reference, supplied, bcg_offset):
    supplied = {} if supplied is None else dict(supplied)
    concentration_rows = [dict(name='reference', concentration=settings['concentration'],
                               c_fixed=settings['c_fixed'], M200c_Msun=reference.get('M200c_Msun'))]
    other_c = float(settings['c_fixed']) * 1.5 if settings['concentration'] == 'fixed' else 4.0
    other_mode = 'fixed'
    alt_settings = dict(settings, concentration=other_mode, c_fixed=other_c)
    # The variant fitter reads settings at call time; swap the closure values via a copied dict
    # by temporarily calling _fit_average through _refit_variant, which uses `settings`.
    alt_mass = None
    saved_settings = dict(settings)
    settings.update(concentration=other_mode, c_fixed=other_c)
    try:
        alt_mass = _refit_variant(states, obs_rows, cosmo, settings)
    finally:
        settings.clear()
        settings.update(saved_settings)
    concentration_rows.append(dict(name='varied', concentration=other_mode, c_fixed=other_c, M200c_Msun=alt_mass))
    centred = _refit_variant(states, obs_rows, cosmo, settings, miscenter=True, f_cen=0.8, sigma_off=0.1)
    inner = [row for row in obs_rows if row['bin'] != obs_rows[-1]['bin']]
    outer = [row for row in obs_rows if row['bin'] != obs_rows[0]['bin']]
    radius_rows = [
        dict(name='drop_outer_bin', M200c_Msun=_refit_variant(states, obs_rows, cosmo, settings, rows=inner)),
        dict(name='drop_inner_bin', M200c_Msun=_refit_variant(states, obs_rows, cosmo, settings, rows=outer))]

    def caller_row(key, fallback):
        block = supplied.get(key)
        if isinstance(block, dict) and str(block.get('note', '')).strip():
            out = dict(block)
            out['supplied'] = True
            out['source'] = 'caller'
            return out
        return dict(supplied=False, source='model', note=fallback)

    return dict(
        concentration=concentration_rows,
        centering=dict(f_cen=0.8, sigma_off_mpc=0.1, M200c_Msun=centred,
                       note='Centring variation is a centred fraction 0.8 and a Rayleigh offset of 0.1 Mpc.'),
        fit_radius=radius_rows,
        non_sphericity=caller_row(
            'non_sphericity', 'The reference model is spherical. No triaxial mass is computed here.'),
        projection=caller_row(
            'projection', 'Line-of-sight structure is not varied here.'),
        bcg_mass_offset_mpc=None if bcg_offset is None else float(bcg_offset))


def _science_gaps(intersection, quality_supplied, fit, sensitivity):
    gaps = []
    if not quality_supplied or intersection is None or intersection <= 0:
        gaps.append('spec_intersection')
    if fit is None or not fit.get('success') or int(fit.get('dof', 0)) <= 0 or fit.get('at_bound'):
        gaps.append('fit_constraint')
    if sensitivity is None:
        gaps.append('sensitivity')
    else:
        if not sensitivity['non_sphericity'].get('supplied'):
            gaps.append('non_sphericity')
        if not sensitivity['projection'].get('supplied'):
            gaps.append('projection')
        if sensitivity.get('bcg_mass_offset_mpc') is None:
            gaps.append('bcg_mass_offset')
    return gaps


def _one_source_plane(states):
    redshifts = np.concatenate([np.asarray(state['sel']['z'], float) for state in states]) if states else np.zeros(0)
    lenses = np.array([state['z_lens'] for state in states], float)
    if redshifts.size == 0 or lenses.size == 0:
        return False
    return float(np.nanmax(redshifts) - np.nanmin(redshifts)) < 1e-4 and float(np.nanmax(lenses) - np.nanmin(lenses)) < 1e-4


def stack_shear_profiles(clusters, cosmo=None, radius='mpc', rmin=0.3, rmax=3.0, nbins=8,
                         z_gap=0.1, sigma_e=0.26, min_per_bin=8, gmax=1.5,
                         inverse_variance_weights=None, boost_uncertainty=None,
                         icl_weights_differ=False, concentration='fixed', c_fixed=4.0,
                         two_halo=False, reduced=True, power='linear', lss_cov=False,
                         count_slope=0.0, samples=0, seed=1, per_cluster_mass=True, fit=True,
                         sensitivity=False, sensitivity_supplied=None, bcg_mass_offset_mpc=None,
                         miscenter=False, f_cen=0.8, sigma_off=0.1):
    """Average tangential-shear profiles, then fit that average with a spherical NFW.

    Each cluster record holds x_arcsec, y_arcsec, g1, g2, z, weight, and z_lens.
    spec_quality and crossmatch_ok, when present, are applied before select_sources.
    radius is 'mpc' or 'r200c'. An r200c stack reads record['mass'] as M_200c.
    The default samples value is 0, so the fit reports the least-squares point.
    Pixel images are rejected.
    """
    _reject_image_input(clusters)
    if not isinstance(clusters, (list, tuple)) or len(clusters) == 0:
        raise ValueError('no clusters to stack')
    if radius not in ('mpc', 'r200c'):
        raise ValueError('radius must be mpc or r200c')
    if power not in ('linear', 'halofit'):
        raise ValueError('power must be linear or halofit')
    cosmo = Cosmology() if cosmo is None else cosmo
    edges = bin_edges(rmin, rmax, nbins)
    states = []
    photometric = []
    intersections = []
    quality_flags = []
    for i, record in enumerate(clusters):
        if isinstance(record, np.ndarray):
            raise ValueError('a stacked image is not a shear catalogue')
        filt = _filter_cluster_rows(record)
        quality_flags.append(filt['spec_quality_supplied'])
        intersections.append(filt['n_spec_intersection'])
        if 'z_grid' in record and 'probability' in record:
            photometric.append(photometric_geometry(
                cosmo, float(record['z_lens']), record['z_grid'], record['probability'], z_gap=z_gap,
                cluster_window=record.get('cluster_window'), background_window=record.get('background_window')))
        if filt['n_kept'] == 0:
            continue
        z_lens = float(record['z_lens'])
        scale = 1.0
        if radius == 'r200c':
            if 'mass' not in record:
                raise ValueError('an r200c stack needs a mass on each cluster')
            scale = r200c_mpc(record['mass'], z_lens, cosmo)
        sel = select_sources(
            filt['x'], filt['y'], filt['g1'], filt['g2'], filt['z'], filt['w'], z_lens, cosmo,
            z_gap=z_gap, rmin=rmin * scale, rmax=rmax * scale, gmax=gmax, member=filt['member'])
        if len(sel['R']) == 0:
            continue
        scaled = dict(sel)
        scaled['R'] = sel['R'] / scale
        random_R = record.get('random_R')
        if random_R is not None:
            random_R = np.asarray(random_R, float) / scale
        rows, bins, _boost = _profile_from_selection(scaled, edges, sigma_e, random_R)
        rows = [row for row in rows if row['n'] >= int(min_per_bin) and np.isfinite(row['gt'])]
        if not rows:
            continue
        state = dict(id=str(record['id']) if record.get('id') is not None else str(i),
                     z_lens=z_lens, sel=sel, rows=rows, bins=bins, scale=float(scale),
                     n_spec_intersection=filt['n_spec_intersection'])
        state['model'] = _halo(state, cosmo, two_halo, power, miscenter, f_cen, sigma_off)
        states.append(state)
    quality_supplied = any(quality_flags)
    if quality_supplied:
        intersection_total = int(sum(value or 0 for value in intersections))
    else:
        intersection_total = None
    population = bool(icl_weights_differ) or any(isinstance(record, dict) and 'icl_weight' in record for record in clusters)
    base = dict(
        product='tangential shear profile stack',
        radius=radius,
        radius_unit='proper Mpc' if radius == 'mpc' else 'R/R_200c',
        deltasigma_unit='Msun/pc^2',
        deltasigma_note=DELTASIGMA_NOTE,
        density_note=DENSITY_NOTE,
        population_note=POPULATION_NOTE if population else None,
        n_clusters=len(clusters),
        n_clusters_used=len(states),
        n_spec_intersection=intersection_total,
        spec_quality_note=None if quality_supplied else 'A finite redshift was not treated as a spectroscopic quality flag.',
        photometric=photometric,
        metacal_scope=METACAL_SCOPE,
        reference='Shirasaki & Takada 2018 MNRAS 478, 4277',
        science_measurement=False)
    if radius == 'r200c':
        base['radius_note'] = (
            'The scaled radius uses the caller mass as R_200c. '
            'That mass is not an error-free alignment coordinate.')
    if quality_supplied and intersection_total == 0:
        base.update(success=False, status='insufficient_sample', profile=[],
                    note='The spectroscopic intersection has no sources.',
                    science_gaps=['spec_intersection'])
        return _plain(base)
    if not states:
        base.update(success=False, status='insufficient_sample', profile=[],
                    note='No cluster retained a source behind the lens.',
                    science_gaps=['spec_intersection'])
        return _plain(base)

    certified = _weights_certified([state['sel'] for state in states], inverse_variance_weights)
    for state in states:
        _annotate_shape_error(state['rows'], sigma_e, certified, boost_uncertainty)
    by_bin = {}
    for index, state in enumerate(states):
        for row in state['rows']:
            err = row.get('gt_err')
            weight = (1.0 / err ** 2) if certified and err not in (None,) and err > 0 else 1.0
            by_bin.setdefault(row['bin'], []).append((index, weight, row))
    profile = []
    for bin_id in sorted(by_bin):
        parts = by_bin[bin_id]
        weights = np.array([part[1] for part in parts], float)
        gt = np.array([part[2]['gt'] for part in parts], float)
        gx = np.array([part[2]['gx'] for part in parts], float)
        inv = np.array([part[2]['sigma_crit_inv'] for part in parts], float)
        mean_gt = float(np.sum(weights * gt) / np.sum(weights))
        mean_gx = float(np.sum(weights * gx) / np.sum(weights))
        mean_inv = float(np.sum(weights * inv) / np.sum(weights))
        if certified and np.all(np.array([part[2]['gt_err'] or 0 for part in parts]) > 0):
            err = float(1.0 / math.sqrt(np.sum(weights)))
            formula = 'sigma_e/sqrt(sum w), clusters combined as independent shape noise'
        else:
            err = None
            formula = 'omitted'
        dS = mean_gt / mean_inv if mean_inv and np.isfinite(mean_inv) and mean_inv > 0 else None
        profile.append(dict(
            bin=int(bin_id),
            R=float((edges[bin_id] * edges[bin_id + 1]) ** 0.5),
            n=int(sum(part[2]['n'] for part in parts)),
            n_clusters=len(parts),
            ids=[states[part[0]]['id'] for part in parts],
            gt=mean_gt, gx=mean_gx, gt_err=err,
            deltasigma=(dS / 1e12) if dS is not None else None,
            deltasigma_err=(err / mean_inv / 1e12) if err is not None and dS is not None else None,
            sigma_crit_inv=mean_inv,
            boost=float(np.mean([part[2]['boost'] for part in parts])),
            contributors=[(part[0], float(part[1])) for part in parts],
            error_formula=formula))
    covariance = 'shape' if certified else 'omitted'
    covariance_note = '' if certified else WEIGHT_NOTE
    fit_covariance = None
    lss_marginal = None
    if lss_cov:
        if _one_source_plane(states) and radius == 'mpc' and certified and len(profile) >= 2:
            distance = float(states[0]['sel']['Dl'])
            theta1 = np.array([edges[row['bin']] / distance for row in profile])
            theta2 = np.array([edges[row['bin'] + 1] / distance for row in profile])
            z_source = float(np.median(states[0]['sel']['z']))
            try:
                matter = lss_bin_covariance(theta1, theta2, cosmo, z_source, power=power)
                shape = np.diag([row['gt_err'] ** 2 for row in profile])
                fit_covariance = matter + shape
                lss_marginal = np.sqrt(np.clip(np.diag(0.5 * (fit_covariance + fit_covariance.T)), 0, None))
                covariance = 'shape+lss'
                covariance_note = LSS_USED_NOTE
            except (ValueError, FloatingPointError) as exc:
                covariance_note = '%s (%s)' % (LSS_SKIP_NOTE, exc)
        else:
            covariance_note = LSS_SKIP_NOTE
    for row in profile:
        row.pop('contributors', None)
    # The fit closure needs the contributors that pop() removed. Rebuild from by_bin.
    for row in profile:
        row['contributors'] = [(part[0], float(part[1])) for part in by_bin[row['bin']]]
    base['profile'] = profile
    base['covariance'] = covariance
    if covariance_note:
        base['covariance_note'] = covariance_note
    settings = dict(concentration=concentration, c_fixed=float(c_fixed), two_halo=bool(two_halo),
                    power=power, reduced=bool(reduced), count_slope=float(count_slope))
    per_cluster = []
    if per_cluster_mass:
        for state in states:
            mpc_edges = edges * state['scale']
            summary, _rows = fit_mass(
                state['sel'], mpc_edges, cosmo, state['z_lens'], sigma_e=sigma_e,
                concentration=concentration, c_fixed=c_fixed, two_halo=two_halo, reduced=reduced,
                power=power, lss_cov=False, samples=int(samples), seed=seed,
                min_per_bin=min_per_bin, count_slope=count_slope, miscenter=miscenter,
                f_cen=f_cen, sigma_off=sigma_off)
            summary['id'] = state['id']
            per_cluster.append(summary)
    base['per_cluster'] = per_cluster
    fitted = None
    saved_model = [row.get('gt_model') for row in profile]
    if fit:
        fitted = _fit_average(
            states, profile, cosmo, concentration, c_fixed, reduced, count_slope, samples, seed,
            covariance=fit_covariance)
        saved_model = [row.get('gt_model') for row in profile]
        fitted['covariance'] = covariance
        if covariance_note:
            fitted['covariance_note'] = covariance_note
        base['fit'] = fitted
        base['success'] = bool(fitted.get('success'))
        base['M200c_Msun'] = fitted.get('M200c_Msun')
        base['c200'] = fitted.get('c200')
        base['at_bound'] = fitted.get('at_bound')
        base['status'] = 'stacked'
    else:
        base['success'] = False
        base['status'] = 'profile_only'
        base['fit'] = None
    table = None
    if sensitivity and fitted is not None and fitted.get('success'):
        table = _sensitivity_table(states, profile, cosmo, settings, fitted, sensitivity_supplied, bcg_mass_offset_mpc)
    elif sensitivity_supplied or bcg_mass_offset_mpc is not None:
        table = _sensitivity_table(states, profile, cosmo, settings, fitted or {}, sensitivity_supplied, bcg_mass_offset_mpc)
    base['sensitivity'] = table
    for row, value in zip(profile, saved_model):
        row['gt_model'] = value
    if lss_marginal is not None:
        for row, value in zip(profile, lss_marginal):
            row['gt_err_shape'] = row['gt_err']
            row['gt_err'] = float(value)
            if row.get('sigma_crit_inv'):
                row['deltasigma_err'] = float(value / row['sigma_crit_inv'] / 1e12)
    gaps = _science_gaps(intersection_total, quality_supplied, fitted, table)
    base['science_gaps'] = gaps
    base['science_measurement'] = len(gaps) == 0
    base['science_note'] = (
        'science_measurement records whether the caller filled the sensitivity and spectroscopic-intersection gates. '
        'A filled gate is not a named-data mass and not a survey calibration.')
    for row in profile:
        row.pop('contributors', None)
    return _plain(base)


# ------------------------------------------------------------------------------- shapes: metacalibration of one stamp
def adaptive_moments(image, niter=30):
    """Gaussian-weighted second moments. Returns e1, e2, T=mxx+myy, flux, or None."""
    img = np.asarray(image, float)
    if img.ndim != 2 or min(img.shape) < 5:
        return None
    n0, n1 = img.shape
    yy, xx = np.mgrid[0:n0, 0:n1]
    pos = np.clip(img, 0, None)
    flux0 = pos.sum()
    if flux0 <= 0 or not np.isfinite(flux0):
        return None
    cx = float((pos * xx).sum() / flux0)
    cy = float((pos * yy).sum() / flux0)
    dx, dy = xx - cx, yy - cy
    sigma2 = float((pos * (dx * dx + dy * dy)).sum() / flux0 / 2.0)
    mxx = myy = max(sigma2, 0.4)
    mxy = 0.0
    for _ in range(niter):
        det = mxx * myy - mxy * mxy
        if det <= 1e-8:
            return None
        ixx, iyy, ixy = myy / det, mxx / det, -mxy / det
        chi2 = ixx * dx * dx + iyy * dy * dy + 2 * ixy * dx * dy
        ww = np.exp(-0.5 * np.clip(chi2, 0, 40)) * pos
        s = float(ww.sum())
        if s <= 0:
            return None
        cx = float((ww * xx).sum() / s)
        cy = float((ww * yy).sum() / s)
        dx, dy = xx - cx, yy - cy
        # For a Gaussian galaxy and a matched Gaussian weight, the weighted
        # second moment is half the covariance. Twice that moment is the fixed point.
        mxx = 2.0 * float((ww * dx * dx).sum() / s)
        myy = 2.0 * float((ww * dy * dy).sum() / s)
        mxy = 2.0 * float((ww * dx * dy).sum() / s)
    T = mxx + myy
    if T <= 0:
        return None
    return dict(e1=(mxx - myy) / T, e2=2 * mxy / T, T=T, flux=s, x=cx, y=cy)


def _psf_correct(mom, mom_psf):
    """Ellipticity after subtracting a rounder PSF in the resolution-factor sense. None if unresolved."""
    if mom is None or mom_psf is None:
        return None
    T, Tp = mom['T'], mom_psf['T']
    if T <= Tp * 1.05:
        return None
    fac = T / (T - Tp)
    return dict(e1=mom['e1'] * fac, e2=mom['e2'] * fac, T=T - Tp, flux=mom['flux'])


def shear_image(image, g1, g2):
    """Apply a reduced shear by I_obs(theta) = I(A theta), A = [[1-g1, -g2], [-g2, 1+g1]], about the stamp centre."""
    img = np.asarray(image, float)
    n0, n1 = img.shape
    cy = (n0 - 1) / 2.0
    cx = (n1 - 1) / 2.0
    yy, xx = np.mgrid[0:n0, 0:n1]
    dx, dy = xx - cx, yy - cy
    xs = (1.0 - g1) * dx - g2 * dy + cx
    ys = -g2 * dx + (1.0 + g1) * dy + cy
    return map_coordinates(img, [ys, xs], order=3, mode='constant', cval=0.0)


def _deconvolve(image, psf, eps=1e-3):
    F = np.fft.fft2(np.fft.ifftshift(image))
    P = np.fft.fft2(np.fft.ifftshift(psf))
    power = np.abs(P) ** 2
    W = np.conj(P) / (power + eps * float(power.max()))
    return np.fft.fftshift(np.fft.ifft2(F * W).real)


def gaussian_psf(size, sigma):
    """Normalised round Gaussian PSF on a stamp of the given size. sigma is in pixels."""
    c = (size - 1) / 2.0
    yy, xx = np.mgrid[0:size, 0:size]
    p = np.exp(-0.5 * ((xx - c) ** 2 + (yy - c) ** 2) / max(sigma, 0.3) ** 2)
    return p / p.sum()


def _ellipticity_of(image, psf):
    mom = adaptive_moments(image)
    if mom is None:
        return None
    if psf is None:
        return dict(e1=mom['e1'], e2=mom['e2'], T=mom['T'], T_obs=mom['T'], T_psf=0.0,
                    resolved=True, flux=mom['flux'])
    mom_p = adaptive_moments(psf)
    corr = _psf_correct(mom, mom_p)
    if corr is None:
        tp = mom_p['T'] if mom_p is not None else np.nan
        return dict(e1=np.nan, e2=np.nan, T=np.nan, T_obs=mom['T'], T_psf=tp,
                    resolved=False, flux=mom['flux'])
    out = dict(corr)
    out['T_obs'] = mom['T']
    out['T_psf'] = mom_p['T']
    out['resolved'] = True
    return out


def dilate_kernel(psf, factor):
    """Widen a centred kernel by a linear factor. P(theta) becomes P(theta/factor) / factor^2."""
    psf = np.asarray(psf, float)
    factor = float(factor)
    if factor <= 0:
        raise ValueError('PSF dilation factor must be positive')
    if abs(factor - 1.0) < 1e-6:
        out = np.clip(psf, 0.0, None)
        return out / out.sum()
    n0, n1 = psf.shape
    cy = (n0 - 1) / 2.0
    cx = (n1 - 1) / 2.0
    yy, xx = np.mgrid[0:n0, 0:n1]
    ys = (yy - cy) / factor + cy
    xs = (xx - cx) / factor + cx
    out = map_coordinates(psf, [ys, xs], order=3, mode='constant', cval=0.0)
    out = np.clip(out, 0.0, None)
    s = float(out.sum())
    if s <= 0:
        raise ValueError('dilated PSF has no flux')
    return out / s


def _side_record(measurement):
    """Compact sheared-stamp record for a selection cut."""
    if measurement is None:
        return None
    return dict(e1=measurement.get('e1', np.nan), e2=measurement.get('e2', np.nan),
                T=measurement.get('T_obs', measurement.get('T', np.nan)),
                T_psf=measurement.get('T_psf', 0.0),
                resolved=bool(measurement.get('resolved', False)))


def metacal_stamp(image, psf=None, dg=0.01, dilate=1.1):
    """Shear response of one stamp.

    With a PSF, the counterfactual image is Wiener-deconvolved, sheared, and reconvolved
    with the PSF widened by ``dilate`` (a linear scale, Sheldon & Huff 2017). The
    ellipticity that is divided by the response is the measurement on that reconvolved
    stamp at zero applied shear, with the widened PSF. The response is
    R = (e(+dg) - e(-dg)) / (2 dg). g = e / R.
    sides holds the sheared measurements used for a selection response. The per-object
    shear does not include that response.
    Flag bits: 1 moments failed, 2 unresolved or response below 0.2, 4 non-finite.
    """
    def _blank(flag):
        return dict(e1=np.nan, e2=np.nan, g1=np.nan, g2=np.nan, R1=np.nan, R2=np.nan, flag=flag, sides={})

    img = np.asarray(image, float)
    img = img - np.nanmedian(img)
    img = np.where(np.isfinite(img), img, 0.0)
    wide = None
    if psf is not None:
        psf = np.clip(np.asarray(psf, float), 0.0, None)
        psf = psf / psf.sum()
        wide = dilate_kernel(psf, dilate)

    def measure(g1, g2):
        if psf is None:
            return _ellipticity_of(shear_image(img, g1, g2), None)
        deconv = _deconvolve(img, psf)
        sheared = shear_image(deconv, g1, g2)
        recon = fftconvolve(sheared, wide, mode='same')
        return _ellipticity_of(recon, wide)

    def usable(m):
        return m is not None and m.get('resolved', False) and np.isfinite(m.get('e1', np.nan))

    if psf is None:
        e = _ellipticity_of(img, None)
        if e is None:
            return _blank(1)
    else:
        mom = adaptive_moments(img)
        if mom is None:
            return _blank(1)
        if _psf_correct(mom, adaptive_moments(psf)) is None:
            return _blank(2)
        e = measure(0.0, 0.0)
        if not usable(e):
            return _blank(2)

    raw = {'+1': measure(dg, 0.0), '-1': measure(-dg, 0.0), '+2': measure(0.0, dg), '-2': measure(0.0, -dg)}
    sides = {key: _side_record(val) for key, val in raw.items()}
    ep, em, fp, fm = raw['+1'], raw['-1'], raw['+2'], raw['-2']
    if not all(usable(m) for m in (ep, em, fp, fm)):
        return dict(e1=float(e['e1']), e2=float(e['e2']), g1=np.nan, g2=np.nan,
                    R1=np.nan, R2=np.nan, flag=2, sides=sides)
    R1 = (ep['e1'] - em['e1']) / (2 * dg)
    R2 = (fp['e2'] - fm['e2']) / (2 * dg)
    flag = 0
    if R1 < 0.2 or R2 < 0.2:
        flag = 2
        g1 = g2 = np.nan
    else:
        g1 = e['e1'] / R1
        g2 = e['e2'] / R2
        if not (np.isfinite(g1) and np.isfinite(g2)):
            flag = 4
    return dict(e1=float(e['e1']), e2=float(e['e2']), g1=float(g1) if np.isfinite(g1) else np.nan,
                g2=float(g2) if np.isfinite(g2) else np.nan, R1=float(R1), R2=float(R2),
                flag=int(flag), sides=sides)


def measure_stamps(image, x_pix, y_pix, psf_fwhm_pix=0.0, stamp=32, dg=0.01):
    """Metacalibration at catalog positions. x_pix, y_pix are 0-based pixel coordinates.

    psf_fwhm_pix <= 0 measures the stamp with no PSF model. Returns a list of metacal dicts.
    The PSF in this function is one round Gaussian, and the response is the per-object
    diagonal response. METACAL_SCOPE names what this step leaves to a survey calibration.
    """
    img = np.asarray(image, float)
    half = max(int(stamp) // 2, 5)
    # Odd stamp: the rounded pixel sits on the centre, which matches the FFT origin.
    size = half * 2 + 1
    psf = None
    if psf_fwhm_pix and psf_fwhm_pix > 0:
        psf = gaussian_psf(size, psf_fwhm_pix / (2.0 * math.sqrt(2.0 * math.log(2.0))))
    out = []
    for x, y in zip(np.asarray(x_pix, float), np.asarray(y_pix, float)):
        ix = int(round(float(x)))
        iy = int(round(float(y)))
        if ix - half < 0 or iy - half < 0 or ix + half >= img.shape[1] or iy + half >= img.shape[0]:
            out.append(dict(e1=np.nan, e2=np.nan, g1=np.nan, g2=np.nan, R1=np.nan, R2=np.nan, flag=1, sides={}))
            continue
        cut = img[iy - half:iy + half + 1, ix - half:ix + half + 1].copy()
        out.append(metacal_stamp(cut, psf, dg=dg))
    return out


def _passes_selection(side, min_resolution, emax):
    """True when a sheared stamp passes the resolution and |e| cuts. A non-positive cut is off."""
    if not side or not side.get('resolved'):
        return False
    T = side.get('T', np.nan)
    Tp = side.get('T_psf', 0.0)
    if np.isfinite(T) and np.isfinite(Tp) and T > Tp:
        resolution = 1.0 - Tp / T
    else:
        resolution = 0.0
    if min_resolution > 0 and resolution < min_resolution:
        return False
    if emax > 0:
        e1, e2 = side.get('e1', np.nan), side.get('e2', np.nan)
        if not (np.isfinite(e1) and np.isfinite(e2)) or math.hypot(float(e1), float(e2)) > emax:
            return False
    return True


def selection_response(records, dg, min_resolution=0.0, emax=0.0):
    """Selection response of a stamp catalogue.

    The cut is applied to each sheared stamp. The mean ellipticity is the zero-shear
    measurement of the objects that pass that cut. The same objects on both sides give
    R_s = 0. The ensemble estimator is <e> / (<R_gamma> + R_s). Per-object g = e / R_gamma
    is left unchanged; the mass step applies <R> / (<R> + R_s) when the caller supplies it.
    """
    dg = float(dg)

    def one(component, plus_key, minus_key):
        def mean_on(key):
            vals = []
            for rec in records:
                side = (rec.get('sides') or {}).get(key)
                val = rec.get(component, np.nan)
                if _passes_selection(side, min_resolution, emax) and np.isfinite(val):
                    vals.append(float(val))
            if not vals:
                return None, 0
            return float(np.mean(vals)), len(vals)

        plus, n_plus = mean_on(plus_key)
        minus, n_minus = mean_on(minus_key)
        if plus is None or minus is None or dg == 0:
            return 0.0, n_plus, n_minus
        return (plus - minus) / (2.0 * dg), n_plus, n_minus

    rs1, n1p, n1m = one('e1', '+1', '-1')
    rs2, n2p, n2m = one('e2', '+2', '-2')
    r1 = [rec['R1'] for rec in records if rec.get('flag') == 0 and np.isfinite(rec.get('R1', np.nan))]
    r2 = [rec['R2'] for rec in records if rec.get('flag') == 0 and np.isfinite(rec.get('R2', np.nan))]
    mean_r1 = float(np.mean(r1)) if r1 else 0.0
    mean_r2 = float(np.mean(r2)) if r2 else 0.0
    return dict(selection_R1=float(rs1), selection_R2=float(rs2),
                mean_R1=mean_r1, mean_R2=mean_r2,
                mean_response=float(0.5 * (mean_r1 + mean_r2)),
                n_pass_g1_plus=int(n1p), n_pass_g1_minus=int(n1m),
                n_pass_g2_plus=int(n2p), n_pass_g2_minus=int(n2m))


def apply_selection_response(g1, g2, mean_response, selection_r1=0.0, selection_r2=0.0):
    """Scale observed shear by <R> / (<R> + R_s). A zero selection response leaves g unchanged."""
    mean_response = float(mean_response)

    def one(g, rs):
        rs = float(rs)
        arr = np.asarray(g, float)
        if mean_response <= 0 or rs == 0.0:
            return arr, 1.0
        den = mean_response + rs
        if den <= 0.05:
            return np.full(np.shape(arr), np.nan), float('nan')
        fac = mean_response / den
        return arr * fac, float(fac)

    o1, f1 = one(g1, selection_r1)
    o2, f2 = one(g2, selection_r2)
    return o1, o2, f1, f2


def _sersic_b(n):
    """Ciotti & Bertin approximation for the Sersic half-light coefficient."""
    n = float(n)
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n * n)


def draw_galaxy(size, g1, g2, kind='gaussian', re=4.0, n_sersic=1.0, psf_sigma=0.0, noise=0.0, seed=1):
    """Round Gaussian or Sersic profile, sheared, with an optional Gaussian PSF and noise.

    This is an internal drawing. It is not a survey image simulation.
    """
    size = int(size)
    if size % 2 == 0:
        size += 1
    kind = kind.lower()
    if kind not in ('gaussian', 'sersic'):
        raise ValueError('galaxy kind must be gaussian or sersic')
    c = (size - 1) / 2.0
    yy, xx = np.mgrid[0:size, 0:size]
    rr = np.maximum(np.hypot(xx - c, yy - c), 1e-3)
    re = max(float(re), 0.5)
    if kind == 'gaussian':
        img = np.exp(-0.5 * (rr / re) ** 2)
    else:
        n = max(float(n_sersic), 0.3)
        img = np.exp(-_sersic_b(n) * np.power(rr / re, 1.0 / n))
    sheared = shear_image(img, float(g1), float(g2))
    psf = None
    if psf_sigma and float(psf_sigma) > 0:
        psf = gaussian_psf(size, float(psf_sigma))
        sheared = fftconvolve(sheared, psf, mode='same')
    if noise and float(noise) > 0:
        rng = np.random.default_rng(seed)
        sheared = sheared + rng.normal(0.0, float(noise), sheared.shape)
    return sheared, psf


def calibrate_shear(kind='gaussian', g_values=None, re=4.0, n_sersic=1.0, size=49,
                    psf_sigma=0.0, noise=0.0, seed=1, dg=0.01):
    """Fit g_meas = (1+m) g_true + c on galaxies drawn by this module.

    The result is an internal check of the shape measurement on simple profiles.
    It is not a survey multiplicative bias, and the mass step does not apply it.
    """
    if g_values is None:
        g_values = (-0.06, -0.03, 0.0, 0.03, 0.06)
    rows = []
    n_fail = 0
    for g in g_values:
        for axis in (1, 2):
            g1t = float(g) if axis == 1 else 0.0
            g2t = float(g) if axis == 2 else 0.0
            image, psf = draw_galaxy(size, g1t, g2t, kind=kind, re=re, n_sersic=n_sersic,
                                     psf_sigma=psf_sigma, noise=noise, seed=seed + len(rows))
            rec = metacal_stamp(image, psf, dg=dg)
            if rec['flag'] != 0 or not np.isfinite(rec['g1']) or not np.isfinite(rec['g2']):
                n_fail += 1
                continue
            rows.append((g1t, g2t, float(rec['g1']), float(rec['g2'])))
    note = 'internal simple-galaxy simulation, not a survey calibration'
    if len(rows) < 4:
        return dict(success=False, m=None, c1=None, c2=None, n=len(rows), n_fail=n_fail,
                    kind=kind, note=note, metacal_scope=METACAL_SCOPE)
    design = []
    observed = []
    for g1t, g2t, g1m, g2m in rows:
        design.append([g1t, 1.0, 0.0])
        observed.append(g1m)
        design.append([g2t, 0.0, 1.0])
        observed.append(g2m)
    coef, _res, _rank, _sv = np.linalg.lstsq(np.asarray(design, float), np.asarray(observed, float), rcond=None)
    one_m, c1, c2 = [float(v) for v in coef]
    pred1 = one_m * np.array([r[0] for r in rows]) + c1
    pred2 = one_m * np.array([r[1] for r in rows]) + c2
    rms = float(np.sqrt(np.mean((pred1 - np.array([r[2] for r in rows])) ** 2
                                + (pred2 - np.array([r[3] for r in rows])) ** 2)))
    return dict(success=True, m=one_m - 1.0, c1=c1, c2=c2, rms=rms, n=len(rows), n_fail=n_fail,
                kind=kind, n_sersic=float(n_sersic) if kind == 'sersic' else None,
                psf_sigma_pix=float(psf_sigma), noise=float(noise), note=note,
                metacal_scope=METACAL_SCOPE)
