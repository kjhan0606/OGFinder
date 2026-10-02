"""Strong-lens modelling engine (numpy / scipy only; lenstronomy is an optional cross-check, see ``lenstronomy_alpha``).

Conventions (all angles in arcsec on the tangent plane, x to the right and y up in the *pixel frame* of the image; no sky rotation is applied):
  * SIE: convergence kappa = theta_E / (2 sqrt(q x'^2 + y'^2 / q)), x' along the mass major axis; ``phi`` is the position angle of the MAJOR axis,
    degrees counter-clockwise from +x.  theta_E is the area-equivalent Einstein radius (the tangential critical curve of the lone SIE has area pi theta_E^2).
  * SHEAR: gamma1 = gamma cos 2phi, gamma2 = gamma sin 2phi (phi = direction of the stretching axis of the shear, from +x, CCW).
  * NFW: spherical, kappa_s and scale radius theta_s (alpha = 4 kappa_s theta_s h(x)/x).
  * PIEMD: spherical dual pseudo-isothermal sphere, kappa = b0/2 (1/sqrt(a^2+R^2) - 1/sqrt(s^2+R^2)), alpha = b0/R (sqrt(a^2+R^2) - a - sqrt(s^2+R^2) + s);
    used for cluster-member perturbers scaled with luminosity (b0, s ~ L^0.5, i.e. sigma ~ L^0.25).
  * Source plane: beta = theta - alpha(theta).   Fermat potential phi = 0.5 |theta - beta|^2 - psi(theta) (arcsec^2).
"""
import math

import numpy as np

ARCSEC = math.pi / 648000.0
C_KMS = 299792.458
G_MSUN = 4.30091e-3 * 1e-6  # G in Mpc (km/s)^2 / Msun -> used with D in Mpc: see mass_E


# ----------------------------------------------------------------------------------------------------------------- components
class Component:
    kind = ''

    def alpha(self, x, y):
        raise NotImplementedError

    def psi(self, x, y):
        raise NotImplementedError

    def kappa(self, x, y, h=1e-3):
        ax1, ay1 = self.alpha(x + h, y)
        ax0, ay0 = self.alpha(x - h, y)
        ax3, ay3 = self.alpha(x, y + h)
        ax2, ay2 = self.alpha(x, y - h)
        return 0.5 * ((ax1 - ax0) / (2 * h) + (ay3 - ay2) / (2 * h))


class SIE(Component):
    kind = 'SIE'

    def __init__(self, theta_E, q=1.0, phi=0.0, x0=0.0, y0=0.0):
        self.theta_E, self.q, self.phi, self.x0, self.y0 = float(theta_E), float(q), float(phi), float(x0), float(y0)

    def _rot(self, x, y):
        c, s = math.cos(math.radians(self.phi)), math.sin(math.radians(self.phi))
        dx, dy = x - self.x0, y - self.y0
        return c * dx + s * dy, -s * dx + c * dy

    def alpha(self, x, y):
        q = min(max(self.q, 0.05), 1.0)
        xp, yp = self._rot(np.asarray(x, float), np.asarray(y, float))
        c, s = math.cos(math.radians(self.phi)), math.sin(math.radians(self.phi))
        if q > 1 - 1e-6:
            r = np.hypot(xp, yp) + 1e-30
            axp, ayp = self.theta_E * xp / r, self.theta_E * yp / r
        else:
            b = self.theta_E * math.sqrt(q)          # Keeton b: kappa = b / (2 sqrt(q^2 x'^2 + y'^2))
            f = math.sqrt(1 - q * q)
            psi = np.sqrt(q * q * xp * xp + yp * yp) + 1e-30
            axp = b / f * np.arctan(f * xp / psi)
            ayp = b / f * np.arctanh(np.clip(f * yp / psi, -1 + 1e-14, 1 - 1e-14))
        return c * axp - s * ayp, s * axp + c * ayp

    def psi(self, x, y):
        ax, ay = self.alpha(x, y)
        return (np.asarray(x) - self.x0) * ax + (np.asarray(y) - self.y0) * ay   # homogeneous degree-1 potential (no core)


class Shear(Component):
    kind = 'SHEAR'

    def __init__(self, gamma=0.0, phi=0.0, x0=0.0, y0=0.0):
        self.gamma, self.phi, self.x0, self.y0 = float(gamma), float(phi), float(x0), float(y0)

    def alpha(self, x, y):
        g1, g2 = self.gamma * math.cos(2 * math.radians(self.phi)), self.gamma * math.sin(2 * math.radians(self.phi))
        dx, dy = np.asarray(x, float) - self.x0, np.asarray(y, float) - self.y0
        return g1 * dx + g2 * dy, g2 * dx - g1 * dy

    def psi(self, x, y):
        g1, g2 = self.gamma * math.cos(2 * math.radians(self.phi)), self.gamma * math.sin(2 * math.radians(self.phi))
        dx, dy = np.asarray(x, float) - self.x0, np.asarray(y, float) - self.y0
        return 0.5 * g1 * (dx * dx - dy * dy) + g2 * dx * dy


class _Circular(Component):
    """circular profile with tabulated potential: subclasses define alpha_r(r)"""
    x0 = y0 = 0.0
    _tab = None

    def alpha_r(self, r):
        raise NotImplementedError

    def alpha(self, x, y):
        dx, dy = np.asarray(x, float) - self.x0, np.asarray(y, float) - self.y0
        r = np.hypot(dx, dy) + 1e-30
        a = self.alpha_r(r)
        return a * dx / r, a * dy / r

    def psi(self, x, y):
        if self._tab is None:
            rmax = 400.0
            rr = np.concatenate([[0.0], np.geomspace(1e-5, rmax, 6000)])
            a = self.alpha_r(np.maximum(rr, 1e-30))
            a[0] = 0.0
            psi = np.concatenate([[0.0], np.cumsum(0.5 * (a[1:] + a[:-1]) * np.diff(rr))])
            self._tab = (rr, psi)
        rr, psi = self._tab
        r = np.hypot(np.asarray(x, float) - self.x0, np.asarray(y, float) - self.y0)
        return np.interp(r, rr, psi)


class NFW(_Circular):
    kind = 'NFW'

    def __init__(self, kappa_s, theta_s, x0=0.0, y0=0.0):
        self.kappa_s, self.theta_s, self.x0, self.y0 = float(kappa_s), float(theta_s), float(x0), float(y0)

    def alpha_r(self, r):
        x = np.maximum(np.asarray(r, float) / self.theta_s, 1e-8)
        h = np.empty_like(x)
        lo, hi = x < 1 - 1e-6, x > 1 + 1e-6
        mid = ~(lo | hi)
        h[lo] = np.log(x[lo] / 2) + np.arccosh(1 / x[lo]) / np.sqrt(1 - x[lo] ** 2)
        h[hi] = np.log(x[hi] / 2) + np.arccos(1 / x[hi]) / np.sqrt(x[hi] ** 2 - 1)
        h[mid] = np.log(0.5) + 1.0
        return 4 * self.kappa_s * self.theta_s * h / x

    def kappa_r(self, r):
        x = np.maximum(np.asarray(r, float) / self.theta_s, 1e-8)
        F = np.empty_like(x)
        lo, hi = x < 1 - 1e-6, x > 1 + 1e-6
        mid = ~(lo | hi)
        F[lo] = np.arccosh(1 / x[lo]) / np.sqrt(1 - x[lo] ** 2)
        F[hi] = np.arccos(1 / x[hi]) / np.sqrt(x[hi] ** 2 - 1)
        out = np.empty_like(x)
        out[~mid] = 2 * self.kappa_s * (1 - F[~mid]) / (x[~mid] ** 2 - 1)
        out[mid] = 2 * self.kappa_s / 3.0
        return out


class PIEMD(_Circular):
    kind = 'PIEMD'

    def __init__(self, b0, a=0.05, s=30.0, x0=0.0, y0=0.0):
        self.b0, self.a, self.s, self.x0, self.y0 = float(b0), float(a), float(s), float(x0), float(y0)

    def alpha_r(self, r):
        r = np.asarray(r, float)
        return self.b0 / r * (np.sqrt(self.a ** 2 + r * r) - self.a - np.sqrt(self.s ** 2 + r * r) + self.s)

    def kappa_r(self, r):
        r = np.asarray(r, float)
        return 0.5 * self.b0 * (1 / np.sqrt(self.a ** 2 + r * r) - 1 / np.sqrt(self.s ** 2 + r * r))


KINDS = {'SIE': SIE, 'SHEAR': Shear, 'NFW': NFW, 'PIEMD': PIEMD}


def make_component(d):
    d = dict(d)
    k = d.pop('type').upper()
    return KINDS[k](**d)


class LensModel:
    def __init__(self, comps):
        self.comps = [c if isinstance(c, Component) else make_component(c) for c in comps]

    def alpha(self, x, y):
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        ax = np.zeros(np.broadcast(x, y).shape)
        ay = np.zeros_like(ax)
        for c in self.comps:
            a, b = c.alpha(x, y)
            ax = ax + a
            ay = ay + b
        return ax, ay

    def psi(self, x, y):
        return sum(c.psi(x, y) for c in self.comps)

    def beta(self, x, y):
        ax, ay = self.alpha(x, y)
        return np.asarray(x) - ax, np.asarray(y) - ay

    def jacobian(self, x, y, h=2e-4):
        """A = d beta / d theta -> (a11, a12, a21, a22), central differences"""
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        ax1, ay1 = self.alpha(x + h, y)
        ax0, ay0 = self.alpha(x - h, y)
        ax3, ay3 = self.alpha(x, y + h)
        ax2, ay2 = self.alpha(x, y - h)
        return 1 - (ax1 - ax0) / (2 * h), -(ax3 - ax2) / (2 * h), -(ay1 - ay0) / (2 * h), 1 - (ay3 - ay2) / (2 * h)

    def detA(self, x, y):
        a11, a12, a21, a22 = self.jacobian(x, y)
        return a11 * a22 - a12 * a21

    def mu(self, x, y):
        return 1.0 / self.detA(x, y)

    def kappa_gamma(self, x, y):
        a11, a12, a21, a22 = self.jacobian(x, y)
        return 1 - 0.5 * (a11 + a22), 0.5 * np.hypot(a22 - a11, a12 + a21)

    # ------------------------------------------------------------------------------------------------ lens equation
    def solve_images(self, bx, by, centre=(0.0, 0.0), half=None, n=401, min_abs_mu=0.0):
        """all images of a point source at (bx, by): grid triangle mapping + Newton refinement.  -> list of dicts (x, y, mu, parity, dt-free)"""
        if half is None:
            half = 4.0 * max(self.einstein_guess(), 0.5)
        cx, cy = centre
        xs = np.linspace(cx - half, cx + half, n)
        ys = np.linspace(cy - half, cy + half, n)
        X, Y = np.meshgrid(xs, ys)
        BX, BY = self.beta(X, Y)
        cand = []
        i00 = (slice(0, -1), slice(0, -1))
        i10 = (slice(0, -1), slice(1, None))
        i01 = (slice(1, None), slice(0, -1))
        i11 = (slice(1, None), slice(1, None))
        for tri in ((i00, i10, i01), (i11, i10, i01)):
            p0 = (BX[tri[0]], BY[tri[0]])
            p1 = (BX[tri[1]], BY[tri[1]])
            p2 = (BX[tri[2]], BY[tri[2]])
            d = (p1[1] - p2[1]) * (p0[0] - p2[0]) + (p2[0] - p1[0]) * (p0[1] - p2[1])
            with np.errstate(divide='ignore', invalid='ignore'):
                l0 = ((p1[1] - p2[1]) * (bx - p2[0]) + (p2[0] - p1[0]) * (by - p2[1])) / d
                l1 = ((p2[1] - p0[1]) * (bx - p2[0]) + (p0[0] - p2[0]) * (by - p2[1])) / d
            l2 = 1 - l0 - l1
            ok = (l0 >= 0) & (l1 >= 0) & (l2 >= 0)
            jj, ii = np.nonzero(ok)
            for j, i in zip(jj, ii):
                t = [(X[tri[k]][j, i], Y[tri[k]][j, i]) for k in range(3)]
                cand.append((l0[j, i] * t[0][0] + l1[j, i] * t[1][0] + l2[j, i] * t[2][0], l0[j, i] * t[0][1] + l1[j, i] * t[1][1] + l2[j, i] * t[2][1]))
        out = []
        for (x, y) in cand:
            r = self.refine(x, y, bx, by)
            if r is None:
                continue
            if any(math.hypot(r[0] - o['x'], r[1] - o['y']) < 2e-3 for o in out):
                continue
            mu = float(self.mu(r[0], r[1]))
            if abs(mu) < min_abs_mu:
                continue
            out.append(dict(x=float(r[0]), y=float(r[1]), mu=mu, parity='min' if mu > 0 and self.detA_trace(r[0], r[1]) > 0 else ('saddle' if mu < 0 else 'max')))
        out.sort(key=lambda o: -abs(o['mu']))
        return out

    def detA_trace(self, x, y):
        a11, a12, a21, a22 = self.jacobian(x, y)
        return a11 + a22

    def refine(self, x, y, bx, by, iters=60, tol=1e-10):
        """Newton iterations for theta - alpha(theta) = beta starting at (x, y); None when it does not converge"""
        x0, y0 = x, y
        for _ in range(iters):
            ax, ay = self.alpha(x, y)
            rx, ry = x - ax - bx, y - ay - by
            if abs(rx) < tol and abs(ry) < tol:
                return float(x), float(y)
            a11, a12, a21, a22 = self.jacobian(x, y)
            d = a11 * a22 - a12 * a21
            if not np.isfinite(d) or abs(d) < 1e-12:
                return None
            dx = (a22 * rx - a12 * ry) / d
            dy = (-a21 * rx + a11 * ry) / d
            st = math.hypot(dx, dy)
            if st > 0.5:
                dx, dy = dx * 0.5 / st, dy * 0.5 / st
            x, y = x - dx, y - dy
            if math.hypot(x - x0, y - y0) > 20:
                return None
        ax, ay = self.alpha(x, y)
        return (float(x), float(y)) if math.hypot(x - ax - bx, y - ay - by) < 1e-6 else None

    def einstein_guess(self):
        t = 0.0
        for c in self.comps:
            if isinstance(c, SIE):
                t += c.theta_E
            elif isinstance(c, PIEMD):
                t += c.b0 * 0.0 + 0
            elif isinstance(c, NFW):
                t += 0.5 * c.theta_s
        return t if t > 0 else 1.0

    # ------------------------------------------------------------------------------------------------ curves
    def critical_curves(self, centre=(0.0, 0.0), half=None, n=801):
        """-> list of dicts: x, y (critical curve polygon), cx, cy (caustic), area (arcsec^2), theta_eff, kind 'tangential'|'radial'"""
        import contourpy
        if half is None:
            half = 4.0 * max(self.einstein_guess(), 0.5)
        xs = np.linspace(centre[0] - half, centre[0] + half, n)
        ys = np.linspace(centre[1] - half, centre[1] + half, n)
        X, Y = np.meshgrid(xs, ys)
        D = self.detA(X, Y)
        D = np.where(np.isfinite(D), D, 1.0)
        gen = contourpy.contour_generator(X, Y, D)
        curves = []
        for line in gen.lines(0.0):
            if len(line) < 8:
                continue
            x, y = line[:, 0], line[:, 1]
            ref = self._refine_curve(x, y)
            x, y = ref
            bx, by = self.beta(x, y)
            area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
            closed = math.hypot(x[0] - x[-1], y[0] - y[-1]) < 3 * half / n
            curves.append(dict(x=x, y=y, cx=bx, cy=by, area=float(area), theta_eff=float(math.sqrt(area / math.pi)), closed=bool(closed)))
        curves.sort(key=lambda c: -c['area'])
        for i, c in enumerate(curves):
            c['kind'] = 'tangential' if i == 0 else 'radial'
        return curves

    def _refine_curve(self, x, y, iters=3):
        """slide contour vertices onto detA = 0 along the local gradient (removes the grid-interpolation error)"""
        x = x.copy()
        y = y.copy()
        h = 1e-3
        for _ in range(iters):
            d = self.detA(x, y)
            gx = (self.detA(x + h, y) - self.detA(x - h, y)) / (2 * h)
            gy = (self.detA(x, y + h) - self.detA(x, y - h)) / (2 * h)
            g2 = gx * gx + gy * gy + 1e-30
            step = d / g2
            ok = np.isfinite(step) & (np.abs(step * np.sqrt(g2)) < 0.05)
            x = np.where(ok, x - step * gx, x)
            y = np.where(ok, y - step * gy, y)
        return x, y

    def einstein_radius(self, centre=(0.0, 0.0), half=None, n=801):
        """area-equivalent radius of the tangential (largest) critical curve, arcsec; nan when there is no closed curve"""
        cc = self.critical_curves(centre, half, n)
        return cc[0]['theta_eff'] if cc else float('nan')

    # ------------------------------------------------------------------------------------------------ time delays
    def fermat(self, x, y, bx, by):
        return 0.5 * ((np.asarray(x) - bx) ** 2 + (np.asarray(y) - by) ** 2) - self.psi(x, y)


# ----------------------------------------------------------------------------------------------------------------- cosmology and physical units
def distances(z_l, z_s, H0=70.0, Om=0.3):
    """angular-diameter distances (Mpc): D_l, D_s, D_ls (flat LCDM via astropy)"""
    from astropy.cosmology import FlatLambdaCDM
    cos = FlatLambdaCDM(H0=H0, Om0=Om)
    try:
        dls = cos.angular_diameter_distance(z_l, z_s).value          # astropy >= 7
    except TypeError:
        dls = cos.angular_diameter_distance_z1z2(z_l, z_s).value
    return cos.angular_diameter_distance(z_l).value, cos.angular_diameter_distance(z_s).value, dls


def sigma_crit(z_l, z_s, H0=70.0, Om=0.3):
    """critical surface density, Msun / Mpc^2"""
    Dl, Ds, Dls = distances(z_l, z_s, H0, Om)
    # c^2 / (4 pi G) with G = 4.30091e-9 Mpc (km/s)^2 / Msun
    return C_KMS ** 2 / (4 * math.pi * 4.30091e-9) * Ds / (Dl * Dls)


def mass_in_radius(theta_arcsec, z_l, z_s, H0=70.0, Om=0.3):
    """mass inside a circle of angular radius theta that has mean convergence 1 (the Einstein mass), Msun"""
    Dl = distances(z_l, z_s, H0, Om)[0]
    R = theta_arcsec * ARCSEC * Dl
    return math.pi * R * R * sigma_crit(z_l, z_s, H0, Om)


def sis_sigma_v(theta_E, z_l, z_s, H0=70.0, Om=0.3):
    """velocity dispersion (km/s) of the SIS with the same Einstein radius"""
    Dl, Ds, Dls = distances(z_l, z_s, H0, Om)
    return C_KMS * math.sqrt(theta_E * ARCSEC * Ds / (4 * math.pi * Dls))


def time_delays(model, images, bx, by, z_l, z_s, H0=70.0, Om=0.3):
    """-> list of delays (days) relative to the first arriving image, same order as ``images`` (dicts with x, y)"""
    Dl, Ds, Dls = distances(z_l, z_s, H0, Om)
    fac = (1 + z_l) * Dl * Ds / Dls * ARCSEC ** 2 * 3.0856775814914e19 / (C_KMS * 86400.0)   # Mpc -> km, seconds -> days
    ph = np.array([float(model.fermat(i['x'], i['y'], bx, by)) for i in images])
    return list((ph - ph.min()) * fac)


# ----------------------------------------------------------------------------------------------------------------- fitting
DEFAULT_PARAMS = {
    'theta_E': (1.0, 0.05, 20.0), 'q': (0.8, 0.2, 1.0), 'phi': (0.0, -360.0, 360.0), 'x0': (0.0, -1e3, 1e3), 'y0': (0.0, -1e3, 1e3),
    'gamma': (0.0, 0.0, 0.5), 'phi_g': (0.0, -360.0, 360.0),
}


def build_sie_shear(p, perturbers=()):
    comps = [SIE(p['theta_E'], p['q'], p['phi'], p['x0'], p['y0']), Shear(p.get('gamma', 0.0), p.get('phi_g', 0.0), p['x0'], p['y0'])]
    comps += list(perturbers)
    return LensModel(comps)


def _local_images(model, obs_xy, bx, by):
    out = np.full((len(obs_xy), 2), np.nan)
    for k, (x, y) in enumerate(obs_xy):
        r = model.refine(x, y, bx, by, iters=40, tol=1e-9)
        if r is not None:
            out[k] = r
    return out


def fit_sie_shear(obs_xy, sigma, centre, perturbers=(), free=None, start=None, fixed=None, fluxes=None, flux_sigma=None, n_starts=8, seed=0, h_extent=None):
    """Image-plane fit of SIE (+ external shear) to observed multiple-image positions (arcsec, 2-D array N x 2).

    free:   names to vary among theta_E q phi gamma phi_g x0 y0 (default: theta_E q phi gamma phi_g; the lens centre is held at ``centre``)
    fixed:  dict of values for non-free parameters
    The source position (bx, by) is always fitted (nuisance).  Residuals are image-plane: observed minus predicted position of the same image, in units of
    ``sigma`` (arcsec, scalar or per image); an image that cannot be reproduced by the model gets a large residual.  Optional ``fluxes`` add a
    magnification-ratio term (ln |mu_i / mu_0| against ln(F_i / F_0), width ``flux_sigma`` dex-free: natural-log units).
    -> dict(params, source, chi2, dof, rms_arcsec, cov, errors, predicted, n_images_predicted, success)
    """
    from scipy.optimize import least_squares
    obs = np.asarray(obs_xy, float).reshape(-1, 2)
    N = len(obs)
    sig = np.broadcast_to(np.asarray(sigma, float), (N,)).copy()
    free = list(free) if free else ['theta_E', 'q', 'phi', 'gamma', 'phi_g']
    base = dict(theta_E=1.0, q=0.8, phi=0.0, gamma=0.0, phi_g=0.0, x0=float(centre[0]), y0=float(centre[1]))
    if fixed:
        base.update(fixed)
    rng = np.random.default_rng(seed)
    lo = {k: DEFAULT_PARAMS[k][1] for k in DEFAULT_PARAMS}
    hi = {k: DEFAULT_PARAMS[k][2] for k in DEFAULT_PARAMS}
    lo['x0'], hi['x0'] = centre[0] - 3, centre[0] + 3
    lo['y0'], hi['y0'] = centre[1] - 3, centre[1] + 3
    rel = np.hypot(obs[:, 0] - centre[0], obs[:, 1] - centre[1])
    th0 = float(np.mean(rel))
    beta0 = (float(np.mean(obs[:, 0] - 0)), float(np.mean(obs[:, 1])))

    def make(pv):
        p = dict(base)
        for k, v in zip(free, pv[:len(free)]):
            p[k] = float(v)
        return p

    def residual_lin(pv):
        """first-order image-plane residual: A_i^-1 (beta_i - beta_mean) / sigma (cheap, used for initial convergence)"""
        p = make(pv)
        m = build_sie_shear(p, perturbers)
        bxs, bys = m.beta(obs[:, 0], obs[:, 1])
        mbx, mby = np.mean(bxs), np.mean(bys)
        a11, a12, a21, a22 = m.jacobian(obs[:, 0], obs[:, 1])
        d = a11 * a22 - a12 * a21
        d = np.where(np.abs(d) < 1e-6, 1e-6, d)
        dbx, dby = bxs - mbx, bys - mby
        rx = (a22 * dbx - a12 * dby) / d
        ry = (-a21 * dbx + a11 * dby) / d
        return np.concatenate([rx / sig, ry / sig])

    def residual_img(pv):
        p = make(pv[:len(free)])
        bx, by = pv[len(free):len(free) + 2]
        m = build_sie_shear(p, perturbers)
        pred = _local_images(m, obs, bx, by)
        r = (obs - pred) / sig[:, None]
        r = np.where(np.isfinite(r), r, 1e3)
        res = [r[:, 0], r[:, 1]]
        if fluxes is not None:
            mu = np.abs(m.mu(obs[:, 0], obs[:, 1]))
            f = np.asarray(fluxes, float)
            fs = flux_sigma if flux_sigma is not None else 0.3
            res.append((np.log(mu / mu[0]) - np.log(f / f[0]))[1:] / fs)
        return np.concatenate(res)

    names = free
    lb = np.array([lo[k] for k in names])
    ub = np.array([hi[k] for k in names])
    starts = []
    for i in range(max(1, n_starts)):
        pv = []
        for k in names:
            if start and k in start:
                v = start[k] if i == 0 else start[k] + rng.normal(0, 0.1 * abs(start[k]) + 0.05)
            elif k == 'theta_E':
                v = th0 * (1 + 0.15 * rng.normal()) if i else th0
            elif k == 'q':
                v = 0.85 if i == 0 else rng.uniform(0.5, 0.98)
            elif k == 'phi':
                v = 30.0 * i + rng.uniform(0, 20)
            elif k == 'gamma':
                v = 0.05 if i == 0 else rng.uniform(0.0, 0.15)
            elif k == 'phi_g':
                v = 45.0 * i
            else:
                v = base[k]
            pv.append(v)
        starts.append(np.clip(pv, lb + 1e-9, ub - 1e-9))
    best = None
    for s0 in starts:
        try:
            r1 = least_squares(residual_lin, s0, bounds=(lb, ub), x_scale=np.maximum(np.abs(s0), 0.05), max_nfev=200)
            pm = make(r1.x)
            m = build_sie_shear(pm, perturbers)
            bxs, bys = m.beta(obs[:, 0], obs[:, 1])
            full0 = np.concatenate([r1.x, [np.mean(bxs), np.mean(bys)]])
            r2 = least_squares(residual_img, full0, bounds=(np.concatenate([lb, [-1e3, -1e3]]), np.concatenate([ub, [1e3, 1e3]])),
                               x_scale=np.maximum(np.abs(full0), 0.05), max_nfev=300, xtol=1e-12, ftol=1e-12)
        except Exception:
            continue
        c = 2 * r2.cost
        if best is None or c < best[0]:
            best = (c, r2)
    if best is None:
        return dict(success=False, message='all starts failed')
    chi2, r2 = best
    p = make(r2.x[:len(free)])
    bx, by = float(r2.x[len(free)]), float(r2.x[len(free) + 1])
    m = build_sie_shear(p, perturbers)
    pred = _local_images(m, obs, bx, by)
    ndata = 2 * N + (N - 1 if fluxes is not None else 0)
    npar = len(free) + 2
    dof = ndata - npar
    cov = None
    errs = {}
    try:
        J = r2.jac
        cov = np.linalg.pinv(J.T @ J)
        sc = max(chi2 / dof, 1.0) if dof > 0 else 1.0
        for i, k in enumerate(free):
            errs[k] = float(math.sqrt(max(cov[i, i], 0.0) * sc))
        errs['bx'] = float(math.sqrt(max(cov[len(free), len(free)], 0) * sc))
        errs['by'] = float(math.sqrt(max(cov[len(free) + 1, len(free) + 1], 0) * sc))
    except Exception:
        pass
    d = np.hypot(obs[:, 0] - pred[:, 0], obs[:, 1] - pred[:, 1])
    return dict(params=p, free=list(free), source=(bx, by), chi2=float(chi2), dof=int(dof), n_data=int(ndata), rms_arcsec=float(np.sqrt(np.nanmean(d ** 2))),
                residuals_arcsec=[float(v) for v in d], predicted=pred.tolist(), errors=errs, success=bool(np.isfinite(chi2)))


def predict(model, bx, by, centre, observed=None, match_tol=0.15, half=None, n=401, min_abs_mu=0.0):
    """full lens-equation solution for the source (bx, by).  Each predicted image gets ``matched`` = index of the observed image within match_tol (or -1)
    -> list (counter-images = the unmatched ones)"""
    ims = model.solve_images(bx, by, centre=centre, half=half, n=n, min_abs_mu=min_abs_mu)
    used = set()
    for im in ims:
        im['matched'] = -1
        if observed is not None and len(observed):
            d = np.hypot(np.asarray(observed)[:, 0] - im['x'], np.asarray(observed)[:, 1] - im['y'])
            k = int(np.argmin(d))
            if d[k] < match_tol and k not in used:
                im['matched'] = k
                used.add(k)
    return ims


# ----------------------------------------------------------------------------------------------------------------- source plane
def ray_trace_image(model, img, pixscale, origin_xy, centre_xy, src_pix=None, src_half=None, mask=None):
    """Map every pixel of ``img`` to the source plane (surface brightness is conserved) and average into source-plane pixels.

    origin_xy: pixel (0-based x, y) that sits at the lens centre (theta = 0 is placed at ``centre_xy`` arcsec); theta = (pix - origin) * pixscale + centre.
    -> dict(src (mean surface brightness), count, x0, y0 (arcsec of src[0,0] centre), pix (arcsec))
    """
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    tx = (xx - origin_xy[0]) * pixscale + centre_xy[0]
    ty = (yy - origin_xy[1]) * pixscale + centre_xy[1]
    bx, by = model.beta(tx, ty)
    ok = np.isfinite(img) & np.isfinite(bx) & np.isfinite(by)
    if mask is not None:
        ok &= mask > 0
    src_pix = src_pix or 0.5 * pixscale
    if src_half is None:
        src_half = float(max(np.percentile(np.abs(bx[ok] - np.median(bx[ok])), 99), np.percentile(np.abs(by[ok] - np.median(by[ok])), 99), 5 * src_pix))
    cx, cy = float(np.median(bx[ok])), float(np.median(by[ok]))
    n = int(2 * src_half / src_pix) + 1
    ix = np.floor((bx[ok] - (cx - src_half)) / src_pix + 0.5).astype(int)
    iy = np.floor((by[ok] - (cy - src_half)) / src_pix + 0.5).astype(int)
    good = (ix >= 0) & (ix < n) & (iy >= 0) & (iy < n)
    s = np.zeros((n, n))
    c = np.zeros((n, n))
    np.add.at(s, (iy[good], ix[good]), img[ok][good])
    np.add.at(c, (iy[good], ix[good]), 1.0)
    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(c > 0, s / c, np.nan)
    return dict(src=mean, count=c, x0=cx - src_half, y0=cy - src_half, pix=src_pix)


def lens_source_map(model, srcmap, shape, pixscale, origin_xy, centre_xy):
    """forward model: surface brightness of the (reconstructed or analytic) source map at beta(theta) for every image pixel"""
    from scipy.ndimage import map_coordinates
    ny, nx = shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    tx = (xx - origin_xy[0]) * pixscale + centre_xy[0]
    ty = (yy - origin_xy[1]) * pixscale + centre_xy[1]
    bx, by = model.beta(tx, ty)
    fx = (bx - srcmap['x0']) / srcmap['pix']
    fy = (by - srcmap['y0']) / srcmap['pix']
    s = np.nan_to_num(srcmap['src'], nan=0.0)
    return map_coordinates(s, [fy, fx], order=1, mode='constant', cval=0.0)


# ----------------------------------------------------------------------------------------------------------------- optional cross-check
def lenstronomy_alpha(model_kind, params, x, y):
    """deflection from lenstronomy for an SIE / SHEAR / NFW cross-check; ImportError when lenstronomy is not installed"""
    from lenstronomy.LensModel.lens_model import LensModel as LM
    if model_kind == 'SIE':
        q = params['q']
        e1 = (1 - q) / (1 + q) * math.cos(2 * math.radians(params['phi']))
        e2 = (1 - q) / (1 + q) * math.sin(2 * math.radians(params['phi']))
        lm = LM(['SIE'])
        kw = [dict(theta_E=params['theta_E'], e1=e1, e2=e2, center_x=params['x0'], center_y=params['y0'])]
    elif model_kind == 'NFW':
        lm = LM(['NFW'])
        kw = [dict(Rs=params['theta_s'], alpha_Rs=float(NFW(params['kappa_s'], params['theta_s']).alpha_r(np.array([params['theta_s']]))[0]),
                   center_x=params.get('x0', 0.0), center_y=params.get('y0', 0.0))]
    else:
        raise ValueError(model_kind)
    return lm.alpha(np.asarray(x), np.asarray(y), kw)
