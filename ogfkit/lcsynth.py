"""Synthetic single-band transient / variable light curves with known class (flux zero point 25, AB).

Classes: SNIa, SNIbc, SNII, fast, AGN, variable, static.  The families are parametric caricatures (Bazin-type rise/decline for the SNe, a
plateau + drop for SN II, a damped random walk for AGN, a harmonic series for variable stars, constant flux for failed subtractions);
a classifier trained on them is only as good as these templates are representative (see docs/lightcurves.md).
"""
import numpy as np

CLASSES = ['SNIa', 'SNIbc', 'SNII', 'fast', 'AGN', 'variable', 'static']
SN_CLASSES = ['SNIa', 'SNIbc', 'SNII']
ZP = 25.0


def mag2flux(m):
    return 10.0 ** (-0.4 * (np.asarray(m, float) - ZP))


def flux2mag(f):
    with np.errstate(divide='ignore', invalid='ignore'):
        return ZP - 2.5 * np.log10(np.asarray(f, float))


def _bazin(t, A, t0, tr, tf):
    return A * np.exp(-(t - t0) / tf) / (1.0 + np.exp(-(t - t0) / tr))


def model_flux(cls, t, rng, peak_mag, z, t0):
    """Noise-free flux at times t (days) for one object; returns (flux, params dict)."""
    s = 1.0 + z                                    # time dilation
    A = mag2flux(peak_mag)
    p = dict(cls=cls, peak_mag=peak_mag, z=z, t0=t0)
    if cls == 'SNIa':
        tr, tf = rng.uniform(2.0, 4.5) * s, rng.uniform(16.0, 28.0) * s
        x = t - t0
        # Bazin shape normalised to its own peak, plus a weak secondary shoulder
        peak = _bazin(np.array([tr * np.log(tf / tr - 1.0 + 1e-9)]), 1.0, 0.0, tr, tf)[0] if tf > tr else 1.0
        f = _bazin(x, 1.0, 0.0, tr, tf) / peak
        f = f + 0.15 * np.exp(-0.5 * ((x - 25 * s) / (6 * s)) ** 2) * rng.uniform(0, 1)
        p.update(tau_r=tr, tau_f=tf)
    elif cls == 'SNIbc':
        tr, tf = rng.uniform(3.5, 7.0) * s, rng.uniform(9.0, 20.0) * s
        x = t - t0
        peak = _bazin(np.array([tr * np.log(tf / tr - 1.0 + 1e-9)]), 1.0, 0.0, tr, tf)[0] if tf > tr else 1.0
        f = _bazin(x, 1.0, 0.0, tr, tf) / peak
        p.update(tau_r=tr, tau_f=tf)
    elif cls == 'SNII':
        tr = rng.uniform(4.0, 11.0) * s
        kind = rng.rand()
        x = t - t0
        rise = 1.0 / (1.0 + np.exp(-x / tr))
        if kind < 0.6:                            # plateau, then drop
            tp = rng.uniform(60, 110) * s; slope = rng.uniform(0.0, 0.012)
            plat = 10 ** (-0.4 * slope * np.clip(x, 0, None))
            drop = 1.0 / (1.0 + np.exp((x - tp) / (4.0 * s)))
            tail = 0.03 * np.exp(-(x - tp) / (50.0 * s)) * (x > tp)
            f = rise * (plat * drop + tail)
            p.update(tau_r=tr, plateau=tp)
        else:                                      # linear decline (II-L)
            sl = rng.uniform(0.02, 0.05) / s
            f = rise * 10 ** (-0.4 * sl * np.clip(x, 0, None))
            p.update(tau_r=tr, decline=sl)
    elif cls == 'fast':
        tr, tf = rng.uniform(0.3, 1.5), rng.uniform(0.8, 4.0)
        x = t - t0
        peak = _bazin(np.array([tr * np.log(tf / tr - 1.0 + 1e-9)]), 1.0, 0.0, tr, tf)[0] if tf > tr else 1.0
        f = _bazin(x, 1.0, 0.0, tr, tf) / peak
        p.update(tau_r=tr, tau_f=tf)
    elif cls == 'AGN':
        tau = rng.uniform(40, 300); sig = rng.uniform(0.08, 0.35)
        ts = np.sort(t); dt = np.diff(ts, prepend=ts[0])
        x = np.zeros(len(ts)); x[0] = rng.randn() * sig
        for i in range(1, len(ts)):
            a = np.exp(-dt[i] / tau)
            x[i] = a * x[i - 1] + np.sqrt(1 - a * a) * sig * rng.randn()
        order = np.argsort(np.argsort(t))
        f = 10 ** (-0.4 * x[order])
        p.update(tau=tau, sigma=sig)
    elif cls == 'variable':
        P = 10 ** rng.uniform(np.log10(0.2), np.log10(12.0)); amp = rng.uniform(0.08, 0.8)
        ph = rng.rand() * 2 * np.pi; h2 = rng.uniform(0, 0.5); ph2 = rng.rand() * 2 * np.pi
        w = 2 * np.pi * t / P
        mm = amp * (np.sin(w + ph) + h2 * np.sin(2 * w + ph2)) / (1 + h2)
        f = 10 ** (-0.4 * mm)
        p.update(period=P, amp=amp)
    else:                                          # static
        f = np.ones_like(t)
    return A * f, p


def make_lc(cls, rng, mlim=None, baseline=None, n=None, host=True):
    """-> dict(t, flux, err, cls, params, host_offset_re (or None)).  Error: sky-limited sigma from the 5-sigma depth plus 3 % of the flux."""
    baseline = baseline or 10 ** rng.uniform(np.log10(15.0), np.log10(220.0))
    n = n or int(rng.randint(6, 41))
    mlim = mlim if mlim is not None else rng.uniform(23.8, 26.3)
    t = np.sort(rng.rand(n) * baseline)
    z = rng.uniform(0.03, 0.8) if cls in ('SNIa', 'SNIbc', 'SNII') else 0.0
    peak = rng.uniform(mlim - 3.8, mlim - 0.8) if cls != 'static' else rng.uniform(mlim - 2.5, mlim - 0.2)
    if cls in ('SNIa', 'SNIbc', 'SNII', 'fast'):
        t0 = rng.uniform(-15.0, baseline * 0.9)
    else:
        t0 = 0.0
    f, p = model_flux(cls, t, rng, peak, z, t0)
    sig_sky = mag2flux(mlim) / 5.0
    err = np.sqrt(sig_sky ** 2 + (0.03 * np.abs(f)) ** 2)
    obs = f + rng.randn(n) * err
    ho = None
    if host:
        if cls == 'AGN':
            ho = abs(rng.randn()) * 0.12
        elif cls in ('SNIa', 'SNIbc', 'SNII'):
            ho = rng.gamma(2.2, 0.45) if cls != 'SNIbc' else rng.gamma(2.8, 0.30)      # disks; Ibc/II more concentrated
        elif cls == 'variable':
            ho = None                                  # foreground star: no host
        else:
            ho = rng.gamma(2.0, 0.6)
    return dict(t=t, flux=obs, err=err, cls=cls, params=p, host_offset_re=ho, mlim=mlim)


def make_set(n_per_class, seed, classes=CLASSES, **kw):
    rng = np.random.RandomState(seed)
    out = []
    for c in classes:
        for _ in range(n_per_class):
            out.append(make_lc(c, rng, **kw))
    return out
