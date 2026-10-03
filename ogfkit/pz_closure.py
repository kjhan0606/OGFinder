"""Photo-z closure: cross-fitted PIT recalibration fed back into the predictive distributions, and a consistency check between two redshift estimates
(e.g. a template SED fit and the empirical / MDN photo-z).

recalibration.  The PIT values u_i = F_i(z_spec,i) of a calibrated predictive are uniform.  PITRecal learns the empirical CDF G of the training PITs (monotone,
piecewise linear through `nknots` quantile knots, with a floor on its slope) and maps new PITs through it.  `recalibrated_predictive` feeds this back into the model:
p'(z) = p(z) g(F(z)) with g = G' (so that F'(z) = G(F(z)) is the recalibrated CDF), returned as a grid-type `photoz_stats.Predictive`.  `crossfit` estimates
the benefit honestly: each fold is recalibrated with a map learnt on the other folds only.
consistency.  `consistency` compares two estimates z_A +- s_A, z_B +- s_B of the same objects.  Their errors are correlated (same photometry), so the normalised
difference d = (z_A - z_B)/sqrt(s_A^2 + s_B^2) is *not* N(0,1): `calibrate` measures its robust width c on objects with spec-z (or any reference sample), and flags
|d| > nsig c.  With spec-z the flagged / unflagged samples are compared (which estimate is right when they disagree; accuracy of the agreeing subsample).
"""
import math

import numpy as np

from . import photoz_stats as ps


class PITRecal:
    def __init__(self, knots_u, knots_g):
        self.u = np.asarray(knots_u, float)
        self.g = np.asarray(knots_g, float)

    @classmethod
    def fit(cls, pit, nknots=40, min_slope=0.05):
        t = np.sort(np.clip(np.asarray(pit, float), 0.0, 1.0))
        t = t[np.isfinite(t)]
        n = len(t)
        q = np.linspace(0.0, 1.0, nknots + 1)
        u = np.quantile(t, q)                        # PIT values at G = q
        u[0], u[-1] = 0.0, 1.0
        u = np.maximum.accumulate(u)
        # slopes dG/du = dq/du, floored; renormalise so that G(1) = 1
        du = np.maximum(np.diff(u), 1e-9)
        dg = np.diff(q)
        slope = np.maximum(dg / du, min_slope)
        dgn = slope * du
        g = np.concatenate([[0.0], np.cumsum(dgn)])
        g = g / g[-1]
        # merge duplicate u knots (ties) for interpolation
        uu, idx = np.unique(u, return_index=True)
        return cls(uu, g[idx])

    def apply(self, pit):
        return np.interp(np.clip(np.asarray(pit, float), 0, 1), self.u, self.g)

    def slope(self, pit):
        """g = dG/du at the PIT values (piecewise constant between knots)."""
        pit = np.clip(np.asarray(pit, float), 0, 1)
        i = np.clip(np.searchsorted(self.u, pit, side='right') - 1, 0, len(self.u) - 2)
        return (self.g[i + 1] - self.g[i]) / np.maximum(self.u[i + 1] - self.u[i], 1e-12)

    def to_dict(self):
        return dict(u=self.u.tolist(), g=self.g.tolist())

    @classmethod
    def from_dict(cls, d):
        return cls(d['u'], d['g'])


def recalibrated_predictive(pred, recal, zgrid=None):
    """Grid predictive with F'(z) = G(F(z)): p'(z) = p(z) G'(F(z))."""
    zg = np.linspace(0.0 if pred.zmin is None else pred.zmin, 4.0, 2001) if zgrid is None else np.asarray(zgrid, float)
    P = pred.grid(zg)
    C = np.concatenate([np.zeros((pred.n, 1)), np.cumsum(0.5 * (P[:, 1:] + P[:, :-1]) * np.diff(zg), axis=1)], axis=1)
    C = C / np.maximum(C[:, -1:], 1e-12)
    P2 = P * recal.slope(C.ravel()).reshape(C.shape)
    return ps.Predictive('grid', zmin=pred.zmin, zgrid=zg, pdf=P2)


def _subset(pred, idx):
    c = ps.Predictive.__new__(ps.Predictive)
    c.__dict__.update(pred.__dict__)
    for k in ('pi', 'mu', 'sig', 'med', 'slo', 'shi', 'p', 'cg'):
        if k in pred.__dict__:
            setattr(c, k, getattr(pred, k)[idx])
    c.n = len(idx)
    return c


def summary(pred, zspec):
    s, pit, cr, ls, zsc = ps.calibration_summary(pred, zspec)
    cov = {str(c['level']): c['observed'] for c in s['pit']['coverage']}
    return dict(ks_p=s['pit']['ks_p'], shape=s['pit']['shape'], crps=s['crps_mean'], logscore=s['logscore_mean'], coverage=cov, zscore_std=s['zscore']['std']), pit


def crossfit(pred, zspec, k=5, seed=1, nknots=40, zgrid=None):
    """K-fold: recalibrate each fold with a map learnt on the others.  -> dict(before, after_oof, pit_after [out of fold], map [learnt on everything])"""
    zspec = np.asarray(zspec, float)
    n = len(zspec)
    fold = np.random.default_rng(seed).permutation(n) % k
    zg = np.linspace(0.0 if pred.zmin is None else pred.zmin, 4.0, 2001) if zgrid is None else zgrid
    before, pit0 = summary(pred, zspec)
    pit_all = pred.cdf(zspec)
    P = np.empty((n, len(zg)))
    for f in range(k):
        te = np.where(fold == f)[0]
        rc = PITRecal.fit(pit_all[fold != f], nknots)
        sub = _subset(pred, te)
        P[te] = recalibrated_predictive(sub, rc, zg).p
    after_pred = ps.Predictive('grid', zmin=pred.zmin, zgrid=zg, pdf=P)
    after, pit1 = summary(after_pred, zspec)
    return dict(before=before, after_oof=after, pit_after=pit1, map=PITRecal.fit(pit_all, nknots))


# ---------------------------------------------------------------------------------------------------------------- consistency
def _d(zA, sA, zB, sB):
    return (np.asarray(zA, float) - np.asarray(zB, float)) / np.sqrt(np.asarray(sA, float) ** 2 + np.asarray(sB, float) ** 2)


def calibrate(zA, sA, zB, sB):
    """Robust width c of the normalised difference d on a reference sample (use objects with spec-z, or the whole catalogue if the estimates are good)."""
    d = _d(zA, sA, zB, sB)
    d = d[np.isfinite(d)]
    return dict(c=float(1.4826 * np.median(np.abs(d - np.median(d)))), median=float(np.median(d)), n=int(len(d)))


def consistency(zA, sA, zB, sB, c=1.0, nsig=3.0, zspec=None, outlier=0.15):
    """Flags |d - median| > nsig c.  -> dict(flag [bool], d, frac_flagged, plus with zspec: closer-estimate fractions of the flagged objects, outlier fractions of
    the unflagged / flagged sets, accuracy (sigma_NMAD) of the agreeing objects and of the inverse-variance mean of A and B)."""
    zA = np.asarray(zA, float); zB = np.asarray(zB, float); sA = np.asarray(sA, float); sB = np.asarray(sB, float)
    d = _d(zA, sA, zB, sB)
    flag = np.abs(d) > nsig * c
    out = dict(flag=flag, d=d, frac_flagged=float(np.mean(flag)))
    if zspec is not None:
        z = np.asarray(zspec, float)
        dzA = (zA - z) / (1 + z); dzB = (zB - z) / (1 + z)
        wm = (zA / sA ** 2 + zB / sB ** 2) / (1 / sA ** 2 + 1 / sB ** 2)
        dzm = (wm - z) / (1 + z)
        ok = ~flag
        out.update(n=int(len(z)), n_flagged=int(flag.sum()),
                   outlier_A=float(np.mean(np.abs(dzA) > outlier)), outlier_B=float(np.mean(np.abs(dzB) > outlier)),
                   outlier_A_unflagged=float(np.mean(np.abs(dzA[ok]) > outlier)), outlier_B_unflagged=float(np.mean(np.abs(dzB[ok]) > outlier)),
                   outlier_A_flagged=float(np.mean(np.abs(dzA[flag]) > outlier)) if flag.any() else float('nan'),
                   outlier_B_flagged=float(np.mean(np.abs(dzB[flag]) > outlier)) if flag.any() else float('nan'),
                   A_closer_when_flagged=float(np.mean(np.abs(dzA[flag]) < np.abs(dzB[flag]))) if flag.any() else float('nan'),
                   nmad_A=ps.nmad(dzA), nmad_B=ps.nmad(dzB), nmad_A_unflagged=ps.nmad(dzA[ok]), nmad_B_unflagged=ps.nmad(dzB[ok]), nmad_mean_unflagged=ps.nmad(dzm[ok]),
                   outlier_mean_unflagged=float(np.mean(np.abs(dzm[ok]) > outlier)))
    return out
