"""Synthetic observations for validating moving.nightlink.  All of this is SYNTHETIC: truth orbits are 2-body (so the dynamics
model of the linker is exact) unless the tracklets come from JPL Horizons ephemerides (n-body truth, see validation script)."""
import numpy as np
from . import kepler as K
from . import nightlink as N


def random_orbits(n, rng, kind="mainbelt"):
    """Heliocentric ecliptic elements -> list of dict(a,e,i,om,w,M0).  kind: mainbelt | neo | mixed."""
    out = []
    for _ in range(n):
        k = kind if kind != "mixed" else rng.choice(["mainbelt", "mainbelt", "mainbelt", "neo", "hilda"])
        if k == "mainbelt":
            a = rng.uniform(2.1, 3.3); e = min(rng.rayleigh(0.1), 0.35); inc = min(rng.rayleigh(8.0), 30.0)
        elif k == "neo":
            a = rng.uniform(1.0, 2.2); e = rng.uniform(0.2, 0.6); inc = min(rng.rayleigh(10.0), 40.0)
        else:
            a = rng.uniform(3.9, 4.1); e = rng.uniform(0.05, 0.3); inc = rng.uniform(0, 20)
        out.append(dict(a=a, e=e, i=inc, om=rng.uniform(0, 360), w=rng.uniform(0, 360), M0=rng.uniform(0, 360), kind=k))
    return out


def orbit_state(el, t_ref=None, t_epoch=None):
    r, v = K.elements_to_state(el["a"], el["e"], el["i"], el["om"], el["w"], el["M0"], ecliptic=False)
    return np.asarray(r, float), np.asarray(v, float)


def random_state(rng, ra_window, dec_window, t_epoch, code=None, kind="mainbelt"):
    """Heliocentric state of a bound object that is inside the sky window at t_epoch: sky position uniform, range from a
    population-like distribution, velocity = prograde near-circular motion about the ecliptic pole with scatter."""
    ra = rng.uniform(*ra_window); de = rng.uniform(*dec_window)
    op, _ = N.observer_helio(t_epoch, code); op = op[0]
    los = N._los(ra, de)
    rho = {"mainbelt": rng.uniform(0.9, 2.6), "neo": rng.uniform(0.05, 0.6), "hilda": rng.uniform(3.0, 4.0)}[kind]
    r = op + rho * los
    zecl = K.R_EQ2ECL.T @ np.array([0.0, 0.0, 1.0])
    tdir = np.cross(zecl, r); tdir /= np.linalg.norm(tdir)
    tilt = np.radians(rng.normal(0, 8.0)); ax = np.cross(tdir, r); ax /= np.linalg.norm(ax)
    tdir = tdir * np.cos(tilt) + ax * np.sin(tilt)
    rn = np.linalg.norm(r)
    vc = np.sqrt(K.MU / rn)
    v = tdir * vc * rng.uniform(0.82, 1.12) + r / rn * vc * rng.normal(0, 0.12)
    if 0.5 * v @ v - K.MU / rn >= 0:
        v *= 0.8
    return r, v


def observe_state(state, t_epoch, times, rng, sig=0.3, code=None):
    op, _ = N.observer_helio(times, code)
    ra, de, rg = N.predict(state[0], state[1], t_epoch, np.asarray(times, float), op)
    cd = np.cos(np.radians(de))
    return ra + rng.normal(0, sig, len(ra)) / 3600.0 / cd, de + rng.normal(0, sig, len(de)) / 3600.0, rg


def make_field(n_obj, nights, rng, sig=0.3, kind="mainbelt", n_per_night=3, spacing_h=0.5, t0=61000.0, p_detect=1.0,
               ra_window=(100.0, 106.0), dec_window=(0.0, 6.0), n_decoy=0, code=None, jitter_h=0.0):
    """Objects placed in the sky window on the first night (bound near-circular orbits, see random_state), observed on the
    given integer nights (days after t0) in `n_per_night` exposures `spacing_h` apart; each tracklet is recovered with
    probability p_detect.  n_decoy random single-night tracklets (same window, rates ~|N(0, 0.2)| deg/day) are added per night.
    Returns Tracklet list with truth = object index (decoys -1).  SYNTHETIC."""
    trks = []
    for oi in range(n_obj):
        st = random_state(rng, ra_window, dec_window, t0 + nights[0], code, kind)
        for nt in nights:
            if rng.random() > p_detect:
                continue
            times = t0 + nt + np.arange(n_per_night) * spacing_h / 24.0 + rng.uniform(-jitter_h, jitter_h) / 24.0
            ra, de, _ = observe_state(st, t0 + nights[0], times, rng, sig, code)
            trks.append(N.Tracklet(times, ra, de, sig, code=code, tid=len(trks), night=int(nt), truth=oi))
    for nt in nights:
        for _ in range(n_decoy):
            times = t0 + nt + np.arange(n_per_night) * spacing_h / 24.0
            ra0 = rng.uniform(*ra_window); de0 = rng.uniform(*dec_window)
            sp = abs(rng.normal(0.0, 0.2)); pa = rng.uniform(0, 2 * np.pi); dt = times - times[0]
            cd = np.cos(np.radians(de0))
            ra = ra0 + sp * np.cos(pa) * dt / cd + rng.normal(0, sig, n_per_night) / 3600 / cd
            de = de0 + sp * np.sin(pa) * dt + rng.normal(0, sig, n_per_night) / 3600
            trks.append(N.Tracklet(times, ra, de, sig, code=code, tid=len(trks), night=int(nt), truth=-1))
    return trks


def _score_groups(trks, groups, min_members):
    by_obj = {}
    for t in trks:
        if t.truth is not None and t.truth >= 0:
            by_obj.setdefault(t.truth, []).append(t.id)
    objs = {o: v for o, v in by_obj.items() if len(v) >= min_members}
    tid2truth = {t.id: t.truth for t in trks}
    pure = wrong = complete = 0; found = set()
    for g in groups:
        tr = [tid2truth[i] for i in g["ids"]]
        if len(set(tr)) == 1 and tr[0] is not None and tr[0] >= 0:
            pure += 1; found.add(tr[0])
            if len(g["ids"]) == len(by_obj.get(tr[0], [])):
                complete += 1
        else:
            wrong += 1
    return dict(n_objects=len(objs), n_groups=len(groups), pure=pure, wrong=wrong, recovered=len(found & set(objs)), complete=complete,
                recall=len(found & set(objs)) / max(len(objs), 1), precision=pure / max(len(groups), 1))


def score(trks, result):
    """Truth scoring (SYNTHETIC/known labels).  `groups` (>= 3 tracklets) are scored against objects with >= 3 tracklets; `pairs` against
    objects with >= 2.  A group is 'pure' when all members are one object, 'wrong' otherwise, 'complete' when it holds all the object's tracklets."""
    sc = _score_groups(trks, result["groups"], 3)
    pr = _score_groups(trks, result.get("pairs", []), 2)
    sc.update(pair_groups=pr["n_groups"], pair_pure=pr["pure"], pair_wrong=pr["wrong"])
    return sc
