"""Child process for an explicit rebound and assist propagation.

The parent orbit module does not import these packages. This file imports
them only when it is executed. A missing package or a missing kernel file
exits 2 and writes no result.
"""
import json
import os
import sys

FORCE_MODEL = "ASSIST(DE440+16 asteroids+GR+J2/J3/J4)"


def main():
    try:
        import rebound
        import assist
    except Exception:
        sys.stderr.write(
            "rebound or assist was not found. "
            "The in-tree orbit model is ogfmeas.\n")
        return 2
    from moving.orbit import ephem_paths
    planets, asteroids = ephem_paths()
    if not (os.path.exists(planets) and os.path.exists(asteroids)):
        sys.stderr.write(
            "ephemeris files were not found. "
            "The in-tree orbit model is ogfmeas.\n")
        return 2
    req = json.load(sys.stdin)
    import numpy as np
    states = np.atleast_2d(np.asarray(req["states"], float))
    jd0 = float(req["jd0"])
    jds = np.asarray(req["jds"], float)
    ep = assist.Ephem(planets, asteroids)
    jd_ref = ep.jd_ref
    order = np.argsort(jds)
    out = np.zeros((len(states), len(jds), 6))
    fw = [k for k in order if jds[k] >= jd0]
    bw = [k for k in order[::-1] if jds[k] < jd0]
    # One fresh simulation per direction. Reusing a simulation after the
    # forward sweep corrupts ASSIST bookkeeping.
    for seq in (fw, bw):
        if not seq:
            continue
        sim = rebound.Simulation()
        sim.t = jd0 - jd_ref
        for s_ in states:
            sim.add(x=s_[0], y=s_[1], z=s_[2], vx=s_[3], vy=s_[4], vz=s_[5])
        ex = assist.Extras(sim, ep)
        for k in seq:
            ex.integrate_or_interpolate(jds[k] - jd_ref)
            for i in range(len(states)):
                p = sim.particles[i]
                out[i, k] = (p.x, p.y, p.z, p.vx, p.vy, p.vz)
        del ex, sim
    json.dump({"force_model": FORCE_MODEL, "states": out.tolist()}, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main() or 0)
