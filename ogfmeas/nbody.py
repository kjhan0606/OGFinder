"""Newtonian point-mass leapfrog.

Kick-drift-kick with mutual 1/r^2 forces. This is a small N-body step for
tests and for a point-mass trajectory. It does not read planetary ephemerides
and does not include relativity or zonal harmonics. Solar-system test
particles, with planet positions from the built-in ephemeris, are integrated
by ``ogfmeas.ssbody``.
"""
import numpy as np


def accelerations(positions, masses, G=1.0):
    r = np.asarray(positions, dtype=np.float64)
    m = np.asarray(masses, dtype=np.float64)
    acc = np.zeros_like(r)
    n = m.size
    for i in range(n):
        for j in range(i + 1, n):
            delta = r[j] - r[i]
            dist2 = float(np.dot(delta, delta))
            dist = dist2 ** 0.5
            factor = G / (dist2 * dist + 1e-30)
            acc[i] += m[j] * factor * delta
            acc[j] -= m[i] * factor * delta
    return acc


def propagate(positions, velocities, masses, dt, steps, G=1.0):
    """Advance point masses. ``positions`` and ``velocities`` have shape (n, 3)."""
    r = np.array(positions, dtype=np.float64, copy=True)
    v = np.array(velocities, dtype=np.float64, copy=True)
    m = np.asarray(masses, dtype=np.float64)
    if r.shape != v.shape or r.ndim != 2 or r.shape[1] != 3:
        raise ValueError("positions and velocities must have shape (n, 3)")
    acc = accelerations(r, m, G=G)
    step = float(dt)
    for _ in range(int(steps)):
        v += 0.5 * step * acc
        r += step * v
        acc = accelerations(r, m, G=G)
        v += 0.5 * step * acc
    return r, v
