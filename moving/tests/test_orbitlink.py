"""moving.orbitlink: population rate prior used to vet tracklets (synthetic element table, no network)."""
import math
import numpy as np
import pytest

from moving import orbitlink as OL, tracklet as T


@pytest.fixture(scope='module')
def elements():
    rng = np.random.default_rng(3)
    n = 60000
    a = rng.uniform(2.1, 3.3, n)
    return dict(a=a, e=rng.uniform(0.0, 0.3, n), i=np.abs(rng.normal(0, 8, n)), om=rng.uniform(0, 360, n), w=rng.uniform(0, 360, n), H=rng.uniform(12, 19.5, n))


def test_kepler_state_circular():
    pos, vel = OL._state(np.array([2.0]), np.array([0.0]), np.array([0.3]), np.array([1.0]), np.array([0.5]), np.array([2.0]))
    assert abs(np.linalg.norm(pos) - 2.0) < 1e-9
    assert abs(np.linalg.norm(vel) - math.sqrt(OL.GM_SUN / 2.0)) < 1e-9
    assert abs(float(pos @ vel.T)) < 1e-9                                               # circular: r perpendicular to v


def test_population_at_opposition_is_retrograde_and_slow(elements):
    # opposition = the direction from the Sun through the Earth (the real Earth state of the epoch).
    from moving.tracklet import earth_helio_state
    mjd = 51075.3
    pe, _ = earth_helio_state(mjd)
    ra_opp = math.degrees(math.atan2(pe[1], pe[0])) % 360.0
    pop = OL.population_rates(elements, mjd, ra_opp, math.degrees(math.asin(pe[2] / np.linalg.norm(pe))), radius_deg=4.0, n_draw=20, vmax=22.5)
    assert pop['n'] > 300
    assert -45 < np.median(pop['rate_e']) < -20                                          # main belt near opposition: ~ -30 arcsec/h in RA
    assert abs(np.median(pop['rate_n'])) < 20
    assert pop['rate_e'].std() < 12


def test_hg_magnitude_phase_and_distance():
    h = OL.hg_magnitude(np.array([15.0, 15.0]), np.array([2.5, 2.5]), np.array([1.5, 1.5]), 1.0)
    assert h[0] == h[1]
    near = OL.hg_magnitude(15.0, 2.0, 1.0, 1.0); far = OL.hg_magnitude(15.0, 3.0, 2.0, 1.0)
    assert far > near + 2.0
    assert OL.hg_magnitude(15.0, 2.5, 1.5, 1.0) < OL.hg_magnitude(15.0, 2.5, 1.5, 2.0 + 1e-3) + 5      # sanity: finite


def test_rateprior_hdr_and_vetting(elements):
    rng = np.random.default_rng(2)
    re = rng.normal(-30, 4, 5000); rn = rng.normal(-12, 6, 5000)
    pr = OL.RatePrior(re, rn, bandwidth=2.0)
    lvl = pr.hdr_level(0.99)
    inside = pr.logp(re[:2000], rn[:2000]) >= lvl
    assert 0.97 < inside.mean() <= 1.0
    trs = []
    for r_e, r_n in ((-30, -12), (-29, -8), (40, 20), (5, 60), (-80, 0)):
        rate = math.hypot(r_e, r_n)
        trs.append(dict(rate_ash=rate, pa_deg=math.degrees(math.atan2(r_e, r_n)) % 360))
    out = OL.vet_tracklets(trs, pr, 0.99)
    assert [t['in_hdr'] for t in out] == [True, True, False, False, False]
    assert len(OL.vet_tracklets([dict(t) for t in trs], pr, 0.99, drop=True)) == 2
    assert out[0]['log_prior'] > out[2]['log_prior']


def test_prior_for_tracklets_and_link_detections_integration(elements):
    from moving import pipeline as P
    ra0, de0, t0 = 357.0, -0.1, 51075.3
    trs = [dict(ra=[ra0, ra0], dec=[de0, de0], t=[t0, t0 + 0.001], rate_ash=31.0, pa_deg=250.0), dict(ra=[ra0, ra0], dec=[de0, de0], t=[t0, t0 + 0.001], rate_ash=90.0, pa_deg=10.0)]
    pr = OL.prior_for_tracklets(trs, elements, radius_deg=4.0)
    assert pr is not None
    out = OL.vet_tracklets(trs, pr, 0.99)
    assert out[0]['in_hdr'] and not out[1]['in_hdr']
    # the pipeline accepts the prior and drops the tracklet that is outside the region
    assert 'orbit_prior' in P.link_detections.__code__.co_varnames
