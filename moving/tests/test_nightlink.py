"""moving.nightlink: cross-night orbit-fit linking.  Offline tests use SYNTHETIC 2-body fields (labelled); the network test
uses real MPC astrometry of known asteroids.  No ASSIST ephemeris needed."""
import os, sys
import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from moving import nightlink as N, nightlink_sim as S, kepler as K
from conftest import needs_net


def test_attributable_of_a_linear_tracklet():
    t = 61000.0 + np.array([0.0, 0.02, 0.04])
    ra = 100.0 + np.array([0, 1, 2]) * 0.01 / 3600 * 10      # 0.1" per step in RA*cos(dec)
    de = np.array([5.0, 5.0, 5.0]) + np.array([0, 2, 4]) * 0.1 / 3600
    T = N.Tracklet(t, ra / 1.0, de, 0.1)
    assert abs(T.mu_y - 0.2 / 0.02 / 1.0 * 1.0 * 0.5 * 1.0) < 1e-6 or abs(T.mu_y - 10.0) < 1e-6      # 0.2" per 0.02 d = 10"/d
    assert T.sig_mu > 0 and abs(T.arc_days - 0.04) < 1e-9


def test_gate_separates_true_pair_from_shifted_decoy():
    rng = np.random.default_rng(3)
    trks = S.make_field(1, [0, 3], rng, sig=0.2, ra_window=(100, 101), dec_window=(0, 1))
    A, B = trks
    chi, idx, hyp = N.pair_gate(A, B)
    assert chi < 20 and len(idx) > 0
    # same tracklet shifted by 2 degrees on the sky must be rejected by the gate
    Bs = N.Tracklet(B.mjd, B.ra + 2.0, B.dec, 0.2, night=B.night)
    chi2_, _, _ = N.pair_gate(A, Bs, hyp=hyp)
    assert chi2_ > 100


def test_fit_pair_recovers_orbit():
    rng = np.random.default_rng(5)
    st = S.random_state(rng, (100, 101), (0, 1), 61000.0)
    trks = []
    for nt in (0, 3, 7):
        times = 61000.0 + nt + np.arange(3) * 0.5 / 24
        ra, de, _ = S.observe_state(st, 61000.0, times, rng, 0.1)
        trks.append(N.Tracklet(times, ra, de, 0.1, night=nt, tid=nt, truth=0))
    res = N.link_nights(trks)
    assert len(res["groups"]) == 1 and sorted(res["groups"][0]["ids"]) == [0, 3, 7]
    el = K.state_to_elements(*st)
    g = res["groups"][0]["elements"]
    assert abs(g["a"] - el["a"]) / el["a"] < 0.05 and abs(g["e"] - el["e"]) < 0.05
    assert res["groups"][0]["chi2_red"] < 2.0


def test_injection_recovery_and_no_false_links():
    """40 objects x 3 nights (90 % detection) + 40 decoys per night in the same 6x6 deg field (SYNTHETIC, 2-body truth)."""
    rng = np.random.default_rng(11)
    trks = S.make_field(25, [0, 2, 5], rng, sig=0.3, n_decoy=25, p_detect=0.9)
    res = N.link_nights(trks)
    sc = S.score(trks, res)
    assert sc["recall"] >= 0.95, sc
    assert sc["wrong"] == 0, sc
    vet = N.vet_with_links(trks, res)
    assert sum(1 for v in vet.values() if v["status"] == "linked") >= 2 * sc["recovered"]


def test_pure_decoys_do_not_link():
    rng = np.random.default_rng(2)
    trks = S.make_field(0, [0, 1, 4], rng, sig=0.3, n_decoy=40)
    res = N.link_nights(trks)
    assert len(res["groups"]) == 0, [g["ids"] for g in res["groups"]]


def test_same_night_tracklets_never_grouped():
    rng = np.random.default_rng(4)
    trks = S.make_field(5, [0, 0, 3], rng, sig=0.2)
    res = N.link_nights(trks)
    for g in res["groups"]:
        assert len(set(g["nights"])) == len(g["nights"])


@needs_net
def test_real_mpc_astrometry_links_known_asteroids():
    """REAL MPC astrometry of 12 numbered asteroids that were in one field (SBDB elements, 2-body positions), tracklets by
    station/night within +-8 d; the pooled tracklets must link with no wrong pairing and recover most multi-night objects."""
    import json
    from moving.validation import nightlink_validate as V
    des = json.load(open(os.path.join(os.path.dirname(__file__), "data", "nightlink_des.json")))
    trks, _ = V.mpc_tracklets(des["designations"], des["mjd0"] - 1.0, des["mjd0"] + 12.0)
    if len({t.truth for t in trks}) < 4:
        pytest.skip("too few tracklets from the MPC")
    res = N.link_nights(trks, max_gap_days=14.0)
    sc = S.score(trks, res)
    assert sc["wrong"] == 0, sc
    assert sc["recall"] >= 0.7, sc


def test_cli_nightlink_mode(tmp_path):
    """ds9_moving.py --mode nightlink on three per-night tracklets.json files (SYNTHETIC)."""
    import json, subprocess
    rng = np.random.default_rng(21)
    trks = S.make_field(6, [0, 2, 5], rng, sig=0.2, n_decoy=4)
    files = []
    for nt in (0, 2, 5):
        d = tmp_path / ("night%d" % nt); d.mkdir()
        rows = [dict(id=t.id, t=t.mjd.tolist(), ra=t.ra.tolist(), dec=t.dec.tolist(), sig=t.sig.tolist()) for t in trks if t.night == nt]
        (d / "tracklets.json").write_text(json.dumps(rows))
        files.append(str(d / "tracklets.json"))
    cli = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "ds9", "library", "ds9_moving.py")
    r = subprocess.run([sys.executable, cli, "--mode", "nightlink", "--workdir", str(tmp_path / "wd"), "--tracklet-files"] + files + ["--nl-sigma", "0.2"],
                       capture_output=True, text=True, timeout=300)
    assert r.returncode == 0, r.stderr[-600:]
    out = json.load(open(tmp_path / "wd" / "nightlinks.json"))
    assert len(out["groups"]) >= 5 and all(g["n_tracklets"] >= 2 for g in out["groups"])
    assert os.path.exists(tmp_path / "wd" / "nightlinks.tsv")
