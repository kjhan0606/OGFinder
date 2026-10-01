"""Coordinate / time / MPC-format utilities (no network)."""
import numpy as np
from moving import util, obs as OBS, kepler as K


def test_tangent_offsets_roundtrip():
    from moving.align import tangent_inverse
    xi, eta = util.tangent_offsets(np.array([150.01]), np.array([2.34]), 150.0, 2.3)
    ra, de = tangent_inverse(xi, eta, 150.0, 2.3)
    assert abs(ra[0] - 150.01) < 1e-9 and abs(de[0] - 2.34) < 1e-9
    assert abs(xi[0] - 0.01 * np.cos(np.radians(2.3)) * 3600) < 0.05   # small-angle check, arcsec


def test_mpc80_roundtrip():
    rec = [dict(mjd_utc=58495.123456, ra=150.123456, dec=-2.654321, stn="F51", cat="L", band="r", mag=21.3)]
    txt = OBS.write_obs80(rec, desig="25153", stn="F51")
    back = OBS.parse_obs80(txt)[0]
    assert abs(back["mjd_utc"] - rec[0]["mjd_utc"]) < 1.2e-7 + 0.5e-5      # 1e-5 day format precision
    assert abs(back["ra"] - rec[0]["ra"]) * 3600 * np.cos(np.radians(rec[0]["dec"])) < 0.16   # 0.01 s of RA = 0.15 arcsec (equator)
    assert abs(back["dec"] - rec[0]["dec"]) * 3600 < 0.06
    assert back["stn"] == "F51"


def test_kepler_element_roundtrip():
    r, v = K.elements_to_state(2.4, 0.1, 9.3, 195.8, 160.5, 92.2)
    el = K.state_to_elements(r, v, ecliptic_input=True)
    for k, x in zip(("a", "e", "i", "om", "w"), (2.4, 0.1, 9.3, 195.8, 160.5)):
        assert abs(el[k] - x) < 1e-9 * max(1, abs(x)) + 1e-8
    assert abs(((el["ma"] - 92.2 + 180) % 360) - 180) < 1e-7


def test_two_body_period():
    r, v = K.elements_to_state(1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    P = 2 * np.pi / np.sqrt(K.MU)
    r2, v2 = K.propagate_2body(r, v, P)
    assert np.linalg.norm(r2 - r) < 1e-9
