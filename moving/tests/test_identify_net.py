"""Live public-service checks: SkyBoT position of a known asteroid agrees with Horizons within SkyBoT's stated accuracy."""
from conftest import needs_net


@needs_net
def test_skybot_vs_horizons():
    import numpy as np
    from moving import skybot, horizons, util
    mjd = 58495.5
    e = horizons.observer_table("433;", [mjd + 2400000.5], center="500")
    try:
        res = skybot.cone(e[0, 0], e[0, 1], 0.05, mjd)
    except Exception:
        import pytest; pytest.skip("SkyBoT not reachable")
    hit = [r for r in res if r["num"] == "433"]
    assert hit, "Eros not returned by SkyBoT"
    assert util.angsep_arcsec(hit[0]["ra"], hit[0]["dec"], e[0, 0], e[0, 1]) < 3.0
