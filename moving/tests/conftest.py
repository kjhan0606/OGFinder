import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
import pytest


def _net():
    try:
        import requests
        requests.get("https://ssd.jpl.nasa.gov/api/horizons.api", timeout=8)
        return True
    except Exception:
        return False


needs_net = pytest.mark.skipif(not _net(), reason="JPL Horizons not reachable")


def _eph():
    try:
        from moving import orbit
        return orbit.assist_available()
    except Exception:
        return False


needs_ephem = pytest.mark.skipif(not _eph(), reason="ASSIST + DE440/sb441 ephemeris files missing (ds9_moving.py --mode setup)")
