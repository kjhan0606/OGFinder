"""MPC public observation API (no credentials): https://data.minorplanetcenter.net/api/get-obs"""
import json
from .util import http_get, CACHE_DIR
import os, hashlib

URL = "https://data.minorplanetcenter.net/api/get-obs"


def get_obs(designation):
    """Return dict with OBS80 text for a number or designation (cached on disk)."""
    import requests
    key = os.path.join(CACHE_DIR, "mpc_obs_%s.json" % hashlib.sha1(str(designation).encode()).hexdigest()[:16])
    if os.path.exists(key):
        return json.load(open(key))
    r = requests.get(URL, json={"desigs": [str(designation)], "output_format": ["OBS80"]}, timeout=120)
    r.raise_for_status()
    j = r.json()[0]
    os.makedirs(CACHE_DIR, exist_ok=True)
    json.dump(j, open(key, "w"))
    return j
