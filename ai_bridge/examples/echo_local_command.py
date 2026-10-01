#!/usr/bin/env python3
"""EXAMPLE local_command service (not a model).  Reads the bridge JSON from stdin and writes one result per record.

Contract (ogf-ai-bridge/1):
  stdin : {"contract": "...", "task": "star_galaxy", "service": "...", "params": {...},
           "records": [{"id": "1", "x":..., "y":..., "ra":..., "dec":..., "mags": {"AUTO": 20.1},
                        "mag_errs": {...}, "cutouts": {"main": {"path": "...png", "format": "png", ...}},
                        "catalog_row": {...}}, ...]}
  stdout: {"results": [{"id": "1", "STAR_PROB": 0.5, "CLASS_LABEL": "GALAXY"}, ...],
           "model": "<name + version of YOUR model>", "request_id": "<optional>"}
  A record that fails: {"id": "7", "error": "reason"}.   Exit status != 0 = the whole batch failed.
This toy rule (brighter than mag 16 => star) only shows the plumbing.
"""
import json
import sys

req = json.load(sys.stdin)
out = []
for r in req['records']:
    mags = [m for m in (r.get('mags') or {}).values() if m is not None]
    if not mags:
        out.append({'id': r['id'], 'error': 'no magnitude'})
        continue
    p = 0.9 if min(mags) < 16 else 0.1
    out.append({'id': r['id'], 'STAR_PROB': p, 'CLASS_LABEL': 'STAR' if p > 0.5 else 'GALAXY'})
json.dump({'results': out, 'model': 'echo-example 0.0 (toy rule, not a model)'}, sys.stdout)
