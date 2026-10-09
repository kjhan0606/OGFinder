#!/usr/bin/env python3
"""Score a single elliptical Sérsic against truth with the in-tree model.

    python compare_pysersic.py FIT.json WORKDIR OUT.json
    python compare_pysersic.py CMP.json FIT.json WORKDIR OUT.json

FIT.json is the list written by pysersic_map.py (ogfmeas.sersic). CMP.json is
optional. A row is scored when it is case A, or when no case is recorded.
Rows do not need a stored external-fitter solution. Chi-squared is the sum of
squared residuals of the ogfmeas model on that cutout. The external binary is
not executed. ``--galfit`` and ``--ld`` are accepted and ignored.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import pysersic_map as PM  # noqa: E402


def _positional(argv):
    pos = []
    i = 0
    while i < len(argv):
        token = argv[i]
        if token in ("--galfit", "--ld", "--n"):
            i += 2
            continue
        if token.startswith("--"):
            i += 2 if i + 1 < len(argv) and not argv[i + 1].startswith("--") else 1
            continue
        pos.append(token)
        i += 1
    return pos


def _load_list(path):
    payload = json.load(open(path))
    if isinstance(payload, dict):
        payload = payload.get("results") or payload.get("rows") or []
    return payload


def _truth_for(cmp_row, work, seed):
    if cmp_row and cmp_row.get("truth"):
        truth = cmp_row["truth"][0]
        return dict(x=truth["x"], y=truth["y"], mag=truth["mag"], re=truth["re"], n=truth["n"],
                    q=truth["q"], pa=truth["pa"], sky=cmp_row.get("truth_sky"))
    feed = os.path.join(work, "A%02d" % seed, "truth.feedme")
    if not os.path.isfile(feed):
        return None
    from ogfkit import galfitio as GI
    from ogfmeas.sersic import galfit_pa_from_theta
    import math
    cfg = GI.parse_feedme(feed, strict=False, base_dir=os.path.dirname(feed))
    comps = cfg.get("components") or []
    if len(comps) != 1 or comps[0].get("galfit_type") != "sersic":
        return None
    component = comps[0]
    return dict(x=component["x"], y=component["y"], mag=component["mag"], re=component["re"], n=component["n"],
                q=component["q"], pa=galfit_pa_from_theta(math.radians(component["pa"])),
                sky=cfg.get("sky_value_galfit"))


def _deltas(fit, truth):
    return dict(
        x=fit["x"] - truth["x"], y=fit["y"] - truth["y"], mag=fit["mag"] - truth["mag"],
        re_rel=fit["re"] / truth["re"] - 1.0, n_rel=fit["n"] / truth["n"] - 1.0,
        q=fit["q"] - truth["q"], pa=(fit["pa"] - truth["pa"] + 90.0) % 180.0 - 90.0,
    )


def _summary(rows):
    usable = [row["ogfmeas"] for row in rows if row.get("ogfmeas")]
    if not usable:
        return {}
    keys = usable[0]
    return {key: dict(
        median=float(np.median([item[key] for item in usable])),
        std=float(np.std([item[key] for item in usable])),
        p90_abs=float(np.percentile(np.abs([item[key] for item in usable]), 90)),
    ) for key in keys}


def compare(cmp_json, fit_json, work, out):
    fits = {row["seed"]: row for row in _load_list(fit_json) if "error" not in row and "seed" in row}
    cmp_rows = {row["seed"]: row for row in _load_list(cmp_json)} if cmp_json else {}
    seeds = sorted(fits) if not cmp_rows else sorted(set(fits) & set(cmp_rows))
    rows = []
    for seed in seeds:
        cmp_row = cmp_rows.get(seed)
        if cmp_row is not None and cmp_row.get("kind") not in (None, "A"):
            continue
        fit = fits[seed]
        if fit.get("fitter") not in (None, "ogfmeas.sersic"):
            continue
        wd = os.path.join(work, "A%02d" % seed)
        chi2, npix = PM.pixel_chi2(wd, fit)
        truth = _truth_for(cmp_row, work, seed)
        row = dict(seed=seed, ogfmeas_chi2=chi2, npix=npix, fitter="ogfmeas.sersic", time=fit.get("time"))
        if truth is not None:
            row["ogfmeas"] = _deltas(fit, truth)
        rows.append(row)
    summary = dict(ogfmeas=_summary(rows), n=len(rows), fitter="ogfmeas.sersic")
    json.dump(dict(rows=rows, summary=summary), open(out, "w"))
    print("n", len(rows), "fitter ogfmeas.sersic")
    if summary["ogfmeas"]:
        print("ogfmeas", " ".join("%s %+.3f/%.3f" % (key, value["median"], value["p90_abs"])
                                   for key, value in summary["ogfmeas"].items()))
    return summary


def main(argv=None):
    pos = _positional(list(sys.argv[1:] if argv is None else argv))
    if len(pos) == 3:
        cmp_json, fit_json, work, out = None, pos[0], pos[1], pos[2]
    elif len(pos) == 4:
        cmp_json, fit_json, work, out = pos
    else:
        sys.exit("usage: compare_pysersic.py [CMP.json] FIT.json WORKDIR OUT.json")
    compare(cmp_json, fit_json, work, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
