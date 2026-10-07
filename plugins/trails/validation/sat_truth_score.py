#!/usr/bin/env python3
"""Match detector segments to the ACS satellite-trail labels.

The rule below was written before ``detect_trails`` ran on these exposures.
Do not change the thresholds, and do not edit ``sat_truth_labels.json``, to
improve the score. Notes taken after the run are in ``trails_validation.md``.
They are not inputs to this script.

A segment is ``{id, chip, x1, y1, x2, y2}`` in 1-based FITS pixels.
A chip with no trail is ``{id, chip, none: true}``.
The detections file is either a list of segments or the detector payload
``{"rows": [...]}``. Rows with ``rejected: true`` are not detections.

Match, per chip, directionless:

* dθ is the smallest angle between the two segment directions, in [0, 90] deg.
* dρ is the median perpendicular distance, in pixels, of 21 equally spaced
  points on the detection segment from the label's supporting line.
* overlap is the length of the intersection of the two segments projected
  onto the label axis.
* A detection matches a label when dθ <= 3 deg, dρ <= 8 px, and
  overlap >= max(100 px, 0.30 * min(length_label, length_detection)).
* One detection matches at most one label (the largest overlap).
  One label needs only one match.

Recall = matched labels / all labels. Precision = matched detections / all
detections. A chip with no label and no detection is a true-negative chip.
The same counts are also printed for chips that have a label.

usage: python3 sat_truth_score.py sat_truth_labels.json sat_truth_detections.json
"""
import json
import sys
import numpy as np


def ang(x1, y1, x2, y2):
    return float(np.degrees(np.arctan2(y2 - y1, x2 - x1)) % 180.0)


def length(s):
    return float(np.hypot(s["x2"] - s["x1"], s["y2"] - s["y1"]))


def dtheta(a, b):
    d = abs(ang(a["x1"], a["y1"], a["x2"], a["y2"]) - ang(b["x1"], b["y1"], b["x2"], b["y2"])) % 180.0
    return min(d, 180.0 - d)


def d_rho(label, det):
    dx, dy = label["x2"] - label["x1"], label["y2"] - label["y1"]
    L = np.hypot(dx, dy)
    nx, ny = -dy / L, dx / L
    ts = np.linspace(0.0, 1.0, 21)
    xs = det["x1"] + ts * (det["x2"] - det["x1"])
    ys = det["y1"] + ts * (det["y2"] - det["y1"])
    off = (xs - label["x1"]) * nx + (ys - label["y1"]) * ny
    return float(np.median(np.abs(off)))


def overlap(label, det):
    dx, dy = label["x2"] - label["x1"], label["y2"] - label["y1"]
    L = np.hypot(dx, dy)
    ux, uy = dx / L, dy / L

    def span(s):
        a = (s["x1"] - label["x1"]) * ux + (s["y1"] - label["y1"]) * uy
        b = (s["x2"] - label["x1"]) * ux + (s["y2"] - label["y1"]) * uy
        return (min(a, b), max(a, b))

    a0, a1 = span(label)
    b0, b1 = span(det)
    return float(max(0.0, min(a1, b1) - max(a0, b0)))


def matches(label, det):
    if dtheta(label, det) > 3.0:
        return False
    if d_rho(label, det) > 8.0:
        return False
    ov = overlap(label, det)
    return ov >= max(100.0, 0.30 * min(length(label), length(det)))


def load_segments(path, drop_rejected):
    data = json.load(open(path))
    if isinstance(data, dict):
        data = data["rows"]
    if drop_rejected:
        data = [s for s in data if not s.get("rejected")]
    return data


def main():
    labels = load_segments(sys.argv[1], drop_rejected=False)
    dets = load_segments(sys.argv[2], drop_rejected=True)
    labs = [s for s in labels if not s.get("none")]
    chips = sorted({(s["id"], s["chip"]) for s in labels})
    used = set()
    tp_lab = 0
    rows = []
    for lab in labs:
        hits = []
        for i, det in enumerate(dets):
            if det["id"] != lab["id"] or det["chip"] != lab["chip"]:
                continue
            if i in used:
                continue
            if matches(lab, det):
                hits.append((overlap(lab, det), i, det))
        if hits:
            hits.sort(reverse=True)
            _, i, det = hits[0]
            used.add(i)
            tp_lab += 1
            rows.append((lab, det, True))
        else:
            rows.append((lab, None, False))
    tp_det = len(used)
    fp = [d for i, d in enumerate(dets) if i not in used]
    fn = len(labs) - tp_lab
    tn_chips = 0
    for id_, chip in chips:
        nlab = sum(1 for s in labs if s["id"] == id_ and s["chip"] == chip)
        ndet = sum(1 for s in dets if s["id"] == id_ and s["chip"] == chip)
        if nlab == 0 and ndet == 0:
            tn_chips += 1
    labelled = {(s["id"], s["chip"]) for s in labs}
    det_on = [d for d in dets if (d["id"], d["chip"]) in labelled]
    fp_on = [d for d in fp if (d["id"], d["chip"]) in labelled]
    print("labels %d  detections %d" % (len(labs), len(dets)))
    print("TP labels %d  FN %d  TP detections %d  FP %d  true-negative chips %d" % (
        tp_lab, fn, tp_det, len(fp), tn_chips))
    if labs:
        print("recall %.3f (%d/%d)" % (tp_lab / len(labs), tp_lab, len(labs)))
    else:
        print("recall undefined (0 labels)")
    if dets:
        print("precision %.3f (%d/%d)" % (tp_det / len(dets), tp_det, len(dets)))
    else:
        print("precision undefined (0 detections)")
    if det_on:
        print("on labelled chips: detections %d  precision %.3f (%d/%d)" % (
            len(det_on), tp_det / len(det_on), tp_det, len(det_on)))
    print("false detections on chips with no label: %d" % (len(fp) - len(fp_on)))
    print("--- per label ---")
    for lab, det, ok in rows:
        if ok:
            print("MATCH %s %s dθ %.2f dρ %.2f overlap %.0f det (%.0f,%.0f)-(%.0f,%.0f)" % (
                lab["id"], lab["chip"], dtheta(lab, det), d_rho(lab, det), overlap(lab, det),
                det["x1"], det["y1"], det["x2"], det["y2"]))
        else:
            print("MISS  %s %s label (%.0f,%.0f)-(%.0f,%.0f)" % (
                lab["id"], lab["chip"], lab["x1"], lab["y1"], lab["x2"], lab["y2"]))
    print("--- false detections ---")
    for d in fp:
        print("FP    %s %s (%.0f,%.0f)-(%.0f,%.0f) theta %.1f" % (
            d["id"], d["chip"], d["x1"], d["y1"], d["x2"], d["y2"], ang(d["x1"], d["y1"], d["x2"], d["y2"])))


if __name__ == "__main__":
    main()
