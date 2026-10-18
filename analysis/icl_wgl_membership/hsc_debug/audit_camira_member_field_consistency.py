#!/usr/bin/env python3
"""Audit whether CAMIRA member columns 7/8 behave like mass and weight.

This is a schema-consistency check on the already downloaded S20A catalog, not
an authoritative column definition and not a science measurement.
"""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

EXPECTED_SHA256 = "bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e"
EXPECTED_COLUMNS = ["ra_cl_deg", "dec_cl_deg", "N_mem", "z_cl",
                    "ra_deg", "dec_deg", "field7", "weight"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quantiles(values):
    return [float(np.quantile(values, q)) for q in (0.01, 0.50, 0.99)]


def run(data_root):
    data_root = Path(data_root)
    member_path = data_root / "catalog/camira_s20a_wide_sm_member.dat.gz"
    summary_path = data_root / "logs/member_summary.json"
    manifest_path = data_root / "sample_800_2026-10-13/sample_manifest.json"
    digest = sha256(member_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if digest != EXPECTED_SHA256 or summary.get("sha256") != digest:
        raise ValueError("CAMIRA source or acquisition-summary SHA-256 mismatch")
    if summary.get("columns") != EXPECTED_COLUMNS:
        raise ValueError("Unexpected acquisition column mapping")

    clusters = json.loads(manifest_path.read_text(encoding="utf-8"))["clusters"]
    by_center = {(round(float(c["ra_deg"]), 5), round(float(c["dec_deg"]), 5)): c
                 for c in clusters}
    if len(clusters) != 800 or len(by_center) != 800:
        raise ValueError("Expected 800 unique selected clusters")
    members = {str(c["id"]): [] for c in clusters}
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)
    field7, field8, radii = [], [], []
    with gzip.open(member_path, "rt", encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            values = [float(v) for v in line.split()]
            if len(values) != 8 or not all(math.isfinite(v) for v in values):
                raise ValueError(f"Malformed CAMIRA row {line_number}")
            center = (round(values[0], 5), round(values[1], 5))
            cluster = by_center.get(center)
            if cluster is None:
                continue
            if abs(values[2] - float(cluster["richness"])) > 0.001:
                raise ValueError(f"N_mem mismatch at row {line_number}")
            z = float(cluster.get("z_cl", cluster["z"]))
            if abs(values[3] - z) > 0.001:
                raise ValueError(f"z_cl mismatch at row {line_number}")
            if not 0.0 <= values[7] <= 1.0:
                raise ValueError(f"Column 8 outside [0, 1] at row {line_number}")
            ident = str(cluster["id"])
            members[ident].append(values)
            field7.append(values[6])
            field8.append(values[7])
            dra = (values[4] - values[0]) * math.cos(math.radians(values[5]))
            ddec = values[5] - values[1]
            angle_rad = math.radians(math.hypot(dra, ddec))
            radii.append(float(cosmo.angular_diameter_distance(z).to_value("Mpc"))
                         * angle_rad)

    ratios, ratios_cut = [], []
    for cluster in clusters:
        rows = members[str(cluster["id"])]
        if not rows:
            raise ValueError(f"No member rows for cluster {cluster['id']}")
        n_mem = float(cluster["richness"])
        ratios.append(sum(r[7] for r in rows) / n_mem)
        ratios_cut.append(sum(r[7] for r in rows if r[6] >= 10.2) / n_mem)
    if len(field7) != 36513 or sum(bool(v) for v in members.values()) != 800:
        raise ValueError("Selected catalog row/cluster totals differ from expected values")

    ratios = np.asarray(ratios)
    ratios_cut = np.asarray(ratios_cut)
    abs_diff = np.abs(ratios - 1.0)
    radii = np.asarray(radii)
    return {
        "schema": "camira-s20a-member-field-consistency-v1",
        "science_measurement": False,
        "schema_verified": False,
        "interpretation": "column 7 likely log10 stellar mass; column 8 behaves like a membership weight; exact S20A member-file dictionary not located",
        "source_sha256": digest,
        "clusters": len(clusters),
        "member_rows": len(field7),
        "field7_quantiles_1_50_99": quantiles(field7),
        "field7_fraction_ge_10_2": float(np.mean(np.asarray(field7) >= 10.2)),
        "field8_quantiles_1_50_99": quantiles(field8),
        "field8_fraction_equal_1": float(np.mean(np.asarray(field8) == 1.0)),
        "sum_field8_over_Nmem_median": float(np.median(ratios)),
        "sum_field8_over_Nmem_p10_p90": [float(np.quantile(ratios, q)) for q in (0.10, 0.90)],
        "median_absolute_relative_difference": float(np.median(abs_diff)),
        "fraction_clusters_within_20_percent_of_Nmem": float(np.mean(abs_diff <= 0.20)),
        "sum_field8_for_field7_ge_10_2_over_Nmem_median": float(np.median(ratios_cut)),
        "projected_member_radius_mpc_q50_q95_q99_max": [
            float(np.quantile(radii, q)) for q in (0.50, 0.95, 0.99, 1.0)],
        "fraction_radius_above_1_h_inverse_mpc": float(np.mean(radii > (1.0 / 0.7))),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("/gpfs/kjhan/Astrafek/data"))
    args = parser.parse_args()
    print(json.dumps(run(args.data_root), indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
