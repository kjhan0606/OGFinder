#!/usr/bin/env python3
"""Characterize Gaia-supported HSC PSF candidates before image review."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quantiles(values: list[float]) -> dict:
    if not values:
        return {"q16": None, "median": None, "q84": None}
    q16, median, q84 = np.percentile(np.asarray(values, dtype=float), [16, 50, 84])
    return {"q16": float(q16), "median": float(median), "q84": float(q84)}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build(plan_dir: Path, gaia_work: Path, prepared: Path,
          output: Path, approved_plan_sha256: str) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite support assessment: {output}")
    plan_path = plan_dir / "plan.json"
    request_path = plan_dir / "candidates.json"
    prep_summary_path = prepared / "summary.json"
    gaia_status_path = gaia_work / "status.json"
    response_path = gaia_work / "gaia_matches.csv"
    evidence_path = gaia_work / "evidence_summary.json"
    if sha256(plan_path) != approved_plan_sha256:
        raise ValueError("plan SHA-256 differs from approved request")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    requests = json.loads(request_path.read_text(encoding="utf-8"))["candidates"]
    prep_summary = json.loads(prep_summary_path.read_text(encoding="utf-8"))
    status = json.loads(gaia_status_path.read_text(encoding="utf-8"))
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if (status.get("status") != "complete"
            or status.get("plan_sha256") != approved_plan_sha256
            or status.get("request_sha256") != sha256(request_path)
            or status.get("response_sha256") != sha256(response_path)
            or plan.get("preparation_summary_sha256") != sha256(prep_summary_path)
            or prep_summary.get("status") != "complete"
            or prep_summary.get("n_prepared") != 800
            or evidence.get("n_candidates") != len(requests)
            or evidence.get("science_measurement") is not False
            or evidence.get("psf_approved") is not False):
        raise ValueError("Gaia evidence does not match its approved candidate set")
    matches: dict[int, list[dict[str, str]]] = {}
    with response_path.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            matches.setdefault(int(row["candidate_id"]), []).append(row)
    selected = [row for row in evidence["records"]
                if row["status"] == "stellar_evidence_requires_epoch_and_psf_review"]
    selected_ids = {int(row["candidate_id"]) for row in selected}
    if len(selected_ids) != len(selected) or any(len(matches.get(i, [])) != 1 for i in selected_ids):
        raise ValueError("stellar-evidence candidates are not unique one-to-one Gaia matches")

    input_hashes = {}
    preparation_records = {str(row["id"]): row for row in prep_summary["records"]}
    results = []
    for record in selected:
        candidate_id = int(record["candidate_id"])
        candidate = requests[candidate_id]
        cluster_id = str(candidate["cluster_id"])
        table_path = prepared / f"cluster_{cluster_id}" / "psf_candidates.tsv"
        prep_path = prepared / f"cluster_{cluster_id}" / "preparation.json"
        if (plan.get("candidate_tables_sha256", {}).get(cluster_id) != sha256(table_path)):
            raise ValueError(f"candidate table changed for HSC frame {cluster_id}")
        input_hashes[str(table_path)] = sha256(table_path)
        input_hashes[str(prep_path)] = sha256(prep_path)
        table_rows = read_tsv(table_path)
        table_row = next((row for row in table_rows
                          if int(row["NUMBER"]) == int(candidate["number"])), None)
        if table_row is None:
            raise ValueError(f"candidate {candidate_id} is absent from {table_path}")
        for field, candidate_key in (("SNR_NOMINAL", "snr_nominal"),
                                     ("FWHM_MOMENT_PIX", "fwhm_moment_pix"),
                                     ("WING_VALID_FRACTION", "wing_valid_fraction")):
            if not np.isclose(float(table_row[field]), float(candidate[candidate_key]), rtol=0, atol=1e-9):
                raise ValueError(f"{field} changed for candidate {candidate_id}")
        prep = json.loads(prep_path.read_text(encoding="utf-8"))
        aggregate = preparation_records.get(cluster_id)
        if (aggregate is None
                or prep.get("psf_fwhm_locus_pix") != aggregate.get("psf_fwhm_locus_pix")
                or prep.get("n_psf_candidates") != aggregate.get("n_psf_candidates")):
            raise ValueError(f"frame preparation record mismatch for HSC frame {cluster_id}")
        locus = prep.get("psf_fwhm_locus_pix")
        fwhm = float(candidate["fwhm_moment_pix"])
        gaia = matches[candidate_id][0]
        ratio = fwhm / float(locus) if locus and float(locus) > 0 else None
        pmra_snr = (abs(float(gaia["pmra"])) / float(gaia["pmra_error"])
                    if gaia.get("pmra") not in (None, "") and gaia.get("pmra_error") not in (None, "")
                    and float(gaia["pmra_error"]) > 0 else None)
        pmdec_snr = (abs(float(gaia["pmdec"])) / float(gaia["pmdec_error"])
                     if gaia.get("pmdec") not in (None, "") and gaia.get("pmdec_error") not in (None, "")
                     and float(gaia["pmdec_error"]) > 0 else None)
        results.append({
            "candidate_id": candidate_id, "cluster_id": cluster_id,
            "detection_number": int(candidate["number"]),
            "ra_deg": float(candidate["ra_deg"]), "dec_deg": float(candidate["dec_deg"]),
            "gaia_source_id": gaia["source_id"],
            "gaia_separation_arcsec": float(gaia["separation_arcsec"]),
            "gaia_g_mag": float(gaia["phot_g_mean_mag"]),
            "gaia_phot_variable_flag": gaia["phot_variable_flag"],
            "gaia_ruwe": float(gaia["ruwe"]),
            "gaia_parallax_over_error": float(gaia["parallax_over_error"]),
            "gaia_pmra_snr": pmra_snr, "gaia_pmdec_snr": pmdec_snr,
            "gaia_star_probability": float(gaia["classprob_dsc_combmod_star"]),
            "hsc_snr_nominal": float(candidate["snr_nominal"]),
            "hsc_fwhm_moment_pix": fwhm,
            "frame_fwhm_locus_pix": float(locus) if locus is not None else None,
            "hsc_fwhm_to_locus_ratio": ratio,
            "hsc_wing_valid_fraction": float(candidate["wing_valid_fraction"]),
            "review_order": "brightest Gaia G first; then higher HSC SNR and wing coverage",
            "classification": "review_priority_only_not_psf_approved",
        })
    results.sort(key=lambda row: (row["gaia_g_mag"], -row["hsc_snr_nominal"],
                                  -row["hsc_wing_valid_fraction"], row["candidate_id"]))
    output.mkdir(parents=True, exist_ok=False)
    tsv_path = output / "review_candidates.tsv"
    with tsv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(results)
    g_mags = [row["gaia_g_mag"] for row in results]
    snrs = [row["hsc_snr_nominal"] for row in results]
    wing = [row["hsc_wing_valid_fraction"] for row in results]
    ratios = [row["hsc_fwhm_to_locus_ratio"] for row in results
              if row["hsc_fwhm_to_locus_ratio"] is not None]
    from collections import Counter
    summary = {
        "schema": "astrafek-gaia-hsc-psf-candidate-support-v1",
        "status": "complete",
        "approved_plan_sha256": approved_plan_sha256,
        "preparation_summary_sha256": sha256(prep_summary_path),
        "gaia_response_sha256": sha256(response_path),
        "gaia_evidence_summary_sha256": sha256(evidence_path),
        "prepared_candidate_tables": len({row["cluster_id"] for row in results}),
        "n_candidates_with_stellar_evidence": len(results),
        "n_cluster_frames": len({row["cluster_id"] for row in results}),
        "gaia_phot_variable_flags": dict(Counter(row["gaia_phot_variable_flag"] for row in results)),
        "gaia_g_magnitude_quantiles": quantiles(g_mags),
        "gaia_g_magnitude_at_or_brighter_than": {
            str(limit): sum(value <= limit for value in g_mags) for limit in (18, 19, 20, 21, 22)
        },
        "hsc_nominal_snr_quantiles": quantiles(snrs),
        "hsc_wing_valid_fraction_quantiles": quantiles(wing),
        "hsc_wing_valid_fraction_ge_0p8": sum(value >= 0.8 for value in wing),
        "hsc_fwhm_to_frame_locus_ratio_quantiles": quantiles(ratios),
        "images_read": False,
        "science_measurement": False,
        "psf_approved": False,
        "interpretation": [
            "The ordering prioritizes review workload only; it is not a calibrated star score.",
            "Gaia NOT_AVAILABLE variability flags do not establish non-variability.",
            "FWHM-to-locus ratios are descriptive morphology comparisons, not PSF fits.",
            "No image stamps, extended wings, HSC coadd epochs, or sky profiles were inspected.",
        ],
        "review_candidates_sha256": sha256(tsv_path),
        "input_candidate_tables_sha256": input_hashes,
    }
    summary_path = output / "summary.json"
    summary["script_sha256"] = sha256(Path(__file__))
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--gaia-work", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--approved-plan-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run on lageunha CPU")
    result = build(args.plan_dir, args.gaia_work, args.prepared,
                   args.output, args.approved_plan_sha256)
    print(json.dumps({key: result[key] for key in (
        "status", "n_candidates_with_stellar_evidence", "n_cluster_frames",
        "gaia_phot_variable_flags", "gaia_g_magnitude_quantiles",
        "hsc_wing_valid_fraction_ge_0p8", "hsc_fwhm_to_frame_locus_ratio_quantiles",
        "science_measurement", "psf_approved")}, indent=2))


if __name__ == "__main__":
    main()
