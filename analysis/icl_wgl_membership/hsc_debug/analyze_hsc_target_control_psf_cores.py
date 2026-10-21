#!/usr/bin/env python3
"""Describe target/control coadd PSF-core mismatch from approved HSC Picker data.

This is a diagnostic of the local model-PSF cores only. It does not estimate
extended PSF wings, alter science images, or approve a PSF correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prepare_hsc_psf_picker_target_control_plan import collect_positions, _load_rows
from summarize_hsc_psf_core_variation import measure_core


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PIXEL_SCALE_ARCSEC = 0.168
TARGET_PLAN = ROOT / "icl_psf_picker_target_control_centers_plan_2026-10-21_01"
TARGET_DOWNLOAD = ROOT / "icl_psf_picker_target_control_centers_download_2026-10-21_01"
INJECTION_PLAN = ROOT / "icl_psf_picker_injection_plan_2026-10-19_02"
INJECTION_DOWNLOAD = ROOT / "icl_psf_picker_injection_2026-10-19_01"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _key(row: dict) -> tuple[float, float]:
    return round(float(row["ra_deg"]), 9), round(float(row["dec_deg"]), 9)


def _read_audit(plan_dir: Path, download_dir: Path, *, count: int,
                schema: str, audit_name: str) -> tuple[dict, dict]:
    plan_path = plan_dir / "psf_picker_plan.json"
    mapping_path = plan_dir / "position_mapping.json"
    coord_path = plan_dir / "psf_picker_coordinates.txt"
    audit_path = download_dir / audit_name
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    if audit.get("status") != "verified" or audit.get("n_returned") != count:
        raise ValueError(f"incomplete PSF audit: {audit_path}")
    if plan.get("schema") != schema:
        raise ValueError(f"unexpected PSF plan schema: {plan_path}")
    if audit.get("science_measurement") is not False or audit.get("psf_approved") is True:
        raise ValueError(f"PSF audit has been promoted beyond a diagnostic: {audit_path}")
    if sha256(mapping_path) != plan.get("position_mapping_sha256"):
        raise ValueError(f"mapping checksum differs from plan: {mapping_path}")
    if sha256(coord_path) != plan.get("coordinate_list_sha256"):
        raise ValueError(f"coordinate checksum differs from plan: {coord_path}")
    if (audit.get("position_mapping_sha256") != sha256(mapping_path)
            or audit.get("coordinate_list_sha256") != sha256(coord_path)):
        raise ValueError(f"audit was generated from different plan inputs: {audit_path}")
    if len(mapping) != count:
        raise ValueError(f"unexpected mapping count in {mapping_path}")
    records = {}
    for row in audit.get("files", []):
        key = (round(float(row["ra_deg_filename"]), 9),
               round(float(row["dec_deg_filename"]), 9))
        if key in records:
            raise ValueError(f"duplicate coordinate in {audit_path}: {key}")
        path = download_dir / row["file"]
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise ValueError(f"PSF file missing or hash mismatch: {path}")
        records[key] = (row, path)
    if len(records) != count:
        raise ValueError(f"audit file record count differs from {count}: {audit_path}")
    return {"plan_path": plan_path, "mapping_path": mapping_path,
            "coordinate_path": coord_path, "audit_path": audit_path,
            "plan": plan, "audit": audit, "mapping": mapping,
            "records": records}, records


def distribution(values: list[float]) -> dict:
    data = np.asarray(values, dtype=float)
    if not len(data) or not np.isfinite(data).all():
        raise ValueError("distribution needs finite values")
    return {"n": int(len(data)), "median": float(np.median(data)),
            "q16": float(np.quantile(data, 0.16)),
            "q84": float(np.quantile(data, 0.84)),
            "q05": float(np.quantile(data, 0.05)),
            "q95": float(np.quantile(data, 0.95)),
            "minimum": float(np.min(data)), "maximum": float(np.max(data))}


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _pair_rows(targets: dict[str, dict], controls: dict[str, dict]) -> list[dict]:
    rows = []
    for target_id in sorted(controls):
        if target_id not in targets:
            raise ValueError(f"matched control has no target PSF: {target_id}")
        target, control = targets[target_id], controls[target_id]
        row = {"target_cluster_id": target_id,
               "control_id": control["control_id"],
               "target_psf_file": target["file"],
               "control_psf_file": control["file"]}
        for name in ("moment_fwhm_major_arcsec", "moment_fwhm_minor_arcsec",
                     "moment_axis_ratio_minor_over_major",
                     "encircled_energy_r50_arcsec", "encircled_energy_r80_arcsec",
                     "encircled_energy_r90_arcsec"):
            t, c = float(target[name]), float(control[name])
            row[f"target_{name}"] = t
            row[f"control_{name}"] = c
            row[f"target_minus_control_{name}"] = t - c
            row[f"absolute_difference_{name}"] = abs(t - c)
        rows.append(row)
    return rows


def analyze(target_plan_dir: Path, target_download_dir: Path,
            injection_plan_dir: Path, injection_download_dir: Path,
            output_dir: Path) -> dict:
    target_plan_dir, target_download_dir = map(Path, (target_plan_dir, target_download_dir))
    injection_plan_dir, injection_download_dir = map(Path, (injection_plan_dir, injection_download_dir))
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_dir}")

    target, _ = _read_audit(target_plan_dir, target_download_dir,
                            count=1542,
                            schema="astrafek-hsc-psf-picker-target-control-center-plan-v1",
                            audit_name="psf_picker_audit_reproducible.json")
    injection, _ = _read_audit(injection_plan_dir, injection_download_dir,
                               count=744,
                               schema="astrafek-hsc-psf-picker-injection-plan-v1",
                               audit_name="psf_picker_audit.json")
    tplan = target["plan"]
    iplan = injection["plan"]
    if (tplan.get("schema") != "astrafek-hsc-psf-picker-target-control-center-plan-v1"
            or tplan.get("n_target_centers") != 800
            or tplan.get("n_control_centers") != 798
            or tplan.get("n_unique_combined_centers") != 1598
            or tplan.get("n_already_covered_by_approved_psfs") != 56):
        raise ValueError("target/control plan does not describe the approved 800/798 sample")
    if (iplan.get("schema") != "astrafek-hsc-psf-picker-injection-plan-v1"
            or iplan.get("user_approved_requests") != 744):
        raise ValueError("prior 744-position injection plan is not the approved plan")
    for doc, label in ((target, "target/control"), (injection, "injection")):
        if (doc["plan"].get("filter") != "HSC-I"
                or doc["plan"].get("rerun") != "pdr3_wide"
                or doc["plan"].get("type") != "coadd"):
            raise ValueError(f"unexpected product in {label} plan")

    target_manifest = Path(tplan["target_manifest"])
    control_manifest = Path(tplan["control_manifest"])
    if (sha256(target_manifest) != tplan["target_manifest_sha256"]
            or sha256(control_manifest) != tplan["control_manifest_sha256"]):
        raise ValueError("target or matched-control manifest hash mismatch")
    targets_input = _load_rows(target_manifest, "target")
    controls_input = _load_rows(control_manifest, "control")
    approved_mapping = json.loads(Path(tplan["previous_approved_mapping"]).read_text(encoding="utf-8"))
    if sha256(Path(tplan["previous_approved_mapping"])) != tplan["previous_approved_mapping_sha256"]:
        raise ValueError("previous approved PSF mapping hash mismatch")
    all_positions, remaining, _ = collect_positions(
        targets_input, controls_input, approved_mapping
    )
    if len(all_positions) != 1598 or len(remaining) != 1542:
        raise ValueError("reconstructed target/control center counts differ from plan")

    target_mapping = {_key(row): row for row in target["mapping"]}
    if len(target_mapping) != 1542 or set(target_mapping) != {_key(row) for row in remaining}:
        raise ValueError("target/control plan mapping differs from reconstructed remaining centers")
    old_mapping = {_key(row): row for row in injection["mapping"]}
    if len(old_mapping) != 744:
        raise ValueError("prior approved mapping does not contain 744 unique positions")
    merged = {}
    for position in all_positions:
        key = _key(position)
        if key in target["records"]:
            source = target["records"][key]
            source_name = "target_control_center_batch"
        elif key in injection["records"]:
            source = injection["records"][key]
            source_name = "approved_injection_batch_overlap"
        else:
            raise ValueError(f"no approved PSF core for sample center {key}")
        row, path = source
        metrics = measure_core(_read_psf(path))
        merged[key] = {"file": path.name, "sha256": row["sha256"],
                       "ra_deg": key[0], "dec_deg": key[1],
                       "source": source_name, **metrics}
    if len(merged) != 1598:
        raise ValueError(f"expected 1,598 center PSFs, measured {len(merged)}")

    target_by_id = {}
    control_by_target = {}
    center_rows = []
    for position in all_positions:
        key = _key(position)
        record = merged[key]
        roles = position["roles"]
        for role in roles:
            name = role["role"]
            if name == "target_cluster_center":
                ident = str(role["cluster_id"])
                if ident in target_by_id:
                    raise ValueError(f"duplicate target center role: {ident}")
                target_by_id[ident] = record
            elif name == "matched_control_center":
                target_id = str(role["target_cluster_id"])
                if target_id in control_by_target:
                    raise ValueError(f"multiple controls mapped to target: {target_id}")
                control_by_target[target_id] = {**record, "control_id": str(role["control_id"])}
            else:
                raise ValueError(f"unexpected sample center role: {name}")
        center_rows.append({
            "ra_deg": key[0], "dec_deg": key[1], "file": record["file"],
            "sha256": record["sha256"], "psf_source": record["source"],
            "roles_json": json.dumps(roles, ensure_ascii=False, sort_keys=True),
            **{k: v for k, v in record.items() if k not in {"file", "sha256", "ra_deg", "dec_deg", "source"}},
        })
    if len(target_by_id) != 800 or len(control_by_target) != 798:
        raise ValueError(f"unexpected sample roles: targets={len(target_by_id)} controls={len(control_by_target)}")
    pairs = _pair_rows(target_by_id, control_by_target)
    if len(pairs) != 798:
        raise ValueError("expected 798 paired target/control PSF-core measurements")

    metrics = ("moment_fwhm_major_arcsec", "moment_fwhm_minor_arcsec",
               "moment_axis_ratio_minor_over_major", "encircled_energy_r50_arcsec",
               "encircled_energy_r80_arcsec", "encircled_energy_r90_arcsec")
    distributions = {}
    for label, records in (("target_800", list(target_by_id.values())),
                           ("matched_control_798", list(control_by_target.values()))):
        distributions[label] = {metric: distribution([float(r[metric]) for r in records])
                                for metric in metrics}
    pair_differences = {}
    for metric in metrics:
        values = [float(row[f"target_minus_control_{metric}"]) for row in pairs]
        pair_differences[metric] = {
            **distribution(values),
            "n_target_wider_or_larger": sum(value > 0 for value in values),
            "n_control_wider_or_larger": sum(value < 0 for value in values),
            "n_equal_at_machine_precision": sum(value == 0 for value in values),
            "median_absolute_difference": float(np.median(np.abs(values))),
        }

    output_dir.mkdir(parents=True, exist_ok=False)
    center_path = output_dir / "center_psf_core_metrics.tsv"
    pair_path = output_dir / "paired_target_control_psf_core_metrics.tsv"
    _write_tsv(center_path, center_rows)
    _write_tsv(pair_path, pairs)
    script_path = Path(__file__).resolve()
    result = {
        "schema": "astrafek-hsc-target-control-psf-core-diagnostic-v1",
        "status": "complete_descriptive_core_mismatch_diagnostic",
        "inputs": {
            "target_control_plan": str(target["plan_path"].resolve()),
            "target_control_plan_sha256": sha256(target["plan_path"]),
            "target_control_audit": str(target["audit_path"].resolve()),
            "target_control_audit_sha256": sha256(target["audit_path"]),
            "injection_plan": str(injection["plan_path"].resolve()),
            "injection_plan_sha256": sha256(injection["plan_path"]),
            "injection_audit": str(injection["audit_path"].resolve()),
            "injection_audit_sha256": sha256(injection["audit_path"]),
            "target_manifest_sha256": sha256(target_manifest),
            "control_manifest_sha256": sha256(control_manifest),
            "pixel_scale_arcsec": PIXEL_SCALE_ARCSEC,
        },
        "support": {"target_centers": len(target_by_id),
                    "matched_control_centers": len(control_by_target),
                    "paired_target_control_centers": len(pairs),
                    "unmatched_targets": len(target_by_id) - len(pairs),
                    "unique_centers": len(center_rows),
                    "from_new_target_control_batch": sum(r["psf_source"] == "target_control_center_batch" for r in center_rows),
                    "from_prior_approved_psf_overlap": sum(r["psf_source"] == "approved_injection_batch_overlap" for r in center_rows)},
        "descriptive_core_distributions": distributions,
        "paired_target_minus_control_differences": pair_differences,
        "tables": {
            "center_psf_core_metrics": {"path": str(center_path), "sha256": sha256(center_path), "rows": len(center_rows)},
            "paired_target_control_psf_core_metrics": {"path": str(pair_path), "sha256": sha256(pair_path), "rows": len(pairs)},
        },
        "science_measurement": False,
        "psf_approved": False,
        "interpretation_limits": [
            "Moment widths and encircled-energy radii describe the position-dependent model PSF cores only.",
            "The 41-43 pixel stamps cover about 6.9-7.2 arcsec and do not measure the extended wings needed for ICL calibration.",
            "Pair differences are descriptive; no hypothesis test, confidence interval, PSF homogenization, or image correction was performed.",
            "The model-core comparison does not establish image-specific sky transfer or validate the ICL measurement.",
        ],
        "script": str(script_path),
        "script_sha256": sha256(script_path),
        "python": platform.python_version(),
        "numpy": np.__version__,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    return result


def _read_psf(path: Path) -> np.ndarray:
    from astropy.io import fits
    with fits.open(path, memmap=False, checksum=True) as hdus:
        if len(hdus) != 1 or hdus[0].data is None:
            raise ValueError(f"unexpected PSF FITS structure: {path}")
        return np.asarray(hdus[0].data, dtype=np.float64)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-plan", type=Path, default=TARGET_PLAN)
    parser.add_argument("--target-download", type=Path, default=TARGET_DOWNLOAD)
    parser.add_argument("--injection-plan", type=Path, default=INJECTION_PLAN)
    parser.add_argument("--injection-download", type=Path, default=INJECTION_DOWNLOAD)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.target_plan, args.target_download,
                     args.injection_plan, args.injection_download,
                     args.output_dir)
    print(json.dumps({"status": result["status"], "support": result["support"],
                      "paired_target_minus_control_differences": result["paired_target_minus_control_differences"],
                      "science_measurement": result["science_measurement"],
                      "psf_approved": result["psf_approved"],
                      "summary_sha256": sha256(Path(args.output_dir) / "summary.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
