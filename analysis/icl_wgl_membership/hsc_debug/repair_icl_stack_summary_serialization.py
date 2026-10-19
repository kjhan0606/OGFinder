#!/usr/bin/env python3
"""Correct a bool-to-int JSON serialization defect without changing stack data."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--executed-script", type=Path, required=True)
    args = parser.parse_args()
    summary_path = args.summary.resolve()
    script_path = args.executed_script.resolve()
    raw = summary_path.read_bytes()
    before_sha = hashlib.sha256(raw).hexdigest()
    summary = json.loads(raw)
    executed_hash = sha256(script_path)
    named_script = str((Path(__file__).resolve().parent / "stack_hsc_icl_target_control_2d.py").resolve())
    expected = summary.get("script_sha256", {}).get(named_script)
    if expected != executed_hash:
        parser.error("saved executed script does not match the run's recorded code hash")
    if type(summary.get("science_measurement")) is not int or summary["science_measurement"] != 0:
        parser.error("summary does not contain the specific false-to-zero serialization defect")
    if "serialization_correction" in summary:
        parser.error("summary already has a serialization correction record")
    summary["science_measurement"] = False
    summary["executed_script_snapshot"] = str(script_path)
    summary["executed_script_snapshot_sha256"] = executed_hash
    summary["serialization_correction"] = {
        "field": "science_measurement", "incorrect_json_value": 0,
        "correct_json_value": False, "data_products_modified": False,
        "pre_correction_summary_sha256": before_sha,
        "executed_script_snapshot": str(script_path),
        "executed_script_snapshot_sha256": executed_hash,
        "repair_script": str(Path(__file__).resolve()),
        "repair_script_sha256": sha256(Path(__file__).resolve()),
    }
    temporary = summary_path.with_name(summary_path.name + ".repair.tmp")
    if temporary.exists():
        parser.error(f"refusing to overwrite temporary file: {temporary}")
    temporary.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(summary_path)
    print(json.dumps({"status": "corrected_metadata_only", "science_measurement": False,
                      "pre_correction_summary_sha256": before_sha,
                      "post_correction_summary_sha256": sha256(summary_path),
                      "executed_script_sha256": executed_hash,
                      "repair_script_sha256": sha256(Path(__file__).resolve())}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
