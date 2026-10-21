#!/usr/bin/env python3
"""Run a restartable, checksum-tracked core/public-wing ICL sensitivity grid."""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import os
import platform
import time
from pathlib import Path

from extended_psf_injection import (
    DEFAULT_DATA, DEFAULT_PLAN_DIR, DEFAULT_PSFS, DEFAULT_REFERENCE_DIR,
    DEFAULT_TEST_CENTER, run_batch,
)
from validate_extended_psf_injection_batch import audit_batch


SPLICE_RADII_ARCSEC = (1.8, 2.4, 3.0)
WING_AMPLITUDES = (0.50, 0.65, 0.80, 1.00)
BLEND_WIDTH_ARCSEC = 0.6
REFERENCE_NAME = "HSC_PSF_i.fits"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def grid_spec(cluster_id: str, filter_name: str) -> dict:
    if filter_name not in ("i", "i2"):
        raise ValueError("filter_name must be i or i2")
    return {
        "schema": "astrafek-core-wing-hybrid-grid-v1",
        "cluster_id": str(cluster_id),
        "filter": filter_name,
        "splice_radii_arcsec": list(SPLICE_RADII_ARCSEC),
        "blend_width_arcsec": BLEND_WIDTH_ARCSEC,
        "wing_amplitudes": list(WING_AMPLITUDES),
        "reference_name": REFERENCE_NAME,
        "n_configurations": len(SPLICE_RADII_ARCSEC) * len(WING_AMPLITUDES),
        "science_measurement": False,
        "psf_approved": False,
    }


def configuration_key(parameters: dict) -> str:
    radius = float(parameters["splice_radius_arcsec"])
    amplitude = float(parameters["wing_amplitude"])
    return f"r{radius:.1f}_a{amplitude:.2f}"


def configurations() -> list[dict]:
    return [
        {"splice_radius_arcsec": radius,
         "blend_width_arcsec": BLEND_WIDTH_ARCSEC,
         "wing_amplitude": amplitude}
        for radius, amplitude in itertools.product(SPLICE_RADII_ARCSEC,
                                                    WING_AMPLITUDES)
    ]


def _save_checkpoint(path: Path, state: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(state, indent=2, sort_keys=True,
                                     allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def run_grid(output_root: Path, *, cluster_id: str, filter_name: str = "i",
             data: Path = DEFAULT_DATA,
             reference_dir: Path = DEFAULT_REFERENCE_DIR,
             plan_dir: Path = DEFAULT_PLAN_DIR,
             psf_dir: Path = DEFAULT_PSFS,
             test_center_psf: Path = DEFAULT_TEST_CENTER) -> dict:
    if (not platform.node().lower().startswith("lageunha")
            and not os.environ.get("SLURM_JOB_ID")):
        raise RuntimeError("run the full-image hybrid grid on lageunha or through Slurm")
    output_root = Path(output_root).resolve()
    spec = grid_spec(cluster_id, filter_name)
    state_path = output_root / "grid_status.json"
    if output_root.exists():
        if not state_path.is_file():
            raise FileExistsError("existing grid root has no checkpoint; refusing to reuse it")
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("spec") != spec:
            raise ValueError("grid checkpoint does not match the requested specification")
    else:
        output_root.mkdir(parents=True, exist_ok=False)
        state = {"spec": spec, "status": "running", "records": [],
                 "started_unix": time.time(), "errors": []}
        _save_checkpoint(state_path, state)

    reference_path = reference_dir / REFERENCE_NAME
    completed = {row["configuration"]: row for row in state["records"]
                 if row.get("status") == "complete"}
    for index, parameters in enumerate(configurations(), 1):
        key = configuration_key(parameters)
        batch_dir = output_root / key
        if key in completed:
            summary_path = batch_dir / "summary.json"
            if (not summary_path.is_file()
                    or sha256(summary_path) != completed[key].get("summary_sha256")):
                raise ValueError(f"completed configuration {key} failed checkpoint hash verification")
            audit_batch(batch_dir, expected_filter=filter_name, expected_frames=1)
            continue
        if batch_dir.exists():
            raise FileExistsError(f"uncheckpointed output exists for {key}: {batch_dir}")

        print(f"[{index}/{spec['n_configurations']}] {key}", flush=True)
        started = time.time()
        summary = run_batch(
            Path(data), batch_dir, reference_path, Path(reference_dir),
            Path(plan_dir), Path(psf_dir), Path(test_center_psf),
            only_cluster=str(cluster_id), filter_name=filter_name,
            hybrid_parameters=parameters,
        )
        if summary.get("status") != "complete" or summary.get("n_errors") != 0:
            message = f"configuration {key} did not complete cleanly"
            state["status"] = "failed"
            state["errors"].append({"configuration": key, "error": message})
            _save_checkpoint(state_path, state)
            raise RuntimeError(message)
        validation = audit_batch(batch_dir, expected_filter=filter_name,
                                 expected_frames=1)
        summary_path = batch_dir / "summary.json"
        state["records"].append({
            "configuration": key,
            "parameters": parameters,
            "status": "complete",
            "summary_path": str(summary_path),
            "summary_sha256": sha256(summary_path),
            "n_frames_validated": validation["n_frames_validated"],
            "max_arithmetic_error": validation["max_stored_vs_recomputed_response_error"],
            "elapsed_seconds": time.time() - started,
        })
        _save_checkpoint(state_path, state)

    if len(state["records"]) != spec["n_configurations"]:
        raise ValueError("grid checkpoint has an unexpected configuration count")
    state["status"] = "complete"
    state["finished_unix"] = time.time()
    _save_checkpoint(state_path, state)
    return state


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--cluster-id", default="1056")
    parser.add_argument("--filter", choices=("i", "i2"), default="i")
    args = parser.parse_args()
    state = run_grid(args.output_root, cluster_id=args.cluster_id,
                     filter_name=args.filter)
    print(json.dumps({"status": state["status"],
                      "configurations": len(state["records"]),
                      "output_root": str(args.output_root.resolve())}, indent=2))


if __name__ == "__main__":
    main()
