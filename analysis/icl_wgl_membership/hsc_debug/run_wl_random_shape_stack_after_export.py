#!/usr/bin/env python3
"""Wait for the approved random-shape export, then run its CPU stack once.

This watcher reads the existing checkpoint only. It never contacts HSC CAS or
submits queries. The stacker performs the full immutable-plan, batch-hash,
count, target-stack, mask, n(z), and radial-bin checks before writing output.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
ROWS = PLAN / "results"
CHECKPOINT = ROWS / "row_status.json"
EXPECTED_PLAN_SHA256 = "fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a"
TARGET_STACK = DATA / "wl_cluster_stack_bmode0_2026-10-19_01/summary.json"
LENS_ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
LENS_PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
LENS_PLAN = DATA / "wl_shape_rows_plan_2026-10-18_04"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
OUTPUT = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
RUNNER = Path(__file__).with_name("run_wl_random_shape_stack.py")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--submitter-pid", type=int, required=True,
                        help="live approved export process; no restart is attempted")
    parser.add_argument("--poll-seconds", type=float, default=30.0)
    parser.add_argument("--log", type=Path, default=ROWS / "random_stack_after_export.log")
    return parser.parse_args()


def checkpoint_state() -> dict:
    state = json.loads(CHECKPOINT.read_text(encoding="utf-8"))
    if state.get("schema") != "hsc-s19a-random-shape-row-run-v1":
        raise ValueError("unexpected random-row checkpoint schema")
    if state.get("plan_sha256") != EXPECTED_PLAN_SHA256:
        raise ValueError("checkpoint does not match the explicitly approved row plan")
    if state.get("count_plan_sha256") != "18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625":
        raise ValueError("checkpoint count-plan hash differs from the approved plan")
    if state.get("credentials_saved") is not False or state.get("photoz_pdf_rows_included") is not False or state.get("images_returned") is not False or state.get("science_measurement") is not False:
        raise ValueError("checkpoint scope flags changed unexpectedly")
    return state


def run() -> int:
    args = parse_args()
    if args.poll_seconds <= 0:
        raise SystemExit("--poll-seconds must be positive")
    if not CHECKPOINT.is_file():
        raise FileNotFoundError(CHECKPOINT)
    if not RUNNER.is_file():
        raise FileNotFoundError(RUNNER)
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite stack output: {OUTPUT}")

    while True:
        state = checkpoint_state()
        status = state.get("status")
        if status == "complete":
            break
        if status != "running":
            raise RuntimeError(f"random-row export terminated with status={status!r}")
        try:
            os.kill(args.submitter_pid, 0)
        except ProcessLookupError as exc:
            raise RuntimeError(
                "approved export is incomplete and its submitter process is gone; "
                "no query or process was restarted"
            ) from exc
        print(f"export_status={status}; waiting {args.poll_seconds:g}s", flush=True)
        time.sleep(args.poll_seconds)

    # Re-read after observing complete to avoid racing a checkpoint replacement.
    state = checkpoint_state()
    export_summary = state.get("summary", {})
    if (state.get("status") != "complete"
            or state.get("source_shape_rows_downloaded") is not True
            or export_summary.get("status") != "complete"
            or export_summary.get("source_shape_rows_downloaded") is not True):
        raise ValueError("export completion flags are inconsistent")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite stack output: {OUTPUT}")
    args.log.parent.mkdir(parents=True, exist_ok=True)
    command = [
        "nice", "-n", "10", sys.executable, str(RUNNER),
        "--rows", str(ROWS), "--plan-dir", str(PLAN), "--nz-sacc", str(NZ_SACC),
        "--b-mode-masked-value", "0",
        "--radial-edges-mpc", "0.1", "0.2", "0.4", "0.7", "1.0", "1.5", "2.0", "3.0",
        "--output", str(OUTPUT),
        "--lens-stack", str(TARGET_STACK), "--lens-rows", str(LENS_ROWS),
        "--lens-preflight", str(LENS_PREFLIGHT), "--lens-plan", str(LENS_PLAN),
        "--block-width-deg", "5.0",
    ]
    env = os.environ.copy()
    env.update({"OPENBLAS_NUM_THREADS": "1", "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1"})
    print(f"starting_cpu_stack={OUTPUT}", flush=True)
    with args.log.open("x", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=RUNNER.parents[3], env=env,
                                stdout=log, stderr=subprocess.STDOUT, check=False)
    print(f"stack_exit_code={result.returncode}; log={args.log}", flush=True)
    if result.returncode != 0:
        return result.returncode
    if not OUTPUT.is_file():
        raise RuntimeError("stack process exited successfully without its summary")
    summary = json.loads(OUTPUT.read_text(encoding="utf-8"))
    print(json.dumps({
        "status": summary.get("status"),
        "science_measurement": summary.get("science_measurement"),
        "source_pairs": summary.get("input", {}).get("n_source_pairs"),
        "random_centers": summary.get("input", {}).get("n_random_centers"),
        "random_subtraction_applied": summary.get("random_lens_subtraction_applied"),
        "spatial_block_covariance_applied": summary.get("spatial_block_covariance_applied"),
        "output": str(OUTPUT),
    }, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
