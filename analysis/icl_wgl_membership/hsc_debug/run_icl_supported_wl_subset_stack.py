#!/usr/bin/env python3
"""Run filter-separated random-subtracted WL stacks on ICL-supported lenses.

Reads only existing, hash-verified HSC exports. No HSC request is made. The
products are exploratory because random-aperture matching and several
photo-z/shear systematics remain open.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from run_wl_random_shape_stack import verify_and_load_pairs as verify_and_load_random_pairs
from run_wl_shape_stack import verify_and_load_pairs as verify_and_load_lens_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    paired_tomographic_endpoint_jackknife,
    stack_s19a_tomographic_pairs,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_AUDIT = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03/summary.json"
DEFAULT_WL_STACK = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
DEFAULT_LENS_STACK = DATA / "wl_cluster_stack_bmode0_2026-10-19_01/summary.json"
DEFAULT_LENS_ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
DEFAULT_LENS_PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
DEFAULT_LENS_PLAN = DATA / "wl_shape_rows_plan_2026-10-18_04"
DEFAULT_RANDOM_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
DEFAULT_RANDOM_ROWS = DEFAULT_RANDOM_PLAN / "results"
DEFAULT_NZ = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "icl_supported_wl_subset_stack_2026-11-05_01"


def sha256_text(values) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, values))).encode("utf-8")).hexdigest()


def subset_pairs(pairs: dict[str, np.ndarray], field: str, values: set[str]) -> dict[str, np.ndarray]:
    mask = np.isin(np.asarray(pairs[field]).astype(str), list(values))
    if not np.any(mask):
        raise ValueError(f"empty pair subset for {field}")
    return {key: np.asarray(value)[mask] for key, value in pairs.items()}


def occupied_blocks(pair_arrays: dict[str, np.ndarray], width: float):
    lens_blocks = equal_area_block_ids(
        pair_arrays["lens_ra_deg"], pair_arrays["lens_dec_deg"], cell_width_deg=width,
    )
    source_blocks = equal_area_block_ids(
        pair_arrays["source_ra_deg"], pair_arrays["source_dec_deg"], cell_width_deg=width,
    )
    return lens_blocks, source_blocks


def run(args: argparse.Namespace) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing subset-stack output: {output}")
    audit = json.loads(args.audit.read_text(encoding="utf-8"))
    wl_saved = json.loads(args.wl_stack.read_text(encoding="utf-8"))
    target_saved = json.loads(args.lens_stack.read_text(encoding="utf-8"))
    if audit.get("science_measurement") is not False or audit.get("hsc_query_submitted") is not False:
        raise ValueError("sample-overlap audit has an unexpected gate/network state")
    wl_input = wl_saved["input"]
    expected_wl_sha = "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8"
    if sha(args.wl_stack) != expected_wl_sha:
        raise ValueError("the reference full random-subtracted stack changed")
    masked_value = int(wl_input["b_mode_masked_value"])
    radial_edges = [float(v) for v in wl_input["radial_edges_mpc"]]
    h0 = float(wl_input["cosmology"]["H0_km_s_Mpc"])
    omega_m = float(wl_input["cosmology"]["Omega_m"])
    block_width = float(args.block_width_deg or wl_input["block_width_deg"])
    if not np.isfinite(block_width) or block_width <= 0:
        raise ValueError("block width must be finite and positive")
    if sha(args.nz_sacc) != wl_input["n_z_sacc_sha256"]:
        raise ValueError("HSC Y3 n(z) input differs from the reference stack")
    if sha(args.lens_stack) != wl_input["target_stack_sha256"]:
        raise ValueError("target cluster stack differs from the random-subtracted reference")

    lens_pairs, lens_prov = verify_and_load_lens_pairs(
        args.lens_rows, args.lens_preflight, args.lens_plan, masked_value,
    )
    random_pairs, random_prov = verify_and_load_random_pairs(
        args.random_rows, args.random_plan, masked_value,
    )
    if (lens_prov["row_plan_sha256"] != wl_input["lens_row_plan_sha256"]
            or random_prov["random_row_plan_sha256"] != wl_input["random_row_plan_sha256"]
            or random_prov["count_plan_sha256"] != wl_input["count_plan_sha256"]):
        raise ValueError("verified shape rows do not match the reference stack provenance")
    if random_prov["matched_lens_ids"] != wl_input["matched_lens_ids"]:
        raise ValueError("random center matched-lens selection differs from the reference stack")
    common_ids = (set(random_prov["matched_lens_ids"])
                  & set(map(str, np.unique(lens_pairs["lens_id"])))
                  & {random_prov["random_center_matched_lens_ids"][center]
                     for center, count in random_prov["random_center_expected_source_counts"].items()
                     if int(count) > 0})
    if (len(common_ids) != 435
            or sha256_text(common_ids) != audit["wl_sample"]["common_lens_ids_sha256"]):
        raise ValueError("recomputed common-support IDs disagree with the audited set")

    # Reproduce the saved 800-lens target estimator before forming subsets.
    z_nz, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=h0, Om=omega_m)
    full_target = stack_s19a_tomographic_pairs(
        lens_pairs, z_nz, nz_by_bin, radial_edges, cosmology,
    )
    saved_target_profile = target_saved.get("stack", {}).get("profile", [])
    if len(saved_target_profile) != len(full_target["profile"]):
        raise ValueError("recomputed target profile has a different bin layout")
    for index, (saved, reproduced) in enumerate(zip(saved_target_profile, full_target["profile"])):
        for key in ("delta_sigma_weak_msun_mpc2", "cross_delta_sigma_weak_msun_mpc2"):
            a, b = saved.get(key), reproduced.get(key)
            if a is None or b is None or not np.isclose(float(a), float(b), rtol=1e-8, atol=1e-8):
                raise ValueError(f"target profile reproduction failed at bin {index}: {key}")

    center_to_lens = random_prov["random_center_matched_lens_ids"]
    centers_by_lens: dict[str, set[str]] = {}
    for center_id, lens_id in center_to_lens.items():
        centers_by_lens.setdefault(str(lens_id), set()).add(str(center_id))
    results = {}
    all_selected_ids = {}
    for filt in ("i", "i2"):
        if args.selection == "icl-profile-supported":
            lens_ids = set(map(str, audit["icl_profile_support_by_filter"][filt][
                "icl_supported_ids_in_common_wl_support"
            ]))
        else:
            lens_ids = set(map(str, audit["joint_icl_satellite_support_by_filter"][filt][
                "cluster_ids_in_common_support"
            ]))
        if not lens_ids or not lens_ids <= common_ids:
            raise ValueError(f"{filt} ICL-supported set is empty or outside WL common support")
        if args.selection == "icl-profile-supported":
            expected_n = audit["icl_profile_support_by_filter"][filt][
                "n_in_verified_random_subtracted_common_support"]
        else:
            expected_n = audit["joint_icl_satellite_support_by_filter"][filt][
                "n_in_verified_random_subtracted_common_support"]
        if len(lens_ids) != expected_n:
            raise ValueError(f"{filt} selected lens list length differs from the overlap audit")
        target_subset = subset_pairs(lens_pairs, "lens_id", lens_ids)
        expected_centers = {
            center for ident in lens_ids for center in centers_by_lens.get(ident, set())
        }
        if len(expected_centers) != 5 * len(lens_ids):
            raise ValueError(f"not every {filt} cluster has exactly five matched random centers")
        random_subset = subset_pairs(random_pairs, "lens_id", expected_centers)
        expected_positive_centers = {
            center for center in expected_centers
            if int(random_prov["random_center_expected_source_counts"][center]) > 0
        }
        actual_centers = set(map(str, np.unique(random_subset["lens_id"])))
        if actual_centers != expected_positive_centers:
            raise ValueError(f"{filt} random shape rows disagree with count-only support")

        target_lens_blocks, target_source_blocks = occupied_blocks(target_subset, block_width)
        random_lens_blocks, random_source_blocks = occupied_blocks(random_subset, block_width)
        jackknife = paired_tomographic_endpoint_jackknife(
            target_subset, random_subset,
            target_lens_blocks, target_source_blocks,
            random_lens_blocks, random_source_blocks,
            z_nz, nz_by_bin, radial_edges, cosmology,
            block_method=f"matched delete-one equal-area sky cells; RA/sin(Dec) width {block_width:g} deg",
        )
        n_target_clusters = len(set(map(str, np.unique(target_subset["lens_id"]))))
        if n_target_clusters != len(lens_ids):
            raise ValueError(f"{filt} selected target lenses lost source rows after subsetting")
        results[filt] = {
            "status": "exploratory",
            "science_measurement": False,
            "n_icl_profile_supported_lenses": len(lens_ids),
            "n_target_lenses_with_shape_rows": n_target_clusters,
            "n_matched_random_centers": len(expected_centers),
            "n_random_centers_with_shape_rows": len(actual_centers),
            "n_target_shape_pairs": int(len(target_subset["lens_id"])),
            "n_random_shape_pairs": int(len(random_subset["lens_id"])),
            "lens_ids_sha256": sha256_text(lens_ids),
            "lens_ids": sorted(lens_ids, key=int),
            "paired_endpoint_jackknife": jackknife,
        }
        all_selected_ids[filt] = sorted(lens_ids, key=int)

    payload = {
        "schema": "astrafek-icl-supported-wl-subset-stack-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "method": "existing HSC S19A calibrated tomographic estimator; matched-random subtraction; paired endpoint spatial jackknife",
        "selection": ("ICL 200-300 kpc profile-supported lens IDs intersected with the exact 435-lens common cluster/random shape support"
                      if args.selection == "icl-profile-supported" else
                      "joint ICL and 50-75 kpc satellite-profile-supported lens IDs intersected with the exact 435-lens common cluster/random shape support"),
        "selection_key": args.selection,
        "n_z_sacc_path": str(args.nz_sacc.resolve()),
        "n_z_sacc_sha256": sha(args.nz_sacc),
        "target_stack_sha256": sha(args.lens_stack),
        "reference_random_subtracted_stack_sha256": sha(args.wl_stack),
        "audit_summary_sha256": sha(args.audit),
        "lens_row_plan_sha256": lens_prov["row_plan_sha256"],
        "lens_row_batch_sha256": hashlib.sha256("\n".join(
            f"{key}:{value}" for key, value in sorted(lens_prov["row_batch_sha256"].items())
        ).encode()).hexdigest(),
        "random_row_plan_sha256": random_prov["random_row_plan_sha256"],
        "random_row_batch_sha256": hashlib.sha256("\n".join(
            f"{key}:{value}" for key, value in sorted(random_prov["random_row_batch_sha256"].items())
        ).encode()).hexdigest(),
        "count_plan_sha256": random_prov["count_plan_sha256"],
        "cosmology": {"H0_km_s_Mpc": h0, "Omega_m": omega_m},
        "radial_edges_mpc": radial_edges,
        "b_mode_masked_value": masked_value,
        "spatial_block_width_deg": block_width,
        "common_support_lenses": len(common_ids),
        "common_support_lens_ids_sha256": sha256_text(common_ids),
        "by_filter": results,
        "limitations": [
            "The 5-degree equal-area cells are not validated contiguous HSC jackknife regions.",
            "Random-aperture depth, mask, PSF, and sky matching are incomplete.",
            "Photo-z calibration uncertainty, reduced shear, and source-density boost correction are not propagated.",
            "This ICL-supported stack does not itself test the satellite-profile correlation.",
            "The point estimates remain exploratory; science_measurement is false.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "science_measurement": False,
        "by_filter": {
            filt: {
                "n_lenses": result["n_icl_profile_supported_lenses"],
                "n_target_pairs": result["n_target_shape_pairs"],
                "n_random_pairs": result["n_random_shape_pairs"],
                "n_blocks": result["paired_endpoint_jackknife"]["n_occupied_blocks"],
            } for filt, result in results.items()
        },
        "output": str(path),
    }, ensure_ascii=False))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--wl-stack", type=Path, default=DEFAULT_WL_STACK)
    parser.add_argument("--lens-stack", type=Path, default=DEFAULT_LENS_STACK)
    parser.add_argument("--lens-rows", type=Path, default=DEFAULT_LENS_ROWS)
    parser.add_argument("--lens-preflight", type=Path, default=DEFAULT_LENS_PREFLIGHT)
    parser.add_argument("--lens-plan", type=Path, default=DEFAULT_LENS_PLAN)
    parser.add_argument("--random-rows", type=Path, default=DEFAULT_RANDOM_ROWS)
    parser.add_argument("--random-plan", type=Path, default=DEFAULT_RANDOM_PLAN)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-width-deg", type=float, default=None,
                        help="spatial jackknife width; defaults to the reference 5-degree value")
    parser.add_argument("--selection", choices=("icl-profile-supported", "joint-satellite-supported"),
                        default="icl-profile-supported")
    return parser


if __name__ == "__main__":
    run(build_parser().parse_args())
