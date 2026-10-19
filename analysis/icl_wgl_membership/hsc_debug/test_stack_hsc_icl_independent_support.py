import importlib.util
import json
from pathlib import Path

import numpy as np


STACK_PATH = Path(__file__).with_name("stack_hsc_icl_target_control_2d.py")
STACK_SPEC = importlib.util.spec_from_file_location("icl_2d_stack", STACK_PATH)
stack = importlib.util.module_from_spec(STACK_SPEC)
STACK_SPEC.loader.exec_module(stack)

SCRIPT = Path(__file__).with_name("stack_hsc_icl_independent_support.py")
SPEC = importlib.util.spec_from_file_location("stack_hsc_icl_independent_support", SCRIPT)
independent = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(independent)


def test_independent_support_means_use_each_frames_own_pixels():
    target = np.asarray([[2.0, 100.0], [2.0, 100.0]])
    control = np.asarray([[100.0, 4.0], [100.0, 4.0]])
    target_weight = np.asarray([[1, 0], [1, 0]])
    control_weight = np.asarray([[0, 1], [0, 1]])
    rows = stack.independent_support_radial_profile(
        target, control, target_weight, control_weight,
        axis_kpc=np.asarray([-1.0, 1.0]), edges_kpc=np.asarray([0.0, 3.0]))
    row = rows[0]
    assert row["n_target_pixels"] == row["n_control_pixels"] == 2
    assert row["n_common_pixels"] == 0
    assert row["target_coverage_fraction"] == row["control_coverage_fraction"] == 0.5
    assert row["delta_mean_njy_arcsec2"] == -2.0


def test_summary_separates_support_gain_from_same_pair_estimator_change():
    independent_profile = []
    for lo, hi in zip(stack.PROFILE_EDGES_KPC[:-1], stack.PROFILE_EDGES_KPC[1:]):
        independent_profile.append({
            "r_min_kpc": float(lo), "r_max_kpc": float(hi),
            "n_target_pixels": 8, "n_control_pixels": 7, "n_common_pixels": 2,
            "target_coverage_fraction": 0.4, "control_coverage_fraction": 0.35,
            "common_coverage_fraction": 0.1, "target_mean_njy_arcsec2": 5.0,
            "control_mean_njy_arcsec2": 3.0, "delta_mean_njy_arcsec2": 2.0,
        })
    pair = {"id": "a", "control_id": "ca", "filter": "i",
            "independent_profile": independent_profile}
    lo, hi = (int(stack.PROFILE_EDGES_KPC[0]), int(stack.PROFILE_EDGES_KPC[1]))
    common_profile = [{"delta_mean_njy_arcsec2": None} for _ in range(20)]
    common_profile[0]["delta_mean_njy_arcsec2"] = 1.0
    common = independent.profile_by_pair([{
        "target_id": "a", "control_id": "ca", "filter": "i",
        "profile_json": json.dumps(common_profile),
    }])
    rows = independent.summarize_pairs([pair], common)
    first = rows[0]
    assert first["n_pairs_independent_delta"] == 1
    assert first["n_pairs_common_support"] == 1
    assert first["n_pairs_added_vs_common"] == 0
    assert first["mean_independent_delta_on_common_pairs_njy_arcsec2"] == 2.0
    assert first["mean_common_delta_on_same_pairs_njy_arcsec2"] == 1.0

    rows = independent.summarize_pairs([pair], {})
    assert rows[0]["n_pairs_independent_delta"] == 1
    assert rows[0]["n_pairs_common_support"] == 0
    assert rows[0]["n_pairs_added_vs_common"] == 1
