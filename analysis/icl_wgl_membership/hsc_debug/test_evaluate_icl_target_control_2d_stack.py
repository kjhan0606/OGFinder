import importlib.util
import json
from pathlib import Path


SCRIPT = Path(__file__).with_name("evaluate_icl_target_control_2d_stack.py")
SPEC = importlib.util.spec_from_file_location("evaluate_icl_target_control_2d_stack", SCRIPT)
evaluate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(evaluate)


def test_broad_annuli_combine_pixels_within_pair_then_weight_pairs_equally():
    bins = []
    for i in range(20):
        lo, hi = 25 * i, 25 * (i + 1)
        value = 2.0 if i < 2 else 100.0
        count = 10 if i == 0 else 30
        bins.append({"r_min_kpc": lo, "r_max_kpc": hi, "n_valid_pixels": count,
                     "geometric_pixels": 100, "target_mean_njy_arcsec2": value + 1,
                     "control_mean_njy_arcsec2": 1.0, "delta_mean_njy_arcsec2": value})
    input_rows = [{"target_id": "a", "control_id": "ca", "filter": "i",
                   "profile_json": json.dumps(bins)},
                  {"target_id": "b", "control_id": "cb", "filter": "i",
                   "profile_json": json.dumps(bins)}]
    pair_rows = evaluate.broad_pair_profiles(input_rows)
    first = next(row for row in pair_rows if row["target_id"] == "a" and row["r_min_kpc"] == 0)
    assert first["delta_njy_arcsec2"] == 2.0
    assert first["common_coverage_fraction"] == 0.2
    combined = evaluate.aggregate_broad_profiles(pair_rows)
    first_combined = next(row for row in combined if row["filter"] == "i" and row["r_min_kpc"] == 0)
    assert first_combined["n_pairs"] == 2
    assert first_combined["mean_delta_njy_arcsec2"] == 2.0
