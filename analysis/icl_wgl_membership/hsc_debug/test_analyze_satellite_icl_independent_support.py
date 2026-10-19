import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).with_name("analyze_satellite_icl_independent_support.py")
SPEC = importlib.util.spec_from_file_location("analyze_satellite_icl_independent_support", SCRIPT)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def test_independent_annulus_uses_pixel_count_weighting_per_frame():
    profile = []
    for index, lo in enumerate((200, 225, 250, 275)):
        count = index + 1
        profile.append({
            "r_min_kpc": lo, "r_max_kpc": lo + 25,
            "n_target_pixels": count, "n_control_pixels": count * 2,
            "geometric_pixels": 100,
            "target_mean_njy_arcsec2": float(index + 1),
            "control_mean_njy_arcsec2": float(2 * (index + 1)),
        })
    result = analysis.independent_annulus(profile)
    assert result["target_mean_njy_arcsec2"] == 3.0
    assert result["control_mean_njy_arcsec2"] == 6.0
    assert result["delta_njy_arcsec2"] == -3.0
    assert result["target_coverage"] == 0.025
    assert result["control_coverage"] == 0.05


def test_block_bootstrap_pairs_baseline_and_alternate_draws():
    rows = []
    for index in range(24):
        x = float(index % 7) + 0.15 * index
        rows.append({
            "icl_delta_200_300": x,
            "icl_independent_delta_200_300": 1.15 * x + (index % 3) * 0.2,
            "satellite_delta_50_75": 0.7 * x + (index % 5) * 0.3,
            "redshift": 0.15 + (index % 6) * 0.04,
            "richness": 10 + (index * 7) % 23,
        })
    labels = np.asarray([f"block-{index % 8}" for index in range(len(rows))])
    result = analysis.paired_block_bootstrap(rows, labels, replicas=100, seed=17)
    assert result["replicas_valid"] > 0
    assert result["n_blocks"] == 8
    assert np.isclose(result["partial_r_difference"],
                      result["independent_support_partial_r"] - result["baseline_partial_r"])
    assert len(result["paired_difference_block_ci95"]) == 2
