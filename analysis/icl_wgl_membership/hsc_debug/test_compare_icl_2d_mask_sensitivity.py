import importlib.util
from pathlib import Path


SCRIPT = Path(__file__).with_name("compare_icl_2d_mask_sensitivity.py")
SPEC = importlib.util.spec_from_file_location("compare_icl_2d_mask_sensitivity", SCRIPT)
compare = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(compare)


def test_comparison_uses_matching_filter_annulus_and_absolute_delta_change():
    bands = {}
    for factor, delta, n in ((2.5, 5.0, 300), (3.5, 4.0, 250), (5.0, 1.0, 100)):
        bands[factor] = [
            {
                "filter": band, "r_min_kpc": str(r_min), "r_max_kpc": str(r_min + 100),
                "n_pairs": str(n), "median_pair_common_coverage": "0.1",
                "mean_delta_njy_arcsec2": str(
                    delta if (band, r_min) == ("i", 200) else delta + (0.5 if band == "i2" else 0)
                ), "pair_sem_delta_njy_arcsec2": "1.0",
            }
            for band in ("i", "i2")
            for r_min in (0, 50, 100, 200, 300, 400)
        ]
    rows = compare.assemble_comparison(bands)
    assert len(rows) == 12
    target = next(row for row in rows if (row["filter"], row["r_min_kpc"]) == ("i", 200))
    assert target["baseline_k2p5_n_pairs"] == 300
    assert target["k3p5_delta_minus_k2p5_njy_arcsec2"] == -1.0
    assert target["k5p0_delta_minus_k2p5_njy_arcsec2"] == -4.0
