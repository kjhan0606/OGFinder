import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("audit_icl_sample_selection.py")
SPEC = importlib.util.spec_from_file_location("audit_icl_sample_selection", MODULE_PATH)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_validate_inputs_keeps_unique_filter_matched_subsets():
    manifest_rows = []
    profile_rows = []
    for i in range(800):
        ident = str(i)
        band = "i" if i < 400 else "i2"
        manifest_rows.append({"id": ident})
        profile_rows.append({
            "id": ident, "filter": band, "z": str(.2 + i / 10000),
            "richness": str(10 + i / 100), "ra_deg": str(i % 360),
            "dec_deg": str(-50 + (i % 100)),
            "auto_unique_center_candidate": str(i % 2 == 0),
            "center_masked_by_kron": '{"2.5": false}',
        })
    strict = []
    for i in range(137):
        profile = profile_rows[i]
        strict.append({"id": str(i), "filter": profile["filter"], "z": profile["z"],
            "richness": profile["richness"], "auto_unique_center_candidate": "True",
            "center_masked_k2p5": "False", **{field: "1.0" for field in audit.STRICT_FIELDS}})
    annulus = []
    for i in range(224):
        profile = profile_rows[200 + i]
        annulus.append({"id": str(200 + i), "filter": profile["filter"],
            "z": profile["z"], "richness": profile["richness"],
            "auto_unique_center_candidate": "True", "center_masked_by_kron": "False",
            "kron_factor": "2.5", "coverage_floor": "0.4",
            "radial_bin_kpc": "200-300", "common_coverage": "0.41"})

    parent, got_strict, got_outer = audit.validate_inputs(
        {"clusters": manifest_rows}, profile_rows, strict, annulus)
    assert len(parent) == 800
    assert len(got_strict) == 137
    assert len(got_outer) == 224
    assert len({row["id"] for row in got_outer}) == 224
    assert {row["filter"] for row in got_outer} == {"i", "i2"}


def test_describe_reports_parent_standardized_mean_difference():
    parent = [
        {"z": 0., "richness": 10., "ra_deg": 1., "dec_deg": 0.,
         "auto_unique_center_candidate": "True", "center_masked_k2p5": "False"},
        {"z": 1., "richness": 20., "ra_deg": 7., "dec_deg": 0.,
         "auto_unique_center_candidate": "False", "center_masked_k2p5": "True"},
    ]
    selected = [parent[1]]
    result = audit.describe(selected, parent)
    assert result["n_clusters"] == 1
    assert result["z"]["standardized_mean_difference"] > 0
    assert result["auto_unique_center_fraction"] == 0
    assert result["center_masked_k2p5_fraction"] == 1
    assert result["unique_5deg_sky_blocks"] == 1
