import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).with_name("audit_s23b_candidate_overlap_association_support.py")
SPEC = importlib.util.spec_from_file_location("audit_s23b_candidate_overlap_association_support", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_classifies_icl_support_and_joint_association_independently():
    matches = []
    effective = []
    preexcluded = []
    for index in range(7):
        target = f"t{index}"
        row = {"target_id": target, "control_id": f"c{index}", "filter": "i",
               "candidate_ids": [f"obj{index}"]}
        matches.append({"target_cluster_id": target, "frame_id": f"c{index}",
                        "filter": "i", "candidates": [{"candidate_id": f"obj{index}"}]})
        (preexcluded if index == 0 else effective).append(row)
    stack_summary = {"s23b_candidate_control_exclusions": effective,
                     "candidate_controls_already_filter_excluded": preexcluded}
    icl_rows = [{"id": "t1", "filter": "i", "radial_bin_kpc": "200-300",
                 "control_id": "c1", "delta_sb_njy_arcsec2": "1.5"}]
    association_rows = [{"cluster_id": "t1", "filter": "i",
                         "satellite_delta_50_75": "0.25"}]
    satellite_rows = [{"cluster_id": "t1", "filter": "i", "radial_bin_kpc": "50-75"}]

    rows = audit.build_records(matches, stack_summary, icl_rows, association_rows, satellite_rows)

    assert len(rows) == 7
    assert rows[0]["baseline_pair_eligible"] is False
    assert rows[1]["has_icl_200_300_support"] is True
    assert rows[1]["has_icl_satellite_association_row"] is True
    assert rows[1]["has_satellite_50_75_profiles"] is True
    assert rows[1]["n_satellite_profiles_50_75"] == 1
    assert rows[2]["has_icl_200_300_support"] is False
    with pytest.raises(ValueError, match="different control"):
        audit.build_records(matches, stack_summary,
                            [{**icl_rows[0], "control_id": "wrong"}], association_rows,
                            satellite_rows)
