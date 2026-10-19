#!/usr/bin/env python3
import pytest

from audit_icl_satellite_association_selection import summarize_groups


def row(cluster_id, filt, z, richness, eligible):
    return {
        "id": str(cluster_id), "filter": filt, "z": str(z),
        "richness": str(richness), "delta_sb_njy_arcsec2": str(z + richness),
        "target_coverage": "0.8", "control_coverage": "0.7",
        "common_coverage": "0.7", "n_member_positions": "3",
        "member_density_mpc2": "20", "auto_unique_center_candidate": str(eligible),
        "center_masked_k2p5": "False",
    }


def test_selection_audit_reports_numeric_shift_and_unique_parent_keys():
    parent = [row(i, "i2", i / 10, 10 + i, i % 2 == 0) for i in range(8)]
    selected = parent[:4]
    result = summarize_groups(parent, selected)
    assert result["n_icl_support_rows"] == 8
    assert result["n_association_rows"] == 4
    assert result["association_fraction"] == 0.5
    assert result["numeric"]["richness"]["standardized_mean_difference"] < 0
    assert result["center_eligible_fraction"]["association"] == 0.5


def test_selection_audit_rejects_duplicate_and_unmatched_ids():
    parent = [row(i, "i", i / 10, i + 1, True) for i in range(6)]
    with pytest.raises(ValueError, match="duplicated"):
        summarize_groups(parent, [parent[0], parent[0]])
    with pytest.raises(ValueError, match="match"):
        summarize_groups(parent, [row(99, "i", 0.4, 100, True)])
