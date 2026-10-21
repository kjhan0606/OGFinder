import pytest

from analyze_icl_satellite_weighted_radial_matrix_familywise import aggregate_member_profiles


def member(cluster, index, delta, weight, radial="50-75", band="i2"):
    return {"cluster_id": cluster, "filter": band, "radial_bin_kpc": radial,
            "member_index": str(index), "delta": str(delta), "camira_weight": str(weight)}


def test_member_profiles_apply_camira_weight_with_equal_cluster_unit():
    rows = [member("A", 0, 1, 0.25), member("A", 1, 3, 0.75),
            member("B", 0, 5, 1.0)]
    profiles = aggregate_member_profiles(rows)
    assert profiles[("A", "i2", "50-75")]["equal_mean"] == pytest.approx(2.0)
    assert profiles[("A", "i2", "50-75")]["weighted_mean"] == pytest.approx(2.5)
    assert profiles[("A", "i2", "50-75")]["n_satellites"] == 2
    assert profiles[("B", "i2", "50-75")]["weighted_mean"] == pytest.approx(5.0)


def test_zero_total_membership_weight_is_rejected():
    with pytest.raises(ValueError, match="zero total CAMIRA weight"):
        aggregate_member_profiles([member("A", 0, 1, 0.0)])


def test_duplicate_member_annulus_is_rejected():
    row = member("A", 0, 1, 1.0)
    with pytest.raises(ValueError, match="duplicate CAMIRA-weighted member profile"):
        aggregate_member_profiles([row, row])
