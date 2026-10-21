import pytest

from audit_icl_satellite_radial_support_matrix import support_matrix


def _icl(cluster, band, radial):
    return {"id": cluster, "filter": band, "radial_bin_kpc": radial}


def _sat(cluster, band, radial):
    return {"cluster_id": cluster, "filter": band, "radial_bin_kpc": radial,
            "cohort": "full_sample"}


def test_support_matrix_counts_exact_cluster_filter_intersections():
    icl = [_icl("A", "i", "200-300"), _icl("B", "i", "200-300"),
           _icl("A", "i2", "200-300")]
    sat = [_sat("A", "i", "50-75"), _sat("C", "i", "50-75"),
           _sat("A", "i2", "50-75")]

    rows = support_matrix(icl, sat)

    i_outer = next(row for row in rows if row["filter"] == "i"
                   and row["icl_radial_bin_kpc"] == "200-300"
                   and row["satellite_radial_bin_kpc"] == "50-75")
    i2_outer = next(row for row in rows if row["filter"] == "i2"
                    and row["icl_radial_bin_kpc"] == "200-300"
                    and row["satellite_radial_bin_kpc"] == "50-75")
    assert i_outer["n_icl_supported"] == 2
    assert i_outer["n_satellite_supported"] == 2
    assert i_outer["n_common"] == 1
    assert i_outer["common_fraction_of_icl"] == 0.5
    assert i2_outer["n_common"] == 1


def test_support_matrix_excludes_non_full_sample_satellite_rows():
    rows = support_matrix([_icl("A", "i", "50-100")],
                          [{**_sat("A", "i", "5-10"), "cohort": "pilot"}])
    cell = next(row for row in rows if row["filter"] == "i"
                and row["icl_radial_bin_kpc"] == "50-100"
                and row["satellite_radial_bin_kpc"] == "5-10")
    assert cell["n_satellite_supported"] == 0
    assert cell["n_common"] == 0


def test_support_matrix_rejects_duplicate_profile_keys():
    with pytest.raises(ValueError, match="duplicate ICL"):
        support_matrix([_icl("A", "i", "50-100"), _icl("A", "i", "50-100")], [])
