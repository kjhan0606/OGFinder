from prepare_s23b_icl_coverage_plan import OUTPUT_COLUMNS, build_probes
from prepare_icl_expansion_coverage_plan import render_query
from ogfkit.wlcluster import Cosmology


def test_build_probes_uses_release_prefixed_id_and_nine_positions():
    cosmology = Cosmology(H0=70.0, Om=0.3)
    rows = build_probes([{
        "id": "S23B:42", "ra_deg": 359.9, "dec_deg": 12.0, "z": 0.4,
    }], cosmology)
    assert len(rows) == 9
    assert {row["candidate_id"] for row in rows} == {"S23B:42"}
    assert {row["probe_id"] for row in rows} == {
        "center", *(f"r1Mpc_pa{angle:03d}" for angle in range(0, 360, 45))
    }
    assert next(row for row in rows if row["probe_id"] == "center")["radius_kpc"] == 0.0
    assert all(row["radius_kpc"] == 1000.0 for row in rows if row["probe_id"] != "center")


def test_query_schema_is_metadata_only_and_hsc_i_limited():
    sql = render_query([{
        "candidate_id": "S23B:42", "probe_id": "center",
        "ra_deg": 12.0, "dec_deg": -1.0,
    }]).upper()
    assert "PDR3_WIDE.MOSAIC" in sql
    assert "PUBLIC.SKYMAP" in sql
    assert "PDR3_WIDE.PATCH_QA" in sql
    assert "M.BAND = 'HSC-I'" in sql
    assert "SELECT *" not in sql
    assert "SHAPE" not in sql and "PHOTOZ" not in sql
    assert OUTPUT_COLUMNS == [
        "candidate_id", "probe_id", "skymap_id", "tract", "patch", "band",
        "naxis1", "naxis2", "seeing", "zeropt", "imag_psf_depth", "iseeing",
        "imeas_exists", "idet_bkgd_exists",
    ]
