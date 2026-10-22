import pytest

from submit_icl_expansion_coverage_metadata import OUTPUT_COLUMNS, validate_result
from validate_icl_expansion_coverage_metadata import aggregate_candidates


def make_probes(candidate_id="10"):
    return [{"candidate_id": candidate_id, "probe_id": "center"}] + [
        {"candidate_id": candidate_id, "probe_id": f"r1Mpc_pa{angle:03d}"}
        for angle in range(0, 360, 45)
    ]


def make_row(probe_id, *, candidate_id="10", depth="26.0", seeing="0.7", imeas="true", idet="true"):
    return {
        "candidate_id": candidate_id, "probe_id": probe_id, "skymap_id": f"map-{probe_id}",
        "tract": "100", "patch": "4,4", "band": "HSC-I", "naxis1": "4200",
        "naxis2": "4200", "seeing": "0.8", "zeropt": "27.0",
        "imag_psf_depth": depth, "iseeing": seeing,
        "imeas_exists": imeas, "idet_bkgd_exists": idet,
    }


def test_aggregate_reports_complete_nine_point_coverage_and_metadata_only_metrics():
    candidate = {"id": "10", "ra_deg": "12.3", "dec_deg": "4.5", "z_cl": "0.4", "richness": "20"}
    probes = make_probes()
    rows = {("10", row["probe_id"]): [make_row(row["probe_id"])] for row in probes}
    result = aggregate_candidates([candidate], probes, rows)[0]
    assert result["coverage_probes"] == 9
    assert result["all_nine_probes_covered"] is True
    assert result["median_best_imag_psf_depth"] == 26.0
    assert result["median_best_iseeing"] == 0.7
    assert result["imeas_exists_true_fraction"] == 1.0
    assert result["science_measurement"] is False


def test_aggregate_counts_missing_probes_as_uncovered():
    probes = make_probes()
    rows = {("10", row["probe_id"]): [make_row(row["probe_id"])] for row in probes[:-1]}
    result = aggregate_candidates(
        [{"id": "10", "ra_deg": "1", "dec_deg": "2", "z_cl": "0.3", "richness": "5"}],
        probes, rows,
    )[0]
    assert result["coverage_probes"] == 8
    assert result["all_nine_probes_covered"] is False
    assert result["probe_coverage_fraction"] == pytest.approx(8 / 9)
    assert result["idet_bkgd_exists_true_fraction"] == 1.0


def test_aggregate_accepts_release_prefixed_s23b_candidate_key():
    candidate_id = "S23B:10"
    probes = make_probes(candidate_id)
    rows = {
        (candidate_id, row["probe_id"]): [make_row(row["probe_id"], candidate_id=candidate_id)]
        for row in probes
    }
    result = aggregate_candidates(
        [{"candidate_id": candidate_id, "ra_deg": "12.3", "dec_deg": "4.5",
          "z_cl": "0.4", "richness": "20"}], probes, rows,
    )[0]
    assert result["candidate_id"] == candidate_id
    assert result["all_nine_probes_covered"] is True


def test_aggregate_rejects_incomplete_or_duplicate_probe_inventory():
    with pytest.raises(ValueError, match="nine probes"):
        aggregate_candidates(
            [{"id": "10", "ra_deg": "1", "dec_deg": "2", "z_cl": "0.3", "richness": "5"}],
            make_probes()[:-1], {},
        )


def test_validate_result_accepts_only_the_approved_metadata_schema(tmp_path):
    path = tmp_path / "metadata.csv"
    row = make_row("center")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = __import__("csv").DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerow(row)
    assert validate_result(path, {"n_probe_points": 9}, {("10", "center")}) == 1


def test_validate_result_rejects_unrequested_probe_and_extra_fields(tmp_path):
    path = tmp_path / "metadata.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = __import__("csv").DictWriter(stream, fieldnames=OUTPUT_COLUMNS + ["shape"])
        writer.writeheader()
        writer.writerow({**make_row("center"), "shape": "not allowed"})
    with pytest.raises(ValueError, match="columns outside"):
        validate_result(path, {"n_probe_points": 9}, {("10", "center")})
