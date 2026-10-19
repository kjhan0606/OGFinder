import csv
import json

import pytest

from summarize_psf_candidate_coverage import sha256, summarize


def _prepared(tmp_path):
    prepared = tmp_path / "prepared"
    prepared.mkdir()
    (prepared / "summary.json").write_text(json.dumps({
        "status": "complete", "science_measurement": False,
    }))
    records = []
    hashes = {}
    for ident, filt, n, fwhm in (("1", "i", 2, 5.0), ("2", "i2", 0, None)):
        folder = prepared / f"cluster_{ident}"
        folder.mkdir()
        path = folder / "psf_candidates.tsv"
        with path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, ["FWHM_MOMENT_PIX", "SNR_NOMINAL",
                "WING_VALID_FRACTION", "CLASSIFICATION"], delimiter="\t")
            writer.writeheader()
            if n:
                for _ in range(n):
                    writer.writerow({"FWHM_MOMENT_PIX": fwhm, "SNR_NOMINAL": 60,
                                     "WING_VALID_FRACTION": 0.9,
                                     "CLASSIFICATION": "provisional_compact_candidate"})
        hashes[str(path)] = sha256(path)
        records.append({"id": ident, "filter": filt, "n_psf_candidates": n})
    (prepared / "artifact_audit.json").write_text(json.dumps({
        "passed": True, "science_measurement": False, "mask_ready_for_science": False,
        "n_frames": 2, "records": records, "table_sha256": hashes,
    }))
    return prepared


def test_summarize_verifies_candidate_tables_and_keeps_psf_unapproved(tmp_path):
    result = summarize(_prepared(tmp_path), expected_frames=2)

    assert result["status"] == "complete_descriptive"
    assert result["science_measurement"] is False
    assert result["psf_approved"] is False
    assert result["n_frames_verified"] == 2
    assert result["n_candidate_rows"] == 2
    assert result["classification_counts"] == {"provisional_compact_candidate": 2}
    assert result["by_filter"]["i"]["n_frames_with_at_least_4_candidates"] == 0
    assert result["by_filter"]["i2"]["n_frames_with_no_candidates"] == 1


def test_summarize_rejects_changed_candidate_table(tmp_path):
    prepared = _prepared(tmp_path)
    candidate = prepared / "cluster_1" / "psf_candidates.tsv"
    candidate.write_text(candidate.read_text() + "# changed\n")
    with pytest.raises(ValueError, match="missing or changed"):
        summarize(prepared, expected_frames=2)
