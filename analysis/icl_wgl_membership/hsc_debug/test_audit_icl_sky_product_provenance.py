import pytest

from audit_icl_sky_product_provenance import validate_collection


def make_collection(tmp_path):
    images = tmp_path / "images"
    images.mkdir()
    image = images / "frame.fits"
    image.write_bytes(b"synthetic fits placeholder")
    frame = {
        "id": "C1", "image": image.name, "sky": "global-sky",
        "product": "coadd/bg", "product_reviewed": True,
        "product_evidence": "HSC manual maps product before aggressive background subtraction",
        "coadd_header_filter": "i",
    }
    manifest = {"clusters": [frame], "downloaded_images": str(images)}
    status = {
        "status": "complete", "n_expected": 1, "downloaded": 1,
        "credentials_saved": False,
        "completed": [{"id": "C1", "file": image.name,
                       "bytes": image.stat().st_size, "sha256": "a" * 64}],
    }
    return manifest, status


def test_collection_accepts_reviewed_global_sky_coadd_bg(tmp_path):
    manifest, status = make_collection(tmp_path)
    result = validate_collection("target", manifest, status, expected_count=1)
    assert result["n_frames"] == 1
    assert result["products"] == {"coadd/bg": 1}
    assert result["sky_types"] == {"global-sky": 1}
    assert result["recorded_image_bytes"] == status["completed"][0]["bytes"]


def test_collection_rejects_local_sky_or_wrong_product(tmp_path):
    manifest, status = make_collection(tmp_path)
    manifest["clusters"][0]["sky"] = "local-sky"
    with pytest.raises(ValueError, match="unexpected sky/product"):
        validate_collection("target", manifest, status, expected_count=1)


def test_collection_rejects_missing_file_or_size_mismatch(tmp_path):
    manifest, status = make_collection(tmp_path)
    status["completed"][0]["bytes"] += 1
    with pytest.raises(ValueError, match="missing or size differs"):
        validate_collection("target", manifest, status, expected_count=1)
