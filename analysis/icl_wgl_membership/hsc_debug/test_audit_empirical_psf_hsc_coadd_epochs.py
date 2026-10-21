from astropy.io import fits

import csv

from audit_empirical_psf_hsc_coadd_epochs import audit_fits_header, run


def test_epoch_audit_finds_standard_date_in_nonprimary_hdu(tmp_path):
    path = tmp_path / "with_epoch.fits"
    primary = fits.PrimaryHDU()
    image = fits.ImageHDU()
    image.header["MJD-MID"] = 59321.25
    fits.HDUList([primary, image]).writeto(path)

    result = audit_fits_header(path, "1")

    assert result["has_standard_epoch"] is True
    assert '"key":"MJD-MID"' in result["epoch_cards"]
    assert '"hdu":1' in result["epoch_cards"]


def test_epoch_audit_does_not_mistake_exposure_time_for_an_epoch(tmp_path):
    path = tmp_path / "without_epoch.fits"
    primary = fits.PrimaryHDU()
    primary.header["EXPTIME"] = 30.0
    fits.HDUList([primary]).writeto(path)

    result = audit_fits_header(path, "2")

    assert result["has_standard_epoch"] is False
    assert result["epoch_cards"] == "[]"


def test_run_counts_unique_frames_and_records_header_only_scope(tmp_path):
    sample_root = tmp_path / "sample"
    image_dir = sample_root / "download_01/images"
    image_dir.mkdir(parents=True)
    ledger = tmp_path / "star_ledger.tsv"
    with ledger.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["cluster_id", "stamp"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"cluster_id": "1", "stamp": "/some/star_1.fits"})
        writer.writerow({"cluster_id": "1", "stamp": "/some/star_2.fits"})
        writer.writerow({"cluster_id": "2", "stamp": "/some/star_3.fits"})
    fits.PrimaryHDU().writeto(image_dir / "cluster_1.fits")
    hdu = fits.PrimaryHDU()
    hdu.header["DATE-OBS"] = "2020-01-01T00:00:00"
    hdu.writeto(image_dir / "cluster_2.fits")

    result = run(sample_root, ledger, tmp_path / "audit")

    assert result["n_gaia_stamps"] == 3
    assert result["n_unique_target_frames"] == 2
    assert result["n_frames_with_standard_epoch"] == 1
    assert result["n_frames_missing_standard_epoch"] == 1
    assert result["image_pixels_read"] is False
    assert result["hsc_requests_made"] is False
    assert result["stamp_center_source"] == "HSC detection X_IMAGE/Y_IMAGE in native coadd pixels"
    assert result["gaia_epoch_propagation_required_to_recenter_stamp"] is False
    assert result["gaia_match_identity_epoch_exactly_verified"] is False
    assert result["science_measurement"] is False
