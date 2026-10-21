from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_public_extended_psf_reference import (
    analyze_image,
    audit_fits,
    physical_aperture_diagnostics,
)


def test_radial_profile_tracks_encircled_flux_and_center():
    image = np.zeros((9, 9), dtype=float)
    image[4, 4] = 10.0
    image[4, 5] = image[4, 3] = 1.0
    profile, metrics = analyze_image(image, pixel_scale_arcsec=0.2,
                                    radial_edges_pix=np.array([0.0, 0.5, 1.5, 3.0, 8.0]))
    assert metrics["peak_xy_pixel"] == [4, 4]
    assert metrics["peak_offset_from_center_pixel"] == 0
    assert metrics["profile_covered_flux_fraction"] == pytest.approx(1.0)
    assert profile[1]["annular_flux_fraction"] == pytest.approx(2 / 12)
    assert profile[-1]["encircled_flux_fraction_at_outer_edge"] == pytest.approx(1.0)
    assert metrics["psf_approved"] is False
    assert metrics["science_measurement"] is False


@pytest.mark.parametrize("image", [np.zeros((5, 5)), np.ones((4, 4)), np.full((5, 5), np.nan)])
def test_rejects_invalid_psf_arrays(image):
    with pytest.raises(ValueError):
        analyze_image(image, radial_edges_pix=np.array([1.0, 1.0]))


def test_audit_reads_second_fits_extension_and_hashes(tmp_path):
    data = np.zeros((9, 9), dtype=np.float32)
    data[4, 4] = 1
    path = tmp_path / "reference.fits"
    fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(data)]).writeto(path)
    metadata, profile = audit_fits("sym", path, 0.168)
    assert metadata["shape"] == [9, 9]
    assert metadata["bytes"] == path.stat().st_size
    assert len(metadata["sha256"]) == 64
    assert metadata["hdu_shapes"] == [None, [9, 9]]
    assert len(profile) == 49


def test_physical_annulus_diagnostic_has_eight_normalized_sectors():
    image = np.ones((1001, 1001), dtype=float)
    result = physical_aperture_diagnostics(image, redshifts=(0.3,))[0]
    annulus = result["annuli"]["icl_200_300_kpc"]
    assert result["proper_kpc_per_pixel"] > 0
    assert len(annulus["eight_sector_flux_shares"]) == 8
    assert sum(annulus["eight_sector_flux_shares"]) == pytest.approx(1.0)
    assert annulus["sector_share_rms_about_equal"] < 0.01
