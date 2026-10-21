import numpy as np

from stack_empirical_psf_star_profiles import frame_weighted_stack, measure_stamp


def test_frame_weighting_avoids_pseudoreplication_and_bootstrap_is_reproducible():
    rows = []
    for frame, values in (("A", [1.0, 3.0, 5.0]), ("B", [9.0])):
        for i, value in enumerate(values):
            rows.append({"candidate_id": f"{frame}{i}", "frame_id": frame, "kron_factor": 2.5,
                         "sky_r_min_pix": 120, "sky_r_max_pix": 150, "r_min_pix": 0,
                         "r_max_pix": 2, "r_mid_arcsec": 0.168, "status": "ok",
                         "profile_per_core_per_pixel": value})
    a = frame_weighted_stack(rows, n_bootstrap=1000, seed=7)[0]
    b = frame_weighted_stack(rows, n_bootstrap=1000, seed=7)[0]
    assert a["n_stars"] == 4 and a["n_frames"] == 2
    assert a["frame_equal_mean"] == 6.0  # frame A mean=3, frame B=9
    assert a["frame_equal_mean"] == b["frame_equal_mean"]
    assert a["bootstrap_mean_q16"] == b["bootstrap_mean_q16"]


def test_stamp_profile_preserves_signed_wings_and_science_gates(tmp_path):
    from astropy.io import fits

    size = 401
    yy, xx = np.indices((size, size), dtype=float)
    radius = np.hypot(xx - 200.25, yy - 199.75)
    image = 100.0 * np.exp(-radius / 2.0) - 0.01 - 0.01 * ((radius >= 64) & (radius < 96))
    variance = np.ones((size, size), dtype=float)
    base = np.zeros((size, size), dtype=np.uint8)
    mask = np.zeros((size, size), dtype=np.uint8)
    header = fits.Header()
    header["SCIVAL"] = False; header["PSFOK"] = False
    header["XCENTER"] = 200.25; header["YCENTER"] = 199.75; header["PIXSCALE"] = 0.168
    path = tmp_path / "star.fits"
    fits.HDUList([fits.PrimaryHDU(header=header), fits.ImageHDU(image, name="SCIENCE"),
                  fits.ImageHDU(variance, name="VARIANCE"), fits.ImageHDU(base, name="BASE_MASK"),
                  fits.ImageHDU(mask, name="MASK_K2P5"), fits.ImageHDU(mask, name="MASK_K3P5")]).writeto(path)
    rows = measure_stamp(path, "frame", "candidate")
    outer = [r for r in rows if r["r_min_pix"] == 64 and r["r_max_pix"] == 96 and r["kron_factor"] == 2.5]
    assert len(outer) == 2 and all(r["status"] == "ok" for r in outer)
    assert all(r["profile_per_core_per_pixel"] < 0 for r in outer)
    assert all(r["science_measurement"] is False and r["psf_approved"] is False for r in outer)
