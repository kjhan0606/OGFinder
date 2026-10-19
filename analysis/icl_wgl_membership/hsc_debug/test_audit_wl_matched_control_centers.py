from audit_wl_matched_control_centers import control_overlap_counts


def test_control_aperture_overlap_uses_great_circle_separation_and_paired_z():
    rows = [
        {"ra_deg": 359.9, "dec_deg": 0.0, "z_lens": 0.3},
        {"ra_deg": 0.1, "dec_deg": 0.0, "z_lens": 0.3},
        {"ra_deg": 2.0, "dec_deg": 0.0, "z_lens": 0.3},
    ]
    counts, radii = control_overlap_counts(rows, 70.0, 0.3)
    assert counts.tolist() == [1, 1, 0]
    assert radii[0] > 0


def test_larger_redshift_increases_angular_overlap_radius():
    rows = [
        {"ra_deg": 0.0, "dec_deg": 0.0, "z_lens": 0.3},
        {"ra_deg": 0.0, "dec_deg": 0.0, "z_lens": 0.8},
    ]
    _, radii = control_overlap_counts(rows, 70.0, 0.3)
    assert radii[0] > radii[1]
