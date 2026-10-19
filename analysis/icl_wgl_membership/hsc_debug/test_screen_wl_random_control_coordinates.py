import numpy as np

from screen_wl_random_control_coordinates import screen_geometry


def test_geometry_screen_handles_ra_wrap_and_flags_duplicate_candidate():
    lenses = [
        {"lens_id": "L1", "z_lens": 0.3},
        {"lens_id": "L2", "z_lens": 0.4},
    ]
    coordinates = [
        {"lens_id": "L1", "ra_deg": 359.99, "dec_deg": 0.0},
        {"lens_id": "L1", "ra_deg": 1.0, "dec_deg": 0.0},
        {"lens_id": "L2", "ra_deg": 1.0, "dec_deg": 0.0},
    ]
    camira = np.array([[1.0, 0.01, 0.0, 0.3]])

    rows = screen_geometry(coordinates, lenses, camira, h0=70, omega_m=0.3)

    assert rows[0]["camira_3mpc_overlap_count"] == 1
    assert not rows[0]["geometry_candidate"]
    assert rows[1]["camira_3mpc_overlap_count"] == 0
    assert rows[1]["coordinate_multiplicity"] == 2
    assert rows[2]["coordinate_multiplicity"] == 2
    assert not rows[1]["geometry_candidate"]
    assert not rows[2]["geometry_candidate"]
    assert all(row["science_ready"] is False for row in rows)


def test_geometry_screen_rejects_unknown_lens_and_bad_chunk_size():
    lenses = [{"lens_id": "L1", "z_lens": 0.3}]
    camira = np.array([[1.0, 100.0, 0.0, 0.3]])
    with np.testing.assert_raises_regex(ValueError, "unknown lens"):
        screen_geometry([{"lens_id": "L2", "ra_deg": 1.0, "dec_deg": 0.0}],
                        lenses, camira, 70, 0.3)
    with np.testing.assert_raises_regex(ValueError, "chunk size"):
        screen_geometry([], lenses, camira, 70, 0.3, chunk_size=0)
