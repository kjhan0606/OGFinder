import numpy as np

import hashlib
import json

import pytest
from astropy.cosmology import FlatLambdaCDM

from screen_wl_random_control_coordinates import load_preflight_cosmology, screen_geometry


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


def test_load_preflight_cosmology_uses_hash_verified_source_manifest(tmp_path):
    sample_path = tmp_path / "sample_manifest.json"
    sample_path.write_text(json.dumps({"h0": 70.0, "omega_m": 0.3}), encoding="utf-8")
    digest = hashlib.sha256(sample_path.read_bytes()).hexdigest()
    preflight_dir = tmp_path / "preflight"
    preflight_dir.mkdir()
    (preflight_dir / "preflight_manifest.json").write_text(json.dumps({
        "source_manifest": str(sample_path), "source_sha256": digest,
    }), encoding="utf-8")

    assert load_preflight_cosmology(preflight_dir) == (70.0, 0.3)

    sample_path.write_text(json.dumps({"h0": 68.0, "omega_m": 0.31}), encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        load_preflight_cosmology(preflight_dir)


def test_overlap_uses_each_redshift_proper_radius_and_is_chunk_invariant():
    h0, omega_m = 70.0, 0.3
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    target_z, neighbor_z = 0.3, 0.65
    target_radius = 3.0 / cosmology.angular_diameter_distance(target_z).value
    neighbor_radius = 3.0 / cosmology.angular_diameter_distance(neighbor_z).value
    radius_sum_deg = np.rad2deg(target_radius + neighbor_radius)
    base_ra = 140.0
    coordinates = [
        {"lens_id": "L1", "ra_deg": base_ra + 0.99 * radius_sum_deg, "dec_deg": 0.0},
        {"lens_id": "L1", "ra_deg": base_ra + 1.01 * radius_sum_deg, "dec_deg": 0.0},
    ]
    lenses = [{"lens_id": "L1", "z_lens": target_z}]
    camira = np.array([[7.0, base_ra, 0.0, neighbor_z]])

    one_row_chunks = screen_geometry(coordinates, lenses, camira, h0, omega_m, chunk_size=1)
    one_chunk = screen_geometry(coordinates, lenses, camira, h0, omega_m, chunk_size=10)

    assert one_row_chunks == one_chunk
    assert one_row_chunks[0]["camira_3mpc_overlap_count"] == 1
    assert not one_row_chunks[0]["geometry_candidate"]
    assert one_row_chunks[1]["camira_3mpc_overlap_count"] == 0
    assert one_row_chunks[1]["geometry_candidate"]
