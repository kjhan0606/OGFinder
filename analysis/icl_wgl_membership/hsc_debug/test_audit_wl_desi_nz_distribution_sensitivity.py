from __future__ import annotations

import json
import hashlib

import numpy as np
import pytest

from audit_wl_desi_nz_distribution_sensitivity import (
    CURVE_KINDS,
    SCALE_CUTS,
    load_verified_baseline,
    load_desi_curves,
    profile_comparison,
)


def _write_variant(directory, variant, *, wrong_scale=False, missing_key=None):
    folder = directory / "results" / "distributions" / variant
    folder.mkdir(parents=True)
    lo, hi = SCALE_CUTS[variant]
    metadata = {"scale_cuts": [lo + 0.1, hi] if wrong_scale else [lo, hi]}
    (folder / f"nz_res_metadata_{variant}.json").write_text(json.dumps(metadata))
    z = np.linspace(0.0, 3.0, 101)
    arrays = {"z": z}
    for source_bin in range(1, 5):
        curve = np.exp(-0.5 * ((z - (0.5 + source_bin * 0.2)) / 0.2) ** 2)
        for kind in CURVE_KINDS:
            if missing_key == (source_bin, kind):
                continue
            arrays[f"{source_bin}/nz_{kind}"] = curve * (1 + 0.01 * CURVE_KINDS.index(kind))
    np.savez(folder / f"distributions_{variant}.npz", **arrays)


def _write_manifest(directory):
    commit = "abc70af29910e49fce17a63ef3083160d515fb6c"
    rows = []
    for variant in SCALE_CUTS:
        for filename in (f"distributions_{variant}.npz", f"nz_res_metadata_{variant}.json"):
            rel = f"results/distributions/{variant}/{filename}"
            path = directory / rel
            data = path.read_bytes()
            rows.append({
                "path": rel,
                "url": f"https://raw.githubusercontent.com/JeanCHDJdev/desi-y3-hsc/{commit}/{rel}",
                "bytes": len(data),
                "sha256": hashlib.sha256(data).hexdigest(),
            })
    (directory / "manifest.json").write_text(json.dumps({"commit": commit, "files": rows}))


def test_loader_checks_metadata_grids_and_normalizes_all_curve_kinds(tmp_path):
    for variant in SCALE_CUTS:
        _write_variant(tmp_path, variant)
    _write_manifest(tmp_path)
    z, curves, hashes = load_desi_curves(tmp_path)
    assert z.shape == (101,)
    assert set(curves) == set(SCALE_CUTS)
    assert set(curves["0.3_3_v_1p1"]) == set(CURVE_KINDS)
    for by_kind in curves.values():
        for array in by_kind.values():
            np.testing.assert_allclose(np.trapezoid(array, z, axis=1), np.ones(4))
    assert len(hashes) == 5  # four approved artifacts plus their provenance manifest


def test_loader_rejects_scale_cut_metadata_mismatch(tmp_path):
    _write_variant(tmp_path, "0.3_3_v_1p1", wrong_scale=True)
    _write_variant(tmp_path, "1_5_v_1p1")
    _write_manifest(tmp_path)
    with pytest.raises(ValueError, match="scale cut"):
        load_desi_curves(tmp_path)


def test_loader_rejects_missing_required_curve(tmp_path):
    _write_variant(tmp_path, "0.3_3_v_1p1", missing_key=(2, "upper_mu"))
    _write_variant(tmp_path, "1_5_v_1p1")
    _write_manifest(tmp_path)
    with pytest.raises(ValueError, match="missing required"):
        load_desi_curves(tmp_path)


def test_profile_comparison_pairs_source_bin_and_radius():
    baseline = {"profile": [{
        "source_zbin": 2, "r_min_mpc": 0.1, "r_max_mpc": 0.2,
        "delta_sigma_weak_msun_mpc2": 10.0,
        "delta_sigma_jackknife_error_msun_mpc2": 2.0,
    }]}
    candidate = {"profile": [{
        "source_zbin": 2, "r_min_mpc": 0.1, "r_max_mpc": 0.2,
        "delta_sigma_weak_msun_mpc2": 12.0,
        "delta_sigma_jackknife_error_msun_mpc2": 3.0,
    }]}
    result = profile_comparison(baseline, candidate)
    assert result[0]["change_msun_mpc2"] == 2.0
    assert result[0]["fractional_change"] == pytest.approx(0.2)
    candidate["profile"][0]["source_zbin"] = 3
    with pytest.raises(ValueError, match="different"):
        profile_comparison(baseline, candidate)


def test_profile_comparison_accepts_envelope_without_jackknife_error():
    base_row = {
        "source_zbin": 1, "r_min_mpc": 0.1, "r_max_mpc": 0.2,
        "delta_sigma_weak_msun_mpc2": 4.0,
        "delta_sigma_jackknife_error_msun_mpc2": 1.0,
    }
    envelope_row = {
        "source_zbin": 1, "r_min_mpc": 0.1, "r_max_mpc": 0.2,
        "delta_sigma_weak_msun_mpc2": 4.5,
    }
    comparison = profile_comparison({"profile": [base_row]}, {"profile": [envelope_row]})
    assert comparison[0]["desi_spatial_jackknife_error_msun_mpc2"] is None


def test_persisted_baseline_requires_exact_provenance(tmp_path):
    sacc = tmp_path / "baseline.sacc"
    sacc.write_bytes(b"approved test sacc")
    baseline_path = tmp_path / "baseline.json"
    provenance = {
        "row_plan_sha256": "plan", "row_batch_sha256": {"batch": "hash"},
        "n_pair_rows": 10, "n_lenses": 800, "masked_value": 0,
    }
    payload = {
        "schema": "astrafek-hsc-y3-tomographic-consistency-jackknife-v1",
        "science_measurement": False, "random_lens_subtraction_applied": False,
        "n_occupied_cells": 52, "profile": [{"source_zbin": 1}],
        "cell_width_equatorial_deg": 5.0,
        "input": {
            **provenance, "n_z_sacc_sha256": __import__("hashlib").sha256(sacc.read_bytes()).hexdigest(),
            "radial_edges_mpc": [0.1, 0.2],
            "cosmology": {"H0_km_s_mpc": 70.0, "Omega_m": 0.3},
        },
    }
    baseline_path.write_text(json.dumps(payload))
    baseline, digest = load_verified_baseline(
        baseline_path, provenance, sacc, np.array([0.1, 0.2]), 5.0, 70.0, 0.3,
    )
    assert baseline["n_occupied_cells"] == 52
    assert len(digest) == 64
    provenance["masked_value"] = 1
    with pytest.raises(ValueError, match="masked_value"):
        load_verified_baseline(
            baseline_path, provenance, sacc, np.array([0.1, 0.2]), 5.0, 70.0, 0.3,
        )
