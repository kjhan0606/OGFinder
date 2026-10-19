import json
import pytest

from plot_wl_tomographic_lens_redshift import load_profiles, plot


def _summary():
    profiles = []
    radial = ((0.1, 0.2), (0.2, 0.4), (0.4, 0.7), (0.7, 1.0),
              (1.0, 1.5), (1.5, 2.0), (2.0, 3.0))
    for cohort in ("lower_z", "upper_z"):
        for source_bin in range(1, 5):
            for index, (r_min, r_max) in enumerate(radial):
                profiles.append({
                    "lens_redshift_cohort": cohort,
                    "source_zbin": source_bin,
                    "r_min_mpc": r_min,
                    "r_max_mpc": r_max,
                    "delta_sigma_weak_msun_mpc2": float(source_bin * (index + 1)),
                    "delta_sigma_jackknife_error_msun_mpc2": 0.5,
                    "cross_delta_sigma_weak_msun_mpc2": 0.1,
                    "cross_jackknife_error_msun_mpc2": 0.4,
                })
    return {
        "schema": "astrafek-hsc-y3-lens-redshift-tomographic-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "lens_redshift_split": {
            "n_lenses_lower_z": 400,
            "n_lenses_upper_z": 400,
            "lower_z_max": 0.43434,
            "upper_z_min": 0.43519,
            "midpoint_between_boundary_lenses": 0.434765,
        },
        "profile": profiles,
    }


def test_lens_redshift_plot_rejects_incomplete_profiles(tmp_path):
    summary = _summary()
    summary["profile"].pop()
    path = tmp_path / "summary.json"
    path.write_text(json.dumps(summary))
    with pytest.raises(ValueError, match="expected 2 lens-redshift cohorts"):
        load_profiles(path)


def test_plot_lens_redshift_tomography_writes_png(tmp_path):
    summary_path = tmp_path / "summary.json"
    summary_path.write_text(json.dumps(_summary()))
    _, profiles = load_profiles(summary_path)
    assert len(profiles["lower_z"]) == 4
    assert len(profiles["upper_z"]) == 4
    output = tmp_path / "profiles.png"
    plot(summary_path, output)
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
