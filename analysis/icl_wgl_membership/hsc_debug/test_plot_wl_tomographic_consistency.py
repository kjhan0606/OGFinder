import json

from plot_wl_tomographic_consistency import plot


def test_plot_tomographic_profiles_writes_png(tmp_path):
    profiles = []
    for source_bin in range(1, 5):
        for index, (r_min, r_max) in enumerate(((0.1, 0.2), (0.2, 0.4), (0.4, 0.7), (0.7, 1.0), (1.0, 1.5), (1.5, 2.0), (2.0, 3.0))):
            profiles.append({
                "source_zbin": source_bin,
                "r_min_mpc": r_min,
                "r_max_mpc": r_max,
                "delta_sigma_weak_msun_mpc2": float(source_bin * (index + 1)),
                "delta_sigma_jackknife_error_msun_mpc2": 0.5,
                "cross_delta_sigma_weak_msun_mpc2": 0.1,
                "cross_jackknife_error_msun_mpc2": 0.4,
            })
    summary = tmp_path / "summary.json"
    summary.write_text(json.dumps({
        "schema": "astrafek-hsc-y3-tomographic-consistency-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "profile": profiles,
    }))
    output = tmp_path / "profiles.png"
    plot(summary, output)
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
