import importlib.util
from pathlib import Path


MODULE_PATH = Path(__file__).with_name("audit_icl_center_quality.py")
SPEC = importlib.util.spec_from_file_location("audit_icl_center_quality", MODULE_PATH)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def test_center_flag_requires_unique_candidate_and_unmasked_center():
    assert audit.center_eligible({"auto_unique_center_candidate": "True",
                                  "center_masked_by_kron": "False"})
    assert not audit.center_eligible({"auto_unique_center_candidate": "False",
                                     "center_masked_by_kron": "False"})
    assert not audit.center_eligible({"auto_unique_center_candidate": "True",
                                     "center_masked_by_kron": "True"})


def test_block_bootstrap_summarizes_paired_spatial_groups():
    rows = []
    for i in range(12):
        ra = i * 10 + 1
        rows.append({"target_ra_deg": str(ra), "target_dec_deg": "0",
                     "control_ra_deg": str(ra + .2), "control_dec_deg": "0",
                     "delta_flux": str(i)})
    result = audit.summarize(rows, replicas=1000, seed=11)
    assert result["n_clusters"] == 12
    assert result["n_joint_blocks"] == 12
    assert result["mean_delta_flux"] == 5.5
    assert result["spatial_block_bootstrap_ci95"][0] < 5.5
    assert result["spatial_block_bootstrap_ci95"][1] > 5.5
    assert result["leave_one_block_out_mean_range"] == [5.0, 6.0]
    assert result["leave_one_block_out_max_abs_shift"] == .5
