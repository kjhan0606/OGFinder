#!/usr/bin/env bash
# CPU-only, approved PSF-core injection followed by matched baseline comparison.
set -euo pipefail
if [[ $# -ne 3 ]]; then
    echo 'Usage: run_s23b_member_psf_injection_and_compare.sh PSF_PLAN_DIR AUDITED_PSF_DIR NEW_RUN_DIR' >&2
    exit 2
fi
case "$(hostname | tr '[:upper:]' '[:lower:]')" in
    lageunha*) ;;
    *) echo 'Run this CPU calculation on lageunha.' >&2; exit 2 ;;
esac
plan_dir=$1
psf_dir=$2
run_root=$3
for path in "$plan_dir" "$psf_dir" "$run_root"; do
    if [[ "$path" != /* ]]; then echo 'All paths must be absolute.' >&2; exit 2; fi
done
plan_hash=03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d
actual_hash=$(sha256sum "$plan_dir/psf_request_plan.json")
if [[ "${actual_hash%% *}" != "$plan_hash" ]]; then
    echo 'Approved PSF plan hash mismatch.' >&2; exit 2
fi
if [[ ! -f "$psf_dir/psf_picker_audit.json" ]]; then
    echo 'The approved PSF receipt must pass its auditor first.' >&2; exit 2
fi
mkdir -- "$run_root"
exec >"$run_root/run.log" 2>&1
printf '%s\n' "$$" >"$run_root/runner.pid"
printf '%s\n' RUNNING >"$run_root/status"
trap 'rc=$?; if [[ $rc -ne 0 ]]; then printf "FAILED exit=%s\n" "$rc" >"$run_root/status"; fi' EXIT
analysis_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_root=$(cd -- "$analysis_dir/../../.." && pwd)
cd -- "$repo_root"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python_bin=/home/kjhan/miniconda3/bin/python3
{
    hostname
    date -u '+%Y-%m-%dT%H:%M:%SZ'
    git rev-parse HEAD
    sha256sum "$plan_dir/psf_request_plan.json" "$psf_dir/psf_picker_audit.json" \
        "$analysis_dir/run_s23b_member_outer_profile_injection_recovery.py" \
        "$analysis_dir/compare_s23b_member_psf_injection.py" "${BASH_SOURCE[0]}"
} >"$run_root/provenance.txt"
"$python_bin" "$analysis_dir/run_s23b_member_outer_profile_injection_recovery.py" \
    --workers 1 --psf-plan-dir "$plan_dir" --psf-download-dir "$psf_dir" \
    --approved-psf-plan-sha256 "$plan_hash" --output "$run_root/injection"
"$python_bin" "$analysis_dir/compare_s23b_member_psf_injection.py" \
    --psf-dir "$run_root/injection" --output "$run_root/comparison"
printf '%s\n' COMPLETE >"$run_root/status"
date -u '+%Y-%m-%dT%H:%M:%SZ'
