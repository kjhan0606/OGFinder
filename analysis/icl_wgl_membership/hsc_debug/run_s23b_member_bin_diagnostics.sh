#!/usr/bin/env bash
# Local CPU rerun of approved images, then offline historical comparison.
set -euo pipefail
if [[ $# -ne 1 ]]; then
    echo 'Usage: run_s23b_member_bin_diagnostics.sh NEW_RUN_DIRECTORY' >&2
    exit 2
fi
case "$(hostname | tr '[:upper:]' '[:lower:]')" in
    lageunha*) ;;
    *) echo 'Run this CPU analysis on lageunha.' >&2; exit 2 ;;
esac
run_root=$1
if [[ "$run_root" != /* ]]; then
    echo 'Use an absolute new run directory.' >&2
    exit 2
fi
mkdir -- "$run_root" # refuse existing paths, including prior/active runs
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
    sha256sum "$analysis_dir/measure_s23b_camira_member_outer_profiles.py" \
        "$analysis_dir/audit_s23b_member_bin_diagnostics.py" \
        "$analysis_dir/measure_satellite_outer_profiles.py" \
        "$analysis_dir/prepare_detection_psf.py" \
        "$analysis_dir/member_profile_checkpoint.py" "${BASH_SOURCE[0]}"
} >"$run_root/provenance.txt"
while [[ ! -d "$run_root/measurement" ]]; do
    resume_args=()
    if [[ -f "$run_root/measurement.checkpoints/configuration.json" ]]; then
        resume_args=(--resume)
    fi
    "$python_bin" "$analysis_dir/measure_s23b_camira_member_outer_profiles.py" \
        --all-eligible --max-new-frames 100 "${resume_args[@]}" --output "$run_root/measurement"
done
"$python_bin" "$analysis_dir/audit_s23b_member_bin_diagnostics.py" \
    --measurement "$run_root/measurement" --output "$run_root/audit"
printf '%s\n' COMPLETE >"$run_root/status"
date -u '+%Y-%m-%dT%H:%M:%SZ'
