#!/usr/bin/env bash
set -euo pipefail

host="$(hostname -s | tr '[:upper:]' '[:lower:]')"
[[ "$host" == lageunha* ]] || { echo "run on lageunha CPU; current host is $host" >&2; exit 2; }

repo=/home/kjhan/BACKUP/Astrafek/OGFinder
data=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
python=/home/kjhan/miniconda3/bin/python3
runner=$repo/analysis/icl_wgl_membership/hsc_debug/run_joint_icl_satellite_batch.py
auditor=$repo/analysis/icl_wgl_membership/hsc_debug/audit_joint_icl_satellite_psf_redetection.py
plan=$data/icl_psf_picker_injection_plan_2026-10-19_02
psfs=$data/icl_psf_picker_injection_2026-10-19_01
test_psf=$data/icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits

export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg PYTHONUNBUFFERED=1
cd "$repo"

for scale in 200 300; do
  run=$data/joint_icl_satellite_psf_redetection_${scale}kpc_2026-10-20_01
  audit=$run/psf_redetection_audit.json
  [[ ! -e "$run" ]] || { echo "refusing to overwrite $run" >&2; exit 3; }
  timeout 30m "$python" "$runner" \
    --output-dir "$run" \
    --diffuse-amplitude 1.0 --diffuse-scale-kpc "$scale" \
    --satellite-amplitude 1.0 --satellite-scale-kpc 15.0 \
    --rerun-detection-masks --unresolved-bcg-policy mask_all \
    --psf-plan-dir "$plan" --psf-download-dir "$psfs"
  timeout 5m "$python" "$auditor" \
    --output-dir "$run" --plan-dir "$plan" --download-dir "$psfs" \
    --test-center-psf "$test_psf" --output "$audit"
done
