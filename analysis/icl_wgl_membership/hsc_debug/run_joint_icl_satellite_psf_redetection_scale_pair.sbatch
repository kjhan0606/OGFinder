#!/usr/bin/env bash
#SBATCH --job-name=icl-scale-pair
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-scale-pair-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-scale-pair-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
runner=analysis/icl_wgl_membership/hsc_debug/run_joint_icl_satellite_batch.py
python=/home/kjhan/miniconda3/bin/python3
data=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
psf_plan="$data/icl_psf_picker_injection_plan_2026-10-19_02"
psf_files="$data/icl_psf_picker_injection_2026-10-19_01"

exec flock /tmp/astrafex_icl_redetection_scale_pair.lock bash -c '
  set -euo pipefail
  "$1" "$2" --output-dir "$3/joint_icl_satellite_psf_redetection_100kpc_audited_2026-10-20_01" \
    --diffuse-amplitude 1.0 --diffuse-scale-kpc 100.0 \
    --satellite-amplitude 1.0 --satellite-scale-kpc 15.0 \
    --rerun-detection-masks --unresolved-bcg-policy mask_all \
    --psf-plan-dir "$4" --psf-download-dir "$5"
  "$1" "$2" --output-dir "$3/joint_icl_satellite_psf_redetection_500kpc_audited_2026-10-20_01" \
    --diffuse-amplitude 1.0 --diffuse-scale-kpc 500.0 \
    --satellite-amplitude 1.0 --satellite-scale-kpc 15.0 \
    --rerun-detection-masks --unresolved-bcg-policy mask_all \
    --psf-plan-dir "$4" --psf-download-dir "$5"
' _ "$python" "$runner" "$data" "$psf_plan" "$psf_files"
