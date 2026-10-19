#!/usr/bin/env bash
#SBATCH --job-name=icl-2984-recheck
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=00:30:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-2984-recheck-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-2984-recheck-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
outdir=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_500kpc_frame_2984_recheck_2026-10-20_01
test ! -e "$outdir"
mkdir "$outdir"
exec flock /tmp/astrafex_icl_redetection_injection.lock /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/joint_icl_satellite_injection.py \
  --cluster-id 2984 \
  --output "$outdir/cluster_2984.json" \
  --diffuse-amplitude 1.0 \
  --diffuse-scale-kpc 500.0 \
  --satellite-amplitude 1.0 \
  --satellite-scale-kpc 15.0 \
  --rerun-detection-masks \
  --unresolved-bcg-policy mask_all \
  --psf-plan-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_injection_plan_2026-10-19_02 \
  --psf-download-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_injection_2026-10-19_01
