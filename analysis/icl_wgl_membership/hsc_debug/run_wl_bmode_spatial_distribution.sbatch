#!/usr/bin/env bash
#SBATCH --job-name=wl-bmode-spatial
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-bmode-spatial-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-bmode-spatial-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_wl_bmode_spatial_distribution.py \
  --rows /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03 \
  --preflight /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01 \
  --plan-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_plan_2026-10-18_04 \
  --grid-width-deg 1.0 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_spatial_diagnostic_2026-10-19_01/summary.json
