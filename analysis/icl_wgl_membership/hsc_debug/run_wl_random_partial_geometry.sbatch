#!/usr/bin/env bash
#SBATCH --job-name=wl-random-geom-pilot
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/wl-random-geom-pilot-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/wl-random-geom-pilot-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/screen_wl_random_completed_batches.py \
  --plan /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_coordinate_plan_2026-10-19_01 \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_geometry_partial_2026-10-20_04
