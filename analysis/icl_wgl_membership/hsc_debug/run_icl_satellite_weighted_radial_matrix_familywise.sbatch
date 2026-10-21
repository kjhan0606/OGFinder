#!/usr/bin/env bash
#SBATCH --job-name=icl-wgt-radial
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_logs_2026-10-20_01/icl-wgt-radial-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_logs_2026-10-20_01/icl-wgt-radial-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_weighted_radial_matrix_familywise.py \
  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_familywise_2026-10-20_01 \
  --replicas 5000 \
  --seed 20261022
