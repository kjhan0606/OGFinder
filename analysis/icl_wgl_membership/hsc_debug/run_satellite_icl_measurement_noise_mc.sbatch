#!/usr/bin/env bash
#SBATCH --job-name=icl-sat-noise-mc
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --gres=gpu:1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-sat-noise-mc-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-sat-noise-mc-%j.err

set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/propagate_satellite_icl_measurement_noise.py \
  --association-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_outer_icl_association_perbin_2026-10-18_03 \
  --covariance-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_target_control_covariance_2026-10-19_03 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_measurement_noise_mc_2026-10-19_01 \
  --replicas 10000 --seed 20261019
