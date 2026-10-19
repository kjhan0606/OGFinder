#!/usr/bin/env bash
#SBATCH --job-name=icl-sat-lag-cov
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-sat-lag-cov-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-sat-lag-cov-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_satellite_icl_shared_pixel_covariance.py \
  --lag-covariance-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_noise_correlation_2026-10-19_04 \
  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lag_covariance_2026-10-19_01
