#!/usr/bin/env bash
#SBATCH --job-name=icl-noise-full
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-full-noise-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-full-noise-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_icl_full_sample_coadd_covariance.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_noise_correlation_full_2026-11-03_01 \
  --patch-size 64 --max-patches 96 --max-lag 4 --min-valid-fraction 0.80
