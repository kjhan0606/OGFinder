#!/bin/bash
# CPU diagnostic over existing HSC images; the GPU allocation satisfies the
# cluster's Slurm-node policy but this analysis does not use CUDA.
#SBATCH --job-name=icl-psf-stars
#SBATCH --partition=a10
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_gaia_full_sample_2026-10-20_01_logs/slurm-%j.out

set -euo pipefail
module purge
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

repo=/home/kjhan/BACKUP/Astrafek/OGFinder
sample=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13
python=/home/kjhan/miniconda3/bin/python3
cd "$repo"
exec "$python" analysis/icl_wgl_membership/hsc_debug/diagnose_empirical_psf.py \
  --sample "$sample" \
  --gaia-work "$sample/psf_gaia_full_sample_match_2026-10-20_03" \
  --detection-dir "$sample/satellite_detection_psf_2026-10-18_01" \
  --work "$sample/empirical_psf_gaia_full_sample_2026-10-20_01"
