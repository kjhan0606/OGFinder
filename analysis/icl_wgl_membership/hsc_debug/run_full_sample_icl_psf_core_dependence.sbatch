#!/usr/bin/env bash
#SBATCH --job-name=icl-psf-mismatch
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-psf-mismatch-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-psf-mismatch-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
output_dir="${OUTPUT_DIR:-/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_core_dependence_full_sample_2026-10-21_01}"
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/analyze_full_sample_icl_psf_core_dependence.py \
  --output-dir "$output_dir" --replicas 5000 --seed 20261021
