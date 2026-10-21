#!/usr/bin/env bash
#SBATCH --job-name=icl-psf-wings
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-wings-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-wings-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
output_dir="${OUTPUT_DIR:-/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_2026-10-20_01}"
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_outer_field_wings.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --output "$output_dir" \
  --bootstrap-draws 10000 \
  --seed 20261020
