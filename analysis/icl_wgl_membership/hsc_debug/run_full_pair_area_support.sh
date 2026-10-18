#!/bin/bash
#SBATCH --job-name=icl-pair-area
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=01:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/pair-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/pair-%j.err

set -euo pipefail
module purge
module load gnu13
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
/home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/analyze_full_pair_area_support.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --extraction /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_extraction_2026-10-15_01 \
  --cluster-masks /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_cluster_masks_2026-10-16_01 \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_full_pair_area_support_2026-10-16_01
