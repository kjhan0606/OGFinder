#!/bin/bash
#SBATCH --job-name=icl-mask-support
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --gres=gpu:1
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/mask-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/mask-%j.err

set -euo pipefail
module purge
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
/home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/prepare_full_cluster_mask_support.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_cluster_masks_2026-10-16_01
