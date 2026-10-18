#!/usr/bin/env bash
#SBATCH --job-name=hsc-wl-bmode-sens
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_sensitivity_%j.log
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_sensitivity_%j.err

set -euo pipefail
module purge
source /home/kjhan/miniconda3/etc/profile.d/conda.sh
conda activate base
cd /home/kjhan/BACKUP/Astrafek/OGFinder

python3 -u analysis/icl_wgl_membership/hsc_debug/run_wl_bmode_polarity_sensitivity.py \
  --rows /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03 \
  --preflight /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01 \
  --nz-sacc /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_polarity_sensitivity_2026-10-19_01/stack.json
