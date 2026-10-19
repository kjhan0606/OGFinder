#!/usr/bin/env bash
#SBATCH --job-name=wl-tract-infl
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --time=03:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-tract-influence-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl-tract-influence-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/run_wl_hsc_tract_region_influence.py \
  --rows /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03 \
  --preflight /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01 \
  --plan-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_plan_2026-10-18_04 \
  --lens-manifest /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json \
  --patch-qa /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/catalogues_2026-10-13_01/patch_qa.csv \
  --nz-sacc /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc \
  --masked-value 0 \
  --adjacency-deg 2.0 \
  --radial-edges-mpc 0.1 0.2 0.4 0.7 1.0 1.5 2.0 3.0 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_hsc_tract_region_influence_2026-10-20_01
