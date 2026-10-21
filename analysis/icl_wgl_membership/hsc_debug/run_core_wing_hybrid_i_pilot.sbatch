#!/usr/bin/env bash
#SBATCH --job-name=icl-psf-hybrid-pilot
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-hybrid-pilot-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_measurement_logs_2026-10-18_01/icl-psf-hybrid-pilot-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
/home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/run_core_wing_hybrid_grid.py \
  --filter i --cluster-id 1056 \
  --output-root /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_pilot_cluster1056_2026-10-21_01
