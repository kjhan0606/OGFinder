#!/usr/bin/env bash
#SBATCH --job-name=icl-pair-sky
#SBATCH --partition=a10
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2-00:00:00
#SBATCH --output=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-pair-sky-%j.out
#SBATCH --error=/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl-pair-sky-%j.err
set -euo pipefail
module purge
module load gnu14
unset PYTHONHOME PYTHONPATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 MPLBACKEND=Agg
cd /home/kjhan/BACKUP/Astrafek/OGFinder
exec /home/kjhan/miniconda3/bin/python3 \
  analysis/icl_wgl_membership/hsc_debug/audit_icl_target_control_outer_sky.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_target_control_outer_sky_2026-11-03_01 \
  --replicas 5000 --seed 20261103
